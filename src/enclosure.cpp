/***************************************************************************
 *            enclosure.cpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Evolution.
 *
 *  Ariadne Evolution is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Evolution is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Evolution. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file src/enclosure.cpp
 *  @brief Implement enclosure evolution operations.
 */

#include "evolution/enclosure.hpp"

#include <iomanip>
#include <stdexcept>

#include "algebra/algebra.hpp"
#include "algebra/differential.hpp"
#include "algebra/expansion.inl.hpp"
#include "algebra/matrix.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/polynomial.hpp"
#include "algebra/sweeper.hpp"
#include "algebra/vector.hpp"
#include "drawing/drawer.hpp"
#include "drawing/drawer_manager.hpp"
#include "enclosure_recorded_draw.hpp"
#include "function/constraint.hpp"
#include "function/formula.hpp"
#include "function/function.hpp"
#include "function/function_model.hpp"
#include "function/functional.hpp"
#include "function/procedure.hpp"
#include "function/taylor_function.hpp"
#include "geometry/box.hpp"
#include "labelled/figure.hpp"
#include "labelled/graphics_interface.hpp"
#include "labelled/set.hpp"
#include "logging/logging.hpp"
#include "logging/progress_indicator.hpp"
#include "numeric/numeric.hpp"
#include "numeric/operators.hpp"
#include "paving/binary_tree_paving.hpp"
#include "paving/grid.hpp"
#include "paving/paver.hpp"
#include "paving/paving_interface.hpp"
#include "rendering/figure.hpp"
#include "rendering/rendering_manager.hpp"
#include "set/affine_set.hpp"
#include "set/function_set.hpp"
#include "symbolic/space.hpp"
#include "utility/exceptions.hpp"
#include "utility/macros.hpp"
#include "utility/stopwatch.hpp"

#include "constraint-solving/constraint_solver.hpp"
#include "expression-set/expression_set.hpp"
#include "nonlinear-programming/nonlinear_programming.hpp"

namespace Ariadne {

template<class T>
inline StringType str(const T& t) {
    StringStream ss;
    ss << t;
    return ss.str();
}

// This overload was removed from the extracted Set module. Keep the original
// Real-interval rounding and error-enclosure algorithm local to dynamics.
static Pair<Interval<FloatDP>, FloatDPError> make_domain(const Interval<Real>& ivl) {
    FloatDPBounds dlb(ivl.lower_bound(), dp);
    FloatDPBounds dub(ivl.upper_bound(), dp);
    FloatDPApproximation dla(ivl.lower_bound(), dp);
    FloatDPApproximation dua(ivl.upper_bound(), dp);
    // Keep the full FloatDP precision of the approximated endpoints.  Reducing
    // them to Float32 here creates an avoidable, step-size-independent error
    // which is later transferred to the enclosure remainder.
    FloatDP dl(dla.raw());
    FloatDP du(dua.raw());
    FloatDPError e = cast_positive(max(max(dub.upper() - du, du - dub.lower()), max(dlb.upper() - dl, dl - dlb.lower())));
    return make_pair(make_interval(dl, du), e);
}

namespace {

ValidatedVectorMultivariateFunctionPatch make_identity(const RealBox& bx, const EnclosureConfiguration& configuration) {
    ExactBoxType dom(bx.dimension());
    Vector<FloatDPError> errs(bx.dimension(), dp);

    for (SizeType i = 0; i != bx.dimension(); ++i) {
        make_lpair(dom[i], errs[i]) = make_domain(bx[i]);
    }

    ValidatedVectorMultivariateFunctionPatch res = configuration.function_factory().create_identity(dom);
    for (SizeType i = 0; i != bx.dimension(); ++i) {
        res[i] = res[i] + FloatDPBounds(-errs[i], +errs[i]);
    }

    return res;
}

} // namespace

OutputStream& operator<<(OutputStream& os, const EnclosureVariableKind& vk) {
    switch (vk) {
        case EnclosureVariableKind::INITIAL:
            return os << "x";
        case EnclosureVariableKind::TEMPORAL:
            return os << "t";
        case EnclosureVariableKind::PARAMETER:
            return os << "a";
        case EnclosureVariableKind::INPUT:
            return os << "u";
        case EnclosureVariableKind::NOISE:
            return os << "v";
        case EnclosureVariableKind::ERROR:
            return os << "e";
        case EnclosureVariableKind::UNKNOWN:
        default:
            return os << "s";
    }
}

List<Identifier> canonical_variable_names(const List<EnclosureVariableKind>& vks) {
    std::map<EnclosureVariableKind, Nat> counts;
    List<Identifier> result;
    for (SizeType i = 0; i != vks.size(); ++i) {
        result.append(str(vks[i]) + str(counts[vks[i]]++));
    }
    return result;
}

inline Pair<ValidatedScalarMultivariateFunctionPatch, ValidatedScalarMultivariateFunctionPatch> split(const ValidatedScalarMultivariateFunctionPatch& f, SizeType k) {
    Pair<ExactBoxType, ExactBoxType> domains = split(f.domain(), k);
    return make_pair(restriction(f, domains.first), restriction(f, domains.second));
}

inline Pair<ValidatedVectorMultivariateFunctionPatch, ValidatedVectorMultivariateFunctionPatch> split(const ValidatedVectorMultivariateFunctionPatch& f, SizeType k) {
    Pair<ExactBoxType, ExactBoxType> domains = split(f.domain(), k);
    return make_pair(restriction(f, domains.first), restriction(f, domains.second));
}

EnclosureConfiguration::EnclosureConfiguration(ValidatedFunctionPatchFactory function_factory, SizeType reconditioning_num_blocks)
    : _function_factory(function_factory), _reconditioning_num_blocks(reconditioning_num_blocks) { ARIADNE_PRECONDITION(reconditioning_num_blocks > 0u); }

OutputStream& operator<<(OutputStream& os, const EnclosureConfiguration& ec) {
    return os << "EnclosureConfiguration( function_factory=" << ec._function_factory
              << ", reconditioning_num_blocks=" << ec._reconditioning_num_blocks << ")";
}

const EnclosureConfiguration&
    Enclosure::configuration() const {
    return this->_configuration;
}

/*

// TODO: Make more efficient
inline Void assign_all_but_last(MultiIndex& r, const MultiIndex& a) {
    for(SizeType i=0; i!=r.size(); ++i) { r[i]=a[i]; }
}

// FIXME: What if solving for constraint leaves domain?
Void Enclosure::_solve_zero_constraints() {
    for(List<ValidatedScalarMultivariateFunctionPatch>::Iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ) {
        const ExactBoxType& domain=this->domain();
        const ValidatedTaylorModelDP& model=iter->model();
        const SizeType k=model.argument_size()-1u;
        ValidatedTaylorModelDP zeroth_order(k,this->sweeper());
        ValidatedTaylorModelDP first_order(k,this->sweeper());
        Bool is_zeroth_order=true;
        Bool is_first_order=true;
        MultiIndex r(k);
        // Try linear approach in last coefficient
        for(ValidatedTaylorModelDP::ConstIterator tmiter=model.begin(); tmiter!=model.end(); ++tmiter) {
            if(tmiter->index()[k]==0) {
                assign_all_but_last(r,tmiter->index());
                zeroth_order.expansion().append(r,tmiter->coefficient());
            } else if(tmiter->index()[k]==1) {
                is_zeroth_order=false;
                assign_all_but_last(r,tmiter->index());
                first_order.expansion().append(r,tmiter->coefficient());
            } else {
                is_first_order=false; break;
            }
        }
        if(is_first_order && !is_zeroth_order) {
            const ExactBoxType new_domain=project(domain,range(0,k));
            ValidatedTaylorModelDP substitution_model=-zeroth_order/first_order;
            this->_state_function=this->configuration().function_factory().create(new_domain,Ariadne::substitute(this->_state_function.models(),k,substitution_model));
            for(List<ValidatedScalarMultivariateFunctionPatch>::Iterator constraint_iter=this->_negative_constraints.begin();
                    constraint_iter!=this->_negative_constraints.end(); ++constraint_iter) {
                ValidatedScalarMultivariateFunctionPatch& constraint=*constraint_iter;
                constraint=this->configuration().function_factory().create(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            for(List<ValidatedScalarMultivariateFunctionPatch>::Iterator constraint_iter=this->_zero_constraints.begin();
                    constraint_iter!=this->_zero_constraints.end(); ++constraint_iter) {
                ValidatedScalarMultivariateFunctionPatch& constraint=*constraint_iter;
                constraint=this->configuration().function_factory().create(new_domain,Ariadne::substitute(constraint.model(),k,substitution_model));
            }
            // Since we are using an std::vector, assign Iterator to next element
            iter=this->_zero_constraints.erase(iter);
        } else {
            ARIADNE_WARN("No method for solving constraint "<<*iter<<" currently implemented.");
            ++iter;
        }
    }
}
*/

Enclosure::Enclosure()
    : _domain(), _auxiliary_mapping(), _state_function(), _time_function(), _dwell_time_function(), _reduced_domain(), _variable_kinds(), _configuration(ValidatedFunctionPatchFactory(nullptr)) {
}

Enclosure::Enclosure(const EnclosureConfiguration& config)
    : _domain(), _auxiliary_mapping(), _state_function(), _time_function(), _dwell_time_function(), _reduced_domain(), _variable_kinds(), _configuration(config) {
}

Enclosure* Enclosure::clone() const {
    return new Enclosure(*this);
}

Enclosure::Enclosure(const RealBox& box, const EnclosureConfiguration& configuration)
    : Enclosure(BoundedConstraintSet(box), configuration) {
}

Enclosure::Enclosure(const BoundedConstraintSet& set, const EnclosureConfiguration& configuration)
    : _configuration(configuration) {
    this->_state_function = make_identity(set.domain(), configuration);
    this->_domain = this->_state_function.domain();
    this->_time_function = this->configuration().function_factory().create_zero(this->domain());
    this->_dwell_time_function = this->configuration().function_factory().create_zero(this->domain());
    this->_auxiliary_mapping = EffectiveVectorMultivariateFunction(0u, EuclideanDomain(this->_state_function.result_size()));
    this->_variable_kinds = List<EnclosureVariableKind>(this->_domain.size(), EnclosureVariableKind::INITIAL);
    for (SizeType i = 0; i != set.number_of_constraints(); ++i) {
        this->new_state_constraint(set.constraint(i));
    }
    this->_reduced_domain = this->_domain;
}

Enclosure::Enclosure(const ExactBoxType& box, const EnclosureConfiguration& configuration)
    : _configuration(configuration) {
    // Ensure domain elements have nonempty radius
    this->_domain = box;

    this->_variable_kinds = List<EnclosureVariableKind>(this->_domain.dimension(), EnclosureVariableKind::INITIAL);

    this->_state_function = this->configuration().function_factory().create_identity(this->_domain);
    this->_time_function = this->configuration().function_factory().create_zero(this->_domain);
    this->_dwell_time_function = this->configuration().function_factory().create_zero(this->domain());
    this->_auxiliary_mapping = EffectiveVectorMultivariateFunction(0u, EuclideanDomain(this->_state_function.result_size()));
    this->_reduced_domain = this->_domain;
}

Enclosure::Enclosure(const ExactBoxType& domain, const ValidatedVectorMultivariateFunction& function, const EnclosureConfiguration& configuration)
    : Enclosure(domain, function, List<ValidatedConstraint>(), configuration) {
}

Enclosure::Enclosure(const ExactBoxType& domain, const ValidatedVectorMultivariateFunction& function, const List<ValidatedConstraint>& constraints, const EnclosureConfiguration& configuration)
    : Enclosure(domain, function, configuration.function_factory().create_zero(domain), constraints, configuration) {
}

Enclosure::Enclosure(const ExactBoxType& domain, const ValidatedVectorMultivariateFunction& state_function, const ValidatedScalarMultivariateFunction& time_function, const List<ValidatedConstraint>& constraints, const EnclosureConfiguration& configuration)
    : _configuration(configuration) {
    ARIADNE_PRECONDITION_MSG(domain.size() == state_function.argument_size(), "domain=" << domain << ", state_function=" << state_function);
    ARIADNE_PRECONDITION_MSG(domain.size() == time_function.argument_size(), "domain=" << domain << ", time_function=" << time_function);
    // Generic function domains are EuclideanDomain(argument_size), checked above.
    this->_domain = domain;
    this->_variable_kinds = List<EnclosureVariableKind>(this->_domain.size(), EnclosureVariableKind::INITIAL);

    this->_state_function = this->configuration().function_factory().create(this->_domain, state_function);
    this->_time_function = this->configuration().function_factory().create(this->_domain, time_function);
    this->_dwell_time_function = this->configuration().function_factory().create_zero(this->domain());
    this->_auxiliary_mapping = EffectiveVectorMultivariateFunction(0u, EuclideanDomain(this->_state_function.result_size()));

    for (SizeType i = 0; i != constraints.size(); ++i) {
        ARIADNE_PRECONDITION_MSG(domain.size() == constraints[i].function().argument_size(), "domain=" << domain << ", constraint=" << constraints[i]);
        this->new_parameter_constraint(constraints[i]);
    }

    this->_reduced_domain = domain;
    this->reduce();
}

// Returns true if the entire set is positive; false if entire set is negative
ValidatedKleenean Enclosure::satisfies(ValidatedScalarMultivariateFunction constraint) const {
    UpperIntervalType constraint_range = apply(constraint, this->codomain());
    if (definitely(constraint_range.upper_bound() < 0)) {
        return false;
    }
    if (definitely(constraint_range.lower_bound() > 0)) {
        return true;
    }
    return ValidatedKleenean(LogicalValue::INDETERMINATE);
}

/*
Void Enclosure::substitute(SizeType j, ValidatedScalarMultivariateFunctionPatch v)
{
    ARIADNE_ASSERT_MSG(v.argument_size()+1u==this->number_of_parameters(),
                       "number_of_parameters="<<this->number_of_parameters()<<", variable="<<v);
                       this->_state_function = Ariadne::substitute(this->_state_function,j,v);
                       for(List<ValidatedScalarMultivariateFunctionPatch>::Iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }
                       for(List<ValidatedScalarMultivariateFunctionPatch>::Iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
                           *iter = Ariadne::substitute(*iter,j,v);
                       }

}

Void Enclosure::substitute(SizeType j, FloatDP c)
{
    this->_state_function = Ariadne::partial_evaluate(this->_state_function,j,c);
    for(List<ValidatedScalarMultivariateFunctionPatch>::Iterator iter=this->_negative_constraints.begin(); iter!=this->_negative_constraints.end(); ++iter) {
        *iter = Ariadne::partial_evaluate(*iter,j,c);
    }
    for(List<ValidatedScalarMultivariateFunctionPatch>::Iterator iter=this->_zero_constraints.begin(); iter!=this->_zero_constraints.end(); ++iter) {
        *iter = Ariadne::partial_evaluate(*iter,j,c);
    }
}
*/

const List<EnclosureVariableKind>& Enclosure::variable_kinds() const {
    return this->_variable_kinds;
}

const EffectiveVectorMultivariateFunction& Enclosure::auxiliary_mapping() const {
    return this->_auxiliary_mapping;
}

EnclosureFunctionArguments Enclosure::auxiliary_arguments() const { return _auxiliary_arguments; }
SizeType Enclosure::auxiliary_parameter_count() const { return _auxiliary_parameter_count; }

Bool Enclosure::has_constant_parameter_binding() const { return _constant_parameters.has_value(); }
RealSpace Enclosure::constant_parameter_space() const {
    return _constant_parameters ? RealSpace(_constant_parameters->names) : RealSpace{};
}
SizeType Enclosure::constant_parameter_count() const {
    return _constant_parameters ? _constant_parameters->names.size() : 0u;
}
const ExactBoxType& Enclosure::constant_parameter_domain() const {
    if (!_constant_parameters)
        throw std::logic_error("Enclosure has no named physical parameter binding");
    return _constant_parameters->domain;
}

Void Enclosure::_constrain_constant_parameters() {
    for (SizeType i = 0u; i != constant_parameter_count(); ++i) {
        const auto interval = _constant_parameters->domain[i];
        _constraints.append(ValidatedConstraintModel(interval.lower_bound(),
                                                     _constant_parameters->function.get(i), interval.upper_bound()));
    }
}

Bool Enclosure::_expose_constant_parameter_errors() {
    if (!_constant_parameters)
        return false;
    List<SizeType> indices;
    for (SizeType i = 0u; i != constant_parameter_count(); ++i)
        if (possibly(_constant_parameters->function.get(i).error() > 0))
            indices.append(i);
    if (indices.empty())
        return false;
    Enclosure extended = *this;
    const SizeType offset = number_of_parameters();
    const ValidatedFunctionPatchFactory persistent_factory{new TaylorFunctionFactory(TrivialSweeper<FloatDP>(dp))};
    extended.adjoin_parameters(ExactBoxType(indices.size(), ExactIntervalType(-1, 1)), EnclosureVariableKind::PARAMETER);
    for (SizeType j = 0u; j != indices.size(); ++j) {
        auto component = extended._constant_parameters->function.get(indices[j]);
        const FloatDP error = cast_exact(component.error().get(double_precision));
        const auto* taylor = dynamic_cast<const ValidatedScalarMultivariateTaylorFunctionModelDP*>(&component.reference());
        if (!taylor)
            throw std::invalid_argument("Persistent physical parameter errors currently require Taylor models");
        auto model = taylor->model();
        model.set_sweeper(TrivialSweeper<FloatDP>(dp));
        model.clobber();
        component = ValidatedScalarMultivariateTaylorFunctionModelDP(extended.parameter_domain(), model);
        component = component + persistent_factory.create_coordinate(extended.parameter_domain(), offset + j) * FloatDPBounds(error);
        if (possibly(component.error() > 0))
            throw std::runtime_error("Function factory cannot expose constant parameter error as an exact persistent coordinate");
        extended._constant_parameters->function.set(indices[j], component);
    }
    extended._constrain_constant_parameters();
    *this = std::move(extended);
    return true;
}

Void Enclosure::bind_constant_parameters(const RealSpace& space, const ExactBoxType& domain,
                                         const ValidatedVectorMultivariateFunctionPatch& function) {
    if (_constant_parameters)
        throw std::invalid_argument("Enclosure already has a physical parameter binding");
    if (space.dimension() == 0u || space.dimension() != domain.dimension() ||
        space.dimension() != function.result_size() ||
        (enclosure_has_parameter_arguments(_auxiliary_arguments) &&
         _auxiliary_parameter_count > space.dimension()))
        throw std::invalid_argument("Physical parameter names, domain and map dimensions do not agree");
    ARIADNE_PRECONDITION(Ariadne::subset(parameter_domain(), function.domain()));
    Enclosure bound = *this;
    bound._constant_parameters = ConstantParameterBinding{space.variable_names(), domain, restriction(function, parameter_domain())};
    if (const auto* taylor = dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP*>(&bound._constant_parameters->function.reference())) {
        auto persistent = *taylor;
        persistent.set_properties(TrivialSweeper<FloatDP>(dp));
        bound._constant_parameters->function = persistent;
    }
    if (!bound._expose_constant_parameter_errors())
        bound._constrain_constant_parameters();
    *this = std::move(bound);
}

Void Enclosure::adjoin_constant_parameters(const RealSpace& space, const ExactBoxType& domain) {
    if (_constant_parameters || space.dimension() != domain.dimension())
        throw std::invalid_argument("Cannot adjoin incompatible or duplicate physical parameters");
    if (space.dimension() == 0u)
        return;
    Enclosure extended = *this;
    extended.adjoin_parameters(domain, EnclosureVariableKind::PARAMETER);
    extended.bind_constant_parameters(space, domain, extended.constant_parameter_function(space.dimension()));
    *this = std::move(extended);
}

List<Bool> Enclosure::_protected_parameter_axes() const {
    List<Bool> result(number_of_parameters(), false);
    for (SizeType i = 0u; i != number_of_parameters(); ++i)
        result[i] = _variable_kinds[i] == EnclosureVariableKind::PARAMETER;
    // An auxiliary a(x,u_current) must retain the exact current INPUT
    // coordinates across Kuhn reconditioning. Historical INPUT axes may
    // instead be summarized as validated state error.
    if (enclosure_has_input_arguments(_auxiliary_arguments)) {
        SizeType remaining = _auxiliary_parameter_count;
        for (SizeType i = number_of_parameters(); i != 0u && remaining != 0u; --i) {
            if (_variable_kinds[i - 1u] == EnclosureVariableKind::INPUT) {
                result[i - 1u] = true;
                --remaining;
            }
        }
    }
    if (_constant_parameters) {
        const auto* models = dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP*>(&_constant_parameters->function.reference());
        if (!models)
            return List<Bool>(number_of_parameters(), true);
        for (SizeType k = 0u; k != models->result_size(); ++k)
            for (auto term = models->models()[k].begin(); term != models->models()[k].end(); ++term)
                for (SizeType i = 0u; i != number_of_parameters(); ++i)
                    if (term->coefficient() != FloatDP(0, dp) && term->index()[i] != 0)
                        result[i] = true;
    }
    return result;
}

SizeType Enclosure::number_of_protected_parameters() const {
    const auto axes = _protected_parameter_axes();
    return static_cast<SizeType>(std::count(axes.begin(), axes.end(), true));
}

ValidatedVectorMultivariateFunctionPatch Enclosure::constant_parameter_function(SizeType count) const {
    if (_constant_parameters) {
        const SizeType parameter_dimension = constant_parameter_count();
        if (count > parameter_dimension)
            throw std::logic_error("Physical parameter binding has insufficient components");
        if (count == parameter_dimension)
            return _constant_parameters->function;
        auto selected = configuration().function_factory().create_zeros(count, parameter_domain());
        for (SizeType i = 0u; i != count; ++i)
            selected.set(i, _constant_parameters->function.get(parameter_dimension - count + i));
        return selected;
    }
    Array<SizeType> indices(count);
    SizeType remaining = count;
    for (SizeType i = number_of_parameters(); i != 0u && remaining != 0u; --i) {
        if (_variable_kinds[i - 1u] == EnclosureVariableKind::PARAMETER)
            indices[--remaining] = i - 1u;
    }
    if (remaining != 0u)
        throw std::logic_error("Enclosure has insufficient trajectory-constant parameter coordinates");
    const auto& factory = configuration().function_factory();
    auto result = factory.create_zeros(count, parameter_domain());
    for (SizeType i = 0u; i != count; ++i)
        result[i] = factory.create_coordinate(parameter_domain(), indices[i]);
    return result;
}

ValidatedVectorMultivariateFunctionPatch Enclosure::input_parameter_function(SizeType count) const {
    Array<SizeType> indices(count);
    SizeType remaining = count;
    for (SizeType i = number_of_parameters(); i != 0u && remaining != 0u; --i) {
        if (_variable_kinds[i - 1u] == EnclosureVariableKind::INPUT)
            indices[--remaining] = i - 1u;
    }
    if (remaining != 0u)
        throw std::invalid_argument("Enclosure has insufficient current input coordinates");
    const auto& factory = configuration().function_factory();
    auto values = factory.create_zeros(count, parameter_domain());
    for (SizeType j = 0u; j != count; ++j)
        values[j] = factory.create_coordinate(parameter_domain(), indices[j]);
    return values;
}

ValidatedVectorMultivariateFunctionPatch Enclosure::argument_function(EnclosureFunctionArguments arguments, SizeType parameter_count) const {
    const bool selected_layout = enclosure_has_parameter_arguments(arguments) ||
                                 enclosure_has_input_arguments(arguments);
    if (selected_layout != (parameter_count != 0u))
        throw std::invalid_argument("Enclosure function parameter count does not match its explicit argument layout");
    switch (arguments) {
        case EnclosureFunctionArguments::STATE:
            return _state_function;
        case EnclosureFunctionArguments::STATE_TIME:
            return join(_state_function, _time_function);
        case EnclosureFunctionArguments::STATE_PARAMETERS:
            return join(_state_function, constant_parameter_function(parameter_count));
        case EnclosureFunctionArguments::STATE_TIME_PARAMETERS:
            return join(join(_state_function, _time_function), constant_parameter_function(parameter_count));
        case EnclosureFunctionArguments::STATE_INPUTS:
            return join(_state_function, input_parameter_function(parameter_count));
        default:
            throw std::invalid_argument("Unknown enclosure function argument layout");
    }
}

Void Enclosure::set_auxiliary_mapping(const EffectiveVectorMultivariateFunction& aux, EnclosureFunctionArguments arguments, SizeType parameter_count) {
    ARIADNE_PRECONDITION(argument_function(arguments, parameter_count).result_size() == aux.argument_size());
    this->_auxiliary_mapping = aux;
    this->_auxiliary_arguments = arguments;
    this->_auxiliary_parameter_count = parameter_count;
}

Void Enclosure::_unchecked_new_variable(ExactIntervalType ivl, EnclosureVariableKind vk) {
    ValidatedScalarMultivariateFunctionPatch variable_function = this->configuration().function_factory().create_coordinate({ivl}, 0);
    this->_domain = product(this->_domain, ivl);
    this->_reduced_domain = product(this->_reduced_domain, ivl);
    this->_state_function = combine(this->_state_function, variable_function);
    this->_time_function = embed(this->_time_function, ivl);
    this->_dwell_time_function = embed(this->_dwell_time_function, ivl);
    if (_constant_parameters)
        _constant_parameters->function = embed(_constant_parameters->function, ivl);

    for (SizeType i = 0; i != this->_constraints.size(); ++i) {
        ValidatedConstraintModel& constraint = this->_constraints[i];
        constraint.set_function(embed(constraint.function(), ivl));
    }
    this->_variable_kinds.append(vk);
}

Void Enclosure::set_time_function(const ValidatedScalarMultivariateFunctionPatch& time) {
    ARIADNE_PRECONDITION(time.argument_size() == this->number_of_parameters());
    ARIADNE_PRECONDITION(Ariadne::subset(this->parameter_domain(), time.domain()));
    this->_time_function = restriction(time, this->parameter_domain());
}

Void Enclosure::adjoin_parameters(const ExactBoxType& parameters, EnclosureVariableKind kind) {
    if (parameters.dimension() == 0u) {
        return;
    }
    Enclosure extended = *this;
    extended._domain = product(this->_domain, parameters);
    extended._reduced_domain = product(this->_reduced_domain, parameters);
    extended._state_function = embed(this->_state_function, parameters);
    extended._time_function = embed(this->_time_function, parameters);
    extended._dwell_time_function = embed(this->_dwell_time_function, parameters);
    if (extended._constant_parameters)
        extended._constant_parameters->function = embed(this->_constant_parameters->function, parameters);
    for (auto& constraint: extended._constraints) {
        constraint.set_function(embed(constraint.function(), parameters));
    }
    extended._variable_kinds.concatenate(
        List<EnclosureVariableKind>(parameters.dimension(), kind));
    *this = std::move(extended);
}

Void Enclosure::clear_time() {
    this->_time_function = 0;
    this->_dwell_time_function = 0;
}

Void Enclosure::apply_map(ValidatedVectorMultivariateFunction map) {
    ARIADNE_PRECONDITION_MSG(map.argument_size() == this->state_dimension(), "state_dimension=" << this->state_dimension() << ", map=" << map);
    const SizeType nt = enclosure_has_time_argument(_auxiliary_arguments) ? 1u : 0u;
    ARIADNE_PRECONDITION(map.result_size() + nt + _auxiliary_parameter_count == this->_auxiliary_mapping.argument_size());
    this->_state_function = compose(map, this->_state_function);
    this->_dwell_time_function = this->configuration().function_factory().create_zero(this->domain());
}

Void Enclosure::apply_map(ValidatedVectorMultivariateFunction map, EffectiveVectorMultivariateFunction aux_map,
                          EnclosureFunctionArguments source_arguments, EnclosureFunctionArguments auxiliary_arguments) {
    const auto arguments = argument_function(source_arguments, enclosure_has_parameter_arguments(source_arguments) ? constant_parameter_count() : 0u);
    ARIADNE_PRECONDITION_MSG(map.argument_size() == arguments.result_size(), "arguments=" << arguments << ", map=" << map);
    const SizeType nt = enclosure_has_time_argument(auxiliary_arguments) ? 1u : 0u;
    const SizeType np = enclosure_has_parameter_arguments(auxiliary_arguments) ? constant_parameter_count() : 0u;
    ARIADNE_PRECONDITION(auxiliary_arguments == EnclosureFunctionArguments::STATE || auxiliary_arguments == EnclosureFunctionArguments::STATE_TIME || enclosure_has_parameter_arguments(auxiliary_arguments));
    if (enclosure_has_parameter_arguments(auxiliary_arguments) && np == 0u)
        throw std::invalid_argument("Parameter-dependent reset outputs require an explicit physical binding");
    ARIADNE_PRECONDITION_MSG(aux_map.argument_size() == map.result_size() + nt + np, "map=" << map << ", aux_map=" << aux_map);
    this->_state_function = compose(map, arguments);
    this->_dwell_time_function = this->configuration().function_factory().create_zero(this->domain());
    this->_auxiliary_mapping = aux_map;
    this->_auxiliary_arguments = auxiliary_arguments;
    this->_auxiliary_parameter_count = np;
}

Void Enclosure::apply_discrete_time_map_step(ValidatedVectorMultivariateFunction map) {
    this->apply_map(map);
    this->_time_function = this->_time_function + StepSizeType(1, 0u);
}

Void Enclosure::apply_fixed_evolve_step(ValidatedVectorMultivariateFunction flow, StepSizeType time) {
    ARIADNE_PRECONDITION(flow.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION_MSG(flow.argument_size() == this->state_dimension() + 1u, "state_dimension=" << this->state_dimension() << ", flow=" << flow);
    try {
        ValidatedVectorMultivariateFunctionPatch flow_model = dynamic_handle_cast<const ValidatedVectorMultivariateFunctionPatch>(flow);
        ValidatedVectorMultivariateFunctionPatch flow_step_model = partial_evaluate(flow_model, flow_model.argument_size() - 1u, time);
        this->_state_function = compose(flow_step_model, this->_state_function);
    } catch (const std::bad_cast&) {
        ValidatedScalarMultivariateFunctionPatch evolve_time_function = this->configuration().function_factory().create_constant(this->domain(), time);
        this->_state_function = compose(flow, join(this->_state_function, evolve_time_function));
    }
    this->_time_function = this->_time_function + time;
    this->_dwell_time_function = this->_dwell_time_function + time;
}

Void Enclosure::apply_parameterised_fixed_evolve_step(
    ValidatedVectorMultivariateFunctionPatch state, StepSizeType time) {
    ARIADNE_PRECONDITION(state.argument_size() == this->number_of_parameters());
    ARIADNE_PRECONDITION(state.domain() == this->parameter_domain());
    ARIADNE_PRECONDITION(state.result_size() == this->state_dimension());
    this->_state_function = std::move(state);
    this->_time_function = this->_time_function + time;
    this->_dwell_time_function = this->_dwell_time_function + time;
}

Void Enclosure::apply_space_evolve_step(ValidatedVectorMultivariateFunction flow, ValidatedScalarMultivariateFunction time) {
    ARIADNE_PRECONDITION(flow.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION_MSG(flow.argument_size() == this->state_dimension() + 1u, "state_dimension=" << this->state_dimension() << ", flow=" << flow);
    ARIADNE_PRECONDITION_MSG(time.argument_size() == this->state_dimension(), "state_dimension=" << this->state_dimension() << ", time=" << time);
    ValidatedScalarMultivariateFunctionPatch evolve_time_function = compose(time, this->_state_function);
    this->_state_function = compose(flow, join(this->_state_function, evolve_time_function));
    this->_time_function = this->_time_function + evolve_time_function;
    this->_dwell_time_function = this->_dwell_time_function + evolve_time_function;
}

Void Enclosure::apply_spacetime_evolve_step(ValidatedVectorMultivariateFunction flow, ValidatedScalarMultivariateFunction time) {
    ARIADNE_PRECONDITION(flow.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION_MSG(flow.argument_size() == this->state_dimension() + 1u, "state_dimension=" << this->state_dimension() << ", flow=" << flow);
    ARIADNE_PRECONDITION_MSG(time.argument_size() == this->state_dimension() + 1u, "state_dimension=" << this->state_dimension() << ", time=" << time);
    ValidatedScalarMultivariateFunctionPatch evolve_time_function = compose(time, join(this->_state_function, this->_time_function));
    this->_state_function = compose(flow, join(this->_state_function, evolve_time_function));
    this->_time_function = this->_time_function + evolve_time_function;
    this->_dwell_time_function = this->_dwell_time_function + evolve_time_function;
}

Void Enclosure::apply_parameter_evolve_step(ValidatedVectorMultivariateFunction flow, ValidatedScalarMultivariateFunction time) {
    ARIADNE_PRECONDITION(flow.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION_MSG(flow.argument_size() == this->state_dimension() + 1u, "state_dimension=" << this->state_dimension() << ", flow=" << flow);
    ARIADNE_PRECONDITION_MSG(time.argument_size() == this->number_of_parameters(), "number_of_parameters=" << this->number_of_parameters() << ", time=" << time << ", enclosure=" << *this);
    this->_state_function = compose(flow, join(this->_state_function, this->configuration().function_factory().create(this->_state_function.domain(), time)));
    this->_time_function = this->_time_function + time;
    this->_dwell_time_function = this->_dwell_time_function + time;
}

Void Enclosure::apply_finishing_parameter_evolve_step(ValidatedVectorMultivariateFunction flow, ValidatedScalarMultivariateFunction finishing_time) {
    ARIADNE_PRECONDITION(flow.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION_MSG(flow.argument_size() == this->state_dimension() + 1u, "state_dimension=" << this->state_dimension() << ", flow=" << flow);
    ARIADNE_PRECONDITION_MSG(finishing_time.argument_size() == this->number_of_parameters(), "number_of_parameters=" << this->number_of_parameters() << ", finishing_time=" << finishing_time);
    ValidatedScalarMultivariateFunctionPatch omega = this->configuration().function_factory().create(this->domain(), finishing_time);
    this->_state_function = compose(flow, join(this->_state_function, omega - this->_time_function));
    this->_dwell_time_function = this->_dwell_time_function + (omega - this->_time_function);
    this->_time_function = omega;
}

Void Enclosure::apply_full_reach_step(ValidatedVectorMultivariateFunctionPatch phi) {
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_PRECONDITION(phi.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION(phi.argument_size() == this->state_dimension() + 1);
    FloatDP h = phi.domain()[phi.result_size()].upper_bound();
    ValidatedScalarMultivariateFunctionPatch elps = this->configuration().function_factory().create_constant(this->domain(), h);
    this->apply_parameter_reach_step(phi, elps);
}

Void Enclosure::apply_parameterised_full_reach_step(
    ValidatedVectorMultivariateFunctionPatch flowpipe) {
    ARIADNE_PRECONDITION(flowpipe.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION(flowpipe.argument_size() == this->number_of_parameters() + 1u);

    const ExactIntervalType time_domain =
        flowpipe.domain()[flowpipe.domain().dimension() - 1u];
    const ExactBoxType old_parameter_domain = this->parameter_domain();
    ARIADNE_PRECONDITION(flowpipe.domain() == product(old_parameter_domain, time_domain));

    this->_unchecked_new_variable(time_domain, EnclosureVariableKind::TEMPORAL);
    this->_state_function = std::move(flowpipe);

    ValidatedScalarMultivariateFunctionPatch time_step_function =
        this->configuration().function_factory().create_coordinate(
            this->domain(), this->domain().size() - 1u);
    this->_time_function = this->_time_function + time_step_function;
    this->_dwell_time_function = this->_dwell_time_function + time_step_function;
}

Void Enclosure::apply_spacetime_reach_step(ValidatedVectorMultivariateFunctionPatch phi, ValidatedScalarMultivariateFunction elps) {
    ARIADNE_PRECONDITION(phi.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION(phi.argument_size() == this->state_dimension() + 1);
    ARIADNE_PRECONDITION(elps.argument_size() == this->state_dimension() + 1);
    this->apply_parameter_reach_step(phi, compose(elps, join(this->state_function(), this->time_function())));
}

Void Enclosure::apply_parameter_reach_step(ValidatedVectorMultivariateFunctionPatch phi, ValidatedScalarMultivariateFunction elps) {
    // xi'(s,t) = phi(xi(s),t) for t in [0,h] with constraint t<=eps(s) where range(eps) in [0,h]
    // tau'(s) = tau(s)+t
    ARIADNE_PRECONDITION(phi.result_size() == this->state_dimension());
    ARIADNE_PRECONDITION(phi.argument_size() == this->state_dimension() + 1);
    ARIADNE_PRECONDITION(elps.argument_size() == this->number_of_parameters());
    FloatDP h = phi.domain()[phi.result_size()].upper_bound();
    ExactBoxType parameter_domain = this->parameter_domain();
    ExactIntervalType time_domain = ExactIntervalType(FloatDP(0, dp), h);
    ValidatedScalarMultivariateFunctionPatch time_function = this->configuration().function_factory().create_coordinate({time_domain}, 0u);
    this->_unchecked_new_variable(time_domain, EnclosureVariableKind::TEMPORAL);
    // The appended temporal coordinate supplies the final flow argument.
    this->apply_map(phi);
    ExactBoxType new_domain = this->parameter_domain();
    ValidatedScalarMultivariateFunctionPatch time_step_function = this->configuration().function_factory().create_coordinate(new_domain, new_domain.size() - 1u);
    this->_time_function = this->_time_function + time_step_function;
    this->_dwell_time_function = this->_dwell_time_function + time_step_function;
    if (phi.domain()[phi.result_size()].lower_bound() < time_domain.upper_bound()) {
        auto restricted_elps = this->configuration().function_factory().create(parameter_domain, elps);
        this->new_negative_parameter_constraint(time_step_function - embed(restricted_elps, time_domain));
    }
}

Void Enclosure::new_state_time_bound(ValidatedScalarMultivariateFunction gamma) {
    ARIADNE_PRECONDITION(gamma.argument_size() == this->state_dimension());
    this->_constraints.append(this->time_function() <= compose(gamma, this->state_function()));
}

Void Enclosure::new_state_constraint(ValidatedConstraint constraint) {
    ARIADNE_PRECONDITION(constraint.function().argument_size() == this->state_dimension());
    ValidatedNumber lower_bound = constraint.lower_bound();
    ValidatedScalarMultivariateFunctionPatch composed_function_model = compose(constraint.function(), this->_state_function);
    ValidatedNumber upper_bound = constraint.upper_bound();
    this->_constraints.append(ValidatedConstraintModel(lower_bound, composed_function_model, upper_bound));
}

Void Enclosure::new_state_time_constraint(ValidatedConstraint constraint) {
    ARIADNE_PRECONDITION(constraint.function().argument_size() == this->state_dimension() + 1u);
    ValidatedNumber lower_bound = constraint.lower_bound();
    ValidatedScalarMultivariateFunctionPatch composed_function_model = compose(constraint.function(), join(this->_state_function, this->_time_function));
    ValidatedNumber upper_bound = constraint.upper_bound();
    this->_constraints.append(ValidatedConstraintModel(lower_bound, composed_function_model, upper_bound));
}

Void Enclosure::new_parameter_constraint(ValidatedConstraint constraint) {
    ARIADNE_PRECONDITION(constraint.function().argument_size() == this->number_of_parameters());
    ValidatedNumber lower_bound = constraint.lower_bound();
    ValidatedScalarMultivariateFunctionPatch function_model = this->configuration().function_factory().create(this->domain(), constraint.function());
    ValidatedNumber upper_bound = constraint.upper_bound();
    this->_constraints.append(ValidatedConstraintModel(lower_bound, function_model, upper_bound));
}

Void Enclosure::new_positive_state_constraint(ValidatedScalarMultivariateFunction constraint_function) {
    ARIADNE_PRECONDITION_MSG(constraint_function.argument_size() == this->state_dimension(), "state_dimension=" << this->state_dimension() << ", constraint_function=" << constraint_function);
    ExactNumber zero = 0;
    this->_constraints.append(compose(constraint_function, this->state_function()) >= zero);
}

Void Enclosure::new_negative_state_constraint(ValidatedScalarMultivariateFunction constraint_function) {
    ARIADNE_PRECONDITION_MSG(constraint_function.argument_size() == this->state_dimension(), "state_dimension=" << this->state_dimension() << ", constraint_function=" << constraint_function);
    ExactNumber zero = 0;
    this->_constraints.append(compose(constraint_function, this->state_function()) <= zero);
}

Void Enclosure::new_zero_state_constraint(ValidatedScalarMultivariateFunction constraint_function) {
    ARIADNE_PRECONDITION_MSG(constraint_function.argument_size() == this->state_dimension(), "state_dimension=" << this->state_dimension() << ", constraint_function=" << constraint_function);
    ExactNumber zero = 0;
    this->_constraints.append(compose(constraint_function, this->state_function()) == zero);
}

Void Enclosure::new_negative_parameter_constraint(ValidatedScalarMultivariateFunction constraint_function) {
    ARIADNE_PRECONDITION_MSG(constraint_function.argument_size() == this->domain().size(), "domain=" << this->domain() << ", constraint_function=" << constraint_function);
    ExactNumber zero = 0;
    this->_constraints.append(this->configuration().function_factory().create(this->domain(), constraint_function) <= zero);
}

Void Enclosure::new_zero_parameter_constraint(ValidatedScalarMultivariateFunction constraint_function) {
    ARIADNE_PRECONDITION_MSG(constraint_function.argument_size() == this->domain().size(), "domain=" << this->domain() << ", constraint_function=" << constraint_function);
    ExactNumber zero = 0;
    this->_constraints.append(this->configuration().function_factory().create(this->domain(), constraint_function) == zero);
}

ExactBoxType Enclosure::domain() const {
    return this->_domain;
}

ExactBoxType Enclosure::parameter_domain() const {
    return this->_domain;
}

ExactBoxType Enclosure::reduced_domain() const {
    return this->_reduced_domain;
}

ExactBoxType Enclosure::codomain() const {
    return cast_exact_box(widen(this->_state_function.range()));
}

const ValidatedVectorMultivariateFunctionPatch Enclosure::state_time_auxiliary_function() const {
    return join(join(this->state_function(), this->time_function()), this->auxiliary_function());
}

const ValidatedVectorMultivariateFunctionPatch Enclosure::state_auxiliary_function() const {
    return join(this->state_function(), this->auxiliary_function());
}

const ValidatedVectorMultivariateFunctionPatch& Enclosure::state_function() const {
    return this->_state_function;
}

const ValidatedScalarMultivariateFunctionPatch& Enclosure::time_function() const {
    return this->_time_function;
}

const ValidatedScalarMultivariateFunctionPatch& Enclosure::dwell_time_function() const {
    return this->_dwell_time_function;
}

const ValidatedVectorMultivariateFunctionPatch Enclosure::auxiliary_function() const {
    static_assert(Same<decltype(this->_auxiliary_mapping), EffectiveVectorMultivariateFunction>);
    static_assert(Same<decltype(this->_state_function), ValidatedVectorMultivariateFunctionPatch>);
    static_assert(not Convertible<ValidatedVectorMultivariateFunction, EffectiveVectorMultivariateFunction>);
    static_assert(not Convertible<ValidatedVectorMultivariateFunctionPatch, EffectiveVectorMultivariateFunction>);
    static_assert(not Convertible<EffectiveVectorMultivariateFunction, ValidatedVectorMultivariateFunctionPatch>);
    return compose(this->_auxiliary_mapping, this->argument_function(this->_auxiliary_arguments, this->_auxiliary_parameter_count));
}

SizeType Enclosure::number_of_constraints() const {
    return this->_constraints.size();
}

const List<ValidatedConstraintModel>& Enclosure::constraint_models() const {
    return this->_constraints;
}

const List<ValidatedConstraint> Enclosure::constraints() const {
    List<ValidatedConstraint> result;
    for (SizeType i = 0; i != this->_constraints.size(); ++i) {
        result.append(ValidatedConstraint(this->_constraints[i].lower_bound(), this->_constraints[i].function(), this->_constraints[i].upper_bound()));
    }
    return result;
}

const ValidatedConstraintModel& Enclosure::constraint(SizeType i) const {
    return this->_constraints[i];
}

const ValidatedVectorMultivariateFunctionPatch Enclosure::constraint_function() const {
    ValidatedVectorMultivariateFunctionPatch g = this->configuration().function_factory().create_zeros(this->number_of_constraints(), this->domain());
    for (SizeType i = 0; i != this->number_of_constraints(); ++i) {
        g.set(i, this->constraint(i).function());
    }
    return g;
}

const ExactBoxType Enclosure::constraint_bounds() const {
    ExactBoxType c(this->number_of_constraints());
    for (SizeType i = 0; i != this->number_of_constraints(); ++i) {
        c[i] = this->constraint(i).bounds();
    }
    return c;
}

DimensionType Enclosure::dimension() const {
    return this->_state_function.result_size() + this->_auxiliary_mapping.result_size();
}

DimensionType Enclosure::state_dimension() const {
    return this->_state_function.result_size();
}

SizeType Enclosure::number_of_parameters() const {
    return this->_state_function.argument_size();
}

UpperBoxType Enclosure::bounding_box() const {
    return this->_state_function.codomain().bounding_box();
}

FloatDPError Enclosure::radius() const {
    return cast_positive(this->bounding_box().radius());
}

FloatDPPoint Enclosure::centre() const {
    return cast_exact(this->bounding_box().centre());
}

ValidatedKleenean
    Enclosure::satisfies(ValidatedConstraint c) const {
    Enclosure copy = *this;
    copy.new_state_constraint(c);
    if (definitely(copy.is_empty())) {
        return false;
    } else {
        return ValidatedKleenean(LogicalValue::INDETERMINATE);
    }
}

ValidatedLowerKleenean Enclosure::is_bounded() const {
    return this->domain().is_bounded() || ValidatedKleenean(LogicalValue::INDETERMINATE);
}

ValidatedLowerKleenean Enclosure::is_empty() const {
    if (definitely(this->_reduced_domain.is_empty())) {
        return true;
    }
    if (this->_constraints.empty()) {
        return this->domain().is_empty();
    }
    // Constraint models can be appended by constructors and reconditioning.
    // Contract every constrained enclosure rather than relying on a stale flag.
    // Repeating reduction can tighten coupled constraints, but once a pass
    // leaves the exact parameter box unchanged, further identical passes
    // cannot contract it. Avoid rebuilding the solver for those passes.
    if (this->number_of_parameters() > 0u) {
        for (SizeType pass = 0u; pass != 3u; ++pass) {
            const ExactBoxType previous_domain = this->_reduced_domain;
            this->reduce();
            if (this->_reduced_domain.is_empty() ||
                same(previous_domain, this->_reduced_domain))
                break;
        }
    }
    // Contraction may make the parameter domain empty. Avoid evaluating models
    // on that empty box before returning the established infeasibility.
    if (this->_reduced_domain.is_empty()) {
        return true;
    }

    for (SizeType i = 0; i != this->_constraints.size(); ++i) {
        // A zero-dimensional domain has one point and constant models. No
        // interval contractor or vector of evaluation arguments is needed.
        UpperIntervalType constraint_range = this->number_of_parameters() == 0u
                                               ? this->_constraints[i].function().range()
                                               : Ariadne::apply(this->_constraints[i].function(), this->_reduced_domain);
        if (definitely(disjoint(constraint_range, this->_constraints[i].bounds()))) {
            if (this->_reduced_domain.size() > 0) {
                this->_reduced_domain[0] = ExactIntervalType(1, -1);
            }
            return true;
        }
    }
    return ValidatedKleenean(LogicalValue::INDETERMINATE);
}

ValidatedLowerKleenean Enclosure::inside(const ExactBoxType& bx) const {
    return Ariadne::subset(Ariadne::apply(this->_state_function, this->_reduced_domain), bx);
}

ValidatedLowerKleenean Enclosure::subset(const ExactBoxType& bx) const {
    this->reduce();

    return ValidatedLowerKleenean(Ariadne::subset(Ariadne::apply(this->_state_function, this->_reduced_domain), bx)) || ValidatedKleenean(LogicalValue::INDETERMINATE);
}

ValidatedLowerKleenean Enclosure::separated(const ExactBoxType& bx) const {
    ARIADNE_PRECONDITION_MSG(this->state_dimension() == bx.dimension(), "Enclosure::subset(ExactBoxType): self=" << *this << ", box=" << bx);
    List<ValidatedConstraint> constraints = this->constraints();
    ConstraintSolver contractor = ConstraintSolver();
    UpperBoxType reduced_domain(this->_reduced_domain);
    contractor.reduce(reduced_domain, constraints);
    this->_reduced_domain = cast_exact_box(reduced_domain);

    if (_reduced_domain.is_empty()) {
        return true;
    }

    const ExactBoxType test_domain = this->_reduced_domain;
    for (SizeType i = 0; i != bx.dimension(); ++i) {
        // FIXME: Conversion should be automatic
        ValidatedScalarMultivariateFunction fi(static_cast<const ValidatedScalarMultivariateFunction::Interface&>(this->_state_function[i]));
        constraints.append(fi >= bx[i].lower_bound());
        constraints.append(fi <= bx[i].upper_bound());
    }
    return !contractor.feasible(test_domain, constraints).first;
}

Void Enclosure::reduce() const {
    List<ValidatedConstraint> constraints = this->constraints();
    ConstraintSolver contractor = ConstraintSolver();
    UpperBoxType reduced_domain(this->_reduced_domain);
    contractor.reduce(reduced_domain, constraints);
    this->_reduced_domain = cast_exact_box(reduced_domain);

    // The hull contractor only intersects each coordinate with candidates via
    // max(old_lower,candidate) / min(old_upper,candidate). With valid input
    // endpoints, std::max/min retain the old endpoint for a NaN candidate too;
    // contraction therefore preserves valid endpoints and domain inclusion.

    /*
        // Remove redundant constraints
        SizeType j=0;
        List<ValidatedScalarMultivariateFunctionPatch>& mutable_constraints=const_cast<List<ValidatedScalarMultivariateFunctionPatch>&>(this->_negative_constraints);
        for(SizeType i=0; i!=mutable_constraints.size(); ++i) {
            if(mutable_constraints[i](this->_reduced_domain).upper_bound()<0.0) { redundant_constraints.append(i); }
            else { if(i>j) { mutable_constraints[j]=mutable_constraints[j]; } ++j; }
        }
        mutable_constraints.resize(j);
    */
}

// The former generic overload lived in geometry/function_set.cpp. Preserve
// its Taylor-model dispatch here while linking to the extracted Set implementation.
Matrix<FloatDPError> nonlinearities_zeroth_order(
    const ValidatedVectorMultivariateTaylorFunctionModelDP&, const ExactBoxType&);

static Matrix<FloatDPError> nonlinearities_zeroth_order(
    const ValidatedVectorMultivariateFunction& f, const ExactBoxType& dom) {
    return nonlinearities_zeroth_order(
        dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP&>(*f.raw_pointer()), dom);
}
Pair<SizeType, FloatDPError> nonlinearity_index_and_error(const ValidatedVectorMultivariateFunction& function, const ExactBoxType& domain);
Pair<SizeType, FloatDPError> lipschitz_index_and_error(const ValidatedVectorMultivariateFunction& function, const ExactBoxType& domain);

Pair<Enclosure, Enclosure>
    Enclosure::split_zeroth_order() const {
    return this->split(this->splitting_index_zeroth_order());
}

List<ExactBoxType>
    Enclosure::splitting_subdomains_zeroth_order() const {
    List<ExactBoxType> result;
    SizeType k = this->splitting_index_zeroth_order();
    if (k == this->number_of_parameters()) {
        result.append(this->_reduced_domain);
    } else {
        Pair<ExactBoxType, ExactBoxType> subdomains = this->_reduced_domain.split(this->splitting_index_zeroth_order());
        result.append(subdomains.first);
        result.append(subdomains.second);
    }
    return result;
}

SizeType
    Enclosure::splitting_index_zeroth_order() const {
    Matrix<UpperIntervalType> jacobian = Ariadne::jacobian_range(this->state_function(), cast_vector(this->reduced_domain()));

    // Compute the column of the matrix which has the norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    SizeType jmax = this->number_of_parameters();
    FloatDPError max_column_norm(0u, dp);
    for (SizeType j = 0; j != this->number_of_parameters(); ++j) {
        FloatDPError column_norm(0u, dp);
        for (SizeType i = 0; i != this->state_dimension(); ++i) {
            column_norm += mag(jacobian[i][j]);
        }
        column_norm *= this->reduced_domain()[j].radius();
        if (column_norm.raw() > max_column_norm.raw()) {
            max_column_norm = column_norm;
            jmax = j;
        }
    }

    return jmax;
}

Pair<Enclosure, Enclosure>
    Enclosure::split_first_order() const {
    Matrix<FloatDPError> nonlinearities = Ariadne::nonlinearities_zeroth_order(this->_state_function, this->_reduced_domain);

    // Compute the row of the nonlinearities Array which has the highest norm
    // i.e. the highest sum of $mag(a_ij)$ where mag([l,u])=max(|l|,|u|)
    SizeType jmax_in_row_imax = nonlinearities.column_size();
    FloatDPError max_row_sum(0u, dp);
    for (SizeType i = 0; i != nonlinearities.row_size(); ++i) {
        SizeType jmax = nonlinearities.column_size();
        FloatDPError row_sum(0u, dp);
        FloatDPError max_mag_j_in_i(0u, dp);
        for (SizeType j = 0; j != nonlinearities.column_size(); ++j) {
            row_sum += mag(nonlinearities[i][j]);
            if (mag(nonlinearities[i][j]).raw() > max_mag_j_in_i.raw()) {
                jmax = j;
                max_mag_j_in_i = mag(nonlinearities[i][j]);
            }
        }
        if (row_sum.raw() > max_row_sum.raw()) {
            max_row_sum = row_sum;
            jmax_in_row_imax = jmax;
        }
    }

    if (jmax_in_row_imax == nonlinearities.column_size()) {
        ARIADNE_THROW(std::runtime_error, "split_first_order", "No need to split");
    }

    return this->split(jmax_in_row_imax);
}

Pair<Enclosure, Enclosure>
    Enclosure::split() const {
    return this->split_zeroth_order();
}

Pair<Enclosure, Enclosure>
    Enclosure::split(SizeType d) const {
    ARIADNE_PRECONDITION(d < this->number_of_parameters());
    ExactBoxType subdomain1, subdomain2;
    make_lpair(subdomain1, subdomain2) = Ariadne::split(this->_state_function.domain(), d);

    ValidatedVectorMultivariateFunctionPatch function1, function2;
    make_lpair(function1, function2) = Ariadne::split(this->_state_function, d);

    Pair<Enclosure, Enclosure>
        result = make_pair(Enclosure(function1.domain(), function1, this->configuration()),
                           Enclosure(function2.domain(), function2, this->configuration()));
    Enclosure& result1 = result.first;
    Enclosure& result2 = result.second;

    ValidatedScalarMultivariateFunctionPatch constraint_function1, constraint_function2;
    for (List<ValidatedConstraintModel>::ConstIterator iter = this->_constraints.begin();
         iter != this->_constraints.end(); ++iter) {
        const ValidatedConstraintModel& constraint = *iter;
        make_lpair(constraint_function1, constraint_function2) = Ariadne::split(constraint.function(), d);
        result1._constraints.append(ValidatedConstraintModel(constraint.lower_bound(), constraint_function1, constraint.upper_bound()));
        result2._constraints.append(ValidatedConstraintModel(constraint.lower_bound(), constraint_function2, constraint.upper_bound()));
    }

    ValidatedScalarMultivariateFunctionPatch time_function1, time_function2;
    make_lpair(time_function1, time_function2) = Ariadne::split(this->_time_function, d);
    result1._time_function = time_function1;
    result2._time_function = time_function2;
    ValidatedScalarMultivariateFunctionPatch dwell_time_function1, dwell_time_function2;
    make_lpair(dwell_time_function1, dwell_time_function2) = Ariadne::split(this->_dwell_time_function, d);
    result1._dwell_time_function = dwell_time_function1;
    result2._dwell_time_function = dwell_time_function2;

    result1._auxiliary_mapping = this->_auxiliary_mapping;
    result1._auxiliary_arguments = this->_auxiliary_arguments;
    result1._auxiliary_parameter_count = this->_auxiliary_parameter_count;
    result2._auxiliary_mapping = this->_auxiliary_mapping;
    result2._auxiliary_arguments = this->_auxiliary_arguments;
    result2._auxiliary_parameter_count = this->_auxiliary_parameter_count;
    result1._variable_kinds = this->_variable_kinds;
    result2._variable_kinds = this->_variable_kinds;
    if (_constant_parameters) {
        result1._constant_parameters = _constant_parameters;
        result2._constant_parameters = _constant_parameters;
        result1._constant_parameters->function = restriction(_constant_parameters->function, result1.parameter_domain());
        result2._constant_parameters->function = restriction(_constant_parameters->function, result2.parameter_domain());
        result1._expose_constant_parameter_errors();
        result2._expose_constant_parameter_errors();
    }

    return result;
}

ValidatedConstrainedImageSet Enclosure::state_set() const {
    return ValidatedConstrainedImageSet(this->domain(), this->state_function(), this->constraints());
}

ValidatedConstrainedImageSet Enclosure::state_auxiliary_set() const {
    return ValidatedConstrainedImageSet(this->domain(), this->state_auxiliary_function(), this->constraints());
}

ValidatedConstrainedImageSet Enclosure::state_time_auxiliary_set() const {
    return ValidatedConstrainedImageSet(this->domain(), this->state_time_auxiliary_function(), this->constraints());
}

TaylorModel<ValidatedTag, FloatDP> recondition(const TaylorModel<ValidatedTag, FloatDP>& tm, Array<SizeType>& discarded_variables, SizeType number_of_error_variables, SizeType index_of_error);
TaylorModel<ValidatedTag, FloatDP> recondition(const TaylorModel<ValidatedTag, FloatDP>& tm, Array<SizeType>& discarded_variables, SizeType number_of_error_variables);

Void Enclosure::recondition() {
    this->kuhn_recondition();
    // Kuhn reduction can turn discarded clock dependencies into uniform error.
    // Expose those errors last, so subsequent time-dependent steps retain them.
    this->uniform_error_recondition();
}

Void Enclosure::
    uniform_error_recondition() {
    this->_expose_constant_parameter_errors();
    const ExactDouble MAXIMUM_ERROR = cast_exact(std::numeric_limits<double>::epsilon() * 1024);
    const SizeType old_number_of_parameters = this->number_of_parameters();
    const SizeType state_dimension = this->state_dimension();
    List<SizeType> large_error_indices;
    for (SizeType i = 0; i != state_dimension; ++i) {
        if (possibly(this->_state_function.get(i).error() > MAXIMUM_ERROR))
            large_error_indices.append(i);
    }
    if (possibly(this->_time_function.error() > MAXIMUM_ERROR))
        large_error_indices.append(state_dimension);
    if (possibly(this->_dwell_time_function.error() > MAXIMUM_ERROR))
        large_error_indices.append(state_dimension + 1u);
    if (large_error_indices.empty())
        return;

    ExactBoxType error_domains(large_error_indices.size(), ExactIntervalType(-1, +1));
    this->_domain = product(this->_domain, error_domains);
    this->_reduced_domain = product(this->_reduced_domain, error_domains);
    this->_state_function = embed(this->_state_function, error_domains);
    this->_time_function = embed(this->_time_function, error_domains);
    this->_dwell_time_function = embed(this->_dwell_time_function, error_domains);
    if (_constant_parameters)
        _constant_parameters->function = embed(_constant_parameters->function, error_domains);
    for (auto& constraint: this->_constraints) {
        constraint.function() = embed(constraint.function(), error_domains);
    }
    SizeType parameter = old_number_of_parameters;
    for (SizeType index: large_error_indices) {
        ValidatedScalarMultivariateFunctionPatch component = index < state_dimension
                                                               ? this->_state_function.get(index)
                                                               : (index == state_dimension ? this->_time_function : this->_dwell_time_function);
        const FloatDP error = cast_exact(component.error().get(double_precision));
        component.clobber();
        component = component + this->configuration().function_factory().create_coordinate(this->_domain, parameter) * FloatDPBounds(error);
        if (index < state_dimension)
            this->_state_function[index] = component;
        else if (index == state_dimension)
            this->_time_function = component;
        else
            this->_dwell_time_function = component;
        ++parameter;
    }
    this->_variable_kinds.concatenate(List<EnclosureVariableKind>(large_error_indices.size(), EnclosureVariableKind::ERROR));
}

TaylorModel<ValidatedTag, FloatDP> recondition(const TaylorModel<ValidatedTag, FloatDP>& tm, Array<SizeType>& discarded_variables, SizeType number_of_error_variables, SizeType index_of_error) {
    // The two local callers supply sorted indices chosen from [0,argument_size)
    // and an error index in [0,number_of_error_variables]. They only reduce
    // when at least one parameter is discarded.
    const SizeType number_of_variables = tm.argument_size();
    const SizeType number_of_discarded_variables = discarded_variables.size();
    const SizeType number_of_kept_variables = number_of_variables - number_of_discarded_variables;

    // Make an Array of the variables to be kept
    Array<SizeType> kept_variables = complement(number_of_variables, discarded_variables);

    // Construct result and reserve memory
    TaylorModel<ValidatedTag, FloatDP> r(number_of_kept_variables + number_of_error_variables, tm.sweeper());
    r.expansion().reserve(tm.number_of_nonzeros() + 1u);
    MultiIndex ra(number_of_kept_variables + number_of_error_variables);

    // Accumulate bounds in an Error object, not through a FloatDP alias.
    const Bool embed_error_in_variable = index_of_error < number_of_error_variables;
    if (embed_error_in_variable) {
        ra[number_of_kept_variables + index_of_error] = 1;
        r.expansion().append(ra, FloatDP(dp));
        ra[number_of_kept_variables + index_of_error] = 0;
    }
    FloatDPError error = tm.error();

    for (TaylorModel<ValidatedTag, FloatDP>::ConstIterator iter = tm.begin(); iter != tm.end(); ++iter) {
        UniformConstReference<MultiIndex> xa = iter->index();
        UniformConstReference<FloatDP> xv = iter->coefficient();
        Bool keep = true;
        for (SizeType k = 0; k != number_of_discarded_variables; ++k) {
            if (xa[discarded_variables[k]] != 0) {
                error += mag(xv);
                keep = false;
                break;
            }
        }
        if (keep) {
            for (SizeType k = 0; k != number_of_kept_variables; ++k) {
                ra[k] = xa[kept_variables[k]];
            }
            r.expansion().append(ra, xv);
        }
    }
    if (embed_error_in_variable) {
        // Store the upward-rounded bound as an exact coefficient.
        r.expansion().begin()->coefficient() = error.raw();
    } else {
        r.error() = error;
    }
    return r;
}

TaylorModel<ValidatedTag, FloatDP> recondition(const TaylorModel<ValidatedTag, FloatDP>& tm, Array<SizeType>& discarded_variables, SizeType number_of_error_variables) {
    return recondition(tm, discarded_variables, number_of_error_variables, number_of_error_variables);
}

Void Enclosure::kuhn_recondition() {
    this->_expose_constant_parameter_errors();
    if (!dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP*>(&this->state_function().reference())) {
        ARIADNE_WARN("Cannot Kuhn reduce an Enclosure which is not given by TaylorFunctions.");
        // Uniform error reconditioning also supports alternative function models.
        this->uniform_error_recondition();
        return;
    }

    const auto protected_axes = _protected_parameter_axes();
    const SizeType fixed_parameters = static_cast<SizeType>(std::count(protected_axes.begin(), protected_axes.end(), true));
    const SizeType number_of_kept_parameters = std::max(
        fixed_parameters, (this->configuration().reconditioning_num_blocks() - 1) * this->state_dimension());
    const SizeType number_of_error_parameters = this->state_dimension();

    if (this->number_of_parameters() <= number_of_kept_parameters) {
        this->uniform_error_recondition();
        return;
    }
    const SizeType number_of_discarded_parameters = this->number_of_parameters() - number_of_kept_parameters;

    const ValidatedVectorMultivariateTaylorFunctionModelDP& function = dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP&>(this->state_function().reference());
    const Vector<ValidatedTaylorModelDP>& models = function.models();
    Matrix<FloatDPApproximation> dependencies(this->state_dimension(), this->number_of_parameters(), dp);
    for (SizeType i = 0; i != dependencies.row_size(); ++i) {
        for (ValidatedTaylorModelDP::ConstIterator iter = models[i].begin(); iter != models[i].end(); ++iter) {
            for (SizeType j = 0; j != dependencies.column_size(); ++j) {
                if (iter->index()[j] != 0) {
                    dependencies[i][j] = dependencies[i][j] + abs(iter->coefficient());
                }
            }
        }
    }
    Array<Pair<FloatDP, SizeType>> column_max_dependencies(this->number_of_parameters(), make_pair(FloatDP(dp), 0u));
    for (SizeType j = 0; j != dependencies.column_size(); ++j) {
        column_max_dependencies[j] = make_pair(FloatDP(0.0_x, dp), SizeType(j));
        for (SizeType i = 0; i != dependencies.row_size(); ++i) {
            column_max_dependencies[j].first = max(column_max_dependencies[j].first, cast_exact(dependencies[i][j]));
        }
    }
    std::sort(column_max_dependencies.begin(), column_max_dependencies.end(),
              [&](const auto& left, const auto& right) {
                  const bool left_fixed = protected_axes[left.second];
                  const bool right_fixed = protected_axes[right.second];
                  return left_fixed != right_fixed ? left_fixed : left > right;
              });

    Array<SizeType> kept_parameters(number_of_kept_parameters);
    Array<SizeType> discarded_parameters(number_of_discarded_parameters);
    for (SizeType j = 0; j != number_of_kept_parameters; ++j) {
        kept_parameters[j] = column_max_dependencies[j].second;
    }
    for (SizeType j = 0; j != number_of_discarded_parameters; ++j) {
        discarded_parameters[j] = column_max_dependencies[number_of_kept_parameters + j].second;
    }
    std::sort(kept_parameters.begin(), kept_parameters.end());
    std::sort(discarded_parameters.begin(), discarded_parameters.end());

    Vector<ValidatedTaylorModelDP> new_models(models.size(), ValidatedTaylorModelDP(number_of_kept_parameters + number_of_error_parameters, function.properties()));
    for (SizeType i = 0; i != this->state_dimension(); ++i) {
        new_models[i] = Ariadne::recondition(models[i], discarded_parameters, number_of_error_parameters, i);
    }

    ExactBoxType new_domain(number_of_kept_parameters + number_of_error_parameters);
    ExactBoxType new_reduced_domain(number_of_kept_parameters + number_of_error_parameters);
    for (SizeType j = 0; j != number_of_kept_parameters; ++j) {
        new_domain[j] = this->parameter_domain()[kept_parameters[j]];
        new_reduced_domain[j] = this->reduced_domain()[kept_parameters[j]];
    }
    for (SizeType j = number_of_kept_parameters; j != number_of_kept_parameters + number_of_error_parameters; ++j) {
        new_domain[j] = ExactIntervalType(-1, +1);
        new_reduced_domain[j] = ExactIntervalType(-1, +1);
    }
    Enclosure new_set(new_domain, ValidatedVectorMultivariateTaylorFunctionModelDP(new_domain, new_models), this->configuration());
    for (SizeType i = 0; i != this->_constraints.size(); ++i) {
        const ValidatedConstraintModel& constraint = this->_constraints[i];
        const ValidatedScalarMultivariateTaylorFunctionModelDP& constraint_function = dynamic_cast<const ValidatedScalarMultivariateTaylorFunctionModelDP&>(constraint.function().reference());
        ValidatedScalarMultivariateFunctionPatch new_constraint_function = ValidatedScalarMultivariateTaylorFunctionModelDP(new_domain, Ariadne::recondition(constraint_function.model(), discarded_parameters, number_of_error_parameters));
        new_set._constraints.append(ValidatedConstraintModel(constraint.lower_bound(), new_constraint_function, constraint.upper_bound()));
    }
    const ValidatedScalarMultivariateTaylorFunctionModelDP& time = dynamic_cast<const ValidatedScalarMultivariateTaylorFunctionModelDP&>(this->_time_function.reference());
    new_set._time_function = ValidatedScalarMultivariateTaylorFunctionModelDP(new_domain, Ariadne::recondition(time.model(), discarded_parameters, number_of_error_parameters));
    const ValidatedScalarMultivariateTaylorFunctionModelDP& dwell_time = dynamic_cast<const ValidatedScalarMultivariateTaylorFunctionModelDP&>(this->_dwell_time_function.reference());
    new_set._dwell_time_function = ValidatedScalarMultivariateTaylorFunctionModelDP(new_domain, Ariadne::recondition(dwell_time.model(), discarded_parameters, number_of_error_parameters));

    new_set._auxiliary_mapping = this->_auxiliary_mapping;
    new_set._auxiliary_arguments = this->_auxiliary_arguments;
    new_set._auxiliary_parameter_count = this->_auxiliary_parameter_count;
    if (_constant_parameters) {
        const auto& parameter_function = dynamic_cast<const ValidatedVectorMultivariateTaylorFunctionModelDP&>(_constant_parameters->function.reference());
        Vector<ValidatedTaylorModelDP> new_parameters(constant_parameter_count(),
                                                      ValidatedTaylorModelDP(new_domain.size(), parameter_function.properties()));
        for (SizeType i = 0u; i != constant_parameter_count(); ++i)
            new_parameters[i] = Ariadne::recondition(parameter_function.models()[i], discarded_parameters, number_of_error_parameters);
        new_set._constant_parameters = _constant_parameters;
        new_set._constant_parameters->function = ValidatedVectorMultivariateTaylorFunctionModelDP(new_domain, new_parameters);
    }
    new_set._reduced_domain = new_reduced_domain;
    new_set._variable_kinds = List<EnclosureVariableKind>(new_domain.size(), EnclosureVariableKind::ERROR);
    for (SizeType j = 0; j < number_of_kept_parameters; ++j) {
        new_set._variable_kinds[j] = this->_variable_kinds[kept_parameters[j]];
    }

    (*this) = new_set;
}

Void Enclosure::restrict(const ExactBoxType& subdomain) {
    ARIADNE_PRECONDITION_MSG(subdomain.size() == this->number_of_parameters(), "set=" << *this << ", subdomain=" << subdomain);
    ARIADNE_PRECONDITION_MSG(Ariadne::subset(subdomain, this->domain()), "set.domain()=" << this->domain() << ", subdomain=" << subdomain);
    Enclosure& result(*this);
    result._domain = subdomain;
    result._reduced_domain = Ariadne::intersection(static_cast<const ExactBoxType&>(result._reduced_domain), subdomain);
    result._state_function = restriction(result._state_function, subdomain);
    result._time_function = restriction(result._time_function, subdomain);
    result._dwell_time_function = restriction(result._dwell_time_function, subdomain);
    if (result._constant_parameters)
        result._constant_parameters->function = restriction(result._constant_parameters->function, subdomain);
    ValidatedScalarMultivariateFunctionPatch new_constraint;
    for (List<ValidatedConstraintModel>::Iterator iter = result._constraints.begin();
         iter != result._constraints.end(); ++iter) {
        ValidatedScalarMultivariateFunctionPatch& constraint_function = iter->function();
        constraint_function = restriction(constraint_function, subdomain);
    }
    this->_expose_constant_parameter_errors();
    this->reduce();
}

Enclosure restriction(const Enclosure& encl, const ExactBoxType& subdomain) {
    Enclosure result(encl);
    result.restrict(subdomain);
    return result;
}

const ValidatedScalarMultivariateFunctionPatch Enclosure::get_function(SizeType i) const {
    if (i < this->state_dimension()) {
        return this->_state_function[i];
    } else if (i == this->state_dimension()) {
        return this->_time_function;
    } else {
        return this->auxiliary_function()[i - this->state_dimension() - 1u];
    }
}

Void Enclosure::draw(CanvasInterface& canvas, const Projection2d& projection) const {
    DrawerManager::instance().drawer().draw(canvas, projection, this->state_time_auxiliary_set());
}

template<class K, class V>
Map<K, V> filter(const Map<K, V>& m, const Set<K>& s) {
    Map<K, V> r;
    for (typename Set<K>::ConstIterator iter = s.begin(); iter != s.end(); ++iter) {
        r.insert(*m.find(*iter));
    }
    return r;
}

OutputStream& Enclosure::_write(OutputStream& os) const {
    // TODO: improve by using writers

    os << "Enclosure"
       << "(\n  domain=" << this->domain()
       << ",\n  range=" << this->bounding_box()
       << ",\n  reduced_domain=" << this->reduced_domain()
       << ",\n  is_empty=" << this->reduced_domain().is_empty()
       << ",\n  state_function=" << this->state_function()
       << ",\n  time_function=" << this->time_function()
       << ",\n  constraints=" << this->constraints()
       << ",\n  auxiliary_mapping = [" << this->auxiliary_mapping().argument_size() << "]" << this->auxiliary_mapping()
       << ",\n  parameters = " << canonical_variable_names(this->variable_kinds())
       << "\n)\n";
    return os;
}

Enclosure product(const Enclosure& set, const ExactIntervalType& ivl) {
    typedef List<ValidatedConstraintModel>::ConstIterator ConstIterator;

    ValidatedVectorMultivariateFunctionPatch new_function = combine(set.state_function(), set.configuration().function_factory().create_coordinate({ivl}, 0));

    Enclosure result(new_function.domain(), new_function, set.configuration());
    for (ConstIterator iter = set._constraints.begin(); iter != set._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(), embed(iter->function(), ivl), iter->upper_bound()));
    }
    result._time_function = embed(set._time_function, ivl);
    result._dwell_time_function = embed(set._dwell_time_function, ivl);
    if (set._constant_parameters) {
        result._constant_parameters = set._constant_parameters;
        result._constant_parameters->function = embed(set._constant_parameters->function, ivl);
        result._variable_kinds = set._variable_kinds;
        result._variable_kinds.append(EnclosureVariableKind::INITIAL);
    }

    return result;
}

Enclosure product(const Enclosure& set, const ExactBoxType& bx) {
    typedef List<ValidatedConstraintModel>::ConstIterator ConstIterator;

    ValidatedVectorMultivariateFunctionPatch new_function = combine(set.state_function(), set.configuration().function_factory().create_identity(bx));

    Enclosure result(new_function.domain(), new_function, set.configuration());
    for (ConstIterator iter = set._constraints.begin(); iter != set._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(), embed(iter->function(), bx), iter->upper_bound()));
    }
    result._time_function = embed(set._time_function, bx);
    result._dwell_time_function = embed(set._dwell_time_function, bx);
    if (set._constant_parameters) {
        result._constant_parameters = set._constant_parameters;
        result._constant_parameters->function = embed(set._constant_parameters->function, bx);
        result._variable_kinds = set._variable_kinds;
        result._variable_kinds.concatenate(List<EnclosureVariableKind>(bx.dimension(), EnclosureVariableKind::INITIAL));
    }

    return result;
}

Enclosure product(const Enclosure& set1, const Enclosure& set2) {
    if (set1.has_constant_parameter_binding() || set2.has_constant_parameter_binding())
        throw std::invalid_argument("Products of enclosures with physical parameter bindings require an explicit reconciliation policy");
    ARIADNE_PRECONDITION(same(set1.time_function().range(), set2.time_function().range()));
    ARIADNE_PRECONDITION(same(set1.dwell_time_function().range(), set2.dwell_time_function().range()));

    typedef List<ValidatedConstraintModel>::ConstIterator ConstIterator;

    ValidatedVectorMultivariateFunctionPatch new_state_function = combine(set1.state_function(), set2.state_function());

    Enclosure result(new_state_function.domain(), new_state_function, set1.configuration());
    for (ConstIterator iter = set1._constraints.begin(); iter != set1._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(), embed(iter->function(), set2.domain()), iter->upper_bound()));
    }
    for (ConstIterator iter = set2._constraints.begin(); iter != set2._constraints.end(); ++iter) {
        result._constraints.append(ValidatedConstraintModel(iter->lower_bound(), embed(set1.domain(), iter->function()), iter->upper_bound()));
    }
    result._time_function = embed(set1.time_function(), set2.time_function().domain());
    result._dwell_time_function = embed(set1.dwell_time_function(), set2.dwell_time_function().domain());

    return result;
}

Enclosure apply(const ValidatedVectorMultivariateFunction& function, const Enclosure& set) {
    Enclosure result(set);
    result.apply_map(function);
    return result;
}

Enclosure unchecked_apply(const ValidatedVectorMultivariateFunctionPatch& function, const Enclosure& set) {
    Enclosure result(set);
    const ValidatedVectorMultivariateFunctionPatch& state_function = result.state_function();
    const_cast<ValidatedVectorMultivariateFunctionPatch&>(state_function) = unchecked_compose(function, set.state_function());
    return result;
}

LabelledEnclosure::LabelledEnclosure(const LabelledExactBoxType& bx, const EnclosureConfiguration& config)
    : Enclosure(bx.euclidean_set(), config), _state_variables(bx.space().variable_names()) {
}

LabelledEnclosure::LabelledEnclosure(const ExactBoxType& bx, const RealSpace& state_space, const EnclosureConfiguration& config)
    : Enclosure(bx, config), _state_variables(state_space.variable_names()) {
}

LabelledEnclosure::LabelledEnclosure(const RealBox& bx, const RealSpace& state_space, const EnclosureConfiguration& config)
    : Enclosure(bx, config), _state_variables(state_space.variable_names()) {
}

LabelledEnclosure::LabelledEnclosure(const RealVariablesBox& bx, const RealSpace& state_space, const EnclosureConfiguration& config)
    : Enclosure(bx.euclidean_set(state_space), config), _state_variables(state_space.variable_names()) {
}

LabelledEnclosure::LabelledEnclosure(const BoundedConstraintSet& set, const RealSpace& state_space, const EnclosureConfiguration& config)
    : Enclosure(set, config), _state_variables(state_space.variable_names()) {
}

LabelledEnclosure::LabelledEnclosure(const Enclosure& set, const RealSpace& state_space)
    : Enclosure(set), _state_variables(state_space.variable_names()) {
}

LabelledEnclosure::LabelledEnclosure(const Enclosure& set, const RealSpace& state_space, const RealSpace& auxiliary_space)
    : Enclosure(set), _state_variables(state_space.variable_names()), _auxiliary_variables(auxiliary_space.variable_names()) {
}

LabelledEnclosure* LabelledEnclosure::clone() const {
    return new LabelledEnclosure(*this);
}

Void LabelledEnclosure::set_state_space(const RealSpace& space) {
    ARIADNE_PRECONDITION(this->state_function().result_size() == space.dimension());
    this->_state_variables = space.variable_names();
}

const RealSpace LabelledEnclosure::state_space() const {
    return RealSpace(this->_state_variables);
}

const RealVariable LabelledEnclosure::time_variable() const {
    return TimeVariable();
}

const RealSpace LabelledEnclosure::auxiliary_space() const {
    return RealSpace(this->_auxiliary_variables);
}

const RealSpace LabelledEnclosure::state_auxiliary_space() const {
    return join(this->state_space(), this->auxiliary_space());
}

const RealSpace LabelledEnclosure::state_time_auxiliary_space() const {
    auto state_time_space = (this->state_space().contains(this->time_variable()) ? this->state_space() : join(this->state_space(), this->time_variable()));
    return join(state_time_space, this->auxiliary_space());
}

const ValidatedVectorMultivariateFunctionPatch LabelledEnclosure::state_time_auxiliary_function() const {
    if (this->state_space().contains(this->time_variable())) {
        return this->state_auxiliary_function();
    }
    return this->Enclosure::state_time_auxiliary_function();
}

const RealSpace LabelledEnclosure::space() const {
    return this->state_time_auxiliary_space();
}

Void LabelledEnclosure::set_auxiliary(const RealSpace& spc, const EffectiveVectorMultivariateFunction& aux, EnclosureFunctionArguments arguments, SizeType parameter_count) {
    ARIADNE_PRECONDITION(spc.size() == aux.result_size());
    this->Enclosure::set_auxiliary_mapping(aux, arguments, parameter_count);
    this->_auxiliary_variables = spc.variable_names();
}

Pair<LabelledEnclosure, LabelledEnclosure> LabelledEnclosure::split() const {
    Pair<Enclosure, Enclosure> split_enclosures = this->Enclosure::split();
    return make_pair(LabelledEnclosure(split_enclosures.first, this->state_space(), this->auxiliary_space()),
                     LabelledEnclosure(split_enclosures.second, this->state_space(), this->auxiliary_space()));
}

LabelledEnclosure product(const LabelledEnclosure& set1, const ExactIntervalType& ivl2) {
    std::stringstream strstr;
    strstr << "x" << to_str(set1.dimension());
    return product(set1, LabelledExactIntervalType(RealVariable(strstr.str()), ivl2));
}

LabelledEnclosure product(const LabelledEnclosure& set1, const LabelledExactIntervalType& ivl2) {
    return LabelledEnclosure(product(set1.euclidean_set(), ivl2.interval()),
                             join(set1.state_space(), ivl2.variable()),
                             set1.auxiliary_space());
}

LabelledEnclosure product(const LabelledEnclosure& set1, const LabelledExactBoxType& bx2) {
    return LabelledEnclosure(product(set1.euclidean_set(), bx2.euclidean_set()),
                             join(set1.state_space(), bx2.space()),
                             set1.auxiliary_space());
}

LabelledEnclosure product(const LabelledEnclosure& set1, const LabelledEnclosure& set2) {
    return LabelledEnclosure(product(static_cast<const Enclosure&>(set1), static_cast<const Enclosure&>(set2)),
                             join(set1.state_space(), set2.state_space()), join(set1.auxiliary_space(), set2.auxiliary_space()));
}

Void LabelledEnclosure::apply_full_reach_step(ValidatedVectorMultivariateFunctionPatch phi) {
    this->Enclosure::apply_full_reach_step(phi);
}

Void LabelledEnclosure::apply_map(const ValidatedVectorMultivariateFunction& f) {
    this->Enclosure::apply_map(f);
}

Void LabelledEnclosure::apply_map(const ValidatedVectorMultivariateFunction& f, const RealSpace& spc) {
    ARIADNE_PRECONDITION(spc.size() == f.result_size());
    this->Enclosure::apply_map(f);
    this->_state_variables = spc.variable_names();
}

Void LabelledEnclosure::apply_map(const ValidatedVectorMultivariateFunction& f, const RealSpace& spc,
                                  const EffectiveVectorMultivariateFunction& aux_f, const RealSpace& aux_spc,
                                  EnclosureFunctionArguments source_arguments, EnclosureFunctionArguments auxiliary_arguments) {
    ARIADNE_PRECONDITION(aux_spc.size() == aux_f.result_size());
    this->Enclosure::apply_map(f, aux_f, source_arguments, auxiliary_arguments);
    this->_state_variables = spc.variable_names();
    this->_auxiliary_variables = aux_spc.variable_names();
}

Void LabelledEnclosure::draw(CanvasInterface& canvas, const Variables2d& axes) const {
    Projection2d proj = projection(this->state_time_auxiliary_space(), axes);
    ValidatedConstrainedImageSet set(this->domain(), this->state_time_auxiliary_function(), this->constraints());
    DrawerManager::instance().drawer().draw(canvas, proj, set);
}

OutputStream& LabelledEnclosure::_write(OutputStream& os) const {
    // TODO: improve by using writers

    os << "LabelledEnclosure"
       << "(\n  domain=" << this->domain()
       << ",\n  range=" << this->bounding_box()
       << ",\n  reduced_domain=" << this->reduced_domain()
       << ",\n  is_empty=" << this->reduced_domain().is_empty()
       << ",\n  state_function=" << this->state_function()
       << ",\n  time_function=" << this->time_function()
       << ",\n  constraints=" << this->constraints()
       << ",\n  auxiliary_mapping = [" << this->auxiliary_mapping().argument_size() << "]" << this->auxiliary_mapping()
       << ",\n  parameters = " << canonical_variable_names(this->variable_kinds())
       << ",\n  space =" << this->space()
       << ",\n  auxiliary_space =" << this->auxiliary_space()
       << "\n)\n";
    return os;
}

LabelledUpperBoxType LabelledEnclosure::bounding_box() const {
    return LabelledBox(this->state_space(), this->euclidean_set().bounding_box());
}

LabelledUpperBoxType
    ListSet<LabelledEnclosure>::bounding_box() const {
    ARIADNE_PRECONDITION(this->size() > 0u);
    auto iter = this->begin();
    auto box = iter->state_auxiliary_set().bounding_box();
    RealSpace space = join(iter->state_space(), iter->auxiliary_space());
    ++iter;
    for (; iter != this->end(); ++iter) {
        box = hull(box, iter->state_auxiliary_set().bounding_box());
    }
    return LabelledUpperBoxType(space, box);
}

const ListSet<LabelledSet<UpperBoxType>> ListSet<LabelledEnclosure>::bounding_boxes() const {
    ListSet<LabelledSet<UpperBoxType>> boxes;
    for (SizeType i = 0; i != this->size(); ++i) {
        boxes.adjoin(LabelledSet((*this)[i].state_space(), (*this)[i].euclidean_set().bounding_box()));
    }
    return boxes;
}

ListSet<LabelledEnclosure>* ListSet<LabelledEnclosure>::clone() const {
    return new ListSet<LabelledEnclosure>(*this);
}

Void ListSet<LabelledEnclosure>::draw(CanvasInterface& cnvs, const Variables2d& prj) const {
    draw_enclosures_recorded(cnvs, *this, [&](const LabelledEnclosure& enclosure, RecordingCanvas& recording) {
        enclosure.draw(recording, prj);
    });
}

} // namespace Ariadne
