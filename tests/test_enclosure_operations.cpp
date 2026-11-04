#include "evolution/enclosure.hpp"
#include "function/taylor_function.hpp"
#include "labelled/figure.hpp"
#include "paving/paver.hpp"
#include "set/list_set.hpp"
#include "macro/testing.hpp"
#include <sstream>
using namespace Ariadne;

class TestEnclosureOperations {
    ExactBoxType domain{{-1,1}};
    EnclosureConfiguration config{TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8))};
    EffectiveVectorMultivariateFunction identity=EffectiveVectorMultivariateFunction::identity(1u);
    EffectiveScalarMultivariateFunction coordinate=identity[0];
    EffectiveScalarMultivariateFunction zero=EffectiveScalarMultivariateFunction::zero(1u);
    EffectiveScalarMultivariateFunction quarter=EffectiveScalarMultivariateFunction::constant(1u,Real(1/4_q));
    auto const& factory() const { return config.function_factory(); }
    Enclosure initial() const { return Enclosure(domain,config); }
    void expect_value(ValidatedScalarMultivariateFunctionPatch const& f,double expected) const {
        auto value=f(Vector<FloatDPApproximation>({0.0_x},dp));
        ARIADNE_TEST_EQUALS(value.raw().get_d(),expected);
    }
public:
    void test() {
        ARIADNE_TEST_CALL(test_constructors_and_configuration());
        ARIADNE_TEST_CALL(test_evolve_and_reach_steps());
        ARIADNE_TEST_CALL(test_map_and_constraint_contracts());
        ARIADNE_TEST_CALL(test_labelled_adapters_and_parameter_names());
    }
    void test_constructors_and_configuration() {
        EnclosureConfiguration copied(static_cast<ValidatedFunctionPatchFactory::Interface const&>(factory()),2u);
        ARIADNE_TEST_EQUALS(copied.reconditioning_num_blocks(),2u);
        copied.set_paver(Paver(new SubdivisionPaver()));
        std::ostringstream description; description<<copied.paver();
        ARIADNE_TEST_ASSERT(description.str().find("SubdivisionPaver")!=std::string::npos);
        Enclosure plain(domain,identity,config);
        Enclosure constrained(domain,identity,List<ValidatedConstraint>{ValidatedConstraint(-1,coordinate,1)},config);
        ARIADNE_TEST_EQUALS(constrained.number_of_constraints(),1u);
        ARIADNE_TEST_EQUALS(plain.state_dimension(),1u);
        ARIADNE_TEST_FAIL(Enclosure(domain,EffectiveVectorMultivariateFunction::identity(2u),zero,List<ValidatedConstraint>{},config));
        ARIADNE_TEST_FAIL(Enclosure(domain,identity,EffectiveScalarMultivariateFunction::zero(2u),List<ValidatedConstraint>{},config));
        ARIADNE_TEST_FAIL(Enclosure(domain,identity,List<ValidatedConstraint>{ValidatedConstraint(-1,EffectiveScalarMultivariateFunction::zero(2u),1)},config));
    }
    void test_evolve_and_reach_steps() {
        auto coordinates=EffectiveVectorMultivariateFunction::identity(2u);
        EffectiveVectorMultivariateFunction flow={coordinates[0]+coordinates[1]};
        auto flow_domain=product(domain,ExactIntervalType(0.0_x,0.5_x));
        auto phi=factory().create(flow_domain,flow);
        auto a=initial(); a.apply_fixed_evolve_step(flow,Dyadic(1/4_q));
        expect_value(a.state_function()[0],0.25); expect_value(a.time_function(),0.25);
        auto b=initial(); b.apply_space_evolve_step(flow,quarter);
        expect_value(b.state_function()[0],0.25); expect_value(b.time_function(),0.25);
        auto c=initial(); c.apply_spacetime_evolve_step(flow,EffectiveScalarMultivariateFunction::constant(2u,Real(1/4_q)));
        expect_value(c.state_function()[0],0.25); expect_value(c.time_function(),0.25);
        auto d=initial(); d.apply_parameter_evolve_step(flow,quarter);
        expect_value(d.state_function()[0],0.25); expect_value(d.time_function(),0.25);
        d.apply_finishing_parameter_evolve_step(flow,EffectiveScalarMultivariateFunction::constant(1u,Real(1/2_q)));
        expect_value(d.state_function()[0],0.5); expect_value(d.time_function(),0.5);
        auto direct=initial();
        direct.apply_parameterised_fixed_evolve_step(factory().create(domain,EffectiveVectorMultivariateFunction{coordinate+1}),Dyadic(1/4_q));
        expect_value(direct.state_function()[0],1.0); expect_value(direct.time_function(),0.25);
        auto reach=initial(); reach.apply_parameterised_full_reach_step(phi);
        ARIADNE_TEST_EQUALS(reach.number_of_parameters(),2u);
        ARIADNE_TEST_EQUALS(reach.variable_kinds().back(),EnclosureVariableKind::TEMPORAL);
        ARIADNE_TEST_ASSERT(definitely(reach.inside(ExactBoxType({{-1.01_x,1.51_x}}))));
        ARIADNE_TEST_FAIL(initial().apply_parameterised_fixed_evolve_step(factory().create_identity(ExactBoxType({{-1,1},{-1,1}})),Dyadic(1/4_q)));
        ARIADNE_TEST_FAIL(initial().apply_parameterised_fixed_evolve_step(factory().create_identity(ExactBoxType({{0,1}})),Dyadic(1/4_q)));
        ARIADNE_TEST_FAIL(initial().apply_parameterised_fixed_evolve_step(factory().create_zeros(2u,domain),Dyadic(1/4_q)));
        ARIADNE_TEST_FAIL(initial().apply_parameterised_full_reach_step(factory().create_zeros(2u,flow_domain)));
        ARIADNE_TEST_FAIL(initial().apply_parameterised_full_reach_step(factory().create_identity(domain)));
        auto incompatible_domain=ExactBoxType({{0,1},{0,1}});
        ARIADNE_TEST_FAIL(initial().apply_parameterised_full_reach_step(factory().create_zeros(1u,incompatible_domain)));
        auto wrong_argument=identity;
        auto wrong_result=EffectiveVectorMultivariateFunction::identity(2u);
        ARIADNE_TEST_FAIL(initial().apply_fixed_evolve_step(wrong_argument,Dyadic(1/4_q)));
        ARIADNE_TEST_FAIL(initial().apply_fixed_evolve_step(wrong_result,Dyadic(1/4_q)));
        ARIADNE_TEST_FAIL(initial().apply_space_evolve_step(wrong_argument,quarter));
        ARIADNE_TEST_FAIL(initial().apply_space_evolve_step(wrong_result,quarter));
        ARIADNE_TEST_FAIL(initial().apply_space_evolve_step(flow,EffectiveScalarMultivariateFunction::zero(2u)));
        ARIADNE_TEST_FAIL(initial().apply_spacetime_evolve_step(wrong_argument,zero));
        ARIADNE_TEST_FAIL(initial().apply_spacetime_evolve_step(wrong_result,zero));
        ARIADNE_TEST_FAIL(initial().apply_spacetime_evolve_step(flow,zero));
        ARIADNE_TEST_FAIL(initial().apply_parameter_evolve_step(wrong_argument,zero));
        ARIADNE_TEST_FAIL(initial().apply_parameter_evolve_step(wrong_result,zero));
        ARIADNE_TEST_FAIL(initial().apply_parameter_evolve_step(flow,EffectiveScalarMultivariateFunction::zero(2u)));
        ARIADNE_TEST_FAIL(initial().apply_finishing_parameter_evolve_step(wrong_argument,zero));
        ARIADNE_TEST_FAIL(initial().apply_finishing_parameter_evolve_step(wrong_result,zero));
        ARIADNE_TEST_FAIL(initial().apply_finishing_parameter_evolve_step(flow,EffectiveScalarMultivariateFunction::zero(2u)));
    }
    void test_map_and_constraint_contracts() {
        EffectiveVectorMultivariateFunction map={2*coordinate};
        auto mapped=apply(map,initial()); expect_value(mapped.state_function()[0],0.0);
        auto unchecked=unchecked_apply(factory().create(domain,map),initial());
        ARIADNE_TEST_ASSERT(definitely(unchecked.inside(ExactBoxType({{-2.01_x,2.01_x}}))));
        auto a=initial(); a.set_auxiliary_mapping(map);
        ARIADNE_TEST_EQUALS(a.get_function(0u).argument_size(),1u);
        expect_value(a.get_function(1u),0.0); expect_value(a.get_function(2u),0.0);
        ARIADNE_TEST_FAIL(initial().apply_map(EffectiveVectorMultivariateFunction::identity(2u)));
        ARIADNE_TEST_FAIL(initial().apply_map(EffectiveVectorMultivariateFunction::zeros(2u,1u)));
        ARIADNE_TEST_FAIL(initial().apply_map(identity,EffectiveVectorMultivariateFunction::zeros(0u,2u)));
        auto wrong=EffectiveScalarMultivariateFunction::zero(2u);
        ARIADNE_TEST_FAIL(initial().new_state_time_bound(wrong));
        ARIADNE_TEST_FAIL(initial().new_state_constraint(ValidatedConstraint(-1,wrong,1)));
        ARIADNE_TEST_FAIL(initial().new_state_time_constraint(ValidatedConstraint(-1,zero,1)));
        ARIADNE_TEST_FAIL(initial().new_parameter_constraint(ValidatedConstraint(-1,wrong,1)));
        ARIADNE_TEST_FAIL(initial().new_positive_state_constraint(wrong));
        ARIADNE_TEST_FAIL(initial().new_negative_state_constraint(wrong));
        ARIADNE_TEST_FAIL(initial().new_zero_state_constraint(wrong));
        ARIADNE_TEST_FAIL(initial().new_negative_parameter_constraint(wrong));
        ARIADNE_TEST_FAIL(initial().new_zero_parameter_constraint(wrong));
        ARIADNE_TEST_FAIL(initial().restrict(ExactBoxType({{0,1},{0,1}})));
        ARIADNE_TEST_FAIL(initial().restrict(ExactBoxType({{-2,2}})));
        ARIADNE_TEST_FAIL(initial().separated(ExactBoxType({{0,1},{0,1}})));
    }
    void test_labelled_adapters_and_parameter_names() {
        RealVariable x("x"),z("z");
        LabelledEnclosure labelled(initial(),RealSpace({x}));
        labelled.apply_map(identity,RealSpace({z}));
        ARIADNE_TEST_EQUALS(labelled.state_space(),RealSpace({z}));
        auto joined=product(labelled,ExactIntervalType(0.0_x,1.0_x));
        ARIADNE_TEST_EQUALS(joined.state_space().dimension(),2u);
        ListSet<LabelledEnclosure> sets; sets.adjoin(labelled); sets.adjoin(labelled);
        ARIADNE_TEST_EQUALS(sets.bounding_boxes().size(),2u);
        ARIADNE_TEST_EQUALS(ListSet<LabelledEnclosure>().bounding_boxes().size(),0u);
        List<EnclosureVariableKind> kinds={EnclosureVariableKind::INITIAL,EnclosureVariableKind::INITIAL,
            EnclosureVariableKind::TEMPORAL,EnclosureVariableKind::PARAMETER,EnclosureVariableKind::INPUT,
            EnclosureVariableKind::NOISE,EnclosureVariableKind::ERROR,EnclosureVariableKind::UNKNOWN};
        auto names=canonical_variable_names(kinds);
        ARIADNE_TEST_EQUALS(names,List<Identifier>({"x0","x1","t0","a0","u0","v0","e0","s0"}));
    }
};
int main() { TestEnclosureOperations().test(); return ARIADNE_TEST_FAILURES; }
