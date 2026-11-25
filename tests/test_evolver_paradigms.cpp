/***************************************************************************
 *            test_evolver_paradigms.cpp
 *
 *  Copyright  2007-20  Alberto Casagrande, Pieter Collins
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

/** @file tests/test_evolver_paradigms.cpp
 *  @brief Verify distinct numerical paradigm contracts and polymorphic evolution.
 */

#include <memory>
#include <type_traits>

#include "logging/test_runner.hpp"

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

using ApproximateContinuous = ContinuousTimeEvolver<ApproximateTag>;
using ValidatedContinuous = ContinuousTimeEvolver<ValidatedTag>;
using ApproximateHybrid = HybridContinuousTimeEvolver<ApproximateTag>;
using ValidatedHybrid = HybridContinuousTimeEvolver<ValidatedTag>;

template<class E>
concept HasSampledOrbit = requires(const E& e, const typename E::InitialType& initial,
                                   const typename E::TerminationType& termination) {
    e.orbit(initial, termination);
};
template<class E>
concept HasSetSemantics = requires(const E& e, const typename E::InitialType& initial,
                                   const typename E::TerminationType& termination) {
    e.orbit(initial, termination, Semantics::UPPER);
};
template<class P, class S>
concept HasEvolutionTraits = requires { typename EvolutionTraits<P, S>::EvolverType; };
template<class T>
concept CompleteType = requires { sizeof(T); };

static_assert(std::is_same_v<ApproximateContinuous::Paradigm, ApproximateTag>);
static_assert(std::is_same_v<ValidatedContinuous::Paradigm, ValidatedTag>);
static_assert(std::is_same_v<ApproximateHybrid::Paradigm, ApproximateTag>);
static_assert(std::is_same_v<ValidatedHybrid::Paradigm, ValidatedTag>);
static_assert(std::is_same_v<DiscreteTimeEvolver<ValidatedTag>::Paradigm, ValidatedTag>);
static_assert(!std::is_same_v<ApproximateContinuous::InitialType, ValidatedContinuous::InitialType>);
static_assert(!std::is_same_v<ApproximateContinuous::ResultType, ValidatedContinuous::ResultType>);
static_assert(!std::is_same_v<ApproximateHybrid::ResultType, ValidatedHybrid::ResultType>);
static_assert(HasSampledOrbit<ApproximateContinuous> && !HasSetSemantics<ApproximateContinuous>);
static_assert(HasSampledOrbit<ApproximateHybrid> && !HasSetSemantics<ApproximateHybrid>);
static_assert(!HasSampledOrbit<ValidatedContinuous> && HasSetSemantics<ValidatedContinuous>);
static_assert(!HasSampledOrbit<ValidatedHybrid> && HasSetSemantics<ValidatedHybrid>);
static_assert(!HasSampledOrbit<DiscreteTimeEvolver<ValidatedTag>>);
static_assert(HasEvolutionTraits<ApproximateTag, OrdinaryDifferentialAlgebraicEquation>);
static_assert(HasEvolutionTraits<ValidatedTag, OrdinaryDifferentialAlgebraicInclusion>);
static_assert(HasEvolutionTraits<ApproximateTag, OrdinaryDifferentialAlgebraicInclusion>);
static_assert(std::is_constructible_v<ApproximateContinuous, const OrdinaryDifferentialAlgebraicInclusion&, const InputSamplingPolicy&>);
static_assert(HasEvolutionTraits<ValidatedTag, DifferenceAlgebraicEquation>);
static_assert(std::is_same_v<DiscreteTimeEvolver<ValidatedTag>::SystemType,
                             DifferenceAlgebraicEquation>);
static_assert(HasEvolutionTraits<ApproximateTag, DifferenceAlgebraicEquation>);
static_assert(HasEvolutionTraits<ApproximateTag, DifferenceAlgebraicInclusion>);
static_assert(std::is_constructible_v<DiscreteTimeEvolver<ApproximateTag>,
    const DifferenceAlgebraicInclusion&, const InputSamplingPolicy&>);
static_assert(!std::is_constructible_v<DiscreteTimeEvolver<ApproximateTag>,
    const DifferenceAlgebraicInclusion&>);
static_assert(CompleteType<DiscreteTimeEvolver<ApproximateTag>>);
static_assert(HasSampledOrbit<DiscreteTimeEvolver<ApproximateTag>>);
static_assert(!HasSetSemantics<DiscreteTimeEvolver<ApproximateTag>>);
static_assert(!CompleteType<ContinuousTimeEvolver<ExactTag>>);
static_assert(!std::is_constructible_v<ApproximateContinuous, const OrdinaryDifferentialAlgebraicInclusion&>);
static_assert(!std::is_constructible_v<ApproximateContinuous, const DifferenceAlgebraicEquation&>);

class TestEvolverParadigms : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_continuous_polymorphic_contracts() {
        RealVariable x("x");
        const TimeVariable t;
        const OrdinaryDifferentialAlgebraicEquation equation({dot(x) = t});
        ApproximateContinuous approximate(equation);
        approximate.configuration().set_step_size(0.03125);
        std::unique_ptr<ApproximateContinuous::Interface> sampled(approximate.clone());
        ARIADNE_TESTING_PROVE(sampled->system().has_explicit_time());
        const ApproximateContinuous::InitialType points{
            ApproximateContinuous::ApproximatePointType(RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp))};
        auto path = sampled->orbit(points, Real(0.125_x));
        ARIADNE_TESTING_EXPECT_WITHIN(path.curve(0u).end()->second[0u].raw().get_d(), 0.0078125, 1e-14);

        TaylorPicardIntegrator integrator(step_maximum_error = 1e-7);
        ValidatedContinuous validated(equation, integrator);
        validated.configuration().set_maximum_step_size(0.03125);
        std::unique_ptr<ValidatedContinuous::Interface> certified(validated.clone());
        const auto result = certified->orbit(validated.enclosure(ExactBoxType({
                                                 {0, 0}
        })),
                                             Real(0.125_x), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(result.final().size() > 0u);
        for (const auto& final: result.final()) {
            ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.0078125_x, 0.0078125_x),
                                         cast_exact_interval(final.state_function()[0u].range())));
        }
    }

    void test_hybrid_approximate_interface() {
        RealVariable x("x");
        const DiscreteLocation q(1);
        ContinuousTimeHybridAutomaton system;
        system.new_mode(q, OrdinaryDifferentialAlgebraicEquation({dot(x) = 1}));
        ApproximateHybrid evolver(system);
        evolver.configuration().set_step_size(0.03125);
        std::unique_ptr<HybridContinuousTimeEvolverInterface<ApproximateTag>> sampled(evolver.clone());
        ARIADNE_TESTING_EXPECT_EQUAL(sampled->system().continuous_state_space(q), RealSpace({x}));
        const ApproximateHybrid::InitialType initial(q, RealSpace({x}), FloatDPApproximationPoint({0.0_x}, dp));
        const auto result = sampled->orbit(initial, HybridTime(0.125_x, 1));
        ARIADNE_TESTING_EXPECT_WITHIN(result.curves()[0u].euclidean_set().end()->second[0u].raw().get_d(), 0.125, 1e-14);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_evolver_paradigms");
    ARIADNE_TESTING_ADD_TEST(TestEvolverParadigms, test_continuous_polymorphic_contracts);
    ARIADNE_TESTING_ADD_TEST(TestEvolverParadigms, test_hybrid_approximate_interface);
    return Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
