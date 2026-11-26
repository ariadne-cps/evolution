/***************************************************************************
 *            test_continuous_inclusions.cpp
 *
 *  Copyright  2026  Luca Geretti
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

/** @file tests/test_continuous_inclusions.cpp
 *  @brief Verify continuous inclusions behavior.
 */

#include <memory>
#include <sstream>

#include "evolution/continuous_time_evolver.hpp"
#include "function/taylor_model.hpp"
#include "integration/integrator.hpp"
#include "integration/validated_integration_step.hpp"
#include "logging/test_runner.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestContinuousInclusions : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_explicit_reach_and_endpoint_choices() {
        RealVariable x("x"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = u}, RealVariablesBox({-1 <= u <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        for (const auto& reach: List<InputApproximation>{ZeroApproximation(), ConstantApproximation()}) {
            for (const auto& evolve: List<InputApproximation>{ZeroApproximation(), ConstantApproximation(),
                                                              AffineApproximation(), SinusoidalApproximation(), PiecewiseApproximation()}) {
                ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator, reach, evolve);
                evolver.configuration().set_maximum_step_size(0.03125);
                const auto result = evolver.orbit(RealVariablesBox({-0.125_dec <= x <= 0.125_dec}),
                                                  Real(0.03125_x), Semantics::UPPER);
                ARIADNE_TESTING_PROVE(result.reach().size() > 0u);
                const auto final = cast_exact_box(result.final().bounding_box().euclidean_set());
                ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                                 {-0.15625_x, 0.15625_x}
                }),
                                             final));
                ARIADNE_TESTING_PROVE(subset(final, ExactBoxType({
                                                        {-0.1572265625_x, 0.1572265625_x}
                })));
            }
        }
        for (const auto& invalid: List<InputApproximation>{AffineApproximation(), SinusoidalApproximation(), PiecewiseApproximation()}) {
            ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                          ContinuousTimeEvolver<ValidatedTag>(system, integrator, invalid, ConstantApproximation()));
        }
    }

    Void test_reconditioning_has_no_winner_statistics() {
        RealVariable x("x"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = u}, RealVariablesBox({-1 <= u <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator, ZeroApproximation(), AffineApproximation());
        evolver.configuration().set_maximum_step_size(0.03125);
        const auto result = evolver.orbit(RealVariablesBox({-0.125_dec <= x <= 0.125_dec}), Real(0.25_x), Semantics::UPPER);
        const auto final = cast_exact_box(result.final().bounding_box().euclidean_set());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.375_x, 0.375_x}
        }),
                                     final));
        ARIADNE_TESTING_PROVE(subset(final, ExactBoxType({
                                                {-0.376953125_x, 0.376953125_x}
        })));
        ARIADNE_TESTING_PROVE(result.final()[0u].number_of_parameters() < 10u);
    }

    Void test_affine_dynamics_with_uncertain_inputs() {
        RealVariable x("x"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = u - x}, RealVariablesBox({0.5_dec <= u <= 1.5_dec}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        for (const auto& method: List<InputApproximation>{ConstantApproximation(), AffineApproximation()}) {
            ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator, ConstantApproximation(), method);
            evolver.configuration().set_maximum_step_size(0.03125);
            const auto result = evolver.orbit(RealVariablesBox({-0.03125_dec <= x <= 0.03125_dec}), Real(0.125_x), Semantics::UPPER);
            const FloatDPBounds decay = exp(FloatDPBounds(-0.125_x, dp));
            const auto lo = -FloatDPBounds(0.03125_x, dp) * decay + FloatDPBounds(0.5_x, dp) * (FloatDPBounds(1, dp) - decay);
            const auto hi = FloatDPBounds(0.03125_x, dp) * decay + FloatDPBounds(1.5_x, dp) * (FloatDPBounds(1, dp) - decay);
            const UpperBoxType analytical({UpperIntervalType(lo.lower(), hi.upper())});
            const auto box = cast_exact_box(result.final().bounding_box().euclidean_set());
            ARIADNE_TESTING_PROVE(subset(cast_exact_box(analytical), box));
            ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                                  {0.025390625_x, 0.2099609375_x}
            })));
        }
    }

    Void test_switching_inputs_are_independent_between_steps() {
        RealVariable x("x"), y("y"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = u, dot(y) = x}, RealVariablesBox({-1 <= u <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator, ConstantApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.0625);
        evolver.configuration().set_enable_reconditioning(false);
        const auto result = evolver.orbit(RealVariablesBox({x == 0, y == 0}), Real(0.125_x), Semantics::UPPER);
        const auto& final = result.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(final.number_of_parameters(), 4u);
        // u=+1 in the first half and -1 in the second gives x(H)=0,
        // y(H)=H^2/4. A single frozen input cannot produce this point.
        const Vector<FloatDPBounds> switching({FloatDPBounds(0, dp), FloatDPBounds(0, dp),
                                               FloatDPBounds(1, dp), FloatDPBounds(-1, dp)});
        const auto values = final.state_function()(switching);
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0, dp), values[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.00390625_x, dp), values[1u]));
    }

    Void test_uniform_initial_error_is_preserved_without_reconditioning() {
        RealVariable x("x"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = u}, RealVariablesBox({-1 <= u <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_enable_reconditioning(false);
        evolver.configuration().set_maximum_step_size(0.03125);
        const ExactBoxType domain({
            {0.0_x, 0.0_x}
        });
        auto model = ValidatedVectorMultivariateTaylorFunctionModelDP::identity(domain, ThresholdSweeper<FloatDP>(dp, 1e-10));
        model[0u].set_error(FloatDPError(0x1p-10_x, dp));
        const LabelledEnclosure initial(Enclosure(domain, model, EnclosureConfiguration(integrator.function_factory())), RealSpace({x}));
        const auto result = evolver.orbit(initial, Real(0.0625_x), Semantics::UPPER);
        const auto box = cast_exact_box(result.final().bounding_box().euclidean_set());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.0634765625_x, 0.0634765625_x}
        }),
                                     box));
        ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                              {-0.064453125_x, 0.064453125_x}
        })));
    }

    Void test_timed_constant_parameters_survive_input_choices() {
        RealVariable x("x"), t("t"), p("p"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = t + p + u},
                                                            RealVariablesBox({-0.125_dec <= u <= 0.125_dec}), RealVariablesBox({0.125_dec <= p <= 0.25_dec}), t);
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator, ZeroApproximation(), ConstantApproximation());
        // Input-approximation errors are certified separately from the ODE
        // solver's Taylor error. Resolve them with a smaller maximum step;
        // retain the same analytical endpoint and tight outer tolerance.
        evolver.configuration().set_maximum_step_size(0.0078125);
        evolver.configuration().set_enable_reconditioning(false);
        const auto result = evolver.orbit(RealVariablesBox({x == 0}), Real(0.09375_x), Semantics::UPPER);
        const auto& final = result.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(final.state_dimension(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(final.number_of_parameters(), 14u);
        const auto box = cast_exact_box(final.euclidean_set().bounding_box());
        ARIADNE_TESTING_NOTE(box);
        ARIADNE_TESTING_NOTE(final.state_function());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {0.00439453125_x, 0.03955078125_x}
        }),
                                     box));
        ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                              {0.00390625_x, 0.0400390625_x}
        })));
        Vector<FloatDPBounds> low(14u, FloatDPBounds(0, dp));
        Vector<FloatDPBounds> high = low;
        low[1u] = FloatDPBounds(0.125_x, dp);
        high[1u] = FloatDPBounds(0.25_x, dp);
        const auto delta = final.state_function()(high) - final.state_function()(low);
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.01171875_x, dp), delta[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.09375_x, dp), final.time_function()(low)));
    }

    Void test_autonomous_parameters_keep_the_declared_state() {
        RealVariable x("x"), p("p"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = p + u},
                                                            RealVariablesBox({-1 <= u <= 1}), RealVariablesBox({-0.5_dec <= p <= 0.5_dec}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.0078125);
        evolver.configuration().set_enable_reconditioning(false);
        const auto result = evolver.orbit(RealVariablesBox({x == 0}), Real(0.09375_x), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(!evolver.system().has_explicit_time());
        ARIADNE_TESTING_EXPECT_EQUAL(result.final()[0u].state_dimension(), 1u);
        const auto box = cast_exact_box(result.final().bounding_box().euclidean_set());
        ARIADNE_TESTING_NOTE(box);
        ARIADNE_TESTING_NOTE(result.final()[0u].state_function());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.140625_x, 0.140625_x}
        }),
                                     box));
        ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                              {-0.1416015625_x, 0.1416015625_x}
        })));
    }

    Void test_first_order_box_paths_and_zero_free_inputs() {
        RealVariable x("x"), t("t"), p("p"), u("u");
        const FirstOrderMultifunctionIntegrator integrator;
        const OrdinaryDifferentialAlgebraicInclusion system({dot(x) = t + p + u},
                                                            RealVariablesBox({-0.125_dec <= u <= 0.125_dec}), RealVariablesBox({0.125_dec <= p <= 0.25_dec}), t);
        ContinuousTimeEvolver<ValidatedTag> evolver(system, integrator);
        evolver.configuration().set_maximum_step_size(0.0625);
        const auto result = evolver.orbit(RealVariablesBox({x == 0}), Real(0.09375_x), Semantics::UPPER);
        const auto box = cast_exact_box(result.final().bounding_box().euclidean_set());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {0.00439453125_x, 0.03955078125_x}
        }),
                                     box));
        ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                              {0.0_x, 0.044921875_x}
        })));
        const OrdinaryDifferentialAlgebraicInclusion deterministic({dot(x) = 1}, RealVariablesBox{});
        ContinuousTimeEvolver<ValidatedTag> deterministic_evolver(deterministic, integrator);
        const auto deterministic_result = deterministic_evolver.orbit(RealVariablesBox({x == 0}), Real(0.0625_x), Semantics::UPPER);
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {0.0625_x, 0.0625_x}
        }),
                                     cast_exact_box(deterministic_result.final().bounding_box().euclidean_set())));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      evolver.orbit(RealVariablesBox({x == 0}), Real(0u), Semantics::LOWER));
    }


    Void test_structured_lower_inclusion_preserves_switching_inputs() {
        const RealVariable x("x"), y("y"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system(
            {dot(x) = u, dot(y) = x}, RealVariablesBox({-1 <= u <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(
            system, integrator, ConstantApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_step_size(0.0625);
        evolver.configuration().set_enable_reconditioning(false);

        const auto initial = RealVariablesBox({x == 0, y == 0});
        const auto lower = evolver.orbit(initial, Real(0.125_x), Semantics::LOWER);
        const auto upper = evolver.orbit(initial, Real(0.125_x), Semantics::UPPER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        ARIADNE_TESTING_REQUIRE(upper.final().size() == 1u);
        ARIADNE_TESTING_PROVE(lower.reach().size() > 0u);
        const auto& lower_final = lower.final()[0u];
        const auto& upper_final = upper.final()[0u];
        ARIADNE_TESTING_EXPECT_EQUAL(lower_final.state_dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(lower_final.number_of_parameters(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(
            lower_final.number_of_parameters(), upper_final.number_of_parameters());
        ARIADNE_TESTING_EXPECT_EQUAL(lower_final.number_of_constraints(), 0u);

        // Two independent full input choices, without freezing u:
        // u[0]=+1, u[1]=-1 gives x[2]=0 and y[2]=1/256.
        // The enclosure must contain the same actual trajectory for
        // LOWER, not merely a numerical input-midpoint witness.
        const Vector<FloatDPBounds> switching({
            FloatDPBounds(0, dp), FloatDPBounds(0, dp),
            FloatDPBounds(1, dp), FloatDPBounds(-1, dp)});
        const auto lower_values = lower_final.state_function()(switching);
        const auto upper_values = upper_final.state_function()(switching);
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0, dp), lower_values[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.00390625_x, dp), lower_values[1u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0, dp), upper_values[0u]));
        ARIADNE_TESTING_PROVE(refines(FloatDPBounds(0.00390625_x, dp), upper_values[1u]));
    }

    Void test_lower_preserves_unsplit_initial_enclosure() {
        const RealVariable x("x"), u("u");
        const OrdinaryDifferentialAlgebraicInclusion system(
            {dot(x) = u}, RealVariablesBox({-1 <= u <= 1}));
        TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> evolver(
            system, integrator, ZeroApproximation(), ConstantApproximation());
        evolver.configuration().set_maximum_enclosure_radius(0.01);
        evolver.configuration().set_enable_subdivisions(true);
        const auto source = RealVariablesBox({-0.25_dec <= x <= 0.25_dec});

        // At zero horizon there is no numerical step. A large
        // numerical LOWER source must remain a single enclosure,
        // not be partitioned into purported independent witnesses.
        const auto lower = evolver.orbit(source, Real(0), Semantics::LOWER);
        ARIADNE_TESTING_REQUIRE(lower.final().size() == 1u);
        const auto original = cast_exact_box(
            lower.final()[0u].bounding_box().euclidean_set());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({{-0.25_x, 0.25_x}}), original));
        ARIADNE_TESTING_EXPECT_EQUAL(lower.final()[0u].number_of_parameters(), 1u);

        const FirstOrderMultifunctionIntegrator box_integrator;
        ContinuousTimeEvolver<ValidatedTag> boxes(system, box_integrator);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            boxes.orbit(RealVariablesBox({x == 0}), Real(0.125_x), Semantics::LOWER));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
            boxes.orbit(RealVariablesBox({x == 0}), Real(0), Semantics::LOWER));
    }

    Void test_pointwise_interval_field_and_domain_rejection() {
        RealVariable x("x");
        const ExactBoxType domain({
            {-1.0_x, 1.0_x}
        });
        const auto uncertain = ValidatedIntervalTaylorFunctionModel<FloatDP>::constant(
            domain, FloatDPUpperInterval(-1, 1, dp), ThresholdSweeper<FloatDP>(dp, 1e-10));
        const OrdinaryDifferentialAlgebraicInclusion system(RealSpace({x}),
                                                            ValidatedVectorIntervalTaylorFunctionModel<FloatDP>({uncertain}));
        const FirstOrderMultifunctionIntegrator boxes;
        ContinuousTimeEvolver<ValidatedTag> evolver(system, boxes);
        evolver.configuration().set_maximum_step_size(0.0625);
        const auto result = evolver.orbit(RealVariablesBox({x == 0}), Real(0.09375_x), Semantics::UPPER);
        const auto box = cast_exact_box(result.final().bounding_box().euclidean_set());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.09375_x, 0.09375_x}
        }),
                                     box));
        ARIADNE_TESTING_PROVE(subset(box, ExactBoxType({
                                              {-0.0947265625_x, 0.0947265625_x}
        })));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      evolver.orbit(RealVariablesBox({x == 2}), Real(0.03125_x), Semantics::UPPER));
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        ContinuousTimeEvolver<ValidatedTag> structured(system, base, ZeroApproximation(), ConstantApproximation());
        structured.configuration().set_maximum_step_size(0.0625);
        const auto structured_result = structured.orbit(RealVariablesBox({x == 0}), Real(0.09375_x), Semantics::UPPER);
        const auto structured_box = cast_exact_box(structured_result.final().bounding_box().euclidean_set());
        ARIADNE_TESTING_PROVE(subset(ExactBoxType({
                                         {-0.09375_x, 0.09375_x}
        }),
                                     structured_box));
        ARIADNE_TESTING_PROVE(subset(structured_box, ExactBoxType({
                                                         {-0.0947265625_x, 0.0947265625_x}
        })));
    }

    Void test_owned_components_and_configuration() {
        RealVariable x("x"), u("u");
        const auto make = [&]() {
            OrdinaryDifferentialAlgebraicInclusion system({dot(x) = u}, RealVariablesBox({-1 <= u <= 1}));
            TaylorPicardIntegrator integrator(step_maximum_error = 1e-6);
            return ContinuousTimeEvolver<ValidatedTag>(system, integrator, ZeroApproximation(), ConstantApproximation());
        };
        const auto evolver = make();
        std::unique_ptr<ContinuousTimeEvolver<ValidatedTag>> clone(evolver.clone());
        clone->configuration().set_maximum_step_size(0.03125);
        const auto result = clone->orbit(RealVariablesBox({x == 0}), Real(0.0625_x), Semantics::UPPER);
        ARIADNE_TESTING_EXPECT_EQUAL(result.reach().size(), 2u);
        std::ostringstream text;
        text << clone->configuration();
        ARIADNE_TESTING_PROVE(text.str().find("maximum_step_size=") != std::string::npos);
        auto source = clone->enclosure(ExactBoxType({
            {0.0_x, 0.0_x}
        }));
        ExactBoxType protected_coordinates(6u);
        for (SizeType i = 0u; i != protected_coordinates.dimension(); ++i) {
            protected_coordinates[i] = ExactIntervalType(-1.0_x, 1.0_x);
        }
        source.adjoin_parameters(protected_coordinates, EnclosureVariableKind::PARAMETER);
        const auto protected_result = clone->orbit(source, Real(0.0625_x), Semantics::UPPER);
        for (SizeType i = 1u; i != 7u; ++i) {
            ARIADNE_TESTING_EXPECT_EQUAL(protected_result.final()[0u].variable_kinds()[i], EnclosureVariableKind::PARAMETER);
        }
    }

    Void test_interval_bound_covers_the_selected_numerical_field() {
        struct Probe : ContinuousTimeEvolver<ValidatedTag> {
            using ContinuousTimeEvolver<ValidatedTag>::_inclusion_step;
            using ContinuousTimeEvolver<ValidatedTag>::ContinuousTimeEvolver;
        };
        RealVariable x("x");
        const ExactBoxType domain({
            {0.5_x, 1.5_x}
        });
        const ThresholdSweeper<FloatDP> sweeper(dp, 1e-10);
        using Model = ValidatedIntervalTaylorFunctionModel<FloatDP>;
        const auto uncertainty = Model::constant(domain, FloatDPUpperInterval(-1, 1, dp), sweeper);
        const auto coordinate = Model::coordinate(domain, 0u, sweeper);
        const auto one = Model::constant(domain, FloatDPUpperInterval(1, 1, dp), sweeper);
        const Model component = uncertainty * (coordinate - one);
        const OrdinaryDifferentialAlgebraicInclusion system(RealSpace({x}),
                                                            ValidatedVectorIntervalTaylorFunctionModel<FloatDP>({component}));
        TaylorPicardIntegrator base(step_maximum_error = 1e-6);
        Probe evolver(system, base, ZeroApproximation(), ConstantApproximation());
        const auto initial = evolver.enclosure(ExactBoxType({
            {1.0_x, 1.0_x}
        }));
        // Original c*(x-1) vanishes at x=1. Expanding it into independent
        // interval coefficients admits x'=x+1, whose tube leaves the model
        // domain before h=1/4. A bound on only the original field is unsound.
        const auto step = evolver._inclusion_step(initial, StepSizeType(0u), StepSizeType(0.25_x));
        ARIADNE_TESTING_PROVE(step.time_step() <= StepSizeType(0.125_x));
        ARIADNE_TESTING_PROVE(step.time_step() > StepSizeType(0u));
        ARIADNE_TESTING_EXPECT_EQUAL(evolver.system().dimension(), 1u);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_continuous_inclusions");
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_explicit_reach_and_endpoint_choices);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_reconditioning_has_no_winner_statistics);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_affine_dynamics_with_uncertain_inputs);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_switching_inputs_are_independent_between_steps);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_uniform_initial_error_is_preserved_without_reconditioning);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_timed_constant_parameters_survive_input_choices);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_autonomous_parameters_keep_the_declared_state);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_first_order_box_paths_and_zero_free_inputs);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_structured_lower_inclusion_preserves_switching_inputs);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_lower_preserves_unsplit_initial_enclosure);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_pointwise_interval_field_and_domain_rejection);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_owned_components_and_configuration);
    ARIADNE_TESTING_ADD_TEST(TestContinuousInclusions, test_interval_bound_covers_the_selected_numerical_field);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
