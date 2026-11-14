/***************************************************************************
 *            test_enclosure_recombination.cpp
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

/** @file tests/test_enclosure_recombination.cpp
 *  @brief Verify continuous union coverage, constraints and clock recombination.
 */

#include <limits>
#include <sstream>
#include <utility>

#include "algebra/sweeper.hpp"
#include "evolution/enclosure_recombiner.hpp"
#include "function/taylor_function.hpp"

#include "ariadne-testing.hpp"

using namespace Ariadne;

namespace {

// A separate implementation proves that the handle dispatches collection
// strategies, rather than always invoking the single-output free operation.
class SeparateEnclosureRecombiner : public EnclosureRecombinerInterface {
  public:
    SeparateEnclosureRecombiner* clone() const override { return new SeparateEnclosureRecombiner(*this); }
    EnclosureRecombinationResult recombine(const List<Enclosure>& enclosures,
                                           const EnclosureRecombinationContext& context) const override {
        EnclosureRecombinationResult result;
        for (SizeType index = 0; index != enclosures.size(); ++index) {
            auto source = enclosures[index];
            if (definitely(source.is_empty())) {
                result.discarded_sources.append(index);
                continue;
            }
            if (context.maximum_enclosure_radius())
                result.radius_limit_satisfied = result.radius_limit_satisfied && definitely(source.radius() <= *context.maximum_enclosure_radius());
            result.outputs.push_back({std::move(source), {index}});
        }
        return result;
    }
    OutputStream& _write(OutputStream& os) const override { return os << "SeparateEnclosureRecombiner()"; }
};

} // namespace

class TestEnclosureRecombination : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

  private:
    EnclosureConfiguration config{ValidatedFunctionPatchFactory(
        make_taylor_function_patch_factory(ThresholdSweeper<FloatDP>(dp, 1e-14)))};
    ExactBoxType domain{
        {-1, 1}
    };
    EffectiveScalarMultivariateFunction s = EffectiveScalarMultivariateFunction::coordinate(1u, 0u);
    const auto& factory() const { return config.function_factory(); }

    void expect_bounds(const UpperIntervalType& interval, double lower, double upper) {
        const auto actual_lower = interval.lower_bound().raw().get_d();
        const auto actual_upper = interval.upper_bound().raw().get_d();
        // Both inclusion and tightness matter: a large enclosing box cannot
        // pass the nonlinear-correlation and constrained-domain checks.
        ARIADNE_TESTING_PROVE(actual_lower <= lower);
        ARIADNE_TESTING_PROVE(actual_upper >= upper);
        ARIADNE_TESTING_PROVE(actual_lower >= lower - 1e-10);
        ARIADNE_TESTING_PROVE(actual_upper <= upper + 1e-10);
    }

    String description(const Enclosure& enclosure) {
        std::ostringstream stream;
        stream << enclosure;
        return stream.str();
    }

  public:
    void test_handle_dispatch_clone_lifetime_and_context() {
        auto handle = [] {
            SeparateEnclosureRecombiner implementation;
            const EnclosureRecombinerInterface& interface = implementation;
            return EnclosureRecombiner(interface);
        }();
        const auto copied = handle;
        Enclosure first(domain, EffectiveVectorMultivariateFunction{s, sqr(s)}, config);
        Enclosure second(domain, EffectiveVectorMultivariateFunction{s, sqr(s) + Real(1 / 4_q)}, config);
        const auto separate = copied.recombine({first, second}, EnclosureRecombinationContext(0.25_x));
        ARIADNE_TESTING_EXPECT_EQUAL(separate.outputs.size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(separate.outputs[0].source_indices, List<SizeType>({0u}));
        ARIADNE_TESTING_EXPECT_EQUAL(separate.outputs[1].source_indices, List<SizeType>({1u}));
        ARIADNE_TESTING_PROVE(!separate.radius_limit_satisfied);
        ARIADNE_TESTING_EXPECT_EQUAL(description(separate.outputs[0].enclosure), description(first));
        std::ostringstream stream;
        stream << handle;
        ARIADNE_TESTING_EXPECT_EQUAL(stream.str(), String("SeparateEnclosureRecombiner()"));
    }

    void test_reference_batch_provenance_and_empty_inputs() {
        Enclosure first(domain, EffectiveVectorMultivariateFunction{s, sqr(s)}, config);
        Enclosure second(domain, EffectiveVectorMultivariateFunction{s, sqr(s) + Real(1 / 4_q)}, config);
        auto empty = first;
        empty.new_zero_parameter_constraint(EffectiveScalarMultivariateFunction::constant(1u, Real(1)));
        ReferenceEnclosureRecombiner implementation;
        const EnclosureRecombiner handle(implementation);
        const auto result = handle.recombine({empty, first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.outputs.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.outputs.front().source_indices, List<SizeType>({1u, 2u}));
        ARIADNE_TESTING_EXPECT_EQUAL(result.discarded_sources, List<SizeType>({0u}));
        ARIADNE_TESTING_PROVE(result.radius_limit_satisfied);
        ARIADNE_TESTING_EXPECT_EQUAL(description(result.outputs.front().enclosure), description(Ariadne::recombine({first, second})));
        const auto none = handle.recombine({});
        ARIADNE_TESTING_EXPECT_EQUAL(none.outputs.size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(none.discarded_sources.size(), 0u);
        const auto exhausted = handle.recombine({empty, empty});
        ARIADNE_TESTING_EXPECT_EQUAL(exhausted.outputs.size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(exhausted.discarded_sources, List<SizeType>({0u, 1u}));
        ARIADNE_TESTING_PROVE(exhausted.radius_limit_satisfied);
        const auto singleton = handle.recombine({empty, first});
        ARIADNE_TESTING_EXPECT_EQUAL(singleton.outputs.front().source_indices, List<SizeType>({1u}));
        ARIADNE_TESTING_EXPECT_EQUAL(description(singleton.outputs.front().enclosure), description(first));
    }

    void test_radius_request_is_reported_without_losing_coverage() {
        Enclosure first(domain, EffectiveVectorMultivariateFunction{s, sqr(s)}, config);
        Enclosure second(domain, EffectiveVectorMultivariateFunction{s, sqr(s) + Real(1 / 4_q)}, config);
        const EnclosureRecombiner handle(ReferenceEnclosureRecombiner{});
        const EnclosureRecombinationContext unrestricted;
        ARIADNE_TESTING_PROVE(!unrestricted.maximum_enclosure_radius().has_value());
        const auto small = handle.recombine({first, second}, EnclosureRecombinationContext(0.25_x));
        ARIADNE_TESTING_PROVE(!small.radius_limit_satisfied);
        ARIADNE_TESTING_EXPECT_EQUAL(small.outputs.size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(small.outputs.front().source_indices, List<SizeType>({0u, 1u}));
        expect_bounds(apply(small.outputs.front().enclosure.state_function().get(1) - sqr(small.outputs.front().enclosure.state_function().get(0)), small.outputs.front().enclosure.domain()), 0, 0.25);
        const auto large = handle.recombine({first, second}, EnclosureRecombinationContext(4.0_x));
        ARIADNE_TESTING_PROVE(large.radius_limit_satisfied);
        ARIADNE_TESTING_EXPECT_EQUAL(description(small.outputs.front().enclosure), description(large.outputs.front().enclosure));
        Enclosure point(ExactBoxType({
                            {1, 1}
        }),
                        config);
        const auto exact = handle.recombine({point}, EnclosureRecombinationContext(0.0_x));
        ARIADNE_TESTING_PROVE(exact.radius_limit_satisfied);
        ARIADNE_TESTING_EXPECT_EXCEPTION(EnclosureRecombinationContext(-1.0_x));
        ARIADNE_TESTING_EXPECT_EXCEPTION(EnclosureRecombinationContext(ExactDouble(std::numeric_limits<double>::infinity())));
        ARIADNE_TESTING_EXPECT_EXCEPTION(EnclosureRecombinationContext(ExactDouble(std::numeric_limits<double>::quiet_NaN())));
    }

    void test_nonlinear_reference_and_full_source_coverage() {
        Enclosure first(domain, EffectiveVectorMultivariateFunction{s, sqr(s)}, config);
        Enclosure second(domain, EffectiveVectorMultivariateFunction{s, sqr(s) + Real(1 / 4_q)}, config);
        const auto first_before = description(first);
        const auto second_before = description(second);
        const auto result = recombine({first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.variable_kinds().back(), EnclosureVariableKind::ERROR);
        expect_bounds(apply(result.state_function().get(1) - sqr(result.state_function().get(0)), result.domain()), 0, 0.25);
        expect_bounds(apply(result.time_function(), result.domain()), 0, 0);
        expect_bounds(apply(result.dwell_time_function(), result.domain()), 0, 0);

        // Certify whole curves, rather than sampling their points: the result
        // is exactly (2*u-1,(2*u-1)^2+e), with u in [0,1], e in [0,1/4].
        const auto u = factory().create_coordinate(result.domain(), 0u);
        const auto e = factory().create_coordinate(result.domain(), 1u);
        expect_bounds(apply(result.state_function().get(0) - (2 * u - 1), result.domain()), 0, 0);
        expect_bounds(apply(result.state_function().get(1) - (sqr(2 * u - 1) + e), result.domain()), 0, 0);
        ARIADNE_TESTING_PROVE(same(result.domain(), ExactBoxType({
                                                        {0,      1},
                                                        {0, 0.25_x}
        })));
        ARIADNE_TESTING_EXPECT_EQUAL(description(first), first_before);
        ARIADNE_TESTING_EXPECT_EQUAL(description(second), second_before);
    }

    void test_disjoint_constraints_are_united_and_domains_contracted() {
        ExactBoxType wide{
            {-10, 10}
        };
        Enclosure left(wide, EffectiveVectorMultivariateFunction{s}, config);
        left.new_parameter_constraint(ValidatedConstraint(-1, s, Dyadic(-1 / 2_q)));
        Enclosure right(wide, EffectiveVectorMultivariateFunction{s}, config);
        right.new_parameter_constraint(ValidatedConstraint(Dyadic(1 / 2_q), s, 1));
        const auto left_domain = left.reduced_domain();
        const auto right_domain = right.reduced_domain();
        const auto result = recombine({left, right});
        ARIADNE_TESTING_PROVE(!definitely(result.is_empty()));
        expect_bounds(apply(result.state_function().get(0), result.reduced_domain()), -1, 1);
        ARIADNE_TESTING_PROVE(same(left.reduced_domain(), left_domain));
        ARIADNE_TESTING_PROVE(same(right.reduced_domain(), right_domain));

        // Reverse input order: the reference may change, coverage may not.
        const auto reversed = recombine({right, left});
        expect_bounds(apply(reversed.state_function().get(0), reversed.reduced_domain()), -1, 1);
    }

    void test_oblique_common_constraint_and_different_constraint_lists() {
        ExactBoxType square{
            {-1, 1},
            {-1, 1}
        };
        auto coordinates = EffectiveVectorMultivariateFunction::identity(2u);
        auto sum = coordinates[0] + coordinates[1];
        Enclosure first(square, coordinates, config);
        first.new_zero_parameter_constraint(sum);
        Enclosure second(square, coordinates, config);
        second.new_parameter_constraint(ValidatedConstraint(-2, coordinates[0], 2));
        second.new_zero_parameter_constraint(sum);
        const auto result = recombine({first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_constraints(), 1u);
        expect_bounds(UpperIntervalType(result.constraint(0).bounds()), 0, 0);
        ARIADNE_TESTING_PROVE(same(result.reduced_domain(), result.domain()));
        ARIADNE_TESTING_PROVE(!definitely(result.is_empty()));

        // An unconstrained alternative must not inherit the first predicate.
        const auto relaxed = recombine({first, Enclosure(square, coordinates, config)});
        expect_bounds(UpperIntervalType(relaxed.constraint(0).bounds()), -2, 2);
        const auto reordered = recombine({second, first});
        ARIADNE_TESTING_PROVE(!definitely(reordered.is_empty()));
        expect_bounds(UpperIntervalType(reordered.constraint(1).bounds()), 0, 0);
    }

    void test_constraints_reduce_residuals_without_reducing_the_box() {
        ExactBoxType square{
            {-1, 1},
            {-1, 1}
        };
        const auto coordinates = EffectiveVectorMultivariateFunction::identity(2u);
        const auto x = coordinates[0];
        const auto sum = coordinates[0] + coordinates[1];
        Enclosure first(square, EffectiveVectorMultivariateFunction{x, sqr(x)}, config);
        first.new_zero_parameter_constraint(sum);
        Enclosure second(square, EffectiveVectorMultivariateFunction{x, sqr(x) + sum}, config);
        second.new_zero_parameter_constraint(sum);
        const auto result = recombine({first, second});
        // Every source point has u+v=0, so the apparently different state maps
        // describe the same whole parabola. Boxing the domains alone would
        // introduce an artificial residual of width four.
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_parameters(), 2u);
        expect_bounds(apply(result.state_function().get(1) - sqr(result.state_function().get(0)), result.domain()), 0, 0);
        expect_bounds(UpperIntervalType(result.constraint(0).bounds()), 0, 0);
    }

    void test_clock_maps_and_continued_evolution() {
        ExactBoxType unit{
            {0, 1}
        };
        Enclosure first(unit, EffectiveVectorMultivariateFunction{s, sqr(s)}, s, {}, config);
        Enclosure second(unit, EffectiveVectorMultivariateFunction{s, sqr(s)}, s + 1, {}, config);
        const auto flow_coordinates = EffectiveVectorMultivariateFunction::identity(3u);
        EffectiveVectorMultivariateFunction stationary_flow{flow_coordinates[0], flow_coordinates[1]};
        first.apply_parameter_evolve_step(stationary_flow, s);
        second.apply_parameter_evolve_step(stationary_flow, s + Real(1 / 4_q));
        auto result = recombine({first, second});
        expect_bounds(apply(result.time_function() - 2 * result.state_function().get(0), result.domain()), 0, 1.25);
        expect_bounds(apply(result.dwell_time_function() - result.state_function().get(0), result.domain()), 0, 0.25);

        // Exact nonlinear flow of dx/dt=1, dy/dt=x^2. The result must remain
        // usable by ordinary Enclosure evolution, preserving both clocks.
        const auto x = flow_coordinates[0];
        const auto h = flow_coordinates[2];
        EffectiveVectorMultivariateFunction nonlinear_flow{x + h,
                                                           flow_coordinates[1] + (sqr(x + h) * (x + h) - sqr(x) * x) / 3};
        const auto flow_domain = product(cast_exact_box(result.bounding_box()), ExactIntervalType(0, 0.5_x));
        result.apply_fixed_evolve_step(factory().create(flow_domain, nonlinear_flow), Dyadic(1 / 2_q));
        expect_bounds(apply(result.time_function() - 2 * result.state_function().get(0), result.domain()), -0.5, 0.75);
        expect_bounds(apply(result.dwell_time_function() - result.state_function().get(0), result.domain()), 0, 0.25);
        const auto initial_x = result.state_function().get(0) - FloatDPBounds(0.5_x, dp);
        expect_bounds(apply(result.state_function().get(1) - (FloatDPBounds(1.5_x, dp) * sqr(initial_x) + FloatDPBounds(0.25_x, dp) * initial_x + FloatDPBounds(Real(1 / 24_q), dp)), result.domain()), 0, 0);

        // Finish all the differently timed sources at one absolute horizon.
        // The same exact flow gives y_final=s^2+(x_final^3-s^3)/3.
        const auto finishing_domain = product(cast_exact_box(result.bounding_box()), ExactIntervalType(0, 4));
        result.apply_finishing_parameter_evolve_step(factory().create(finishing_domain, nonlinear_flow),
                                                     EffectiveScalarMultivariateFunction::constant(result.number_of_parameters(), Real(4)));
        expect_bounds(apply(result.time_function(), result.domain()), 4, 4);
        expect_bounds(apply(result.dwell_time_function() - result.state_function().get(0), result.domain()), 0, 0.25);
        const auto original_x = factory().create_coordinate(result.domain(), 0u);
        const auto final_x = result.state_function().get(0);
        expect_bounds(apply(result.state_function().get(1) - sqr(original_x) - (sqr(final_x) * final_x - sqr(original_x) * original_x) * FloatDPBounds(Real(1 / 3_q), dp), result.domain()), 0, 0);
    }

    // Full-domain algebraic witness for joint state/time/dwell coverage.
    // Sources share u in [0,1] and differ by exact constant offsets. The
    // recombined enclosure has independent residual coordinates containing
    // BOTH offset triples, so each whole source curve embeds in its graph.
    // No sampled points, tolerance-based point matching or source boxing.
    void test_exact_joint_clock_and_state_coverage_witness() {
        const ExactBoxType unit{{0, 1}};
        Enclosure first(unit, EffectiveVectorMultivariateFunction{s, sqr(s)}, s, {}, config);
        Enclosure second(unit, EffectiveVectorMultivariateFunction{s, sqr(s) + Real(1 / 4_q)},
                         s + Real(1 / 2_q), {}, config);
        const auto coords = EffectiveVectorMultivariateFunction::identity(3u);
        const EffectiveVectorMultivariateFunction stationary{coords[0], coords[1]};
        first.apply_parameter_evolve_step(stationary, s);
        second.apply_parameter_evolve_step(stationary, s + Real(1 / 4_q));
        // Joint source maps:
        // first (u,u^2,2u,u), second (u,u^2+1/4,2u+3/4,u+1/4).
        const auto result = recombine({first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_parameters(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_constraints(), 0u);
        const auto u = factory().create_coordinate(result.domain(), 0u);
        const auto ey = factory().create_coordinate(result.domain(), 1u);
        const auto et = factory().create_coordinate(result.domain(), 2u);
        const auto ed = factory().create_coordinate(result.domain(), 3u);
        expect_bounds(apply(result.state_function().get(0) - u, result.domain()), 0, 0);
        expect_bounds(apply(result.state_function().get(1) - (sqr(u) + ey), result.domain()), 0, 0);
        expect_bounds(apply(result.time_function() - (2*u + et), result.domain()), 0, 0);
        expect_bounds(apply(result.dwell_time_function() - (u + ed), result.domain()), 0, 0);
        // For every u, either constant residual triple (0,0,0) or
        // (1/4,3/4,1/4) is in this independent residual box.
        expect_bounds(UpperIntervalType(result.domain()[1]), 0, 0.25);
        expect_bounds(UpperIntervalType(result.domain()[2]), 0, 0.75);
        expect_bounds(UpperIntervalType(result.domain()[3]), 0, 0.25);
    }

    // An algebraic whole-domain coverage witness for two urgent-guard curves.
    // The two source maps satisfy x+y=0 identically while their absolute
    // birth times differ. One shared reference parameter and independent
    // residuals preserve both entire curves and their clock correlations.
    void test_guard_constrained_joint_coverage_witness() {
        const ExactBoxType unit{{0, 1}};
        Enclosure first(unit, EffectiveVectorMultivariateFunction{-s, s}, s, {}, config);
        Enclosure second(unit, EffectiveVectorMultivariateFunction{
                             -s - Real(1 / 4_q), s + Real(1 / 4_q)},
                         s + Real(1 / 4_q), {}, config);
        const auto coordinates = EffectiveVectorMultivariateFunction::identity(2u);
        first.new_zero_state_constraint(coordinates[0] + coordinates[1]);
        second.new_zero_state_constraint(coordinates[0] + coordinates[1]);

        const auto result = recombine({first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_parameters(), 4u);
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_constraints(), 1u);
        const auto u = factory().create_coordinate(result.domain(), 0u);
        const auto ex = factory().create_coordinate(result.domain(), 1u);
        const auto ey = factory().create_coordinate(result.domain(), 2u);
        const auto et = factory().create_coordinate(result.domain(), 3u);
        expect_bounds(apply(result.state_function().get(0) + u - ex, result.domain()), 0, 0);
        expect_bounds(apply(result.state_function().get(1) - u - ey, result.domain()), 0, 0);
        expect_bounds(apply(result.time_function() - u - et, result.domain()), 0, 0);
        expect_bounds(apply(result.dwell_time_function(), result.domain()), 0, 0);

        // For every source parameter u in [0,1] choose the constant triple
        // (ex,ey,et)=(0,0,0) or (-1/4,1/4,1/4) respectively.
        // Both choices satisfy the output's guard equality and its bounds.
        expect_bounds(UpperIntervalType(result.domain()[1]), -0.25, 0);
        expect_bounds(UpperIntervalType(result.domain()[2]), 0, 0.25);
        expect_bounds(UpperIntervalType(result.domain()[3]), 0, 0.25);
        expect_bounds(UpperIntervalType(result.constraint(0).bounds()), 0, 0);
    }

    void test_unequal_parameter_dimensions() {
        Enclosure first(domain, EffectiveVectorMultivariateFunction{s, sqr(s)}, config);
        auto coordinates = EffectiveVectorMultivariateFunction::identity(2u);
        ExactBoxType rectangle{
            {     -1,      1},
            {-0.25_x, 0.25_x}
        };
        Enclosure second(rectangle,
                         EffectiveVectorMultivariateFunction{coordinates[0] + coordinates[1], sqr(coordinates[0]) + coordinates[1]}, config);
        const auto result = recombine({first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_parameters(), 4u);
        const auto u = factory().create_coordinate(result.domain(), 0u);
        const auto ex = factory().create_coordinate(result.domain(), 2u);
        const auto ey = factory().create_coordinate(result.domain(), 3u);
        expect_bounds(apply(result.state_function().get(0) - (2 * u - 1 + ex), result.domain()), 0, 0);
        expect_bounds(apply(result.state_function().get(1) - (sqr(2 * u - 1) + ey), result.domain()), 0, 0);
        expect_bounds(UpperIntervalType(result.domain()[2]), -0.25, 0.25);
        expect_bounds(UpperIntervalType(result.domain()[3]), -0.25, 0.25);
    }

    void test_points_and_zero_dimensional_parameters() {
        ExactBoxType point_domain(0u);
        Enclosure first(point_domain, EffectiveVectorMultivariateFunction::constant(0u, Vector<Real>(1u, Real(1))), config);
        Enclosure second(point_domain, EffectiveVectorMultivariateFunction::constant(0u, Vector<Real>(1u, Real(2))), config);
        const auto result = recombine({first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.number_of_parameters(), 1u);
        expect_bounds(apply(result.state_function().get(0), result.domain()), 1, 2);
        expect_bounds(apply(result.time_function(), result.domain()), 0, 0);
        expect_bounds(apply(result.dwell_time_function(), result.domain()), 0, 0);

        Enclosure degenerate(ExactBoxType({
                                 {0.25_x, 0.25_x}
        }),
                             config);
        Enclosure other(ExactBoxType({
                            {0.5_x, 0.5_x}
        }),
                        config);
        const auto degenerate_result = recombine({degenerate, other});
        expect_bounds(apply(degenerate_result.state_function().get(0), degenerate_result.domain()), 0.25, 0.5);
    }

    void test_uniform_errors_and_rounded_parameter_maps() {
        auto uncertain = factory().create(domain, EffectiveVectorMultivariateFunction{s, sqr(s)});
        uncertain[1] = uncertain.get(1) + FloatDPBounds(-0.125_x, 0.125_x, dp);
        Enclosure first(domain, uncertain, config);
        Enclosure second(domain, EffectiveVectorMultivariateFunction{s, sqr(s) + Real(1 / 4_q)}, config);
        const auto result = recombine({first, second});
        const auto difference = apply(result.state_function().get(1) - sqr(result.state_function().get(0)), result.domain());
        ARIADNE_TESTING_PROVE(difference.lower_bound().raw().get_d() <= -0.125);
        ARIADNE_TESTING_PROVE(difference.upper_bound().raw().get_d() >= 0.25);

        Enclosure rounded_first(ExactBoxType({
                                    {ExactDouble(0.1), ExactDouble(0.3)}
        }),
                                EffectiveVectorMultivariateFunction{sqr(s)}, config);
        Enclosure rounded_second(ExactBoxType({
                                     {ExactDouble(0.2), ExactDouble(0.4)}
        }),
                                 EffectiveVectorMultivariateFunction{sqr(s)}, config);
        const auto rounded_result = recombine({rounded_first, rounded_second});
        ARIADNE_TESTING_PROVE(!definitely(rounded_result.is_empty()));
        const auto bounds = apply(rounded_result.state_function().get(0), rounded_result.domain());
        ARIADNE_TESTING_PROVE(bounds.lower_bound().raw().get_d() <= 0.1 * 0.1);
        ARIADNE_TESTING_PROVE(bounds.upper_bound().raw().get_d() >= 0.4 * 0.4);
    }

    void test_empty_singleton_and_duplicate_sources() {
        Enclosure first(domain, EffectiveVectorMultivariateFunction{s, sqr(s)}, config);
        Enclosure empty = first;
        empty.new_zero_parameter_constraint(EffectiveScalarMultivariateFunction::constant(1u, Real(1)));
        ARIADNE_TESTING_PROVE(definitely(recombine({empty, empty}).is_empty()));
        ARIADNE_TESTING_EXPECT_EQUAL(description(recombine({first})), description(first));
        ARIADNE_TESTING_EXPECT_EQUAL(description(recombine({empty, first})), description(first));
        const auto duplicates = recombine({first, first});
        ARIADNE_TESTING_EXPECT_EQUAL(duplicates.number_of_parameters(), 1u);
        expect_bounds(apply(duplicates.state_function().get(1) - sqr(duplicates.state_function().get(0)), duplicates.domain()), 0, 0);
    }

    void test_input_contracts_and_auxiliary_mapping() {
        Enclosure first(domain, config);
        ARIADNE_TESTING_EXPECT_EXCEPTION(recombine({}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(recombine({Enclosure()}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(recombine({first, Enclosure(config)}));
        ARIADNE_TESTING_EXPECT_EXCEPTION(recombine({
            first, Enclosure(ExactBoxType({{-1, 1}, {-1, 1}}
             ), config)
        }));
        first.set_auxiliary_mapping(EffectiveVectorMultivariateFunction{sqr(s)});
        Enclosure second(ExactBoxType({
                             {0, 2}
        }),
                         config);
        second.set_auxiliary_mapping(EffectiveVectorMultivariateFunction{sqr(s)});
        const auto result = recombine({first, second});
        ARIADNE_TESTING_EXPECT_EQUAL(result.auxiliary_mapping().result_size(), 1u);
        expect_bounds(apply(result.auxiliary_function().get(0) - sqr(result.state_function().get(0)), result.domain()), 0, 0);
        second.set_auxiliary_mapping(EffectiveVectorMultivariateFunction{s});
        ARIADNE_TESTING_EXPECT_EXCEPTION(recombine({first, second}));
    }
};

int main(int argc, char* argv[]) {
    ARIADNE_TESTING_SUITE("test_enclosure_recombination");
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_handle_dispatch_clone_lifetime_and_context);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_reference_batch_provenance_and_empty_inputs);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_radius_request_is_reported_without_losing_coverage);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_nonlinear_reference_and_full_source_coverage);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_disjoint_constraints_are_united_and_domains_contracted);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_oblique_common_constraint_and_different_constraint_lists);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_constraints_reduce_residuals_without_reducing_the_box);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_clock_maps_and_continued_evolution);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_exact_joint_clock_and_state_coverage_witness);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_guard_constrained_joint_coverage_witness);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_unequal_parameter_dimensions);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_points_and_zero_dimensional_parameters);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_uniform_errors_and_rounded_parameter_maps);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_empty_singleton_and_duplicate_sources);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecombination, test_input_contracts_and_auxiliary_mapping);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
