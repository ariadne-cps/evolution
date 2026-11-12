/***************************************************************************
 *            test_hybrid_enclosure.cpp
 *
 *  Copyright  2008-21  Luca Geretti
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

/** @file tests/test_hybrid_enclosure.cpp
 *  @brief Test hybrid enclosure behavior and numerical invariants.
 */

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>

#include "algebra/algebra.hpp"
#include "drawing/drawer.hpp"
#include "drawing/drawer_manager.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_grid.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_space.hpp"
#include "labelled/figure.hpp"
#include "recording_canvas.hpp"
#include "test_logging.hpp"

#include "ariadne-testing.hpp"

using namespace std;
using namespace Ariadne;

class TestHybridEnclosure : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
    void test_time_maps_and_steps() {
        RealVariable x("x"), y("y");
        DiscreteLocation q(1);
        ExactBoxType domain({
            {-1, 1},
            {-1, 1}
        });
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
        const auto& factory = config.function_factory();
        HybridEnclosure initial(HybridExactBox(q, RealSpace({x, y}), domain), config);
        auto coordinates = EffectiveVectorMultivariateFunction::identity(3u);
        EffectiveVectorMultivariateFunction flow = {coordinates[0] + coordinates[2], coordinates[1] + 2 * coordinates[2]};
        ExactBoxType flow_domain = product(domain, ExactIntervalType(0.0_x, 0.5_x));
        auto phi = factory.create(flow_domain, flow);
        auto constant_space = EffectiveScalarMultivariateFunction::constant(2u, Real(1 / 4_q));
        auto constant_spacetime = EffectiveScalarMultivariateFunction::constant(3u, Real(1 / 4_q));
        auto space = initial;
        space.apply_space_evolve_step(phi, factory.create(domain, constant_space));
        ARIADNE_TESTING_PROVE(definitely(space.state_bounding_box().inside(ExactBoxType({
            {-ExactDouble(0.76), ExactDouble(1.26)},
            {-ExactDouble(0.51), ExactDouble(1.51)}
        }))));
        ARIADNE_TESTING_EXPECT_EQUAL(space.time_range().upper_bound().raw().get_d(), 0.25);
        auto spacetime = initial;
        spacetime.apply_spacetime_evolve_step(phi, factory.create(flow_domain, constant_spacetime));
        ARIADNE_TESTING_EXPECT_EQUAL(spacetime.time_range().upper_bound().raw().get_d(), 0.25);
        auto reach = initial;
        reach.apply_spacetime_reach_step(phi, factory.create(flow_domain, constant_spacetime));
        ARIADNE_TESTING_EXPECT_EQUAL(reach.number_of_parameters(), 3u);
        ARIADNE_TESTING_PROVE(reach.time_range().upper_bound().raw().get_d() >= 0.25);
        reach.reduce();
        ARIADNE_TESTING_PROVE(!definitely(reach.is_empty()));
        auto changed = space;
        changed.set_time_function(factory.create_constant(domain, Real(2)));
        ARIADNE_TESTING_EXPECT_EQUAL(changed.time_range().lower_bound().raw().get_d(), 2.0);
        ARIADNE_TESTING_EXPECT_EQUAL(changed.time_range().upper_bound().raw().get_d(), 2.0);
        ARIADNE_TESTING_EXPECT_EQUAL(changed.dwell_time_range().upper_bound().raw().get_d(), 0.25);
        ARIADNE_TESTING_EXPECT_EXCEPTION(changed.set_time_function(factory.create_zero(ExactBoxType({
            {0, 1}
        }))));
        ARIADNE_TESTING_EXPECT_EXCEPTION(changed.set_time_function(factory.create_zero(ExactBoxType({
            {0, 1},
            {0, 1}
        }))));
        changed.clear_time();
        ARIADNE_TESTING_EXPECT_EQUAL(changed.time_range().upper_bound().raw().get_d(), 0.0);
        ARIADNE_TESTING_EXPECT_EQUAL(changed.dwell_time_range().upper_bound().raw().get_d(), 0.0);
        RealVariable t("t");
        HybridEnclosure clock(HybridExactBox(q, RealSpace({t, x}), domain), config);
        ARIADNE_TESTING_EXPECT_EQUAL(clock.state_time_set().dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(clock.state_time_auxiliary_function().result_size(), 2u);
        RecordingCanvas canvas;
        clock.draw(canvas, {}, Variables2d(t, x));
        ARIADNE_TESTING_PROVE(canvas.boundaries > 0u);
    }

    void test_constraints_and_collections() {
        RealVariable x("x"), y("y"), z("z");
        RealSpace space({x, y});
        DiscreteLocation q(1), other(2);
        DiscreteEvent event("event");
        ExactBoxType domain({
            {-1, 1},
            {-1, 1}
        });
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
        HybridEnclosure initial(HybridBoxSet(q, {-1 <= x <= 1, -1 <= y <= 1}), space, config);
        HybridEnclosure exact(HybridExactBox(q, space, domain), config);
        std::unique_ptr<HybridEnclosure> cloned(exact.clone());
        ARIADNE_TESTING_EXPECT_EQUAL(cloned->parameter_domain(), domain);
        auto coordinate = EffectiveVectorMultivariateFunction::identity(2u);
        auto range = initial.range_of(coordinate[0] + 2 * coordinate[1]);
        ARIADNE_TESTING_PROVE(range.lower_bound().raw().get_d() <= -3.0);
        ARIADNE_TESTING_PROVE(range.upper_bound().raw().get_d() >= 3.0);
        auto bounded = initial;
        bounded.new_state_time_bound(event, EffectiveScalarMultivariateFunction::constant(2u, Real(1)));
        bounded.reduce();
        ARIADNE_TESTING_PROVE(!definitely(bounded.is_empty()));
        auto impossible = initial;
        impossible.new_state_time_bound(event, EffectiveScalarMultivariateFunction::constant(2u, Real(-1)));
        impossible.reduce();
        ARIADNE_TESTING_PROVE(definitely(impossible.is_empty()));
        auto constrained = initial;
        constrained.new_guard(event, coordinate[0]);
        constrained.new_state_constraint(event, ValidatedConstraint(-1, coordinate[1], 1));
        auto spacetime = EffectiveVectorMultivariateFunction::identity(3u);
        constrained.new_state_time_constraint(event, ValidatedConstraint(-2, spacetime[0] + spacetime[2], 2));
        ARIADNE_TESTING_EXPECT_EQUAL(constrained.number_of_constraints(), 3u);
        constrained.bound_time(ValidatedScalarMultivariateFunction(EffectiveScalarMultivariateFunction::constant(2u, Real(1))));
        constrained.set_time(ValidatedScalarMultivariateFunction(EffectiveScalarMultivariateFunction::zero(2u)));
        constrained.reduce();
        ARIADNE_TESTING_PROVE(!definitely(constrained.is_empty()));
        ARIADNE_TESTING_EXPECT_EXCEPTION(initial.set_auxiliary({z}, EffectiveVectorMultivariateFunction::zeros(1u, 1u)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(initial.set_auxiliary({}, EffectiveVectorMultivariateFunction::zeros(1u, 2u)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(initial.apply_reset(event, other, space, EffectiveVectorMultivariateFunction::identity(1u), RealSpace(), EffectiveVectorMultivariateFunction::zeros(0u, 2u)));
        auto reset = initial;
        reset.apply_reset(event, other, space, coordinate, RealSpace(), EffectiveVectorMultivariateFunction::zeros(0u, 2u));
        ARIADNE_TESTING_EXPECT_EQUAL(reset.previous_events().size(), 1u);
        reset.clear_events();
        ARIADNE_TESTING_EXPECT_EQUAL(reset.previous_events().size(), 0u);
        HybridExactBox enclosing(q, space, ExactBoxType({
                                               {-2, 2},
                                               {-2, 2}
        }));
        HybridExactBox disjoint(q, space, ExactBoxType({
                                              {3, 4},
                                              {3, 4}
        }));
        HybridExactBox elsewhere(other, space, domain);
        ARIADNE_TESTING_PROVE(definitely(initial.inside(enclosing)));
        ARIADNE_TESTING_PROVE(!definitely(initial.inside(elsewhere)));
        ARIADNE_TESTING_PROVE(definitely(separated(initial, disjoint)));
        ARIADNE_TESTING_PROVE(definitely(separated(initial, elsewhere)));
        ARIADNE_TESTING_PROVE(!definitely(separated(initial, enclosing)));
        auto split = initial.split();
        ARIADNE_TESTING_EXPECT_EQUAL(split.size(), 2u);
        for (auto& part: split) {
            ARIADNE_TESTING_PROVE(definitely(part.inside(enclosing)));
            part.uniform_error_recondition();
            part.kuhn_recondition();
            part.recondition();
            part._check();
        }
        ListSet<HybridEnclosure> sets(initial), empty;
        sets.append(reset);
        sets.adjoin(ListSet<HybridEnclosure>(List<HybridEnclosure>{constrained}));
        sets.adjoin(empty);
        ARIADNE_TESTING_EXPECT_EQUAL(sets.size(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(sets.front().location(), q);
        ARIADNE_TESTING_EXPECT_EQUAL(sets.back().location(), q);
        const auto& const_sets = sets;
        ARIADNE_TESTING_EXPECT_EQUAL(const_sets.front().location(), q);
        ARIADNE_TESTING_EXPECT_EQUAL(const_sets.back().location(), q);
        ARIADNE_TESTING_EXPECT_EQUAL(sets[q].size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(sets[other].size(), 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(sets[DiscreteLocation(9)].size(), 0u);
        List<HybridEnclosure>& entries = sets;
        const List<HybridEnclosure>& const_entries = const_sets;
        ARIADNE_TESTING_EXPECT_EQUAL(entries.size(), const_entries.size());
        SizeType count = 0u;
        for (auto iter = sets.begin(); iter != sets.end(); ++iter) {
            ++count;
        }
        ARIADNE_TESTING_EXPECT_EQUAL(count, 3u);
        sets.reduce();
        empty.reduce();
        ARIADNE_TESTING_EXPECT_EQUAL(sets.bounding_box().euclidean_set(space).dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(empty.bounding_box());
        struct DrawerScope {
            Drawer previous{DrawerManager::instance().drawer()};
            DrawerScope() { DrawerManager::instance().set_drawer(BoxDrawer()); }
            ~DrawerScope() { DrawerManager::instance().set_drawer(previous); }
        } drawer_scope;
        RecordingCanvas canvas;
        sets.draw(canvas, {DiscreteLocation(9)}, Variables2d(x, y));
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.boundaries, 0u);
        sets.draw(canvas, {q}, Variables2d(x, y));
        ARIADNE_TESTING_PROVE(canvas.boundaries > 0u);
        sets.draw(canvas, {}, Variables2d(x, y));
        std::ostringstream brief, full;
        initial.print(brief);
        constrained.repr(full);
        ARIADNE_TESTING_PROVE(brief.str().find("time_range") != std::string::npos);
        ARIADNE_TESTING_PROVE(full.str().find("HybridEnclosure") != std::string::npos);
        MonolithicHybridSpace hybrid_space;
        hybrid_space.new_location(q, space);
        HybridPaving paving{HybridGrid(HybridSpace(hybrid_space))};
        initial.adjoin_outer_approximation_to(paving, 0u);
        ARIADNE_TESTING_PROVE(paving[q].size() > 0u);
        MonolithicHybridSpace incompatible;
        incompatible.new_location(q, RealSpace({y, x}));
        HybridPaving wrong{HybridGrid(HybridSpace(incompatible))};
        ARIADNE_TESTING_EXPECT_EXCEPTION(initial.adjoin_outer_approximation_to(wrong, 0u));
        sets.clear();
        ARIADNE_TESTING_EXPECT_EQUAL(sets.size(), 0u);
    }

    void test_alternate_factory_formatting() {
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        DiscreteLocation q(1);
        EnclosureConfiguration config(ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelDP>(ThresholdSweeper<FloatDP>(dp, 1e-8)));
        config.set_reconditioning_num_blocks(1u);
        HybridEnclosure enclosure(HybridExactBox(q, space, ExactBoxType({
                                                               {-1, 1},
                                                               {-1, 1}
        })),
                                  config);
        auto coordinate = EffectiveVectorMultivariateFunction::identity(2u);
        enclosure.new_state_constraint(DiscreteEvent("bound"), ValidatedConstraint(-1, coordinate[0], 1));
        enclosure._check();
        enclosure.kuhn_recondition();
        ARIADNE_TESTING_EXPECT_EQUAL(enclosure.number_of_parameters(), 2u);
        std::ostringstream representation;
        enclosure.repr(representation);
        ARIADNE_TESTING_PROVE(representation.str().find("HybridEnclosure") != std::string::npos);
        ARIADNE_TESTING_PROVE(representation.str().find("<=") != std::string::npos);
        ARIADNE_TESTING_EXPECT_EQUAL(enclosure.bounding_box().location(), q);
        ARIADNE_TESTING_PROVE(definitely(inside(enclosure, HybridRealBox(q, {-2 <= x <= 2, -2 <= y <= 2}))));
        HybridExactBox exact(q, space, ExactBoxType({
                                           {-2, 2},
                                           {-2, 2}
        }));
        ARIADNE_TESTING_PROVE(definitely(inside(enclosure, exact)));
        ARIADNE_TESTING_PROVE(definitely(inside(ListSet<HybridEnclosure>(enclosure), exact)));
        ARIADNE_TESTING_PROVE(definitely(inside(ListSet<HybridEnclosure>(), exact)));
    }

    void test_bounds_factory_reach_steps() {
        // Regression for ariadne#811: interval coefficients must support the
        // extra reach parameter, as well as a fixed-time endpoint.
        RealVariable x("x"), y("y");
        RealSpace space({x, y});
        DiscreteLocation q(1);
        EnclosureConfiguration settings(ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelDP>(ThresholdSweeper<FloatDP>(dp, 1e-8)));
        auto factory = settings.function_factory();
        ExactBoxType domain({
            {-1, 1},
            {-1, 1}
        });
        auto coordinates = EffectiveVectorMultivariateFunction::identity(3u);
        EffectiveVectorMultivariateFunction affine = {coordinates[0] + coordinates[2], coordinates[1] - 2 * coordinates[2]};
        auto flow = factory.create(product(domain, ExactIntervalType(0, 0.25_x)), affine);
        HybridEnclosure initial(HybridExactBox(q, space, domain), settings);
        auto fixed = initial, reach = initial;
        fixed.apply_fixed_evolve_step(flow, Dyadic(1 / 4_q));
        reach.apply_full_reach_step(flow);
        ARIADNE_TESTING_EXPECT_EQUAL(fixed.number_of_parameters(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(reach.number_of_parameters(), 3u);
        ARIADNE_TESTING_EXPECT_EQUAL(reach.state_dimension(), 2u);
        ARIADNE_TESTING_PROVE(definitely(inside(fixed, HybridExactBox(q, space,
                                                                      ExactBoxType({
                                                                          {-ExactDouble(0.751), ExactDouble(1.251)},
                                                                          {-ExactDouble(1.501), ExactDouble(0.501)}
        })))));
        ARIADNE_TESTING_PROVE(definitely(inside(reach, HybridExactBox(q, space,
                                                                      ExactBoxType({
                                                                          {-ExactDouble(1.001), ExactDouble(1.251)},
                                                                          {-ExactDouble(1.501), ExactDouble(1.001)}
        })))));
        auto at_end = reach.state_function()(Vector<FloatDPApproximation>({0, 0, 0.25_x}, dp));
        ARIADNE_TESTING_PROVE(std::abs(at_end[0].raw().get_d() - 0.25) < 1e-10);
        ARIADNE_TESTING_PROVE(std::abs(at_end[1].raw().get_d() + 0.5) < 1e-10);
        ARIADNE_TESTING_PROVE(std::abs(fixed.time_range().lower_bound().raw().get_d() - 0.25) < 1e-10);
        ARIADNE_TESTING_PROVE(std::abs(fixed.time_range().upper_bound().raw().get_d() - 0.25) < 1e-10);
        ARIADNE_TESTING_PROVE(std::abs(reach.time_range().lower_bound().raw().get_d()) < 1e-10);
        ARIADNE_TESTING_PROVE(std::abs(reach.time_range().upper_bound().raw().get_d() - 0.25) < 1e-10);
        // A curved flow exercises interval-coefficient Taylor multiplication
        // and composition, with the independent analytic solution x(t)=x0*exp(t).
        ExactBoxType positive_domain({
            { 1, 2},
            {-1, 1}
        });
        EffectiveVectorMultivariateFunction exponential = {coordinates[0] * exp(coordinates[2]), coordinates[1]};
        auto curved_flow = factory.create(product(positive_domain, ExactIntervalType(0, 0.25_x)), exponential);
        HybridEnclosure curved(HybridExactBox(q, space, positive_domain), settings);
        curved.apply_full_reach_step(curved_flow);
        ARIADNE_TESTING_PROVE(definitely(inside(curved, HybridExactBox(q, space,
                                                                       ExactBoxType({
                                                                           {                  0,                  3},
                                                                           {-ExactDouble(1.001), ExactDouble(1.001)}
        })))));
        auto value = curved.state_function()(Vector<FloatDPBounds>({1.5_x, 0.0_x, 0.25_x}, dp));
        double expected = 1.5 * std::exp(0.25);
        ARIADNE_TESTING_PROVE(value[0].lower().raw().get_d() <= expected + 1e-12);
        ARIADNE_TESTING_PROVE(value[0].upper().raw().get_d() >= expected - 1e-12);
        ARIADNE_TESTING_PROVE(value[1].lower().raw().get_d() <= 0.0);
        ARIADNE_TESTING_PROVE(value[1].upper().raw().get_d() >= 0.0);
    }

    void test_construct() {
        HybridEnclosure encl;

        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp, 1e-8));
        EnclosureConfiguration config(function_factory);
        HybridEnclosure encl2(config);
        ARIADNE_TESTING_NOTE(encl2.configuration());

        DiscreteLocation location;
        RealVariable x("x"), y("y");
        HybridRealBox box(location, {1 <= x <= 2, 0 <= y <= 1});
        HybridEnclosure encl3(box, config);
        ARIADNE_TESTING_EXPECT_EQUAL(encl3.number_of_parameters(), 2);
        ARIADNE_TESTING_EXPECT_EQUAL(encl3.number_of_constraints(), 0);

        Enclosure cencl(config);
        HybridEnclosure(location, RealSpace({x, y}), cencl);
    }

    void test_space_having_time() {
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp, 1e-8));
        EnclosureConfiguration config(function_factory);
        RealVariable x("x"), t("t");
        HybridRealBox box(DiscreteLocation(), {1 <= t <= 2, 0 <= x <= 1});
        HybridEnclosure encl(box, config);
        ARIADNE_TESTING_PROVE(encl.state_time_auxiliary_space() == RealSpace({t, x}));
    }

    void test_auxiliary() {
        ExactBoxType dom({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        Enclosure encl(dom, EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8))));
        ARIADNE_TESTING_NOTE(encl.auxiliary_function());
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2, 0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2, 1);
        EffectiveVectorMultivariateFunction auxiliary(1, x0 + sqr(x1));
        RealVariable x("x"), y("y"), z("z");
        HybridEnclosure he(DiscreteLocation(), RealSpace({x, y}), encl);
        he.set_auxiliary({z}, auxiliary);
        ARIADNE_TESTING_EXPECT_EQUAL(he.dimension(), 3);
        ARIADNE_TESTING_EXPECT_NO_THROW(he.function(x));
        ARIADNE_TESTING_EXPECT_NO_THROW(he.function(y));
        ARIADNE_TESTING_EXPECT_NO_THROW(he.function(z));
        ARIADNE_TESTING_EXPECT_NO_THROW(he.function(TimeVariable()));
        ARIADNE_TESTING_EXPECT_EXCEPTION(he.function(RealVariable("u")));
        ARIADNE_TESTING_EXPECT_EQUAL(he.auxiliary_function().result_size(), 1);
        ARIADNE_TESTING_EXPECT_EQUAL(he.state_auxiliary_function().result_size(), 3);
        ARIADNE_TESTING_EXPECT_EQUAL(he.state_time_auxiliary_function().result_size(), 4);
        ARIADNE_TESTING_EXPECT_EQUAL(he.state_set().dimension(), 2);
        ARIADNE_TESTING_EXPECT_EQUAL(he.state_time_set().dimension(), 3);
        ARIADNE_TESTING_EXPECT_EQUAL(he.state_auxiliary_set().dimension(), 3);

        ARIADNE_TESTING_PROVE(he.state_auxiliary_space() == RealSpace({x, y, z}));

        ARIADNE_TESTING_EXPECT_EQUAL(project(he, RealSpace({x})).dimension(), 1);
        ARIADNE_TESTING_EXPECT_EQUAL(project(he, RealSpace({y, z})).dimension(), 2);
        ARIADNE_TESTING_EXPECT_EXCEPTION(project(he, RealSpace({RealVariable("u")})));
    }
};

int main(int argc, char* argv[]) {
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_hybrid_enclosure");
    ARIADNE_TESTING_ADD_TEST(TestHybridEnclosure, test_time_maps_and_steps);
    ARIADNE_TESTING_ADD_TEST(TestHybridEnclosure, test_constraints_and_collections);
    ARIADNE_TESTING_ADD_TEST(TestHybridEnclosure, test_alternate_factory_formatting);
    ARIADNE_TESTING_ADD_TEST(TestHybridEnclosure, test_bounds_factory_reach_steps);
    ARIADNE_TESTING_ADD_TEST(TestHybridEnclosure, test_construct);
    ARIADNE_TESTING_ADD_TEST(TestHybridEnclosure, test_space_having_time);
    ARIADNE_TESTING_ADD_TEST(TestHybridEnclosure, test_auxiliary);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
