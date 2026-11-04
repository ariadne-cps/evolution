/***************************************************************************
 *            test_hybrid_enclosure.cpp
 *
 *  Copyright  2008-21  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>

#include "macro/testing.hpp"

#include "algebra/algebra.hpp"
#include "function/taylor_function.hpp"
#include "hybrid/hybrid_set.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "hybrid/hybrid_grid.hpp"
#include "hybrid/hybrid_space.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "labelled/figure.hpp"
#include "recording_canvas.hpp"
#include "drawing/drawer.hpp"
#include "drawing/drawer_manager.hpp"
#include <memory>
#include <sstream>

using namespace std;
using namespace Ariadne;

class TestHybridEnclosure {
  public:
    Void test() {
        ARIADNE_TEST_CALL(test_time_maps_and_steps())
        ARIADNE_TEST_CALL(test_constraints_and_collections())
        ARIADNE_TEST_CALL(test_alternate_factory_formatting())
        ARIADNE_TEST_CALL(test_construct())
        ARIADNE_TEST_CALL(test_space_having_time())
        ARIADNE_TEST_CALL(test_auxiliary())
    }

    void test_time_maps_and_steps() {
        RealVariable x("x"),y("y");
        DiscreteLocation q(1);
        ExactBoxType domain({{-1,1},{-1,1}});
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8)));
        auto const& factory=config.function_factory();
        HybridEnclosure initial(HybridExactBox(q,RealSpace({x,y}),domain),config);
        auto coordinates=EffectiveVectorMultivariateFunction::identity(3u);
        EffectiveVectorMultivariateFunction flow={coordinates[0]+coordinates[2],coordinates[1]+2*coordinates[2]};
        ExactBoxType flow_domain=product(domain,ExactIntervalType(0.0_x,0.5_x));
        auto phi=factory.create(flow_domain,flow);
        auto constant_space=EffectiveScalarMultivariateFunction::constant(2u,Real(1/4_q));
        auto constant_spacetime=EffectiveScalarMultivariateFunction::constant(3u,Real(1/4_q));
        auto space=initial;
        space.apply_space_evolve_step(phi,factory.create(domain,constant_space));
        ARIADNE_TEST_ASSERT(definitely(space.state_bounding_box().inside(ExactBoxType({{-0.76_x,1.26_x},{-0.51_x,1.51_x}}))));
        ARIADNE_TEST_EQUALS(space.time_range().upper_bound().raw().get_d(),0.25);
        auto spacetime=initial;
        spacetime.apply_spacetime_evolve_step(phi,factory.create(flow_domain,constant_spacetime));
        ARIADNE_TEST_EQUALS(spacetime.time_range().upper_bound().raw().get_d(),0.25);
        auto reach=initial;
        reach.apply_spacetime_reach_step(phi,factory.create(flow_domain,constant_spacetime));
        ARIADNE_TEST_EQUALS(reach.number_of_parameters(),3u);
        ARIADNE_TEST_ASSERT(reach.time_range().upper_bound().raw().get_d()>=0.25);
        reach.reduce();
        ARIADNE_TEST_ASSERT(!definitely(reach.is_empty()));
        auto changed=space;
        changed.set_time_function(factory.create_constant(domain,Real(2)));
        ARIADNE_TEST_EQUALS(changed.time_range().lower_bound().raw().get_d(),2.0);
        ARIADNE_TEST_EQUALS(changed.time_range().upper_bound().raw().get_d(),2.0);
        ARIADNE_TEST_EQUALS(changed.dwell_time_range().upper_bound().raw().get_d(),0.25);
        ARIADNE_TEST_FAIL(changed.set_time_function(factory.create_zero(ExactBoxType({{0,1}}))));
        ARIADNE_TEST_FAIL(changed.set_time_function(factory.create_zero(ExactBoxType({{0,1},{0,1}}))));
        changed.clear_time();
        ARIADNE_TEST_EQUALS(changed.time_range().upper_bound().raw().get_d(),0.0);
        ARIADNE_TEST_EQUALS(changed.dwell_time_range().upper_bound().raw().get_d(),0.0);
        RealVariable t("t");
        HybridEnclosure clock(HybridExactBox(q,RealSpace({t,x}),domain),config);
        ARIADNE_TEST_EQUALS(clock.state_time_set().dimension(),2u);
        ARIADNE_TEST_EQUALS(clock.state_time_auxiliary_function().result_size(),2u);
        RecordingCanvas canvas; clock.draw(canvas,{},Variables2d(t,x));
        ARIADNE_TEST_ASSERT(canvas.boundaries>0u);
    }

    void test_constraints_and_collections() {
        RealVariable x("x"),y("y"),z("z");
        RealSpace space({x,y});
        DiscreteLocation q(1),other(2);
        DiscreteEvent event("event");
        ExactBoxType domain({{-1,1},{-1,1}});
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8)));
        HybridEnclosure initial(HybridBoxSet(q,{-1<=x<=1,-1<=y<=1}),space,config);
        HybridEnclosure exact(HybridExactBox(q,space,domain),config);
        std::unique_ptr<HybridEnclosure> cloned(exact.clone());
        ARIADNE_TEST_EQUALS(cloned->parameter_domain(),domain);
        auto coordinate=EffectiveVectorMultivariateFunction::identity(2u);
        auto range=initial.range_of(coordinate[0]+2*coordinate[1]);
        ARIADNE_TEST_ASSERT(range.lower_bound().raw().get_d()<=-3.0);
        ARIADNE_TEST_ASSERT(range.upper_bound().raw().get_d()>=3.0);
        auto bounded=initial;
        bounded.new_state_time_bound(event,EffectiveScalarMultivariateFunction::constant(2u,Real(1)));
        bounded.reduce();
        ARIADNE_TEST_ASSERT(!definitely(bounded.is_empty()));
        auto impossible=initial;
        impossible.new_state_time_bound(event,EffectiveScalarMultivariateFunction::constant(2u,Real(-1)));
        impossible.reduce();
        ARIADNE_TEST_ASSERT(definitely(impossible.is_empty()));
        auto constrained=initial;
        constrained.new_guard(event,coordinate[0]);
        constrained.new_state_constraint(event,ValidatedConstraint(-1,coordinate[1],1));
        auto spacetime=EffectiveVectorMultivariateFunction::identity(3u);
        constrained.new_state_time_constraint(event,ValidatedConstraint(-2,spacetime[0]+spacetime[2],2));
        ARIADNE_TEST_EQUALS(constrained.number_of_constraints(),3u);
        constrained.bound_time(ValidatedScalarMultivariateFunction(EffectiveScalarMultivariateFunction::constant(2u,Real(1))));
        constrained.set_time(ValidatedScalarMultivariateFunction(EffectiveScalarMultivariateFunction::zero(2u)));
        constrained.reduce();
        ARIADNE_TEST_ASSERT(!definitely(constrained.is_empty()));
        ARIADNE_TEST_FAIL(initial.set_auxiliary({z},EffectiveVectorMultivariateFunction::zeros(1u,1u)));
        ARIADNE_TEST_FAIL(initial.set_auxiliary({},EffectiveVectorMultivariateFunction::zeros(1u,2u)));
        ARIADNE_TEST_FAIL(initial.apply_reset(event,other,space,EffectiveVectorMultivariateFunction::identity(1u),RealSpace(),EffectiveVectorMultivariateFunction::zeros(0u,2u)));
        auto reset=initial;
        reset.apply_reset(event,other,space,coordinate,RealSpace(),EffectiveVectorMultivariateFunction::zeros(0u,2u));
        ARIADNE_TEST_EQUALS(reset.previous_events().size(),1u);
        reset.clear_events();
        ARIADNE_TEST_EQUALS(reset.previous_events().size(),0u);
        HybridExactBox enclosing(q,space,ExactBoxType({{-2,2},{-2,2}}));
        HybridExactBox disjoint(q,space,ExactBoxType({{3,4},{3,4}}));
        HybridExactBox elsewhere(other,space,domain);
        ARIADNE_TEST_ASSERT(definitely(initial.inside(enclosing)));
        ARIADNE_TEST_ASSERT(!definitely(initial.inside(elsewhere)));
        ARIADNE_TEST_ASSERT(definitely(separated(initial,disjoint)));
        ARIADNE_TEST_ASSERT(definitely(separated(initial,elsewhere)));
        ARIADNE_TEST_ASSERT(!definitely(separated(initial,enclosing)));
        auto split=initial.split();
        ARIADNE_TEST_EQUALS(split.size(),2u);
        for(auto& part : split) {
            ARIADNE_TEST_ASSERT(definitely(part.inside(enclosing)));
            part.uniform_error_recondition(); part.kuhn_recondition(); part.recondition(); part._check();
        }
        ListSet<HybridEnclosure> sets(initial),empty;
        sets.append(reset);
        sets.adjoin(ListSet<HybridEnclosure>(List<HybridEnclosure>{constrained}));
        sets.adjoin(empty);
        ARIADNE_TEST_EQUALS(sets.size(),3u);
        ARIADNE_TEST_EQUALS(sets.front().location(),q);
        ARIADNE_TEST_EQUALS(sets.back().location(),q);
        auto const& const_sets=sets;
        ARIADNE_TEST_EQUALS(const_sets.front().location(),q);
        ARIADNE_TEST_EQUALS(const_sets.back().location(),q);
        ARIADNE_TEST_EQUALS(sets[q].size(),2u);
        ARIADNE_TEST_EQUALS(sets[other].size(),1u);
        ARIADNE_TEST_EQUALS(sets[DiscreteLocation(9)].size(),0u);
        List<HybridEnclosure>& entries=sets;
        List<HybridEnclosure> const& const_entries=const_sets;
        ARIADNE_TEST_EQUALS(entries.size(),const_entries.size());
        SizeType count=0u;
        for(auto iter=sets.begin();iter!=sets.end();++iter) { ++count; }
        ARIADNE_TEST_EQUALS(count,3u);
        sets.reduce(); empty.reduce();
        ARIADNE_TEST_EQUALS(sets.bounding_box().euclidean_set(space).dimension(),2u);
        ARIADNE_TEST_FAIL(empty.bounding_box());
        struct DrawerScope {
            Drawer previous{DrawerManager::instance().drawer()};
            DrawerScope() { DrawerManager::instance().set_drawer(BoxDrawer()); }
            ~DrawerScope() { DrawerManager::instance().set_drawer(previous); }
        } drawer_scope;
        RecordingCanvas canvas;
        sets.draw(canvas,{DiscreteLocation(9)},Variables2d(x,y));
        ARIADNE_TEST_EQUALS(canvas.boundaries,0u);
        sets.draw(canvas,{q},Variables2d(x,y));
        ARIADNE_TEST_ASSERT(canvas.boundaries>0u);
        sets.draw(canvas,{},Variables2d(x,y));
        std::ostringstream brief,full;
        initial.print(brief); constrained.repr(full);
        ARIADNE_TEST_ASSERT(brief.str().find("time_range")!=std::string::npos);
        ARIADNE_TEST_ASSERT(full.str().find("HybridEnclosure")!=std::string::npos);
        MonolithicHybridSpace hybrid_space;
        hybrid_space.new_location(q,space);
        HybridGridTreePaving paving{HybridGrid(HybridSpace(hybrid_space))};
        initial.adjoin_outer_approximation_to(paving,0u);
        ARIADNE_TEST_ASSERT(paving[q].size()>0u);
        MonolithicHybridSpace incompatible;
        incompatible.new_location(q,RealSpace({y,x}));
        HybridGridTreePaving wrong{HybridGrid(HybridSpace(incompatible))};
        ARIADNE_TEST_FAIL(initial.adjoin_outer_approximation_to(wrong,0u));
        sets.clear();
        ARIADNE_TEST_EQUALS(sets.size(),0u);
    }

    void test_alternate_factory_formatting() {
        RealVariable x("x"),y("y");
        RealSpace space({x,y});
        DiscreteLocation q(1);
        EnclosureConfiguration config(ScaledFunctionPatchFactory<ValidatedBoundsTaylorModelDP>(ThresholdSweeper<FloatDP>(dp,1e-8)));
        HybridEnclosure enclosure(HybridExactBox(q,space,ExactBoxType({{-1,1},{-1,1}})),config);
        auto coordinate=EffectiveVectorMultivariateFunction::identity(2u);
        enclosure.new_state_constraint(DiscreteEvent("bound"),ValidatedConstraint(-1,coordinate[0],1));
        enclosure._check();
        std::ostringstream representation; enclosure.repr(representation);
        ARIADNE_TEST_ASSERT(representation.str().find("HybridEnclosure")!=std::string::npos);
        ARIADNE_TEST_ASSERT(representation.str().find("<=")!=std::string::npos);
        ARIADNE_TEST_EQUALS(enclosure.bounding_box().location(),q);
        ARIADNE_TEST_ASSERT(definitely(inside(enclosure,HybridRealBox(q,{-2<=x<=2,-2<=y<=2}))));
        HybridExactBox exact(q,space,ExactBoxType({{-2,2},{-2,2}}));
        ARIADNE_TEST_ASSERT(definitely(inside(enclosure,exact)));
        ARIADNE_TEST_ASSERT(definitely(inside(ListSet<HybridEnclosure>(enclosure),exact)));
        ARIADNE_TEST_ASSERT(definitely(inside(ListSet<HybridEnclosure>(),exact)));
    }

    void test_construct() {
        HybridEnclosure encl;

        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
        EnclosureConfiguration config(function_factory);
        HybridEnclosure encl2(config);
        ARIADNE_TEST_PRINT(encl2.configuration())

        DiscreteLocation location;
        RealVariable x("x"), y("y");
        HybridRealBox box(location,{1<=x<=2,0<=y<=1});
        HybridEnclosure encl3(box,config);
        ARIADNE_TEST_EQUALS(encl3.number_of_parameters(),2)
        ARIADNE_TEST_EQUALS(encl3.number_of_constraints(),0)

        Enclosure cencl(config);
        HybridEnclosure(location,RealSpace({x,y}),cencl);
    }

    void test_space_having_time() {
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-8));
        EnclosureConfiguration config(function_factory);
        RealVariable x("x"), t("t");
        HybridRealBox box(DiscreteLocation(),{1<=t<=2,0<=x<=1});
        HybridEnclosure encl(box,config);
        ARIADNE_TEST_ASSERT(encl.state_time_auxiliary_space() == RealSpace({t,x}))
    }

    void test_auxiliary() {
        ExactBoxType dom({{0.0_x,2.0_x},{1.0_x,3.0_x}});
        Enclosure encl(dom,EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8))));
        ARIADNE_TEST_PRINT(encl.auxiliary_function());
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2,0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2,1);
        EffectiveVectorMultivariateFunction auxiliary(1,x0+sqr(x1));
        RealVariable x("x"), y("y"), z("z");
        HybridEnclosure he(DiscreteLocation(),RealSpace({x,y}),encl);
        he.set_auxiliary({z},auxiliary);
        ARIADNE_TEST_EQUALS(he.dimension(),3)
        ARIADNE_TEST_EXECUTE(he.function(x))
        ARIADNE_TEST_EXECUTE(he.function(y))
        ARIADNE_TEST_EXECUTE(he.function(z))
        ARIADNE_TEST_EXECUTE(he.function(TimeVariable()))
        ARIADNE_TEST_FAIL(he.function(RealVariable("u")))
        ARIADNE_TEST_EQUALS(he.auxiliary_function().result_size(),1)
        ARIADNE_TEST_EQUALS(he.state_auxiliary_function().result_size(),3)
        ARIADNE_TEST_EQUALS(he.state_time_auxiliary_function().result_size(),4)
        ARIADNE_TEST_EQUALS(he.state_set().dimension(),2)
        ARIADNE_TEST_EQUALS(he.state_time_set().dimension(),3)
        ARIADNE_TEST_EQUALS(he.state_auxiliary_set().dimension(),3)

        ARIADNE_TEST_ASSERT(he.state_auxiliary_space() == RealSpace({x,y,z}))

        ARIADNE_TEST_EQUALS(project(he,RealSpace({x})).dimension(),1)
        ARIADNE_TEST_EQUALS(project(he,RealSpace({y,z})).dimension(),2)
        ARIADNE_TEST_FAIL(project(he,RealSpace({RealVariable("u")})))
    }
};

Int main() {
    TestHybridEnclosure().test();
    return ARIADNE_TEST_FAILURES;
}

