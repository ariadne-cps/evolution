#include "evolution/hybrid_evolver.hpp"
#include "system/hybrid_automata.hpp"
#include "integration/integrator.hpp"
#include "ariadne-testing.hpp"
#include "log_verbosity_scope.hpp"
#include "test_logging.hpp"

using namespace Ariadne;

class StepProbe : public GeneralHybridEvolver {
public:
    using GeneralHybridEvolver::GeneralHybridEvolver;
    using HybridEvolverBase::_compute_flow;
    using HybridEvolverBase::_apply_evolve_step;
    using HybridEvolverBase::_apply_guard_step;
    using HybridEvolverBase::_apply_guard;
    using GeneralHybridEvolver::_estimate_timing;
};

class TestHybridStepContracts : public Testing::ContextFixture {
public:
    using ContextFixture::ContextFixture;

    void test_rejected_step_and_guard_inputs() {
        LogVerbosityScope logging(100u);
        DiscreteLocation q(1);
        DiscreteEvent event("guard");
        RealVariable x("x"),y("y");
        RealSpace space({x,y});
        HybridAutomaton system;
        system.new_mode(q,{dot(x)=1,dot(y)=0});
        StepProbe probe(system);
        auto initial=probe.enclosure(HybridExactBox(q,space,ExactBoxType({{-1,1},{-1,1}})));
        auto dynamic=system.dynamic_function(q);
        FlowFunctionModel flow=probe._compute_flow(dynamic,cast_exact_box(initial.state_bounding_box()),Dyadic(0.5_x));
        Set<DiscreteEvent> active;
        Map<DiscreteEvent,CrossingData> crossings;
        Map<DiscreteEvent,TransitionData> transitions;
        auto timing=probe._estimate_timing(active,Real(2),flow,crossings,transitions,initial);
        auto invalid_timing=timing;
        invalid_timing.step_kind=static_cast<StepKind>(255u);
        auto target=initial;
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_evolve_step(target,flow,invalid_timing));
        ARIADNE_TESTING_EXPECT_EQUAL(target.number_of_parameters(),initial.number_of_parameters());

        TransitionData transition;
        transition.event=event;
        transition.event_kind=static_cast<EventKind>(255u);
        transition.guard_function=make_function(space,x);
        CrossingData crossing(CrossingKind::DEGENERATE);
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard_step(
            target,dynamic,flow,timing,transition,crossing,Semantics::UPPER));
        transition.event_kind=EventKind::URGENT;
        crossing.crossing_kind=static_cast<CrossingKind>(255u);
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard_step(
            target,dynamic,flow,timing,transition,crossing,Semantics::UPPER));

        auto elapsed=probe.function_factory().create_constant(initial.parameter_domain(),0.5_x);
        List<HybridEnclosure> sets{initial};
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard(
            sets,elapsed,initial,flow,transition,crossing,Semantics::UPPER));
        crossing.crossing_kind=CrossingKind::DEGENERATE;
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard(
            sets,elapsed,initial,flow,transition,crossing,static_cast<Semantics>(255u)));
        for(auto kind:{CrossingKind::NEGATIVE,CrossingKind::POSITIVE,CrossingKind::DECREASING}) {
            crossing.crossing_kind=kind;
            auto untouched=sets;
            probe._apply_guard(untouched,elapsed,initial,flow,transition,crossing,Semantics::UPPER);
            ARIADNE_TESTING_EXPECT_EQUAL(untouched[0].number_of_constraints(),initial.number_of_constraints());
        }
        auto wrong_elapsed=probe.function_factory().create_constant(ExactBoxType({{-1,1}}),0.5_x);
        ARIADNE_TESTING_EXPECT_EXCEPTION(probe._apply_guard(
            sets,wrong_elapsed,initial,flow,transition,crossing,Semantics::UPPER));
    }

    void test_grazing_whose_maximum_is_in_the_past() {
        LogVerbosityScope logging(100u);
        DiscreteLocation q(1);
        DiscreteEvent event("past maximum");
        RealVariable x("x"),y("y");
        RealSpace space({x,y});
        HybridAutomaton system;
        system.new_mode(q,{dot(x)=1,dot(y)=0});
        StepProbe probe(system);
        auto initial=probe.enclosure(HybridExactBox(q,space,
            ExactBoxType({{1.0_x,1.25_x},{0.5_x,0.75_x}})));
        auto factory=probe.function_factory();
        auto dynamic=system.dynamic_function(q);
        FlowFunctionModel flow=probe._compute_flow(dynamic,cast_exact_box(initial.state_bounding_box()),Dyadic(0.5_x));
        auto elapsed=factory.create_constant(initial.parameter_domain(),0.5_x);
        auto final_set=initial;
        final_set.apply_parameter_evolve_step(flow,elapsed);
        List<HybridEnclosure> sets{final_set};
        TransitionData transition;
        transition.event=event;
        transition.event_kind=EventKind::INVARIANT;
        transition.guard_function=make_function(space,y-sqr(x));
        CrossingData crossing(CrossingKind::GRAZING);
        crossing.critical_time=-factory.create_coordinate(flow.space_domain(),0u);
        crossing.guard_range_at_critical_time=factory.create_coordinate(flow.space_domain(),1u).range();
        probe._apply_guard(sets,elapsed,initial,flow,transition,crossing,Semantics::UPPER);
        // x only increases from x>=1, so y-x^2 stays strictly negative.
        ARIADNE_TESTING_EXPECT_EQUAL(sets.size(),1u);
        ARIADNE_TESTING_PROVE(!definitely(sets[0].is_empty()));
        ARIADNE_TESTING_PROVE(definitely(sets[0].inside(HybridExactBox(q,space,
            ExactBoxType({{1.49_x,1.76_x},{0.49_x,0.76_x}})))));
    }
};

int main(int argc,char** argv) {
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_hybrid_step_contracts");
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts,test_rejected_step_and_guard_inputs);
    ARIADNE_TESTING_ADD_TEST(TestHybridStepContracts,test_grazing_whose_maximum_is_in_the_past);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc,argv);
}
