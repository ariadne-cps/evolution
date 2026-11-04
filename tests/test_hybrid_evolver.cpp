#include "evolution/hybrid_evolver.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "system/hybrid_automata.hpp"
#include "function/function.hpp"
#include "integration/integrator.hpp"
#include "macro/testing.hpp"
#include "function/taylor_function.decl.hpp"
#include "function/taylor_model.hpp"
#include "equation-solving/interval_newton_method.hpp"
#include "log_verbosity_scope.hpp"
#include <sstream>
#include <type_traits>
#include <cmath>
#include <memory>
#include "test_logging.hpp"
using namespace Ariadne;

// Exercise the default timing strategy provided for derived evolvers.
class FixedTimingEvolver : public GeneralHybridEvolver {
public:
    using GeneralHybridEvolver::GeneralHybridEvolver;
    FixedTimingEvolver* clone() const override { return new FixedTimingEvolver(*this); }
protected:
    TimingData _estimate_timing(Set<DiscreteEvent>& active,Real time,FlowFunctionModel const& flow,
        Map<DiscreteEvent,CrossingData>& crossings,Map<DiscreteEvent,TransitionData> const& transitions,
        HybridEnclosure const& initial) const override {
        return HybridEvolverBase::_estimate_timing(active,time,flow,crossings,transitions,initial);
    }
};

// Expose the documented crossing operation for a correlated near-miss test.
class CrossingProbeEvolver : public GeneralHybridEvolver {
public:
    using GeneralHybridEvolver::GeneralHybridEvolver;
    using HybridEvolverBase::_compute_crossings;
    using HybridEvolverBase::_compute_flow;
};

class CountingEvolver : public GeneralHybridEvolver {
public:
    using GeneralHybridEvolver::GeneralHybridEvolver;
    mutable SizeType reconditions=0u;
    CountingEvolver* clone() const override { return new CountingEvolver(*this); }
protected:
    Void _recondition(HybridEnclosure& enclosure) const override {
        ++reconditions;
        HybridEvolverBase::_recondition(enclosure);
    }
};

// Valid timing policies with genuinely varying step lengths. Their common
// linear flow has an exact solution, independently of the policy chosen.
class DependentStepEvolver : public GeneralHybridEvolver {
    StepKind _kind;
public:
    DependentStepEvolver(HybridAutomatonInterface const& system,StepKind kind)
        : GeneralHybridEvolver(system),_kind(kind) { }
    DependentStepEvolver* clone() const override { return new DependentStepEvolver(*this); }
protected:
    TimingData _estimate_timing(Set<DiscreteEvent>& active,Real time,FlowFunctionModel const& flow,
        Map<DiscreteEvent,CrossingData>& crossings,Map<DiscreteEvent,TransitionData> const& transitions,
        HybridEnclosure const& initial) const override {
        auto result=GeneralHybridEvolver::_estimate_timing(active,time,flow,crossings,transitions,initial);
        if(result.finishing_kind==FinishingKind::BEFORE_FINAL_TIME && result.step_kind==StepKind::CONSTANT_EVOLUTION_TIME) {
            auto variable=function_factory().create_coordinate(initial.parameter_domain(),0);
            if(_kind==StepKind::SPACE_DEPENDENT_EVOLUTION_TIME) variable=initial.state_function()[0];
            if(_kind==StepKind::TIME_DEPENDENT_EVOLUTION_TIME) variable=initial.time_function();
            result.step_kind=_kind;
            result.parameter_dependent_evolution_time=hlf(result.step_size)+(hlf(hlf(result.step_size)))*variable;
            result.parameter_dependent_finishing_time=initial.time_function()+result.parameter_dependent_evolution_time;
            if(_kind!=StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME) {
                auto domain=result.spacetime_dependent_evolution_time.domain();
                auto state_time=function_factory().create_coordinate(domain,
                    _kind==StepKind::SPACE_DEPENDENT_EVOLUTION_TIME ? 0u : initial.state_dimension());
                result.spacetime_dependent_evolution_time=hlf(result.step_size)+(hlf(hlf(result.step_size)))*state_time;
                result.spacetime_dependent_finishing_time=result.spacetime_dependent_evolution_time
                    +function_factory().create_coordinate(domain,initial.state_dimension());
            }
        }
        return result;
    }
};

class OvershootingIntegrator : public GradedTaylorSeriesIntegrator {
public:
    OvershootingIntegrator() : GradedTaylorSeriesIntegrator(1e-8) { }
    OvershootingIntegrator* clone() const override { return new OvershootingIntegrator(*this); }
    using GradedTaylorSeriesIntegrator::flow_step;
    FlowStepModelType flow_step(ValidatedVectorMultivariateFunction const& field,
        ExactBoxType const& domain,Suggestion<StepSizeType> const& hint) const override {
        // Suggestions need not be upper bounds in IntegratorInterface.
        return GradedTaylorSeriesIntegrator::flow_step(field,domain,suggest(2*static_cast<StepSizeType const&>(hint)));
    }
};

// An unexpected numerical backend failure must retain its original diagnostic.
class RuntimeFailureSolver : public IntervalNewtonSolver {
public:
    RuntimeFailureSolver() : IntervalNewtonSolver(1e-8,16u) { }
    RuntimeFailureSolver* clone() const override { return new RuntimeFailureSolver(*this); }
    using IntervalNewtonSolver::implicit;
    ValidatedScalarMultivariateFunctionPatch implicit(
        ValidatedScalarMultivariateFunction const&,ExactBoxType const&,ExactIntervalType const&) const override {
        throw std::runtime_error("implicit backend unavailable");
    }
};

class ProgressAutomaton : public HybridAutomaton {
public:
    ProgressAutomaton* clone() const override { return new ProgressAutomaton(*this); }
    EventKind event_kind(DiscreteLocation location,DiscreteEvent event) const override {
        auto kind=HybridAutomaton::event_kind(location,event);
        return kind==EventKind::INVARIANT ? EventKind::PROGRESS : kind;
    }
};

class TestHybridEvolver {
    inline static RealVariable x{"x"},y{"y"};
    inline static DiscreteLocation q{1};
    inline static DiscreteEvent e{"e"};
    inline static Real one{1};
    const Dyadic eighth{0.125_x},quarter{0.25_x};
    mutable std::unique_ptr<GeneralHybridEvolver> evolver_ptr;
    Void _set_evolver(HybridAutomatonInterface const& system) const {
        evolver_ptr=std::make_unique<GeneralHybridEvolver>(system);
        evolver_ptr->set_integrator(GradedTaylorSeriesIntegrator(1e-3));
    }
    Void expect_time(HybridEnclosure const& enclosure,double expected) const {
        auto range=enclosure.time_range();
        ARIADNE_TEST_ASSERT(std::abs(range.lower_bound().raw().get_d()-expected)<1e-8);
        ARIADNE_TEST_ASSERT(std::abs(range.upper_bound().raw().get_d()-expected)<1e-8);
    }
    Void expect_final(Orbit<HybridEnclosure> const& orbit,ExactBoxType const& expected,double time,SizeType events) const {
        ARIADNE_TEST_ASSERT(orbit.final().size()>0u);
        for(auto const& enclosure : orbit.final()) {
            ARIADNE_TEST_EQUALS(enclosure.location(),q);
            ARIADNE_TEST_EQUALS(enclosure.previous_events().size(),events);
            ARIADNE_TEST_ASSERT(definitely(enclosure.inside(HybridExactBox(q,RealSpace({x,y}),expected))));
            expect_time(enclosure,time);
        }
    }
public:
    Void test_all() const {
        ARIADNE_TEST_CALL(test_components_and_copy_ownership());
        ARIADNE_TEST_CALL(test_default_timing_and_policies());
        ARIADNE_TEST_CALL(test_impact_and_optional_transition());
        ARIADNE_TEST_CALL(test_grazing_and_terminating_event());
        ARIADNE_TEST_CALL(test_invariants_and_curvature());
        ARIADNE_TEST_CALL(test_correlated_near_miss());
        ARIADNE_TEST_CALL(test_uncertain_initial_time());
        ARIADNE_TEST_CALL(test_initial_constraints_and_contracts());
        ARIADNE_TEST_CALL(test_solver_budget_fallback());
        ARIADNE_TEST_CALL(test_integrator_step_limit());
        ARIADNE_TEST_CALL(test_dependent_step_policies());
        ARIADNE_TEST_CALL(test_nonlinear_grazing());
        ARIADNE_TEST_CALL(test_partial_crossing_without_creep());
        ARIADNE_TEST_CALL(test_straddling_initial_times());
        ARIADNE_TEST_CALL(test_formatting());
        ARIADNE_TEST_CALL(test_flow());
        ARIADNE_TEST_CALL(test_exact_final_time());
        ARIADNE_TEST_CALL(test_maximum_steps());
        ARIADNE_TEST_CALL(test_urgent_event());
        ARIADNE_TEST_CALL(test_step_size_event());
        ARIADNE_TEST_CALL(test_partial_event());
        ARIADNE_TEST_CALL(test_initially_active_event());
        ARIADNE_TEST_CALL(test_simultaneous_events());
        ARIADNE_TEST_CALL(test_transverse_linear_crossing());
    }
    Void test_components_and_copy_ownership() const;
    Void test_default_timing_and_policies() const;
    Void test_impact_and_optional_transition() const;
    Void test_grazing_and_terminating_event() const;
    Void test_invariants_and_curvature() const;
    Void test_correlated_near_miss() const;
    Void test_uncertain_initial_time() const;
    Void test_initial_constraints_and_contracts() const;
    Void test_solver_budget_fallback() const;
    Void test_integrator_step_limit() const;
    Void test_dependent_step_policies() const;
    Void test_nonlinear_grazing() const;
    Void test_partial_crossing_without_creep() const;
    Void test_straddling_initial_times() const;
    Void test_formatting() const;
    Void test_flow() const;
    Void test_exact_final_time() const;
    Void test_maximum_steps() const;
    Void test_urgent_event() const;
    Void test_step_size_event() const;
    Void test_partial_event() const;
    Void test_initially_active_event() const;
    Void test_simultaneous_events() const;
    Void test_transverse_linear_crossing() const;
};

Void TestHybridEvolver::test_components_and_copy_ownership() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    EnclosureConfiguration enclosure_config{ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(ThresholdSweeper<FloatDP>(dp,1e-8)))};
    auto factory=enclosure_config.function_factory();
    using FactoryInterface=EvolverFactoryInterface<HybridAutomatonInterface,HybridEnclosure,HybridTerminationCriterion>;
    static_assert(std::has_virtual_destructor_v<FactoryInterface>);
    GeneralHybridEvolverFactory standard;
    std::unique_ptr<GeneralHybridEvolver> made(standard.create(system));
    std::unique_ptr<FactoryInterface> polymorphic_factory(standard.clone());
    std::unique_ptr<EvolverInterface<HybridAutomatonInterface,HybridEnclosure,HybridTerminationCriterion>> owned(polymorphic_factory->create(system));
    ARIADNE_TEST_EQUALS(owned->system().continuous_state_space(q),RealSpace({x,y}));
    GeneralHybridEvolverFactory customized(factory);
    std::unique_ptr<GeneralHybridEvolverFactory> copied_factory(customized.clone());
    std::unique_ptr<GeneralHybridEvolver> created(copied_factory->create(system));
    created->set_function_factory(factory);
    created->set_solver(IntervalNewtonSolver(1e-9,16u));
    created->configuration().set_enable_subdivisions(false);
    ARIADNE_TEST_ASSERT(!created->configuration().enable_subdivisions());
    HybridExactBox box(q,RealSpace({x,y}),ExactBoxType({{-ExactDouble(0.01),ExactDouble(0.01)},{-ExactDouble(0.01),ExactDouble(0.01)}}));
    ARIADNE_TEST_EQUALS(created->enclosure(box).state_dimension(),2u);
    HybridBoundedConstraintSet bounded(q,{-0.01_dec<=x<=0.01_dec,-0.01_dec<=y<=0.01_dec});
    ARIADNE_TEST_EQUALS(created->enclosure(bounded).state_dimension(),2u);
    auto computed=created->orbit(box,HybridTime(0.25_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(computed.final().size(),1u);
    std::ostringstream description; description<<*created;
    ARIADNE_TEST_ASSERT(description.str().find("GeneralHybridEvolver")!=std::string::npos);
    ARIADNE_TEST_FAIL(created->reconfigure(HybridExactBoxes(),HybridExactFloatVector()));
    std::unique_ptr<GeneralHybridEvolver> surviving;
    {
        GeneralHybridEvolver original(system,factory);
        original.configuration().set_flow_accuracy(1e-7);
        original.configuration().set_maximum_step_size(0.25);
        original.ALLOW_UNWIND=true;
        surviving.reset(original.clone());
        surviving->configuration().set_maximum_step_size(0.125);
        ARIADNE_TEST_EQUALS(original.configuration().maximum_step_size(),0.25_x);
        ARIADNE_TEST_ASSERT(surviving->ALLOW_UNWIND);
        *made=original;
        made->configuration().set_flow_accuracy(1e-9);
        ARIADNE_TEST_EQUALS(original.configuration().flow_accuracy(),cast_exact(ApproximateDouble(1e-7)));
        auto& alias=*made; *made=alias;
        ARIADNE_TEST_EQUALS(made->configuration().flow_accuracy(),cast_exact(ApproximateDouble(1e-9)));
    }
    surviving->configuration().set_flow_accuracy(1e-6);
    auto copied_orbit=surviving->orbit(box,HybridTime(0.25_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(copied_orbit.final().size(),1u);
    expect_time(copied_orbit.final()[0],0.25);
    ARIADNE_TEST_ASSERT(definitely(copied_orbit.final()[0].inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{ExactDouble(0.23),ExactDouble(0.27)},{-ExactDouble(0.02),ExactDouble(0.02)}})))));
    auto coordinates=EffectiveVectorMultivariateFunction::identity(3u);
    FlowFunctionModel flow(factory.create(ExactBoxType({{-1,1},{-1,1},{0,1}}),EffectiveVectorMultivariateFunction{coordinates[0]+coordinates[2],coordinates[1]}));
    ARIADNE_TEST_EQUALS(flow.codomain().dimension(),2u);
    ARIADNE_TEST_EQUALS(flow.space_domain(),ExactBoxType({{-1,1},{-1,1}}));
    ARIADNE_TEST_EQUALS(flow.step_size(),Dyadic(1));
}

Void TestHybridEvolver::test_default_timing_and_policies() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    FixedTimingEvolver fixed(system);
    fixed.configuration().set_maximum_step_size(0.25);
    HybridExactBox box(q,RealSpace({x,y}),ExactBoxType({{-0.25_x,0.25_x},{-ExactDouble(0.01),ExactDouble(0.01)}}));
    auto fixed_orbit=fixed.orbit(box,HybridTime(0.25_x,3),Semantics::UPPER);
    expect_final(fixed_orbit,ExactBoxType({{-ExactDouble(0.01),ExactDouble(0.51)},{-ExactDouble(0.02),ExactDouble(0.02)}}),0.25,0u);
    GeneralHybridEvolver evolver(system);
    evolver.configuration().set_maximum_enclosure_radius(0.2);
    auto upper=evolver.orbit(box,HybridTime(0.125_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(upper.final().size(),2u);
    auto lower=evolver.orbit(box,HybridTime(0.125_x,3),Semantics::LOWER);
    ARIADNE_TEST_EQUALS(lower.final().size(),0u);
    evolver.configuration().set_enable_subdivisions(false);
    auto unsplit=evolver.orbit(box,HybridTime(0.125_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(unsplit.final().size(),1u);
    evolver.configuration().set_maximum_enclosure_radius(1.0);
    auto lower_flow=evolver.orbit(box,HybridTime(0.125_x,3),Semantics::LOWER);
    expect_final(lower_flow,ExactBoxType({{-ExactDouble(0.126),ExactDouble(0.376)},{-ExactDouble(0.02),ExactDouble(0.02)}}),0.125,0u);
    auto stopped=evolver.orbit(box,HybridTime(1.0_x,0),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(stopped.final().size(),1u);
    ARIADNE_TEST_EQUALS(stopped.reach().size(),0u);
    expect_time(stopped.final()[0],0.0);
    CountingEvolver counted(system);
    auto initial=counted.enclosure(box);
    auto factory=counted.function_factory();
    auto state=factory.create_identity(box.euclidean_set())+Vector<FloatDPBounds>({{-ExactDouble(0.1),ExactDouble(0.1)},{-ExactDouble(0.1),ExactDouble(0.1)}},dp);
    HybridEnclosure uncertain(q,RealSpace({x,y}),Enclosure(box.euclidean_set(),state,EnclosureConfiguration(factory)));
    counted.configuration().set_enable_reconditioning(false);
    auto without=counted.orbit(uncertain,HybridTime(0.125_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(counted.reconditions,0u);
    counted.configuration().set_enable_reconditioning(true);
    auto with=counted.orbit(uncertain,HybridTime(0.125_x,3),Semantics::UPPER);
    ARIADNE_TEST_ASSERT(counted.reconditions>0u);
    ARIADNE_TEST_EQUALS(with.final().size(),without.final().size());
}

Void TestHybridEvolver::test_impact_and_optional_transition() const {
    HybridAutomaton ball;
    ball.new_mode(q,{dot(x)=y,dot(y)=-1});
    ball.new_transition(q,e,q,{next(x)=x,next(y)=-y/2},x<=0,EventKind::IMPACT);
    _set_evolver(ball);
    evolver_ptr->configuration().set_flow_accuracy(1e-8);
    evolver_ptr->configuration().set_maximum_step_size(0.25);
    auto impact=evolver_ptr->orbit(HybridBoxSet(q,{0.999_dec<=x<=1.001_dec,-0.001_dec<=y<=0.001_dec}),HybridTime(1.5_x,3),Semantics::UPPER);
    expect_final(impact,ExactBoxType({{ExactDouble(0.03),ExactDouble(0.08)},{ExactDouble(0.58),ExactDouble(0.67)}}),1.5,1u);
    DiscreteLocation target(2);
    HybridAutomaton optional;
    optional.new_mode(q,{dot(x)=1,dot(y)=0});
    optional.new_mode(target,{dot(x)=0,dot(y)=0});
    optional.new_transition(q,e,target,{next(x)=0,next(y)=y},x>=0,EventKind::PERMISSIVE);
    _set_evolver(optional);
    evolver_ptr->configuration().set_maximum_step_size(1.0);
    auto orbit=evolver_ptr->orbit(HybridBoxSet(q,{-0.2_dec<=x<=-0.1_dec,-0.01_dec<=y<=0.01_dec}),HybridTime(0.5_x,1),Semantics::UPPER);
    Bool flowed=false,jumped=false;
    for(auto const& final : orbit.final()) {
        if(final.location()==q) {
            flowed=true;
            ARIADNE_TEST_EQUALS(final.previous_events().size(),0u);
            expect_time(final,0.5);
            ARIADNE_TEST_ASSERT(definitely(final.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{ExactDouble(0.29),ExactDouble(0.41)},{-ExactDouble(0.02),ExactDouble(0.02)}})))));
        } else {
            jumped=true;
            ARIADNE_TEST_EQUALS(final.location(),target);
            ARIADNE_TEST_EQUALS(final.previous_events().size(),1u);
            ARIADNE_TEST_ASSERT(final.time_range().lower_bound().raw().get_d()>=-1e-8);
            ARIADNE_TEST_ASSERT(final.time_range().upper_bound().raw().get_d()<=0.5+1e-8);
        }
    }
    ARIADNE_TEST_ASSERT(flowed && jumped);
}

Void TestHybridEvolver::test_grazing_and_terminating_event() const {
    HybridAutomaton tangent;
    tangent.new_mode(q,{dot(x)=1,dot(y)=0});
    tangent.new_transition(q,e,q,{next(x)=x,next(y)=y-1},y-sqr(x)>=0,EventKind::URGENT);
    _set_evolver(tangent);
    evolver_ptr->configuration().set_maximum_step_size(2.0);
    auto grazing=evolver_ptr->orbit(HybridBoxSet(q,{-1-eighth<=x<=-1+eighth,-quarter<=y<=quarter}),HybridTime(2.0_x,3),Semantics::UPPER);
    Bool jumped=false,flowed=false;
    for(auto const& final : grazing.final()) {
        auto count=final.previous_events().size();
        ARIADNE_TEST_ASSERT(count<=1u);
        if(count==0u) flowed=true; else jumped=true;
        ARIADNE_TEST_ASSERT(final.time_range().lower_bound().raw().get_d()>=1.9);
        ARIADNE_TEST_ASSERT(final.time_range().upper_bound().raw().get_d()<=2.1);
        ExactBoxType expected=count==0u ? ExactBoxType({{ExactDouble(0.874),ExactDouble(1.126)},{-ExactDouble(0.26),ExactDouble(0.26)}}) : ExactBoxType({{ExactDouble(0.874),ExactDouble(1.126)},{-ExactDouble(1.26),-ExactDouble(0.74)}});
        ARIADNE_TEST_ASSERT(definitely(final.inside(HybridExactBox(q,RealSpace({x,y}),expected))));
    }
    ARIADNE_TEST_ASSERT(jumped && flowed);
    DiscreteLocation target(2);
    HybridAutomaton urgent;
    urgent.new_mode(q,{dot(x)=1,dot(y)=0});
    urgent.new_mode(target,{dot(x)=0,dot(y)=0});
    urgent.new_transition(q,e,target,{next(x)=0,next(y)=y},x>=Real(1/4_q),EventKind::URGENT);
    _set_evolver(urgent);
    HybridTerminationCriterion stop(Real(1),Integer(5),Set<DiscreteEvent>{e});
    for(auto semantics : {Semantics::UPPER,Semantics::LOWER}) {
        auto orbit=evolver_ptr->orbit(HybridBoxSet(q,{-0.01_dec<=x<=0.01_dec,-0.01_dec<=y<=0.01_dec}),stop,semantics);
        ARIADNE_TEST_ASSERT(orbit.final().size()>0u);
        for(auto const& final : orbit.final()) {
            ARIADNE_TEST_EQUALS(final.location(),target);
            ARIADNE_TEST_EQUALS(final.previous_events().size(),1u);
            ARIADNE_TEST_EQUALS(final.previous_events()[0],e);
            ARIADNE_TEST_ASSERT(final.time_range().lower_bound().raw().get_d()>0.2);
            ARIADNE_TEST_ASSERT(final.time_range().upper_bound().raw().get_d()<0.3);
        }
    }
}

Void TestHybridEvolver::test_invariants_and_curvature() const {
    for(Bool progress : {false,true}) {
        HybridAutomaton invariant;
        ProgressAutomaton progressing;
        for(HybridAutomaton* system : {&invariant,static_cast<HybridAutomaton*>(&progressing)}) {
            system->new_mode(q,{dot(x)=1,dot(y)=0});
            system->new_invariant(q,x<=0.25_dec,e);
        }
        _set_evolver(progress ? static_cast<HybridAutomatonInterface const&>(progressing) : invariant);
        evolver_ptr->configuration().set_maximum_step_size(0.25);
        HybridBoxSet initial(q,{-0.01_dec<=x<=0.01_dec,-0.01_dec<=y<=0.01_dec});
        auto before=evolver_ptr->orbit(initial,HybridTime(0.125_x,3),Semantics::UPPER);
        expect_final(before,ExactBoxType({{ExactDouble(0.114),ExactDouble(0.136)},{-ExactDouble(0.011),ExactDouble(0.011)}}),0.125,0u);
        auto blocked=evolver_ptr->orbit(initial,HybridTime(0.5_x,3),Semantics::UPPER);
        ARIADNE_TEST_EQUALS(blocked.final().size(),0u);
        ARIADNE_TEST_ASSERT(blocked.reach().size()>0u);
        for(auto const& reach : blocked.reach()) {
            ARIADNE_TEST_EQUALS(reach.previous_events().size(),0u);
        }
    }
    DiscreteLocation target(2);
    // Under x'=1,y'=0 these guards are respectively convex, concave,
    // and degenerate in time. A jump resets to a stationary origin.
    for(unsigned shape=0;shape!=3;++shape) {
        HybridAutomaton system;
        system.new_mode(q,{dot(x)=1,dot(y)=0});
        system.new_mode(target,{dot(x)=0,dot(y)=0});
        auto guard=shape==0 ? sqr(x)+y : (shape==1 ? y-sqr(x) : x*x*x+y);
        system.new_transition(q,e,target,{next(x)=0,next(y)=0},guard>=0,EventKind::URGENT);
        _set_evolver(system);
        evolver_ptr->ALLOW_CREEP=false;
        evolver_ptr->configuration().set_maximum_step_size(1.0);
        HybridBoxSet initial=shape==0 ? HybridBoxSet(q,{-0.01_dec<=x<=0.01_dec,-0.063_dec<=y<=-0.062_dec})
            : shape==1 ? HybridBoxSet(q,{-0.5_dec<=x<=-0.49_dec,0.062_dec<=y<=0.063_dec})
            : HybridBoxSet(q,{-0.01_dec<=x<=0.01_dec,-0.001_dec<=y<=-0.0005_dec});
        auto orbit=evolver_ptr->orbit(initial,HybridTime(0.75_x,1),Semantics::UPPER);
        ARIADNE_TEST_ASSERT(orbit.final().size()>0u);
        for(auto const& final : orbit.final()) {
            ARIADNE_TEST_EQUALS(final.location(),target);
            ARIADNE_TEST_EQUALS(final.previous_events().size(),1u);
            ARIADNE_TEST_ASSERT(definitely(final.inside(HybridExactBox(target,RealSpace({x,y}),ExactBoxType({{-ExactDouble(0.001),ExactDouble(0.001)},{-ExactDouble(0.001),ExactDouble(0.001)}})))));
            ARIADNE_TEST_ASSERT(final.time_range().lower_bound().raw().get_d()>=-1e-8);
            ARIADNE_TEST_ASSERT(final.time_range().upper_bound().raw().get_d()<=0.75+1e-8);
        }
        if(shape==2) {
            auto lower=evolver_ptr->orbit(initial,HybridTime(0.75_x,1),Semantics::LOWER);
            ARIADNE_TEST_EQUALS(lower.final().size(),0u);
        }
    }
}

Void TestHybridEvolver::test_correlated_near_miss() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=1});
    DiscreteEvent second("second");
    system.new_transition(q,e,q,{next(x)=x,next(y)=y},y-sqr(x)>=0,EventKind::PERMISSIVE);
    system.new_transition(q,second,q,{next(x)=x,next(y)=y},y-sqr(x)+0.0001_dec>=0,EventKind::PERMISSIVE);
    CrossingProbeEvolver probe(system);
    auto factory=probe.function_factory();
    ExactBoxType parameters({{-0.5_x,-ExactDouble(0.49)},{-ExactDouble(0.0001),ExactDouble(0.0001)}});
    auto coordinates=factory.create_identity(parameters);
    auto mapping=coordinates;
    mapping[1]=coordinates[0]-factory.create_constant(parameters,Real(0.251_dec))+coordinates[1];
    HybridEnclosure initial(q,RealSpace({x,y}),Enclosure(parameters,mapping,EnclosureConfiguration(factory)));
    auto dynamic=system.dynamic_function(q);
    auto flow=probe._compute_flow(dynamic,cast_exact_box(initial.state_bounding_box()),Dyadic(2));
    const Set<DiscreteEvent> active={e,second};
    Map<DiscreteEvent,EffectiveScalarMultivariateFunction> guards;
    guards.insert(e,system.guard_function(q,e));
    guards.insert(second,system.guard_function(q,second));
    auto crossings=probe._compute_crossings(active,dynamic,guards,flow,initial);
    // Along every trajectory y-x is constant. Both guard maxima at x=1/2
    // are strictly negative, despite the coarse flow box meeting the guards.
    ARIADNE_TEST_ASSERT(crossings.empty());
    ARIADNE_TEST_EQUALS(active.size(),2u);
}

Void TestHybridEvolver::test_uncertain_initial_time() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    for(double width : {0.0625,0.125,0.5}) {
        GeneralHybridEvolver evolver(system);
        evolver.ALLOW_UNWIND=true;
        evolver.configuration().set_maximum_step_size(0.25);
        ExactBoxType domain({{FloatDP(0,dp),FloatDP(ExactDouble(width),dp)},{-ExactDouble(0.01),ExactDouble(0.01)}});
        HybridEnclosure initial=evolver.enclosure(HybridExactBox(q,RealSpace({x,y}),domain));
        initial.set_time_function(evolver.function_factory().create_coordinate(domain,0));
        auto orbit=evolver.orbit(initial,HybridTime(1.25_x,3),Semantics::UPPER);
        // x(t0)=t0 gives x(T)=T for every initial time in the enclosure.
        expect_final(orbit,ExactBoxType({{ExactDouble(1.249),ExactDouble(1.251)},{-ExactDouble(0.011),ExactDouble(0.011)}}),1.25,0u);
    }
}

Void TestHybridEvolver::test_initial_constraints_and_contracts() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    system.new_invariant(q,x<=0.25_dec,e);
    GeneralHybridEvolver evolver(system);
    evolver.configuration().set_maximum_step_size(0.25);
    HybridBoxSet partial(q,{0.125_dec<=x<=0.375_dec,-0.01_dec<=y<=0.01_dec});
    auto orbit=evolver.orbit(partial,HybridTime(0.0625_x,3),Semantics::UPPER);
    expect_final(orbit,ExactBoxType({{ExactDouble(0.1874),ExactDouble(0.2501)},{-ExactDouble(0.011),ExactDouble(0.011)}}),0.0625,0u);
    auto stopped=evolver.orbit(partial,HybridTime(1.0_x,0),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(stopped.final().size(),1u);
    auto stopped_final=stopped.final()[0];
    stopped_final.reduce();
    ARIADNE_TEST_ASSERT(definitely(stopped_final.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{ExactDouble(0.1249),ExactDouble(0.2501)},{-ExactDouble(0.011),ExactDouble(0.011)}})))));
    auto outside=evolver.orbit(HybridBoxSet(q,{0.5_dec<=x<=0.75_dec,0<=y<=0.125_dec}),HybridTime(0.125_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(outside.reach().size(),0u);
    ARIADNE_TEST_EQUALS(outside.final().size(),0u);
    ExactBoxType domain({{0,0.125_x},{0,0.125_x}});
    auto config=EnclosureConfiguration(evolver.function_factory());
    DiscreteLocation unknown(404);
    HybridEnclosure unknown_location(unknown,RealSpace({x,y}),Enclosure(domain,config));
    ARIADNE_TEST_FAIL(evolver.orbit(unknown_location,HybridTime(0.25_x,3),Semantics::UPPER));
    ARIADNE_TEST_FAIL(evolver.enclosure(HybridExactBox(unknown,RealSpace({x,y}),domain)));
    ARIADNE_TEST_FAIL(evolver.orbit(HybridExactBox(unknown,RealSpace({x,y}),domain),HybridTime(0.25_x,3),Semantics::UPPER));
    RealVariable z("z");
    HybridEnclosure wrong_variables(q,RealSpace({x,z}),Enclosure(domain,config));
    ARIADNE_TEST_FAIL(evolver.orbit(wrong_variables,HybridTime(0.25_x,3),Semantics::UPPER));
}

Void TestHybridEvolver::test_solver_budget_fallback() const {
    DiscreteLocation target(2);
    for(Bool concave : {false,true}) {
        HybridAutomaton system;
        system.new_mode(q,{dot(x)=1,dot(y)=0});
        system.new_mode(target,{dot(x)=0,dot(y)=0});
        auto guard=concave ? y-sqr(x) : x-0.25_dec;
        system.new_transition(q,e,target,{next(x)=0,next(y)=0},guard>=0,EventKind::URGENT);
        GeneralHybridEvolver evolver(system);
        // The solver is allowed to exhaust its proof budget. Evolution must
        // retain a validated event enclosure using guard constraints instead.
        evolver.set_solver(IntervalNewtonSolver(1e-8,0u));
        evolver.ALLOW_CREEP=false;
        evolver.configuration().set_maximum_step_size(1.0);
        HybridBoxSet initial=concave
            ? HybridBoxSet(q,{-0.5_dec<=x<=-0.49_dec,0.062_dec<=y<=0.063_dec})
            : HybridBoxSet(q,{-0.01_dec<=x<=0.01_dec,-0.01_dec<=y<=0.01_dec});
        for(auto semantics : {Semantics::UPPER,Semantics::LOWER}) {
            auto orbit=evolver.orbit(initial,HybridTime(0.75_x,1),semantics);
                ARIADNE_TEST_ASSERT(orbit.final().size()>0u);
                for(auto const& final : orbit.final()) {
                    ARIADNE_TEST_EQUALS(final.location(),target);
                    ARIADNE_TEST_EQUALS(final.previous_events().size(),1u);
                    ARIADNE_TEST_EQUALS(final.previous_events()[0],e);
                    ARIADNE_TEST_ASSERT(definitely(final.inside(HybridExactBox(target,RealSpace({x,y}),ExactBoxType({{-ExactDouble(0.001),ExactDouble(0.001)},{-ExactDouble(0.001),ExactDouble(0.001)}})))));
                    ARIADNE_TEST_ASSERT(final.time_range().lower_bound().raw().get_d()>=-1e-8);
                    ARIADNE_TEST_ASSERT(final.time_range().upper_bound().raw().get_d()<=0.75+1e-8);
                }
        }
    }
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    system.new_transition(q,e,q,{next(x)=0,next(y)=y},x>=0.25_dec,EventKind::URGENT);
    CrossingProbeEvolver probe(system);
    probe.set_solver(RuntimeFailureSolver());
    auto initial=probe.enclosure(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{-ExactDouble(0.01),ExactDouble(0.01)},{-ExactDouble(0.01),ExactDouble(0.01)}})));
    auto dynamic=system.dynamic_function(q);
    auto flow=probe._compute_flow(dynamic,cast_exact_box(initial.state_bounding_box()),Dyadic(1));
    Map<DiscreteEvent,EffectiveScalarMultivariateFunction> guards;
    guards.insert(e,system.guard_function(q,e));
    Bool received=false;
    try { (void)probe._compute_crossings(Set<DiscreteEvent>{e},dynamic,guards,flow,initial); }
    catch(std::runtime_error const& error) {
        received=true;
        ARIADNE_TEST_EQUALS(std::string(error.what()),std::string("implicit backend unavailable"));
    }
    ARIADNE_TEST_ASSERT(received);
}

Void TestHybridEvolver::test_integrator_step_limit() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    CrossingProbeEvolver evolver(system);
    OvershootingIntegrator integrator;
    evolver.set_integrator(integrator);
    evolver.configuration().set_maximum_step_size(0.25);
    ExactBoxType domain({{-ExactDouble(0.01),ExactDouble(0.01)},{-ExactDouble(0.01),ExactDouble(0.01)}});
    auto field=system.dynamic_function(q);
    auto available=integrator.flow_step(field,domain,suggest(Dyadic(1/4_q)));
    ARIADNE_TEST_EQUALS(available.domain()[2].upper_bound(),FloatDP(0.5_x,dp));
    FlowFunctionModel capped(evolver._compute_flow(field,domain,Dyadic(1/4_q)));
    ARIADNE_TEST_EQUALS(capped.step_size(),Dyadic(1/4_q));
    auto orbit=evolver.orbit(HybridExactBox(q,RealSpace({x,y}),domain),HybridTime(0.5_x,3),Semantics::UPPER);
    expect_final(orbit,ExactBoxType({{ExactDouble(0.489),ExactDouble(0.511)},{-ExactDouble(0.011),ExactDouble(0.011)}}),0.5,0u);
    ARIADNE_TEST_EQUALS(orbit.reach().size(),2u);
}

Void TestHybridEvolver::test_dependent_step_policies() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    for(auto kind : {StepKind::SPACE_DEPENDENT_EVOLUTION_TIME,StepKind::TIME_DEPENDENT_EVOLUTION_TIME,
                    StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME}) {
        DependentStepEvolver evolver(system,kind);
        evolver.ALLOW_CREEP=false;
        evolver.configuration().set_maximum_step_size(0.125);
        auto initial=evolver.enclosure(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{-ExactDouble(0.01),ExactDouble(0.01)},{-ExactDouble(0.01),ExactDouble(0.01)}})));
        initial.set_time_function(evolver.function_factory().create_constant(initial.parameter_domain(),Real(1/8_q)));
        auto orbit=evolver.orbit(initial,HybridTime(0.5_x,3),Semantics::UPPER);
        expect_final(orbit,ExactBoxType({{ExactDouble(0.364),ExactDouble(0.386)},
            {-ExactDouble(0.011),ExactDouble(0.011)}}),0.5,0u);
        ARIADNE_TEST_ASSERT(orbit.intermediate().size()>=2u);
        // These clocks distinguish all three policies and verify two steps,
        // rather than merely accepting an unchanged final solution.
        double first=kind==StepKind::TIME_DEPENDENT_EVOLUTION_TIME ? 0.19140625 : 0.1875;
        double second=kind==StepKind::TIME_DEPENDENT_EVOLUTION_TIME ? 0.2598876953125
            : (kind==StepKind::SPACE_DEPENDENT_EVOLUTION_TIME ? 0.251953125 : 0.25);
        for(SizeType i=0;i<2u;++i) {
            auto const& intermediate=orbit.intermediate()[i];
            auto point=Vector<FloatDPApproximation>(intermediate.parameter_domain().midpoint(),dp);
            auto clock=intermediate.time_function()(point);
            ARIADNE_TEST_ASSERT(std::abs(clock.raw().get_d()-(i==0u ? first : second))<1e-10);
        }
    }
}

Void TestHybridEvolver::test_nonlinear_grazing() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    // The guard reaches its unique maximum at x=0; negative y misses it,
    // positive y crosses it, and y=0 is tangent.
    system.new_transition(q,e,q,{next(x)=x,next(y)=y-1},y-exp(x)-exp(-x)+2>=0,EventKind::URGENT);
    GeneralHybridEvolver evolver(system);
    evolver.set_solver(IntervalNewtonSolver(1e-3,16u));
    evolver.ALLOW_CREEP=false;
    evolver.configuration().set_maximum_step_size(1.0);
    auto orbit=evolver.orbit(HybridBoxSet(q,{-0.5_dec<=x<=-0.49_dec,-0.01_dec<=y<=0.01_dec}),
        HybridTime(1.0_x,2),Semantics::UPPER);
    Bool jumped=false,flowed=false;
    for(auto const& final : orbit.final()) {
        auto count=final.previous_events().size();
        ARIADNE_TEST_ASSERT(count<=1u);
        if(count==0u) flowed=true; else { jumped=true; ARIADNE_TEST_EQUALS(final.previous_events()[0],e); }
        expect_time(final,1.0);
        auto expected=count==0u
            ? ExactBoxType({{ExactDouble(0.48),ExactDouble(0.53)},{-ExactDouble(0.02),ExactDouble(0.02)}})
            : ExactBoxType({{ExactDouble(0.48),ExactDouble(0.53)},{-ExactDouble(1.02),-ExactDouble(0.98)}});
        ARIADNE_TEST_ASSERT(definitely(final.inside(HybridExactBox(q,RealSpace({x,y}),expected))));
    }
    ARIADNE_TEST_ASSERT(jumped && flowed);
}

Void TestHybridEvolver::test_partial_crossing_without_creep() const {
    HybridAutomaton system;
    DiscreteLocation target(2);
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    system.new_mode(target,{dot(x)=0,dot(y)=0});
    system.new_transition(q,e,target,{next(x)=0,next(y)=y},x>=0.4_dec,EventKind::URGENT);
    GeneralHybridEvolver evolver(system);
    evolver.ALLOW_CREEP=false;
    evolver.configuration().set_maximum_step_size(0.25);
    auto orbit=evolver.orbit(HybridBoxSet(q,{-0.25_dec<=x<=0.25_dec,-0.01_dec<=y<=0.01_dec}),
        HybridTime(1.0_x,1),Semantics::UPPER);
    ARIADNE_TEST_ASSERT(orbit.final().size()>0u);
    for(auto const& final : orbit.final()) {
        ARIADNE_TEST_EQUALS(final.location(),target);
        ARIADNE_TEST_EQUALS(final.previous_events().size(),1u);
        ARIADNE_TEST_EQUALS(final.previous_events()[0],e);
        ARIADNE_TEST_ASSERT(final.time_range().lower_bound().raw().get_d()>0.149);
        ARIADNE_TEST_ASSERT(final.time_range().upper_bound().raw().get_d()<0.651);
        ARIADNE_TEST_ASSERT(definitely(final.inside(HybridExactBox(target,RealSpace({x,y}),
            ExactBoxType({{-ExactDouble(0.001),ExactDouble(0.001)},{-ExactDouble(0.011),ExactDouble(0.011)}})))));
    }
}

Void TestHybridEvolver::test_straddling_initial_times() const {
    HybridAutomaton system;
    system.new_mode(q,{dot(x)=1,dot(y)=0});
    GeneralHybridEvolver evolver(system);
    evolver.ALLOW_CREEP=false;
    evolver.configuration().set_maximum_step_size(0.125);
    ExactBoxType domain({{0,0.5_x},{-ExactDouble(0.01),ExactDouble(0.01)}});
    auto initial=evolver.enclosure(HybridExactBox(q,RealSpace({x,y}),domain));
    initial.set_time_function(evolver.function_factory().create_coordinate(domain,0)-Dyadic(1/8_q));
    auto orbit=evolver.orbit(initial,HybridTime(0.25_x,3),Semantics::UPPER);
    // x(t0)=t0+1/8, hence x(1/4)=3/8 even when the initial
    // time range extends past the horizon for part of the parameter domain.
    ARIADNE_TEST_ASSERT(orbit.final().size()>0u);
    for(auto const& final : orbit.final()) {
        ARIADNE_TEST_EQUALS(final.previous_events().size(),0u);
        // A constrained slice can retain a non-point clock range. Check the
        // final-time equality and the conserved relation x-t=1/8 numerically.
        auto near_zero=[](UpperIntervalType const& range) {
            return range.lower_bound().raw().get_d()>-1e-10 && range.upper_bound().raw().get_d()<1e-10;
        };
        auto clock_residual=final.time_function()-Dyadic(1/4_q);
        Bool sliced=near_zero(clock_residual.range());
        for(auto const& constraint : final.continuous_set().constraint_models()) {
            if(definitely(constraint.lower_bound()==0) && definitely(constraint.upper_bound()==0)) {
                sliced=sliced || near_zero((constraint.function()-clock_residual).range())
                    || near_zero((constraint.function()+clock_residual).range());
            }
        }
        ARIADNE_TEST_ASSERT(sliced);
        ARIADNE_TEST_ASSERT(near_zero((final.state_function()[0]-final.time_function()-Dyadic(1/8_q)).range()));
        auto point=final.continuous_set().reduced_domain().midpoint();
        auto state=final.state_function()(Vector<FloatDPApproximation>(point,dp));
        auto clock=final.time_function()(Vector<FloatDPApproximation>(point,dp));
        ARIADNE_TEST_ASSERT(std::abs(state[0].raw().get_d()-0.375)<1e-10);
        ARIADNE_TEST_ASSERT(std::abs(clock.raw().get_d()-0.25)<1e-10);
    }
    initial.set_time_function(evolver.function_factory().create_constant(domain,Real(1)));
    auto past=evolver.orbit(initial,HybridTime(0.25_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(past.reach().size(),0u);
    ARIADNE_TEST_EQUALS(past.final().size(),0u);
    initial.set_time_function(evolver.function_factory().create_zero(domain));
    auto zero_horizon=evolver.orbit(initial,HybridTime(0.0_x,3),Semantics::UPPER);
    ARIADNE_TEST_EQUALS(zero_horizon.reach().size(),0u);
    ARIADNE_TEST_ASSERT(zero_horizon.final().size()>0u);
    for(auto const& final : zero_horizon.final()) expect_time(final,0.0);
}

Void TestHybridEvolver::test_formatting() const {
    auto check=[](auto value,const char* expected) {
        std::ostringstream text; text<<value;
        ARIADNE_TEST_EQUALS(text.str(),std::string(expected));
    };
    ARIADNE_TEST_ASSERT(is_blocking(EventKind::INVARIANT));
    ARIADNE_TEST_ASSERT(is_blocking(EventKind::PROGRESS));
    ARIADNE_TEST_ASSERT(!is_activating(EventKind::INVARIANT));
    ARIADNE_TEST_ASSERT(!is_activating(EventKind::PROGRESS));
    ARIADNE_TEST_FAIL(is_blocking(static_cast<EventKind>(255)));
    ARIADNE_TEST_FAIL(is_activating(static_cast<EventKind>(255)));
    check(CrossingKind::DEGENERATE,"DEGENERATE"); check(CrossingKind::POSITIVE,"POSITIVE");
    check(CrossingKind::NEGATIVE,"NEGATIVE"); check(CrossingKind::MONOTONE,"MONOTONE");
    check(CrossingKind::INCREASING,"INCREASING"); check(CrossingKind::DECREASING,"DECREASING");
    check(CrossingKind::CONVEX,"CONVEX"); check(CrossingKind::CONCAVE,"CONCAVE");
    check(CrossingKind::TRANSVERSE,"TRANSVERSE"); check(CrossingKind::GRAZING,"GRAZING");
    check(static_cast<CrossingKind>(255),"unknown");
    check(StepKind::CONSTANT_EVOLUTION_TIME,"CONSTANT_EVOLUTION_TIME");
    check(StepKind::SPACE_DEPENDENT_EVOLUTION_TIME,"SPACE_DEPENDENT_EVOLUTION_TIME");
    check(StepKind::TIME_DEPENDENT_EVOLUTION_TIME,"TIME_DEPENDENT_EVOLUTION_TIME");
    check(StepKind::SPACETIME_DEPENDENT_EVOLUTION_TIME,"SPACETIME_DEPENDENT_EVOLUTION_TIME");
    check(StepKind::PARAMETER_DEPENDENT_EVOLUTION_TIME,"PARAMETER_DEPENDENT_EVOLUTION_TIME");
    check(StepKind::PARAMETER_DEPENDENT_FINISHING_TIME,"PARAMETER_DEPENDENT_FINISHING_TIME");
    check(StepKind::SPACETIME_DEPENDENT_FINISHING_TIME,"SPACETIME_DEPENDENT_FINISHING_TIME");
    check(StepKind::CONSTANT_FINISHING_TIME,"CONSTANT_FINISHING_TIME");
    check(static_cast<StepKind>(255),"UNKNOWN");
    check(FinishingKind::BEFORE_FINAL_TIME,"BEFORE_FINAL_TIME"); check(FinishingKind::AT_FINAL_TIME,"AT_FINAL_TIME");
    check(FinishingKind::AFTER_FINAL_TIME,"AFTER_FINAL_TIME"); check(FinishingKind::STRADDLE_FINAL_TIME,"STRADDLE_FINAL_TIME");
    check(static_cast<FinishingKind>(255),"UNKNOWN");
    EnclosureConfiguration config{ValidatedFunctionPatchFactory(make_taylor_function_patch_factory(ThresholdSweeper<FloatDP>(dp,1e-8)))};
    auto time=config.function_factory().create_constant(ExactBoxType({{-1,1},{-1,1}}),Real(1/4_q));
    CrossingData transverse(CrossingKind::TRANSVERSE,time),grazing(CrossingKind::GRAZING),negative(CrossingKind::NEGATIVE);
    grazing.critical_time=time;
    for(auto const& crossing : {transverse,grazing,negative}) {
        std::ostringstream text; text<<crossing;
        ARIADNE_TEST_ASSERT(text.str().find("kind=")!=std::string::npos);
    }
}

Void TestHybridEvolver::test_flow() const {
    HybridAutomaton automaton;
    automaton.new_mode(q,{dot(x)=one,dot(y)=one/2});
    RealSpace space=automaton.continuous_state_space(q);

    HybridBoxSet initial(q,{-eighth<=x<=eighth,-eighth<=y<=eighth});
    HybridTime time(2.5_x,3);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(1.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);
    ARIADNE_TEST_EQUALS(orbit.final().size(),1u);
    ARIADNE_TEST_EQUALS(orbit.reach().size(),3u);
    ARIADNE_TEST_EQUALS(orbit.intermediate().size(),2u);


    expect_final(orbit,ExactBoxType({{ExactDouble(2.374),ExactDouble(2.626)},{ExactDouble(1.124),ExactDouble(1.376)}}),2.5,0u);
}

Void TestHybridEvolver::test_exact_final_time() const {
    HybridAutomaton automaton;
    automaton.new_mode(q,{dot(x)=one,dot(y)=one/2});
    RealSpace space=automaton.continuous_state_space(q);

    HybridBoxSet initial(q,{-eighth<=x<=eighth,-eighth<=y<=eighth});
    HybridTime time(2.0_x,1);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(1.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);

    ARIADNE_TEST_EQUALS(orbit.final().size(),1u);
    ARIADNE_TEST_EQUALS(orbit.reach().size(),2u);

    expect_final(orbit,ExactBoxType({{ExactDouble(1.874),ExactDouble(2.126)},{ExactDouble(0.874),ExactDouble(1.126)}}),2.0,0u);
}

Void TestHybridEvolver::test_maximum_steps() const {
    HybridAutomaton automaton;
    automaton.new_mode(q,{dot(x)=one,dot(y)=one/2});
    const Dyadic c1(1.5); const Dyadic c2(1.0/16); const Dyadic c3(0.5);
    automaton.new_transition(q,e,q,{next(x)=x-c1,next(y)=y},x+c2*y-c3>=0,EventKind::URGENT);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBoxSet initial(q,{-eighth<=x<=eighth,-eighth<=y<=eighth});
    HybridTime time(1.0_x,1);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(1.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);

    ARIADNE_TEST_EQUALS(orbit.final().size(),1u);
    ARIADNE_TEST_EQUALS(orbit.reach().size(),1u);

    for(auto const& enclosure : orbit.final()) {
        ARIADNE_TEST_EQUALS(enclosure.previous_events().size(),1u);
        ARIADNE_TEST_EQUALS(enclosure.previous_events()[0],e);
        ARIADNE_TEST_ASSERT(enclosure.time_range().upper_bound().raw().get_d()<1.0);
        ARIADNE_TEST_ASSERT(definitely(enclosure.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{-ExactDouble(1.1),-ExactDouble(0.9)},{-ExactDouble(0.2),ExactDouble(0.6)}})))));
    }
}

Void TestHybridEvolver::test_urgent_event() const {
    HybridAutomaton automaton;
    automaton.new_mode(q,{dot(x)=one,dot(y)=one/2});
    const Dyadic c1(1.5); const Dyadic c2(1.0/16); const Dyadic c3(0.5);
    automaton.new_transition(q,e,q,{next(x)=x-c1,next(y)=y},x+c2*y-c3>=0,EventKind::URGENT);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBoxSet initial(q,{-eighth<=x<=eighth,-eighth<=y<=eighth});
    HybridTime time(1.0_x,2);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(1.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);

    ARIADNE_TEST_EQUALS(orbit.final().size(),1u);
    ARIADNE_TEST_EQUALS(orbit.reach().size(),2u);

    expect_final(orbit,ExactBoxType({{-ExactDouble(0.626),-ExactDouble(0.374)},{ExactDouble(0.374),ExactDouble(0.626)}}),1.0,1u);
}

Void TestHybridEvolver::test_step_size_event() const {
    HybridAutomaton automaton;
    automaton.new_mode(q,{dot(x)=one,dot(y)=one/2});
    automaton.new_transition(q,e,q,{next(x)=x-2,next(y)=y},x-2>=0,EventKind::URGENT);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBoxSet initial(q,{-eighth<=x<=eighth,-eighth<=y<=eighth});
    HybridTime time(1.0_x,3);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(2.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);

    ARIADNE_TEST_EQUALS(orbit.final().size(),1u);
    ARIADNE_TEST_EQUALS(orbit.reach().size(),1u);

    expect_final(orbit,ExactBoxType({{ExactDouble(0.874),ExactDouble(1.126)},{ExactDouble(0.374),ExactDouble(0.626)}}),1.0,0u);
}

Void TestHybridEvolver::test_partial_event() const {
    HybridAutomaton automaton;
    automaton.new_mode(q,{dot(x)=one,dot(y)=one/2});
    automaton.new_transition(q,e,q,{next(x)=x-2,next(y)=y},x-y/16-2>=0,EventKind::URGENT);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBoxSet initial(q,{-eighth<=x<=quarter,-eighth<=y<=eighth});
    HybridTime time(2.0_x,3);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(2.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);

    ARIADNE_TEST_EQUALS(orbit.final().size(),2u);
    ARIADNE_TEST_EQUALS(orbit.reach().size(),2u);

    Bool jumped=false,flowed=false;
    for(auto const& enclosure : orbit.final()) {
        expect_time(enclosure,2.0);
        auto count=enclosure.previous_events().size();
        ARIADNE_TEST_ASSERT(count<=1u);
        if(count==0u) {
            flowed=true;
            ARIADNE_TEST_ASSERT(definitely(enclosure.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{ExactDouble(1.8),ExactDouble(2.3)},{ExactDouble(0.8),ExactDouble(1.2)}})))));
        } else {
            jumped=true;
            ARIADNE_TEST_EQUALS(enclosure.previous_events()[0],e);
            ARIADNE_TEST_ASSERT(definitely(enclosure.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{-ExactDouble(0.2),ExactDouble(0.3)},{ExactDouble(0.8),ExactDouble(1.2)}})))));
        }
    }
    ARIADNE_TEST_ASSERT(jumped && flowed);
}

Void TestHybridEvolver::test_initially_active_event() const {
    HybridAutomaton automaton;
    automaton.new_mode(q,{dot(x)=one,dot(y)=one});
    automaton.new_transition(q,e,q,{next(x)=x+1,next(y)=y},-x>=0,EventKind::URGENT);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBoxSet initial(q,{-1.625_dy<=x<=-1.375_dy,-0.125_dy<=y<=0.125_dy});
    HybridTime time(1.0_x,4);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(2.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);
    ARIADNE_TEST_EQUALS(orbit.final().size(),1u);

    ARIADNE_TEST_EQUALS(orbit.reach().size(),3u);



    expect_final(orbit,ExactBoxType({{ExactDouble(1.374),ExactDouble(1.626)},{ExactDouble(0.874),ExactDouble(1.126)}}),1.0,2u);
}

Void TestHybridEvolver::test_simultaneous_events() const {
    HybridAutomaton automaton;
    DiscreteEvent e1("e1");
    DiscreteEvent e2("e2");
    automaton.new_mode(q,{dot(x)=one,dot(y)=one});
    automaton.new_transition(q,e1,q,{next(x)=x-1,next(y)=y-2},x-1>=0,EventKind::URGENT);
    automaton.new_transition(q,e2,q,{next(x)=x-2,next(y)=y-1},y-1>=0,EventKind::URGENT);

    RealSpace space=automaton.continuous_state_space(q);
    HybridBoxSet initial(q,{-quarter<=x<=eighth,-eighth<=y<=quarter});
    HybridTime time(2.5_x,4);

    _set_evolver(automaton);
    evolver_ptr->configuration().set_maximum_step_size(4.0_x);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial,time,Semantics::UPPER);
    ARIADNE_TEST_EQUALS(orbit.final().size(),2u);

    Bool first_order=false,second_order=false;
    for(auto const& enclosure : orbit.final()) {
        expect_time(enclosure,2.5);
        ARIADNE_TEST_ASSERT(enclosure.previous_events().size()>0u);
        first_order=first_order || enclosure.previous_events()[0]==e1;
        second_order=second_order || enclosure.previous_events()[0]==e2;
    }
    ARIADNE_TEST_ASSERT(first_order && second_order);
}

Void TestHybridEvolver::test_transverse_linear_crossing() const
{
    DiscreteLocation q1(1), q2(2);
    Real r(0.125_x);
    ExactDouble tol=1e-5_pr;
    RealExpression guard=x+y/2-1;

    HybridAutomaton system;
    system.new_mode(q1,{dot(x)=1,dot(y)=0});
    system.new_mode(q2,{dot(x)=0,dot(y)=1});
    system.new_transition(q1,e,q2,{prime(x)=x,prime(y)=y},guard>=0,EventKind::URGENT);
    HybridSet initial_set(q1,{-r<=x<=r,-r<=y<=r});
    HybridTime evolution_time(2.0_x,3);

    _set_evolver(system);

    Orbit<HybridEnclosure> orbit=evolver_ptr->orbit(initial_set,evolution_time,Semantics::UPPER);
    Enclosure final_enclosure=orbit.final()[0].continuous_set();

    RealSpace final_space=system.continuous_state_space(q2);
    RealExpression ct=-guard; // Crossing time
    EffectiveVectorMultivariateFunction function=join(make_function(final_space,x+ct),make_function(final_space,y+2-ct));
    LabelledEnclosure expected_final_enclosure(initial_set.euclidean_set(q1,final_space),system.continuous_state_space(q1),EnclosureConfiguration(evolver_ptr->function_factory()));
    expected_final_enclosure.apply_map(function);

    Vector<FloatDPBounds> tolerance(2,FloatDPBounds(-tol,+tol,dp));
    ARIADNE_TEST_BINARY_PREDICATE(refines,expected_final_enclosure.state_function()[0],final_enclosure.state_function()[0]);
    ARIADNE_TEST_BINARY_PREDICATE(refines,final_enclosure.state_function()[0],expected_final_enclosure.state_function()[0]+tolerance[0]);

    Real xl(-1.0_x), xu(3.0_x), yl(-1.0_x), yu(3.0_x);
    ARIADNE_TEST_EQUALS(orbit.final().size(),1u);
    ARIADNE_TEST_EQUALS(orbit.final()[0].location(),q2);
    ARIADNE_TEST_EQUALS(orbit.final()[0].previous_events().size(),1u);
    ARIADNE_TEST_EQUALS(orbit.final()[0].previous_events()[0],e);
    expect_time(orbit.final()[0],2.0);
    ARIADNE_TEST_BINARY_PREDICATE(refines,expected_final_enclosure.state_function()[1],final_enclosure.state_function()[1]);
    ARIADNE_TEST_BINARY_PREDICATE(refines,final_enclosure.state_function()[1],expected_final_enclosure.state_function()[1]+tolerance[1]);
}

int main() {
    redirect_test_logging_to_null();
    TestHybridEvolver().test_all();
    { LogVerbosityScope diagnostics(100u); TestHybridEvolver().test_all(); }
    return ARIADNE_TEST_FAILURES;
}
