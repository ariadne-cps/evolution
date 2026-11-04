#include "evolution/hybrid_evolver.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "system/hybrid_automata.hpp"
#include "function/function.hpp"
#include "integration/integrator.hpp"
#include "macro/testing.hpp"
#include <cmath>
#include <memory>
using namespace Ariadne;

class TestHybridEvolver {
    inline static RealVariable x{"x"},y{"y"};
    inline static DiscreteLocation q;
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


    expect_final(orbit,ExactBoxType({{2.374_x,2.626_x},{1.124_x,1.376_x}}),2.5,0u);
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

    expect_final(orbit,ExactBoxType({{1.874_x,2.126_x},{0.874_x,1.126_x}}),2.0,0u);
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
        ARIADNE_TEST_ASSERT(definitely(enclosure.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{-1.1_x,-0.9_x},{-0.2_x,0.6_x}})))));
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

    expect_final(orbit,ExactBoxType({{-0.626_x,-0.374_x},{0.374_x,0.626_x}}),1.0,1u);
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

    expect_final(orbit,ExactBoxType({{0.874_x,1.126_x},{0.374_x,0.626_x}}),1.0,0u);
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
            ARIADNE_TEST_ASSERT(definitely(enclosure.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{1.8_x,2.3_x},{0.8_x,1.2_x}})))));
        } else {
            jumped=true;
            ARIADNE_TEST_EQUALS(enclosure.previous_events()[0],e);
            ARIADNE_TEST_ASSERT(definitely(enclosure.inside(HybridExactBox(q,RealSpace({x,y}),ExactBoxType({{-0.2_x,0.3_x},{0.8_x,1.2_x}})))));
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



    expect_final(orbit,ExactBoxType({{1.374_x,1.626_x},{0.874_x,1.126_x}}),1.0,2u);
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

int main() { TestHybridEvolver().test_all(); return ARIADNE_TEST_FAILURES; }
