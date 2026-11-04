#include "evolution/hybrid_time.hpp"
#include "macro/testing.hpp"
#include <sstream>
#include <string>
#include "test_logging.hpp"
using namespace Ariadne;

class TestHybridTime {
public:
    Void check_continuous(const Real& value,double expected) {
        auto bounds=value.get(dp);
        ARIADNE_TEST_EQUALS(bounds.lower_raw().get_d(),expected);
        ARIADNE_TEST_EQUALS(bounds.upper_raw().get_d(),expected);
    }
    Void test() {
        ARIADNE_TEST_CALL(test_construction_and_arithmetic());
        ARIADNE_TEST_CALL(test_comparison_semantics());
        ARIADNE_TEST_CALL(test_termination());
    }
    Void test_construction_and_arithmetic() {
        HybridTime unlimited;
        check_continuous(unlimited.continuous_time(),0.0);
        ARIADNE_TEST_EQUALS(unlimited.discrete_time(),Integer(-1));
        HybridTime continuous_only(Real(2));
        check_continuous(continuous_only.continuous_time(),2.0);
        ARIADNE_TEST_EQUALS(continuous_only.discrete_time(),Integer(-1));
        HybridTime real(Real(2),Integer(3));
        HybridTime raw(FloatDP(2,dp),Integer(3));
        HybridTime signed_steps(2.0_x,Int(3));
        HybridTime natural_steps(2.0_x,Nat(3));
        for(const auto& time : {real,raw,signed_steps,natural_steps}) {
            check_continuous(time.continuous_time(),2.0);
            ARIADNE_TEST_EQUALS(time.discrete_time(),Integer(3));
        }
        HybridTime scaled=Integer(2)*real;
        check_continuous(scaled.continuous_time(),4.0);
        ARIADNE_TEST_EQUALS(scaled.discrete_time(),Integer(6));
        HybridTime difference=scaled-real;
        check_continuous(difference.continuous_time(),2.0);
        ARIADNE_TEST_EQUALS(difference.discrete_time(),Integer(3));
        HybridTime clamped=real-scaled;
        check_continuous(clamped.continuous_time(),-2.0);
        ARIADNE_TEST_EQUALS(clamped.discrete_time(),Integer(0));
        std::ostringstream text; text<<real;
        ARIADNE_TEST_ASSERT(text.str().find("steps=3")!=std::string::npos);
    }
    Void test_comparison_semantics() {
        HybridTime earlier(Real(1),Integer(1)),later(Real(2),Integer(2));
        HybridTime same_continuous=earlier; same_continuous._discrete_time=Integer(3);
        ARIADNE_TEST_ASSERT(same(earlier,earlier));
        ARIADNE_TEST_ASSERT(!same(earlier,later));
        ARIADNE_TEST_ASSERT(!same(earlier,same_continuous));
        // Equality of arbitrary computable reals is not decidable. Test
        // robust falsification and inequality, and strict ordered cases.
        ARIADNE_TEST_ASSERT(!possibly(check(earlier==later,Effort(1))));
        ARIADNE_TEST_ASSERT(definitely(check(earlier!=later,Effort(1))));
        ARIADNE_TEST_ASSERT(definitely(check(earlier<later,Effort(1))));
        ARIADNE_TEST_ASSERT(!possibly(check(later<earlier,Effort(1))));
        ARIADNE_TEST_ASSERT(definitely(check(earlier<=later,Effort(1))));
        ARIADNE_TEST_ASSERT(!possibly(check(later<=earlier,Effort(1))));
        ARIADNE_TEST_ASSERT(definitely(check(later>earlier,Effort(1))));
        ARIADNE_TEST_ASSERT(!possibly(check(earlier>later,Effort(1))));
        ARIADNE_TEST_ASSERT(definitely(check(later>Real(1),Effort(1))));
        ARIADNE_TEST_ASSERT(!possibly(check(earlier>Real(2),Effort(1))));
        HybridTime fewer_steps(Real(2),Integer(0));
        ARIADNE_TEST_ASSERT(!possibly(check(earlier<fewer_steps,Effort(1))));
        ARIADNE_TEST_ASSERT(!possibly(check(earlier<=fewer_steps,Effort(1))));
    }
    Void test_termination() {
        DiscreteEvent stop("stop");
        HybridTerminationCriterion with_events(Real(4),Integer(2),Set<DiscreteEvent>{stop});
        HybridTerminationCriterion plain(Real(4),Integer(2));
        HybridTerminationCriterion from_time(HybridTime(Real(4),Integer(2)));
        check_continuous(with_events.maximum_time(),4.0);
        ARIADNE_TEST_EQUALS(with_events.maximum_steps(),Integer(2));
        ARIADNE_TEST_ASSERT(with_events.terminating_events().contains(stop));
        ARIADNE_TEST_ASSERT(plain.terminating_events().empty());
        ARIADNE_TEST_ASSERT(from_time.terminating_events().empty());
        std::ostringstream text; text<<with_events;
        ARIADNE_TEST_ASSERT(text.str().find("stop")!=std::string::npos);
        ARIADNE_TEST_ASSERT(text.str().find("maximum_steps=2")!=std::string::npos);
    }
};
Int main() {
    redirect_test_logging_to_null(); TestHybridTime().test(); return ARIADNE_TEST_FAILURES; }
