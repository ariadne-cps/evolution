// The policy is independently usable through its public header.
#include "evolution/scalings.hpp"
#include <sstream>
#include <string>
#include "macro/testing.hpp"
using namespace Ariadne;

class TestScalings {
public:
    Void test() {
        ARIADNE_TEST_CALL(test_named_and_default_grid());
        ARIADNE_TEST_CALL(test_positive_scaling_contracts());
    }
    Void test_named_and_default_grid() {
        RealVariable x("x"),y("y");
        Scalings uniform(ApproximateDouble(2.0));
        ARIADNE_TEST_EQUALS(uniform.scaling(x),2.0_x);
        uniform.set_scaling(x,ApproximateDouble(4.0));
        ARIADNE_TEST_EQUALS(uniform.scaling(x),4.0_x);
        ARIADNE_TEST_EQUALS(uniform.scaling(y),2.0_x);
        Map<RealVariable,ExactDouble> entries{{x,2.0_x}};
        Scalings named(3.0_x,entries);
        ARIADNE_TEST_EQUALS(named.scaling(x),2.0_x);
        ARIADNE_TEST_EQUALS(named.scaling(y),3.0_x);
        Grid grid=named.grid(RealSpace({y,x}));
        ARIADNE_TEST_EQUALS(grid.dimension(),2u);
        ARIADNE_TEST_EQUALS(grid.lengths()[0].get_d(),3.0);
        ARIADNE_TEST_EQUALS(grid.lengths()[1].get_d(),2.0);
        Scalings empty(1.0_x,Map<RealVariable,ExactDouble>{});
        ARIADNE_TEST_EQUALS(empty.grid(RealSpace()).dimension(),0u);
        std::ostringstream description; description<<named;
        ARIADNE_TEST_ASSERT(description.str().find("Scalings")!=std::string::npos);
        ARIADNE_TEST_ASSERT(description.str().find("x")!=std::string::npos);
    }
    Void test_positive_scaling_contracts() {
        RealVariable x("x");
        ARIADNE_TEST_FAIL(Scalings(ApproximateDouble(0.0)));
        ARIADNE_TEST_FAIL(Scalings(ApproximateDouble(-1.0)));
        const Map<RealVariable,ExactDouble> empty;
        ARIADNE_TEST_FAIL(Scalings(0.0_x,empty));
        Map<RealVariable,ExactDouble> zero{{x,0.0_x}};
        Map<RealVariable,ExactDouble> negative{{x,-1.0_x}};
        ARIADNE_TEST_FAIL(Scalings(1.0_x,zero));
        ARIADNE_TEST_FAIL(Scalings(1.0_x,negative));
        Scalings policy(ApproximateDouble(2.0));
        ARIADNE_TEST_FAIL(policy.set_scaling(x,ApproximateDouble(0.0)));
        ARIADNE_TEST_FAIL(policy.set_scaling(x,ApproximateDouble(-1.0)));
        ARIADNE_TEST_EQUALS(policy.scaling(x),2.0_x);
    }
};
Int main() { TestScalings().test(); return ARIADNE_TEST_FAILURES; }
