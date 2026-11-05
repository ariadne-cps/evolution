// The policy is independently usable through its public header.
#include "evolution/scalings.hpp"
#include <sstream>
#include <string>
#include "ariadne-testing.hpp"
#include "test_logging.hpp"
using namespace Ariadne;

class TestScalings : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

public:
    
    Void test_named_and_default_grid() {
        RealVariable x("x"),y("y");
        Scalings uniform(ApproximateDouble(2.0));
        ARIADNE_TESTING_EXPECT_EQUAL(uniform.scaling(x),2.0_x);
        uniform.set_scaling(x,ApproximateDouble(4.0));
        ARIADNE_TESTING_EXPECT_EQUAL(uniform.scaling(x),4.0_x);
        ARIADNE_TESTING_EXPECT_EQUAL(uniform.scaling(y),2.0_x);
        Map<RealVariable,ExactDouble> entries{{x,2.0_x}};
        Scalings named(3.0_x,entries);
        ARIADNE_TESTING_EXPECT_EQUAL(named.scaling(x),2.0_x);
        ARIADNE_TESTING_EXPECT_EQUAL(named.scaling(y),3.0_x);
        Grid grid=named.grid(RealSpace({y,x}));
        ARIADNE_TESTING_EXPECT_EQUAL(grid.dimension(),2u);
        ARIADNE_TESTING_EXPECT_EQUAL(grid.lengths()[0].get_d(),3.0);
        ARIADNE_TESTING_EXPECT_EQUAL(grid.lengths()[1].get_d(),2.0);
        Scalings empty(1.0_x,Map<RealVariable,ExactDouble>{});
        ARIADNE_TESTING_EXPECT_EQUAL(empty.grid(RealSpace()).dimension(),0u);
        std::ostringstream description; description<<named;
        ARIADNE_TESTING_PROVE(description.str().find("Scalings")!=std::string::npos);
        ARIADNE_TESTING_PROVE(description.str().find("x")!=std::string::npos);
    }
    Void test_positive_scaling_contracts() {
        RealVariable x("x");
        ARIADNE_TESTING_EXPECT_EXCEPTION(Scalings(ApproximateDouble(0.0)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(Scalings(ApproximateDouble(-1.0)));
        const Map<RealVariable,ExactDouble> empty;
        ARIADNE_TESTING_EXPECT_EXCEPTION(Scalings(0.0_x,empty));
        Map<RealVariable,ExactDouble> zero{{x,0.0_x}};
        Map<RealVariable,ExactDouble> negative{{x,-1.0_x}};
        ARIADNE_TESTING_EXPECT_EXCEPTION(Scalings(1.0_x,zero));
        ARIADNE_TESTING_EXPECT_EXCEPTION(Scalings(1.0_x,negative));
        Scalings policy(ApproximateDouble(2.0));
        ARIADNE_TESTING_EXPECT_EXCEPTION(policy.set_scaling(x,ApproximateDouble(0.0)));
        ARIADNE_TESTING_EXPECT_EXCEPTION(policy.set_scaling(x,ApproximateDouble(-1.0)));
        ARIADNE_TESTING_EXPECT_EQUAL(policy.scaling(x),2.0_x);
    }
};

int main(int argc, char* argv[]) {
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_scalings");
    ARIADNE_TESTING_ADD_TEST(TestScalings, test_named_and_default_grid);
    ARIADNE_TESTING_ADD_TEST(TestScalings, test_positive_scaling_contracts);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
