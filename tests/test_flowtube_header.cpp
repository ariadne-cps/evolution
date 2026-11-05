#include "evolution/flowtube.hpp"
#include "function/taylor_function.hpp"
#include "ariadne-testing.hpp"

#include <type_traits>
#include <utility>
#include "test_logging.hpp"

// The interval returned by time_range() may be temporary: endpoints must not
// escape as references to that temporary.
using Ariadne::FlowTube;
static_assert(!std::is_reference_v<decltype(std::declval<FlowTube const&>().initial_time())>);
static_assert(!std::is_reference_v<decltype(std::declval<FlowTube const&>().final_time())>);

using namespace Ariadne;

class TestFlowtubeHeader : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    void test_contracts() {
    
    TaylorFunctionFactory factory(ThresholdSweeper<FloatDP>(dp,1e-8));
    ExactBoxType domain({{-1.0_x,1.0_x},{0.0_x,2.0_x}});
    auto coordinates=EffectiveVectorMultivariateFunction::identity(2);
    EffectiveVectorMultivariateFunction flow={coordinates[0]+coordinates[1]};
    FlowTube tube(factory.create(domain,flow));
    ARIADNE_TESTING_EXPECT_EQUAL(tube.number_of_state_variables(),1u);
    ARIADNE_TESTING_EXPECT_EQUAL(tube.space_domain(),ExactBoxType({{-1.0_x,1.0_x}}));
    ARIADNE_TESTING_EXPECT_EQUAL(tube.initial_time().get_d(),0.0);
    ARIADNE_TESTING_EXPECT_EQUAL(tube.final_time().get_d(),2.0);
    ARIADNE_TESTING_EXPECT_EQUAL(tube.reach_function().argument_size(),2u);
    ARIADNE_TESTING_EXPECT_EQUAL(tube.final_function().argument_size(),1u);
    Vector<FloatDPApproximation> point({0.5_x},dp);
    auto image=tube(point,FloatDPApproximation(1.0_x,dp));
    ARIADNE_TESTING_EXPECT_EQUAL(image[0].raw().get_d(),1.5);
    ARIADNE_TESTING_EXPECT_EQUAL(tube.initial_set().domain(),tube.space_domain());
    ARIADNE_TESTING_EXPECT_EQUAL(tube.reach_set().domain(),domain);
    ARIADNE_TESTING_EXPECT_EQUAL(tube.final_set().domain(),tube.space_domain());
    ARIADNE_TESTING_PROVE(definitely(tube.initial_set().bounding_box().inside(ExactBoxType({{-ExactDouble(1.1),ExactDouble(1.1)}}))));
    ARIADNE_TESTING_PROVE(definitely(tube.final_set().bounding_box().inside(ExactBoxType({{ExactDouble(0.9),ExactDouble(3.1)}}))));
    ARIADNE_TESTING_PROVE(definitely(tube.reach_set().bounding_box().inside(ExactBoxType({{-ExactDouble(1.1),ExactDouble(3.1)}}))));
    ARIADNE_TESTING_EXPECT_EXCEPTION(FlowTube(factory.create(domain,coordinates)));
    

    }
};

int main(int argc, char** argv) {
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_flowtube_header");
    ARIADNE_TESTING_ADD_TEST(TestFlowtubeHeader, test_contracts);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
