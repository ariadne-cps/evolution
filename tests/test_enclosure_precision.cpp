/***************************************************************************
 *            test_enclosure_precision.cpp
 *
 *  Regression test for precision loss in RealBox enclosure construction.
 ****************************************************************************/

#include "algebra/vector.hpp"
#include "function/taylor_function.hpp"
#include "function/constraint.hpp"
#include "geometry/box.hpp"
#include "evolution/enclosure.hpp"

#include "ariadne-testing.hpp"
#include "test_logging.hpp"

using namespace Ariadne;


class TestEnclosurePrecision : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    void test_contracts() {
    
    RealBox box({{Decimal("1.25"),Decimal("1.55")},
                 {Decimal("2.35"),Decimal("2.45")}});

    TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp,1e-14));
    Enclosure enclosure(box,EnclosureConfiguration(function_factory));

    // The RealBox-to-Enclosure conversion must not quantise the endpoints to
    // Float32.  The old conversion left a remainder of about 9.5e-8 here.
    auto const& state_model =
        dynamic_cast<ValidatedVectorMultivariateTaylorFunctionModelDP const&>(
            enclosure.state_function().reference());
    for(SizeType i=0; i!=state_model.size(); ++i) {
        ARIADNE_TESTING_EXPECT_LESS(state_model[i].error().raw(),1e-14_pr);
    }

    

    }
};

int main(int argc, char** argv) {
    redirect_test_logging_to_null();
    ARIADNE_TESTING_SUITE("test_enclosure_precision");
    ARIADNE_TESTING_ADD_TEST(TestEnclosurePrecision, test_contracts);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
