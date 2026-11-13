/***************************************************************************
 *            test_enclosure.cpp
 *
 *  Copyright  2006-21  Luca Geretti
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

/** @file tests/test_enclosure.cpp
 *  @brief Test enclosure behavior and numerical invariants.
 */

#include <fstream>
#include <iostream>

#include "algebra/algebra.hpp"
#include "algebra/matrix.hpp"
#include "algebra/vector.hpp"
#include "drawing/drawer.hpp"
#include "drawing/drawer_manager.hpp"
#include "evolution/enclosure.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "geometry/box.hpp"
#include "integration/integrator.hpp"
#include "labelled/figure.hpp"
#include "logging/logging.hpp"
#include "numeric/rational.hpp"
#include "set/function_set.hpp"
#include "set/list_set.hpp"
#include "utility/tuple.hpp"

#include "ariadne-testing.hpp"
#include "expression-set/expression_set.hpp"

using namespace Ariadne;
using namespace std;

class TestEnclosure : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

  public:
    Void test_construct_without_domain() {
        Enclosure encl;
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp, 1e-8));
        EnclosureConfiguration configuration(function_factory);
        Enclosure encl_config(configuration);
    }

    Void test_construct_with_domain() {
        ExactBoxType dom({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp, 1e-8));
        EnclosureConfiguration configuration(function_factory);
        Enclosure encl(dom, configuration);
        ARIADNE_TESTING_NOTE(configuration);
        ARIADNE_TESTING_NOTE(encl);
        ARIADNE_TESTING_EXPECT_EQUAL(encl.dimension(), 2);
        ARIADNE_TESTING_EXPECT_EQUAL(encl.centre(), FloatDPPoint({FloatDP(1.0_x, DoublePrecision()), FloatDP(2.0_x, DoublePrecision())}));
        ARIADNE_TESTING_EXPECT_EQUAL(encl.radius().raw(), FloatDP(1.0_x, DoublePrecision()));
        ARIADNE_TESTING_PROVE(dom.inside(encl.codomain()));
    }

    Void test_is_bounded() {
        ExactBoxType bounded_dom({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        EnclosureConfiguration configuration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
        Enclosure encl(bounded_dom, configuration);
        ARIADNE_TESTING_PROVE(encl.is_bounded());
        ExactBoxType unbounded_dom({
            {0.0_x, infty},
            {1.0_x, 3.0_x}
        });
        Enclosure encl2(unbounded_dom, configuration);
        ARIADNE_TESTING_PROVE(is_indeterminate(encl2.is_bounded()));
    }

    Void test_subset_separated() {
        ExactBoxType box1({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        Enclosure encl(box1, EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8))));
        ARIADNE_TESTING_PROVE(definitely(encl.subset(box1)));
        ARIADNE_TESTING_PROVE(possibly(not encl.separated(box1)));
        ExactBoxType box2({
            { 2.0_x, 3.0_x},
            {-1.0_x, 0.0_x}
        });
        ARIADNE_TESTING_PROVE(possibly(not encl.subset(box2)));
        ARIADNE_TESTING_PROVE(definitely(encl.separated(box2)));
    }

    Void test_restriction() {
        ExactBoxType box1({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        Enclosure encl(box1, EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8))));
        auto codomain1 = encl.codomain();
        ExactBoxType box2({
            {0.5_x, 1.5_x},
            {1.5_x, 2.5_x}
        });
        auto encl2 = restriction(encl, box2);
        ARIADNE_TESTING_PROVE(definitely(encl2.subset(codomain1)));
    }

    // Discard all original parameters and preserve the represented set.
    Void test_kuhn_reconditioning_encloses_original() {
        ExactBoxType original({
            {-1.0_x, 1.0_x},
            {-2.0_x, 2.0_x}
        });
        auto swp = ThresholdSweeper<FloatDP>(dp, 1e-8);
        EnclosureConfiguration configuration{TaylorFunctionFactory(swp)};
        configuration.set_reconditioning_num_blocks(1u);
        Enclosure encl(original, configuration);
        ExactBoxType before = encl.codomain();
        encl.kuhn_recondition();
        ARIADNE_TESTING_EXPECT_EQUAL(encl.number_of_parameters(), 2u);
        // Non-strict inclusion is required: equal boundaries are permitted.
        ARIADNE_TESTING_NOTE(before);
        ARIADNE_TESTING_NOTE(encl.codomain());
        ARIADNE_TESTING_PROVE(before.subset(encl.codomain()));
    }

    Void test_auxiliary_map() {
        ExactBoxType dom({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        Enclosure encl(dom, EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8))));
        ARIADNE_TESTING_NOTE(encl.auxiliary_function());
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2, 0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2, 1);
        EffectiveVectorMultivariateFunction auxiliary(1, x0 + sqr(x1));
        encl.set_auxiliary_mapping(auxiliary);
        ARIADNE_TESTING_NOTE(encl.state_auxiliary_function());
        ARIADNE_TESTING_NOTE(encl.state_auxiliary_set());
        ARIADNE_TESTING_EXPECT_EQUAL(encl.dimension(), 3);
        for (SizeType i = 0; i < 4; ++i)
            ARIADNE_TESTING_NOTE(encl.get_function(i));
    }

    Void test_constraints() {
        ExactBoxType dom({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        Enclosure encl(dom, EnclosureConfiguration(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8))));
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2, 0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2, 1);
        encl.new_zero_state_constraint(x0 + sqr(x1));
        encl.new_positive_state_constraint(x0);
        ARIADNE_TESTING_EXPECT_EQUAL(encl.number_of_constraints(), 2);
        ARIADNE_TESTING_NOTE(encl.constraint_function());
        ARIADNE_TESTING_NOTE(encl.constraint_models());
        ARIADNE_TESTING_NOTE(encl.constraint(0));
        ARIADNE_TESTING_NOTE(encl.constraint(1));
        ARIADNE_TESTING_NOTE(encl.constraint_bounds());
        ARIADNE_TESTING_PROVE(definitely(encl.satisfies(x1 + x0)));
        ARIADNE_TESTING_PROVE(definitely(not encl.satisfies(-sqr(x1))));
        ARIADNE_TESTING_PROVE(is_indeterminate(encl.satisfies(-sqr(x0))));
    }

    Void test_split() {
        ExactBoxType dom1({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        auto swp = ThresholdSweeper<FloatDP>(dp, 1e-8);
        auto fnc = antiderivative(ValidatedVectorMultivariateTaylorFunctionModelDP::identity(dom1, swp), 0);
        Enclosure encl(dom1, fnc, EnclosureConfiguration(TaylorFunctionFactory(swp)));
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2, 0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2, 1);
        EffectiveVectorMultivariateFunction auxiliary(1, x0 + sqr(x1));
        encl.set_auxiliary_mapping(auxiliary);
        auto split1 = encl.split();
        ARIADNE_TESTING_EXPECT_EQUAL(split1.first.auxiliary_function().result_size(), 1);
        auto split2 = encl.split(1);
        ARIADNE_TESTING_EXPECT_EQUAL(split2.first.auxiliary_function().result_size(), 1);
        auto split3 = encl.split_first_order();
        ARIADNE_TESTING_EXPECT_EQUAL(split3.first.auxiliary_function().result_size(), 1);
        ARIADNE_TESTING_NOTE(encl.splitting_subdomains_zeroth_order());
    }

    Void test_labelled_construction() {
        RealVariable x("x"), y("y"), z("z"), t("t");
        RealSpace spc({x, y});
        ExactBoxType bx({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        LabelledExactBoxType dom(spc, bx);
        TaylorFunctionFactory function_factory(ThresholdSweeper<FloatDP>(dp, 1e-8));
        EnclosureConfiguration configuration(function_factory);
        LabelledEnclosure labelled_encl1(dom, configuration);
        ARIADNE_TESTING_NOTE(labelled_encl1);
        labelled_encl1.set_state_space(RealSpace({y, z}));
        ARIADNE_TESTING_NOTE(labelled_encl1);
        Enclosure encl2(bx, configuration);
        LabelledEnclosure labelled_encl2(encl2, spc);
        ARIADNE_TESTING_NOTE(labelled_encl2);
        Enclosure encl3(bx, configuration);
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2, 0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2, 1);
        EffectiveVectorMultivariateFunction auxiliary(1, x0 + sqr(x1));
        encl3.set_auxiliary_mapping(auxiliary);
        LabelledEnclosure labelled_encl3(encl3, spc, RealSpace({z}));
        ARIADNE_TESTING_NOTE(labelled_encl3);
        ARIADNE_TESTING_NOTE(labelled_encl3.state_auxiliary_space());
        LabelledEnclosure labelled_encl4(encl3, RealSpace({x, t}), RealSpace({z}));
        ARIADNE_TESTING_NOTE(labelled_encl4.state_time_auxiliary_space());
    }

    Void test_labelled_product() {
        RealVariable x("x"), y("y"), z("z");
        ExactBoxType dom1({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        auto swp = ThresholdSweeper<FloatDP>(dp, 1e-8);
        auto fnc = antiderivative(ValidatedVectorMultivariateTaylorFunctionModelDP::identity(dom1, swp), 0);
        Enclosure e1(dom1, fnc, EnclosureConfiguration(TaylorFunctionFactory(swp)));
        LabelledEnclosure le1(e1, RealSpace({x, y}), RealSpace({z}));
        ARIADNE_TESTING_NOTE(product(le1, LabelledExactIntervalType(z, ExactIntervalType(-1.0_x, 1.0_x))));
        ExactBoxType dom2({
            {-1.0_x, 1.0_x}
        });
        LabelledExactBoxType box2(RealSpace({z}), dom2);
        ARIADNE_TESTING_NOTE(product(le1, box2));
        Enclosure e2(dom2, EnclosureConfiguration(TaylorFunctionFactory(swp)));
        LabelledEnclosure le2(e2, RealSpace({z}));
        ARIADNE_TESTING_NOTE(product(le1, le2));
    }

    Void _draw(const Drawer& drawer, String suffix) {
        ARIADNE_TESTING_NOTE("Drawing with " + suffix + " method");
        ARIADNE_TESTING_NOTE(drawer);
        DrawerManager::instance().set_drawer(drawer);
        ExactBoxType dom({
            {0.0_x, 1.0_x},
            {1.0_x, 2.0_x}
        });
        auto swp = ThresholdSweeper<FloatDP>(dp, 1e-8);
        auto fnc = antiderivative(ValidatedVectorMultivariateTaylorFunctionModelDP::identity(dom, swp), 0);
        TaylorFunctionFactory function_factory(swp);
        EnclosureConfiguration configuration(function_factory);
        Enclosure encl(dom, fnc, configuration);
        RealVariable x("x"), y("y");
        LabelledEnclosure lencl(encl, RealSpace({x, y}));
        LabelledFigure fig(Axes2d(-1.0_dec <= x <= 1.0_dec, -1.5_dec <= y <= 1.5_dec));
        fig << line_style(true) << lencl;
        fig.write(("test_enclosure_" + suffix).c_str());
    }

    Void test_labelled_drawers() {
        _draw(BoxDrawer(), "box");
        _draw(GridDrawer(4), "grid");
        _draw(AffineDrawer(1), "affine");
    }

    void test_labelled_with_auxiliary_draw() {
        ExactBoxType dom({
            {0.0_x, 2.0_x},
            {1.0_x, 3.0_x}
        });
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
        Enclosure encl(dom, config);
        auto x0 = EffectiveScalarMultivariateFunction::coordinate(2, 0);
        auto x1 = EffectiveScalarMultivariateFunction::coordinate(2, 1);
        EffectiveVectorMultivariateFunction auxiliary(1, x0 + sqr(x1));
        encl.set_auxiliary_mapping(auxiliary);
        RealVariable x("x"), y("y"), z("z");
        LabelledEnclosure lencl(encl, RealSpace({x, y}), RealSpace({z}));
        LabelledFigure fig(Axes2d(-1.0_dec <= x <= 3, -1 <= z <= 11));
        fig << lencl;
        DrawerManager::instance().set_drawer(AffineDrawer(1));
        fig.write("test_enclosure_auxiliary");
    }
};

int main(int argc, char* argv[]) {
    ARIADNE_TESTING_SUITE("test_enclosure");
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_construct_without_domain);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_construct_with_domain);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_is_bounded);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_subset_separated);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_restriction);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_kuhn_reconditioning_encloses_original);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_auxiliary_map);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_constraints);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_split);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_labelled_construction);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_labelled_product);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_labelled_drawers);
    ARIADNE_TESTING_ADD_TEST(TestEnclosure, test_labelled_with_auxiliary_draw);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
