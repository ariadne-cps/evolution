/***************************************************************************
 *            test_orbit_collections.cpp
 *
 *  Copyright  2026  Luca Geretti
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

/** @file tests/test_orbit_collections.cpp
 *  @brief Test orbit collections behavior and numerical invariants.
 */

#include <memory>
#include <sstream>

#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "evolution/orbit.hpp"
#include "function/taylor_function.hpp"
#include "labelled/figure.hpp"
#include "recording_canvas.hpp"
#include "rendering/figure.hpp"

#include "ariadne-testing.hpp"
using namespace Ariadne;

class TestOrbitCollections : public Ariadne::Testing::ContextFixture {
  public:
    using Ariadne::Testing::ContextFixture::ContextFixture;

    RealVariable x{"x"}, y{"y"}, z{"z"};
    RealSpace space{x, y};
    FloatDPApproximationPoint p{
        {1.0_x, 2.0_x},
        dp
    };
    FloatDPApproximationPoint next{
        {3.0_x, 4.0_x},
        dp
    };
    FloatDP time{1.0_x, dp};

  public:
    Void check_endpoint(const InterpolatedCurve& curve, SizeType expected_size) {
        ARIADNE_TESTING_EXPECT_EQUAL(curve.size(), expected_size);
        ARIADNE_TESTING_EXPECT_EQUAL(curve.end()->second[0].raw().get_d(), 3.0);
        ARIADNE_TESTING_EXPECT_EQUAL(curve.end()->second[1].raw().get_d(), 4.0);
    }
    Void test_point_collections() {
        Orbit<FloatDPApproximationPoint> single(p);
        single.insert(time, next);
        check_endpoint(single.curve(), 2u);
        Orbit<List<FloatDPApproximationPoint>> empty(List<FloatDPApproximationPoint>{});
        ARIADNE_TESTING_EXPECT_EQUAL(empty.curves().size(), 0u);
        Orbit<List<FloatDPApproximationPoint>> multiple(List<FloatDPApproximationPoint>{p, next});
        multiple.insert(time, next, 0u);
        check_endpoint(multiple.curve(0), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(multiple.curves().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(multiple.curve(1).size(), 1u);
    }
    Void test_labelled_collections() {
        LabelledPoint<FloatDPApproximation> labelled(space, p), end(space, next);
        Orbit<LabelledPoint<FloatDPApproximation>> single(labelled);
        single.insert(time, end);
        check_endpoint(single.curves(), 2u);
        RecordingCanvas canvas;
        single.draw(canvas, Variables2d(x, y));
        ARIADNE_TESTING_PROVE(canvas.primitives > 0u);
        std::unique_ptr<LabelledDrawable2dInterface> copy(single.labelled_clone());
        RecordingCanvas copied;
        copy->draw(copied, Variables2d(x, y));
        ARIADNE_TESTING_EXPECT_EQUAL(copied.primitives, canvas.primitives);
        Orbit<List<LabelledPoint<FloatDPApproximation>>> empty(List<LabelledPoint<FloatDPApproximation>>{});
        ARIADNE_TESTING_EXPECT_EQUAL(empty.curves().size(), 0u);
        empty.draw(canvas, Variables2d(x, y));
        Orbit<List<LabelledPoint<FloatDPApproximation>>> multiple(List<LabelledPoint<FloatDPApproximation>>{labelled, end});
        multiple.insert(time, end, 0u);
        check_endpoint(multiple.curve(0u), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(multiple.curves().size(), 2u);
        RecordingCanvas drawn;
        multiple.draw(drawn, Variables2d(x, y));
        std::unique_ptr<LabelledDrawable2dInterface> cloned(multiple.labelled_clone());
        RecordingCanvas redrawn;
        cloned->draw(redrawn, Variables2d(x, y));
        ARIADNE_TESTING_PROVE(drawn.primitives > 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(redrawn.primitives, drawn.primitives);
    }
    template<class Shape>
    Void check_filter(const Shape& shape) {
        RecordingCanvas canvas;
        shape.draw(canvas, {DiscreteLocation(9)}, Variables2d(x, y));
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.primitives, 0u);
        shape.draw(canvas, {}, Variables2d(x, z));
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.primitives, 0u);
        shape.draw(canvas, {DiscreteLocation(1)}, Variables2d(x, y));
        ARIADNE_TESTING_PROVE(canvas.primitives > 0u);
        auto count = canvas.primitives;
        shape.draw(canvas, {}, Variables2d(x, y));
        ARIADNE_TESTING_PROVE(canvas.primitives > count);
    }
    Void test_hybrid_collections() {
        HybridApproximatePoint start(DiscreteLocation(1), space, p), end(DiscreteLocation(1), space, next);
        HybridApproximatePoint jump(DiscreteLocation(2), space, p);
        Orbit<HybridApproximatePoint> single(start);
        single.insert(HybridTime(1.0_x, 0), end);
        single.insert(HybridTime(1.0_x, 1), jump);
        ARIADNE_TESTING_EXPECT_EQUAL(single.size(), 2u);
        check_endpoint(single.curve(0u), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(single.curves()[1].location(), DiscreteLocation(2));
        ARIADNE_TESTING_EXPECT_EXCEPTION(single.insert(HybridTime(1.0_x, -1), end));
        ARIADNE_TESTING_EXPECT_EXCEPTION(single.insert(HybridTime(1.0_x, 3), end));
        check_filter(single);
        Orbit<HybridApproximatePoint> copied(single.curves());
        ARIADNE_TESTING_EXPECT_EQUAL(copied.size(), 2u);
        Orbit<HybridApproximatePoint> empty(List<HybridInterpolatedCurve>{});
        ARIADNE_TESTING_EXPECT_EQUAL(empty.size(), 0u);
        RecordingCanvas blank;
        empty.draw(blank, {}, Variables2d(x, y));
        ARIADNE_TESTING_EXPECT_EQUAL(blank.primitives, 0u);
        std::ostringstream text;
        text << single;
        ARIADNE_TESTING_PROVE(text.str().find("InterpolatedCurve") != std::string::npos);
        Orbit<Vector<HybridApproximatePoint>> multi(Vector<HybridApproximatePoint>{start, end});
        multi.insert(HybridTime(1.0_x, 0), end, 0u);
        multi.insert(HybridTime(1.0_x, 1), jump, 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(multi.size(0u), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(multi.size(1u), 1u);
        check_endpoint(multi.curve(0u, 0u), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(multi.curves(0u)[1].location(), DiscreteLocation(2));
        ARIADNE_TESTING_EXPECT_EXCEPTION(multi.insert(HybridTime(1.0_x, -1), end, 0u));
        ARIADNE_TESTING_EXPECT_EXCEPTION(multi.insert(HybridTime(1.0_x, 3), end, 0u));
        check_filter(multi);
        Orbit<Vector<HybridApproximatePoint>> no_points(Vector<HybridApproximatePoint>{});
        no_points.draw(blank, {}, Variables2d(x, y));
        ARIADNE_TESTING_EXPECT_EQUAL(blank.primitives, 0u);
    }
    template<class ES>
    Void check_adjoin(ES const& enclosure) {
        Orbit<ES> orbit(enclosure);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.initial().dimension(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.intermediate().size(), 0u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 0u);
        orbit.adjoin_reach(enclosure);
        orbit.adjoin_intermediate(enclosure);
        orbit.adjoin_final(enclosure);
        ListSet<ES> sets;
        sets.adjoin(enclosure);
        orbit.adjoin_reach(sets);
        orbit.adjoin_intermediate(sets);
        orbit.adjoin_final(sets);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.reach().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.intermediate().size(), 2u);
        ARIADNE_TESTING_EXPECT_EQUAL(orbit.final().size(), 2u);
        std::ostringstream description;
        description << orbit;
        ARIADNE_TESTING_PROVE(description.str().find("intermediate=") != std::string::npos);
    }
    Void test_enclosure_collections() {
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp, 1e-8)));
        Enclosure enclosure(ExactBoxType({
                                {0, 1},
                                {0, 1}
        }),
                            config);
        check_adjoin(enclosure);
        LabelledEnclosure labelled(enclosure, space);
        check_adjoin(labelled);
        Orbit<LabelledEnclosure> labelled_orbit(labelled);
        labelled_orbit.adjoin_reach(labelled);
        RecordingCanvas canvas;
        labelled_orbit.draw(canvas, Variables2d(x, y));
        ARIADNE_TESTING_PROVE(canvas.boundaries > 0u);
        std::unique_ptr<LabelledDrawable2dInterface> clone(labelled_orbit.labelled_clone());
        RecordingCanvas copy;
        clone->draw(copy, Variables2d(x, y));
        ARIADNE_TESTING_EXPECT_EQUAL(copy.boundaries, canvas.boundaries);
        HybridEnclosure hybrid(DiscreteLocation(1), labelled);
        check_adjoin(hybrid);
        Orbit<HybridEnclosure> hybrid_orbit(hybrid);
        hybrid_orbit.adjoin_reach(hybrid);
        RecordingCanvas hybrid_canvas;
        hybrid_orbit.draw(hybrid_canvas, {}, Variables2d(x, y));
        ARIADNE_TESTING_PROVE(hybrid_canvas.boundaries > 0u);
        Orbit<HybridEnclosure> no_initial;
        ARIADNE_TESTING_EXPECT_EQUAL(no_initial.reach().size(), 0u);
        Figure figure(Box<ApproximateIntervalType>({
                          {-1, 2},
                          {-1, 2}
        }),
                      Projection2d(2, 0, 1));
        Orbit<Enclosure> drawable(enclosure);
        drawable.adjoin_reach(enclosure);
        drawable.adjoin_final(enclosure);
        figure << drawable;
        auto recorded = std::make_shared<RecordingCanvas>();
        BackendScope backend{RenderingBackend(RecordingBackend(recorded))};
        figure.write("orbit-memory");
        ARIADNE_TESTING_PROVE(recorded->boundaries > 0u);
    }
};

int main(int argc, char* argv[]) {
    ARIADNE_TESTING_SUITE("test_orbit_collections");
    ARIADNE_TESTING_ADD_TEST(TestOrbitCollections, test_point_collections);
    ARIADNE_TESTING_ADD_TEST(TestOrbitCollections, test_labelled_collections);
    ARIADNE_TESTING_ADD_TEST(TestOrbitCollections, test_hybrid_collections);
    ARIADNE_TESTING_ADD_TEST(TestOrbitCollections, test_enclosure_collections);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
