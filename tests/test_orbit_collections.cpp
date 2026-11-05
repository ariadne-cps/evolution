#include "evolution/orbit.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "function/taylor_function.hpp"
#include "labelled/figure.hpp"
#include "rendering/figure.hpp"
#include "recording_canvas.hpp"
#include "macro/testing.hpp"
#include <sstream>
#include <memory>
#include "test_logging.hpp"
using namespace Ariadne;

class TestOrbitCollections {
    RealVariable x{"x"},y{"y"},z{"z"};
    RealSpace space{x,y};
    FloatDPApproximationPoint p{{1.0_x,2.0_x},dp};
    FloatDPApproximationPoint next{{3.0_x,4.0_x},dp};
    FloatDP time{1.0_x,dp};
public:
    Void test() {
        ARIADNE_TEST_CALL(test_point_collections());
        ARIADNE_TEST_CALL(test_labelled_collections());
        ARIADNE_TEST_CALL(test_hybrid_collections());
        ARIADNE_TEST_CALL(test_enclosure_collections());
    }
    Void check_endpoint(InterpolatedCurve const& curve,SizeType size) {
        ARIADNE_TEST_EQUALS(curve.size(),size);
        ARIADNE_TEST_EQUALS(curve.end()->second[0].raw().get_d(),3.0);
        ARIADNE_TEST_EQUALS(curve.end()->second[1].raw().get_d(),4.0);
    }
    Void test_point_collections() {
        Orbit<FloatDPApproximationPoint> single(p);
        single.insert(time,next); check_endpoint(single.curve(),2u);
        Orbit<List<FloatDPApproximationPoint>> empty(List<FloatDPApproximationPoint>{});
        ARIADNE_TEST_EQUALS(empty.curves().size(),0u);
        Orbit<List<FloatDPApproximationPoint>> multiple(List<FloatDPApproximationPoint>{p,next});
        multiple.insert(time,next,0u); check_endpoint(multiple.curve(0),2u);
        ARIADNE_TEST_EQUALS(multiple.curves().size(),2u);
        ARIADNE_TEST_EQUALS(multiple.curve(1).size(),1u);
    }
    Void test_labelled_collections() {
        LabelledPoint<FloatDPApproximation> labelled(space,p),end(space,next);
        Orbit<LabelledPoint<FloatDPApproximation>> single(labelled);
        single.insert(time,end); check_endpoint(single.curves(),2u);
        RecordingCanvas canvas; single.draw(canvas,Variables2d(x,y));
        ARIADNE_TEST_ASSERT(canvas.primitives>0u);
        std::unique_ptr<LabelledDrawable2dInterface> copy(single.labelled_clone());
        RecordingCanvas copied; copy->draw(copied,Variables2d(x,y));
        ARIADNE_TEST_EQUALS(copied.primitives,canvas.primitives);
        Orbit<List<LabelledPoint<FloatDPApproximation>>> empty(List<LabelledPoint<FloatDPApproximation>>{});
        ARIADNE_TEST_EQUALS(empty.curves().size(),0u);
        empty.draw(canvas,Variables2d(x,y));
        Orbit<List<LabelledPoint<FloatDPApproximation>>> multiple(List<LabelledPoint<FloatDPApproximation>>{labelled,end});
        multiple.insert(time,end,0u); check_endpoint(multiple.curve(0u),2u);
        ARIADNE_TEST_EQUALS(multiple.curves().size(),2u);
        RecordingCanvas drawn; multiple.draw(drawn,Variables2d(x,y));
        std::unique_ptr<LabelledDrawable2dInterface> cloned(multiple.labelled_clone());
        RecordingCanvas redrawn; cloned->draw(redrawn,Variables2d(x,y));
        ARIADNE_TEST_ASSERT(drawn.primitives>0u);
        ARIADNE_TEST_EQUALS(redrawn.primitives,drawn.primitives);
    }
    template<class Shape> Void check_filter(Shape const& shape) {
        RecordingCanvas canvas;
        shape.draw(canvas,{DiscreteLocation(9)},Variables2d(x,y));
        ARIADNE_TEST_EQUALS(canvas.primitives,0u);
        shape.draw(canvas,{},Variables2d(x,z));
        ARIADNE_TEST_EQUALS(canvas.primitives,0u);
        shape.draw(canvas,{DiscreteLocation(1)},Variables2d(x,y));
        ARIADNE_TEST_ASSERT(canvas.primitives>0u);
        auto count=canvas.primitives;
        shape.draw(canvas,{},Variables2d(x,y));
        ARIADNE_TEST_ASSERT(canvas.primitives>count);
    }
    Void test_hybrid_collections() {
        HybridApproximatePoint start(DiscreteLocation(1),space,p),end(DiscreteLocation(1),space,next);
        HybridApproximatePoint jump(DiscreteLocation(2),space,p);
        Orbit<HybridApproximatePoint> single(start);
        single.insert(HybridTime(1.0_x,0),end);
        single.insert(HybridTime(1.0_x,1),jump);
        ARIADNE_TEST_EQUALS(single.size(),2u);
        check_endpoint(single.curve(0u),2u);
        ARIADNE_TEST_EQUALS(single.curves()[1].location(),DiscreteLocation(2));
        ARIADNE_TEST_FAIL(single.insert(HybridTime(1.0_x,-1),end));
        ARIADNE_TEST_FAIL(single.insert(HybridTime(1.0_x,3),end));
        check_filter(single);
        Orbit<HybridApproximatePoint> copied(single.curves());
        ARIADNE_TEST_EQUALS(copied.size(),2u);
        Orbit<HybridApproximatePoint> empty(List<HybridInterpolatedCurve>{});
        ARIADNE_TEST_EQUALS(empty.size(),0u);
        RecordingCanvas blank; empty.draw(blank,{},Variables2d(x,y));
        ARIADNE_TEST_EQUALS(blank.primitives,0u);
        std::ostringstream text; text<<single;
        ARIADNE_TEST_ASSERT(text.str().find("InterpolatedCurve")!=std::string::npos);
        Orbit<Vector<HybridApproximatePoint>> multi(Vector<HybridApproximatePoint>{start,end});
        multi.insert(HybridTime(1.0_x,0),end,0u);
        multi.insert(HybridTime(1.0_x,1),jump,0u);
        ARIADNE_TEST_EQUALS(multi.size(0u),2u);
        ARIADNE_TEST_EQUALS(multi.size(1u),1u);
        check_endpoint(multi.curve(0u,0u),2u);
        ARIADNE_TEST_EQUALS(multi.curves(0u)[1].location(),DiscreteLocation(2));
        ARIADNE_TEST_FAIL(multi.insert(HybridTime(1.0_x,-1),end,0u));
        ARIADNE_TEST_FAIL(multi.insert(HybridTime(1.0_x,3),end,0u));
        check_filter(multi);
        Orbit<Vector<HybridApproximatePoint>> no_points(Vector<HybridApproximatePoint>{});
        no_points.draw(blank,{},Variables2d(x,y));
        ARIADNE_TEST_EQUALS(blank.primitives,0u);
    }
    template<class ES> Void check_adjoin(ES const& enclosure) {
        Orbit<ES> orbit(enclosure);
        ARIADNE_TEST_EQUALS(orbit.initial().dimension(),2u);
        ARIADNE_TEST_EQUALS(orbit.reach().size(),0u);
        ARIADNE_TEST_EQUALS(orbit.intermediate().size(),0u);
        ARIADNE_TEST_EQUALS(orbit.final().size(),0u);
        orbit.adjoin_reach(enclosure); orbit.adjoin_intermediate(enclosure); orbit.adjoin_final(enclosure);
        ListSet<ES> sets; sets.adjoin(enclosure);
        orbit.adjoin_reach(sets); orbit.adjoin_intermediate(sets); orbit.adjoin_final(sets);
        ARIADNE_TEST_EQUALS(orbit.reach().size(),2u);
        ARIADNE_TEST_EQUALS(orbit.intermediate().size(),2u);
        ARIADNE_TEST_EQUALS(orbit.final().size(),2u);
        std::ostringstream description; description<<orbit;
        ARIADNE_TEST_ASSERT(description.str().find("intermediate=")!=std::string::npos);
    }
    Void test_enclosure_collections() {
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8)));
        Enclosure enclosure(ExactBoxType({{0,1},{0,1}}),config);
        check_adjoin(enclosure);
        LabelledEnclosure labelled(enclosure,space); check_adjoin(labelled);
        Orbit<LabelledEnclosure> labelled_orbit(labelled); labelled_orbit.adjoin_reach(labelled);
        RecordingCanvas canvas; labelled_orbit.draw(canvas,Variables2d(x,y));
        ARIADNE_TEST_ASSERT(canvas.boundaries>0u);
        std::unique_ptr<LabelledDrawable2dInterface> clone(labelled_orbit.labelled_clone());
        RecordingCanvas copy; clone->draw(copy,Variables2d(x,y));
        ARIADNE_TEST_EQUALS(copy.boundaries,canvas.boundaries);
        HybridEnclosure hybrid(DiscreteLocation(1),labelled); check_adjoin(hybrid);
        Orbit<HybridEnclosure> hybrid_orbit(hybrid); hybrid_orbit.adjoin_reach(hybrid);
        RecordingCanvas hybrid_canvas; hybrid_orbit.draw(hybrid_canvas,{},Variables2d(x,y));
        ARIADNE_TEST_ASSERT(hybrid_canvas.boundaries>0u);
        Orbit<HybridEnclosure> no_initial;
        ARIADNE_TEST_EQUALS(no_initial.reach().size(),0u);
        Figure figure(Box<ApproximateIntervalType>({{-1,2},{-1,2}}),Projection2d(2,0,1));
        Orbit<Enclosure> drawable(enclosure); drawable.adjoin_reach(enclosure); drawable.adjoin_final(enclosure);
        figure<<drawable;
        auto recorded=std::make_shared<RecordingCanvas>();
        BackendScope backend{RenderingBackend(RecordingBackend(recorded))};
        figure.write("orbit-memory");
        ARIADNE_TEST_ASSERT(recorded->boundaries>0u);
    }
};
int main() {
    redirect_test_logging_to_null(); TestOrbitCollections().test(); return ARIADNE_TEST_FAILURES; }
