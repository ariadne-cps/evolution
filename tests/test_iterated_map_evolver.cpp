/***************************************************************************
 *            test_iterated_map_evolver.cpp
 *
 *  Copyright  2006-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>

#include "utility/tuple.hpp"
#include "numeric/rational.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/algebra.hpp"
#include "function/taylor_model.hpp"
#include "algebra/differential.hpp"
#include "function/constraint.hpp"
#include "function/function.hpp"
#include "function/taylor_function.hpp"
#include "function/formula.hpp"
#include "evolution/enclosure.hpp"
#include "evolution/orbit.hpp"
#include "geometry/box.hpp"
#include "set/list_set.hpp"
#include "system/iterated_map.hpp"
#include "evolution/iterated_map_evolver.hpp"
#include "labelled/figure.hpp"
#include "drawing/drawing_command_line_interface.hpp"

#include "macro/testing.hpp"
#include "log_verbosity_scope.hpp"
#include <memory>

using namespace Ariadne;
using namespace std;

class TestIteratedMapEvolver
{
  public:
    Void test() const {
        ARIADNE_TEST_CALL(test_stopping_policies());
        ARIADNE_TEST_CALL(test_evolve());
        ARIADNE_TEST_CALL(test_subdivide());
    }

    Void test_stopping_policies() const {
        RealVariable x("x"),other("other");
        IteratedMap map(List<PrimedRealAssignment>{next(x)=x+1});
        IteratedMapEvolver evolver(map);
        auto initial=evolver.enclosure(ExactBoxType({{-1.0_x,1.0_x}}));
        std::unique_ptr<IteratedMapEvolver> clone(evolver.clone());
        auto const& const_evolver=*clone;
        ARIADNE_TEST_EQUALS(const_evolver.configuration().maximum_enclosure_radius(),100.0);
        ARIADNE_TEST_EQUALS(clone->system().state_space(),RealSpace({x}));
        auto zero=evolver.orbit(initial,Integer(0),Semantics::LOWER);
        ARIADNE_TEST_EQUALS(zero.final().size(),1u);
        ARIADNE_TEST_EQUALS(zero.reach().size(),0u);
        ARIADNE_TEST_EQUALS(zero.intermediate().size(),0u);
        evolver.configuration().set_maximum_enclosure_radius(0.5);
        evolver.configuration().set_enable_premature_termination(true);
        auto stopped=evolver.orbit(initial,Integer(1),Semantics::LOWER);
        ARIADNE_TEST_EQUALS(stopped.final().size(),0u);
        ARIADNE_TEST_EQUALS(stopped.reach().size(),0u);
        evolver.configuration().set_maximum_enclosure_radius(2.0);
        auto small=evolver.orbit(initial,Integer(1),Semantics::LOWER);
        ARIADNE_TEST_EQUALS(small.final().size(),1u);
        evolver.configuration().set_maximum_enclosure_radius(0.5);
        evolver.configuration().set_enable_premature_termination(false);
        {
            LogVerbosityScope diagnostics(100u);
            auto continued=evolver.orbit(initial,Integer(2),Semantics::LOWER);
            ARIADNE_TEST_EQUALS(continued.final().size(),1u);
            ARIADNE_TEST_EQUALS(continued.reach().size(),2u);
            ARIADNE_TEST_EQUALS(continued.intermediate().size(),2u);
            ARIADNE_TEST_ASSERT(definitely(continued.final()[0].euclidean_set().subset(ExactBoxType({{0.99_x,3.01_x}}))));
        }
        auto unsplit=evolver.orbit(initial,Integer(1),Semantics::UPPER);
        ARIADNE_TEST_EQUALS(unsplit.final().size(),1u);
        auto incompatible=LabelledEnclosure(initial.euclidean_set(),RealSpace({other}));
        ARIADNE_TEST_FAIL(evolver.orbit(incompatible,Integer(1),Semantics::UPPER));
    }

    Void test_evolve() const
    {
        typedef DoublePrecision PR;
        DoublePrecision pr;

        // Set up the map field
        // The Henon map \f$(x,y)\mapsto(a-x^2+by,x)
        Real a=1.3_dec; Real b=0.3_dec;
        RealVariable x("x"), y("y");
        IteratedMap henon({ next(x)=a-x*x+b*y, next(y)=x });
        ARIADNE_TEST_PRINT(henon);

        // Define the initial set
        RealExpressionBoundedConstraintSet initial_set({0.4_dec<=x<=0.7_dec,0.9_dec<=y<=1.1_dec});

        ARIADNE_TEST_PRINT(initial_set);

        // Function evaluation sanity check
        Vector<FloatApproximation<PR>> p({a,b},pr);
        Vector<FloatApproximation<PR>> xa({0.5_x,0.25_x},pr);
        Vector<FloatApproximation<PR>> hxa={p[0]-xa[0]*xa[0]+xa[1]*p[1], xa[0]};
        ARIADNE_TEST_EQUAL(henon.update_function().evaluate(xa),hxa);
        Matrix<FloatApproximation<PR>> dhxa={{-2*xa[0],p[1]},{{1.0,pr},{0.0,pr}}};
        ARIADNE_TEST_EQUAL(henon.update_function().jacobian(xa),dhxa);

        // Set up the evolution parameters and grid
        IteratedMap::TimeType time(6);

        // Set up the evaluators
        IteratedMapEvolver evolver(henon);
        evolver.configuration().set_maximum_enclosure_radius(0.25);

        ARIADNE_TEST_PRINT(evolver.configuration());

        // Compute the orbit
        auto orbit = evolver.orbit(initial_set,time,Semantics::UPPER);

        auto initial_enclosure = LabelledEnclosure(cast_exact_box(initial_set.euclidean_set(henon.state_space()).bounding_box()),
                                                   henon.state_space(),EnclosureConfiguration(evolver.function_factory()));

        // Print the initial, evolve and reach sets
        LabelledFigure fig(Axes2d(-10<=x<=10, -10<=y<=10));
        fig << line_style(true) << fill_colour(cyan) << orbit.reach();
        fig << fill_colour(yellow) << orbit.final();
        fig << fill_colour(blue) << initial_enclosure;
        fig.write("test_iterated_map_evolver_xy");

        LabelledFigure fig2(Axes2d(0<=TimeVariable()<=10, -10<=y<=10));
        fig2 << orbit.reach();
        fig2.write("test_iterated_map_evolver_ty");
    }

    Void test_subdivide() const
    {
        Real a=1.3_dec; Real b=0.3_dec;
        RealVariable x("x"), y("y");
        IteratedMap henon({ next(x)=a-x*x+b*y, next(y)=x });

        RealExpressionBoundedConstraintSet initial_set({0.4_dec<=x<=0.7_dec,0.9_dec<=y<=1.1_dec});

        IteratedMap::TimeType time(6);

        IteratedMapEvolver evolver(henon);
        evolver.configuration().set_maximum_enclosure_radius(0.1);
        evolver.configuration().set_enable_subdivisions(true);

        auto orbit = evolver.orbit(initial_set,time,Semantics::UPPER);

        ARIADNE_TEST_ASSERT(orbit.final().size()>1);

        LabelledFigure fig(Axes2d(-10<=x<=10, -10<=y<=10));
        fig << line_style(true) << fill_colour(cyan) << orbit.reach();
        fig << fill_colour(yellow) << orbit.final();
        fig.write("test_iterated_map_evolver_xy_subdivisions");
    }
};

Int main(Int argc, const char* argv[])
{
    if (not DrawingCommandLineInterface().acquire(argc,argv)) return -1;

    TestIteratedMapEvolver().test();
    return ARIADNE_TEST_FAILURES;
}
