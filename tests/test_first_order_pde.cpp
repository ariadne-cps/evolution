/***************************************************************************
 *            test_first_order_pde.cpp
 *
 *  Copyright  2008-21  Mirko Albanese
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

/** @file tests/test_first_order_pde.cpp
 *  @brief Test first order pde behavior and numerical invariants.
 */

#include <cmath>
#include "concise_logging.hpp"
#include <vector>

#include "algebra/tensor.hpp"
#include "evolution/first_order_pde.hpp"
#include "function/function.hpp"
#include "numeric/float_bounds.hpp"
#include "numeric/numeric.hpp"
#include "system/first_order_pde.hpp"

#include "ariadne-evolution.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;
using namespace std;

class TestFirstOrderPDE : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;
    void test_double_precision_contracts() { test_constant_and_contracts(double_precision); }
    void test_double_precision_evolution() { test_first_order_pde(double_precision); }
    void test_multiple_precision_contracts() { test_constant_and_contracts(MultiplePrecision(80u)); }

    template<class PR>
    void test_constant_and_contracts(PR pr) {
        Matrix<Real> identity = Matrix<Real>::identity(1u);
        DiagonalMatrix<Real> speed(Array<Real>{Real(1)});
        FirstOrderPDE pde = FirstOrderPDE::homogeneous(identity, 2u);
        pde.add_direction(identity, speed, identity).add_direction(identity, speed, identity);
        auto constant = EffectiveVectorMultivariateFunction::constant(2u, Vector<Real>{Real(2)});
        auto solution = first_order_pde(pde, constant, pr);
        ARIADNE_TESTING_PROVE(solution.uts.size(0) >= 4u);
        ARIADNE_TESTING_PROVE(solution.uts.size(2) > 1u);
        ARIADNE_TESTING_EXPECT_EQUAL(solution.tau.get_d() * static_cast<double>(solution.uts.size(2) - 1u), 0.25);
        for (SizeType i = 0; i < solution.uts.size(0); ++i) {
            for (SizeType j = 0; j < solution.uts.size(1); ++j) {
                auto value = solution.uts[{i, j, solution.uts.size(2) - 1u}][0];
                ARIADNE_TESTING_PROVE(value.lower().get_d() <= 2.0 && value.upper().get_d() >= 2.0);
                ARIADNE_TESTING_PROVE(value.upper().get_d() - value.lower().get_d() < 1e-10);
            }
        }
        FirstOrderPDE reverse = FirstOrderPDE::homogeneous(identity, 2u);
        DiagonalMatrix<Real> negative_speed(Array<Real>{Real(-1)});
        reverse.add_direction(-identity, negative_speed, identity).add_direction(-identity, negative_speed, identity);
        auto reversed = first_order_pde(reverse, constant, pr);
        auto endpoint = reversed.uts[{0u, 0u, reversed.uts.size(2) - 1u}][0];
        ARIADNE_TESTING_PROVE(endpoint.lower().get_d() <= 2.0 && endpoint.upper().get_d() >= 2.0);
        ARIADNE_TESTING_PROVE(endpoint.upper().get_d() - endpoint.lower().get_d() < 1e-10);
        // At Courant number one, upwind transport is an exact cyclic cell shift.
        auto coordinate = EffectiveVectorMultivariateFunction::identity(2u);
        EffectiveVectorMultivariateFunction ramp = {coordinate[0]};
        Matrix<Real> zero_matrix = 0 * identity;
        DiagonalMatrix<Real> zero_speed(Array<Real>{Real(0)});
        for (int direction: {-1, 1}) {
            FirstOrderPDE transport = FirstOrderPDE::homogeneous(identity, 2u);
            auto velocity = direction * identity;
            DiagonalMatrix<Real> characteristic(Array<Real>{Real(direction)});
            transport.add_direction(velocity, characteristic, identity).add_direction(zero_matrix, zero_speed, identity);
            auto translated = first_order_pde(transport, ramp, pr);
            SizeType cells = translated.uts.size(0), steps = translated.uts.size(2) - 1u;
            bool exact_shift = true;
            for (SizeType i = 0; i < cells; ++i) {
                auto source = direction > 0 ? (i + cells - steps % cells) % cells : (i + steps) % cells;
                double expected = (static_cast<double>(source) + 0.5) * translated.h.get_d();
                for (SizeType j = 0; j < translated.uts.size(1); ++j) {
                    auto value = translated.uts[{i, j, steps}][0];
                    exact_shift = exact_shift && value.lower().get_d() <= expected + 1e-12 && value.upper().get_d() >= expected - 1e-12 && value.upper().get_d() - value.lower().get_d() < 1e-10;
                }
            }
            ARIADNE_TESTING_PROVE(exact_shift);
        }
        // Valid dimensions can still request a grid too large to represent.
        auto coordinates = EffectiveVectorMultivariateFunction::identity(2u);
        EffectiveVectorMultivariateFunction huge = {Real(1e60_dec) * sqr(coordinates[0])};
        EffectiveVectorMultivariateFunction overflowing = {Real(1e308_dec) * sqr(coordinates[0])};
        ARIADNE_TESTING_EXPECT_EXCEPTION(first_order_pde(pde, huge, pr));
        ARIADNE_TESTING_EXPECT_EXCEPTION(first_order_pde(pde, overflowing, pr));
        ARIADNE_TESTING_EXPECT_EXCEPTION(first_order_pde(pde, EffectiveVectorMultivariateFunction::zeros(1u, 1u), pr));
        ARIADNE_TESTING_EXPECT_EXCEPTION(first_order_pde(pde, EffectiveVectorMultivariateFunction::zeros(2u, 2u), pr));
        FirstOrderPDE one_dimension = FirstOrderPDE::homogeneous(identity, 1u);
        one_dimension.add_direction(identity, speed, identity);
        ARIADNE_TESTING_EXPECT_EXCEPTION(first_order_pde(one_dimension, EffectiveVectorMultivariateFunction::zeros(1u, 1u), pr));
    }

    template<class PR>
    void test_first_order_pde(PR pr) {
        // DimensionType m=2;
        DimensionType n = 3;

        Real rho = 1 / 2_q;
        rho = 1 / 2_q;
        Real c = 1;
        c = 1;

        Matrix<Real> I = Matrix<Real>::identity(n);

        Matrix<Real> A = {
            {rho,   0, 0},
            {  0, rho, 0},
            {  0,   0, 1}
        };
        Matrix<Real> B0 = {
            {          0, 0, 1},
            {          0, 0, 0},
            {rho * c * c, 0, 0}
        };
        Matrix<Real> B1 = {
            {0,           0, 0},
            {0,           0, 1},
            {0, rho * c * c, 0}
        };
        ;

        Array<Matrix<Real>> Bs = {B0, B1};

        auto D0 = DiagonalMatrix<Real>(Array<Real>{c, 0, -c});
        auto D1 = DiagonalMatrix<Real>(Array<Real>{c, 0, -c});
        Array<DiagonalMatrix<Real>> Ds = {D0, D1};

        Real det = sqrt(1 + sqr(rho) * sqr(c));
        Real c1 = 1 / det;
        Real cr = (rho * c) / det;
        Matrix<Real> T0 = {
            {c1, 0,  c1},
            { 0, 1,   0},
            {cr, 0, -cr}
        };
        Matrix<Real> T1 = {
            { 0, 1,   0},
            {c1, 0,  c1},
            {cr, 0, -cr}
        };
        Array<Matrix<Real>> Ts = {T0, T1};

        Vector<Real> v00 = column(T0, 0);
        Vector<Real> v01 = column(T0, 1);
        Vector<Real> v02 = column(T0, 2);
        Vector<Real> v10 = column(T1, 0);
        Vector<Real> v11 = column(T1, 1);
        Vector<Real> v12 = column(T1, 2);

        auto z = EffectiveScalarMultivariateFunction::zero(EuclideanDomain(2));
        auto x = EffectiveVectorMultivariateFunction::identity(EuclideanDomain(2));

        EffectiveVectorMultivariateFunction f = EffectiveVectorMultivariateFunction::zeros(3, EuclideanDomain(2 + 1));
        EffectiveVectorMultivariateFunction phi0 = {sin(4 * x[0]), sin(6 * x[1]), z};

        FirstOrderPDE pde(A, f);
        pde.add_direction(B0, D0, T0).add_direction(B1, D1, T1);
        pde.validate();

        auto solution = first_order_pde(pde, phi0, pr);
        auto uts = solution.uts;

        // The periodic finite-volume update preserves each component's spatial mean.
        bool finite = true, initial_values_correct = true, conserved = true;
        std::vector<double> expected_sum(n, 0.0);
        for (SizeType i = 0; i < uts.size(0); ++i) {
            for (SizeType j = 0; j < uts.size(1); ++j) {
                double expected[3] = {std::sin(4 * (static_cast<double>(i) + 0.5) * solution.h.get_d()),
                                      std::sin(6 * (static_cast<double>(j) + 0.5) * solution.h.get_d()), 0.0};
                for (SizeType component = 0; component < n; ++component) {
                    auto value = uts[{i, j, 0u}][component];
                    expected_sum[component] += expected[component];
                    initial_values_correct = initial_values_correct && value.lower().get_d() <= expected[component] + 1e-12 && value.upper().get_d() >= expected[component] - 1e-12 && value.upper().get_d() - value.lower().get_d() < 1e-10;
                }
            }
        }
        for (SizeType step = 0; step < uts.size(2); ++step) {
            std::vector<double> lower(n, 0.0), upper(n, 0.0);
            for (SizeType i = 0; i < uts.size(0); ++i) {
                for (SizeType j = 0; j < uts.size(1); ++j) {
                    for (SizeType component = 0; component < n; ++component) {
                        auto value = uts[{i, j, step}][component];
                        finite = finite && std::isfinite(value.lower().get_d()) && std::isfinite(value.upper().get_d());
                        lower[component] += value.lower().get_d();
                        upper[component] += value.upper().get_d();
                    }
                }
            }
            for (SizeType component = 0; component < n; ++component) {
                conserved = conserved && lower[component] <= expected_sum[component] + 1e-8 && upper[component] >= expected_sum[component] - 1e-8;
            }
        }
        ARIADNE_TESTING_PROVE(finite);
        ARIADNE_TESTING_PROVE(initial_values_correct);
        ARIADNE_TESTING_PROVE(conserved);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_first_order_pde");
    ARIADNE_TESTING_ADD_TEST(TestFirstOrderPDE, test_double_precision_contracts);
    ARIADNE_TESTING_ADD_TEST(TestFirstOrderPDE, test_double_precision_evolution);
    ARIADNE_TESTING_ADD_TEST(TestFirstOrderPDE, test_multiple_precision_contracts);
    return run_evolution_test_suite(ariadne_test_suite, argc, argv);
}
