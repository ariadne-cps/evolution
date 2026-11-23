/***************************************************************************
 *            test_hybrid_crossing.cpp
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

/** @file tests/test_hybrid_crossing.cpp
 *  @brief Verify hybrid crossing behavior.
 */

#include <type_traits>

#include "../src/hybrid_crossing_analysis.hpp"
#include "evolution/hybrid_crossing.hpp"
#include "evolution/hybrid_evolution_backend.hpp"
#include "integration/integrator.hpp"
#include "logging/test_runner.hpp"
#include "system/hybrid_automata.hpp"

#include "ariadne-testing.hpp"
#include "equation-solving/interval_newton_method.hpp"

using namespace Ariadne;
using namespace Ariadne::EvolutionDetail;

static_assert(!std::is_default_constructible_v<CrossingTimeCertificate>);
static_assert(!std::is_default_constructible_v<GuardMaximumCertificate>);

namespace {
// Exact flow x(t)=x0+t, y(t)=y0: roots and maxima are known independently
// of the interval solver used by crossing analysis.
struct AnalyticFlow {
    const RealVariable x{"x"}, y{"y"};
    const RealSpace space{x, y};
    const DiscreteEvent event{"guard"};
    TaylorPicardIntegrator integrator{StepMaximumError(1e-6)};
    ValidatedFunctionPatchFactory factory{integrator.function_factory()};
    HybridEnclosure initial;
    ValidatedVectorMultivariateFunctionPatch flow;
    EffectiveVectorMultivariateFunction dynamic{make_function(space, Vector<RealExpression>{1, 0})};

    explicit AnalyticFlow(const ExactBoxType& box)
        : initial(DiscreteLocation(1), space, Enclosure(box, EnclosureConfiguration(factory))),
          flow(factory.create_zeros(2u, product(box, ExactIntervalType(0, 2)))) {
        flow[0u] = factory.create_coordinate(flow.domain(), 0u) + factory.create_coordinate(flow.domain(), 2u);
        flow[1u] = factory.create_coordinate(flow.domain(), 1u);
    }

    Map<DiscreteEvent, CrossingData> analyze(const RealExpression& guard, const EquationSolverInterface& solver) const {
        const Map<DiscreteEvent, EffectiveScalarMultivariateFunction> guards{
            {event, make_function(space, guard)}
        };
        return analyze_crossings(Set<DiscreteEvent>{event}, dynamic, guards, flow, initial, solver);
    }
};

// Deterministic failure injection distinguishes a failed critical-time solve
// from a failed first-entry solve after the maximum has already been certified.
class FailingSolver : public IntervalNewtonSolver {
    const unsigned int _fail_on_call;
    const bool _unexpected;
    mutable unsigned int _calls = 0u;

  public:
    FailingSolver(unsigned int fail_on_call, bool unexpected = false)
        : IntervalNewtonSolver(1e-8, 32u), _fail_on_call(fail_on_call), _unexpected(unexpected) {}
    FailingSolver* clone() const override { return new FailingSolver(*this); }
    using IntervalNewtonSolver::implicit;
    ValidatedScalarMultivariateFunctionPatch implicit(
        const ValidatedScalarMultivariateFunction& function, const ExactBoxType& domain,
        const ExactIntervalType& time) const override {
        if (++_calls == _fail_on_call) {
            if (_unexpected)
                throw std::runtime_error("unexpected analytic solver failure");
            throw EquationSolverException("uncertified analytic root");
        }
        return IntervalNewtonSolver::implicit(function, domain, time);
    }
};
} // namespace

class TestHybridCrossing : public Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_certificate_invariants() {
        AnalyticFlow scenario(ExactBoxType({
            {-1, -0.75_x},
            { 0,       0}
        }));
        const auto time = -scenario.factory.create_coordinate(scenario.initial.parameter_domain(), 0u);
        const CrossingData unknown;
        ARIADNE_TESTING_EXPECT_EQUAL(unknown.geometry(), CrossingGeometry::UNRESOLVED);
        ARIADNE_TESTING_PROVE(unknown.crossing_time() == nullptr && unknown.guard_maximum() == nullptr);
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument, CrossingData(static_cast<CrossingGeometry>(255)));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      (void)CrossingTimeCertificate(ValidatedScalarMultivariateFunctionPatch{}));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      CrossingData(CrossingGeometry::DECREASING, CrossingTimeCertificate(time)));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      GuardMaximumCertificate(time, UpperIntervalType(-1, 1), UpperBoxType({
                                                                                                  {0, 0}
        })));
        ARIADNE_TESTING_EXPECT_THROWS(std::invalid_argument,
                                      GuardMaximumCertificate(time, UpperIntervalType(1, -1), UpperBoxType({
                                                                                                  {0, 0},
                                                                                                  {0, 0}
        })));
        auto owned_map = time;
        const CrossingData crossing(CrossingGeometry::INCREASING, CrossingTimeCertificate(owned_map));
        owned_map = owned_map + 10;
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(crossing.crossing_time()->time().range()),
                                     ExactIntervalType(0.75_x, 1)));
    }

    void test_increasing_root_and_decreasing_geometry() {
        AnalyticFlow scenario(ExactBoxType({
            {-1, -0.75_x},
            { 0,       0}
        }));
        IntervalNewtonSolver solver(1e-8, 32u);
        const auto increasing = scenario.analyze(scenario.x, solver)[scenario.event];
        ARIADNE_TESTING_EXPECT_EQUAL(increasing.geometry(), CrossingGeometry::INCREASING);
        ARIADNE_TESTING_PROVE(increasing.crossing_time() != nullptr);
        const auto root = increasing.crossing_time()->time();
        const auto residual = root + scenario.factory.create_coordinate(root.domain(), 0u);
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-30_x, 0x1p-30_x)));
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(0.75_x, 1), cast_exact_interval(root.range())));
        const auto decreasing = scenario.analyze(-scenario.x, solver)[scenario.event];
        ARIADNE_TESTING_EXPECT_EQUAL(decreasing.geometry(), CrossingGeometry::DECREASING);
        ARIADNE_TESTING_PROVE(decreasing.crossing_time() == nullptr && decreasing.guard_maximum() == nullptr);
    }

    void test_concave_maximum_and_exclusion() {
        IntervalNewtonSolver solver(1e-8, 32u);
        AnalyticFlow uncertain(ExactBoxType({
            {      -1, -0.75_x},
            {-0.125_x, 0.125_x}
        }));
        const auto crossing = uncertain.analyze(uncertain.y - sqr(uncertain.x), solver)[uncertain.event];
        ARIADNE_TESTING_EXPECT_EQUAL(crossing.geometry(), CrossingGeometry::CONCAVE);
        ARIADNE_TESTING_PROVE(crossing.crossing_time() == nullptr && crossing.guard_maximum() != nullptr);
        const auto& maximum = *crossing.guard_maximum();
        const auto residual = maximum.time() + uncertain.factory.create_coordinate(maximum.time().domain(), 0u);
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-30_x, 0x1p-30_x)));
        ARIADNE_TESTING_PROVE(subset(ExactIntervalType(-0.125_x, 0.125_x), cast_exact_interval(maximum.guard_range())));
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(maximum.state_bounds()[0u]), ExactIntervalType(-0x1p-30_x, 0x1p-30_x)));
        AnalyticFlow negative(ExactBoxType({
            {     -1,  -0.75_x},
            {-0.25_x, -0.125_x}
        }));
        ARIADNE_TESTING_PROVE(negative.analyze(negative.y - sqr(negative.x), solver).empty());
        const auto convex = uncertain.analyze(sqr(uncertain.x) + uncertain.y, solver)[uncertain.event];
        ARIADNE_TESTING_EXPECT_EQUAL(convex.geometry(), CrossingGeometry::CONVEX);
        const auto unresolved = uncertain.analyze(pow(uncertain.x, 3), solver)[uncertain.event];
        ARIADNE_TESTING_EXPECT_EQUAL(unresolved.geometry(), CrossingGeometry::UNRESOLVED);
    }

    void test_solver_fallbacks_preserve_available_evidence() {
        AnalyticFlow scenario(ExactBoxType({
            {    -1, -0.75_x},
            {0.25_x,  0.25_x}
        }));
        const auto increasing = scenario.analyze(scenario.x, FailingSolver(1u))[scenario.event];
        ARIADNE_TESTING_EXPECT_EQUAL(increasing.geometry(), CrossingGeometry::INCREASING);
        ARIADNE_TESTING_PROVE(increasing.crossing_time() == nullptr);
        ARIADNE_TESTING_EXPECT_THROWS(std::runtime_error, scenario.analyze(scenario.x, FailingSolver(1u, true)));
        for (bool unexpected: {false, true}) {
            const auto concave = scenario.analyze(scenario.y - sqr(scenario.x), FailingSolver(1u, unexpected))[scenario.event];
            ARIADNE_TESTING_EXPECT_EQUAL(concave.geometry(), CrossingGeometry::CONCAVE);
            ARIADNE_TESTING_PROVE(concave.crossing_time() == nullptr && concave.guard_maximum() == nullptr);
        }
        const auto unexpected_entry = scenario.analyze(scenario.y - sqr(scenario.x), FailingSolver(2u, true))[scenario.event];
        ARIADNE_TESTING_EXPECT_EQUAL(unexpected_entry.geometry(), CrossingGeometry::CONCAVE);
        ARIADNE_TESTING_PROVE(unexpected_entry.crossing_time() == nullptr && unexpected_entry.guard_maximum() == nullptr);
        const auto maximum = scenario.analyze(scenario.y - sqr(scenario.x), FailingSolver(2u))[scenario.event];
        ARIADNE_TESTING_EXPECT_EQUAL(maximum.geometry(), CrossingGeometry::CONCAVE);
        ARIADNE_TESTING_PROVE(maximum.guard_maximum() != nullptr);
        ARIADNE_TESTING_PROVE(definitely(maximum.guard_maximum()->guard_range().lower_bound() > 0));
        // A positive maximum after failed root solving is not a tangency.
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(maximum.guard_maximum()->guard_range()),
                                     ExactIntervalType(0.25_x, 0.25_x)));
    }

    void test_concave_first_entry_and_provenance() {
        AnalyticFlow scenario(ExactBoxType({
            {    -1, -0.984375_x},
            {0.25_x,      0.25_x}
        }));
        IntervalNewtonSolver solver(1e-8, 32u);
        const auto crossing = scenario.analyze(scenario.y - sqr(scenario.x), solver)[scenario.event];
        ARIADNE_TESTING_EXPECT_EQUAL(crossing.geometry(), CrossingGeometry::CONCAVE);
        ARIADNE_TESTING_PROVE(crossing.crossing_time() != nullptr);
        if (crossing.crossing_time() == nullptr)
            throw std::runtime_error("Expected concave first-entry certificate");
        const auto& time = crossing.crossing_time()->time();
        // First entry is -x0-1/2; the other root -x0+1/2 must be excluded.
        const auto residual = time + scenario.factory.create_coordinate(time.domain(), 0u) + 0.5_x;
        ARIADNE_TESTING_PROVE(subset(cast_exact_interval(residual.range()), ExactIntervalType(-0x1p-24_x, 0x1p-24_x)));
        ARIADNE_TESTING_PROVE(cast_exact_interval(time.range()).upper_bound() < FloatDP(0.75_x, dp));
        HybridEvolutionStepResult result;
        result.transition_to(scenario.initial, scenario.event, EventKind::URGENT, crossing.summary());
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.at(0u).crossing->geometry, CrossingGeometry::CONCAVE);
        ARIADNE_TESTING_EXPECT_EQUAL(result.successors.at(0u).crossing->evidence, CrossingEvidenceKind::CROSSING_TIME);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_hybrid_crossing");
    ARIADNE_TESTING_ADD_TEST(TestHybridCrossing, test_certificate_invariants);
    ARIADNE_TESTING_ADD_TEST(TestHybridCrossing, test_increasing_root_and_decreasing_geometry);
    ARIADNE_TESTING_ADD_TEST(TestHybridCrossing, test_concave_maximum_and_exclusion);
    ARIADNE_TESTING_ADD_TEST(TestHybridCrossing, test_solver_fallbacks_preserve_available_evidence);
    ARIADNE_TESTING_ADD_TEST(TestHybridCrossing, test_concave_first_entry_and_provenance);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
