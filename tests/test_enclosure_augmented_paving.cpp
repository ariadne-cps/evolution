/***************************************************************************
 *            test_enclosure_augmented_paving.cpp
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

/** @file tests/test_enclosure_augmented_paving.cpp
 *  @brief Test enclosure augmented paving behavior and numerical invariants.
 */

// Verify the continuous and hybrid enclosure-to-augmented-paving adapters.
#include "evolution/enclosure_augmented_paving.hpp"
#include "function/taylor_function.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class TestEnclosureAugmentedPaving : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    Void test_continuous() {
        ExactBoxType domain({{-1,1},{-1,1}});
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8)));
        Enclosure enclosure(domain,config);
        enclosure.set_auxiliary_mapping(EffectiveVectorMultivariateFunction::identity(2u));

        Grid grid(2u);
        AugmentedPaving generated=outer_approximation(enclosure,grid,2u);
        ARIADNE_TESTING_PROVE(generated.size()>0u);
        ARIADNE_TESTING_EXPECT_EQUAL(generated.dimension(),2u);
        ARIADNE_TESTING_EXPECT_EQUAL(generated.state_auxiliary_paving().dimension(),4u);

        AugmentedPaving explicit_paving(grid,enclosure.auxiliary_mapping());
        adjoin_outer_approximation_to(enclosure,explicit_paving,2u);
        ARIADNE_TESTING_PROVE(explicit_paving.size()>0u);
        ARIADNE_TESTING_PROVE(subset(generated.state_paving(),explicit_paving.state_paving()));
        ARIADNE_TESTING_PROVE(subset(explicit_paving.state_paving(),generated.state_paving()));
        AugmentedPaving independent(grid,EffectiveVectorMultivariateFunction::identity(2u));
        adjoin_outer_approximation_to(enclosure,independent,2u);
        ARIADNE_TESTING_PROVE(subset(independent.state_paving(),generated.state_paving()));
        Enclosure without_auxiliary(domain,config);
        auto state_only=outer_approximation(without_auxiliary,grid,2u);
        ARIADNE_TESTING_PROVE(subset(state_only.state_paving(),generated.state_paving()));
    }

    Void test_hybrid() {
        DiscreteLocation location(1);
        RealVariable x("x"), y("y"), a("a");
        RealSpace state_space({x,y});
        RealSpace auxiliary_space({a});
        ExactBoxType domain({{-1,1},{-1,1}});
        EnclosureConfiguration config(TaylorFunctionFactory(ThresholdSweeper<FloatDP>(dp,1e-8)));
        HybridEnclosure enclosure(HybridExactBox(location,state_space,domain),config);

        MonolithicHybridSpace state;
        state.new_location(location,state_space);
        MonolithicHybridSpace extended;
        extended.new_location(location,join(state_space,auxiliary_space));
        LabelledMapping labelled(state_space,
            EffectiveVectorMultivariateFunction({EffectiveScalarMultivariateFunction::coordinate(2u,0u)}),
            auxiliary_space);
        HybridAuxiliaryMapping mapping(HybridSpace(extended),
            [labelled](DiscreteLocation const&) { return labelled; });
        HybridAugmentedPaving paving(HybridGrid(state),mapping);
        adjoin_outer_approximation_to(enclosure,paving,2u);
        ARIADNE_TESTING_PROVE(paving.size()>0u);
        ARIADNE_TESTING_EXPECT_EQUAL(paving.auxiliary_mapping(location).argument_space(),state_space);
        ARIADNE_TESTING_EXPECT_EQUAL(paving.state_auxiliary_paving().space(location),join(state_space,auxiliary_space));

        // The same spatial precondition as the underlying HybridEnclosure API.
        MonolithicHybridSpace wrong_state;
        wrong_state.new_location(location,auxiliary_space);
        HybridAugmentedPaving wrong(HybridGrid(wrong_state),mapping);
        ARIADNE_TESTING_EXPECT_EXCEPTION(adjoin_outer_approximation_to(enclosure,wrong,2u));
    }
};

int main(int argc,char** argv) {
    ARIADNE_TESTING_SUITE("test_enclosure_augmented_paving");
    ARIADNE_TESTING_ADD_TEST(TestEnclosureAugmentedPaving,test_continuous);
    ARIADNE_TESTING_ADD_TEST(TestEnclosureAugmentedPaving,test_hybrid);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc,argv);
}
