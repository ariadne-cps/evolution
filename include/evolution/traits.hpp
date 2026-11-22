/***************************************************************************
 *            traits.hpp
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

/** @file evolution/traits.hpp
 *  @brief Compile-time correspondence between systems, evolvers and enclosures.
 */

#ifndef ARIADNE_EVOLUTION_TRAITS_HPP
#define ARIADNE_EVOLUTION_TRAITS_HPP

// Associations between systems and their evolution algorithms belong to evolution,
// never to the system-description classes.
#include "evolution/evolver.decl.hpp"

namespace Ariadne {

/** @brief IteratedMap: compile-time correspondence between systems, evolvers and enclosures. */
class IteratedMap;
/** @brief OrdinaryDifferentialAlgebraicEquation: compile-time correspondence between systems, evolvers and enclosures. */
class OrdinaryDifferentialAlgebraicEquation;
/** @brief OrdinaryDifferentialAlgebraicInclusion: compile-time correspondence between systems, evolvers and enclosures. */
class OrdinaryDifferentialAlgebraicInclusion;
/** @brief HybridAutomatonInterface: compile-time correspondence between systems, evolvers and enclosures. */
class HybridAutomatonInterface;

/** @brief DiscreteTimeEvolver<ValidatedTag>: compile-time correspondence between systems, evolvers and enclosures. */
template<>
class DiscreteTimeEvolver<ValidatedTag>;
/** @brief ContinuousTimeEvolver<ValidatedTag>: compile-time correspondence between systems, evolvers and enclosures. */
template<>
class ContinuousTimeEvolver<ValidatedTag>;
class ContinuousDynamics;
/** @brief HybridContinuousTimeEvolverInterface<ValidatedTag>: compile-time correspondence between systems, evolvers and enclosures. */
template<>
class HybridContinuousTimeEvolverInterface<ValidatedTag>;

/** @brief LabelledEnclosure: compile-time correspondence between systems, evolvers and enclosures. */
class LabelledEnclosure;
/** @brief HybridEnclosure: compile-time correspondence between systems, evolvers and enclosures. */
class HybridEnclosure;

template<class P, class SYS>
struct EvolutionTraits;

template<>
struct EvolutionTraits<ValidatedTag, IteratedMap> {
    using EvolverType = DiscreteTimeEvolver<ValidatedTag>;
    using EnclosureType = LabelledEnclosure;
};

template<>
struct EvolutionTraits<ValidatedTag, OrdinaryDifferentialAlgebraicEquation> {
    using EvolverType = ContinuousTimeEvolver<ValidatedTag>;
    using EnclosureType = LabelledEnclosure;
};

template<>
struct EvolutionTraits<ValidatedTag, OrdinaryDifferentialAlgebraicInclusion> {
    using EvolverType = ContinuousTimeEvolver<ValidatedTag>;
    using EnclosureType = LabelledEnclosure;
};

template<>
struct EvolutionTraits<ValidatedTag, ContinuousDynamics> {
    using EvolverType = ContinuousTimeEvolver<ValidatedTag>;
    using EnclosureType = LabelledEnclosure;
};

template<>
struct EvolutionTraits<ValidatedTag, HybridAutomatonInterface> {
    using EvolverType = HybridContinuousTimeEvolverInterface<ValidatedTag>;
    using EnclosureType = HybridEnclosure;
};

template<>
struct EvolutionTraits<ApproximateTag, OrdinaryDifferentialAlgebraicEquation> {
    using EvolverType = ContinuousTimeEvolver<ApproximateTag>;
};

template<>
struct EvolutionTraits<ApproximateTag, HybridAutomatonInterface> {
    using EvolverType = HybridContinuousTimeEvolverInterface<ApproximateTag>;
};

} // namespace Ariadne
#endif // ARIADNE_EVOLUTION_TRAITS_HPP
