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
namespace Ariadne {

/** @brief IteratedMap: compile-time correspondence between systems, evolvers and enclosures. */
class IteratedMap;
/** @brief VectorField: compile-time correspondence between systems, evolvers and enclosures. */
class VectorField;
/** @brief DifferentialInclusion: compile-time correspondence between systems, evolvers and enclosures. */
class DifferentialInclusion;
/** @brief HybridAutomatonInterface: compile-time correspondence between systems, evolvers and enclosures. */
class HybridAutomatonInterface;

/** @brief IteratedMapEvolver: compile-time correspondence between systems, evolvers and enclosures. */
class IteratedMapEvolver;
/** @brief VectorFieldEvolver: compile-time correspondence between systems, evolvers and enclosures. */
class VectorFieldEvolver;
/** @brief DifferentialInclusionEvolver: compile-time correspondence between systems, evolvers and enclosures. */
class DifferentialInclusionEvolver;
/** @brief HybridEvolverInterface: compile-time correspondence between systems, evolvers and enclosures. */
class HybridEvolverInterface;

/** @brief LabelledEnclosure: compile-time correspondence between systems, evolvers and enclosures. */
class LabelledEnclosure;
/** @brief HybridEnclosure: compile-time correspondence between systems, evolvers and enclosures. */
class HybridEnclosure;

template<class SYS> struct EvolutionTraits;

template<> struct EvolutionTraits<IteratedMap> {
    using EvolverType = IteratedMapEvolver;
    using EnclosureType = LabelledEnclosure;
};

template<> struct EvolutionTraits<VectorField> {
    using EvolverType = VectorFieldEvolver;
    using EnclosureType = LabelledEnclosure;
};

template<> struct EvolutionTraits<DifferentialInclusion> {
    using EvolverType = DifferentialInclusionEvolver;
    using EnclosureType = LabelledEnclosure;
};

template<> struct EvolutionTraits<HybridAutomatonInterface> {
    using EvolverType = HybridEvolverInterface;
    using EnclosureType = HybridEnclosure;
};

} // namespace Ariadne
#endif // ARIADNE_EVOLUTION_TRAITS_HPP
