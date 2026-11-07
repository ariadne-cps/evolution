/********************************************************************************
 *            evolution/enclosure_augmented_paving.hpp
 *
 * Outer-paving adapters for continuous and hybrid evolution enclosures.
 * These definitions belong to evolution, not to the reasoning analyser.
 ********************************************************************************/
#ifndef ARIADNE_EVOLUTION_ENCLOSURE_AUGMENTED_PAVING_HPP
#define ARIADNE_EVOLUTION_ENCLOSURE_AUGMENTED_PAVING_HPP

#include "evolution/enclosure.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "paving/augmented_paving.hpp"
#include "hybrid/hybrid_augmented_paving.hpp"

namespace Ariadne {

// Add a validated outer approximation of an enclosure to an augmented
// state paving. The enclosure's configured paver controls approximation.
Void adjoin_outer_approximation_to(Enclosure const& enclosure, AugmentedPaving& paving, Nat fineness);

// Construct an augmented paving with the enclosure's auxiliary mapping
// and a validated outer approximation of its state set.
AugmentedPaving outer_approximation(Enclosure const& enclosure, Grid const& grid, Nat fineness);

// Reuse HybridEnclosure's location-aware outer approximation, without
// requiring access to an automaton or other reasoning-layer objects.
Void adjoin_outer_approximation_to(HybridEnclosure const& enclosure, HybridAugmentedPaving& paving, Nat fineness);

} // namespace Ariadne

#endif // ARIADNE_EVOLUTION_ENCLOSURE_AUGMENTED_PAVING_HPP
