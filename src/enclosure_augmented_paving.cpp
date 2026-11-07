// Adapt evolved enclosures into drawing-independent augmented pavings.
#include "function/functional.hpp"
#include "function/constraint.hpp"
#include "set/function_set.hpp"
#include "evolution/enclosure_augmented_paving.hpp"
#include "hybrid/hybrid_augmented_paving.hpp"
#include "evolution/hybrid_enclosure.hpp"
#include "paving/binary_tree_paving.hpp"
#include "symbolic/space.hpp"
#include "labelled/figure.hpp"
#include "logging/logging.hpp"

namespace Ariadne {

Void adjoin_outer_approximation_to(Enclosure const& enclosure, AugmentedPaving& paving, Nat fineness)
{
    if (enclosure.auxiliary_mapping().result_size()!=0 &&
            enclosure.auxiliary_mapping().managed_pointer()==paving.auxiliary_mapping().managed_pointer())
    {
        ARIADNE_WARN("enclosure="<<enclosure<<", "
                    "enclosure.auxiliary_mapping()="<<enclosure.auxiliary_mapping()<<", "
                    "paving.auxiliary_mapping()="<<paving.auxiliary_mapping());
    }
    enclosure.configuration().paver().adjoin_outer_approximation(paving.state_paving(),enclosure.state_set(),fineness);
}


AugmentedPaving outer_approximation(Enclosure const& enclosure, Grid const& grid, Nat fineness)
{
    AugmentedPaving result(grid,enclosure.auxiliary_mapping());
    adjoin_outer_approximation_to(enclosure,result,fineness);
    return result;
}

Void adjoin_outer_approximation_to(HybridEnclosure const& enclosure, HybridAugmentedPaving& paving, Nat fineness)
{
    enclosure.adjoin_outer_approximation_to(paving.state_paving(),fineness);
}

} // namespace Ariadne
