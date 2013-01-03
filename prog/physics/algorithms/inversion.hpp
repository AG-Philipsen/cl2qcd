/** @file
 * Declaration of the flavour doubles algorithms
 */

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"

namespace physics {

namespace algorithms {

/**
 * Perform the inversion and store result to solution_buffer
 *
 * @param[out] result Spinorfield in which to store the inversion result
 * @param[in] gaugefield Gaugefield on which to base the inversion
 * @param[in] sources Spinorfields from which to start the inversion
 */
void perform_inversion(const std::vector<const physics::lattices::Spinorfield*> * result, const physics::lattices::Gaugefield& gaugefield, const std::vector<const physics::lattices::Spinorfield*>& sources);

}

}
