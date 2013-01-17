/** @file
 * Declaration of the inversion algorithms
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
 * @param[in] params The inputparameters of the application
 *
 * TODO make gaugefield a const-ref
 */
void perform_inversion(const std::vector<physics::lattices::Spinorfield*> * result, physics::lattices::Gaugefield* gaugefield, const std::vector<physics::lattices::Spinorfield*>& sources, const hardware::System& system);

}

}
