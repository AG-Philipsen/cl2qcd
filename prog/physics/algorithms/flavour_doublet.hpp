/** @file
 * Declaration of the flavour doublets algorithms
 */

#ifndef _PHYSICS_ALGORITHMS_FLAVOUR_DOUBLET_
#define _PHYSICS_ALGORITHMS_FLAVOUR_DOUBLET_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"

namespace physics {

namespace algorithms {

/**
 * Calculate flavour multiplet correlators from private solution_buffer and store them to a file
 * @param[in] gaugefield The gaugefield
 * @param[in] result The Spinorfields
 * @param[in] corr_fn filename
 */
void flavour_doublet_correlators(const std::vector<physics::lattices::Spinorfield*>& result, const std::vector<physics::lattices::Spinorfield*>& sources, std::string corr_fn, const meta::Inputparameters& params);

std::vector<hmc_float> calculate_correlator(const std::string& type, const std::vector<physics::lattices::Spinorfield*>& corr, const std::vector<physics::lattices::Spinorfield*>& sources, const meta::Inputparameters& params);

/**
 * Calculate 2 flavour chiral condesate from private solution_buffer and store it to a file
 * @param[in] solved The solution spinorfields
 * @param[in] sources The source spinorfields
 * @param[in] corr_fn filename
 * @param[in] number number of gaugefield configuration
 * @param[in] system The system to operate on
 */
void flavour_doublet_chiral_condensate(const std::vector<physics::lattices::Spinorfield*>& solved, const std::vector<physics::lattices::Spinorfield*>& sources, std::string pbp_fn, int number, const hardware::System& system);
}

}

#endif /* _PHYSICS_ALGORITHMS_FLAVOUR_DOUBLET_ */
