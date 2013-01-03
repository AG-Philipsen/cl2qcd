/** @file
 * Declaration of the flavour doublets algorithms
 */

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
void flavour_doublet_correlators(const std::vector<const physics::lattices::Spinorfield*>& result, const std::vector<const physics::lattices::Spinorfield*>& sources, std::string corr_fn, const meta::Inputparameters& params);

/**
 * Calculate 2 flavour chiral condesate from private solution_buffer and store it to a file
 * @param[in] gaugefield The gaugefield
 * @param[in] result The Spinorfields
 * @param[in] corr_fn filename
 * @param[in] number number of gaugefield configuration
 */
void flavour_doublet_chiral_condensate(const physics::lattices::Gaugefield& gaugefield, const std::vector<const physics::lattices::Spinorfield*>& result, std::string pbp_fn, int number);
}

}
