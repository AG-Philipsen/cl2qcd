/** @file
 * Implementation of the heatbath algorithm
 */

#include "kappa_clover.hpp"

#include <cassert>
#include "../hardware/code/kappa.hpp"

hmc_float physics::algorithms::kappa_clover(physics::lattices::Gaugefield& gf, hmc_float beta)
{
	assert(gf.get_buffers().size() == 1);

	auto gf_dev = gf.get_buffers()[0];

	auto device = gf_dev->get_device();
	hardware::buffers::Plain<hmc_float> kappa_clover_dev(1, device);
	gf_dev->get_device()->get_kappa_code()->run_kappa_clover(&kappa_clover_dev, gf_dev, beta);

	hmc_float kappa_clover_host;
	kappa_clover_dev.dump(&kappa_clover_host);
	return kappa_clover_host;
}
