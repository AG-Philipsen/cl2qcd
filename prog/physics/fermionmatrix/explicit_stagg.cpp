/** @file
 * Implementation of explicit fermionamtrix operations
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "fermionmatrix_stagg.hpp"


void physics::fermionmatrix::DKS_eo(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& in, int evenodd)
{
	auto out_bufs = out->get_buffers();
	auto gf_bufs = gf.get_buffers();
	auto in_bufs = in.get_buffers();

	size_t num_bufs = out_bufs.size();
	if(num_bufs != gf_bufs.size() || num_bufs != in_bufs.size()) {
		throw std::invalid_argument("Given lattices do not use the same devices");
	}

	for(size_t i = 0; i < num_bufs; ++i) {
		auto fermion_code = out_bufs[i]->get_device()->get_fermion_staggered_code();
		fermion_code->D_KS_eo_device(in_bufs[i], out_bufs[i], gf_bufs[i], evenodd);
	}
	if(num_bufs!=1)
	  out->update_halo();
}

