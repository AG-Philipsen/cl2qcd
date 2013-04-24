__kernel void convertGaugefieldToSOA(__global Matrixsu3StorageType * const restrict out, __global const Matrixsu3 * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(site_idx s = get_global_id(0); s < VOL4D_MEM; s += get_global_size(0)) {
			const st_idx site = get_st_idx_from_site_idx(s);
			Matrixsu3 tmp = in[get_link_idx_AOS(d, site)];

			putSU3(out, get_link_idx_SOA(d, site), tmp);
		}
	}
}
__kernel void convertGaugefieldFromSOA(__global Matrixsu3 * const restrict out, __global const Matrixsu3StorageType * const restrict in)
{
	// we need to take care of index converion. in the AOS storage the dimension is the continious index
	// in the soa storage we want the space indices to be continuous and have the dimension as outermost.
	for(uint d = 0; d < NDIM; ++d) {
		for(site_idx s = get_global_id(0); s < VOL4D_MEM; s += get_global_size(0)) {
			const st_idx site = get_st_idx_from_site_idx(s);
			Matrixsu3 tmp = getSU3(in, get_link_idx_SOA(d, site));

			out[get_link_idx_AOS(d, site)] = tmp;
		}
	}
}
