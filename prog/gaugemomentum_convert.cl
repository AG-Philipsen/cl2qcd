__kernel void gaugemomentum_convert_to_soa(__global aeStorageType * const restrict dest, __global const ae * const restrict src)
{
	PARALLEL_FOR(s, VOL4D) {
		for(uint d = 0; d < 4; ++d) {
			const st_idx site = get_st_idx_from_site_idx(s);
			const link_idx aos_i = get_link_idx_AOS(d, site);
			const link_idx soa_i = get_link_idx(d, site);
			putAe(dest, soa_i, src[aos_i]);
		}
	}
}

__kernel void gaugemomentum_convert_from_soa(__global ae * const restrict dest, __global const aeStorageType * const restrict src)
{
	PARALLEL_FOR(s, VOL4D) {
		for(uint d = 0; d < 4; ++d) {
			const st_idx site = get_st_idx_from_site_idx(s);
			const link_idx aos_i = get_link_idx_AOS(d, site);
			const link_idx soa_i = get_link_idx(d, site);
			dest[aos_i] = getAe(src, soa_i);
		}
	}
}
