/** @todo memset... */
__kernel void set_zero_gaugemomentum(__global aeStorageType * const restrict out)
{
	PARALLEL_FOR(i, GAUGEMOMENTASIZE_MEM) {
		const ae zero = (ae) {
			0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f
		};
		putAe(out, i, zero);
	}
}


