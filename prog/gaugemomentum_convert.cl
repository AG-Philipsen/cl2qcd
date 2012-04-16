__kernel void gaugemomentum_convert_to_soa(__global aeStorageType * const restrict dest, __global const ae * const restrict src)
{
	PARALLEL_FOR(i, GAUGEMOMENTASIZE) {
		putAe(dest, i, src[i]);
	}
}

__kernel void gaugemomentum_convert_from_soa(__global ae * const restrict dest, __global const aeStorageType * const restrict src)
{
	PARALLEL_FOR(i, GAUGEMOMENTASIZE) {
		dest[i] = getAe(src, i);
	}
}
