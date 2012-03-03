__kernel void gauge_force_tlsym(__global ocl_s_gaugefield * field, __global ae * out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int id_tmp = id; id_tmp < VOL4D*NDIM; id_tmp += global_size) {

	}
}
