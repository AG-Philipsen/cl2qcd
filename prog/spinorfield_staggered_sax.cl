// Description of variables of sax:
//  - x: The input staggered field (an su3vec per each site => vector of VOL4D
//       components that are su3vec varibles)
//  - aplha: The complex number by which x has to be multiplied
//  - out: The output staggered field: alpha*x (site by site)

__kernel void sax_staggered(__global const su3vec * const restrict x, __global const hmc_complex * const restrict alpha, __global su3vec * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	hmc_complex alpha_tmp = complexLoadHack(alpha);
	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		su3vec x_tmp = x[id_mem];
		out[id_mem] = su3vec_times_complex(x_tmp, alpha_tmp);
	}
}
