hmc_float ae_squarenorm(ae in)
{
	hmc_float result =
	  (in).e0 * (in).e0 +
	  (in).e1 * (in).e1 +
	  (in).e2 * (in).e2 +
	  (in).e3 * (in).e3 +
	  (in).e4 * (in).e4 +
	  (in).e5 * (in).e5 +
	  (in).e6 * (in).e6 +
	  (in).e7 * (in).e7;
	return result;
}

/** @todo add args for reduction... */
__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out)
{
	int id = get_global_id(0);
	if(id == 0) {
		hmc_float result = 0.;
		for(int i = 0; i < GAUGEMOMENTASIZE; i++) {
			result += ae_squarenorm(in[i]);
		}

		/** @todo add reduction.. */
		//CP: Does this work??
		*out = result;
	}
}
