/** @todo add args for reduction... */
__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out)
{
	int id = get_global_id(0);
	if(id == 0) {
		hmc_float result = 0.;
		for(int i = 0; i < GAUGEMOMENTASIZE; i++) {
		  /*
		  int mu_check = 0;
		  if(i%4 != mu_check) continue;
		  //without mu_check (0..NDIM * 256)
		  int j = i;
		  //get rid of NDIM
		  j = j / NDIM;
		  if(j%4 != 2) continue;
		  */
		  result += ae_squarenorm(in[i]);
		}

		/** @todo add reduction.. */
		//CP: Does this work??
		*out = result;
	}
}
