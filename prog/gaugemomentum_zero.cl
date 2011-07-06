ae set_zero_ae()
{
	ae tmp;
	tmp.e0 = 0.;
	tmp.e1 = 0.;
	tmp.e2 = 0.;
	tmp.e3 = 0.;
	tmp.e4 = 0.;
	tmp.e5 = 0.;
	tmp.e6 = 0.;
	tmp.e7 = 0.;
	return tmp;
}

/** @todo memset... */
__kernel void set_zero_gaugemomentum(__global ae * in)
{
	int id = get_global_id(0);
	if(id == 0) {
		for(int i = 0; i < GAUGEMOMENTASIZE; i++) {
			in[i] = set_zero_ae();
		}
	}
}


