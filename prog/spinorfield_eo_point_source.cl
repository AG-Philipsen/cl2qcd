__kernel void create_point_source_eoprec(__global spinorfield_eoprec* b, int i, int n)
{
	int id = get_global_id(0);
	if(id == 0) {
	  //LZ: Note that in the kappa-format the source has norm^2 = 1/(2*kappa)
	  //we assume everything on the device to be in kappa-normalization
	  //going back is achieved by multiplying all fields by sqrt(2*kappa)
		hmc_float tmp = 1.;
		int color     = spinor_color(i);
		int spin      = spinor_spin(i,color);
		int pos       = n;

		switch (color){
			
			case 0:
				switch (spin){
					case 0:
						(b[pos].e0).e0.re = tmp;
						break;
					case 1:
						(b[pos].e1).e0.re = tmp;
						break;
					case 2:
						(b[pos].e2).e0.re = tmp;
						break;
					case 3:
						(b[pos].e3).e0.re = tmp;
						break;
				}
				break;
			case 1:
				switch (spin){
					case 0:
						(b[pos].e0).e1.re = tmp;
						break;
					case 1:
						(b[pos].e1).e1.re = tmp;
						break;
					case 2:
						(b[pos].e2).e1.re = tmp;
						break;
					case 3:
						(b[pos].e3).e1.re = tmp;
						break;
				}
				break;
			case 2:
				switch (spin){
					case 0:
						(b[pos].e0).e2.re = tmp;
						break;
					case 1:
						(b[pos].e1).e2.re = tmp;
						break;
					case 2:
						(b[pos].e2).e2.re = tmp;
						break;
					case 3:
						(b[pos].e3).e2.re = tmp;
						break;
				}
				break;
}
	}
	return;
}
