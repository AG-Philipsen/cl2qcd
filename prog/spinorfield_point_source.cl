__kernel void create_point_source(__global spinor * const restrict b, int i, int spacepos, int timepos)
{
	int id = get_global_id(0);
	if(id == 0) {
		//LZ: note that the conversion from m to kappa works as
		//   M(m) * phi = eta <=> 1/(2kappa) * M(k) * phi = eta
		//thus we can keep everything as in the orginal basis and only have
		//to multiply the resulting field phi by 2kappa
		hmc_float tmp = 1.;
		int color = spinor_color(i);
		int spin = spinor_spin(i, color);
		int pos = get_pos(spacepos, timepos);
		b[pos] = set_spinor_zero();
		switch (color) {

			case 0:
				switch (spin) {
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
				switch (spin) {
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
				switch (spin) {
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
