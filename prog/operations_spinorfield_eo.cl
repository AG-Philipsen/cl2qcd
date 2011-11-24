/** @file
 * Device code for operations on the even-odd preconditioned spinor field
 */

// TODO remove legacy code starting here

//CP: if this works it can also be done directly in the code!!
spinor get_spinor_from_eoprec_field(__global const spinorfield_eoprec * const restrict in, const int n_eoprec)
{
	return in[n_eoprec];
}

void put_spinor_to_eoprec_field(const spinor in, __global spinorfield_eoprec * const restrict out, int n_eoprec)
{
	out[n_eoprec] = in;
}

// TODO remove legacy code ending here

// TODO document
spinor getSpinorSOA_eo(__global const hmc_float * const restrict in, const uint idx)
{
	return (spinor) {
		{
			// su3vec = 3 * cplx
			{
				in[ 0 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 1 * EOPREC_SPINORFIELD_STRIDE + idx]
			},
			{ in[ 2 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 3 * EOPREC_SPINORFIELD_STRIDE + idx] },
			{ in[ 4 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 5 * EOPREC_SPINORFIELD_STRIDE + idx] },
		}, {
			// su3vec = 3 * cplx
			{ in[ 6 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 7 * EOPREC_SPINORFIELD_STRIDE + idx] },
			{ in[ 8 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 9 * EOPREC_SPINORFIELD_STRIDE + idx] },
			{ in[10 * EOPREC_SPINORFIELD_STRIDE + idx], in[11 * EOPREC_SPINORFIELD_STRIDE + idx] },
		}, {
			// su3vec = 3 * cplx
			{ in[12 * EOPREC_SPINORFIELD_STRIDE + idx], in[13 * EOPREC_SPINORFIELD_STRIDE + idx] },
			{ in[14 * EOPREC_SPINORFIELD_STRIDE + idx], in[15 * EOPREC_SPINORFIELD_STRIDE + idx] },
			{ in[16 * EOPREC_SPINORFIELD_STRIDE + idx], in[17 * EOPREC_SPINORFIELD_STRIDE + idx] },
		}, {
			// su3vec = 3 * cplx
			{ in[18 * EOPREC_SPINORFIELD_STRIDE + idx], in[19 * EOPREC_SPINORFIELD_STRIDE + idx] },
			{ in[20 * EOPREC_SPINORFIELD_STRIDE + idx], in[21 * EOPREC_SPINORFIELD_STRIDE + idx] },
			{ in[22 * EOPREC_SPINORFIELD_STRIDE + idx], in[23 * EOPREC_SPINORFIELD_STRIDE + idx] },
		}
	};
}

// TODO document
void putSpinorSOA_eo(__global hmc_float * const restrict out, const uint idx, const spinor val)
{
	// su3vec = 3 * cplx
	out[ 0 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e0.re;
	out[ 1 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e0.im;
	out[ 2 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e1.re;
	out[ 3 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e1.im;
	out[ 4 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e2.re;
	out[ 5 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e2.im;

	// su3vec = 3 * cplx
	out[ 6 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e0.re;
	out[ 7 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e0.im;
	out[ 8 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e1.re;
	out[ 9 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e1.im;
	out[10 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e2.re;
	out[11 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e2.im;

	// su3vec = 3 * cplx
	out[12 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e0.re;
	out[13 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e0.im;
	out[14 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e1.re;
	out[15 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e1.im;
	out[16 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e2.re;
	out[17 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e2.im;

	// su3vec = 3 * cplx
	out[18 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e0.re;
	out[19 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e0.im;
	out[20 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e1.re;
	out[21 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e1.im;
	out[22 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e2.re;
	out[23 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e2.im;
}
