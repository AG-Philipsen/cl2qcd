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
spinor getSpinorSOA_eo(__global const hmc_complex * const restrict in, const uint idx)
{
	return (spinor) {
		{
			// su3vec = 3 * cplx
			in[ 0 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 1 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 2 * EOPREC_SPINORFIELD_STRIDE + idx]
		}, {
			// su3vec = 3 * cplx
			in[ 3 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 4 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 5 * EOPREC_SPINORFIELD_STRIDE + idx]
		}, {
			// su3vec = 3 * cplx
			in[ 6 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 7 * EOPREC_SPINORFIELD_STRIDE + idx], in[ 8 * EOPREC_SPINORFIELD_STRIDE + idx]
		}, {
			// su3vec = 3 * cplx
			in[ 9 * EOPREC_SPINORFIELD_STRIDE + idx], in[10 * EOPREC_SPINORFIELD_STRIDE + idx], in[11 * EOPREC_SPINORFIELD_STRIDE + idx]
		}
	};
}

// TODO document
void putSpinorSOA_eo(__global hmc_complex * const restrict out, const uint idx, const spinor val)
{
	// su3vec = 3 * cplx
	out[ 0 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e0;
	out[ 1 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e1;
	out[ 2 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0.e2;

	// su3vec = 3 * cplx
	out[ 3 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e0;
	out[ 4 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e1;
	out[ 5 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1.e2;

	// su3vec = 3 * cplx
	out[ 6 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e0;
	out[ 7 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e1;
	out[ 8 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2.e2;

	// su3vec = 3 * cplx
	out[ 9 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e0;
	out[10 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e1;
	out[11 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e3.e2;
}
