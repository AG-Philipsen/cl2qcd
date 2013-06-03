/** @file
 * Device code for operations on the even-odd preconditioned staggered field
 */

/**
 * This function returns the value of the staggered field _@b in in the site with superindex_eo idx.
 *
 * @NOTE If _USE_SOA_ is true the output vector has to be reconstructed using properly
 * the option EOPREC_SPINORFIELD_STRIDE.
 */
inline su3vec get_su3vec_from_field_eo(__global const staggeredStorageType * const restrict in, const uint idx)
{
#ifdef _USE_SOA_
	return (su3vec) {
		// su3vec = 3 * cplx
		in[ 0 * EOPREC_SPINORFIELD_STRIDE + idx],
		in[ 1 * EOPREC_SPINORFIELD_STRIDE + idx],
		in[ 2 * EOPREC_SPINORFIELD_STRIDE + idx]
	};
#else
	return in[idx];
#endif
}

/**
 * This function puts the su3vec _@b val into the staggered field _@b out in the site with superindex_eo idx.
 *
 * @NOTE If _USE_SOA_ is true, each component of the su3vec _@b val has to be placed into _@b out properly
 * using the option EOPREC_SPINORFIELD_STRIDE.
 */
inline void put_su3vec_to_field_eo(__global staggeredStorageType * const restrict out, const uint idx, const su3vec val)
{
#ifdef _USE_SOA_
	// su3vec = 3 * cplx
	out[ 0 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e0;
	out[ 1 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e1;
	out[ 2 * EOPREC_SPINORFIELD_STRIDE + idx] = val.e2;
#else
	out[idx] = val;
#endif
}
