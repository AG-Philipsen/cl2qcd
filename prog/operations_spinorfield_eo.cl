/** @file
 * Device code for operations on the even-odd preconditioned spinor field
 */


//CP: if this works it can also be done directly in the code!!
spinor get_spinor_from_eoprec_field(__global const spinorfield_eoprec * const restrict in, const int n_eoprec)
{
	return in[n_eoprec];
}

void put_spinor_to_eoprec_field(const spinor in, __global spinorfield_eoprec * const restrict out, int n_eoprec)
{
	out[n_eoprec] = in;
}
