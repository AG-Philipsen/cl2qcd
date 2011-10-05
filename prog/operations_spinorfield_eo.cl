/** @file
 * Device code for operations on the even-odd preconditioned spinor field
 */


//CP: if this works it can also be done directly in the code!!
spinor get_spinor_from_eoprec_field(__global spinorfield_eoprec* in, int n_eoprec)
{
	spinor out;
	out = in[n_eoprec];
	return out;
}

void put_spinor_to_eoprec_field(spinor in, __global spinorfield_eoprec * out, int n_eoprec)
{
	out[n_eoprec] = in;
}
