/** @file
 * Device code for operations on the spinor field
 */

/** @todo CP: this must be taken out of this code again! */
#define EOPREC_SPINORFIELDSIZE2 EOPREC_SPINORFIELDSIZE

spinor get_spinor_from_field(__global spinorfield* in, int n, int t)
{
	int pos = get_global_pos(n,t);
	spinor out;
	out = in[pos];
	return out;
}

void put_spinor_to_field(spinor in, __global spinorfield* out, int n, int t)
{
	int pos = get_global_pos(n,t);
	out[pos] = in;
}
