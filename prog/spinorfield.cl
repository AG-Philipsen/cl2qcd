/** @file
 * Device code for operations on the spinor field
 */

spinor get_spinor_from_field(__global const spinor * const restrict in, const int n, const int t)
{
	int pos = get_pos(n, t);
	spinor out;
	out = in[pos];
	return out;
}

void put_spinor_to_field(const spinor in, __global spinor * const restrict out, const int n, const int t)
{
	int pos = get_pos(n, t);
	out[pos] = in;
}
