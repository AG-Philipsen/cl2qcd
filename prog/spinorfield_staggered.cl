/** @file
 * Device code for operations on the spinor_staggered field
 */

/**
 * This functions returns the value of the field "in" at the site (n,t)
 */
su3vec get_su3vec_from_field(__global const su3vec * const restrict in, const int n, const int t)
{
	int pos = get_pos(n, t);
	su3vec out;
	out = in[pos];
	return out;
}

/**
 * This functions uploads the value of the field "out" at the site (n,t) with the value "in"
 */
void put_su3vec_to_field(const su3vec in, __global su3vec * const restrict out, const int n, const int t)
{
	int pos = get_pos(n, t);
	out[pos] = in;
}
