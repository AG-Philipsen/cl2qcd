/** @file
 * Device code for operations on the spinor field
 */


#ifdef ENABLE_PRINTF
void print_su3vec(su3vec in)
{
	printf("(%f,%f)\t(%f,%f)\t(%f,%f)", in.e0.re, in.e0.im, in.e1.re, in.e1.im, in.e2.re, in.e2.im);
}

void print_spinor(spinor in)
{
	print_su3vec(in.e0);
	printf("\n");
	print_su3vec(in.e1);
	printf("\n");
	print_su3vec(in.e2);
	printf("\n");
	print_su3vec(in.e3);
	printf("\n");
}
#endif


spinor get_spinor_from_field(__global const spinorfield * const restrict in, const int n, const int t)
{
	int pos = get_global_pos(n, t);
	spinor out;
	out = in[pos];
	return out;
}

void put_spinor_to_field(const spinor in, __global spinorfield * const restrict out, const int n, const int t)
{
	int pos = get_global_pos(n, t);
	out[pos] = in;
}
