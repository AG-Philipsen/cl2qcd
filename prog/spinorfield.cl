/** @file
 * Device code for operations on the spinor field
 */


void print_su3vec(su3vec in){
     printf("(%f,%f)\t(%f,%f)\t(%f,%f)\t", in.e0.re, in.e0.im, in.e1.re, in.e1.im, in.e2.re, in.e2.im);
}

void print_spinor(spinor in){
     print_su3vec(in.e0);
     print_su3vec(in.e1);
     print_su3vec(in.e2);
     print_su3vec(in.e3);
     printf("\n");
}


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
