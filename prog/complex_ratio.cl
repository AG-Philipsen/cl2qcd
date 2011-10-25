__kernel void ratio(__global hmc_complex * a, __global hmc_complex * b, __global hmc_complex * out)
{
	//!!CP: complexdivide cannot handle __global
	if(get_global_id(0) == 0) {
		hmc_complex tmp1 = complexLoadHack(a);
		hmc_complex tmp2 = complexLoadHack(b);
		(*out) =  complexdivide(tmp1, tmp2);
	}
	return;
}
