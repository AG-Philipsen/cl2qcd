__kernel void sum(__global hmc_complex * a, __global hmc_complex * b, __global hmc_complex * out)
{
	if(get_global_id(0) == 0) {
		hmc_complex tmp1 = complexLoadHack(a);
		hmc_complex tmp2 = complexLoadHack(b);
		(*out) =  complexadd(tmp1, tmp2);
	}
	return;
}
