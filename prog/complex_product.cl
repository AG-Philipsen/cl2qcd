__kernel void product(__global hmc_complex * a, __global hmc_complex * b, __global hmc_complex * out)
{
	//!!CP: complexdivide cannot handle __global
	hmc_complex tmp1 = (*a);
	hmc_complex tmp2 = (*b);
	(*out) =  complexmult(tmp1, tmp2);
	return;
}
