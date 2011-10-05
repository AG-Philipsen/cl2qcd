/**
 * @file operations on algebraelements
 */

void print_ae(ae in)
{
	printf("%f, %f, %f, %f, %f, %f, %f, %f\n", in.e0, in.e1, in.e2, in.e3, in.e4, in.e5, in.e6, in.e7);
}

ae acc_factor_times_algebraelement(ae in, hmc_float factor, ae force_in)
{
	ae tmp;
	tmp.e0 = in.e0 + factor * force_in.e0;
	tmp.e1 = in.e1 + factor * force_in.e1;
	tmp.e2 = in.e2 + factor * force_in.e2;
	tmp.e3 = in.e3 + factor * force_in.e3;
	tmp.e4 = in.e4 + factor * force_in.e4;
	tmp.e5 = in.e5 + factor * force_in.e5;
	tmp.e6 = in.e6 + factor * force_in.e6;
	tmp.e7 = in.e7 + factor * force_in.e7;
	return tmp;
}

ae acc_algebraelement(ae in, ae force_in)
{
	ae tmp;
	tmp.e0 = in.e0 + force_in.e0;
	tmp.e1 = in.e1 + force_in.e1;
	tmp.e2 = in.e2 + force_in.e2;
	tmp.e3 = in.e3 + force_in.e3;
	tmp.e4 = in.e4 + force_in.e4;
	tmp.e5 = in.e5 + force_in.e5;
	tmp.e6 = in.e6 + force_in.e6;
	tmp.e7 = in.e7 + force_in.e7;
	return tmp;
}

void update_gaugemomentum(ae in, hmc_float factor, int global_link_pos, __global ae * out)
{
	ae tmp = out[global_link_pos];
	tmp = acc_factor_times_algebraelement(tmp, factor, in);
	out[global_link_pos] = tmp;
}

//calculates the trace of i times generator times 3x3-matrix and stores this in a su3-algebraelement
ae tr_lambda_u(Matrix3x3 in)
{
	ae tmp;
	tmp.e0 = ( -in.e10.im - in.e01.im);
	tmp.e1 = (+in.e10.re - in.e01.re);
	tmp.e2 = (-in.e00.im + in.e11.im);
	tmp.e3 = (-in.e20.im - in.e02.im);
	tmp.e4 = (+in.e20.re - in.e02.re);
	tmp.e5 = (-in.e21.im - in.e12.im);
	tmp.e6 = (+in.e21.re - in.e12.re);
	tmp.e7 = (-in.e00.im - in.e11.im + 2.0 * in.e22.im) * 0.577350269189625;
	return tmp;
}

ae ae_times_factor(ae in, hmc_float factor)
{
	ae tmp;

	tmp.e0 = in.e0 * factor;
	tmp.e1 = in.e1 * factor;
	tmp.e2 = in.e2 * factor;
	tmp.e3 = in.e3 * factor;
	tmp.e4 = in.e4 * factor;
	tmp.e5 = in.e5 * factor;
	tmp.e6 = in.e6 * factor;
	tmp.e7 = in.e7 * factor;

	return tmp;
}
