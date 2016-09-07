/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "host_operations_gaugefield.h"

//taken from the opencl-file

Matrixsu3 multiply_matrixsu3(const Matrixsu3 p, const Matrixsu3 q)
{
	Matrixsu3 out;

	out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e10.re + p.e02.re * q.e20.re
	             - p.e00.im * q.e00.im - p.e01.im * q.e10.im - p.e02.im * q.e20.im;
	out.e00.im = p.e00.re * q.e00.im + p.e01.re * q.e10.im + p.e02.re * q.e20.im
	             + p.e00.im * q.e00.re + p.e01.im * q.e10.re + p.e02.im * q.e20.re;

	out.e01.re = p.e00.re * q.e01.re + p.e01.re * q.e11.re + p.e02.re * q.e21.re
	             - p.e00.im * q.e01.im - p.e01.im * q.e11.im - p.e02.im * q.e21.im;
	out.e01.im = p.e00.re * q.e01.im + p.e01.re * q.e11.im + p.e02.re * q.e21.im
	             + p.e00.im * q.e01.re + p.e01.im * q.e11.re + p.e02.im * q.e21.re;

	out.e02.re = p.e00.re * q.e02.re + p.e01.re * q.e12.re + p.e02.re * q.e22.re
	             - p.e00.im * q.e02.im - p.e01.im * q.e12.im - p.e02.im * q.e22.im;
	out.e02.im = p.e00.re * q.e02.im + p.e01.re * q.e12.im + p.e02.re * q.e22.im
	             + p.e00.im * q.e02.re + p.e01.im * q.e12.re + p.e02.im * q.e22.re;

	out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e10.re + p.e12.re * q.e20.re
	             - p.e10.im * q.e00.im - p.e11.im * q.e10.im - p.e12.im * q.e20.im;
	out.e10.im = p.e10.re * q.e00.im + p.e11.re * q.e10.im + p.e12.re * q.e20.im
	             + p.e10.im * q.e00.re + p.e11.im * q.e10.re + p.e12.im * q.e20.re;

	out.e11.re = p.e10.re * q.e01.re + p.e11.re * q.e11.re + p.e12.re * q.e21.re
	             - p.e10.im * q.e01.im - p.e11.im * q.e11.im - p.e12.im * q.e21.im;
	out.e11.im = p.e10.re * q.e01.im + p.e11.re * q.e11.im + p.e12.re * q.e21.im
	             + p.e10.im * q.e01.re + p.e11.im * q.e11.re + p.e12.im * q.e21.re;

	out.e12.re = p.e10.re * q.e02.re + p.e11.re * q.e12.re + p.e12.re * q.e22.re
	             - p.e10.im * q.e02.im - p.e11.im * q.e12.im - p.e12.im * q.e22.im;
	out.e12.im = p.e10.re * q.e02.im + p.e11.re * q.e12.im + p.e12.re * q.e22.im
	             + p.e10.im * q.e02.re + p.e11.im * q.e12.re + p.e12.im * q.e22.re;

	out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e10.re + p.e22.re * q.e20.re
	             - p.e20.im * q.e00.im - p.e21.im * q.e10.im - p.e22.im * q.e20.im;
	out.e20.im = p.e20.re * q.e00.im + p.e21.re * q.e10.im + p.e22.re * q.e20.im
	             + p.e20.im * q.e00.re + p.e21.im * q.e10.re + p.e22.im * q.e20.re;

	out.e21.re = p.e20.re * q.e01.re + p.e21.re * q.e11.re + p.e22.re * q.e21.re
	             - p.e20.im * q.e01.im - p.e21.im * q.e11.im - p.e22.im * q.e21.im;
	out.e21.im = p.e20.re * q.e01.im + p.e21.re * q.e11.im + p.e22.re * q.e21.im
	             + p.e20.im * q.e01.re + p.e21.im * q.e11.re + p.e22.im * q.e21.re;

	out.e22.re = p.e20.re * q.e02.re + p.e21.re * q.e12.re + p.e22.re * q.e22.re
	             - p.e20.im * q.e02.im - p.e21.im * q.e12.im - p.e22.im * q.e22.im;
	out.e22.im = p.e20.re * q.e02.im + p.e21.re * q.e12.im + p.e22.re * q.e22.im
	             + p.e20.im * q.e02.re + p.e21.im * q.e12.re + p.e22.im * q.e22.re;

	return out;
}

Matrixsu3 multiply_matrixsu3_dagger(const Matrixsu3 p, const Matrixsu3 q)
{
	Matrixsu3 out;

	out.e00.re = p.e00.re * q.e00.re + p.e01.re * q.e01.re + p.e02.re * q.e02.re
	             + p.e00.im * q.e00.im + p.e01.im * q.e01.im + p.e02.im * q.e02.im;
	out.e00.im = -p.e00.re * q.e00.im - p.e01.re * q.e01.im - p.e02.re * q.e02.im
	             + p.e00.im * q.e00.re + p.e01.im * q.e01.re + p.e02.im * q.e02.re;

	out.e01.re = p.e00.re * q.e10.re + p.e01.re * q.e11.re + p.e02.re * q.e12.re
	             + p.e00.im * q.e10.im + p.e01.im * q.e11.im + p.e02.im * q.e12.im;
	out.e01.im = -p.e00.re * q.e10.im - p.e01.re * q.e11.im - p.e02.re * q.e12.im
	             + p.e00.im * q.e10.re + p.e01.im * q.e11.re + p.e02.im * q.e12.re;

	out.e02.re = p.e00.re * q.e20.re + p.e01.re * q.e21.re + p.e02.re * q.e22.re
	             + p.e00.im * q.e20.im + p.e01.im * q.e21.im + p.e02.im * q.e22.im;
	out.e02.im = -p.e00.re * q.e20.im - p.e01.re * q.e21.im - p.e02.re * q.e22.im
	             + p.e00.im * q.e20.re + p.e01.im * q.e21.re + p.e02.im * q.e22.re;

	out.e10.re = p.e10.re * q.e00.re + p.e11.re * q.e01.re + p.e12.re * q.e02.re
	             + p.e10.im * q.e00.im + p.e11.im * q.e01.im + p.e12.im * q.e02.im;
	out.e10.im = -p.e10.re * q.e00.im - p.e11.re * q.e01.im - p.e12.re * q.e02.im
	             + p.e10.im * q.e00.re + p.e11.im * q.e01.re + p.e12.im * q.e02.re;

	out.e11.re = p.e10.re * q.e10.re + p.e11.re * q.e11.re + p.e12.re * q.e12.re
	             + p.e10.im * q.e10.im + p.e11.im * q.e11.im + p.e12.im * q.e12.im;
	out.e11.im = -p.e10.re * q.e10.im - p.e11.re * q.e11.im - p.e12.re * q.e12.im
	             + p.e10.im * q.e10.re + p.e11.im * q.e11.re + p.e12.im * q.e12.re;

	out.e12.re = p.e10.re * q.e20.re + p.e11.re * q.e21.re + p.e12.re * q.e22.re
	             + p.e10.im * q.e20.im + p.e11.im * q.e21.im + p.e12.im * q.e22.im;
	out.e12.im = -p.e10.re * q.e20.im - p.e11.re * q.e21.im - p.e12.re * q.e22.im
	             + p.e10.im * q.e20.re + p.e11.im * q.e21.re + p.e12.im * q.e22.re;

	out.e20.re = p.e20.re * q.e00.re + p.e21.re * q.e01.re + p.e22.re * q.e02.re
	             + p.e20.im * q.e00.im + p.e21.im * q.e01.im + p.e22.im * q.e02.im;
	out.e20.im = -p.e20.re * q.e00.im - p.e21.re * q.e01.im - p.e22.re * q.e02.im
	             + p.e20.im * q.e00.re + p.e21.im * q.e01.re + p.e22.im * q.e02.re;

	out.e21.re = p.e20.re * q.e10.re + p.e21.re * q.e11.re + p.e22.re * q.e12.re
	             + p.e20.im * q.e10.im + p.e21.im * q.e11.im + p.e22.im * q.e12.im;
	out.e21.im = -p.e20.re * q.e10.im - p.e21.re * q.e11.im - p.e22.re * q.e12.im
	             + p.e20.im * q.e10.re + p.e21.im * q.e11.re + p.e22.im * q.e12.re;

	out.e22.re = p.e20.re * q.e20.re + p.e21.re * q.e21.re + p.e22.re * q.e22.re
	             + p.e20.im * q.e20.im + p.e21.im * q.e21.im + p.e22.im * q.e22.im;
	out.e22.im = -p.e20.re * q.e20.im - p.e21.re * q.e21.im - p.e22.re * q.e22.im
	             + p.e20.im * q.e20.re + p.e21.im * q.e21.re + p.e22.im * q.e22.re;

	return out;
}

Matrixsu3 unit_matrixsu3()
{
	Matrixsu3 out;
	out.e00.re = 1.;
	out.e00.im = 0.;
	out.e01.re = 0.;
	out.e01.im = 0.;
	out.e02.re = 0.;
	out.e02.im = 0.;

	out.e10.re = 0.;
	out.e10.im = 0.;
	out.e11.re = 1.;
	out.e11.im = 0.;
	out.e12.re = 0.;
	out.e12.im = 0.;

	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 1.;
	out.e22.im = 0.;

	return out;
}

Matrixsu3 nonTrivialSu3Matrix()
{
	Matrixsu3 out;
	out.e00.re = .130189;
	out.e00.im = .260378;
	out.e01.re = .260378;
	out.e01.im = .390567;
	out.e02.re = .520756;
	out.e02.im = .650945;

	out.e10.re = .572742;
	out.e10.im = .403041;
	out.e11.re = .371222;
	out.e11.im = .321726;
	out.e12.re = -.449002;
	out.e12.im = -.258088;

	out.e20.re = 1.11022e-16;
	out.e20.im = .651751;
	out.e21.re = .0271563;
	out.e21.im = -.733219;
	out.e22.re = -.0271563;
	out.e22.im = .190094;

	return out;
}

//CP: I introduced explicit calculations of the neighbors because this does not rely on any geometric conventions!
Matrixsu3 local_plaquette(Matrixsu3 * field, int coord_in[NDIM], int mu, int nu, const meta::Inputparameters& parameters)
{
	Matrixsu3 res, tmp;
	//coordinates of neighbors
	int coord[NDIM];
	coord[TDIR] = coord_in[TDIR];
	coord[XDIR] = coord_in[XDIR];
	coord[YDIR] = coord_in[YDIR];
	coord[ZDIR] = coord_in[ZDIR];
	//spatial index

	const size_t NTIME = parameters.get_ntime();
	const size_t NSPACE = parameters.get_nspace();
	//u_mu(x)
	res = field[uint(LinkIndex(Index(coord[1], coord[2], coord[3], coord[TDIR], LatticeExtents(parameters.get_nspace(),parameters.get_ntime())),static_cast<Direction>(mu)))];
	//u_nu(x+mu)
	if(mu == TDIR) {
		coord[mu] = (coord_in[mu] + 1) % NTIME;
		tmp = field[uint(LinkIndex(Index(coord[1], coord[2], coord[3], coord[mu], LatticeExtents(parameters.get_nspace(),parameters.get_ntime())),static_cast<Direction>(nu)))];
		coord[mu] = coord_in[mu];
	} else {
		coord[mu] = (coord_in[mu] + 1) % NSPACE;
		tmp = field[uint(LinkIndex(Index(coord[1], coord[2], coord[3], coord[TDIR], LatticeExtents(parameters.get_nspace(),parameters.get_ntime())),static_cast<Direction>(nu)))];
		coord[mu] = coord_in[mu];
	}
	//accumulate_su3matrix_prod(&prod, &tmp);
	res = multiply_matrixsu3(res, tmp);
	//adjoint(u_mu(x+nu))
	if(nu == TDIR) {
		coord[nu] = (coord_in[nu] + 1) % NTIME;
		tmp = field[uint(LinkIndex(Index(coord[1], coord[2], coord[3], coord[nu], LatticeExtents(parameters.get_nspace(),parameters.get_ntime())),static_cast<Direction>(mu)))];
		coord[nu] = coord_in[nu];
	} else {
		coord[nu] = (coord_in[nu] + 1) % NSPACE;
		tmp = field[uint(LinkIndex(Index(coord[1], coord[2], coord[3], coord_in[TDIR], LatticeExtents(parameters.get_nspace(),parameters.get_ntime())),static_cast<Direction>(mu)))];
		coord[nu] = coord_in[nu];
	}
	res = multiply_matrixsu3_dagger(res, tmp);

	//adjoint(u_nu(x))
	tmp = field[uint(LinkIndex(Index(coord[1], coord[2], coord[3], coord[TDIR], LatticeExtents(parameters.get_nspace(),parameters.get_ntime())),static_cast<Direction>(nu)))];
	res = multiply_matrixsu3_dagger(res, tmp);

	return res;
}
