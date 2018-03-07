/*
 * Copyright (c) 2014,2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */


inline hmc_float update_zeta_cgm_alg(hmc_float a, hmc_float b, hmc_float c, hmc_float d, hmc_float e, hmc_float f)
{
	hmc_float out;
	out = (a * b * c) / (d * e * (b - a) + b * c * (1 - f * d));
	return out;
}

inline hmc_float update_beta_cgm_alg(hmc_float a, hmc_float b, hmc_float c)
{
	hmc_float out;
	out = - a * b / c;
	return out;
}

inline hmc_float update_alpha_cgm_alg(hmc_float a, hmc_float b, hmc_float c, hmc_float d, hmc_float e)
{
	hmc_float out;
	out = - a * b * c / (d * e);
	return out;
}
