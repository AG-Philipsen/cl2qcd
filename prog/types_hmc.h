/**
 @file types needed in the HMC-algorithm
 *
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

#ifndef _TYPES_HMCH_
#define _TYPES_HMCH_

#ifndef _INKERNEL_
//CP: define struct for observables
struct hmc_observables {
	hmc_float plaq;
	hmc_float tplaq;
	hmc_float splaq;
	hmc_complex poly;
	hmc_float deltaH;
	hmc_float prob;
	int accept;
	hmc_float rectangles;
};
#endif

#endif // _TYPES_HMC

