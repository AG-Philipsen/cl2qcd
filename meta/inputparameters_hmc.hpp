/** @file
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#ifndef _META_INPUTPARAMETERS_HMC_HPP_
#define _META_INPUTPARAMETERS_HMC_HPP_

#include "inputparameters.hpp"

namespace meta{
	po::options_description Inputparameters::getParameters_hmc()
	{
		po::options_description options("HMC options");
		options.add_options()
		("tau", po::value<double>(&tau)->default_value(0.5))
		("reversibility_check", po::value<bool>(&reversibility_check)->default_value(false))
		("integrationsteps0", po::value<int>(&integrationsteps0)->default_value(10))
		("integrationsteps1", po::value<int>(&integrationsteps1)->default_value(10))
		("integrationsteps2", po::value<int>(&integrationsteps2)->default_value(10))
		("hmcsteps", po::value<int>(&hmcsteps)->default_value(10))
		("benchmarksteps", po::value<int>(&benchmarksteps)->default_value(500))
		("num_timescales", po::value<int>(&num_timescales)->default_value(1))
		("integrator0", po::value<std::string>()->default_value("leapfrog"))
		("integrator1", po::value<std::string>()->default_value("leapfrog"))
		("integrator2", po::value<std::string>()->default_value("leapfrog"))
		// this is the optimal value...
		("lambda0", po::value<double>(&lambda0)->default_value(0.1931833275037836))
		("lambda1", po::value<double>(&lambda1)->default_value(0.1931833275037836))
		("lambda2", po::value<double>(&lambda2)->default_value(0.1931833275037836))
		("use_gauge_only", po::value<bool>(&use_gauge_only)->default_value(false))
		("use_mp", po::value<bool>(&use_mp)->default_value(false));

		return options;
	}

}

#endif
