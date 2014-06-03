/** @file
 * Input file handling
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#ifndef _META_INPUTPARAMETERS_HPP_
#define _META_INPUTPARAMETERS_HPP_

#include "../host_functionality/logger.hpp"

#include <boost/algorithm/string.hpp>

#include "parametersBasic.hpp"
#include "parametersConfig.hpp"
#include "parametersIo.hpp"
#include "parametersObs.hpp"
#include "parametersGauge.hpp"
#include "parametersFermion.hpp"
#include "parametersSources.hpp"
#include "parametersSolver.hpp"
#include "parametersHmc.hpp"
#include "parametersHeatbath.hpp"
#include "parametersRhmc.hpp"
#include "parametersTest.hpp"

/**
 * This namespace contains generic utility code required by the other packages.
 */
namespace meta {

/**
 * Parser and representation of an input file.
 *
 * This class is copyable and assignable, but should
 * be used as a const value after initialization.
 */
	class Inputparameters : 
		public ParametersConfig,
		public ParametersIo,
		public ParametersObs, 
		public ParametersGauge,
		public ParametersFermion,
		public ParametersSources,
		public ParametersSolver,
		public ParametersHmc, 
		public ParametersHeatbath,
		public ParametersRhmc,
		public ParametersTest
	{

public:
		/**
		* The parsing of the input parameters aborted for some reason.
		* Could either be invalid input or the specification of --help.
		*/
		struct parse_aborted {};

		/**
		* Construct from command line and config file.
		*
		* Config file will be retrieved from command line.
		*
		* @throws parse_aborted
		*/
		Inputparameters(int argc, const char** argv, std::string parameterSet = "allParameters");
	};
}

#endif /* _META_INPUTPARAMETERS_H_ */
