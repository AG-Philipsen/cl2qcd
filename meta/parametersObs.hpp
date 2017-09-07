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

#ifndef _META_PARAMETERS_OBS_HPP_
#define _META_PARAMETERS_OBS_HPP_

#include "parametersBasic.hpp"

namespace meta {
class ParametersObs {
public:

	bool get_measure_transportcoefficient_kappa() const noexcept;
	bool get_measure_rectangles() const noexcept;
	bool get_measure_correlators() const noexcept;
	bool get_measure_pbp() const noexcept;
	common::pbp_version get_pbp_version() const noexcept;
	int get_corr_dir() const noexcept;
	int get_pbp_measurements() const noexcept;

private:
	po::options_description options;

	bool measure_transportcoefficient_kappa;
	bool measure_rectangles;
	bool measure_correlators;
	bool measure_pbp;
	int corr_dir;
	int pbp_measurements;

protected:
	ParametersObs();
	virtual ~ParametersObs();
	ParametersObs(ParametersObs const&) = delete;
	ParametersObs & operator=(ParametersObs const&) = delete;
	po::options_description & getOptions();

	common::pbp_version pbp_version_;
};

}

#endif
