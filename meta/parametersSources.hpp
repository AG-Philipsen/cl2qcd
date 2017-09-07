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

#ifndef _META_PARAMETERS_SOURCES_HPP_
#define _META_PARAMETERS_SOURCES_HPP_

#include "parametersBasic.hpp"

namespace meta {
class ParametersSources {
public:

	int get_num_sources() const noexcept;
	int get_source_x() const noexcept;
	int get_source_y() const noexcept;
	int get_source_z() const noexcept;
	int get_source_t() const noexcept;
	bool get_place_sources_on_host() const noexcept;
	common::sourcetypes get_sourcetype() const noexcept;
	common::sourcecontents get_sourcecontent() const noexcept;

private:
	po::options_description options;

	int num_sources;
	int source_x;
	int source_y;
	int source_z;
	int source_t;
	bool place_sources_on_host;

protected:
	ParametersSources();
	virtual ~ParametersSources();
	ParametersSources(ParametersSources const&) = delete;
	ParametersSources & operator=(ParametersSources const&) = delete;
	po::options_description & getOptions();

	common::sourcetypes sourcetype;
	common::sourcecontents sourcecontent;
};

}

#endif
