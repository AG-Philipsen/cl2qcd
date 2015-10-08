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

#include "parametersSources.hpp"

int meta::ParametersSources::get_num_sources() const noexcept
{
	return num_sources;
}
int meta::ParametersSources::get_source_x() const noexcept
{
	return source_x;
}
int meta::ParametersSources::get_source_y() const noexcept
{
	return source_y;
}
int meta::ParametersSources::get_source_z() const noexcept
{
	return source_z;
}
int meta::ParametersSources::get_source_t() const noexcept
{
	return source_t;
}
bool meta::ParametersSources::get_place_sources_on_host() const noexcept
{
	return place_sources_on_host;
}

common::sourcetypes meta::ParametersSources::get_sourcetype() const noexcept
{
	return sourcetype;
}
common::sourcecontents meta::ParametersSources::get_sourcecontent() const noexcept
{
	return sourcecontent;
}

meta::ParametersSources::ParametersSources()
	: options("Source options")
{
	options.add_options()
	("sourcetype",  po::value<std::string>()->default_value("point"), "Type of source to use for inverter")
	("sourcecontent",  po::value<std::string>()->default_value("one"), "Type of content to use with inverter sources")
	("num_sources", po::value<int>(&num_sources)->default_value(12))
	("source_x", po::value<int>(&source_x)->default_value(0))
	("source_y", po::value<int>(&source_y)->default_value(0))
	("source_z", po::value<int>(&source_z)->default_value(0))
	("source_t", po::value<int>(&source_t)->default_value(0))
	("place_sources_on_host", po::value<bool>(&place_sources_on_host)->default_value(false));
}

meta::ParametersSources::~ParametersSources() = default;

po::options_description & meta::ParametersSources::getOptions()
{
	return options;
}
