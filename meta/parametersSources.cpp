/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
 * Copyright (c) 2018 Alessandro Sciarra
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

#include "parametersSources.hpp"

#include "../executables/exceptions.hpp"

#include <boost/algorithm/string.hpp>

static common::sourcetypes translateSourceTypeToEnum(std::string);
static common::sourcecontents translateSourceContentToEnum(std::string);

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
    return translateSourceTypeToEnum(sourcetype);
}
common::sourcecontents meta::ParametersSources::get_sourcecontent() const noexcept
{
    return translateSourceContentToEnum(sourcecontent);
}

meta::ParametersSources::ParametersSources() : options("Source options")
{
    // clang-format off
	options.add_options()
	("sourceType",  po::value<std::string>(&sourcetype)->default_value("point"), "Which type of source to be used in the inverter (one among 'point', 'volume', 'timeslice', 'zslice').")
	("sourceContent",  po::value<std::string>(&sourcecontent)->default_value("one"), "Which ype of content to be used with sources in the inverter (one among 'one', 'z4', 'gaussian' and 'z2').")
	("nSources", po::value<int>(&num_sources)->default_value(12), "The number of sources to be used in the inverter.")
	("sourceX", po::value<int>(&source_x)->default_value(0), "The x coordinate for the position of a point source.")
	("sourceY", po::value<int>(&source_y)->default_value(0), "The y coordinate for the position of a point source.")
	("sourceZ", po::value<int>(&source_z)->default_value(0), "The z coordinate for the position of either a point source or a zslice source.")
	("sourceT", po::value<int>(&source_t)->default_value(0), "The t coordinate for the position of either a point source or a timeslice source.")
	("placeSourcesOnHost", po::value<bool>(&place_sources_on_host)->default_value(false), "Whether the buffer for sources is to remain on the host (this allows to avoid memory limits if many buffers are required).");
    // clang-format on
}

static common::sourcetypes translateSourceTypeToEnum(std::string s)
{
    boost::algorithm::to_lower(s);
    std::map<std::string, common::sourcetypes> m;
    m["point"]     = common::point;
    m["volume"]    = common::volume;
    m["timeslice"] = common::timeslice;
    m["zslice"]    = common::zslice;

    common::sourcetypes a = m[s];
    if (a) {  // map returns 0 if element is not found
        return a;
    } else {
        throw Invalid_Parameters("Invalid source type!", "point, volume, timeslice, zslice", s);
    }
}

static common::sourcecontents translateSourceContentToEnum(std::string s)
{
    boost::algorithm::to_lower(s);
    std::map<std::string, common::sourcecontents> m;
    m["one"]      = common::one;
    m["z4"]       = common::z4;
    m["gaussian"] = common::gaussian;
    m["z2"]       = common::z2;

    common::sourcecontents a = m[s];
    if (a) {  // map returns 0 if element is not found
        return a;
    } else {
        throw Invalid_Parameters("Invalid source content!", "one, z4, gaussian, z2", s);
    }
}
