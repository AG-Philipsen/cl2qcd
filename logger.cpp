/**
 * @file
 *
 * Utility file to generate singleton logging object.
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

#include "logger.hpp"

#include <string>
#include <algorithm>
#include <cctype>
#include <map>

einhard::Logger<einhard::LOG_LEVEL> logger(einhard::ALL);

std::string toLower(std::string in);

void switchLogLevel(const std::string log_level)
{
	std::map<std::string, einhard::LogLevel> log_level_map;
	for(einhard::LogLevel level = einhard::ALL; level <= einhard::OFF; level = static_cast<einhard::LogLevel>(level + 1)) {
		log_level_map[toLower(einhard::getLogLevelString(level))] = level;
	}
	logger.setVerbosity(log_level_map[toLower(log_level)]);
}

std::string toLower(std::string in)
{
	std::transform(in.begin(), in.end(), in.begin(), ::tolower);
	return in;
}
