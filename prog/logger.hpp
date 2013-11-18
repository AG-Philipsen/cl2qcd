/**
 * @file
 *
 * Singleton for logging to screen
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

#ifndef _LOGGERHPP_
#define _LOGGERHPP_

#include "einhard/einhard.hpp"

/**
 * Singleton logger instance for logging to stdout.
 */
extern einhard::Logger<einhard::LOG_LEVEL> logger;

/**
 * Set the log level to the one described by the string.
 *
 * @arg log_level A string describing the log level. Not case sensetive, should be one of
 *                all, trace, debug, info, warn, error, fatal, off
 */
void switchLogLevel(const std::string log_level);

#endif //_LOGGERHPP_
