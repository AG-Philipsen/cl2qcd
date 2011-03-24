/**
 * @file
 *
 * This is the main include file for Klepsydra.
 *
 * Copyright 2011 Matthias Bach <marix@marix.org>
 *
 * This file is part of Klepsydra.
 *
 * Klepsydra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Klepsydra is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Klepsydra.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \mainpage Klepsydra Timing Library
 *
 * \section intro_sec Introduction
 *
 * Klepsydra aims at being a lightweight timing library with the following features:
 *  - Timing using monotonic timers
 *  - Timing on cpu clock level
 * Other features might be added and not all features will be added from the start.
 *
 * To use Klepsydra all you need to do is create an klepsydra::Timer object and measure the time.
 * streams.
 *
 * \code
 * klepsydra::MonotonicTimer timer;
 * // do foo
 * std::cout << "Elapsed Time: " << timer.getTime() << " mus." << std::endl;
 * \endcode
 *
 * \section install_sec Installation
 *
 * Klepsydra is build using cmake. You can install it using the usual cmake triplet:
 *
 * \code
 * cmake
 * make
 * sudo make install
 * \endcode
 *
 * If you want to build a static lib or install to a custom path you can use the usual cmake
 * utilities to adjust the configuration.
 */

// only import all headers

#include "timer.hpp"
#include "monotonic.hpp"

