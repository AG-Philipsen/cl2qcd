/*
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

# pragma once

#include <iostream>

void printKernelInformation(std::string name);

const int ns4 = 4;
const int nt4 = 4;
const int ns8 = 8;
const int nt8 = 8;
const int ns12 = 12;
const int nt12 = 12;
const int ns16 = 16;
const int nt16 = 16;

double sumOfIntegers(const int start, const int end, const int increment) noexcept;
double sumOfIntegersSquared(const int end) noexcept;

