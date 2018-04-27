/*
 * Copyright (c) 2014,2016 Christopher Pinke
 * Copyright (c) 2014 Francesca Cuteri
 * Copyright (c) 2014 Matthias Bach
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

#include "testUtilities.hpp"

#include "../../host_functionality/logger.hpp"

void printKernelInformation(std::string name)
{
    logger.info() << "Test kernel\t\"" << name << "\"\tagainst reference value";
}

double sumOfIntegers(const int start, const int end, const int increment) noexcept
{
    // One could also implement some variant of Faulhabers Formula here to save the loop
    double sum = 0.;
    for (int iteration = start; iteration <= end; iteration += increment) {
        sum += iteration;
    }
    return sum;
}

double sumOfIntegersSquared(const int end) noexcept
{
    return (2 * end * end * end + 3 * end * end + end) / 6.;  // Faulhaber`s formula
}
