/** @file
 * ildg IO utilities
 *
 * Copyright (c) 2014,2015 Christopher Pinke
 * Copyright (c) 2015,2016,2018 Alessandro Sciarra
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

#ifndef ILDGIO_HPP_
#define ILDGIO_HPP_

#include "../common_header_files/types.hpp"
#include "../physics/lattices/latticesInterfaces.hpp"

#include <string>
#include <vector>

namespace ildgIo {
    Matrixsu3* readGaugefieldFromSourcefile(std::string, const physics::lattices::GaugefieldParametersInterface*, int&);
    void writeGaugefieldToFile(std::string, std::vector<Matrixsu3>&,
                               const physics::lattices::GaugefieldParametersInterface*, int);
}  // namespace ildgIo

#endif
