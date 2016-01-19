/** @file
 * additional parameters definition
 *
 * Copyright 2016 Alessandro Sciarra
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

#pragma once

#include "../common_header_files/types.h"
#include "../executables/exceptions.h"

namespace physics{

    class AdditionalParameters {
        public:
            virtual ~AdditionalParameters() = 0;
            virtual hmc_float getKappa() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
            virtual hmc_float getMubar() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
            virtual hmc_float getMass() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
            virtual bool getConservative() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
    };

    inline AdditionalParameters::~AdditionalParameters(){}

}

