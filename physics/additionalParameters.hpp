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
#include "../meta/inputparameters.hpp"
#include "../meta/util.hpp"

namespace physics{

    class AdditionalParameters {
        public:
            virtual ~AdditionalParameters() = 0;
            virtual hmc_float getKappa() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
            virtual hmc_float getMubar() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
            virtual hmc_float getMass() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
            virtual bool getConservative() const {throw Print_Error_Message("Generic AdditionalParameter object used!");}
    };

    inline AdditionalParameters::~AdditionalParameters()
    {
    }

    class WilsonAdditionalParameters final : public AdditionalParameters {
        public:
            WilsonAdditionalParameters() = delete;
            WilsonAdditionalParameters(const meta::Inputparameters& paramsIn, const bool withMassPreconditioning)
            : parameters(paramsIn), withMassPreconditioning(withMassPreconditioning)
            {
            }
            virtual ~WilsonAdditionalParameters()
            {
            }
            hmc_float getKappa() const override
            {
                return withMassPreconditioning ? parameters.get_kappa_mp() : parameters.get_kappa();
            }
            hmc_float getMubar() const override
            {
                return withMassPreconditioning ? meta::get_mubar_mp(parameters) : meta::get_mubar(parameters);
            }

        private:
            const meta::Inputparameters& parameters;
            bool withMassPreconditioning;
    };

    class StaggeredAdditionalParameters final : public AdditionalParameters {
        public:
            StaggeredAdditionalParameters() = delete;
            StaggeredAdditionalParameters(const meta::Inputparameters& paramsIn) : parameters(paramsIn)
            {
            }
            virtual ~StaggeredAdditionalParameters()
            {
            }
            hmc_float getMass() const override
            {
                return parameters.get_mass();
            }
            bool getConservative() const override
            {
                return parameters.get_conservative();
            }

        private:
            const meta::Inputparameters& parameters;
    };

}

