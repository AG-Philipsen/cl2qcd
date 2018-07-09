/** @file
 *
 * Copyright (c) 2018 Francesca Cuteri
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

#ifndef _META_PARAMETERS_INTEGRATOR_HPP_
#define _META_PARAMETERS_INTEGRATOR_HPP_

#include "parametersBasic.hpp"

namespace meta {
    class ParametersIntegrator {
      public:
        double get_tau() const noexcept;
        int get_integrationsteps(size_t timescale) const noexcept;
        int get_num_timescales() const noexcept;
        common::integrator get_integrator(size_t timescale) const noexcept;
        double get_lambda(size_t timescale) const noexcept;

      private:
        double tau;
        int integrationsteps0;
        int integrationsteps1;
        int integrationsteps2;
        int num_timescales;
        double lambda0;
        double lambda1;
        double lambda2;

      protected:
        ParametersIntegrator();
        virtual ~ParametersIntegrator()                   = default;
        ParametersIntegrator(ParametersIntegrator const&) = delete;
        ParametersIntegrator& operator=(ParametersIntegrator const&) = delete;
        void makeNeededTranslations();

        InputparametersOptions options;
        std::string integrator0String;
        std::string integrator1String;
        std::string integrator2String;
        common::integrator integrator0;
        common::integrator integrator1;
        common::integrator integrator2;
    };

}  // namespace meta

#endif
