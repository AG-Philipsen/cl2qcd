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

#ifndef _META_PARAMETERS_GAUGE_HPP_
#define _META_PARAMETERS_GAUGE_HPP_

#include "parametersBasic.hpp"

namespace meta {
    class ParametersGauge {
      public:
        double get_beta() const noexcept;
        double get_rho() const noexcept;
        int get_rho_iter() const noexcept;
        common::action get_gaugeact() const noexcept;
        bool get_use_smearing() const noexcept;

      private:
        double beta;
        double rho;
        int rho_iter;
        bool use_smearing;

      protected:
        ParametersGauge();
        virtual ~ParametersGauge()              = default;
        ParametersGauge(ParametersGauge const&) = delete;
        ParametersGauge& operator=(ParametersGauge const&) = delete;

        InputparametersOptions options;
        common::action gaugeact;
    };

}  // namespace meta

#endif
