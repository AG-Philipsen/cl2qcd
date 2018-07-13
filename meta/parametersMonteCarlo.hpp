/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
 * Copyright (c) 2015,2017,2018 Alessandro Sciarra
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

#ifndef _META_PARAMETERS_MONTE_CARLO_HPP_
#define _META_PARAMETERS_MONTE_CARLO_HPP_

#include "parametersBasic.hpp"

namespace meta {

    class ParametersMonteCarlo {
      public:
        // Heat bath
        int get_thermalizationsteps() const noexcept;
        int get_heatbathsteps() const noexcept;
        int get_overrelaxsteps() const noexcept;
        int get_xi() const noexcept;
        bool get_use_aniso() const noexcept;
        // Hmc
        int get_hmcsteps() const noexcept;
        bool get_use_gauge_only() const noexcept;
        bool get_use_mp() const noexcept;
        // Rhmc
        double get_num_tastes() const noexcept;
        unsigned int get_num_tastes_decimal_digits() const noexcept;
        unsigned int get_num_pseudofermions() const noexcept;
        unsigned int get_rhmcsteps() const noexcept;

      private:
        int thermalizationsteps;
        // TODO: Unify the following three parameters
        int heatbathsteps;
        int hmcsteps;
        unsigned int numberOfRhmcSteps;
        int overrelaxsteps;
        int xi;
        bool use_aniso;
        bool use_gauge_only;
        bool use_mp;
        /**
         * @internal The variable num_tastes should naturally be an integer. But there is no reason that
         *           forbids to run the RHMC with a rational number of tastes. So we decide to declare it
         *           as double. Then, since we know that the power of the fermionic determinant is num_tastes/4
         *           we must then deduce the correct fraction in order to instantiate the Rational Approximation
         *           objects. This is easily done using the command line parameter num_tastes_decimal_digits,
         *           that tells how many digits after the comma are considered as valid. To avoid any numeric
         *           problem one could have dealing with big numbers, we will limit this to be not bigger than 6
         *           (fair enough limit if one thinks to physical applications).
         * @enditernal
         */
        double numberOfTastes;
        unsigned int numberOfDecimalDigitsInNumberOfTastes;
        unsigned int numberOfPseudofermions;

      protected:
        ParametersMonteCarlo();
        virtual ~ParametersMonteCarlo()                   = default;
        ParametersMonteCarlo(ParametersMonteCarlo const&) = delete;
        ParametersMonteCarlo& operator=(ParametersMonteCarlo const&) = delete;
        InputparametersOptions options;
    };

}  // namespace meta

#endif
