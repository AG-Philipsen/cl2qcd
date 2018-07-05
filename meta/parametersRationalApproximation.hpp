/** @file
 *
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

#ifndef _META_PARAMETERS_RATIONAL_APPROXIMATION_HPP_
#define _META_PARAMETERS_RATIONAL_APPROXIMATION_HPP_

#include "parametersBasic.hpp"

namespace meta {

    class ParametersRationalApproximation {
      public:
        unsigned int get_md_approx_ord() const noexcept;
        unsigned int get_metro_approx_ord() const noexcept;
        unsigned int get_findminmax_iteration_block_size() const noexcept;
        unsigned int get_findminmax_max() const noexcept;
        double get_findminmax_prec() const noexcept;
        bool get_conservative() const noexcept;
        double get_approx_lower() const noexcept;
        double get_approx_upper() const noexcept;
        std::string get_approx_heatbath_file() const noexcept;
        std::string get_approx_md_file() const noexcept;
        std::string get_approx_metropolis_file() const noexcept;
        bool get_read_rational_approximations_from_file() const noexcept;

      protected:
        /** @TODO If the rational approximation is read from file than its parameters could differ
         *        from the following! This means, for example, that one could use get_md_approx_ord()
         *        to get a value that is not that loaded from the file!
         *  @TODO If read_rational_approximations_from_file is false it makes no sense to have the
         *        approx_*_file variables, but this is similar to the gauge configuration.
         */
        unsigned int molecularDynamicsRationalApproximationOrder;
        unsigned int metropolisTestRationalApproximationOrder;
        unsigned int findMinMaxEigenvalueIterationBlockSize;
        unsigned int findMinMaxEigenvalueMaxNumberOfIterations;
        double findMinMaxEigenvaluePrecision;
        bool beConservativeInFindMinMaxEigenvalue;  // this is for the strategy in findminmax_eigenvalues
        double lowerBoundForRationalApproximationRange;
        double upperBoundForRationalApproximationRange;  // range of validity of the Rational Approximation
        bool readRationalApproximationsFromFile;
        std::string heatbathRationalApproximationFilename;
        std::string molecularDynamicsRationalApproximationFilename;
        std::string metropolisTestRationalApproximationFilename;

      protected:
        ParametersRationalApproximation();
        virtual ~ParametersRationalApproximation()                              = default;
        ParametersRationalApproximation(ParametersRationalApproximation const&) = delete;
        ParametersRationalApproximation& operator=(ParametersRationalApproximation const&) = delete;

        InputparametersOptions options;
    };

}  // namespace meta

#endif
