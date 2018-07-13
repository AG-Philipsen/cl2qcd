/** @file
 * Input file handling
 *
 * Copyright (c) 2012,2013 Matthias Bach
 * Copyright (c) 2012,2014 Christopher Pinke
 * Copyright (c) 2018 Alessandro Sciarra
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

#ifndef _META_INPUTPARAMETERS_HPP_
#define _META_INPUTPARAMETERS_HPP_

#include "../host_functionality/logger.hpp"
#include "parametersBasic.hpp"
#include "parametersConfig.hpp"
#include "parametersFermion.hpp"
#include "parametersGauge.hpp"
#include "parametersIntegrator.hpp"
#include "parametersIo.hpp"
#include "parametersMonteCarlo.hpp"
#include "parametersObs.hpp"
#include "parametersRationalApproximation.hpp"
#include "parametersSolver.hpp"
#include "parametersSources.hpp"
#include "parametersTest.hpp"

#include <boost/algorithm/string.hpp>

/**
 * This namespace contains generic utility code required by the other packages.
 */
namespace meta {

    /**
     * Parser and representation of an input file.
     *
     * This class is copyable and assignable, but should
     * be used as a const value after initialization.
     *
     * ATTENTION: This class exposes ONLY part of all existing command line options to the user, in the sense
     *            that not all parameters are considered, depending on the label parameterSet. More in detail,
     *            this class shows different options to the user depending on the executable which creates it
     *            (and which, hence, passes the "correct" label as parameterSet).
     *            Note that at the level of C++ code, this class inherits from all the parameters parents and,
     *            then, it always has all the public getters of the parents. This has the kind of drawback
     *            that the interfaces between meta and the rest of the code base may (and they do indeed) use
     *            a getter of a parameter that has not been exposed to/given by the user. In such a case the
     *            default is retrieved by the interface and is used around. Despite it might not be ideal,
     *            it is simply important to be aware of this feature. We decided it is better than giving to
     *            the user options which are not relevant to that executable and which may also be misused.
     */
    class Inputparameters : public ParametersConfig,
                            public ParametersIo,
                            public ParametersObs,
                            public ParametersGauge,
                            public ParametersFermion,
                            public ParametersSources,
                            public ParametersSolver,
                            public ParametersMonteCarlo,
                            public ParametersIntegrator,
                            public ParametersRationalApproximation,
                            public ParametersTest {
      public:
        /**
         * The parsing of the input parameters aborted for some reason.
         * Could either be invalid input or the specification of --help.
         */
        struct parse_aborted {
        };
        struct help_required {
        };

        /**
         * Construct from command line and config file.
         *
         * Config file will be retrieved from command line.
         *
         * @throws parse_aborted
         */
        Inputparameters(int argc, const char** argv, std::string parameterSet = "allParameters");

      private:
        void ChecksStringOptionsAndMapToEnum();
    };

}  // namespace meta

#endif /* _META_INPUTPARAMETERS_H_ */
