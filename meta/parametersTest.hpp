/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
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

#ifndef _META_PARAMETERS_TEST_HPP_
#define _META_PARAMETERS_TEST_HPP_

#include "parametersBasic.hpp"

namespace meta {
    class ParametersTest {
      public:
        double get_test_ref_value() const noexcept;
        double get_test_ref_value2() const noexcept;

      private:
        double test_ref_value;
        double test_ref_value2;

      protected:
        ParametersTest();
        virtual ~ParametersTest()             = default;
        ParametersTest(ParametersTest const&) = delete;
        ParametersTest& operator=(ParametersTest const&) = delete;

        InputparametersOptions options;
    };

}  // namespace meta

#endif
