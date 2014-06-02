/** @file
 *
 * Copyright (c) 2014 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#include "parametersTest.hpp"

double meta::ParametersTest::get_test_ref_value() const noexcept
{
  return test_ref_value;
}
double meta::ParametersTest::get_test_ref_value2() const noexcept
{
  return test_ref_value2;
}

po::options_description meta::ParametersTest::getOptions()
{
	po::options_description options("Test options");
	options.add_options()
		("test_ref_val", po::value<double>(&test_ref_value)->default_value(0.))
		("test_ref_val2", po::value<double>(&test_ref_value2)->default_value(0.));
		
	return options;
}