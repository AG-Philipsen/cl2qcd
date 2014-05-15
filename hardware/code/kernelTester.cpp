/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "kernelTester.hpp"
#include <boost/test/unit_test.hpp>

KernelTester::KernelTester(std::string kernelNameIn, std::string inputfileIn, int numberOfValuesIn):
	kernelResult(numberOfValuesIn, 0), referenceValue(numberOfValuesIn, 0)
{
	printKernelInformation(kernelNameIn);
	parameters = new meta::Inputparameters( createParameters(inputfileIn) );

	system = new hardware::System(*parameters);
	device = system->get_devices()[0];
	
	testPrecision = parameters->get_solver_prec();

	for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++) {
		if(iteration == 0) {
			referenceValue[iteration] = parameters->get_test_ref_value();
		} else if(iteration == 1) {
			referenceValue[iteration] = parameters->get_test_ref_value2();
		} else {
			throw( std::invalid_argument("Can only set 2 reference values at the moment. Aborting...") );
		}
	}
}

KernelTester::~KernelTester()
{
	for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++) {
		BOOST_CHECK_CLOSE(kernelResult[iteration], referenceValue[iteration], testPrecision);
	}
	delete parameters;
	delete system;
}
