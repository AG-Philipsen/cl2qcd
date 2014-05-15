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

KernelTester::KernelTester(std::string kernelNameIn, std::string inputfileIn, int numberOfValuesIn, int typeOfComparisionIn):
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

	if ( (typeOfComparisionIn == 1) || (typeOfComparisionIn == 2)  )
	  {
	    typeOfComparision = typeOfComparisionIn;
	  } else
	  {
	    throw( std::invalid_argument("Do not recognize type of comparision. Aborting...") );
	  }

}

#include <boost/test/floating_point_comparison.hpp>
KernelTester::~KernelTester()
{
  //NOTE: Using "require" in boost throws an exception here, which should not happen in a destructor.
	for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++) {
	  if (typeOfComparision == 1)
	    {
	      BOOST_CHECK_CLOSE(referenceValue[iteration], kernelResult[iteration], testPrecision);
	    }
	  else if (typeOfComparision == 2)
	    {
	      BOOST_CHECK_SMALL(kernelResult[iteration], referenceValue[iteration]);
	    }
	}
	delete parameters;
	delete system;
}
