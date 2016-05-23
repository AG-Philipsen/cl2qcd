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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_CORRELATORS_STAGGERED

#include "correlator_staggered.hpp"
#include "SpinorStaggeredTester.hpp"
#include "PrngSpinorTester.hpp"

struct StaggeredFermionsCorrelatorsTestParameters : public SpinorStaggeredTestParameters
{
	StaggeredFermionsCorrelatorsTestParameters(const LatticeExtents lE):
		TestParameters(lE), SpinorStaggeredTestParameters(lE), iterations(100), sourcecontent(common::sourcecontents::one) {};
	const int iterations;
	common::sourcecontents sourcecontent;
};

//todo: make this virtual inheritance so that calls to device are clear
struct CorrelatorsStaggeredTester : public SpinorStaggeredTester, PrngSpinorTester
{
	CorrelatorsStaggeredTester(const std::string kernelName, const ParameterCollection pC, const StaggeredFermionsCorrelatorsTestParameters tP, const ReferenceValues rV):
	     SpinorStaggeredTester(kernelName, pC, tP, rV),
	     PrngSpinorTester(kernelName, pC, PrngSpinorTestParameters(tP.latticeExtents, tP.iterations), calculateEvenOddSpinorfieldSize(tP.latticeExtents), rV),
	     elements(calculateEvenOddSpinorfieldSize(tP.latticeExtents))
	{
		code = SpinorStaggeredTester::device->getCorrelatorStaggeredCode();
		sourcecontent = tP.sourcecontent;
		outBuffer = new hardware::buffers::SU3vec(tP.latticeExtents, SpinorStaggeredTester::device);
		outHost = new su3vec[elements * tP.iterations];
	}
	
	virtual ~CorrelatorsStaggeredTester(){
		delete outBuffer;
		delete[] outHost;
		code = NULL;
	}
	
   protected:
	const hardware::code::Correlator_staggered * code;
	const hardware::buffers::SU3vec *outBuffer;
	su3vec * outHost;
	common::sourcecontents sourcecontent;
	const int elements;
};

struct VolumeSourceTester : public CorrelatorsStaggeredTester
{
	VolumeSourceTester(const ParameterCollection pC, const StaggeredFermionsCorrelatorsTestParameters tP) :
		CorrelatorsStaggeredTester("Volume_source", pC, tP, defaultReferenceValues())
	{
		hmc_float sum = 0;
		for (int i = 0; i< tP.iterations; i++){
		  if(i%400==0)logger.info() << "Run kernel for the " << i << "th time";
		  outBuffer->clear();
		  code->create_volume_source_stagg_eoprec_device(outBuffer, prngStates);
		  outBuffer->dump(&outHost[i*elements]);
		  //Here we sum the entries to calculate the mean later
		  sum += count_sf(&outHost[i*elements], elements);
		}
		logger.info() << "result: mean";
		//sum is the sum of iterations*spinorfieldEvenOddElements*6 real numbers
		if(sourcecontent == common::z2)
		{
			sum = sum/tP.iterations/elements/3; //because immaginary part is not randomly drawn, it is 0.0 always
		}else{
			sum = sum/tP.iterations/elements/6;
		}

		SpinorStaggeredTester::kernelResult[0] = sum;
		logger.info() << sum;

		hmc_float var=0.;
			for (int i=0; i<tP.iterations; i++){
				var += calc_var_sf(&outHost[i*elements], elements, sum);
			}
			//var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
			if(sourcecontent == common::z2){
				//because immaginary part is not randomly drawn, it is 0.0 always
				var=var/tP.iterations/elements/3;
			}else{
				var=var/tP.iterations/elements/6;
			}
			SpinorStaggeredTester::kernelResult[0] = sqrt(var);
			logger.info() << "result: variance";
			logger.info() << sqrt(var);
		
		//todo: introduce into testParameters
//		if(sourcecontent == common::one || sourcecontent == common::z2 )
//		{
//			SpinorStaggeredTester::typeOfComparison=1;
//		}
//		else
//		{
//			SpinorStaggeredTester::typeOfComparison=2;
//		}
	}
};

void testVolumeSource(const LatticeExtents lE)
{
	StaggeredFermionsCorrelatorsTestParameters parametersForThisTest(lE);
	//todo: Work over these!
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, true);
	hardware::code::OpenClKernelParametersMockupForSpinorStaggered kernelParameters(parametersForThisTest.latticeExtents, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	VolumeSourceTester(parameterCollection, parametersForThisTest);
}

BOOST_AUTO_TEST_SUITE(SRC_VOLUME)

	//@todo: add more tests like in "src_volume_staggered_eo_input_{1-7}"
	BOOST_AUTO_TEST_CASE( SRC_VOLUME_1 )
	{
		testVolumeSource(LatticeExtents{ns4, nt4});
	}
	
BOOST_AUTO_TEST_SUITE_END()


