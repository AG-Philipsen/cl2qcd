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

class CorrelatorsStaggeredTester : public SpinorStaggeredTester{
   public:
	CorrelatorsStaggeredTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
	     SpinorStaggeredTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues){
		
		code = device->get_correlator_staggered_code();
		sourcecontent = parameters->get_sourcecontent();
		outBuffer = new hardware::buffers::SU3vec(spinorfieldEvenOddElements, device);
		outHost = new su3vec[spinorfieldEvenOddElements * iterations];
		
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
	meta::Inputparameters::sourcecontents sourcecontent;
	
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "../correlatorStaggered/" + inputfileIn;
	}

};

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
	  BOOST_CHECK_NO_THROW(CorrelatorsStaggeredTester("build", "correlators_staggered_build_input_1", 0));
	}
	
	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
	  BOOST_CHECK_NO_THROW(CorrelatorsStaggeredTester("build", "correlators_staggered_build_input_2", 0));
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SRC_VOLUME)

	class VolumeSourceTester : public CorrelatorsStaggeredTester{
	   public:
		VolumeSourceTester(std::string inputfile) : CorrelatorsStaggeredTester("Volume_source", inputfile){
			
			hmc_float sum = 0;
			for (int i = 0; i< iterations; i++){
			  if(i%400==0)logger.info() << "Run kernel for the " << i << "th time";
			  outBuffer->clear();
			  code->create_volume_source_stagg_eoprec_device(outBuffer, prng->get_buffers().at(0));
			  outBuffer->dump(&outHost[i*spinorfieldEvenOddElements]);
			  //Here we sum the entries to calculate the mean later
			  sum += count_sf(&outHost[i*spinorfieldEvenOddElements], spinorfieldEvenOddElements);
			}
			logger.info() << "result: mean";
			//sum is the sum of iterations*spinorfieldEvenOddElements*6 real numbers
			if(sourcecontent == meta::Inputparameters::z2){
				//because immaginary part is not randomly drawn, it is 0.0 always
				sum = sum/iterations/spinorfieldEvenOddElements/3;
			}else{
				sum = sum/iterations/spinorfieldEvenOddElements/6;
			}
			if(calcVariance == false){
				kernelResult[0] = sum;
				logger.info() << sum;
			}else{
				hmc_float var=0.;
				for (int i=0; i<iterations; i++){
					var += calc_var_sf(&outHost[i*spinorfieldEvenOddElements], spinorfieldEvenOddElements, sum);
				}
				//var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
				if(sourcecontent == meta::Inputparameters::z2){
					//because immaginary part is not randomly drawn, it is 0.0 always
					var=var/iterations/spinorfieldEvenOddElements/3;
				}else{
					var=var/iterations/spinorfieldEvenOddElements/6;
				}
				kernelResult[0] = sqrt(var);
				logger.info() << "result: variance";
				logger.info() << sqrt(var);
			}
			
			if(sourcecontent == meta::Inputparameters::one ||
			  (sourcecontent == meta::Inputparameters::z2 && calcVariance)){
				typeOfComparison=1;
			}else{
				typeOfComparison=2;
			}
		}
		

	};

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_1 )
	{
	    VolumeSourceTester("/src_volume_staggered_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SRC_VOLUME_2 )
	{
	    VolumeSourceTester("/src_volume_staggered_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SRC_VOLUME_3 )
	{
	    VolumeSourceTester("/src_volume_staggered_eo_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SRC_VOLUME_4 )
	{
	    VolumeSourceTester("/src_volume_staggered_eo_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( SRC_VOLUME_5 )
	{
	    VolumeSourceTester("/src_volume_staggered_eo_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( SRC_VOLUME_6 )
	{
	    VolumeSourceTester("/src_volume_staggered_eo_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( SRC_VOLUME_7 )
	{
	    VolumeSourceTester("/src_volume_staggered_eo_input_7");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

