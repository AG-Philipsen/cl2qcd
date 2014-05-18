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

// #include "../meta/util.hpp"
// #include "../host_functionality/host_random.h"
// #include "../physics/prng.hpp"
// #include "../hardware/device.hpp"
// #include "../hardware/code/spinors_staggered.hpp"
// //spinors.hpp needed for get_spinorfieldsize and get_eoprec_spinorfieldsize
// #include "../hardware/code/spinors.hpp" 
// #include <numeric>      // std::accumulate
/*
// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_SPINORS_STAGGERED
#include <boost/test/unit_test.hpp>*/
// 
// //some functionality
// #include "../../tests/test_util.h"
// #include "../../tests/test_util_staggered.h"


// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_SPINORS_STAGGERED

#include "SpinorStaggeredTester.hpp"
#include "Kolmogorov_Smirnov.h"
#include "Normal_RNG_tests.h"

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
	   BOOST_CHECK_NO_THROW(SpinorStaggeredTester spinorStaggeredTester("build all kernels",
									     "spinors_staggered_build_input_1", 0));
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
	   BOOST_CHECK_NO_THROW(SpinorStaggeredTester spinorStaggeredTester("build all kernels", 
									      "spinors_staggered_build_input_2", 0));
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SQUARENORM)

	class SquarenormTester: public SpinorStaggeredTester{
	  public:
		SquarenormTester(std::string inputfile) : SpinorStaggeredTester("squarenorm", inputfile, 1){
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			in.load(createSpinorfield(spinorfieldElements));
			calcSquarenormAndStoreAsKernelResult(&in);
			
	/*
        print_staggeredfield_to_textfile("ref_vec_sq", createSpinorfield(spinorfieldElements)); 
        logger.info() << "Produced the ref_vec_sq text file with the staggered field for the ref. code.";
	*/
		}
	};

	BOOST_AUTO_TEST_CASE( SQUARENORM_1 )
	{
	   SquarenormTester("squarenorm_input_1");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_2 )
	{
	   SquarenormTester("squarenorm_input_2");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_1 )
	{
	   SquarenormTester("squarenorm_reduction_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_2 )
	{
	   SquarenormTester("squarenorm_reduction_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_REDUCTION_3 )
	{
	   SquarenormTester("squarenorm_reduction_input_3");
	}	
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT)

	class ScalarProductTester: public SpinorStaggeredTester{
	   public:
		ScalarProductTester(std::string inputfile) : SpinorStaggeredTester("scalar_product", inputfile, 2){
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			const hardware::buffers::Plain<su3vec> in2(spinorfieldElements, device);
			in.load(createSpinorfield(spinorfieldElements, 123));
			in2.load(createSpinorfield(spinorfieldElements, 456));
			hardware::buffers::Plain<hmc_complex> result(1, device);
			
			code->set_complex_to_scalar_product_device(&in, &in2, &result);
			hmc_complex resultHost;
			result.dump(&resultHost);

			kernelResult[0] = resultHost.re;
			kernelResult[1] = resultHost.im;
	
	/*
	print_staggeredfield_to_textfile("ref_vec_sp1", createSpinorfield(spinorfieldElements, 123)); 
        logger.info() << "Produced the ref_vec_sp1 text file with the staggered field for the ref. code.";   
        print_staggeredfield_to_textfile("ref_vec_sp2", createSpinorfield(spinorfieldElements, 456)); 
        logger.info() << "Produced the ref_vec_sp2 text file with the staggered field for the ref. code."; 
	*/
			}
	};

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_1 )
	{
	    ScalarProductTester("scalar_product_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_2 )
	{
	    ScalarProductTester("scalar_product_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_1 )
	{
	    ScalarProductTester("scalar_product_reduction_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_2 )
	{
	    ScalarProductTester("scalar_product_reduction_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_3 )
	{
	    ScalarProductTester("scalar_product_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(COLD_AND_ZERO)

	class ColdAndZeroTester: public SpinorStaggeredTester{
	   public:
		ColdAndZeroTester(std::string inputfile, bool switcher) : SpinorStaggeredTester("cold or zero", inputfile)
			{
				const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
				in.load(createSpinorfield(spinorfieldElements));
				(switcher) ? code->set_cold_spinorfield_device(&in) : 	code->set_zero_spinorfield_device(&in);
				calcSquarenormAndStoreAsKernelResult(&in);
			}
	};

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
	    ColdAndZeroTester("set_zero_input_1", false);
	}
	
	BOOST_AUTO_TEST_CASE( COLD_1 )
	{
	    ColdAndZeroTester("set_cold_input_1", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAX)

	class SaxTester: public SpinorStaggeredTester{
	   public:
		SaxTester(std::string inputfile) : SpinorStaggeredTester("sax", inputfile){
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			const hardware::buffers::Plain<su3vec> out(spinorfieldElements, device);
			hardware::buffers::Plain<hmc_complex> alpha(1, device);

			in.load(createSpinorfield(spinorfieldElements, 123));
			alpha.load(&alpha_host);

			code->sax_device(&in, &alpha, &out);
			calcSquarenormAndStoreAsKernelResult(&out);

        /*
        print_staggeredfield_to_textfile("ref_vec_sax", createSpinorfield(spinorfieldElements, 123)); 
        logger.info() << "Produced the ref_vec_sax text file with the staggered field for the ref. code.";
	*/
		}
	};

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
	    SaxTester("sax_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
	    SaxTester("sax_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
	    SaxTester("sax_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
	    SaxTester("sax_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_5 )
	{
	    SaxTester("sax_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_6 )
	{
	    SaxTester("sax_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_7 )
	{
	    SaxTester("sax_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_8 )
	{
	    SaxTester("sax_input_8");
	}
	
	BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPY)

	class SaxpyTester: public SpinorStaggeredTester{
	   public:
		SaxpyTester(std::string inputfile, bool switcher) : SpinorStaggeredTester("saxpy", inputfile){
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			const hardware::buffers::Plain<su3vec> in2(spinorfieldElements, device);
			const hardware::buffers::Plain<su3vec> out(spinorfieldElements, device);
			hardware::buffers::Plain<hmc_complex> alpha(1, device);

			in.load(createSpinorfield(spinorfieldElements, 123));
			in2.load(createSpinorfield(spinorfieldElements, 456));
			alpha.load(&alpha_host);

			(switcher) ? code->saxpy_device(&in, &in2, &alpha, &out) : (throw(std::invalid_argument("The kernel saxpy_stagg_arg doesn't exist yet, this is a meaningless test!")));
			calcSquarenormAndStoreAsKernelResult(&out);
	
	/*
        print_staggeredfield_to_textfile("ref_vec_saxpy1", createSpinorfield(spinorfieldElements, 123)); 
        logger.info() << "Produced the ref_vec_saxpy1 text file with the staggered field for the ref. code.";   
        print_staggeredfield_to_textfile("ref_vec_saxpy2", createSpinorfield(spinorfieldElements, 456)); 
        logger.info() << "Produced the ref_vec_saxpy2 text file with the staggered field for the ref. code."; 
	*/
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
	    SaxpyTester("saxpy_input_1", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
	    SaxpyTester("saxpy_input_2", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
	    SaxpyTester("saxpy_input_3", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
	    SaxpyTester("saxpy_input_4", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_5 )
	{
	    SaxpyTester("saxpy_input_5", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_6 )
	{
	    SaxpyTester("saxpy_input_6", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_7 )
	{
	    SaxpyTester("saxpy_input_7", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_8 )
	{
	    SaxpyTester("saxpy_input_8", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_9 )
	{
	    SaxpyTester("saxpy_input_9", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_10 )
	{
	    SaxpyTester("saxpy_input_10", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_11 )
	{
	    SaxpyTester("saxpy_input_11", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_12 )
	{
	    SaxpyTester("saxpy_input_12", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_13 )
	{
	    SaxpyTester("saxpy_input_13", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_14 )
	{
	    SaxpyTester("saxpy_input_14", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_15 )
	{
	    SaxpyTester("saxpy_input_15", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_16 )
	{
	    SaxpyTester("saxpy_input_16", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_17 )
	{
	    SaxpyTester("saxpy_input_17", true);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_18 )
	{
	    SaxpyTester("saxpy_input_18", true);
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBYPZ)

	class SaxpbypzTester: public SpinorStaggeredTester{
	   public:
		SaxpbypzTester(std::string inputfile) : SpinorStaggeredTester("saxpbypz", inputfile){
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			const hardware::buffers::Plain<su3vec> in2(spinorfieldElements, device);
			const hardware::buffers::Plain<su3vec> in3(spinorfieldElements, device);
			const hardware::buffers::Plain<su3vec> out(spinorfieldElements, device);
			hardware::buffers::Plain<hmc_complex> alpha(1, device);
			hardware::buffers::Plain<hmc_complex> beta(1, device);

			in.load(createSpinorfield(spinorfieldElements, 123));
			in2.load(createSpinorfield(spinorfieldElements, 456));
			in3.load(createSpinorfield(spinorfieldElements, 789));
			alpha.load(&alpha_host);
			beta.load(&beta_host);

			code->saxpbypz_device(&in, &in2, &in3, &alpha, &beta, &out);
			calcSquarenormAndStoreAsKernelResult(&out);
			
	/*
        print_staggeredfield_to_textfile("ref_vec_saxpbypz1", createSpinorfield(spinorfieldElements, 123)); 
        logger.info() << "Produced the ref_vec_saxpbypz1 text file with the staggered field for the ref. code."; 
	print_staggeredfield_to_textfile("ref_vec_saxpbypz2", createSpinorfield(spinorfieldElements, 456)); 
        logger.info() << "Produced the ref_vec_saxpbypz2 text file with the staggered field for the ref. code.";  
        print_staggeredfield_to_textfile("ref_vec_saxpbypz3", createSpinorfield(spinorfieldElements, 789)); 
        logger.info() << "Produced the ref_vec_saxpbypz3 text file with the staggered field for the ref. code.";
	*/
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_1 )
	{
	    SaxpbypzTester("/saxpbypz_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_2 )
	{
	    SaxpbypzTester("/saxpbypz_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_3 )
	{
	    SaxpbypzTester("/saxpbypz_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_4 )
	{
	    SaxpbypzTester("/saxpbypz_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_5 )
	{
	    SaxpbypzTester("/saxpbypz_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_6 )
	{
	    SaxpbypzTester("/saxpbypz_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_7 )
	{
	    SaxpbypzTester("/saxpbypz_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_8 )
	{
	    SaxpbypzTester("/saxpbypz_input_8");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_9 )
	{
	    SaxpbypzTester("/saxpbypz_input_9");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_10 )
	{
	    SaxpbypzTester("/saxpbypz_input_10");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_11 )
	{
	    SaxpbypzTester("/saxpbypz_input_11");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_12 )
	{
	    SaxpbypzTester("/saxpbypz_input_12");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_13 )
	{
	    SaxpbypzTester("/saxpbypz_input_13");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_14 )
	{
	    SaxpbypzTester("/saxpbypz_input_14");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_15 )
	{
	    SaxpbypzTester("/saxpbypz_input_15");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_16 )
	{
	    SaxpbypzTester("/saxpbypz_input_16");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_17 )
	{
	    SaxpbypzTester("/saxpbypz_input_17");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_18 )
	{
	    SaxpbypzTester("/saxpbypz_input_18");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_19 )
	{
	    SaxpbypzTester("/saxpbypz_input_19");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_20 )
	{
	    SaxpbypzTester("/saxpbypz_input_20");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_21 )
	{
	    SaxpbypzTester("/saxpbypz_input_21");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_22 )
	{
	    SaxpbypzTester("/saxpbypz_input_22");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_23 )
	{
	    SaxpbypzTester("/saxpbypz_input_23");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_24 )
	{
	    SaxpbypzTester("/saxpbypz_input_24");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_25 )
	{
	    SaxpbypzTester("/saxpbypz_input_25");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_26 )
	{
	    SaxpbypzTester("/saxpbypz_input_26");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_27 )
	{
	    SaxpbypzTester("/saxpbypz_input_27");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_28 )
	{
	    SaxpbypzTester("/saxpbypz_input_28");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_29 )
	{
	    SaxpbypzTester("/saxpbypz_input_29");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_30 )
	{
	    SaxpbypzTester("/saxpbypz_input_30");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_31 )
	{
	    SaxpbypzTester("/saxpbypz_input_31");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_32 )
	{
	    SaxpbypzTester("/saxpbypz_input_32");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_33 )
	{
	    SaxpbypzTester("/saxpbypz_input_33");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_34 )
	{
	    SaxpbypzTester("/saxpbypz_input_34");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	class GaussianTester: public SpinorStaggeredTester{
	   public:
		GaussianTester(std::string inputfile) : SpinorStaggeredTester("gaussian_spinorfield", inputfile, 1, 2){
			const hardware::buffers::Plain<su3vec> out(spinorfieldElements, device);
			hardware::buffers::Plain<hmc_float> sqnorm(1, device);

			su3vec * outHost;
			outHost = new su3vec[spinorfieldElements * iterations];
			BOOST_REQUIRE(out);
				
			auto prng_buf = prng->get_buffers().at(0);
		
			double sum = 0;
			for (int i = 0; i < iterations; i++) {
				if(i%100==0) logger.info() << "Run kernel for the " << i << "th time";
				code->set_gaussian_spinorfield_device(&out, prng_buf);
				out.dump(&outHost[i * spinorfieldElements]);
				sum += count_sf(&outHost[i * spinorfieldElements], spinorfieldElements);
			}
			//sum is the sum of iterations*NUM_ELEMENTS_SF*6 real gaussian numbers
			sum /= (iterations * spinorfieldElements * 6);
			kernelResult[0] = sum;

			if(calcVariance){
				double var = 0.;
				for (int i = 0; i < iterations; i++) {
				   var += calc_var_sf(&outHost[i * spinorfieldElements], spinorfieldElements, sum);
				}
				//var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
				var /= (iterations * spinorfieldElements * 6);
				kernelResult[0] = sqrt(var);
			}
			
			/** 
			 * @TODO This piece of code contains actually tests for the RNG itself 
			 *       and should be moved elsewhere. 
			 */
			//The Kolmogorov_Smirnov test requires set of n samples with n around 1000
			//(to big n and to small n are not good choices for this test)
			//So we can consider each kernel result as a set (NUM_ELEMENTS_SF*6=1536 for a 4^4 lattice)
			vector<vector<hmc_float>> samples;
			vector<hmc_float> tmp;
			vector<hmc_float> tmp2;
			for(int i=0; i<iterations; i++){
			  vector<hmc_float> tmp;
			  for(uint j=0; j<spinorfieldElements; j++){
			    tmp2=reals_from_su3vec(outHost[i*spinorfieldElements+j]);
			    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
			    tmp2.clear();
			  }
			  samples.push_back(tmp);
			  tmp.clear();
			}
			logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
			logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";
			
			if(!calcVariance){
			  //Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
			}else{
			  //Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  variance_test_multiple_set(samples,2.,sqrt(0.5));
			  variance_test_multiple_set(samples,3.,sqrt(0.5));
			  variance_test_multiple_set(samples,4.,sqrt(0.5));
			}
		}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_1 )
	{
	    GaussianTester("gaussian_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_2 )
	{
	    GaussianTester("gaussian_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_3 )
	{
	    GaussianTester("gaussian_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_4 )
	{
	    GaussianTester("gaussian_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////           EVEN-ODD PRECONDITIONING             //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE(CONVERT_EO)

	class ConvertEvenOddTester: public SpinorStaggeredTester{
	   public:
		ConvertEvenOddTester(std::string inputfile, bool from_eo) : SpinorStaggeredTester("convert fro/to eo", inputfile, 2){
			const hardware::buffers::Plain<su3vec> in(spinorfieldElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in3(spinorfieldEvenOddElements, device);

			if(from_eo){
				in2.load(createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity());
				in3.load(createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity());
				code->convert_from_eoprec_device(&in2, &in3, &in);
				code->set_float_to_global_squarenorm_device(&in, doubleBuffer);
				doubleBuffer->dump(&kernelResult[0]);
				referenceValue[0] = spinorfieldEvenOddElements*3;
				referenceValue[1] = 0.; //this has not actually a meaning
			}else{
				in.load(createSpinorfieldWithOnesAndZerosDependingOnSiteParity());
				code->convert_to_eoprec_device(&in2, &in3, &in);
				code->set_float_to_global_squarenorm_eoprec_device(&in2, doubleBuffer);
				doubleBuffer->dump(&kernelResult[0]);
				code->set_float_to_global_squarenorm_eoprec_device(&in3, doubleBuffer);
				doubleBuffer->dump(&kernelResult[1]);
				
				if (evenOrOdd){
				  referenceValue[0] = spinorfieldEvenOddElements*3; 
				  referenceValue[1] = 0.;
				}else{ 
				  referenceValue[1] = spinorfieldEvenOddElements*3; 
				  referenceValue[0] = 0.;
				}
			}
		}
	};


	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
	    ConvertEvenOddTester("convert_eo_input_1", true);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
	{
	    ConvertEvenOddTester("convert_eo_input_2", true);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
	    ConvertEvenOddTester("convert_eo_input_1", false);
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
	{
	    ConvertEvenOddTester("convert_eo_input_2", false);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SQUARENORM_EO)

	class SquarenormEvenOddTester: public SpinorStaggeredTester{
	   public:
		SquarenormEvenOddTester(std::string inputfile) : SpinorStaggeredTester("squarenorm_eo", inputfile){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements));
			calcSquarenormEvenOddAndStoreAsKernelResult(&in);
			
        /*
        print_staggeredfield_eo_to_textfile("ref_vec_sq_eo", createSpinorfield(spinorfieldEvenOddElements)); 
        logger.info() << "Produced the ref_vec_sq_eo text file with the staggered field for the ref. code.";
        */
		}
	};

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_1 )
	{
	    SquarenormEvenOddTester("squarenorm_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_2 )
	{
	    SquarenormEvenOddTester("squarenorm_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_1 )
	{
	    SquarenormEvenOddTester("squarenorm_eo_reduction_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_2 )
	{
	    SquarenormEvenOddTester("squarenorm_eo_reduction_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_3 )
	{
	    SquarenormEvenOddTester("squarenorm_eo_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT_EO)
      
	class ScalarProductEvenOddTester: public SpinorStaggeredTester{
	   public:
		ScalarProductEvenOddTester(std::string inputfile, bool real=false) :SpinorStaggeredTester("scalar product eo", inputfile, 2){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));

			if(!real){
				hardware::buffers::Plain<hmc_complex> sqnorm(1, device);
				code->set_complex_to_scalar_product_eoprec_device(&in, &in2, &sqnorm);
				hmc_complex resultTmp;
				sqnorm.dump(&resultTmp);
				kernelResult[0] = resultTmp.re;
				kernelResult[1] = resultTmp.im;
			}else{
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);
				code->set_float_to_scalar_product_real_part_eoprec_device(&in, &in2, &sqnorm);
				hmc_float resultTmp;
				sqnorm.dump(&resultTmp);
				kernelResult[0] = resultTmp;
				kernelResult[1] = 0;
			}
			
        /*
        print_staggeredfield_eo_to_textfile("ref_vec_sp1", createSpinorfield(spinorfieldEvenOddElements, 123)); 
        logger.info() << "Produced the ref_vec_sp1 text file with the staggered field for the ref. code.";   
        print_staggeredfield_eo_to_textfile("ref_vec_sp2", createSpinorfield(spinorfieldEvenOddElements, 456)); 
        logger.info() << "Produced the ref_vec_sp2 text file with the staggered field for the ref. code.";
	*/
			}
	};

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_1 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_2 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_1 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_input_1", true);
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_2 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_reduction_input_1", true);
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_3 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_reduction_input_2", true);
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REAL_EO_4 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_reduction_input_3", true);
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_1 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_reduction_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_2 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_reduction_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_3 )
	{
	    ScalarProductEvenOddTester("scalar_product_eo_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(COLD_AND_ZERO_EO)

	class ColdAndZeroEvenOddTester: public SpinorStaggeredTester{
	   public:
		ColdAndZeroEvenOddTester(std::string inputfile, bool switcher): SpinorStaggeredTester("cold or zero eo", inputfile){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements));
			(switcher) ? code->set_cold_spinorfield_eoprec_device(&in) : 
			             code->set_zero_spinorfield_eoprec_device(&in);
			calcSquarenormEvenOddAndStoreAsKernelResult(&in);
		}
	};

	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
	    ColdAndZeroEvenOddTester("set_cold_eo_input_1", true);
	}
	
	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
	    ColdAndZeroEvenOddTester("set_zero_eo_input_1",  false);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAX_EO)

	class SaxEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxEvenOddTester(std::string inputfile, int switcher):SpinorStaggeredTester("sax_eo",inputfile, 1){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			
			if(switcher==0 || switcher==1){
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				alpha.load(&alpha_host);
				if(switcher==0)
				   code->sax_eoprec_device(&in, &alpha, &out);
				if(switcher==1)
				   code->sax_eoprec_device(&in, alpha_host, &out);
			}else if(switcher==2 || switcher==3){
				hardware::buffers::Plain<hmc_float> alpha_real(1, device);
				alpha_real.load(&alpha_host.re);
				if(switcher==2)
				  code->sax_eoprec_device(&in, &alpha_real, &out);
				if(switcher==3)
				  code->sax_eoprec_device(&in, alpha_host.re, &out);
			}else{
				hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
				std::vector<hmc_float> alpha_host_real_vec(5, alpha_host.re);
				const int index_alpha = 3;
				alpha_real_vec.load(&alpha_host_real_vec[0]); 
				code->sax_eoprec_device(&in, &alpha_real_vec, index_alpha, &out);
			}
			
			calcSquarenormEvenOddAndStoreAsKernelResult(&out);
		
	/*
        print_staggeredfield_eo_to_textfile("ref_vec_sax_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
        logger.info() << "Produced the ref_vec_sax_eo text file with the staggered field for the ref. code.";
	*/
		}
	};

	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_1 )
	{
	    SaxEvenOddTester("sax_eo_input_1", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_2 )
	{
	    SaxEvenOddTester("sax_eo_input_2", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_3 )
	{
	    SaxEvenOddTester("sax_eo_input_3", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_4 )
	{
	    SaxEvenOddTester("sax_eo_input_4", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_5 )
	{
	    SaxEvenOddTester("sax_eo_input_5", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_6 )
	{
	    SaxEvenOddTester("sax_eo_input_6", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_7 )
	{
	    SaxEvenOddTester("sax_eo_input_7", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_EO_8 )
	{
	    SaxEvenOddTester("sax_eo_input_8", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_1 )
	{
	    SaxEvenOddTester("sax_eo_input_1", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_2 )
	{
	    SaxEvenOddTester("sax_eo_input_2", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_3 )
	{
	    SaxEvenOddTester("sax_eo_input_3", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_4 )
	{
	    SaxEvenOddTester("sax_eo_input_4", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_5 )
	{
	    SaxEvenOddTester("sax_eo_input_5", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_6 )
	{
	    SaxEvenOddTester("sax_eo_input_6", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_7 )
	{
	    SaxEvenOddTester("sax_eo_input_7", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_CPLX_ARG_EO_8 )
	{
	    SaxEvenOddTester("sax_eo_input_8", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_1 )
	{
	    SaxEvenOddTester("sax_eo_input_1", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_2 )
	{
	    SaxEvenOddTester("sax_eo_input_2", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_3 )
	{
	    SaxEvenOddTester("sax_eo_input_5", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_EO_4 )
	{
	    SaxEvenOddTester("sax_eo_input_6", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_ARG_EO_1 )
	{
	    SaxEvenOddTester("sax_eo_input_1", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_ARG_EO_2 )
	{
	    SaxEvenOddTester("sax_eo_input_2", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_ARG_EO_3 )
	{
	    SaxEvenOddTester("sax_eo_input_5", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_ARG_EO_4 )
	{
	    SaxEvenOddTester("sax_eo_input_6", 3);
	}

	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_1 )
	{
	    SaxEvenOddTester("sax_eo_input_1", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_2 )
	{
	    SaxEvenOddTester("sax_eo_input_2", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_3 )
	{
	    SaxEvenOddTester("sax_eo_input_5", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAX_REAL_VEC_EO_4 )
	{
	    SaxEvenOddTester("sax_eo_input_6", 4);
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPY_EO)

	class SaxpyEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxpyEvenOddTester(std::string inputfile, int switcher):SpinorStaggeredTester("saxpy_eo or saxpy_arg_eo", inputfile, 1){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
			
			if(switcher==0 || switcher==1){
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				alpha.load(&alpha_host);
				if(switcher==0)
				    code->saxpy_eoprec_device(&in, &in2, &alpha, &out);
				if(switcher==1)
				    code->saxpy_eoprec_device(&in, &in2, alpha_host, &out);
			}else if(switcher==2 || switcher==3){
				hardware::buffers::Plain<hmc_float> alpha_real(1, device);
				alpha_real.load(&alpha_host.re);
				if(switcher==2)
				    code->saxpy_eoprec_device(&in, &in2, &alpha_real, &out);
				if(switcher==3)
				    code->saxpy_eoprec_device(&in, &in2, alpha_host.re, &out);
			}else{
				hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
				std::vector<hmc_float> alpha_host_real_vec(5, alpha_host.re);
				const int index_alpha = 3;
				alpha_real_vec.load(&alpha_host_real_vec[0]); 
				code->saxpy_eoprec_device(&in,  &in2, &alpha_real_vec, index_alpha, &out);
			}
				
			calcSquarenormEvenOddAndStoreAsKernelResult(&out);

     /*
     print_staggeredfield_eo_to_textfile("ref_vec_saxpy1_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
     logger.info() << "Produced the ref_vec_saxpy1_eo text file with the staggered field for the ref. code.";   
     print_staggeredfield_eo_to_textfile("ref_vec_saxpy2_eo", createSpinorfield(spinorfieldEvenOddElements, 456)); 
     logger.info() << "Produced the ref_vec_saxpy2_eo text file with the staggered field for the ref. code.";
     */
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_3", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_4", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_5", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_7 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_7", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_8 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_8", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_9 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_9", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_10 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_11 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_12 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_12", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_13 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_13", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_14 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_14", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_15 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_16 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_16", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_17 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_17", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_EO_18 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_18", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_3", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_4", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_5", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_7 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_7", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_8 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_8", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_9 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_9", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_10 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_11 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_12 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_12", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_13 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_13", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_14 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_14", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_15 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_16 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_16", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_17 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_17", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_CPLX_ARG_EO_18 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_18", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_3 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_6", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_4 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_10", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_5 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_11", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_ARG_EO_6 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_15", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_1 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_1", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_2 )
	{
	    SaxpyEvenOddTester("saxpy_eo_input_2", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_3 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_6", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_4 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_10", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_5 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_11", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPY_REAL_VEC_EO_6 )
	{
	      SaxpyEvenOddTester("saxpy_eo_input_15", 4);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBY_EO)

	class SaxpbyEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxpbyEvenOddTester(std::string inputfile, int switcher):SpinorStaggeredTester("saxpy_eo or saxpy_arg_eo", inputfile, 1){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
			
			if(switcher==0 || switcher==1){
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				hardware::buffers::Plain<hmc_complex> beta(1, device);
				alpha.load(&alpha_host);
				beta.load(&beta_host);
				if(switcher==0)
				    code->saxpby_eoprec_device(&in, &in2, &alpha, &beta, &out);
				if(switcher==1)
				    code->saxpby_eoprec_device(&in, &in2, alpha_host, beta_host, &out);
			}else if(switcher==2 || switcher==3){
				hardware::buffers::Plain<hmc_float> alpha_real(1, device);
				hardware::buffers::Plain<hmc_float> beta_real(1, device);
				alpha_real.load(&alpha_host.re);
				beta_real.load(&beta_host.re);
				if(switcher==2)
				    code->saxpby_eoprec_device(&in, &in2, &alpha_real, &beta_real, &out);
				if(switcher==3)
				    code->saxpby_eoprec_device(&in, &in2, alpha_host.re, beta_host.re, &out);
			}else{
				hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device);
				std::vector<hmc_float> alpha_host_real_vec(5, alpha_host.re);
				const int index_alpha = 3;
				hardware::buffers::Plain<hmc_float> beta_real_vec(5, device);
				std::vector<hmc_float> beta_host_real_vec(5, beta_host.re);
				const int index_beta = 2;
				alpha_real_vec.load(&alpha_host_real_vec[0]); 
				beta_real_vec.load(&beta_host_real_vec[0]); 
				code->saxpby_eoprec_device(&in, &in2, &alpha_real_vec, &beta_real_vec, index_alpha, index_beta, &out);
			}
			
			calcSquarenormEvenOddAndStoreAsKernelResult(&out);

     /*
     print_staggeredfield_eo_to_textfile("ref_vec_saxpby1_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
     logger.info() << "Produced the ref_vec_saxpby1_eo text file with the staggered field for the ref. code.";   
     print_staggeredfield_eo_to_textfile("ref_vec_saxpby2_eo", createSpinorfield(spinorfieldEvenOddElements, 456)); 
     logger.info() << "Produced the ref_vec_saxpby2_eo text file with the staggered field for the ref. code.";
     */
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPBY_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_2", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_4", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_5", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_6", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_7", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_8", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_9 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_10 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_10", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_11 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_12 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_12", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_13 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_13", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_14 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_14", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_15 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_15", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_16 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_16", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_17 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_17", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_18 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_19 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_19", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_20 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_21 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_21", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_22 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_22", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_23 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_23", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_24 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_24", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_25 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_25", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_26 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_27 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_27", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_28 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_29 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_29", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_30 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_30", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_31 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_31", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_32 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_32", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_33 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_33", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_EO_34 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_34", 0);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_2", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_4", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_5", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_6", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_7", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_8", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_9 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_10 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_10", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_11 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_12 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_12", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_13 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_13", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_14 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_14", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_15 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_15", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_16 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_16", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_17 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_17", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_18 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_19 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_19", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_20 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_21 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_21", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_22 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_22", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_23 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_23", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_24 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_24", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_25 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_25", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_26 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_27 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_27", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_28 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_29 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_29", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_30 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_30", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_31 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_31", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_32 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_32", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_33 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_33", 1);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_ARG_EO_34 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_34", 1);
	}

	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 2);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_ARG_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 3);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_1 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_1", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_2 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_3", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_3 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_9", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_4 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_11", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_5 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_18", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_6 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_20", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_7 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_26", 4);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBY_REAL_VEC_EO_8 )
	{
	    SaxpbyEvenOddTester("saxpby_eo_input_28", 4);
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAXPBYPZ_EO)

	class SaxpbypzEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxpbypzEvenOddTester(std::string inputfile, bool switcher=true) : 
		                               SpinorStaggeredTester("saxsbypz_eo", inputfile){
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec in3(spinorfieldEvenOddElements, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			hardware::buffers::Plain<hmc_complex> alpha(1, device);
			hardware::buffers::Plain<hmc_complex> beta(1, device);

			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			in2.load(createSpinorfield(spinorfieldEvenOddElements, 456));
			in3.load(createSpinorfield(spinorfieldEvenOddElements, 789));
			alpha.load(&alpha_host);
			beta.load(&beta_host);

			(switcher) ? code->saxpbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out) :
			             code->saxpbypz_eoprec_device(&in, &in2, &in3, alpha_host, beta_host, &out);

			calcSquarenormEvenOddAndStoreAsKernelResult(&out);
			
  /*
  print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz1_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
  logger.info() << "Produced the ref_vec_saxpbypz1_eo text file with the staggered field for the ref. code."; 
  print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz2_eo", createSpinorfield(spinorfieldEvenOddElements, 456)); 
  logger.info() << "Produced the ref_vec_saxpbypz2_eo text file with the staggered field for the ref. code.";  
  print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz3_eo", createSpinorfield(spinorfieldEvenOddElements, 789)); 
  logger.info() << "Produced the ref_vec_saxpbypz3_eo text file with the staggered field for the ref. code.";
  */
		}
	};

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_1 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_2 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_3 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_4 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_5 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_6 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_7 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_8 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_8");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_9 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_9");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_10 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_10");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_11 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_11");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_12 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_12");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_13 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_13");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_14 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_14");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_15 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_15");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_16 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_16");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_17 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_17");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_18 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_18");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_19 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_19");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_20 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_20");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_21 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_21");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_22 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_22");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_23 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_23");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_24 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_24");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_25 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_25");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_26 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_26");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_27 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_27");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_28 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_28");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_29 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_29");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_30 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_30");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_31 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_31");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_32 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_32");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_33 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_33");
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_EO_34 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_34");
	}

	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_1 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_1", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_2 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_2", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_3 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_3", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_4 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_4", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_5 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_5", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_6 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_6", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_7 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_7", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_8 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_8", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_9 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_9", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_10 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_10", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_11 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_11", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_12 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_12", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_13 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_13", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_14 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_14", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_15 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_15", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_16 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_16", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_17 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_17", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_18 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_18", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_19 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_19", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_20 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_20", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_21 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_21", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_22 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_22", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_23 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_23", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_24 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_24", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_25 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_25", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_26 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_26", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_27 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_27", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_28 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_28", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_29 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_29", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_30 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_30", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_31 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_31", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_32 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_32", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_33 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_33", false);
	}
	
	BOOST_AUTO_TEST_CASE( SAXPBYPZ_ARG_EO_34 )
	{
	    SaxpbypzEvenOddTester("saxpbypz_eo_input_34", false);
	}
	
BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(GAUSSIAN_EO)

	class GaussianEvenOddTester: public SpinorStaggeredTester{
	   public:
		GaussianEvenOddTester(std::string inputfile) : SpinorStaggeredTester("gaussian_spinorfield_eo", inputfile, 1, 2){
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements, device);
			hardware::buffers::Plain<hmc_float> sqnorm(1, device);

			su3vec * outHost;
			outHost = new su3vec[spinorfieldEvenOddElements * iterations];
			BOOST_REQUIRE(out);
				
			auto prng_buf = prng->get_buffers().at(0);
		
			double sum = 0;
			for (int i = 0; i < iterations; i++) {
				if(i%200==0) logger.info() << "Run kernel for the " << i << "th time";
				code->set_gaussian_spinorfield_eoprec_device(&out, prng_buf);
				out.dump(&outHost[i * spinorfieldEvenOddElements]);
				sum += count_sf(&outHost[i * spinorfieldEvenOddElements], spinorfieldEvenOddElements);
			}
			//sum is the sum of iterations*NUM_ELEMENTS_SF*6 real gaussian numbers
			sum /= (iterations * spinorfieldEvenOddElements * 6);
			kernelResult[0] = sum;

			if(calcVariance){
				double var = 0.;
				for (int i = 0; i < iterations; i++) {
				   var += calc_var_sf(&outHost[i * spinorfieldEvenOddElements], spinorfieldEvenOddElements, sum);
				}
				//var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
				var /= (iterations * spinorfieldEvenOddElements * 6);
				kernelResult[0] = sqrt(var);
			}
			
			/** 
			 * @TODO This piece of code contains actually tests for the RNG itself 
			 *       and should be moved elsewhere. 
			 */
			//The Kolmogorov_Smirnov test requires set of n samples with n around 1000
			//(to big n and to small n are not good choices for this test)
			//So we can consider each kernel result as a set (NUM_ELEMENTS_SF*6=1536 for a 4^4 lattice)
			vector<vector<hmc_float>> samples;
			vector<hmc_float> tmp;
			vector<hmc_float> tmp2;
			for(int i=0; i<iterations; i++){
			  vector<hmc_float> tmp;
			  for(uint j=0; j<spinorfieldEvenOddElements; j++){
			    tmp2=reals_from_su3vec(outHost[i*spinorfieldEvenOddElements+j]);
			    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
			    tmp2.clear();
			  }
			  samples.push_back(tmp);
			  tmp.clear();
			}
			logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
			logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";
			
			if(!calcVariance){
			  //Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
			  mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
			}else{
			  //Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
			  //Note that in the test BOOST_CHECK is used.
			  variance_test_multiple_set(samples,2.,sqrt(0.5));
			  variance_test_multiple_set(samples,3.,sqrt(0.5));
			  variance_test_multiple_set(samples,4.,sqrt(0.5));
			}
		}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_1 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_2 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_3 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_4 )
	{
	    GaussianEvenOddTester("gaussian_eo_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SAX_VEC_AND_SQNORM)

	class SaxVecAndSqnormEvenOddTester: public SpinorStaggeredTester{
	   public:
		SaxVecAndSqnormEvenOddTester(std::string inputfile) : SpinorStaggeredTester("sax_vectorized_and_squarenorm_eoprec",inputfile, 1){
			int NUM_EQS = parameters->get_md_approx_ord();
			const hardware::buffers::SU3vec in(spinorfieldEvenOddElements*NUM_EQS, device);
			const hardware::buffers::SU3vec out(spinorfieldEvenOddElements*NUM_EQS, device);
			in.load(createSpinorfield(spinorfieldEvenOddElements, 123));
			
			hardware::buffers::Plain<hmc_float> sqnorm(NUM_EQS, device);
			hardware::buffers::Plain<hmc_float> alpha(NUM_EQS, device);
			std::vector<hmc_float> alpha_vec_host(NUM_EQS, alpha_host.re);
			for(uint i=0; i<alpha_vec_host.size(); i++)
			    alpha_vec_host[i] += i*alpha_host.im;
			alpha.load(&alpha_vec_host[0]);
			
			code->sax_vectorized_and_squarenorm_eoprec_device(&in, &alpha, NUM_EQS, &sqnorm);
			std::vector<hmc_float> cpu_res(NUM_EQS);
			sqnorm.dump(&cpu_res[0]);

			kernelResult[0] = std::accumulate(cpu_res.begin(), cpu_res.end(), 0.0);

  /*
  print_staggeredfield_eo_to_textfile("ref_vec_sax_and_sq_eo", createSpinorfield(spinorfieldEvenOddElements, 123)); 
  logger.info() << "Produced the ref_vec_sax_and_sq_eo text file with the staggered field for the ref. code.";
  */
		}
	};

	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_1 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_2 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_3 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_4 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_5 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_6 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_7 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_8 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_8");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_9 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_9");
	}
	
	BOOST_AUTO_TEST_CASE( SAX_VEC_AND_SQNORM_10 )
	{
	    SaxVecAndSqnormEvenOddTester("/sax_vec_sqnorm_eo_input_10");
	}

BOOST_AUTO_TEST_SUITE_END()
