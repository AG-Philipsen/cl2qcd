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
#define BOOST_TEST_MODULE OPENCL_MODULE_REAL

#include "real.hpp"
#include "kernelTester.hpp"

class RealTester : public KernelTester {
   public:
	RealTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1,
	    int typeOfComparision = 1) : 
	    KernelTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues, typeOfComparision) {
		code = device->getRealCode();
		hmc_float alpha_host = parameters->get_beta();
		hmc_float beta_host = parameters->get_kappa();
		alpha = new hardware::buffers::Plain<hmc_float>(1, device);
		beta = new hardware::buffers::Plain<hmc_float>(1, device);
		result = new hardware::buffers::Plain<hmc_float>(1, device);
		alpha->load(&alpha_host);
		beta->load(&beta_host);
	}
	
	virtual ~RealTester(){
		delete alpha;
		delete beta;
		delete result;
		code = NULL;
	}
    
   protected:
	const hardware::code::Real * code;
	hardware::buffers::Plain<hmc_float> *alpha;
	hardware::buffers::Plain<hmc_float> *beta;
	hardware::buffers::Plain<hmc_float> *result;

	std::string getSpecificInputfile(std::string inputfileIn){
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "real/" + inputfileIn;
	}
};

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
	    BOOST_CHECK_NO_THROW(RealTester("build", "real_build_input_1", 0));
	}
	
	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
	    BOOST_CHECK_NO_THROW(RealTester("build", "real_build_input_2", 0));
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(PRODUCT)

	class RealProductTester: public RealTester{
	  public:
		RealProductTester(std::string inputfile, bool multiple_operation = false) :
		   RealTester("product", inputfile){
			code->set_real_to_product_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_product_device(alpha, result, result);
			
			result->dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( PRODUCT_1 )
	{
	    RealProductTester("product_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_2 )
	{
	    RealProductTester("product_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_3 )
	{
	    RealProductTester("product_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_4 )
	{
	    RealProductTester("product_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_5 )
	{
	    RealProductTester("product_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_6 )
	{
	    RealProductTester("product_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_7 )
	{
	    RealProductTester("product_input_7", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_8 )
	{
	    RealProductTester("product_input_8", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_9 )
	{
	    RealProductTester("product_input_9", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_10 )
	{
	    RealProductTester("product_input_10", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_11 )
	{
	    RealProductTester("product_input_11", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_12 )
	{
	    RealProductTester("product_input_12", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(RATIO)

	class RealRatioTester: public RealTester{
	  public:
		RealRatioTester(std::string inputfile, bool multiple_operation = false) :
		   RealTester("ratio", inputfile){
			code->set_real_to_ratio_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_ratio_device(result, beta, result);
			
			result->dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
	    RealRatioTester("ratio_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_2 )
	{
	    RealRatioTester("ratio_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_3 )
	{
	    RealRatioTester("ratio_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_4 )
	{
	    RealRatioTester("ratio_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_5 )
	{
	    RealRatioTester("ratio_input_5", true);
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_6 )
	{
	    RealRatioTester("ratio_input_6", true);
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_7 )
	{
	    RealRatioTester("ratio_input_7", true);
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_8 )
	{
	    RealRatioTester("ratio_input_8", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SUM)

	class RealSumTester: public RealTester{
	  public:
		RealSumTester(std::string inputfile, bool multiple_operation = false) :
		   RealTester("sum", inputfile){
			code->set_real_to_sum_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_sum_device(alpha, result, result);
			
			result->dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( SUM_1 )
	{
	    RealSumTester("sum_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SUM_2 )
	{
	    RealSumTester("sum_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SUM_3 )
	{
	    RealSumTester("sum_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SUM_4 )
	{
	    RealSumTester("sum_input_4", true);
	}
	
	BOOST_AUTO_TEST_CASE( SUM_5 )
	{
	    RealSumTester("sum_input_5", true);
	}
	
	BOOST_AUTO_TEST_CASE( SUM_6 )
	{
	    RealSumTester("sum_input_6", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(DIFFERENCE)

	class RealDifferenceTester: public RealTester{
	  public:
		RealDifferenceTester(std::string inputfile, bool multiple_operation = false) :
		   RealTester("difference", inputfile){
			code->set_real_to_difference_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_difference_device(result, beta, result);
			
			result->dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( DIFFERENCE_1 )
	{
	    RealDifferenceTester("difference_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_2 )
	{
	    RealDifferenceTester("difference_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_3 )
	{
	    RealDifferenceTester("difference_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_4 )
	{
	    RealDifferenceTester("difference_input_4", true);
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_5 )
	{
	    RealDifferenceTester("difference_input_5", true);
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_6 )
	{
	    RealDifferenceTester("difference_input_6", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(ACCESS_ELEMENT)

	//This test is quite different from the others and we inherit from KernelTester but we
	//use BOOST directly here. Probably one could think about a better object.
	class RealAccessVectorElementTester: public KernelTester{
	  public:
		RealAccessVectorElementTester(std::string inputfile, bool get_element) :
		   KernelTester("Access_vector_element", ("real/" + inputfile), 0){
		   
			const hardware::code::Real * code = device->getRealCode();
			hardware::buffers::Plain<hmc_float> scalar_buf(1, device);
			hardware::buffers::Plain<hmc_float> vector_buf(4, device);
			std::vector<hmc_float> vector_host(4);
			hmc_float scalar_host;

			if(get_element){
			   vector_host[0] = parameters->get_beta();
			   vector_host[1] = parameters->get_kappa();
			   vector_host[2] = parameters->get_rho();
			   vector_host[3] = parameters->get_mu();
			   vector_buf.load(&vector_host[0]);
			   for(uint i=0; i<vector_host.size(); i++){
			      code->set_real_to_vector_element_device(&vector_buf, i, &scalar_buf);
			      hmc_float cpu_res;
			      scalar_buf.dump(&cpu_res);
			      BOOST_REQUIRE_CLOSE(cpu_res, vector_host[i], 1.e-8);
			   }
			}else{
			   scalar_host = parameters->get_beta() + parameters->get_kappa() +
			                 parameters->get_rho() + parameters->get_mu();
			   scalar_buf.load(&scalar_host);
			   for(uint i=0; i<vector_host.size(); i++){
			      code->set_vector_element_to_real_device(&scalar_buf, i, &vector_buf);
			      std::vector<hmc_float> cpu_res(4);
			      vector_buf.dump(&cpu_res[0]);
			      BOOST_REQUIRE_CLOSE(cpu_res[i], scalar_host, 1.e-8);
			   }
			}
		}
	};

	BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_1 )
	{
	    RealAccessVectorElementTester("access_element_vector_input_1", true);
	}
	
	BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_2 )
	{
	    RealAccessVectorElementTester("access_element_vector_input_1", false);
	} 

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(REAL_UPDATE)

	class RealUpdateTester: public RealTester{
	  public:
		RealUpdateTester(std::string inputfile, int which_update) : RealTester("update", inputfile){
			
		   //Variable 1 and 3 are in the base class as alpha and beta Plain<hmc_float> buffers
		   hardware::buffers::Plain<hmc_float> variable_2(1, device);
		   hardware::buffers::Plain<hmc_float> variable_4(1, device);
		   hardware::buffers::Plain<hmc_float> variable_5(1, device);
		   hardware::buffers::Plain<hmc_float> variable_6(1, device);
		   
		   hmc_float variable_2_host = parameters->get_rho();
		   hmc_float variable_4_host = parameters->get_mu();
		   hmc_float variable_5_host = parameters->get_mass(); 
		   hmc_float variable_6_host = parameters->get_approx_lower(); 
		   
		   variable_2.load(&variable_2_host);
		   variable_4.load(&variable_4_host);
		   variable_5.load(&variable_5_host);
		   variable_6.load(&variable_6_host);
		   
		   if(which_update == 0)
		     code->update_alpha_cgm_device(alpha, &variable_2, beta, &variable_4, &variable_5, 1, result);
		   else if (which_update ==1)
		     code->update_beta_cgm_device(alpha, &variable_2, beta, 1, result);
		   else if (which_update ==2)
		     code->update_zeta_cgm_device(alpha, &variable_2, beta, &variable_4, &variable_5, &variable_6, 1, result);

		   result->dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( ALPHA_1 )
	{
	    RealUpdateTester("update_alpha_input_1", 0);
	}
	
	BOOST_AUTO_TEST_CASE( BETA_1 )
	{
	    RealUpdateTester("update_beta_input_1", 1);
	}
	
	BOOST_AUTO_TEST_CASE( ZETA_1 )
	{
	    RealUpdateTester("update_zeta_input_1", 2);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////
