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
#define BOOST_TEST_MODULE OPENCL_MODULE_COMPLEX

#include "kernelTester.hpp"
#include "complex.hpp"

class ComplexTester : public KernelTester {
   public:
	ComplexTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 2,
	    int typeOfComparision = 1) : 
	    KernelTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues, typeOfComparision) {
		code = device->getComplexCode();
		hmc_complex alpha_host = {parameters->get_beta(), parameters->get_rho()};
		hmc_complex beta_host = {parameters->get_kappa(), parameters->get_mu()};
		alpha = new hardware::buffers::Plain<hmc_complex>(1, device);
		beta = new hardware::buffers::Plain<hmc_complex>(1, device);
		result = new hardware::buffers::Plain<hmc_complex>(1, device);
		alpha->load(&alpha_host);
		beta->load(&beta_host);
	}
	
	void storeResultAsComplex(){
		hmc_complex tmp;
		result->dump(&tmp);
		kernelResult[0] = tmp.re;
		kernelResult[1] = tmp.im;
	}
	
	virtual ~ComplexTester(){
		delete alpha;
		delete beta;
		delete result;
		code = NULL;
	}
    
   protected:
	const hardware::code::Complex * code;
	hardware::buffers::Plain<hmc_complex> *alpha;
	hardware::buffers::Plain<hmc_complex> *beta;
	hardware::buffers::Plain<hmc_complex> *result;

	std::string getSpecificInputfile(std::string inputfileIn){
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "complex/" + inputfileIn;
	}
};

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
	    BOOST_CHECK_NO_THROW(ComplexTester("build", "complex_build_input_1", 0));
	}
	
	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
	    BOOST_CHECK_NO_THROW(ComplexTester("build", "complex_build_input_2", 0));
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(PRODUCT)

	class ComplexProductTester: public ComplexTester{
	  public:
		ComplexProductTester(std::string inputfile, bool multiple_operation = false) :
		   ComplexTester("product", inputfile){
			code->set_complex_to_product_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_complex_to_product_device(alpha, result, result);
			
			storeResultAsComplex();
		}
	};

	BOOST_AUTO_TEST_CASE( PRODUCT_1 )
	{
	    ComplexProductTester("product_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_2 )
	{
	    ComplexProductTester("product_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_3 )
	{
	    ComplexProductTester("product_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_4 )
	{
	    ComplexProductTester("product_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_5 )
	{
	    ComplexProductTester("product_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_6 )
	{
	    ComplexProductTester("product_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_7 )
	{
	    ComplexProductTester("product_input_7", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_8 )
	{
	    ComplexProductTester("product_input_8", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_9 )
	{
	    ComplexProductTester("product_input_9", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_10 )
	{
	    ComplexProductTester("product_input_10", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_11 )
	{
	    ComplexProductTester("product_input_11", true);
	}
	
	BOOST_AUTO_TEST_CASE( PRODUCT_12 )
	{
	    ComplexProductTester("product_input_12", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(RATIO)

	class ComplexRatioTester: public ComplexTester{
	  public:
		ComplexRatioTester(std::string inputfile, bool multiple_operation = false) :
		   ComplexTester("ratio", inputfile){
			code->set_complex_to_ratio_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_complex_to_ratio_device(result, beta, result);
			
			storeResultAsComplex();
		}
	};

	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
	    ComplexRatioTester("ratio_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_2 )
	{
	    ComplexRatioTester("ratio_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_3 )
	{
	    ComplexRatioTester("ratio_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_4 )
	{
	    ComplexRatioTester("ratio_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_5 )
	{
	    ComplexRatioTester("ratio_input_5", true);
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_6 )
	{
	    ComplexRatioTester("ratio_input_6", true);
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_7 )
	{
	    ComplexRatioTester("ratio_input_7", true);
	}
	
	BOOST_AUTO_TEST_CASE( RATIO_8 )
	{
	    ComplexRatioTester("ratio_input_8", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SUM)

	class ComplexSumTester: public ComplexTester{
	  public:
		ComplexSumTester(std::string inputfile, bool multiple_operation = false) :
		   ComplexTester("sum", inputfile){
			code->set_complex_to_sum_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_complex_to_sum_device(alpha, result, result);
			
			storeResultAsComplex();
		}
	};

	BOOST_AUTO_TEST_CASE( SUM_1 )
	{
	    ComplexSumTester("sum_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( SUM_2 )
	{
	    ComplexSumTester("sum_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( SUM_3 )
	{
	    ComplexSumTester("sum_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( SUM_4 )
	{
	    ComplexSumTester("sum_input_4", true);
	}
	
	BOOST_AUTO_TEST_CASE( SUM_5 )
	{
	    ComplexSumTester("sum_input_5", true);
	}
	
	BOOST_AUTO_TEST_CASE( SUM_6 )
	{
	    ComplexSumTester("sum_input_6", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(DIFFERENCE)

	class ComplexDifferenceTester: public ComplexTester{
	  public:
		ComplexDifferenceTester(std::string inputfile, bool multiple_operation = false) :
		   ComplexTester("difference", inputfile){
			code->set_complex_to_difference_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_complex_to_difference_device(result, beta, result);
			
			storeResultAsComplex();
		}
	};

	BOOST_AUTO_TEST_CASE( DIFFERENCE_1 )
	{
	    ComplexDifferenceTester("difference_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_2 )
	{
	    ComplexDifferenceTester("difference_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_3 )
	{
	    ComplexDifferenceTester("difference_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_4 )
	{
	    ComplexDifferenceTester("difference_input_4", true);
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_5 )
	{
	    ComplexDifferenceTester("difference_input_5", true);
	}
	
	BOOST_AUTO_TEST_CASE( DIFFERENCE_6 )
	{
	    ComplexDifferenceTester("difference_input_6", true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(CONVERT)

	class ComplexConvertTester: public ComplexTester{
	  public:
		ComplexConvertTester(std::string inputfile) : ComplexTester("convert", inputfile){
			hardware::buffers::Plain<hmc_float> gamma(1, device);
			hmc_float tmp = (parameters->	get_beta());
			gamma.load(&tmp);
			code->set_complex_to_float_device(&gamma, result);
			
			storeResultAsComplex();
		}
	};

	BOOST_AUTO_TEST_CASE( CONVERT_1 )
	{
	    ComplexConvertTester("convert_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( CONVERT_2 )
	{
	    ComplexConvertTester("convert_input_2");
	}

BOOST_AUTO_TEST_SUITE_END()
