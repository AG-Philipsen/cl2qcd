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
#include "mockups.hpp"

struct TestParametersBuild : public TestParameters
{
	bool useEo;
	bool useMergeKernelsFermion;
	bool useChemPotRe;
	bool useChemPotIm;

	TestParametersBuild(std::vector<double> referenceValueIn, int nsIn, int ntIn, bool useEoIn, bool useMergeKernelsFermionIn, bool useChemPotReIn, bool useChemPotImIn):
		TestParameters(referenceValueIn, nsIn, ntIn) {
		useEo = useEoIn;
		useMergeKernelsFermion = useMergeKernelsFermionIn;
		useChemPotRe = useChemPotReIn;
		useChemPotIm = useChemPotImIn;
	}
};

struct TestParametersReal : public TestParameters
{
	double beta;
	double kappa;

	TestParametersReal(std::vector<double> referenceValueIn, int nsIn, int ntIn, double betaIn, double kappaIn):
		TestParameters(referenceValueIn, nsIn, ntIn) {
		beta = betaIn;
		kappa = kappaIn;
	}
};

struct TestParametersRealAccess : public TestParametersReal
{
	double rho;
	double mu;

	TestParametersRealAccess(std::vector<double> referenceValueIn, int nsIn, int ntIn ,double betaIn, double kappaIn, double rhoIn, double muIn):
		TestParametersReal(referenceValueIn, nsIn, ntIn, betaIn, kappaIn) {
		rho = rhoIn;
		mu = muIn;
	}
};

struct TestParametersRealUpdate : public TestParametersRealAccess
{
	double mass;

	TestParametersRealUpdate(std::vector<double> referenceValueIn, int nsIn, int ntIn ,double betaIn, double kappaIn, double rhoIn, double muIn, double massIn):
		TestParametersRealAccess(referenceValueIn, nsIn, ntIn, betaIn, kappaIn, rhoIn, muIn) {
		mass = massIn;
	}
};

class RealTester : public KernelTester {
   public:
	RealTester(std::string kernelName, const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup kernelParameters,
			struct TestParametersReal testParams):
//	RealTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1,
//	    int typeOfComparision = 1) :
	    KernelTester(kernelName,  hardwareParameters, kernelParameters, testParams) {
		code = device->getRealCode();
		hmc_float alpha_host = kernelParameters.getBeta();
		hmc_float beta_host = kernelParameters.getKappa();
		alpha = new hardware::buffers::Plain<hmc_float>(1, device);
		beta = new hardware::buffers::Plain<hmc_float>(1, device);
		result = new hardware::buffers::Plain<hmc_float>(1, device);
		alpha->load(&alpha_host);
		beta->load(&beta_host);
	}
	RealTester(std::string kernelName, const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup kernelParameters,
			struct TestParametersBuild testParams):
		KernelTester(kernelName,  hardwareParameters, kernelParameters, testParams) {
		code = device->getRealCode();
		hmc_float alpha_host = kernelParameters.getBeta();
		hmc_float beta_host = kernelParameters.getKappa();
		alpha = new hardware::buffers::Plain<hmc_float>(1, device);
		beta = new hardware::buffers::Plain<hmc_float>(1, device);
		result = new hardware::buffers::Plain<hmc_float>(1, device);
		alpha->load(&alpha_host);
		beta->load(&beta_host);

	};
	
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

	void instantiateMockupsAndCallTester(struct TestParametersBuild testParams)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt);
		BOOST_CHECK_NO_THROW(RealTester("build", hardwareParameters, kernelParameters, testParams));
	}

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		instantiateMockupsAndCallTester(TestParametersBuild {ReferenceValues{0.,0.}, ns4, nt4, true, true, true, true});
//	    BOOST_CHECK_NO_THROW(RealTester("build", "real_build_input_1", 0));
	}

//	BOOST_AUTO_TEST_CASE( BUILD_2 )
//	{
//	    BOOST_CHECK_NO_THROW(RealTester("build", "real_build_input_2", 0));
//	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(PRODUCT)

	class RealProductTester: public RealTester{
	  public:
		RealProductTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
				struct TestParametersReal testParams, bool multiple_operation = false) :
				RealTester("product",hardwareParameters,kernelParameters,testParams){
			code->set_real_to_product_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_product_device(alpha, result, result);

			result->dump(&kernelResult[0]);
		}
	};

	void instantiateMockupsAndCallTester(struct TestParametersReal testParams, bool multiple_operation = false)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt,testParams.beta,testParams.kappa);
		RealProductTester(hardwareParameters, kernelParameters, testParams, multiple_operation);
	}
	BOOST_AUTO_TEST_CASE( PRODUCT_1 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 0.});
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_2 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 1., 0.});
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_3 )//same as 2
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 1., 0.});
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_4 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{1.}, ns4, nt4, 1., 1.});
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_5 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-1.}, ns4, nt4, -1., 1.});
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_6 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-1.}, ns4, nt4, 1., -1.});
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_7 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 0.}, true);
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_8 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 1., 0.}, true);
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_9 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{1.}, ns4, nt4, -1., 1.}, true);
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_10 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{1.}, ns4, nt4, 1., 1.}, true);
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_11 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{63.}, ns4, nt4, -3., 7.}, true);
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_12 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-100.}, ns4, nt4, 10., -1.}, true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(RATIO)

	class RealRatioTester: public RealTester{
	  public:
		RealRatioTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
				struct TestParametersReal testParams, bool multiple_operation = false) :
		   RealTester("ratio",hardwareParameters,kernelParameters,testParams){
			code->set_real_to_ratio_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_ratio_device(result, beta, result);

			result->dump(&kernelResult[0]);
		}
	};

	void instantiateMockupsAndCallTester(struct TestParametersReal testParams, bool multiple_operation = false)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt,testParams.beta,testParams.kappa);
		RealRatioTester(hardwareParameters, kernelParameters, testParams, multiple_operation);
	}

	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 1.});
	}

	BOOST_AUTO_TEST_CASE( RATIO_2 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{1.}, ns4, nt4, 1., 1.});
	}

	BOOST_AUTO_TEST_CASE( RATIO_3 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{1.}, ns4, nt4, 5., 5.});
	}

	BOOST_AUTO_TEST_CASE( RATIO_4 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-1.}, ns4, nt4, 1., -1.});
	}

	BOOST_AUTO_TEST_CASE( RATIO_5 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 1.}, true);
	}

	BOOST_AUTO_TEST_CASE( RATIO_6 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{1.}, ns4, nt4, 1., 1.}, true);
	}

	BOOST_AUTO_TEST_CASE( RATIO_7 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-1.}, ns4, nt4, -4., 2.}, true);
	}

	BOOST_AUTO_TEST_CASE( RATIO_8 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{1.}, ns4, nt4, 1., -1.}, true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SUM)

	class RealSumTester: public RealTester{
	  public:
		RealSumTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
				struct TestParametersReal testParams, bool multiple_operation = false) :
		   RealTester("sum",hardwareParameters,kernelParameters,testParams){
			code->set_real_to_sum_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_sum_device(alpha, result, result);

			result->dump(&kernelResult[0]);
		}
	};

	void instantiateMockupsAndCallTester(struct TestParametersReal testParams, bool multiple_operation = false)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt,testParams.beta,testParams.kappa);
		RealSumTester(hardwareParameters, kernelParameters, testParams, multiple_operation);
	}

	BOOST_AUTO_TEST_CASE( SUM_1 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 0.});
	}

	BOOST_AUTO_TEST_CASE( SUM_2 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 3., -3.});
	}

	BOOST_AUTO_TEST_CASE( SUM_3 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-6.66}, ns4, nt4, 1.23, -7.89});
	}

	BOOST_AUTO_TEST_CASE( SUM_4 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 0.}, true);
	}

	BOOST_AUTO_TEST_CASE( SUM_5 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{3.}, ns4, nt4, 3., -3.}, true);
	}

	BOOST_AUTO_TEST_CASE( SUM_6 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-5.43}, ns4, nt4, 1.23, -7.89}, true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(DIFFERENCE)

	class RealDifferenceTester: public RealTester{
	  public:
		RealDifferenceTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
				struct TestParametersReal testParams, bool multiple_operation = false) :
		   RealTester("difference",hardwareParameters,kernelParameters,testParams){
			code->set_real_to_difference_device(alpha, beta, result);
			if(multiple_operation)
			  code->set_real_to_difference_device(result, beta, result);

			result->dump(&kernelResult[0]);
		}
	};

	void instantiateMockupsAndCallTester(struct TestParametersReal testParams, bool multiple_operation = false)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt,testParams.beta,testParams.kappa);
		RealDifferenceTester(hardwareParameters, kernelParameters, testParams, multiple_operation);
	}

	BOOST_AUTO_TEST_CASE( DIFFERENCE_1 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 0.});
	}

	BOOST_AUTO_TEST_CASE( DIFFERENCE_2 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 1., 1.});
	}

	BOOST_AUTO_TEST_CASE( DIFFERENCE_3 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-2.2}, ns4, nt4, 1.5, 3.7});
	}

	BOOST_AUTO_TEST_CASE( DIFFERENCE_4 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{0.}, ns4, nt4, 0., 0.}, true);
	}

	BOOST_AUTO_TEST_CASE( DIFFERENCE_5 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-1.}, ns4, nt4, 1., 1.}, true);
	}

	BOOST_AUTO_TEST_CASE( DIFFERENCE_6 )
	{
		instantiateMockupsAndCallTester(TestParametersReal {ReferenceValues{-5.9}, ns4, nt4, 1.5, 3.7}, true);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(ACCESS_ELEMENT)

	//This test is quite different from the others and we inherit from KernelTester but we
	//use BOOST directly here. Probably one could think about a better object.
	class RealAccessVectorElementTester: public KernelTester{
	  public:
		RealAccessVectorElementTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
				struct TestParametersRealAccess testParams, bool get_element) :
//		   KernelTester("Access_vector_element", ("real/" + inputfile), 0){
			KernelTester("Access_vector_element", hardwareParameters, kernelParameters, testParams){

			const hardware::code::Real * code = device->getRealCode();
			hardware::buffers::Plain<hmc_float> scalar_buf(1, device);
			hardware::buffers::Plain<hmc_float> vector_buf(4, device);
			std::vector<hmc_float> vector_host(4);
			hmc_float scalar_host;

			if(get_element){
			   vector_host[0] = kernelParameters.getBeta();
			   vector_host[1] = kernelParameters.getKappa();
			   vector_host[2] = kernelParameters.getRho();
			   vector_host[3] = kernelParameters.getMu();
			   vector_buf.load(&vector_host[0]);
			   for(uint i=0; i<vector_host.size(); i++){
			      code->set_real_to_vector_element_device(&vector_buf, i, &scalar_buf);
			      hmc_float cpu_res;
			      scalar_buf.dump(&cpu_res);
			      BOOST_REQUIRE_CLOSE(cpu_res, vector_host[i], 1.e-8);
			   }
			}else{
			   scalar_host = kernelParameters.getBeta() + kernelParameters.getKappa() +
			                 kernelParameters.getRho() + kernelParameters.getMu();
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

	void instantiateMockupsAndCallTester(struct TestParametersRealAccess testParams, bool get_element)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt,testParams.rho,testParams.beta,testParams.kappa,testParams.mu);
		RealAccessVectorElementTester(hardwareParameters, kernelParameters, testParams, get_element);
	}

	BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_1 )
	{
		instantiateMockupsAndCallTester(TestParametersRealAccess {ReferenceValues{}, ns4, nt4, 1.1, 2.2, 3.3, 4.4}, true);
	}

	BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_2 )
	{
		instantiateMockupsAndCallTester(TestParametersRealAccess {ReferenceValues{}, ns4, nt4, 1.1, 2.2, 3.3, 4.4}, false);
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(REAL_UPDATE)

	class RealUpdateTester: public RealTester{
	  public:
		RealUpdateTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
				struct TestParametersRealUpdate testParams, int which_update) :
			RealTester("update", hardwareParameters, kernelParameters, testParams){

		   //Variable 1 and 3 are in the base class as alpha and beta Plain<hmc_float> buffers
		   hardware::buffers::Plain<hmc_float> variable_2(1, device);
		   hardware::buffers::Plain<hmc_float> variable_4(1, device);
		   hardware::buffers::Plain<hmc_float> variable_5(1, device);
		   hardware::buffers::Plain<hmc_float> variable_6(1, device);

		   hmc_float variable_2_host = kernelParameters.getRho();
		   hmc_float variable_4_host = kernelParameters.getMu();
		   hmc_float variable_5_host = kernelParameters.getMass();
		   hmc_float variable_6_host = kernelParameters.getApproxLower();

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

	void instantiateMockupsAndCallTester(struct TestParametersRealUpdate testParams, int which_update)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt,testParams.rho,testParams.beta,testParams.kappa,testParams.mu,testParams.mass);
		RealUpdateTester(hardwareParameters, kernelParameters, testParams, which_update);
	}

	BOOST_AUTO_TEST_CASE( ALPHA_1 )
	{
		instantiateMockupsAndCallTester(TestParametersRealUpdate {ReferenceValues{-0.6375}, ns4, nt4, 1.35, -1.25, 3.4, -2., 4.5}, 0);
	}

	BOOST_AUTO_TEST_CASE( BETA_1 )//needs beta, kappa, rho (hence a different struct)
	{
		instantiateMockupsAndCallTester(TestParametersRealUpdate {ReferenceValues{-6.4752}, ns4, nt4, -1.52, -2.5, 10.65, 0.006, 0.1}, 1); //check these value!! (rho, mass)
//	    RealUpdateTester("update_beta_input_1", 1);
	}

//	BOOST_AUTO_TEST_CASE( ZETA_1 )
//	{
//	    RealUpdateTester("update_zeta_input_1", 2);
//	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////
