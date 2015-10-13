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

#include "testUtilities.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/lattices/gaugefield.hpp"
#include "../device.hpp"
#include "spinors.hpp"
#include "fermions.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

class TestGaugefield {

public:
	TestGaugefield(const hardware::System * system) : system(system), params(&system->get_inputparameters()), prngParameters(&system->get_inputparameters()), prng(*system, &prngParameters), gf(*system, &params, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		const auto & inputfile = system->get_inputparameters();
		meta::print_info_hmc(inputfile);
	};

	const hardware::code::Fermions * get_device();
	const hardware::buffers::SU3 * get_gaugefield();

private:
	const hardware::System * const system;
	const physics::lattices::GaugefieldParametersImplementation params;
	const physics::ParametersPrng_fromMetaInputparameters prngParameters;
	physics::PRNG prng;
	const physics::lattices::Gaugefield gf;
};

void fill_sf_with_one(spinor * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0 = hmc_complex_one;
		sf_in[i].e0.e1 = hmc_complex_one;
		sf_in[i].e0.e2 = hmc_complex_one;
		sf_in[i].e1.e0 = hmc_complex_one;
		sf_in[i].e1.e1 = hmc_complex_one;
		sf_in[i].e1.e2 = hmc_complex_one;
		sf_in[i].e2.e0 = hmc_complex_one;
		sf_in[i].e2.e1 = hmc_complex_one;
		sf_in[i].e2.e2 = hmc_complex_one;
		sf_in[i].e3.e0 = hmc_complex_one;
		sf_in[i].e3.e1 = hmc_complex_one;
		sf_in[i].e3.e2 = hmc_complex_one;
	}
	return;
}

void fill_sf_with_random(spinor * sf_in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0.re = prng_double();
		sf_in[i].e0.e1.re = prng_double();
		sf_in[i].e0.e2.re = prng_double();
		sf_in[i].e1.e0.re = prng_double();
		sf_in[i].e1.e1.re = prng_double();
		sf_in[i].e1.e2.re = prng_double();
		sf_in[i].e2.e0.re = prng_double();
		sf_in[i].e2.e1.re = prng_double();
		sf_in[i].e2.e2.re = prng_double();
		sf_in[i].e3.e0.re = prng_double();
		sf_in[i].e3.e1.re = prng_double();
		sf_in[i].e3.e2.re = prng_double();

		sf_in[i].e0.e0.im = prng_double();
		sf_in[i].e0.e1.im = prng_double();
		sf_in[i].e0.e2.im = prng_double();
		sf_in[i].e1.e0.im = prng_double();
		sf_in[i].e1.e1.im = prng_double();
		sf_in[i].e1.e2.im = prng_double();
		sf_in[i].e2.e0.im = prng_double();
		sf_in[i].e2.e1.im = prng_double();
		sf_in[i].e2.e2.im = prng_double();
		sf_in[i].e3.e0.im = prng_double();
		sf_in[i].e3.e1.im = prng_double();
		sf_in[i].e3.e2.im = prng_double();
	}
	return;
}

void fill_sf_with_random(spinor * sf_in, int size)
{
	fill_sf_with_random(sf_in, size, 123456);
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_hmc";
	logger.info() << "Init device";
	auto params = createParameters("fermionsMerged/" + inputfile);
	hardware::System system(*params);
	TestGaugefield cpu(&system);
	BOOST_MESSAGE("Test done");
}

hmc_float calc_sf_sum(size_t NUM_ELEMS, spinor * in)
{
	hmc_float res = 0.;
	for(size_t i = 0; i < NUM_ELEMS; i++) {
		spinor tmp = in[i];
		res +=
		  tmp.e0.e0.re + tmp.e0.e0.im +
		  tmp.e0.e1.re + tmp.e0.e1.im +
		  tmp.e0.e2.re + tmp.e0.e2.im +
		  tmp.e1.e0.re + tmp.e1.e0.im +
		  tmp.e1.e1.re + tmp.e1.e1.im +
		  tmp.e1.e2.re + tmp.e1.e2.im +
		  tmp.e2.e0.re + tmp.e2.e0.im +
		  tmp.e2.e1.re + tmp.e2.e1.im +
		  tmp.e3.e2.re + tmp.e2.e2.im +
		  tmp.e3.e0.re + tmp.e3.e0.im +
		  tmp.e3.e1.re + tmp.e3.e1.im +
		  tmp.e3.e2.re + tmp.e3.e2.im ;
	}
	return res;
}

const hardware::code::Fermions* TestGaugefield::get_device()
{
	return system->get_devices()[0]->getFermionCode();
}

const hardware::buffers::SU3 * TestGaugefield::get_gaugefield()
{
	return gf.get_buffers().at(0);
}


void test_dslash_and_m_tm_inverse_sitediagonal_plus_minus(std::string inputfile, bool switcher)
{
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher)
		kernelName = "dslash_AND_m_tm_inverse_sitediagonal";
	else
		kernelName = "dslash_AND_m_tm_inverse_sitediagonal_minus";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	auto params = createParameters("fermionsMerged/" + inputfile);
	hardware::System system(*params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	auto * device = cpu.get_device();
	spinor * sf_in;
	spinor * sf_out;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(*params);

	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params->get_solver() == common::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	out.load(sf_in);

	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->getSpinorCode();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);

	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	if(params->get_read_multiple_configs()) {
		if(switcher)
			device->dslash_AND_M_tm_inverse_sitediagonal_eo_device(&in, &out, cpu.get_gaugefield(), EVEN, params->get_kappa(), meta::get_mubar(*params));
		else
			device->dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(&in, &out, cpu.get_gaugefield(), EVEN, params->get_kappa(), meta::get_mubar(*params));
	} else {
		if(switcher)
			device->dslash_AND_M_tm_inverse_sitediagonal_eo_device(&in, &out, cpu.get_gaugefield(), ODD, params->get_kappa(), meta::get_mubar(*params));
		else
			device->dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(&in, &out, cpu.get_gaugefield(), ODD, params->get_kappa(), meta::get_mubar(*params));
	}
	out.dump(sf_out);
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, *params);
	BOOST_MESSAGE("Test done");
}

void test_dslash_and_m_tm_inverse_sitediagonal(std::string inputfile)
{
	test_dslash_and_m_tm_inverse_sitediagonal_plus_minus(inputfile, true);
}

void test_dslash_and_m_tm_inverse_sitediagonal_minus(std::string inputfile)
{
	test_dslash_and_m_tm_inverse_sitediagonal_plus_minus(inputfile, false);
}

void test_m_tm_sitediagonal_plus_minus_and_gamma5_eo(std::string inputfile, bool switcher)
{
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher) kernelName = "m_tm_sitedigaonal_AND_gamma5_eo";
	else kernelName = "m_tm_sitediagonal_minus_AND_gamma5_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	auto params = createParameters("fermionsMerged/" + inputfile);
	hardware::System system(*params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	auto * device = cpu.get_device();
	spinor * sf_in;
	spinor * sf_out;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(*params);

	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params->get_solver() == common::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	out.load(sf_in);
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->getSpinorCode();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);

	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	if(switcher)
		device->M_tm_sitediagonal_AND_gamma5_eo_device(&in, &out, meta::get_mubar(*params));
	else
		device->M_tm_sitediagonal_minus_AND_gamma5_eo_device(&in, &out, meta::get_mubar(*params));
	out.dump(sf_out);
	logger.info() << "result:";
	hmc_float cpu_res;
	cpu_res = calc_sf_sum(NUM_ELEMENTS_SF, sf_out);
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, *params);
	BOOST_MESSAGE("Test done");
}

void test_m_tm_sitediagonal_and_gamma5_eo(std::string inputfile)
{
	test_m_tm_sitediagonal_plus_minus_and_gamma5_eo(inputfile, true);
}

void test_m_tm_sitediagonal_minus_and_gamma5_eo(std::string inputfile)
{
	test_m_tm_sitediagonal_plus_minus_and_gamma5_eo(inputfile, false);
}

//CP: Note: this is the same test as in the "normal" opencl_module_fermions test, I left it here, too.
BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO )

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_1)
{
	test_dslash_and_m_tm_inverse_sitediagonal("/dslash_and_m_tm_inverse_sitediagonal_input_1");
}

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_2)
{
	test_dslash_and_m_tm_inverse_sitediagonal("/dslash_and_m_tm_inverse_sitediagonal_input_2");
}

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_3)
{
	test_dslash_and_m_tm_inverse_sitediagonal("/dslash_and_m_tm_inverse_sitediagonal_input_3");
}

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_4)
{
	test_dslash_and_m_tm_inverse_sitediagonal("/dslash_and_m_tm_inverse_sitediagonal_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO )

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_1)
{
	test_dslash_and_m_tm_inverse_sitediagonal_minus("/dslash_and_m_tm_inverse_sitediagonal_minus_input_1");
}

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_2)
{
	test_dslash_and_m_tm_inverse_sitediagonal_minus("/dslash_and_m_tm_inverse_sitediagonal_minus_input_2");
}

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_3)
{
	test_dslash_and_m_tm_inverse_sitediagonal_minus("/dslash_and_m_tm_inverse_sitediagonal_minus_input_3");
}

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_4)
{
	test_dslash_and_m_tm_inverse_sitediagonal_minus("/dslash_and_m_tm_inverse_sitediagonal_minus_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO )

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_1)
{
	test_m_tm_sitediagonal_and_gamma5_eo("/m_tm_sitediagonal_and_gamma5_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_2)
{
	test_m_tm_sitediagonal_and_gamma5_eo("/m_tm_sitediagonal_and_gamma5_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_3)
{
	test_m_tm_sitediagonal_and_gamma5_eo("/m_tm_sitediagonal_and_gamma5_eo_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO )

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_1)
{
	test_m_tm_sitediagonal_minus_and_gamma5_eo("/m_tm_sitediagonal_minus_and_gamma5_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_2)
{
	test_m_tm_sitediagonal_minus_and_gamma5_eo("/m_tm_sitediagonal_minus_and_gamma5_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_3)
{
	test_m_tm_sitediagonal_minus_and_gamma5_eo("/m_tm_sitediagonal_minus_and_gamma5_eo_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

#include "FermionTester.hpp"

BOOST_AUTO_TEST_SUITE(SAXPY_AND_GAMMA5_EO )

	class SaxpyAndGamma5EvenOddTester : public FermionTester
	{
	public:
		SaxpyAndGamma5EvenOddTester(std::vector<std::string> parameterStrings, uint num) :
			FermionTester("saxpy_AND_gamma5", parameterStrings, 1)
		{
			const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
			const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
			const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);
			hardware::buffers::Plain<hmc_complex> alpha(1, device);

			if (num == 1){
				in.load(createSpinorfield(spinorfieldEvenOddElements));
				in2.load(createSpinorfield(spinorfieldEvenOddElements));
			}
			else{
				in.load(createSpinorfieldWithOnesAndMinusOneForGamma5Use(spinorfieldEvenOddElements));
				in2.load(createSpinorfieldWithOnesAndMinusOneForGamma5Use(spinorfieldEvenOddElements));
			}
			alpha.load(&alpha_host);
			
			spinor * sf_in;
			sf_in = new spinor[spinorfieldEvenOddElements];
			
			code->saxpy_AND_gamma5_eo_device(&in, &in2, alpha_host, &out);
			out.dump(sf_in);
			kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);
	
			delete sf_in;
		}
	};

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_1)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=0", "--rho=0", "--test_ref_val=0", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,1);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_2)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=0", "--rho=1", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,1);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_3)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=0", "--rho=-1", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,1);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_4)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=1", "--rho=0", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,1);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_5)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=-1", "--rho=0", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,1);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_6)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=1", "--rho=-1", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,1);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_7)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=-1", "--rho=-1", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,1);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_8)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=0", "--rho=0", "--test_ref_val=1536", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,2);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_9)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=0", "--rho=1", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,2);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_10)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=0", "--rho=-1", "--test_ref_val=3072.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,2);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_11)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=1", "--rho=0", "--test_ref_val=0.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,2);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_12)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=-1", "--rho=0", "--test_ref_val=3072.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,2);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_13)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=1", "--rho=-1", "--test_ref_val=1536.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,2);
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_14)
	{
		std::vector<std::string> parameterStrings {"--nspace=4", "--ntime=4", "--solver=cg", "--use_merge_kernels_fermion=true" , "--beta=-1", "--rho=-1", "--test_ref_val=4608.", "--test_ref_val2=0"};
		SaxpyAndGamma5EvenOddTester tester(parameterStrings,2);
	}
	
BOOST_AUTO_TEST_SUITE_END()
