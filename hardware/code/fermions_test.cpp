/*
 * Copyright 2012, 2013, 2014 Christopher Pinke, Matthias Bach
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_FERMIONS

#include "SpinorTester.hpp"
#include "gaugefield.hpp"
#include "fermions.hpp"
#include "../../physics/lattices/gaugefield.hpp"

class FermionTester : public SpinorTester
{
public:
	FermionTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
	SpinorTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues)
	{
			code = device->get_fermion_code();
			gaugefield = new physics::lattices::Gaugefield(*system, *prng);
	}
	~FermionTester()
	{
		delete gaugefield;
	}
	
protected:
	const hardware::code::Fermions * code;
	physics::lattices::Gaugefield * gaugefield;
	
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "../fermions/" + inputfileIn;
	}
	
	const hardware::buffers::SU3* getGaugefieldBuffer() {
		return gaugefield->get_buffers()[0];
	}
};

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		FermionTester tester("build", "build_input_1");
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
		FermionTester tester("build", "build_input_2");
	}

BOOST_AUTO_TEST_SUITE_END()

class FermionmatrixTester : public FermionTester
{
public:
	FermionmatrixTester(std::string kernelName, std::string inputfile) :
	FermionTester(kernelName, inputfile, 1)
	{
		in = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		out = new const hardware::buffers::Plain<spinor>(spinorfieldElements, device);
		in->load(createSpinorfield(spinorfieldElements));
		out->load(createSpinorfield(spinorfieldElements));
	}
	~FermionmatrixTester()
	{
		calcSquarenormAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Plain<spinor> * in;
	const hardware::buffers::Plain<spinor> * out;
};

BOOST_AUTO_TEST_SUITE( M_WILSON )

	class MWilsonTester : public FermionmatrixTester
	{
public:
		MWilsonTester(std::string inputfile) :
		FermionmatrixTester("m_wilson", inputfile)
		{
			code->M_wilson_device(in, out,  this->getGaugefieldBuffer(), parameters->get_kappa());
		}
	};

	BOOST_AUTO_TEST_CASE( M_WILSON_1)
	{
		MWilsonTester tester("m_wilson_input_1");
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_2)
	{
		MWilsonTester tester("m_wilson_input_2");
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_3)
	{
		MWilsonTester tester("m_wilson_input_3");
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_4)
	{
		MWilsonTester tester("m_wilson_input_4");
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_5)
	{
		MWilsonTester tester("m_wilson_input_5");
	}

	BOOST_AUTO_TEST_CASE( M_WILSON_6)
	{
		MWilsonTester tester("m_wilson_input_6");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_MINUS  )

	class MTmMinusTester : public FermionmatrixTester
	{
public:
		MTmMinusTester(std::string inputfile) :
		FermionmatrixTester("m_tm_minus", inputfile)
		{
			code->M_tm_minus_device(in, out,  this->getGaugefieldBuffer(), parameters->get_kappa(), meta::get_mubar(*parameters));
		}
	};

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_1 )
	{
		MTmMinusTester tester("m_tm_minus_input_1");
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_2 )
	{
		MTmMinusTester tester("m_tm_minus_input_2");
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_3 )
	{
		MTmMinusTester tester("m_tm_minus_input_3");
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_4 )
	{
		MTmMinusTester tester("m_tm_minus_input_4");
	}

	BOOST_AUTO_TEST_CASE( M_TM_MINUS_5 )
	{
		MTmMinusTester tester("m_tm_minus_input_5");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_PLUS )

	class MTmPlusTester : public FermionmatrixTester
	{
public:
		MTmPlusTester(std::string inputfile) :
		FermionmatrixTester("m_tm_plus", inputfile)
		{
			code->M_tm_plus_device(in, out,  this->getGaugefieldBuffer(), parameters->get_kappa(), meta::get_mubar(*parameters));
		}
	};

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_1 )
	{
		MTmPlusTester tester("m_tm_plus_input_1");
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_2 )
	{
		MTmPlusTester tester("m_tm_plus_input_2");
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_3 )
	{
		MTmPlusTester tester("m_tm_plus_input_3");
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_4 )
	{
		MTmPlusTester tester("m_tm_plus_input_4");
	}

	BOOST_AUTO_TEST_CASE( M_TM_PLUS_5 )
	{
		MTmPlusTester tester("m_tm_plus_input_5");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5 )

	class Gamma5Tester : public FermionTester
	{
public:
		Gamma5Tester(std::string inputfile) :
			FermionTester("gamma5", inputfile, 1)
		{
			const hardware::buffers::Plain<spinor> in(spinorfieldElements, device);
			spinor * sf_in;
			sf_in = new spinor[spinorfieldElements];
			
			in.load( createSpinorfield(spinorfieldElements) );
			code->gamma5_device(&in);
			in.dump(sf_in);
			kernelResult[0] = count_sf(sf_in, spinorfieldElements);
	
			delete sf_in;
		}
	};

	BOOST_AUTO_TEST_CASE( GAMMA5_1)
	{
		Gamma5Tester tester("gamma5_input_1");
	}

	BOOST_AUTO_TEST_CASE( GAMMA5_2 )
	{
		Gamma5Tester tester("gamma5_input_2");
	}

BOOST_AUTO_TEST_SUITE_END()

	class Gamma5EvenOddTester : public FermionTester
	{
public:
		Gamma5EvenOddTester(std::string inputfile) :
			FermionTester("gamma5_eo", inputfile, 1)
		{
			const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
			spinor * sf_in;
			sf_in = new spinor[spinorfieldEvenOddElements];
			
			in.load( createSpinorfield(spinorfieldEvenOddElements) );
			code->gamma5_eo_device(&in);
			in.dump(sf_in);
			kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);
	
			delete sf_in;
		}
	};

	BOOST_AUTO_TEST_SUITE( GAMMA5_EO)

	BOOST_AUTO_TEST_CASE( GAMMA5_EO_1)
	{
		Gamma5EvenOddTester tester("gamma5_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( GAMMA5_EO_2 )
	{
		Gamma5EvenOddTester tester("gamma5_eo_input_2");
	}

BOOST_AUTO_TEST_SUITE_END()

class FermionmatrixEvenOddTester : public FermionTester
{
public:
	FermionmatrixEvenOddTester(std::string kernelName, std::string inputfile) :
	FermionTester(kernelName, inputfile, 1)
	{
		in = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		out = new const hardware::buffers::Spinor(spinorfieldEvenOddElements, device);
		in->load(createSpinorfield(spinorfieldEvenOddElements));
		out->load(createSpinorfield(spinorfieldEvenOddElements));
	}
	~FermionmatrixEvenOddTester()
	{
		calcSquarenormEvenOddAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Spinor * in;
	const hardware::buffers::Spinor * out;
};

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL )

	class MTmSitediagonalTester: public FermionmatrixEvenOddTester
	{
	public:
		MTmSitediagonalTester(std::string inputfile):
			FermionmatrixEvenOddTester("m_tm_sitediagonal", inputfile)
			{
				code->M_tm_sitediagonal_device( in, out);
			}
	};

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_1)
	{
		MTmSitediagonalTester tester("m_tm_sitediagonal_input_1");
	}

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_2)
	{
		MTmSitediagonalTester tester("m_tm_sitediagonal_input_2");
	}

	BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_3)
	{
		MTmSitediagonalTester tester("m_tm_sitediagonal_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()


#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"

#include "fermions.hpp"

//some functionality
#include "../../tests/test_util.h"

class TestGaugefield {

public:
	TestGaugefield(const hardware::System * system) : system(system), prng(*system), gf(*system, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		auto inputfile = system->get_inputparameters();
		meta::print_info_hmc(inputfile);
	};

	const hardware::code::Fermions * get_device();
	const hardware::buffers::SU3 * get_gaugefield();

private:
	const hardware::System * const system;
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

const hardware::code::Fermions* TestGaugefield::get_device()
{
	return system->get_devices()[0]->get_fermion_code();
}

const hardware::buffers::SU3 * TestGaugefield::get_gaugefield()
{
	return gf.get_buffers().at(0);
}

// void test_m_tm_sitediagonal(std::string inputfile)
// {
// 	test_m_tm_sitediagonal_plus_minus(inputfile, true);
// }
// 
// void test_m_tm_sitediagonal_minus(std::string inputfile)
// {
// 	test_m_tm_sitediagonal_plus_minus(inputfile, false);
// }
// 
// void test_m_tm_inverse_sitediagonal_plus_minus(std::string inputfile, bool switcher)
// {
// 	using namespace hardware::buffers;
// 
// 	std::string kernelName;
// 	if(switcher)
// 		kernelName = "m_tm_inverse_sitediagonal";
// 	else
// 		kernelName = "m_tm_inverse_sitediagonal_minus";
// 	printKernelInfo(kernelName);
// 	logger.info() << "Init device";
// 	meta::Inputparameters params = create_parameters(inputfile);
// 	hardware::System system(params);
// 	TestGaugefield cpu(&system);
// 	auto * device = cpu.get_device();
// 	spinor * sf_in;
// 	spinor * sf_out;
// 
// 	logger.info() << "Fill buffers...";
// 	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
// 
// 	sf_in = new spinor[NUM_ELEMENTS_SF];
// 	sf_out = new spinor[NUM_ELEMENTS_SF];
// 
// 	//use the variable use_cg to switch between cold and random input sf
// 	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
// 	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
// 	BOOST_REQUIRE(sf_in);
// 
// 	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
// 	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
// 	in.load(sf_in);
// 	out.load(sf_in);
// 	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
// 
// 	auto spinor_code = device->get_device()->get_spinor_code();
// 
// 	logger.info() << "|phi|^2:";
// 	hmc_float cpu_back;
// 	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
// 	sqnorm.dump(&cpu_back);
// 	logger.info() << cpu_back;
// 
// 	hmc_float cpu_res;
// 	if(switcher) {
// 		device->M_tm_inverse_sitediagonal_device( &in, &out);
// 	} else {
// 		device->M_tm_inverse_sitediagonal_minus_device( &in, &out);
// 	}
// 	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
// 	sqnorm.dump(&cpu_res);
// 	logger.info() << "result:";
// 	logger.info() << cpu_res;
// 
// 	logger.info() << "Clear buffers";
// 	delete[] sf_in;
// 	delete[] sf_out;
// 
// 	testFloatAgainstInputparameters(cpu_res, params);
// 	BOOST_MESSAGE("Test done");
// }
// 
// void test_m_tm_inverse_sitediagonal(std::string inputfile)
// {
// 	test_m_tm_inverse_sitediagonal_plus_minus(inputfile, true);
// }
// 
// void test_m_tm_inverse_sitediagonal_minus(std::string inputfile)
// {
// 	test_m_tm_inverse_sitediagonal_plus_minus(inputfile, false);
// }

void test_dslash_eo(std::string inputfile)
{
	using namespace hardware::buffers;
	std::string kernelName = "dslash_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	size_t NUM_ELEMENTS_SF_EO = hardware::code::get_eoprec_spinorfieldsize(params);
	spinor * sf_in_eo;
	sf_in_eo = new spinor[NUM_ELEMENTS_SF_EO];
	const Spinor in_eo_even(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor out_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in_eo, NUM_ELEMENTS_SF_EO);
	else fill_sf_with_random(sf_in_eo, NUM_ELEMENTS_SF_EO);
	in_eo_even.load(sf_in_eo);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in_eo_even, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	hmc_float cpu_res;
	if(params.get_read_multiple_configs()) {
		device->dslash_eo_device( &in_eo_even, &out_eo, cpu.get_gaugefield(), EVEN, params.get_kappa() );
	} else {
		device->dslash_eo_device( &in_eo_even, &out_eo, cpu.get_gaugefield(), ODD, params.get_kappa() );
	}
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out_eo, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << "result:";
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in_eo;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}


/*
BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL )

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_1)
{
	test_m_tm_inverse_sitediagonal("/m_tm_inverse_sitediagonal_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_2)
{
	test_m_tm_inverse_sitediagonal("/m_tm_inverse_sitediagonal_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_3)
{
	test_m_tm_inverse_sitediagonal("/m_tm_inverse_sitediagonal_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS )

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_1)
{
	test_m_tm_sitediagonal_minus("/m_tm_sitediagonal_minus_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_2)
{
	test_m_tm_sitediagonal_minus("/m_tm_sitediagonal_minus_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_3)
{
	test_m_tm_sitediagonal_minus("/m_tm_sitediagonal_minus_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL_MINUS)

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_1)
{
	test_m_tm_inverse_sitediagonal_minus("/m_tm_inverse_sitediagonal_minus_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_2)
{
	test_m_tm_inverse_sitediagonal_minus("/m_tm_inverse_sitediagonal_minus_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_3)
{
	test_m_tm_inverse_sitediagonal_minus("/m_tm_inverse_sitediagonal_minus_input_3");
}

BOOST_AUTO_TEST_SUITE_END()*/

BOOST_AUTO_TEST_SUITE(DSLASH_EO )

BOOST_AUTO_TEST_CASE( DSLASH_EO_1)
{
	test_dslash_eo("/dslash_eo_input_1");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_2)
{
	test_dslash_eo("/dslash_eo_input_2");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_3)
{
	test_dslash_eo("/dslash_eo_input_3");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_4)
{
	test_dslash_eo("/dslash_eo_input_4");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_5)
{
	test_dslash_eo("/dslash_eo_input_5");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_6)
{
	test_dslash_eo("/dslash_eo_input_6");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_7)
{
	test_dslash_eo("/dslash_eo_input_7");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_8)
{
	test_dslash_eo("/dslash_eo_input_8");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_9)
{
	test_dslash_eo("/dslash_eo_input_9");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_10)
{
	test_dslash_eo("/dslash_eo_input_10");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_11)
{
	test_dslash_eo("/dslash_eo_input_11");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_12)
{
	test_dslash_eo("/dslash_eo_input_12");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_13)
{
	test_dslash_eo("/dslash_eo_input_13");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_14)
{
	test_dslash_eo("/dslash_eo_input_14");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_15)
{
	test_dslash_eo("/dslash_eo_input_15");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_16)
{
	test_dslash_eo("/dslash_eo_input_16");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_17)
{
	test_dslash_eo("/dslash_eo_input_17");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_18)
{
	test_dslash_eo("/dslash_eo_input_18");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_19)
{
	test_dslash_eo("/dslash_eo_input_19");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_20)
{
	test_dslash_eo("/dslash_eo_input_20");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_21)
{
	test_dslash_eo("/dslash_eo_input_21");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_22)
{
	test_dslash_eo("/dslash_eo_input_22");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_23)
{
	test_dslash_eo("/dslash_eo_input_23");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_24)
{
	test_dslash_eo("/dslash_eo_input_24");
}

BOOST_AUTO_TEST_SUITE_END()

void test_m_fermion_compare_noneo_eo(std::string inputfile, int switcher)
{
	//switcher switches between similar functions
	//0: m_wilson (pure wilson)
	//1: m_tm_plus (twisted mass, upper flavour)
	//2: m_tm_minus (twisted mass, lower flavour)
	using namespace hardware::buffers;

	std::string kernelName = "Test equivalence of ";
	if(switcher == 0) {
		kernelName += "m_wilson";
	} else if(switcher == 1) {
		kernelName += "m_tm_plus";
	} else if(switcher == 2) {
		kernelName += "m_tm_minus";
	} else {
		logger.fatal() << "wrong parameter in test_m_fermion";
	}
	kernelName += " in eo- and non-eo formulation";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	auto * device = cpu.get_device();
	spinor * sf_in_noneo;
	spinor * sf_out_noneo;
	spinor * sf_in_eo1;
	spinor * sf_in_eo2;
	spinor * sf_out_eo;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF_EO = hardware::code::get_eoprec_spinorfieldsize(params);

	sf_in_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_out_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_in_eo1 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in_eo2 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out_eo = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in_eo1, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in_eo2, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in_noneo, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random(sf_in_eo1, NUM_ELEMENTS_SF_EO, 123456);
		fill_sf_with_random(sf_in_eo2, NUM_ELEMENTS_SF_EO, 78910);
		fill_sf_with_random(sf_in_noneo, NUM_ELEMENTS_SF, 123456);
	}
	BOOST_REQUIRE(sf_in_eo1);
	BOOST_REQUIRE(sf_in_eo2);
	BOOST_REQUIRE(sf_in_noneo);
	fill_sf_with_one(sf_out_noneo, NUM_ELEMENTS_SF);
	fill_sf_with_one(sf_out_eo, NUM_ELEMENTS_SF);

	const Spinor in_eo1(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in_eo2(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<spinor> out_eo(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in_noneo(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out_noneo(NUM_ELEMENTS_SF, device->get_device());
	//create 3 buffers for intermediate results
	const Spinor out_tmp_eo1(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor out_tmp_eo2(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor tmp_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	in_eo1.load(sf_in_eo1);
	in_eo2.load(sf_in_eo2);
	out_eo.load(sf_out_eo);
	in_noneo.load(sf_in_noneo);
	out_noneo.load(sf_out_noneo);

	auto spinor_code = device->get_device()->get_spinor_code();

	if(params.get_solver() != meta::Inputparameters::cg) {
		//use read_multiple_configs to choose whether to copy the eo rnd vectors to noneo or vise versa
	  if(params.get_read_multiple_configs())
			spinor_code->convert_from_eoprec_device(&in_eo1, &in_eo2, &in_noneo);
		else
			spinor_code->convert_to_eoprec_device(&in_eo1, &in_eo2, &in_noneo);
	}

	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	//create -1 on device
	hmc_complex minusone_tmp = { -1., 0.};
	hardware::buffers::Plain<hmc_complex> minusone(1, device->get_device());
	minusone.load(&minusone_tmp);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	logger.info() << "|phi_noneo|^2:";
	hmc_float cpu_back_noneo;
	spinor_code->set_float_to_global_squarenorm_device(&in_noneo, &sqnorm);
	sqnorm.dump(&cpu_back_noneo);
	logger.info() << cpu_back_noneo;
	logger.info() << "Run kernel";
	if(switcher == 0) {
		device->M_wilson_device(&in_noneo, &out_noneo,  cpu.get_gaugefield(), params.get_kappa());
	} else if(switcher == 1) {
		device->M_tm_plus_device(&in_noneo, &out_noneo,  cpu.get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
	} else if(switcher == 2) {
		device->M_tm_minus_device(&in_noneo, &out_noneo,  cpu.get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
	} else {
		logger.fatal() << "wrong parameter in test_m_fermion";
	}
	logger.info() << "result:";
	hmc_float cpu_res_noneo;
	spinor_code->set_float_to_global_squarenorm_device(&out_noneo, &sqnorm);
	sqnorm.dump(&cpu_res_noneo);
	logger.info() << cpu_res_noneo;

	logger.info() << "|phi_eo1|^2:";
	hmc_float cpu_back_eo1;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in_eo1, &sqnorm);
	sqnorm.dump(&cpu_back_eo1);
	logger.info() << cpu_back_eo1;
	logger.info() << "|phi_eo2|^2:";
	hmc_float cpu_back_eo2;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in_eo2, &sqnorm);
	sqnorm.dump(&cpu_back_eo2);
	logger.info() << cpu_back_eo2;
	logger.info() << "Run kernel";
	if(switcher == 0) {
		//suppose in1 is the even, in2 the odd input vector
		//now calc out_tmp_eo1 = (1 in1 + D_eo in2)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo1);

		device->dslash_eo_device(&in_eo2, &out_tmp_eo1, cpu.get_gaugefield(), EO, params.get_kappa());

		spinor_code->saxpy_eoprec_device(&out_tmp_eo1, &in_eo1, &minusone, &out_tmp_eo1);

		//now calc out_tmp_eo2 = ( 1 in2 + D_oe in1)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo2);

		device->dslash_eo_device(&in_eo1, &out_tmp_eo2, cpu.get_gaugefield(), OE, params.get_kappa());

		spinor_code->saxpy_eoprec_device(&out_tmp_eo2, &in_eo2, &minusone, &out_tmp_eo2);

		//now, both output vectors have to be converted back to noneo
		spinor_code->convert_from_eoprec_device(&out_tmp_eo1, &out_tmp_eo2, &out_eo);
	} else if(switcher == 1) {
		//suppose in1 is the even, in2 the odd input vector
		//now calc out_tmp_eo1 = (R_even in1 + D_eo in2)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo1);

		device->dslash_eo_device(&in_eo2, &out_tmp_eo1, cpu.get_gaugefield(), EO, params.get_kappa());
		device->M_tm_sitediagonal_device(&in_eo1, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo1, &tmp_eo, &minusone, &out_tmp_eo1);

		//now calc out_tmp_eo2 = ( R_odd in2 + D_oe in1)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo2);

		device->dslash_eo_device(&in_eo1, &out_tmp_eo2, cpu.get_gaugefield(), OE, params.get_kappa());
		device->M_tm_sitediagonal_device(&in_eo2, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo2, &tmp_eo, &minusone, &out_tmp_eo2);

		//now, both output vectors have to be converted back to noneo
		spinor_code->convert_from_eoprec_device(&out_tmp_eo1, &out_tmp_eo2, &out_eo);
	} else if(switcher == 2) {
		//suppose in1 is the even, in2 the odd input vector
		//now calc out_tmp_eo1 = (R_even in1 + D_eo in2)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo1);

		device->dslash_eo_device(&in_eo2, &out_tmp_eo1, cpu.get_gaugefield(), EO, params.get_kappa());
		device->M_tm_sitediagonal_minus_device(&in_eo1, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo1, &tmp_eo, &minusone, &out_tmp_eo1);

		//now calc out_tmp_eo2 = ( R_odd in2 + D_oe in1)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo2);

		device->dslash_eo_device(&in_eo1, &out_tmp_eo2, cpu.get_gaugefield(), OE, params.get_kappa());
		device->M_tm_sitediagonal_minus_device(&in_eo2, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo2, &tmp_eo, &minusone, &out_tmp_eo2);

		//now, both output vectors have to be converted back to noneo
		spinor_code->convert_from_eoprec_device(&out_tmp_eo1, &out_tmp_eo2, &out_eo);
	}
	logger.info() << "result:";
	hmc_float cpu_res_eo;
	spinor_code->set_float_to_global_squarenorm_device(&out_eo, &sqnorm);
	sqnorm.dump(&cpu_res_eo);
	logger.info() << cpu_res_eo;

	logger.info() << "Clear buffers";
	delete[] sf_in_noneo;
	delete[] sf_out_noneo;
	delete[] sf_in_eo1;
	delete[] sf_in_eo2;
	delete[] sf_out_eo;

	logger.info() << "Compare eo and non-eo results";
	BOOST_REQUIRE_CLOSE(cpu_res_eo, cpu_res_noneo, params.get_solver_prec() );
	testFloatAgainstInputparameters(cpu_res_noneo, params);
	testFloatAgainstInputparameters(cpu_res_eo, params);
	BOOST_MESSAGE("Test done");
}

void test_m_wilson_compare_noneo_eo(std::string inputfile)
{
	test_m_fermion_compare_noneo_eo(inputfile, 0);
}

void test_m_tm_plus_compare_noneo_eo(std::string inputfile)
{
	test_m_fermion_compare_noneo_eo(inputfile, 1);
}

void test_m_tm_minus_compare_noneo_eo(std::string inputfile)
{
	test_m_fermion_compare_noneo_eo(inputfile, 2);
}

BOOST_AUTO_TEST_SUITE(M_WILSON_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_1)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_2)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_3)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_4)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_5)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_1)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_2)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_3)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_4)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_5)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_MINUS_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_1)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_2)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_3)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_4)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_5)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_SUITE_END()
