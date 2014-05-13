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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_SPINORS
#include <boost/test/unit_test.hpp>

#include "../hardware/code/kernelTester.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/prng.hpp"
#include "../device.hpp"
#include "spinors.hpp"
#include "complex.hpp"

class SpinorTester : public KernelTester {
public:
	SpinorTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
		inputfile(getSpecificInputfile(inputfileIn)), KernelTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues)
		{
		//todo: this object should be a member of KernelTester!
		meta::Inputparameters parameters = createParameters(inputfile);

		system = new hardware::System(parameters);
		device = system->get_devices()[0];

		code = device->get_spinor_code();
		
		NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(parameters);
		NUM_ELEMENTS_EO = hardware::code::get_eoprec_spinorfieldsize(parameters);
		(parameters.get_solver() == meta::Inputparameters::cg) ? useRandom = false : useRandom =true;
		alpha_host = {parameters.get_beta(), parameters.get_rho()};
		beta_host = {parameters.get_kappa(), parameters.get_mu()};
	}
protected:
	std::string inputfile;
	
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		return "spinors/" + inputfileIn;
	}
	
	spinor * createSpinorfield(size_t numberOfElements, int seed = 123456)
	{
		spinor * in;
		in = new spinor[numberOfElements];
		if(useRandom)
		{
			fill_with_random(in, numberOfElements, seed);
		}
		else
		{
			fill_with_one(in, numberOfElements);
		}
		BOOST_REQUIRE(in);
		return in;
	}
	
	void fill_with_one(spinor * in, int size)
	{
		for(int i = 0; i < size; ++i) {
			in[i].e0.e0 = hmc_complex_one;
			in[i].e0.e1 = hmc_complex_one;
			in[i].e0.e2 = hmc_complex_one;
			in[i].e1.e0 = hmc_complex_one;
			in[i].e1.e1 = hmc_complex_one;
			in[i].e1.e2 = hmc_complex_one;
			in[i].e2.e0 = hmc_complex_one;
			in[i].e2.e1 = hmc_complex_one;
			in[i].e2.e2 = hmc_complex_one;
			in[i].e3.e0 = hmc_complex_one;
			in[i].e3.e1 = hmc_complex_one;
			in[i].e3.e2 = hmc_complex_one;
		}
		return;
	}
	
	void fill_with_random(spinor * in, int size, int seed)
	{
		prng_init(seed);
		for(int i = 0; i < size; ++i) {
			in[i].e0.e0.re = prng_double();
			in[i].e0.e1.re = prng_double();
			in[i].e0.e2.re = prng_double();
			in[i].e1.e0.re = prng_double();
			in[i].e1.e1.re = prng_double();
			in[i].e1.e2.re = prng_double();
			in[i].e2.e0.re = prng_double();
			in[i].e2.e1.re = prng_double();
			in[i].e2.e2.re = prng_double();
			in[i].e3.e0.re = prng_double();
			in[i].e3.e1.re = prng_double();
			in[i].e3.e2.re = prng_double();

			in[i].e0.e0.im = prng_double();
			in[i].e0.e1.im = prng_double();
			in[i].e0.e2.im = prng_double();
			in[i].e1.e0.im = prng_double();
			in[i].e1.e1.im = prng_double();
			in[i].e1.e2.im = prng_double();
			in[i].e2.e0.im = prng_double();
			in[i].e2.e1.im = prng_double();
			in[i].e2.e2.im = prng_double();
			in[i].e3.e0.im = prng_double();
			in[i].e3.e1.im = prng_double();
			in[i].e3.e2.im = prng_double();
		}
		return;
	}
	
	const hardware::System * system;
	hardware::Device * device;

	const hardware::code::Spinors * code;
	
	size_t NUM_ELEMENTS_SF;
	size_t NUM_ELEMENTS_EO;
	bool useRandom;
	hmc_complex alpha_host;
	hmc_complex beta_host;
};

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	BOOST_CHECK_NO_THROW(   SpinorTester spinorTester("build all kernels", "spinors_build_input_1") );
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	BOOST_CHECK_NO_THROW(   SpinorTester spinorTester("build all kernels", "spinors_build_input_2") );
}

BOOST_AUTO_TEST_SUITE_END()


//some functionality
#include "../../tests/test_util.h"

void fill_with_one(spinor * in, int size)
{
	for(int i = 0; i < size; ++i) {
		in[i].e0.e0 = hmc_complex_one;
		in[i].e0.e1 = hmc_complex_one;
		in[i].e0.e2 = hmc_complex_one;
		in[i].e1.e0 = hmc_complex_one;
		in[i].e1.e1 = hmc_complex_one;
		in[i].e1.e2 = hmc_complex_one;
		in[i].e2.e0 = hmc_complex_one;
		in[i].e2.e1 = hmc_complex_one;
		in[i].e2.e2 = hmc_complex_one;
		in[i].e3.e0 = hmc_complex_one;
		in[i].e3.e1 = hmc_complex_one;
		in[i].e3.e2 = hmc_complex_one;
	}
	return;
}

void fill_with_zero(spinor * in, int size)
{
	for(int i = 0; i < size; ++i) {
		in[i].e0.e0 = hmc_complex_zero;
		in[i].e0.e1 = hmc_complex_zero;
		in[i].e0.e2 = hmc_complex_zero;
		in[i].e1.e0 = hmc_complex_zero;
		in[i].e1.e1 = hmc_complex_zero;
		in[i].e1.e2 = hmc_complex_zero;
		in[i].e2.e0 = hmc_complex_zero;
		in[i].e2.e1 = hmc_complex_zero;
		in[i].e2.e2 = hmc_complex_zero;
		in[i].e3.e0 = hmc_complex_zero;
		in[i].e3.e1 = hmc_complex_zero;
		in[i].e3.e2 = hmc_complex_zero;
	}
	return;
}

void fill_with_one_eo(spinor * in, int size, bool eo, meta::Inputparameters &params)
{
	int ns = params.get_nspace();
	int nt = params.get_ntime();
	int x, y, z, t;
	for (x = 0; x < ns; x++) {
		for (y = 0; y < ns; y++) {
			for (z = 0; z < ns; z++) {
				for (t = 0; t < nt; t++) {
					int coord[4];
					coord[0] = t;
					coord[1] = x;
					coord[2] = y;
					coord[3] = z;
					int nspace =  get_nspace(coord, params);
					int global_pos = get_global_pos(nspace, t, params);
					if (global_pos > size)
						break;
					hmc_complex content;
					if ((x + y + z + t) % 2 == 0) {
						if (eo)
							content = hmc_complex_one;
						else
							content = hmc_complex_zero;
					} else {
						if (eo)
							content = hmc_complex_zero;
						else
							content = hmc_complex_one;
					}

					in[global_pos].e0.e0 = content;
					in[global_pos].e0.e1 = content;
					in[global_pos].e0.e2 = content;
					in[global_pos].e1.e0 = content;
					in[global_pos].e1.e1 = content;
					in[global_pos].e1.e2 = content;
					in[global_pos].e2.e0 = content;
					in[global_pos].e2.e1 = content;
					in[global_pos].e2.e2 = content;
					in[global_pos].e3.e0 = content;
					in[global_pos].e3.e1 = content;
					in[global_pos].e3.e2 = content;
				}
			}
		}
	}
	return;
}

hmc_float count_eo(spinor * in, int size, bool eo, meta::Inputparameters &params)
{
	int ns = params.get_nspace();
	int nt = params.get_ntime();
	int x, y, z, t;
	hmc_float sum = 0.;
	for (x = 0; x < ns; x++) {
		for (y = 0; y < ns; y++) {
			for (z = 0; z < ns; z++) {
				for (t = 0; t < nt; t++) {
					int coord[4];
					coord[0] = t;
					coord[1] = x;
					coord[2] = y;
					coord[3] = z;
					int nspace =  get_nspace(coord, params);
					int global_pos = get_global_pos(nspace, t, params);
					if (global_pos > size)
						break;
					if (
					  ( eo == true && (x + y + z + t) % 2 == 0) ||
					  ( eo == false &&  (x + y + z + t) % 2 == 1 )
					) {
						int i = global_pos;
						sum +=
						  in[i].e0.e0.re + in[i].e0.e0.im
						  + in[i].e0.e1.re + in[i].e0.e1.im
						  + in[i].e0.e2.re + in[i].e0.e2.im
						  + in[i].e1.e0.re + in[i].e0.e0.im
						  + in[i].e1.e1.re + in[i].e0.e1.im
						  + in[i].e1.e2.re + in[i].e0.e2.im
						  + in[i].e2.e0.re + in[i].e0.e0.im
						  + in[i].e2.e1.re + in[i].e0.e1.im
						  + in[i].e2.e2.re + in[i].e0.e1.im
						  + in[i].e3.e0.re + in[i].e0.e0.im
						  + in[i].e3.e1.re + in[i].e0.e1.im
						  + in[i].e3.e2.re + in[i].e0.e1.im;
					} else {
						continue;
					}
				}
			}
		}
	}
	return sum;
}

hmc_float count_sf(spinor * in, int size)
{
	hmc_float sum = 0.;
	for (int i = 0; i < size; i++) {
		sum +=
		  in[i].e0.e0.re + in[i].e0.e0.im
		  + in[i].e0.e1.re + in[i].e0.e1.im
		  + in[i].e0.e2.re + in[i].e0.e2.im
		  + in[i].e1.e0.re + in[i].e1.e0.im
		  + in[i].e1.e1.re + in[i].e1.e1.im
		  + in[i].e1.e2.re + in[i].e1.e2.im
		  + in[i].e2.e0.re + in[i].e2.e0.im
		  + in[i].e2.e1.re + in[i].e2.e1.im
		  + in[i].e2.e2.re + in[i].e2.e2.im
		  + in[i].e3.e0.re + in[i].e3.e0.im
		  + in[i].e3.e1.re + in[i].e3.e1.im
		  + in[i].e3.e2.re + in[i].e3.e2.im;
	}
	return sum;
}

hmc_float calc_var(hmc_float in, hmc_float mean)
{
	return (in - mean) * (in - mean);
}

hmc_float calc_var_sf(spinor * in, int size, hmc_float sum)
{
	hmc_float var = 0.;
	for(int k = 0; k < size; k++) {
		var +=
		  calc_var( in[k].e0.e0.re , sum)
		  + calc_var( in[k].e0.e0.im , sum)
		  + calc_var( in[k].e0.e1.re , sum)
		  + calc_var( in[k].e0.e1.im , sum)
		  + calc_var( in[k].e0.e2.re , sum)
		  + calc_var( in[k].e0.e2.im , sum)
		  + calc_var( in[k].e1.e0.re , sum)
		  + calc_var( in[k].e1.e0.im , sum)
		  + calc_var( in[k].e1.e1.re , sum)
		  + calc_var( in[k].e1.e1.im , sum)
		  + calc_var( in[k].e1.e2.re , sum)
		  + calc_var( in[k].e1.e2.im , sum)
		  + calc_var( in[k].e2.e0.re , sum)
		  + calc_var( in[k].e2.e0.im , sum)
		  + calc_var( in[k].e2.e1.re , sum)
		  + calc_var( in[k].e2.e1.im , sum)
		  + calc_var( in[k].e2.e2.re , sum)
		  + calc_var( in[k].e2.e2.im , sum)
		  + calc_var( in[k].e3.e0.re , sum)
		  + calc_var( in[k].e3.e0.im , sum)
		  + calc_var( in[k].e3.e1.re , sum)
		  + calc_var( in[k].e3.e1.im , sum)
		  + calc_var( in[k].e3.e2.re , sum)
		  + calc_var( in[k].e3.e2.im , sum);
	}
	return var;
}


void fill_with_random(spinor * in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		in[i].e0.e0.re = prng_double();
		in[i].e0.e1.re = prng_double();
		in[i].e0.e2.re = prng_double();
		in[i].e1.e0.re = prng_double();
		in[i].e1.e1.re = prng_double();
		in[i].e1.e2.re = prng_double();
		in[i].e2.e0.re = prng_double();
		in[i].e2.e1.re = prng_double();
		in[i].e2.e2.re = prng_double();
		in[i].e3.e0.re = prng_double();
		in[i].e3.e1.re = prng_double();
		in[i].e3.e2.re = prng_double();

		in[i].e0.e0.im = prng_double();
		in[i].e0.e1.im = prng_double();
		in[i].e0.e2.im = prng_double();
		in[i].e1.e0.im = prng_double();
		in[i].e1.e1.im = prng_double();
		in[i].e1.e2.im = prng_double();
		in[i].e2.e0.im = prng_double();
		in[i].e2.e1.im = prng_double();
		in[i].e2.e2.im = prng_double();
		in[i].e3.e0.im = prng_double();
		in[i].e3.e1.im = prng_double();
		in[i].e3.e2.im = prng_double();
	}
	return;
}

void fill_with_random(spinor * in, int size)
{
	fill_with_random(in, size, 123456);
}

BOOST_AUTO_TEST_SUITE(GLOBAL_SQUARENORM)

	class SquarenormTester: public SpinorTester
	{
	public:
		SquarenormTester(std::string inputfile):
			SpinorTester("global_squarenorm", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(NUM_ELEMENTS_SF, device);
				in.load(createSpinorfield(NUM_ELEMENTS_SF));
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);

				code->set_float_to_global_squarenorm_device(&in, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_1 )
	{
		SquarenormTester("/global_squarenorm_input_1");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_2 )
	{
		SquarenormTester("/global_squarenorm_input_2");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_1 )
	{
		SquarenormTester("/global_squarenorm_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_2 )
	{
		SquarenormTester("/global_squarenorm_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( GLOBAL_SQUARENORM_REDUCTION_3 )
	{
		SquarenormTester("/global_squarenorm_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GLOBAL_SQUARENORM_EO)

	class SquarenormEvenOddTester: public SpinorTester
	{
	public:
		SquarenormEvenOddTester(std::string inputfile):
			SpinorTester("global_squarenorm_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(NUM_ELEMENTS_EO, device);
				in.load(createSpinorfield(NUM_ELEMENTS_EO));
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);

				code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};
	
	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_1 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_2 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_1 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_2 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( SQUARENORM_EO_REDUCTION_3 )
	{
		SquarenormEvenOddTester squarenormEoTester("global_squarenorm_eo_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT)

	class ScalarProductTester: public SpinorTester
	{
	public:
		ScalarProductTester(std::string inputfile):
			SpinorTester("scalar_product", inputfile, 2)
			{
				const hardware::buffers::Plain<spinor> in(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Plain<spinor> in2(NUM_ELEMENTS_SF, device);
				in.load(createSpinorfield(NUM_ELEMENTS_SF, 123));
				in2.load(createSpinorfield(NUM_ELEMENTS_SF, 456));
				hardware::buffers::Plain<hmc_complex> sqnorm(1, device);

				code->set_complex_to_scalar_product_device(&in, &in2, &sqnorm);
				hmc_complex resultTmp;
				sqnorm.dump(&resultTmp);
				
				kernelResult[0] = resultTmp.re;
				kernelResult[1] = resultTmp.im;
			}
	};

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_1 )
	{
		ScalarProductTester scalarProductTester("scalar_product_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_2 )
	{
		ScalarProductTester scalarProductTester("scalar_product_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_1 )
	{
		ScalarProductTester scalarProductTester("scalar_product_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_2 )
	{
		ScalarProductTester scalarProductTester("scalar_product_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_REDUCTION_3 )
	{
		ScalarProductTester scalarProductTester("scalar_product_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SCALAR_PRODUCT_EO)

	class ScalarProductEvenOddTester: public SpinorTester
	{
	public:
		ScalarProductEvenOddTester(std::string inputfile):
			SpinorTester("scalar_product_eo", inputfile, 2)
			{
				const hardware::buffers::Spinor in(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor in2(NUM_ELEMENTS_EO, device);
				in.load(createSpinorfield(NUM_ELEMENTS_EO, 123));
				in2.load(createSpinorfield(NUM_ELEMENTS_EO, 456));
				hardware::buffers::Plain<hmc_complex> sqnorm(1, device);

				code->set_complex_to_scalar_product_eoprec_device(&in, &in2, &sqnorm);
				hmc_complex resultTmp;
				sqnorm.dump(&resultTmp);
				
				kernelResult[0] = resultTmp.re;
				kernelResult[1] = resultTmp.im;
			}
	};
	
	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_1 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_2 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_1 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_reduction_input_1");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_2 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_reduction_input_2");
	}

	BOOST_AUTO_TEST_CASE( SCALAR_PRODUCT_EO_REDUCTION_3 )
	{
		ScalarProductEvenOddTester scalarProductEoTester("scalar_product_eo_reduction_input_3");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_AND_ZERO)

	class ColdAndZeroTester: public SpinorTester
	{
	public:
		ColdAndZeroTester(std::string inputfile, bool switcher):
			SpinorTester("cold or zero", inputfile, 2)
			{
				const hardware::buffers::Plain<spinor> in(NUM_ELEMENTS_SF, device);
				in.load(createSpinorfield(NUM_ELEMENTS_SF));
				hardware::buffers::Plain<double> sqnorm(1, device);

				(switcher) ? code->set_spinorfield_cold_device(&in) : 	code->set_zero_spinorfield_device(&in);
				code->set_float_to_global_squarenorm_device(&in, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( COLD_1 )
	{
		ColdAndZeroTester tester("cold_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( ZERO_1 )
	{
		ColdAndZeroTester tester("zero_input_1", false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(COLD_AND_ZERO_EO)

	class ColdAndZeroEvenOddTester: public SpinorTester
	{
	public:
		ColdAndZeroEvenOddTester(std::string inputfile, bool switcher):
			SpinorTester("cold or zero eo", inputfile, 2)
			{
				const hardware::buffers::Spinor in(NUM_ELEMENTS_EO, device);
				in.load(createSpinorfield(NUM_ELEMENTS_EO));
				hardware::buffers::Plain<double> sqnorm(1, device);

				(switcher) ? code->set_eoprec_spinorfield_cold_device(&in) : 	code->set_zero_spinorfield_eoprec_device(&in);
				code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};
	
	BOOST_AUTO_TEST_CASE( COLD_EO_1 )
	{
		ColdAndZeroEvenOddTester tester("cold_eo_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( ZERO_EO_1 )
	{
		ColdAndZeroEvenOddTester tester("zero_eo_input_1",  false);
	}
	
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX)

	class SaxTester: public SpinorTester
	{
	public:
		SaxTester(std::string inputfile):
			SpinorTester("sax", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Plain<spinor> out(NUM_ELEMENTS_SF, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(NUM_ELEMENTS_SF, 123));
				alpha.load(&alpha_host);

				code->sax_device(&in, &alpha, &out);
				code->set_float_to_global_squarenorm_device(&out, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_1 )
	{
		SaxTester tester ("sax_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAX_2 )
	{
		SaxTester tester ("sax_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAX_3 )
	{
		SaxTester tester ("sax_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAX_4 )
	{
		SaxTester tester ("sax_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAX_5 )
	{
		SaxTester tester ("sax_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAX_6 )
	{
		SaxTester tester ("sax_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAX_7 )
	{
		SaxTester tester ("sax_input_7");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAX_EO)

	class SaxEvenOddTester: public SpinorTester
	{
	public:
		SaxEvenOddTester(std::string inputfile):
			SpinorTester("sax_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor out(NUM_ELEMENTS_EO, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(NUM_ELEMENTS_EO, 123));
				alpha.load(&alpha_host);

				code->sax_eoprec_device(&in, &alpha, &out);
				code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( SAX_EO_1 )
	{
		SaxEvenOddTester tester("sax_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_2 )
	{
		SaxEvenOddTester tester("sax_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_3 )
	{
		SaxEvenOddTester tester("sax_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_4 )
	{
		SaxEvenOddTester tester("sax_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_5 )
	{
		SaxEvenOddTester tester("sax_eo_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_6 )
	{
		SaxEvenOddTester tester("sax_eo_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAX_EO_7 )
	{
		SaxEvenOddTester tester("sax_eo_input_7");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY)

	class SaxpyTester: public SpinorTester
	{
	public:
		SaxpyTester(std::string inputfile, bool switcher):
			SpinorTester("saxpy", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Plain<spinor> in2(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Plain<spinor> out(NUM_ELEMENTS_SF, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(NUM_ELEMENTS_SF, 123));
				in2.load(createSpinorfield(NUM_ELEMENTS_SF, 456));
				alpha.load(&alpha_host);

				(switcher) ? code->saxpy_device(&in, &in2, &alpha, &out) : code->saxpy_device(&in, &in2, alpha_host, &out);
				code->set_float_to_global_squarenorm_device(&out, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		SaxpyTester tester("saxpy_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		SaxpyTester tester("saxpy_input_2", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_3 )
	{
		SaxpyTester tester("saxpy_input_3", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_4 )
	{
		SaxpyTester tester("saxpy_input_4", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_5 )
	{
		SaxpyTester tester("/saxpy_input_5", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_6 )
	{
		SaxpyTester tester("saxpy_input_6", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_7 )
	{
		SaxpyTester tester("saxpy_input_7", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_8 )
	{
		SaxpyTester tester("saxpy_input_8", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_9 )
	{
		SaxpyTester tester("saxpy_input_9", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_10 )
	{
		SaxpyTester tester("saxpy_input_10", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_11 )
	{
		SaxpyTester tester("saxpy_input_11", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_12 )
	{
		SaxpyTester tester("saxpy_input_12", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_13 )
	{
		SaxpyTester tester("saxpy_input_13", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_14 )
	{
		SaxpyTester tester("saxpy_input_14", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_1 )
	{
		SaxpyTester tester("saxpy_input_1", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_2 )
	{
		SaxpyTester tester("saxpy_input_2", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_3 )
	{
		SaxpyTester tester("saxpy_input_3", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_4 )
	{
		SaxpyTester tester("saxpy_input_4", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_5 )
	{
		SaxpyTester tester("saxpy_input_5", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_6 )
	{
		SaxpyTester tester("saxpy_input_6", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_7 )
	{
		SaxpyTester tester("saxpy_input_7", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_8 )
	{
		SaxpyTester tester("saxpy_input_8", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_9 )
	{
		SaxpyTester tester("saxpy_input_9", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_10 )
	{
		SaxpyTester tester("saxpy_input_10", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_11 )
	{
		SaxpyTester tester("saxpy_input_11", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_12 )
	{
		SaxpyTester tester("saxpy_input_12", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_13 )
	{
		SaxpyTester tester("saxpy_input_13", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_14 )
	{
		SaxpyTester tester("saxpy_input_14", false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_EO)

	class SaxpyEvenOddTester: public SpinorTester
	{
	public:
		SaxpyEvenOddTester(std::string inputfile, bool switcher):
			SpinorTester("saxpy_eo or saxpy_arg_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor in2(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor out(NUM_ELEMENTS_EO, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);

				in.load(createSpinorfield(NUM_ELEMENTS_EO, 123));
				in2.load(createSpinorfield(NUM_ELEMENTS_EO, 456));
				alpha.load(&alpha_host);

				(switcher) ? code->saxpy_eoprec_device(&in, &in2, &alpha, &out) : code->saxpy_eoprec_device(&in, &in2, alpha_host, &out);
				code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXPY_EO_1 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_1", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_2 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_2", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_3 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_3", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_4 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_4", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_5 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_5", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_6 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_6", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_7 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_7", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_8 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_8", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_9 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_9", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_10 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_10", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_11 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_11", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_12 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_12", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_13 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_13", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_EO_14 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_14", true);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_1 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_1", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_2 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_2", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_3 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_3", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_4 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_4", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_5 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_5", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_6 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_6", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_7 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_7", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_8 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_8", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_9 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_9", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_10 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_10", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_11 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_11", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_12 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_12", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_13 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_13", false);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_ARG_EO_14 )
	{
		SaxpyEvenOddTester tester("saxpy_eo_input_14", false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXSBYPZ)

	class SaxsbypzTester: public SpinorTester
	{
	public:
		SaxsbypzTester(std::string inputfile):
			SpinorTester("saxsbypz", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> in(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Plain<spinor> in2(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Plain<spinor> in3(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Plain<spinor> out(NUM_ELEMENTS_SF, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				hardware::buffers::Plain<hmc_complex> beta(1, device);

				in.load(createSpinorfield(NUM_ELEMENTS_SF, 123));
				in2.load(createSpinorfield(NUM_ELEMENTS_SF, 456));
				in3.load(createSpinorfield(NUM_ELEMENTS_SF, 789));
				alpha.load(&alpha_host);
				beta.load(&beta_host);

				code->saxsbypz_device(&in, &in2, &in3, &alpha, &beta, &out);
				code->set_float_to_global_squarenorm_device(&out, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_1 )
	{
		SaxsbypzTester tester("saxsbypz_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_2 )
	{
		SaxsbypzTester tester("saxsbypz_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_3 )
	{
		SaxsbypzTester tester("saxsbypz_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_4 )
	{
		SaxsbypzTester tester("saxsbypz_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_5 )
	{
		SaxsbypzTester tester("saxsbypz_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_6 )
	{
		SaxsbypzTester tester("saxsbypz_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_7 )
	{
		SaxsbypzTester tester("saxsbypz_input_7");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_8 )
	{
		SaxsbypzTester tester("saxsbypz_input_8");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_9 )
	{
		SaxsbypzTester tester("saxsbypz_input_9");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_10 )
	{
		SaxsbypzTester tester("saxsbypz_input_10");
	}

BOOST_AUTO_TEST_SUITE_END()

// void test_saxsbypz_eo(std::string inputfile)
// {
// 	using namespace hardware::buffers;
// 
// 	std::string kernelName;
// 	kernelName = "saxsbypz_eo";
// 	printKernelInfo(kernelName);
// 	logger.info() << "Init device";
// 	meta::Inputparameters params = create_parameters(inputfile);
// 	hardware::System system(params);
// 	auto * device = system.get_devices().at(0)->get_spinor_code();
// 
// 	logger.info() << "Fill buffers...";
// 	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
// 	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
// 	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
// 	const Spinor in3(NUM_ELEMENTS_SF, device->get_device());
// 	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
// 	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
// 	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
// 	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());
// 
// 	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
// 	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};
// 	logger.info() << "Use alpha = (" << alpha_host.re << "," << alpha_host.im << ")";
// 	logger.info() << "Use beta = (" << beta_host.re << "," << beta_host.im << ")";
// 
// 	spinor * in;
// 	spinor * in2;
// 	spinor * in3;
// 	in = new spinor[NUM_ELEMENTS_SF];
// 	in2 = new spinor[NUM_ELEMENTS_SF];
// 	in3 = new spinor[NUM_ELEMENTS_SF];
// 	//use the variable use_cg to switch between cold and random input sf
// 	if(params.get_solver() == meta::Inputparameters::cg) {
// 		fill_with_one(in, NUM_ELEMENTS_SF);
// 		fill_with_one(in2, NUM_ELEMENTS_SF);
// 		fill_with_one(in3, NUM_ELEMENTS_SF);
// 	} else {
// 		fill_with_random(in, NUM_ELEMENTS_SF, 123);
// 		fill_with_random(in2, NUM_ELEMENTS_SF, 456);
// 		fill_with_random(in3, NUM_ELEMENTS_SF, 789);
// 	}
// 	BOOST_REQUIRE(in);
// 	BOOST_REQUIRE(in2);
// 	BOOST_REQUIRE(in3);
// 
// 	in.load(in);
// 	in2.load(in2);
// 	in3.load(in3);
// 	alpha.load(&alpha_host);
// 	beta.load(&beta_host);
// 
// 	auto spinor_code = device->get_device()->get_spinor_code();
// 
// 	logger.info() << "Run kernel";
// 	device->saxsbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out);
// 
// 	logger.info() << "result:";
// 	hmc_float cpu_res;
// 	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
// 	sqnorm.dump(&cpu_res);
// 	logger.info() << cpu_res;
// 
// 	testFloatAgainstInputparameters(cpu_res, params);
// 	BOOST_MESSAGE("Test done");
// }
// 
// void test_convert_to_eo(std::string inputfile)
// {
// 	using namespace hardware::buffers;
// 
// 	std::string kernelName;
// 	kernelName = "convert_to_eo";
// 	printKernelInfo(kernelName);
// 	logger.info() << "Init device";
// 	meta::Inputparameters params = create_parameters(inputfile);
// 	hardware::System system(params);
// 	auto * device = system.get_devices().at(0)->get_spinor_code();
// 
// 	logger.info() << "Fill buffers...";
// 	size_t NUM_ELEMENTS_EO = hardware::code::get_eoprec_spinorfieldsize(params);
// 	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
// 	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
// 	const Spinor in2(NUM_ELEMENTS_EO, device->get_device());
// 	const Spinor in3(NUM_ELEMENTS_EO, device->get_device());
// 	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
// 	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
// 
// 	spinor * in;
// 	in = new spinor[NUM_ELEMENTS_SF];
// 	if(params.get_read_multiple_configs() )
// 		fill_with_one_eo(in, NUM_ELEMENTS_SF, true, params);
// 	else
// 		fill_with_one_eo(in, NUM_ELEMENTS_SF, false, params);
// 	BOOST_REQUIRE(in);
// 
// 	in.load(in);
// 
// 	auto spinor_code = device->get_device()->get_spinor_code();
// 
// 	logger.info() << "Run kernel";
// 	device->convert_to_eoprec_device(&in2, &in3, &in);
// 
// 	logger.info() << "result:";
// 	hmc_float cpu_res;
// 	if(params.get_read_multiple_configs() ) {
// 		spinor_code->set_float_to_global_squarenorm_eoprec_device(&in3, &sqnorm);
// 		sqnorm.dump(&cpu_res);
// 		logger.info() << cpu_res;
// 		//CP: this must be zero since only the even sites should be filled!
// 		BOOST_REQUIRE_CLOSE(cpu_res, 0., params.get_solver_prec());
// 		spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
// 		sqnorm.dump(&cpu_res);
// 		logger.info() << cpu_res;
// 	} else {
// 		spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
// 		sqnorm.dump(&cpu_res);
// 		logger.info() << cpu_res;
// 		//CP: this must be zero since only the odd sites should be filled!
// 		BOOST_REQUIRE_CLOSE(cpu_res, 0., params.get_solver_prec());
// 		spinor_code->set_float_to_global_squarenorm_eoprec_device(&in3, &sqnorm);
// 		sqnorm.dump(&cpu_res);
// 		logger.info() << cpu_res;
// 	}
// 
// 	testFloatAgainstInputparameters(cpu_res, params);
// 	BOOST_MESSAGE("Test done");
// }
// 
// void test_convert_from_eo(std::string inputfile)
// {
// 	using namespace hardware::buffers;
// 
// 	std::string kernelName;
// 	kernelName = "convert_from_eo";
// 	printKernelInfo(kernelName);
// 	logger.info() << "Init device";
// 	meta::Inputparameters params = create_parameters(inputfile);
// 	hardware::System system(params);
// 	auto * device = system.get_devices().at(0)->get_spinor_code();
// 
// 	logger.info() << "Fill buffers...";
// 	size_t NUM_ELEMENTS_EO = hardware::code::get_eoprec_spinorfieldsize(params);
// 	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
// 	const Spinor in2(NUM_ELEMENTS_EO, device->get_device());
// 	const Spinor in3(NUM_ELEMENTS_EO, device->get_device());
// 	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
// 	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
// 
// 	spinor * eo1;
// 	spinor * eo2;
// 	spinor * out;
// 	eo1 = new spinor[NUM_ELEMENTS_EO];
// 	eo2 = new spinor[NUM_ELEMENTS_EO];
// 	out = new spinor[NUM_ELEMENTS_SF];
// 	if(params.get_read_multiple_configs() ) {
// 		fill_with_one(eo1, NUM_ELEMENTS_EO);
// 		fill_with_zero(eo2, NUM_ELEMENTS_EO);
// 	} else {
// 		fill_with_zero(eo1, NUM_ELEMENTS_EO);
// 		fill_with_one(eo2, NUM_ELEMENTS_EO);
// 	}
// 	BOOST_REQUIRE(eo1);
// 	BOOST_REQUIRE(eo2);
// 
// 	in2.load(eo1);
// 	in3.load(eo2);
// 
// 	auto spinor_code = device->get_device()->get_spinor_code();
// 
// 	logger.info() << "Run kernel";
// 	device->convert_from_eoprec_device(&in2, &in3, &out);
// 
// 	out.dump(out);
// 
// 	logger.info() << "result:";
// 	hmc_float cpu_res;
// 	if(params.get_read_multiple_configs() ) {
// 		cpu_res = count_eo(out, NUM_ELEMENTS_SF, true, params);
// 		logger.info() << cpu_res;
// 	} else {
// 		cpu_res = count_eo(out, NUM_ELEMENTS_SF, false, params);
// 		logger.info() << cpu_res;
// 	}
// 
// 	testFloatAgainstInputparameters(cpu_res, params);
// 	BOOST_MESSAGE("Test done");
// }
// 
// void test_gaussian(std::string inputfile)
// {
// 	using namespace hardware::buffers;
// 
// 	std::string kernelName;
// 	kernelName = "generate_gaussian_spinorfield";
// 	printKernelInfo(kernelName);
// 	logger.info() << "Init device";
// 	meta::Inputparameters params = create_parameters(inputfile);
// 	hardware::System system(params);
// 
// 	physics::PRNG prng(system);
// 	auto * device = system.get_devices().at(0)->get_spinor_code();
// 
// 	logger.info() << "Fill buffers...";
// 	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
// 	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
// 	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
// 
// 	//CP: run the kernel a couple of times times
// 	int iterations = params.get_integrationsteps(0);
// 
// 	spinor * out;
// 	out = new spinor[NUM_ELEMENTS_SF * iterations];
// 	BOOST_REQUIRE(out);
// 
// 	auto spinor_code = device->get_device()->get_spinor_code();
// 	auto prng_buf = prng.get_buffers().at(0);
// 
// 	hmc_float sum = 0;
// 	for (int i = 0; i < iterations; i++) {
// 		logger.info() << "Run kernel";
// 		device->generate_gaussian_spinorfield_device(&out, prng_buf);
// 		out.dump(&out[i * NUM_ELEMENTS_SF]);
// 		sum += count_sf(&out[i * NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
// 	}
// 	logger.info() << "result: mean";
// 	hmc_float cpu_res = 0.;
// 	sum = sum / iterations / NUM_ELEMENTS_SF / 24;
// 	cpu_res = sum;
// 	logger.info() << cpu_res;
// 
// 	if(params.get_read_multiple_configs()  == false) {
// 		//CP: calc std derivation
// 		hmc_float var = 0.;
// 		for (int i = 0; i < iterations; i++) {
// 			var += calc_var_sf(&out[i * NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
// 		}
// 		var = var / iterations / NUM_ELEMENTS_SF / 24;
// 
// 		cpu_res = sqrt(var);
// 		logger.info() << "result: variance";
// 		logger.info() << cpu_res;
// 	}
// 
// 
// 	testFloatSizeAgainstInputparameters(cpu_res, params);
// 	BOOST_MESSAGE("Test done");
// 
// }
// 
// void test_gaussian_eo(std::string inputfile)
// {
// 	using namespace hardware::buffers;
// 
// 	std::string kernelName;
// 	kernelName = "generate_gaussian_spinorfield_eo";
// 	printKernelInfo(kernelName);
// 	logger.info() << "Init device";
// 	meta::Inputparameters params = create_parameters(inputfile);
// 	hardware::System system(params);
// 
// 	physics::PRNG prng(system);
// 	auto * device = system.get_devices().at(0)->get_spinor_code();
// 
// 	logger.info() << "Fill buffers...";
// 	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
// 	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
// 
// 	//CP: run the kernel a couple of times times
// 	int iterations = params.get_integrationsteps(0);
// 
// 	spinor * out;
// 	out = new spinor[NUM_ELEMENTS_SF * iterations];
// 	BOOST_REQUIRE(out);
// 
// 	auto spinor_code = device->get_device()->get_spinor_code();
// 	auto prng_buf = prng.get_buffers().at(0);
// 
// 	hmc_float sum = 0;
// 	for (int i = 0; i < iterations; i++) {
// 		logger.info() << "Run kernel";
// 		device->generate_gaussian_spinorfield_eo_device(&out, prng_buf);
// 		out.dump(&out[i * NUM_ELEMENTS_SF]);
// 		sum += count_sf(&out[i * NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
// 	}
// 	logger.info() << "result: mean";
// 	hmc_float cpu_res = 0.;
// 	sum = sum / iterations / NUM_ELEMENTS_SF / 24;
// 	cpu_res = sum;
// 	logger.info() << cpu_res;
// 
// 	if(params.get_read_multiple_configs()  == false) {
// 		//CP: calc std derivation
// 		hmc_float var = 0.;
// 		for (int i = 0; i < iterations; i++) {
// 			var += calc_var_sf(&out[i * NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
// 		}
// 		var = var / iterations / NUM_ELEMENTS_SF / 24;
// 
// 		cpu_res = sqrt(var);
// 		logger.info() << "result: variance";
// 		logger.info() << cpu_res;
// 	}
// 
// 	testFloatSizeAgainstInputparameters(cpu_res, params);
// 	BOOST_MESSAGE("Test done");
// }
// 
// 
// BOOST_AUTO_TEST_SUITE(SAXSBYPZ_EO)
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_1 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_1");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_2 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_2");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_3 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_3");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_4 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_4");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_5 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_5");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_6 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_6");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_7 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_7");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_8 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_8");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_9 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_9");
// }
// 
// BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_10 )
// {
// 	test_saxsbypz_eo("/saxsbypz_eo_input_10");
// }
// 
// BOOST_AUTO_TEST_SUITE_END()
// 
// BOOST_AUTO_TEST_SUITE(CONVERT_EO)
// 
// BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
// {
// 	test_convert_to_eo("/convert_eo_input_1");
// }
// 
// BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
// {
// 	test_convert_to_eo("/convert_eo_input_2");
// }
// 
// BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
// {
// 	test_convert_from_eo("/convert_eo_input_1");
// }
// 
// BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
// {
// 	test_convert_from_eo("/convert_eo_input_2");
// }
// 
// BOOST_AUTO_TEST_SUITE_END()
// 
// BOOST_AUTO_TEST_SUITE(GAUSSIAN)
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_1 )
// {
// 	test_gaussian("/gaussian_input_1");
// }
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_2 )
// {
// 	test_gaussian("/gaussian_input_2");
// }
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_3 )
// {
// 	test_gaussian("/gaussian_input_3");
// }
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_4 )
// {
// 	test_gaussian("/gaussian_input_4");
// }
// 
// BOOST_AUTO_TEST_SUITE_END()
// 
// BOOST_AUTO_TEST_SUITE(GAUSSIAN_EO)
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_1 )
// {
// 	test_gaussian_eo("/gaussian_eo_input_1");
// }
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_2 )
// {
// 	test_gaussian_eo("/gaussian_eo_input_2");
// }
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_3 )
// {
// 	test_gaussian_eo("/gaussian_eo_input_3");
// }
// 
// BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_4 )
// {
// 	test_gaussian_eo("/gaussian_eo_input_4");
// }
// 
// BOOST_AUTO_TEST_SUITE_END()
