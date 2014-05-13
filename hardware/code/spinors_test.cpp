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
		(parameters.get_read_multiple_configs() ) ? evenOrOdd = true : evenOrOdd = false;
		alpha_host = {parameters.get_beta(), parameters.get_rho()};
		beta_host = {parameters.get_kappa(), parameters.get_mu()};
		ns = parameters.get_nspace();
		nt = parameters.get_ntime();
		iterations = parameters.get_integrationsteps(0);
		parameters.get_read_multiple_configs() ? calcVariance=false : calcVariance = true;
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
	
	spinor * createSpinorfieldWithOnesAndZerosDependingOnSiteParity()
	{
		spinor * in;
		in = new spinor[NUM_ELEMENTS_SF];
		fill_with_one_eo(in, NUM_ELEMENTS_SF, evenOrOdd);
		return in;
	}
	
	//todo: use the fct. from geometry.h here!
	int get_nspace(int* coord)
	{
		int n = 0;
		n = ns * ns * coord[3] + ns * coord[2] + coord[1];
		return n;
	}
	
	int get_global_pos(int spacepos, int t)
	{
		return spacepos + ns*ns*ns * t;
	}
	
	void fill_with_one_eo(spinor * in, int size, bool eo)
	{
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
						int nspace =  get_nspace(coord);
						int global_pos = get_global_pos(nspace, t);
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
	
	const hardware::System * system;
	hardware::Device * device;

	const hardware::code::Spinors * code;
	
	size_t NUM_ELEMENTS_SF;
	size_t NUM_ELEMENTS_EO;
	bool useRandom;
	bool evenOrOdd;
	bool calcVariance;
	hmc_complex alpha_host;
	hmc_complex beta_host;
	int ns;
	int nt;
	int iterations;
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

BOOST_AUTO_TEST_SUITE(SAXSBYPZ_EO)

	class SaxsbypzEvenOddTester: public SpinorTester
	{
	public:
		SaxsbypzEvenOddTester(std::string inputfile):
			SpinorTester("saxsbypz_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor in(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor in2(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor in3(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor out(NUM_ELEMENTS_EO, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);
				hardware::buffers::Plain<hmc_complex> alpha(1, device);
				hardware::buffers::Plain<hmc_complex> beta(1, device);

				in.load(createSpinorfield(NUM_ELEMENTS_EO, 123));
				in2.load(createSpinorfield(NUM_ELEMENTS_EO, 456));
				in3.load(createSpinorfield(NUM_ELEMENTS_EO, 789));
				alpha.load(&alpha_host);
				beta.load(&beta_host);

				code->saxsbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out);
				code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
			}
	};

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_1 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_2 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_3 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_4 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_4");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_5 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_5");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_6 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_6");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_7 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_7");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_8 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_8");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_9 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_9");
	}

	BOOST_AUTO_TEST_CASE( SAXSBYPZ_EO_10 )
	{
		SaxsbypzEvenOddTester tester("saxsbypz_eo_input_10");
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(CONVERT_EO)

	class ConvertEvenOddTester: public SpinorTester
	{
	public:
		ConvertEvenOddTester(std::string inputfile):
			SpinorTester("convert_eo", inputfile, 2)
			{
				const hardware::buffers::Plain<spinor> in(NUM_ELEMENTS_SF, device);
				const hardware::buffers::Spinor in2(NUM_ELEMENTS_EO, device);
				const hardware::buffers::Spinor in3(NUM_ELEMENTS_EO, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);

				in.load( createSpinorfieldWithOnesAndZerosDependingOnSiteParity() );

				code->convert_to_eoprec_device(&in2, &in3, &in);
				code->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
				sqnorm.dump(&kernelResult[0]);
				code->set_float_to_global_squarenorm_eoprec_device(&in3, &sqnorm);
				sqnorm.dump(&kernelResult[1]);
				
				if (evenOrOdd)
				{
					referenceValue[0] = NUM_ELEMENTS_EO*12.; 
					referenceValue[1] = 0.;
				}
				else
				{ 
					referenceValue[1] = NUM_ELEMENTS_EO*12.; 
					referenceValue[0] = 0.;
				}
			}
	};

	BOOST_AUTO_TEST_CASE( CONVERT_EO_1 )
	{
		ConvertEvenOddTester tester("convert_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_2 )
	{
		ConvertEvenOddTester tester("convert_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_3 )
	{
		ConvertEvenOddTester tester("convert_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( CONVERT_EO_4 )
	{
		ConvertEvenOddTester tester("convert_eo_input_2");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN)

	class GaussianTester: public SpinorTester
	{
	public:
		GaussianTester(std::string inputfile):
			SpinorTester("generate_gaussian_spinorfield", inputfile, 1)
			{
				const hardware::buffers::Plain<spinor> out(NUM_ELEMENTS_SF, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);

				spinor * outHost;
				outHost = new spinor[NUM_ELEMENTS_SF * iterations];
				BOOST_REQUIRE(out);
				
// 				physics::PRNG prng(*system);
// 				auto prng_buf = prng.get_buffers().at(0);
		
				double sum = 0;
				for (int i = 0; i < iterations; i++) {
// 					code->generate_gaussian_spinorfield_device(&out, prng_buf);
					out.dump(&outHost[i * NUM_ELEMENTS_SF]);
					sum += count_sf(&outHost[i * NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
				}
				kernelResult[0] = sum / iterations / NUM_ELEMENTS_SF / 24;

				if(calcVariance) 
				{
					double var = 0.;
					for (int i = 0; i < iterations; i++) {
						var += calc_var_sf(&outHost[i * NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
					}
					var = var / iterations / NUM_ELEMENTS_SF / 24;

					kernelResult[0] = sqrt(var);
				}
			}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_1 )
	{
		GaussianTester tester("gaussian_input_1");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_2 )
	{
		GaussianTester tester("gaussian_input_2");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_3 )
	{
		GaussianTester tester("gaussian_input_3");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_4 )
	{
		GaussianTester tester("gaussian_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GAUSSIAN_EO)

	class GaussianEvenOddTester: public SpinorTester
	{
	public:
		GaussianEvenOddTester(std::string inputfile):
			SpinorTester("generate_gaussian_spinorfield_eo", inputfile, 1)
			{
				const hardware::buffers::Spinor out(NUM_ELEMENTS_EO, device);
				hardware::buffers::Plain<hmc_float> sqnorm(1, device);

				spinor * outHost;
				outHost = new spinor[NUM_ELEMENTS_EO * iterations];
				BOOST_REQUIRE(out);
				
// 				physics::PRNG prng(*system);
// 				auto prng_buf = prng.get_buffers().at(0);

				double sum = 0;
				for (int i = 0; i < iterations; i++) {
// 					code->generate_gaussian_spinorfield_eo_device(&out, prng_buf);
					out.dump(&outHost[i * NUM_ELEMENTS_EO]);
					sum += count_sf(&outHost[i * NUM_ELEMENTS_EO], NUM_ELEMENTS_EO);
				}
				kernelResult[0] = sum / iterations / NUM_ELEMENTS_SF / 24;

				if(calcVariance) 
				{
					double var = 0.;
					for (int i = 0; i < iterations; i++) {
						var += calc_var_sf(&outHost[i * NUM_ELEMENTS_EO], NUM_ELEMENTS_EO, sum);
					}
					var = var / iterations / NUM_ELEMENTS_EO / 24;

					logger.info() << "variance:\t" << sqrt(var);
				}
			}
	};

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_1 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_1");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_2 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_2");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_3 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_3");
	}

	BOOST_AUTO_TEST_CASE( GAUSSIAN_EO_4 )
	{
		GaussianEvenOddTester tester("gaussian_eo_input_4");
	}

BOOST_AUTO_TEST_SUITE_END()
