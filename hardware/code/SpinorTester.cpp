#include "SpinorTester.hpp"

void SpinorTester::fill_with_one(spinor * in, int size)
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

void SpinorTester::fill_with_random(spinor * in, int size, int seed)
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

spinor * SpinorTester::createSpinorfield(size_t numberOfElements, int seed)
{
  spinor * in;
  in = new spinor[numberOfElements];
  useRandom ? fill_with_random(in, numberOfElements, seed) : fill_with_one(in, numberOfElements);
  BOOST_REQUIRE(in);
  return in;
}

spinor * SpinorTester::createSpinorfieldWithOnesAndZerosDependingOnSiteParity()
{
  spinor * in;
  in = new spinor[NUM_ELEMENTS_SF];
  fill_with_one_eo(in, NUM_ELEMENTS_SF, evenOrOdd);
  return in;
}

std::string SpinorTester::getSpecificInputfile(std::string inputfileIn)
{
  return "spinors/" + inputfileIn;
}

void SpinorTester::fill_with_one_eo(spinor * in, int size, bool eo)
	{
		int x, y, z, t;
		hmc_complex content;
		int coord[4];
		bool parityOfSite;
		int nspace;
		int global_pos;

		for (x = 0; x < ns; x++) {
			for (y = 0; y < ns; y++) {
				for (z = 0; z < ns; z++) {
					for (t = 0; t < nt; t++) {
						coord[0] = t;
						coord[1] = x;
						coord[2] = y;
						coord[3] = z;
						nspace =  get_nspace(coord);
						global_pos = get_global_pos(nspace, t);
						if (global_pos > size)
							break;

						parityOfSite = (x + y + z + t) % 2 == 0;
						content = (parityOfSite) ?
							(eo ? hmc_complex_one : hmc_complex_zero) :
							(eo ? hmc_complex_zero : hmc_complex_one);

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

	hmc_float SpinorTester::count_sf(spinor * in, int size)
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

	hmc_float SpinorTester::calc_var(hmc_float in, hmc_float mean)
	{
		return (in - mean) * (in - mean);
	}

	hmc_float SpinorTester::calc_var_sf(spinor * in, int size, hmc_float sum)
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

void SpinorTester::calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<spinor> * in)
{
  code->set_float_to_global_squarenorm_device(in, doubleBuffer);
  doubleBuffer->dump(&kernelResult[0]);
}

void SpinorTester::calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::Spinor * in)
{
  code->set_float_to_global_squarenorm_eoprec_device(in, doubleBuffer);
  doubleBuffer->dump(&kernelResult[0]);
}
