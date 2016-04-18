/*
 * Copyright 2014 Christopher Pinke
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

#include "SpinorTester.hpp"
#include "../../geometry/index.hpp"

SpinorTester::SpinorTester(std::string kernelName, const ParameterCollection parameterCollection, const SpinorTestParameters & testParameters, const ReferenceValues rV):
		KernelTester(kernelName, parameterCollection.hardwareParameters, parameterCollection.kernelParameters, testParameters, rV)
{
	code = device->getSpinorCode();
	doubleBuffer = new hardware::buffers::Plain<double> (1, device);
}

void fill_with_zero(spinor * in, int size);
void fill_with_one(spinor * in, int size);
void fill_with_zero_one(spinor * in, int size);
void fill_with_one_zero(spinor * in, int size);
void fill_with_ascending(spinor * in, int size);
void fillWithAscendingComplex(spinor * in, int size);
void fill_with_one_eo(spinor * in, const int, const bool, const int ns, const int nt);

spinor * SpinorfieldCreator::createSpinorfield(SpinorFillType fillTypeIn)
{
  spinor * in;
  in = new spinor[numberOfElements];
  switch (fillTypeIn) {
	case SpinorFillType::zero :
		fill_with_zero(in, numberOfElements);
		break;
	case SpinorFillType::one :
		fill_with_one(in, numberOfElements);
		break;
	case SpinorFillType::zeroOne :
		fill_with_zero_one(in, numberOfElements);
		break;
	case SpinorFillType::oneZero :
		fill_with_one_zero(in, numberOfElements);
		break;
	case SpinorFillType::ascendingReal :
		fill_with_ascending(in, numberOfElements);
		break;
	case SpinorFillType::ascendingComplex :
		fillWithAscendingComplex(in, numberOfElements);
		break;
	default:
		logger.fatal() << "do not know fill type!";
  }
  BOOST_REQUIRE(in);
  return in;
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

int calculateSpinorfieldSize(const LatticeExtents latticeExtendsIn) noexcept
{
	return 	latticeExtendsIn.getLatticeVolume();
}

int calculateEvenOddSpinorfieldSize(const LatticeExtents latticeExtendsIn) noexcept
{
	return 	calculateSpinorfieldSize(latticeExtendsIn) / 2;
}

//@todo: should this take the same fillType for both spinors?
void EvenOddSpinorfieldCreator::fillTwoSpinorBuffers(const hardware::buffers::Spinor * in1, const SpinorFillType fillTypeIn1, const hardware::buffers::Spinor * in2, const SpinorFillType fillTypeIn2)
{
	spinor * sf_in1;
	spinor * sf_in2;
	sf_in1 = new spinor[numberOfElements];
	sf_in2 = new spinor[numberOfElements];
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);

	in1->load(SpinorfieldCreator::createSpinorfield(fillTypeIn1));
	in2->load(SpinorfieldCreator::createSpinorfield(fillTypeIn2));

	delete sf_in1;
	delete sf_in2;
}

void EvenOddSpinorfieldCreator::fillTwoSpinorBuffersDependingOnParity(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2)
{
	spinor * sf_in1;
	spinor * sf_in2;
	sf_in1 = new spinor[numberOfElements];
	sf_in2 = new spinor[numberOfElements];
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);

	fillTwoSpinorfieldsDependingOnParity(sf_in1, sf_in2, numberOfElements);

	in1->load(sf_in1);
	in2->load(sf_in2);

	delete sf_in1;
	delete sf_in2;
}

static spinor fillSpinorWithNumber(hmc_complex content)
{
	spinor in;
	in.e0.e0 = content;
	in.e0.e1 = content;
	in.e0.e2 = content;
	in.e1.e0 = content;
	in.e1.e1 = content;
	in.e1.e2 = content;
	in.e2.e0 = content;
	in.e2.e1 = content;
	in.e2.e2 = content;
	in.e3.e0 = content;
	in.e3.e1 = content;
	in.e3.e2 = content;
	return in;
}

void EvenOddSpinorfieldCreator::fillTwoSpinorfieldsDependingOnParity(spinor * in1, spinor * in2, int size)
{
		int x, y, z, t;
		bool parityOfSite;
		int global_pos;
		int ns, nt;

		ns = latticeExtents.getNs();
		nt = latticeExtents.getNt();

		for (x = 0; x < ns; x++) {
			for (y = 0; y < ns; y++) {
				for (z = 0; z < ns; z++) {
					for (t = 0; t < nt; t++) {
						global_pos = uint(Index(x,y,z,t,latticeExtents));
						if (global_pos >= size)
							break;

						parityOfSite = (x + y + z + t) % 2 == 0;
						if (parityOfSite)
						{
							in1[global_pos] =fillSpinorWithNumber( hmc_complex_one );
							in2[global_pos] =fillSpinorWithNumber( hmc_complex_zero );
						}
						else
						{
							in1[global_pos] =fillSpinorWithNumber( hmc_complex_zero );
							in2[global_pos] =fillSpinorWithNumber( hmc_complex_one );
						}

					}
				}
			}
		}
		return;
	}

spinor * NonEvenOddSpinorfieldCreator::createSpinorfieldWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites)
{
  spinor * in;
  in = new spinor[numberOfElements];
  fill_with_one_eo(in, numberOfElements, fillEvenSites, latticeExtents.getNs(), latticeExtents.getNt());
  return in;
}



void fill_with_zero(spinor * in, int size)
{
  for(int i = 0; i < size; ++i) {
    in[i].e0.e0 = {0., 0.};
    in[i].e0.e1 = {0., 0.};
    in[i].e0.e2 = {0., 0.};
    in[i].e1.e0 = {0., 0.};
    in[i].e1.e1 = {0., 0.};
    in[i].e1.e2 = {0., 0.};
    in[i].e2.e0 = {0., 0.};
    in[i].e2.e1 = {0., 0.};
    in[i].e2.e2 = {0., 0.};
    in[i].e3.e0 = {0., 0.};
    in[i].e3.e1 = {0., 0.};
    in[i].e3.e2 = {0., 0.};
  }
  return;
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

void fill_with_zero_one(spinor * in, int size)
{
	for(int i= 0; i < size; i++) {
	    in[i].e0.e0 = {0., 0.};
	    in[i].e0.e1 = {0., 0.};
	    in[i].e0.e2 = {0., 0.};
	    in[i].e1.e0 = hmc_complex_one;
	    in[i].e1.e1 = hmc_complex_one;
	    in[i].e1.e2 = hmc_complex_one;
	    in[i].e2.e0 = {0., 0.};
	    in[i].e2.e1 = {0., 0.};
	    in[i].e2.e2 = {0., 0.};
	    in[i].e3.e0 = hmc_complex_one;
	    in[i].e3.e1 = hmc_complex_one;
	    in[i].e3.e2 = hmc_complex_one;
	}
}

void fill_with_one_zero(spinor * in, int size)
{
	for(int i= 0; i < size; i++) {
	    in[i].e0.e0 = hmc_complex_one;
	    in[i].e0.e1 = hmc_complex_one;
	    in[i].e0.e2 = hmc_complex_one;
	    in[i].e1.e0 = {0., 0.};
	    in[i].e1.e1 = {0., 0.};
	    in[i].e1.e2 = {0., 0.};
	    in[i].e2.e0 = hmc_complex_one;
	    in[i].e2.e1 = hmc_complex_one;
	    in[i].e2.e2 = hmc_complex_one;
	    in[i].e3.e0 = {0., 0.};
	    in[i].e3.e1 = {0., 0.};
	    in[i].e3.e2 = {0., 0.};
	}
}

void fill_with_ascending(spinor * in, int size)
{
	for(int i = 0; i < size; i++) {
	    in[i].e0.e0 = {1., 0.};
	    in[i].e0.e1 = {2., 0.};
	    in[i].e0.e2 = {3., 0.};
	    in[i].e1.e0 = {4., 0.};
	    in[i].e1.e1 = {5., 0.};
	    in[i].e1.e2 = {6., 0.};
	    in[i].e2.e0 = {7., 0.};
	    in[i].e2.e1 = {8., 0.};
	    in[i].e2.e2 = {9., 0.};
	    in[i].e3.e0 = {10., 0.};
	    in[i].e3.e1 = {11., 0.};
	    in[i].e3.e2 = {12., 0.};
	}
}

void fillWithAscendingComplex(spinor * in, int size)
{
	for(int i = 0; i < size; i++) {
	    in[i].e0.e0 = {1., 2.};
	    in[i].e0.e1 = {3., 4.};
	    in[i].e0.e2 = {5., 6.};
	    in[i].e1.e0 = {7., 8.};
	    in[i].e1.e1 = {9., 10.};
	    in[i].e1.e2 = {11., 12.};
	    in[i].e2.e0 = {13., 14.};
	    in[i].e2.e1 = {15., 16.};
	    in[i].e2.e2 = {17., 18.};
	    in[i].e3.e0 = {19., 20.};
	    in[i].e3.e1 = {21., 22.};
	    in[i].e3.e2 = {23., 24.};
	}
}

void fill_with_one_eo(spinor * in, const int size, const bool fillEvenSites, const int ns, const int nt)
{
	int x, y, z, t;
	hmc_complex content;
	bool parityOfSite;
	int global_pos;
	LatticeExtents lE(ns,nt);

	for (x = 0; x < ns; x++) {
		for (y = 0; y < ns; y++) {
			for (z = 0; z < ns; z++) {
				for (t = 0; t < nt; t++) {
					global_pos = uint(Index(x,y,z,t,lE));
					if (global_pos >= size)
						break;

					parityOfSite = (x + y + z + t) % 2 == 0;
					content = (parityOfSite) ?
						(fillEvenSites ? hmc_complex_one : hmc_complex_zero) :
						(fillEvenSites ? hmc_complex_zero : hmc_complex_one);

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
