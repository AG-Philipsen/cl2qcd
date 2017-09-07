/* @file
 * Implementation of the hardware::code::SpinorStaggeredTester class
 *
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

#include "SpinorStaggeredTester.hpp"
#include "../../geometry/index.hpp"

SpinorStaggeredTester::SpinorStaggeredTester(std::string kN, const ParameterCollection pC, const SpinorStaggeredTestParameters & tP, const ReferenceValues rV):
		KernelTester(kN, pC.hardwareParameters, pC.kernelParameters, tP, rV)
{
	code = device->getSpinorStaggeredCode();
	doubleBuffer = new hardware::buffers::Plain<double> (1,device);
}

static void fill_with_one(su3vec * sf_in, int size)
{
  for(int i = 0; i < size; ++i) {
    sf_in[i].e0 = hmc_complex_one;
    sf_in[i].e1 = hmc_complex_one;
    sf_in[i].e2 = hmc_complex_one;
  }
  return;
}

void fill_with_zero(su3vec * sf_in, int size)
{
  for(int i = 0; i < size; ++i) {
    sf_in[i].e0 = hmc_complex_zero;
    sf_in[i].e1 = hmc_complex_zero;
    sf_in[i].e2 = hmc_complex_zero;
  }
  return;
}

void fill_with_ascending(su3vec * sf_in, int size)
{
  for(int i = 0; i < size; ++i) {
    sf_in[i].e0 = {1., 0.};
    sf_in[i].e1 = {2., 0.};
    sf_in[i].e2 = {3., 0.};
  }
  return;
}

void fillWithAscendingComplex(su3vec * sf_in, int size)
{
  for(int i = 0; i < size; ++i) {
    sf_in[i].e0 = {1., 2.};
    sf_in[i].e1 = {3., 4.};
    sf_in[i].e2 = {5., 6.};
  }
  return;
}

void fill_with_one_eo(su3vec * sf_in, int size, bool eo, const int ns, const int nt)
{
  int x,y,z,t;
  hmc_complex content;
  bool parityOfSite;
  int global_pos;

  LatticeExtents lE(ns,nt);

  for (x = 0; x<ns; x++){
    for (y = 0; y<ns; y++){
      for (z = 0; z<ns; z++){
		for (t = 0; t<nt; t++){
		  global_pos  = uint(Index(x,y,z,t,lE));
		  if (global_pos >= size)
			break;
		  parityOfSite = (x + y + z + t) % 2 == 0;
		  content = (parityOfSite) ? (eo ? hmc_complex_one : hmc_complex_zero) :
						  (eo ? hmc_complex_zero : hmc_complex_one);
		  sf_in[global_pos].e0 = content;
		  sf_in[global_pos].e1 = content;
		  sf_in[global_pos].e2 = content;
		}
      }
    }
  }
  return;
}

su3vec * SpinorStaggeredfieldCreator::createSpinorfield(SpinorFillType fillTypeIn)
{
	su3vec * in;
  in = new su3vec[numberOfElements];
  switch (fillTypeIn) {
	case SpinorFillType::zero :
		fill_with_zero(in, numberOfElements);
		break;
	case SpinorFillType::one :
		fill_with_one(in, numberOfElements);
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

su3vec * NonEvenOddSpinorStaggeredfieldCreator::createSpinorfieldWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites)
{
  su3vec * in;
  in = new su3vec[numberOfElements];
  fill_with_one_eo(in, numberOfElements, fillEvenSites, latticeExtents.getNs(), latticeExtents.getNt());
  return in;
}

su3vec * EvenOddSpinorStaggeredfieldCreator::createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites)
{
  su3vec * in;
  in = new su3vec[numberOfElements];
  fill_with_one_eo(in, numberOfElements, fillEvenSites, latticeExtents.getNs(), latticeExtents.getNt());
  return in;
}

void SpinorStaggeredTester::calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<su3vec> * in)
{
  code->set_float_to_global_squarenorm_device(in, doubleBuffer);
  doubleBuffer->dump(&kernelResult[0]);
}

void SpinorStaggeredTester::calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::SU3vec * in)
{
  code->set_float_to_global_squarenorm_eoprec_device(in, doubleBuffer);
  doubleBuffer->dump(&kernelResult[0]);
}

//This function sums the real and imaginary parts of all su3vec contained in sf_in
hmc_float count_sf(su3vec * sf_in, int size)
{
  hmc_float sum = 0.;
  for (int i=0; i<size; i++){
    sum +=
        sf_in[i].e0.re + sf_in[i].e0.im
      + sf_in[i].e1.re + sf_in[i].e1.im
      + sf_in[i].e2.re + sf_in[i].e2.im;
  }
  return sum;
}

hmc_float squareNorm(su3vec * sf_in, int size)
{
  hmc_float sum = 0.;
  for (int i=0; i<size; i++){
    sum +=
        sf_in[i].e0.re * sf_in[i].e0.re + sf_in[i].e0.im * sf_in[i].e0.im
      + sf_in[i].e1.re * sf_in[i].e1.re + sf_in[i].e1.im * sf_in[i].e1.im
      + sf_in[i].e2.re * sf_in[i].e2.re + sf_in[i].e2.im * sf_in[i].e2.im;
  }
  return sum;
}

hmc_float calc_var_sf(su3vec * sf_in, int size, hmc_float sum){
  hmc_float var = 0.;
  for(int k=0; k<size; k++){
    var +=
        calc_var(sf_in[k].e0.re, sum)
      + calc_var(sf_in[k].e0.im, sum)
      + calc_var(sf_in[k].e1.re, sum)
      + calc_var(sf_in[k].e1.im, sum)
      + calc_var(sf_in[k].e2.re, sum)
      + calc_var(sf_in[k].e2.im, sum);
  }
  return var;
}

/**
 * Function that returns a vector with the 6 real number contained in an su3vec (used in gaussian tests)
 */
std::vector<hmc_float> reals_from_su3vec(su3vec v){
  std::vector<hmc_float> out;
  out.push_back(v.e0.re);
  out.push_back(v.e0.im);
  out.push_back(v.e1.re);
  out.push_back(v.e1.im);
  out.push_back(v.e2.re);
  out.push_back(v.e2.im);
  return out;
}
