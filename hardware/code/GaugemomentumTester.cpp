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

#include "GaugemomentumTester.hpp"

GaugemomentumTester::GaugemomentumTester(std::string kernelName, std::string inputfile, int numberOfValues, int typeOfComparision) :
  KernelTester(kernelName, getSpecificInputfile(inputfile), numberOfValues, typeOfComparision)
{
  numberOfAlgebraElements = meta::get_vol4d(*parameters) * NDIM * meta::get_su3algebrasize();
  numberOfGaugemomentumElements = meta::get_vol4d(*parameters) * NDIM;
  useRandom = (parameters->get_solver() == meta::Inputparameters::cg)  ? false : true;
  
  code = device->get_gaugemomentum_code();
  doubleBuffer = new hardware::buffers::Plain<double> (1, device);
  gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(numberOfGaugemomentumElements, device);
}

GaugemomentumTester::~GaugemomentumTester()
{
  delete doubleBuffer;
  delete gaugemomentumBuffer;
}


std::string GaugemomentumTester::getSpecificInputfile(std::string inputfileIn)
{
  return "gaugemomentum/" + inputfileIn;
}

double * GaugemomentumTester::createGaugemomentum(int seed)
{
  double * gm_in;
  gm_in = new double[numberOfAlgebraElements];
  useRandom ? fill_with_random(gm_in, seed) : fill_with_one(gm_in);
  BOOST_REQUIRE(gm_in);
  return gm_in;    
}

double * GaugemomentumTester::createGaugemomentumBasedOnFilltype(Filltype filltype)
{
  double * gm_in;
  gm_in = new double[numberOfAlgebraElements];
  switch(filltype)
	{
		case one:
			fill_with_one(gm_in);
			break;
		case zero:
			fill_with_zero(gm_in);
			break;
	}
  BOOST_REQUIRE(gm_in);
  return gm_in;    
}

void GaugemomentumTester::fill_with_one(double * sf_in)
{
  for(int i = 0; i < (int) numberOfAlgebraElements; ++i) {
    sf_in[i] = 1.;
  }
  return;
}

void GaugemomentumTester::fill_with_zero(double * sf_in)
{
  for(int i = 0; i < (int) numberOfAlgebraElements; ++i) {
    sf_in[i] = 0.;
  }
  return;
}

void GaugemomentumTester::fill_with_random(double * sf_in, int seed)
{
  prng_init(seed);
  for(int i = 0; i < (int) numberOfAlgebraElements; ++i) {
    sf_in[i] = prng_double();
  }
  return;
}

void GaugemomentumTester::calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Gaugemomentum * in)
{
  code->set_float_to_gaugemomentum_squarenorm_device(in, doubleBuffer);
  doubleBuffer->dump(&kernelResult[0]);
}

double GaugemomentumTester::count_gm(ae * ae_in, int size)
{
  double sum = 0.;
  for (int i = 0; i<size;i++){
    sum +=
      ae_in[i].e0
      + ae_in[i].e1
      + ae_in[i].e2
      + ae_in[i].e3
      + ae_in[i].e4
      + ae_in[i].e5
      + ae_in[i].e6
      + ae_in[i].e7;
    }
  return sum;
}

double GaugemomentumTester::calc_var(double in, double mean){
  return (in - mean) * (in - mean);
}
  
double GaugemomentumTester::calc_var_gm(ae * ae_in, int size, double sum){
  double var = 0.;
  for(int k = 0; k<size; k++){
    var +=
      calc_var(   ae_in[k].e0 , sum) 
      + calc_var( ae_in[k].e1 , sum) 
      + calc_var( ae_in[k].e2 , sum)
      + calc_var( ae_in[k].e3 , sum) 
      + calc_var( ae_in[k].e4 , sum) 
      + calc_var( ae_in[k].e5 , sum) 
	+ calc_var( ae_in[k].e6 , sum) 
      + calc_var( ae_in[k].e7 , sum) ;
  }
  return var;
}
