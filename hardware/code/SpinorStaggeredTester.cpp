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
#include "../../host_functionality/host_geometry.h"
#include "spinors.hpp" //this is for get_spinorfieldsize, get_eoprec_spinorfieldsize


void SpinorStaggeredTester::setMembers(){
	//todo: some of these could also be put into the specific child-classes where they are actually used.
	inputfield = NULL;
	spinorfieldElements = parameters->get_nspace() * parameters->get_nspace() * parameters->get_nspace() * parameters->get_ntime(); //todo: make proper
	spinorfieldEvenOddElements = parameters->get_nspace() * parameters->get_nspace() * parameters->get_nspace() * parameters->get_ntime() / 2; //todo: make proper
	(parameters->get_solver() == common::cg) ? useRandom = false : useRandom =true;
	parameters->get_read_multiple_configs() ? evenOrOdd = true : evenOrOdd = false;
	parameters->get_read_multiple_configs() ? calcVariance=false : calcVariance = true;
	alpha_host = {parameters->get_beta(), parameters->get_rho()};
	beta_host = {parameters->get_kappa(), parameters->get_mu()};
	iterations = parameters->get_integrationsteps(0);
}

SpinorStaggeredTester::SpinorStaggeredTester(std::string kernelName, std::string inputfileIn,
					      int numberOfValues, int typeOfComparision) :
	prngParameters( parameters ), KernelTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues, typeOfComparision)
{
	code = device->getSpinorStaggeredCode();
	prng = new physics::PRNG(*system, &prngParameters);
	doubleBuffer = new hardware::buffers::Plain<double> (1, device);
	allocatedObjects = true;
	setMembers();
}

SpinorStaggeredTester::SpinorStaggeredTester(meta::Inputparameters * parameters, const hardware::System * system,
					     hardware::Device * device) : 
      KernelTester(parameters, system, device), allocatedObjects(false), prngParameters( parameters )
{
	setMembers();
	code = device->getSpinorStaggeredCode();
}

SpinorStaggeredTester::~SpinorStaggeredTester()
{
	if(allocatedObjects){
		delete doubleBuffer;
 		delete prng;
	}
	doubleBuffer = NULL;
	prng = NULL;
	code = NULL;
	if(inputfield != NULL)
	  delete[] inputfield;
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

void SpinorStaggeredTester::fill_with_zero(su3vec * sf_in, int size)
{
  for(int i = 0; i < size; ++i) {
    sf_in[i].e0 = hmc_complex_zero;
    sf_in[i].e1 = hmc_complex_zero;
    sf_in[i].e2 = hmc_complex_zero;
  }
  return;
}

static void fill_with_random(su3vec * sf_in, int size, int seed)
{
  prng_init(seed);
  for(int i = 0; i < size; ++i) {
    sf_in[i].e0.re = prng_double();
    sf_in[i].e1.re = prng_double();
    sf_in[i].e2.re = prng_double();
    
    sf_in[i].e0.im = prng_double();
    sf_in[i].e1.im = prng_double();
    sf_in[i].e2.im = prng_double();
  } 
  return;
}

//This function sums the real and imaginary parts of all su3vec contained in sf_in
hmc_float SpinorStaggeredTester::count_sf(su3vec * sf_in, int size)
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

//The following two function return the sum of the square deviation (frome the mean) of the numbers
//in sf_in. To get the variance, i.e. the mean square deviation, the result
//must be divided by the number of numbers summed.
hmc_float SpinorStaggeredTester::calc_var(hmc_float in, hmc_float mean){
  return (in - mean) * (in - mean);
}

hmc_float SpinorStaggeredTester::calc_var_sf(su3vec * sf_in, int size, hmc_float sum){
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

//This function fills the field sf_in in the following way
// eo==true  ---> sf_in[even]=ONE  and sf_in[odd]=ZERO
// eo==false ---> sf_in[even]=ZERO and sf_in[odd]=ONE
static void fill_with_one_eo(su3vec * sf_in, int size, bool eo, const meta::Inputparameters * parameters)
{
  int ns = parameters->get_nspace();
  int nt = parameters->get_ntime();
  int x,y,z,t;
  hmc_complex content;
  int coord[4];
  bool parityOfSite;
  int nspace;
  int global_pos;
  
  for (x = 0; x<ns; x++){
    for (y = 0; y<ns; y++){
      for (z = 0; z<ns; z++){
	for (t = 0; t<nt; t++){
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  nspace =  get_nspace(coord, *parameters);
	  global_pos = get_global_pos(nspace, t, *parameters);
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

//This function works in the following way
// eo==true  ---> sum of all components of sf_in[even] is returned
// eo==false ---> sum of all components of sf_in[odd]  is returned
hmc_float SpinorStaggeredTester::count_sf_eo(su3vec * sf_in, int size, bool eo)
{
  int ns = parameters->get_nspace();
  int nt = parameters->get_ntime();
  int x,y,z,t;
  hmc_float sum = 0.;
  for (x = 0; x<ns; x++){
    for (y = 0; y<ns; y++){
      for (z = 0; z<ns; z++){
	for (t = 0; t<nt; t++){
	  int coord[4];
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  int nspace =  get_nspace(coord, *parameters);
	  int global_pos = get_global_pos(nspace, t, *parameters);
	  if (global_pos >= size)
	    break;
	  if (
	      ( eo ==true && (x+y+z+t) %2 == 0) ||
	      ( eo ==false &&  (x+y+z+t) %2 == 1 )
	      )
	    {
	      int i = global_pos;
	      sum +=
		sf_in[i].e0.re+ sf_in[i].e0.im 
		+sf_in[i].e1.re+ sf_in[i].e1.im 
		+sf_in[i].e2.re+ sf_in[i].e2.im;
	    }
	  else{
	    continue;
	  }
	}
      }
    }
  }
  return sum;
}


su3vec * SpinorStaggeredTester::createSpinorfield(size_t numberOfElements, int seed)
{
  if(inputfield !=NULL)
    delete[] inputfield;
  inputfield = new su3vec[numberOfElements];
  useRandom ? fill_with_random(inputfield, numberOfElements, seed) : fill_with_one(inputfield, numberOfElements);
  BOOST_REQUIRE(inputfield);
  return inputfield;
}

su3vec * SpinorStaggeredTester::createSpinorfieldWithOnesAndZerosDependingOnSiteParity()
{
  if(inputfield !=NULL)
    delete[] inputfield;
  inputfield = new su3vec[spinorfieldElements];
  fill_with_one_eo(inputfield, spinorfieldElements, evenOrOdd, parameters);
  return inputfield;
}

su3vec * SpinorStaggeredTester::createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity()
{
  if(inputfield !=NULL)
    delete[] inputfield;
  inputfield = new su3vec[spinorfieldEvenOddElements];
  fill_with_one_eo(inputfield, spinorfieldEvenOddElements, evenOrOdd, parameters);
  return inputfield;
}

std::string SpinorStaggeredTester::getSpecificInputfile(std::string inputfileIn)
{
  return "spinors_staggered/" + inputfileIn;
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

/////////////////////////////////////////////////////////////////////////

/**
 * Fuction that "convert" a su3vec to a string with a proper structure to be   
 * written to the text file that will be later used for the reference code     
 */
static std::string su3vec_to_string(su3vec m)
{
  std::ostringstream os;
  os.precision(16);
  os << "(" << m.e0.re << "," << m.e0.im << ") (" << m.e1.re << "," << m.e1.im << ") (" << m.e2.re << "," << m.e2.im << ")\n\n";
  return os.str();
}

//Tool to be used in the function print_gaugefield_to_textfile or print_staggeredfield_to_textfile
static void get_full_coord_from_site_idx(int site_idx, int &x, int &y, int &z, int &t, const int ns)
{
  int volspace=ns*ns*ns;
  int space=site_idx%volspace;
  t=site_idx/volspace;
  z=space/ns/ns;
  int acc=z;
  y=space/ns-ns*acc;
  acc=ns*acc+y;
  x=space-ns*acc;
}

/**
 *  In the reference code the lattice is reorganized in the following way: 
 * 
 *   0             size
 *   |-------|-------|   
 *       e       o 
 * 
 *  where e=even, o=odd, whereas size=VOL4D. 
 *  Hence, in order to use the same random staggered field in tests  
 *  I have to print it to a text file according this scheme.
 */
void SpinorStaggeredTester::print_staggeredfield_to_textfile(std::string outputfile, su3vec * sf)
{
  int nt=parameters->get_ntime();
  int ns=parameters->get_nspace();
  if(ns!=nt){
    logger.fatal() << "The lattice must be isotropic to call the function print_staggeredfield_to_textfile(...)!";
    abort();
  }
  //sf     is the su3vec array ordered with the "superindex scheme"                                                                     
  //sf_new is the su3vec array in the right order (ref. code scheme) to be written to the file                                          
  su3vec *sf_new = new su3vec[ns*ns*ns*nt];
  //Now I have conf_old and I have to fill properly conf_new                                                                            
  int x,y,z,t,num,even,size;
  size=ns*ns*ns*nt;
  for(int i=0; i<ns*ns*ns*nt; i++){
    get_full_coord_from_site_idx(i,x,y,z,t,ns);
    even = (x+y+z+t)%2;
    // even=0 for even sites                                                                                                            
    // even=1 for odd sites                                                                                                             
    num = even*size/2 + (x+y*ns+z*ns*ns+t*ns*ns*ns)/2;
    // num is where, in conf_new, conf_old[...] is to be written                                                                        
    sf_new[num]=sf[i];
  }
  //Now we can write sf_new to the file 
  std::ofstream file(outputfile.c_str());
  file << ns << " " << ns << " " << ns << " " << nt << std::endl;
  for(int i=0; i<ns*ns*ns*nt; i++){
    get_full_coord_from_site_idx(i,x,y,z,t,ns);
    file << su3vec_to_string(sf_new[i]);
  }
  file.close();
}

void SpinorStaggeredTester::print_staggeredfield_eo_to_textfile(std::string outputfile, su3vec * sf)
{
  int nt=parameters->get_ntime();
  int ns=parameters->get_nspace();
  if(ns!=nt){
    logger.fatal() << "The lattice must be isotropic to call the function print_staggeredfield_to_textfile(...)!";
    abort();
  }
  //sf     is the su3vec array ordered with the "even-odd superindex scheme"                                                                     
  //sf_new is the su3vec array in the right order (ref. code scheme) to be written to the file
  // ======> hence sf_new is in this case equal to sf that contain the values of the field only in
  //         even (or odd) sites
  //We can write sf directly to the file 
  std::ofstream file(outputfile.c_str());
  file << ns << " " << ns << " " << ns << " " << nt << std::endl;
  for(int i=0; i<ns*ns*ns*nt/2; i++)
    file << su3vec_to_string(sf[i]);
  file.close();
}

/**
 * Function that returns a vector with the 6 real number contained in an su3vec (used in gaussian tests)
 */
std::vector<hmc_float> SpinorStaggeredTester::reals_from_su3vec(su3vec v){
  std::vector<hmc_float> out;
  out.push_back(v.e0.re);
  out.push_back(v.e0.im);
  out.push_back(v.e1.re);
  out.push_back(v.e1.im);
  out.push_back(v.e2.re);
  out.push_back(v.e2.im);
  return out;
}


