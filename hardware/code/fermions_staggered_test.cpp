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

// #include "testUtilities.hpp"
// 
// #include "../meta/util.hpp"
// #include "../host_functionality/host_random.h"
// #include "../hardware/code/spinors_staggered.hpp"
// #include "../hardware/code/spinors.hpp"
// 
// // use the boost test framework
// #define BOOST_TEST_DYN_LINK
// #define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_STAGGERED
// #include <boost/test/unit_test.hpp>
// 
// //some functionality
// #include "test_util.h"
// #include "test_util_staggered.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_STAGGERED

#include "SpinorStaggeredTester.hpp"
#include "fermions_staggered.hpp"
#include "../../physics/lattices/gaugefield.hpp"

class FermionStaggeredTester : public SpinorStaggeredTester{
   public:
	FermionStaggeredTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
	SpinorStaggeredTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues){
		code = device->get_fermion_staggered_code();
		gaugefield = new physics::lattices::Gaugefield(*system, *prng);
		/*
		print_gaugefield_to_textfile("ref_conf");
		logger.info() << "Produced the ref_conf text file with the links for the Ref.Code.";
		*/
	}
	
	virtual ~FermionStaggeredTester(){
		delete gaugefield;
		code = NULL;
	}
	
   protected:
	const hardware::code::Fermions_staggered * code;
	physics::lattices::Gaugefield * gaugefield;
	
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		//todo: this is ugly, find a better solution.
		// The problem is that the parent class calls a similar fct.
		return "../fermionsStaggered/" + inputfileIn;
	}
	
	const hardware::buffers::SU3* getGaugefieldBuffer() {
		return gaugefield->get_buffers()[0];
	}
	
	//This method is used to produce file for the Reference Code (D'Elia et al) -> implemented at the end
	void print_gaugefield_to_textfile(std::string outputfile);
};

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE(BUILD)

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
	    BOOST_CHECK_NO_THROW(FermionStaggeredTester("build", "fermions_staggered_build_input_1"));
	}
	
	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
	    BOOST_CHECK_NO_THROW(FermionStaggeredTester("build", "fermions_staggered_build_input_2"));
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

class FermionmatrixStaggeredTester : public FermionStaggeredTester{
   public:
	FermionmatrixStaggeredTester(std::string kernelName, std::string inputfile) : 
	    FermionStaggeredTester(kernelName, inputfile){
		in = new const hardware::buffers::Plain<su3vec>(spinorfieldElements, device);
		out = new const hardware::buffers::Plain<su3vec>(spinorfieldElements, device);
		in->load(createSpinorfield(spinorfieldElements));
		out->load(createSpinorfield(spinorfieldElements));
		/*
		print_staggeredfield_to_textfile("ref_vec_Mmatrix", createSpinorfield(spinorfieldElements));
		logger.info() << "Produced the ref_vec_Mmatrix file with the staggered field for the Ref.Code.";
		*/
	}
	virtual ~FermionmatrixStaggeredTester(){
		calcSquarenormAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Plain<su3vec> * in;
	const hardware::buffers::Plain<su3vec> * out;
};

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE( M_MATRIX )

	class MTester : public FermionmatrixStaggeredTester{
	   public:
		MTester(std::string inputfile) : FermionmatrixStaggeredTester("M_staggered", inputfile){
 			code->M_staggered_device(in, out,  getGaugefieldBuffer()); //mass is the def.value
		}
	};

	BOOST_AUTO_TEST_CASE( M_MATRIX_1)
	{
	    MTester("m_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_2)
	{
	    MTester("m_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_3)
	{
	    MTester("m_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_4)
	{
	    MTester("m_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_5)
	{
	    MTester("m_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_6)
	{
	    MTester("m_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_7)
	{
	    MTester("m_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_8)
	{
	    MTester("m_input_8");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_9)
	{
	    MTester("m_input_9");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_10)
	{
	    MTester("m_input_10");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_11)
	{
	    MTester("m_input_11");
	}
	
	BOOST_AUTO_TEST_CASE( M_MATRIX_12)
	{
	    MTester("m_input_12");
	}

BOOST_AUTO_TEST_SUITE_END()

///////////////////////////////////////

class FermionmatrixStaggeredEvenOddTester : public FermionStaggeredTester{
   public:
	FermionmatrixStaggeredEvenOddTester(std::string kernelName, std::string inputfile) :
	   FermionStaggeredTester(kernelName, inputfile){
		in = new const hardware::buffers::SU3vec(spinorfieldEvenOddElements, device);
		out = new const hardware::buffers::SU3vec(spinorfieldEvenOddElements, device);
		in->load(createSpinorfield(spinorfieldEvenOddElements));
		out->load(createSpinorfield(spinorfieldEvenOddElements));
	}
	~FermionmatrixStaggeredEvenOddTester(){
		calcSquarenormEvenOddAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
   protected:
	const hardware::buffers::SU3vec * in;
	const hardware::buffers::SU3vec * out;
};

///////////////////////////////////////

BOOST_AUTO_TEST_SUITE( DKS_EO )

	class DksEvenOddTester : public FermionmatrixStaggeredEvenOddTester{
	   public:
		DksEvenOddTester(std::string inputfile) : FermionmatrixStaggeredEvenOddTester("DKS_eo", inputfile){
 			evenOrOdd ? code->D_KS_eo_device(in, out,  getGaugefieldBuffer(), EVEN)
				  : code->D_KS_eo_device(in, out,  getGaugefieldBuffer(), ODD);
		}
	};

	BOOST_AUTO_TEST_CASE( DKS_EO_1)
	{
	    DksEvenOddTester("dks_input_1");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_2)
	{
	    DksEvenOddTester("dks_input_2");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_3)
	{
	    DksEvenOddTester("dks_input_3");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_4)
	{
	    DksEvenOddTester("dks_input_4");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_5)
	{
	    DksEvenOddTester("dks_input_5");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_6)
	{
	    DksEvenOddTester("dks_input_6");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_7)
	{
	    DksEvenOddTester("dks_input_7");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_8)
	{
	    DksEvenOddTester("dks_input_8");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_9)
	{
	    DksEvenOddTester("dks_input_9");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_10)
	{
	    DksEvenOddTester("dks_input_10");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_11)
	{
	    DksEvenOddTester("dks_input_11");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_12)
	{
	    DksEvenOddTester("dks_input_12");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_13)
	{
	    DksEvenOddTester("dks_input_13");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_14)
	{
	    DksEvenOddTester("dks_input_14");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_15)
	{
	    DksEvenOddTester("dks_input_15");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_16)
	{
	    DksEvenOddTester("dks_input_16");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_17)
	{
	    DksEvenOddTester("dks_input_17");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_18)
	{
	    DksEvenOddTester("dks_input_18");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_19)
	{
	    DksEvenOddTester("dks_input_19");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_20)
	{
	    DksEvenOddTester("dks_input_20");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_21)
	{
	    DksEvenOddTester("dks_input_21");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_22)
	{
	    DksEvenOddTester("dks_input_22");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_23)
	{
	    DksEvenOddTester("dks_input_23");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_24)
	{
	    DksEvenOddTester("dks_input_24");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_25)
	{
	    DksEvenOddTester("dks_input_25");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_26)
	{
	    DksEvenOddTester("dks_input_26");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_27)
	{
	    DksEvenOddTester("dks_input_27");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_28)
	{
	    DksEvenOddTester("dks_input_28");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_29)
	{
	    DksEvenOddTester("dks_input_29");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_30)
	{
	    DksEvenOddTester("dks_input_30");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_31)
	{
	    DksEvenOddTester("dks_input_31");
	}
	
	BOOST_AUTO_TEST_CASE( DKS_EO_32)
	{
	    DksEvenOddTester("dks_input_32");
	}
	
BOOST_AUTO_TEST_SUITE_END()

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Tool to be used in the function print_gaugefield_to_textfile 
/**
 * Fuction that "convert" a matrix to a string with a proper structure to be 
 * written to the text file that will be later used for the reference code 
 */
static std::string matrix_to_string(Matrixsu3 m)
{
  std::ostringstream os;
  os.precision(16);
  os << "(" << m.e00.re << "," << m.e00.im << ") (" << m.e01.re << "," << m.e01.im << ") (" << m.e02.re << "," << m.e02.im << ")\n";
  os << "(" << m.e10.re << "," << m.e10.im << ") (" << m.e11.re << "," << m.e11.im << ") (" << m.e12.re << "," << m.e12.im << ")\n";
  os << "(" << m.e20.re << "," << m.e20.im << ") (" << m.e21.re << "," << m.e21.im << ") (" << m.e22.re << "," << m.e22.im << ")\n\n";
  return os.str();
}

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

static inline Matrixsu3 multiply_matrixsu3_by_complex (Matrixsu3 in, hmc_complex factor)
{
  Matrixsu3 out;
  out.e00 = complexmult(in.e00, factor);
  out.e01 = complexmult(in.e01, factor);
  out.e02 = complexmult(in.e02, factor);
  out.e10 = complexmult(in.e10, factor);
  out.e11 = complexmult(in.e11, factor);
  out.e12 = complexmult(in.e12, factor);
  out.e20 = complexmult(in.e20, factor);
  out.e21 = complexmult(in.e21, factor);
  out.e22 = complexmult(in.e22, factor);
  return out;
}

/**
 *  In the reference code the lattice is reorganized in the following way:    
 *
 *  links used according to this scheme
 *   0             size            size2           size3         no_links  
 *   |-------|-------|-------|-------|-------|-------|-------|-------|
 *      e        o       e       o       e       o       e       o  
 *        x-dir         y-dir         z-dir         t-dir  
 *
 *  where e=even, o=odd, whereas size=VOL4D.  
 *  Hence, in order to use the same random configuration in tests 
 *  I have to print all links to a text file according this scheme.  
 *
 *  @note: In our program mu=0 is the TIME direction and mu=1,2,3 are the x,y,z direction!!!   
 * 
 */
void FermionStaggeredTester::print_gaugefield_to_textfile(std::string outputfile)
{
  int nt=parameters->get_ntime();
  int ns=parameters->get_nspace();
  if(ns!=nt){
    logger.fatal() << "The lattice must be isotropic to call the function print_gaugefield_to_textfile(...)!";
    abort();
  }
  //conf_old is the Matrixsu3 array with the links in the standard order (standard for this code)                                     
  //conf_new is the Matrixsu3 array in the right order (ref. code scheme) to be written to the file                                   
  Matrixsu3 *conf_old=new Matrixsu3[ns*ns*ns*nt*4];
  Matrixsu3 *conf_new=new Matrixsu3[ns*ns*ns*nt*4];
  device->get_gaugefield_code()->exportGaugefield(conf_old, gaugefield->get_buffers()[0]);
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
    conf_new[num       ]=conf_old[4*i+1]; //x-dir                                                                                     
    conf_new[num+size  ]=conf_old[4*i+2]; //y-dir                                                                                     
    conf_new[num+size*2]=conf_old[4*i+3]; //z-dir                                                                                     
    conf_new[num+size*3]=conf_old[4*i  ]; //t-dir      
  }
  //Now we can write conf_new to the file                                                                                             
  std::ofstream file(outputfile.c_str());
  file << ns << " " << ns << " " << ns << " " << nt << " ";
  file << parameters->get_beta() << " " << parameters->get_mass() << " 12345" << std::endl;
  //The last number that I set to 12345 should be the hmc iteration; here it is not relevant                                          
  for(int i=0; i<ns*ns*ns*nt*4; i++)
    file << matrix_to_string(conf_new[i]);
  file.close();
}




