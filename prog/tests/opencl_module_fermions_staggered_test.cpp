#include "../meta/util.hpp"
#include "../host_random.h"
#include "../physics/lattices/gaugefield.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/spinors_staggered.hpp"
#include "../hardware/code/fermions_staggered.hpp"
#include </usr/include/c++/4.7/fstream>

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_STAGGERED
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

class TestGaugefield {

public:
	TestGaugefield(const hardware::System * system) : system(system), prng(*system), gf(*system, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		auto inputfile = system->get_inputparameters();
		meta::print_info_hmc("test program", inputfile);
	};

	const hardware::code::Fermions_staggered * get_device();
	const hardware::buffers::SU3 * get_gaugefield();
	const hardware::code::Gaugefield* get_gf_code();
        void save_conf();

private:
	const hardware::System * const system;
	physics::PRNG prng;
        physics::lattices::Gaugefield gf; //I changed this variable from const to not const to be able to save the conf to a lime file!
};

void fill_sf_with_one(su3vec * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0 = hmc_complex_one;
		sf_in[i].e1 = hmc_complex_one;
		sf_in[i].e2 = hmc_complex_one;
	}
	return;
}

void fill_sf_with_random(su3vec * sf_in, int size, int seed)
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

void fill_sf_with_random(su3vec * sf_in, int size)
{
	fill_sf_with_random(sf_in, size, 123456);
}

const hardware::code::Fermions_staggered* TestGaugefield::get_device()
{
	return system->get_devices()[0]->get_fermion_staggered_code();
}

const hardware::code::Gaugefield* TestGaugefield::get_gf_code()
{
	return system->get_devices()[0]->get_gaugefield_code();
}

const hardware::buffers::SU3 * TestGaugefield::get_gaugefield()
{
	return gf.get_buffers().at(0);
}

void TestGaugefield::save_conf()
{
        gf.save("conf_lime",13);
}  

/**
 * Fuction that "convert" a matrix to a string with a proper structure to be
 * written to the text file that will be later used for the reference code
 */
std::string matrix_to_string(Matrixsu3 m)
{
       std::ostringstream os;
       os.precision(16);
       os << "(" << m.e00.re << "," << m.e00.im << ") (" << m.e01.re << "," << m.e01.im << ") (" << m.e02.re << "," << m.e02.im << ")\n";
       os << "(" << m.e10.re << "," << m.e10.im << ") (" << m.e11.re << "," << m.e11.im << ") (" << m.e12.re << "," << m.e12.im << ")\n";
       os << "(" << m.e20.re << "," << m.e20.im << ") (" << m.e21.re << "," << m.e21.im << ") (" << m.e22.re << "," << m.e22.im << ")\n\n";
       return os.str();
}

/**
 * Fuction that "convert" a su3vec to a string with a proper structure to be
 * written to the text file that will be later used for the reference code
 */
std::string su3vec_to_string(su3vec m)
{
       std::ostringstream os;
       os.precision(16);
       os << "(" << m.e0.re << "," << m.e0.im << ") (" << m.e1.re << "," << m.e1.im << ") (" << m.e2.re << "," << m.e2.im << ")\n\n";
       return os.str();
}

//Tool to be used in the function print_gaugefield_to_textfile
void get_full_coord_from_site_idx(int site_idx, int &x, int &y, int &z, int &t, const int ns)
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

//Tool to be used in the function print_gaugefield_to_textfile
void copy_Matrixsu3(Matrixsu3 &a, const Matrixsu3 b)
{
  //a=b
  a.e00.re=b.e00.re;
  a.e00.im=b.e00.im;
  a.e01.re=b.e01.re;
  a.e01.im=b.e01.im;
  a.e02.re=b.e02.re;
  a.e02.im=b.e02.im;
  a.e10.re=b.e10.re;
  a.e10.im=b.e10.im;
  a.e11.re=b.e11.re;
  a.e11.im=b.e11.im;
  a.e12.re=b.e12.re;
  a.e12.im=b.e12.im;
  a.e20.re=b.e20.re;
  a.e20.im=b.e20.im;
  a.e21.re=b.e21.re;
  a.e21.im=b.e21.im;
  a.e22.re=b.e22.re;
  a.e22.im=b.e22.im;
}

void copy_su3vec(su3vec &a, const su3vec b)
{
  //a=b
  a.e0.re=b.e0.re;
  a.e0.im=b.e0.im;
  a.e1.re=b.e1.re;
  a.e1.im=b.e1.im;
  a.e2.re=b.e2.re;
  a.e2.im=b.e2.im;
}

/**
 *  In the reference code the lattice is reorganized in the following way:
 * 
 *  links used according to this scheme
 *   0            size         size2         size3         no_links
 *   |------|------|------|------|------|------|------|------|
 *      e      o      e       o     e      o      e      o
 *        x-dir         y-dir         z-dir         t-dir
 * 
 *  where e=even, o=odd, whereas size=VOL4D.
 *  Hence, in order to use the same random configuration in tests
 *  I have to print all links to a text file according this scheme. 
 * 
 *  @note: In our program mu=0 is the TIME direction and mu=1,2,3 are the x,y,z direction!!! 
 * 
 */
void print_gaugefield_to_textfile(std::string outputfile, TestGaugefield * cpu, meta::Inputparameters params)
{
	int nt=params.get_ntime();
	int ns=params.get_nspace();
	if(ns!=nt){
	  logger.fatal() << "The lattice must be isotropic to call the function print_gaugefield_to_textfile(...)!";
	  abort();
	}
	//conf_old is the Matrixsu3 array with the links in the standard order (standard for this code)
	//conf_new is the Matrixsu3 array in the right order (ref. code scheme) to be written to the file
	Matrixsu3 *conf_old=new Matrixsu3[ns*ns*ns*nt*4];
	Matrixsu3 *conf_new=new Matrixsu3[ns*ns*ns*nt*4];
	cpu->get_gf_code()->exportGaugefield(conf_old,cpu->get_gaugefield());
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
	  copy_Matrixsu3(conf_new[num       ],conf_old[4*i+1]); //x-dir
	  copy_Matrixsu3(conf_new[num+size  ],conf_old[4*i+2]); //y-dir
	  copy_Matrixsu3(conf_new[num+size*2],conf_old[4*i+3]); //z-dir
	  copy_Matrixsu3(conf_new[num+size*3],conf_old[4*i  ]); //t-dir
	  //Only to check if the function works correctly ---> Comment out these line if all is fine
	  /*
	  conf_new[num].e00={x,0.};
	  conf_new[num].e01={y,0.};
	  conf_new[num].e02={z,0.};
	  conf_new[num].e10={t,0.};
	  conf_new[num].e11={0.,0.};
	  conf_new[num].e12={0.,0.};
	  conf_new[num].e20={0.,0.};
	  conf_new[num].e21={num,0.};
	  conf_new[num].e22={even,0.};
	  conf_new[num+size].e00={x,0.};
	  conf_new[num+size].e01={y,0.};
	  conf_new[num+size].e02={z,0.};
	  conf_new[num+size].e10={t,0.};
	  conf_new[num+size].e11={0.,0.};
	  conf_new[num+size].e12={0.,0.};
	  conf_new[num+size].e20={0.,0.};
	  conf_new[num+size].e21={num,0.};
	  conf_new[num+size].e22={even,0.};
	  conf_new[num+size*2].e00={x,0.};
	  conf_new[num+size*2].e01={y,0.};
	  conf_new[num+size*2].e02={z,0.};
	  conf_new[num+size*2].e10={t,0.};
	  conf_new[num+size*2].e11={0.,0.};
	  conf_new[num+size*2].e12={0.,0.};
	  conf_new[num+size*2].e20={0.,0.};
	  conf_new[num+size*2].e21={num,0.};
	  conf_new[num+size*2].e22={even,0.};
	  conf_new[num+size*3].e00={x,0.};
	  conf_new[num+size*3].e01={y,0.};
	  conf_new[num+size*3].e02={z,0.};
	  conf_new[num+size*3].e10={t,0.};
	  conf_new[num+size*3].e11={0.,0.};
	  conf_new[num+size*3].e12={0.,0.};
	  conf_new[num+size*3].e20={0.,0.};
	  conf_new[num+size*3].e21={num,0.};
	  conf_new[num+size*3].e22={even,0.};
	  */
	}
	//Now we can write conf_new to the file
	std::ofstream file(outputfile.c_str());
	file << ns << " " << ns << " " << ns << " " << nt << " ";
	file << params.get_beta() << " " << params.get_kappa() << " 12345" << std::endl;
	//The last number that I set to 12345 should be the hmc iteration; here it is not relevant
	for(int i=0; i<ns*ns*ns*nt*4; i++)
	  file << matrix_to_string(conf_new[i]);
	file.close();
	/*
        //To keep only links in one direction: it is usefull only for debugging
	for(int i=0; i<ns*ns*ns*nt*4; i++){
	  if((i%4)!=0){
	    conf_old[i].e00.re=0.;
	    conf_old[i].e00.im=0.;
	    conf_old[i].e01.re=0.;
	    conf_old[i].e01.im=0.;
	    conf_old[i].e02.re=0.;
	    conf_old[i].e02.im=0.;
	    conf_old[i].e10.re=0.;
	    conf_old[i].e10.im=0.;
	    conf_old[i].e11.re=0.;
	    conf_old[i].e11.im=0.;
	    conf_old[i].e12.re=0.;
	    conf_old[i].e12.im=0.;
	    conf_old[i].e20.re=0.;
	    conf_old[i].e20.im=0.;
	    conf_old[i].e21.re=0.;
	    conf_old[i].e21.im=0.;
	    conf_old[i].e22.re=0.;
	    conf_old[i].e22.im=0.;
	  }
	}
	cpu->get_gf_code()->importGaugefield(cpu->get_gaugefield(),conf_old);
	cpu->save_conf();
	*/
}

/**
 *  In the reference code the lattice is reorganized in the following way:
 * 
 *   0            size
 *   |------|------|
 *      e      o    
 * 
 *  where e=even, o=odd, whereas size=VOL4D.
 *  Hence, in order to use the same random staggered field in tests
 *  I have to print it to a text file according this scheme. 
 * 
 */
void print_staggeredfield_to_textfile(std::string outputfile, su3vec * sf, meta::Inputparameters params)
{
	int nt=params.get_ntime();
	int ns=params.get_nspace();
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
	  copy_su3vec(sf_new[num],sf[i]);
	}	
	//Now we can write sf_new to the file
	std::ofstream file(outputfile.c_str());
	file << ns << " " << ns << " " << ns << " " << nt << std::endl;
	for(int i=0; i<ns*ns*ns*nt; i++)
	  file << su3vec_to_string(sf_new[i]);
	file.close();
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_fermions_staggered";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	BOOST_MESSAGE("Test done");
}

void test_m_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "m_staggered";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	//The following three lines are to be used to produce the ref_conf file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	/*
	print_gaugefield_to_textfile("ref_conf",&cpu,params);
	logger.info() << "Produced the ref_conf text file with the links for the ref. code. Returning...";
	return;
	*/
	cl_int err = CL_SUCCESS;
	const hardware::code::Fermions_staggered * device = cpu.get_device();
	su3vec * sf_in;
	su3vec * sf_out;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);

	sf_in = new su3vec[NUM_ELEMENTS_SF];
	sf_out = new su3vec[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);
	//The following three lines are to be used to produce the ref_vec file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	 /*
	print_staggeredfield_to_textfile("ref_vec",sf_in,params);
	logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";
	return;
	 */

	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());
	out.load(sf_in);
	const hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_staggered_code();
		
	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	device->M_staggered_device(&in, &out,  cpu.get_gaugefield(), params.get_kappa());
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	
	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("/opencl_module_fermions_staggered_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/opencl_module_fermions_staggered_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_STAGGERED )

BOOST_AUTO_TEST_CASE( M_STAGGERED_1)
{
	test_m_staggered("/m_staggered_input_1");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_2)
{
	test_m_staggered("/m_staggered_input_2");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_3)
{
	test_m_staggered("/m_staggered_input_3");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_4)
{
	test_m_staggered("/m_staggered_input_4");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_5)
{
	test_m_staggered("/m_staggered_input_5");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_6)
{
	test_m_staggered("/m_staggered_input_6");
}

BOOST_AUTO_TEST_SUITE_END()

