#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE m_compare_noneo_eo
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "M_compare_noneo_eo";

//CP: this test is related to issue#282

class Device : public Opencl_Module_Hmc {

	meta::Counter counter1, counter2, counter3, counter4;
public:
	Device(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module_Hmc(params, device, &counter1, &counter2, &counter3, &counter4) {
		Opencl_Module_Hmc::init(); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel2(const hardware::buffers::Plain<spinor> * , const hardware::buffers::Spinor * , const hardware::buffers::Spinor * , const hardware::buffers::SU3 * gf, int gs, int ls, hmc_float kappa, hmc_float);
	void runTestKernel(const hardware::buffers::Plain<spinor> * out, const hardware::buffers::Plain<spinor> * in1, const hardware::buffers::SU3 * gf, int gs, int ls, hmc_float kappa, hmc_float);
	void fill_kernels();
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc(exec_name.c_str(), inputfile);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm_eo(int which);
	hmc_float get_squarenorm_noneo(int which);
	void verify(hmc_float, hmc_float);
	void verify_converted_vectors();
	void verify_output_vectors();
	void runTestKernel();
	void runTestKernel2();
	void runTestKernelwithconvertedfields();
	void runTestKernel2withconvertedfields();
	void reset_outfield_eo();
	void reset_outfield_eo_converted();
	void reset_outfield_noneo();
	void reset_outfield_noneo_converted();
	void convert_to_noneo(spinor * sf_noneo, spinor* eo1, spinor * eo2);
	void convert_to_eo(spinor * sf_noneo, spinor* eo1, spinor * eo2);
	void convert_input();

private:
	void fill_buffers();
	void clear_buffers();
	const hardware::buffers::Spinor * in_eo1, * in_eo2;
	const hardware::buffers::Spinor * in_eo1_converted, * in_eo2_converted;
	const hardware::buffers::Plain<spinor> * in_noneo, * in_noneo_converted;
	const hardware::buffers::Plain<spinor> * out_noneo, * out_noneo_converted, * out_eo, * out_eo_converted;
	spinor * sf_in_noneo;
	spinor * sf_in_noneo_converted;
	spinor * sf_out_noneo;
	spinor * sf_out_noneo_converted;

	spinor * sf_in_eo1;
	spinor * sf_in_eo1_converted;
	spinor * sf_in_eo2;
	spinor * sf_in_eo2_converted;

	spinor * sf_out1_eo;
	spinor * sf_out1_eo_converted;
	spinor * sf_out2_eo;
	spinor * sf_out2_eo_converted;

	spinor * sf_out_eo;
	spinor * sf_out_eo_converted;

	const hardware::buffers::Plain<hmc_float> * sqnorm;
};

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(get_parameters(), get_device_for_task(0));

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

//it is assumed that idx iterates only over half the number of sites
void get_even_site(int idx, int * out_space, int * out_t, const meta::Inputparameters & params)
{
	const size_t NSPACE = params.get_nspace();
	const size_t VOLSPACE = meta::get_volspace(params);
	int x, y, z, t;
	x = idx;
	t = (int)(idx / (VOLSPACE / 2));
	x -= t * VOLSPACE / 2;
	z = (int)(x / (NSPACE * NSPACE / 2));
	x -= z * NSPACE * NSPACE / 2;
	y = (int)(x / NSPACE);
	x -= y * NSPACE;
	(*out_space) =  (int)((z + t) % 2) * (1 + 2 * x - (int) (2 * x / NSPACE)) + (int)((t + z + 1) % 2) * (2 * x + (int) (2 * x / NSPACE)) + 2 * NSPACE * y + NSPACE * NSPACE * z;
	(*out_t) = t;
}

//it is assumed that idx iterates only over half the number of sites

void get_odd_site(int idx, int * out_space, int * out_t, const meta::Inputparameters& params)
{
	const size_t NSPACE = params.get_nspace();
	const size_t VOLSPACE = meta::get_volspace(params);
	int x, y, z, t;
	x = idx;
	t = (int)(idx / (VOLSPACE / 2));
	x -= t * VOLSPACE / 2;
	z = (int)(x / (NSPACE * NSPACE / 2));
	x -= z * NSPACE * NSPACE / 2;
	y = (int)(x / NSPACE);
	x -= y * NSPACE;

	(*out_space) =  (int)((z + t + 1) % 2) * (1 + 2 * x - (int) (2 * x / NSPACE)) + (int)((t + z) % 2) * (2 * x + (int) (2 * x / NSPACE)) + 2 * NSPACE * y + NSPACE * NSPACE * z;
	(*out_t) = t;
}


void fill_sf_with_one(spinor * sf_in1, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.e0 = hmc_complex_one;
		sf_in1[i].e0.e1 = hmc_complex_one;
		sf_in1[i].e0.e2 = hmc_complex_one;
		sf_in1[i].e1.e0 = hmc_complex_one;
		sf_in1[i].e1.e1 = hmc_complex_one;
		sf_in1[i].e1.e2 = hmc_complex_one;
		sf_in1[i].e2.e0 = hmc_complex_one;
		sf_in1[i].e2.e1 = hmc_complex_one;
		sf_in1[i].e2.e2 = hmc_complex_one;
		sf_in1[i].e3.e0 = hmc_complex_one;
		sf_in1[i].e3.e1 = hmc_complex_one;
		sf_in1[i].e3.e2 = hmc_complex_one;
	}
	return;
}

void fill_sf_with_float(spinor * sf_in, int size, hmc_float val)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0.re = val;
		sf_in[i].e0.e1.re = val;
		sf_in[i].e0.e2.re = val;
		sf_in[i].e1.e0.re = val;
		sf_in[i].e1.e1.re = val;
		sf_in[i].e1.e2.re = val;
		sf_in[i].e2.e0.re = val;
		sf_in[i].e2.e1.re = val;
		sf_in[i].e2.e2.re = val;
		sf_in[i].e3.e0.re = val;
		sf_in[i].e3.e1.re = val;
		sf_in[i].e3.e2.re = val;

		sf_in[i].e0.e0.im = val;
		sf_in[i].e0.e1.im = val;
		sf_in[i].e0.e2.im = val;
		sf_in[i].e1.e0.im = val;
		sf_in[i].e1.e1.im = val;
		sf_in[i].e1.e2.im = val;
		sf_in[i].e2.e0.im = val;
		sf_in[i].e2.e1.im = val;
		sf_in[i].e2.e2.im = val;
		sf_in[i].e3.e0.im = val;
		sf_in[i].e3.e1.im = val;
		sf_in[i].e3.e2.im = val;
	}
	return;
}

spinor fill_spinor_with_specific_float(hmc_float val)
{
	spinor tmp;

	tmp.e0.e0.re =  1. / 2. * 1. / 2. * val;
	tmp.e0.e1.re =  1. / 2. * 1. / 4. *val;
	tmp.e0.e2.re =  1. / 2. * 1. / 6. *val;
	tmp.e0.e0.im =  1. / 2. * 1. / 3. *val;
	tmp.e0.e1.im =  1. / 2. * 1. / 5. *val;
	tmp.e0.e2.im =  1. / 2. * 1. / 7. *val;

	tmp.e1.e0.re =  1. / 8. * 1. / 2. *val;
	tmp.e1.e1.re =  1. / 8. * 1. / 4. *val;
	tmp.e1.e2.re =  1. / 8. * 1. / 6. *val;
	tmp.e1.e0.im =  1. / 8. * 1. / 3. * val;
	tmp.e1.e1.im =  1. / 8. * 1. / 5. *val;
	tmp.e1.e2.im =  1. / 8. * 1. / 7. *val;

	tmp.e2.e0.re =  1. / 9. * 1. / 2. *val;
	tmp.e2.e1.re =  1. / 9. * 1. / 4. *val;
	tmp.e2.e2.re =  1. / 9. * 1. / 6. *val;
	tmp.e2.e0.im =  1. / 9. * 1. / 3. *val;
	tmp.e2.e1.im =  1. / 9. * 1. / 5. *val;
	tmp.e2.e2.im =  1. / 9. * 1. / 7. *val;

	tmp.e3.e0.re =  1. / 11. * 1. / 2. *val;
	tmp.e3.e1.re =  1. / 11. * 1. / 4. *val;
	tmp.e3.e2.re =  1. / 11. * 1. / 6. *val;
	tmp.e3.e0.im =  1. / 11. * 1. / 3. *val;
	tmp.e3.e1.im =  1. / 11. * 1. / 5. *val;
	tmp.e3.e2.im =  1. / 11. * 1. / 7. *val;

	return tmp;
}

//this function fills every lattice site with a specific value depending
//on its site index.
void fill_noneo_sf_with_specific_float(spinor * sf_in, const meta::Inputparameters& params)
{
	int x, y, z, t;
	int ns = params.get_nspace();
	for(x = 0;  x < params.get_nspace(); x++) {
		for(y = 0;  y < params.get_nspace();  y++) {
			for(z = 0;  z < params.get_nspace();  z++) {
				for(t = 0; t < params.get_ntime(); ++t) {
					//this has to match the conventions in operations_geometry.cl!!!!
					int i = t * ns * ns * ns + x * ns * ns + y * ns + z;

					hmc_float val = 17. / (i + 1);

					sf_in[i] = fill_spinor_with_specific_float(val);

				}
			}
		}
	}
	return;
}

void fill_eo_sf_with_specific_float(spinor * sf_even, spinor * sf_odd, const meta::Inputparameters& params)
{
	int x, y, z, t;
	int ns = params.get_nspace();
	for(x = 0;  x < params.get_nspace(); x++) {
		for(y = 0;  y < params.get_nspace();  y++) {
			for(z = 0;  z < params.get_nspace();  z++) {
				for(t = 0; t < params.get_ntime(); ++t) {
					//this has to match the conventions in operations_geometry.cl!!!!
					int i = t * ns * ns * ns + x * ns * ns + y * ns + z;
					hmc_float val = 17. / (i + 1);

					//distinguish between even and odd fields
					if( (t + x + y + z) % 2 == 0)
						sf_even[i / 2] = fill_spinor_with_specific_float(val);
					else
						sf_odd[i / 2] = fill_spinor_with_specific_float(val);
				}
			}
		}
	}
	return;
}


void fill_with_one(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 1.;
	}
	return;
}

hmc_complex make_complex(const hmc_float re, const hmc_float im)
{
	const hmc_complex tmp = {re, im};
	return tmp;
}

void fill_with_zero(spinor * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0 = make_complex(0., 0.);
		sf_in[i].e0.e1 = make_complex(0., 0.);
		sf_in[i].e0.e2 = make_complex(0., 0.);
		sf_in[i].e1.e0 = make_complex(0., 0.);
		sf_in[i].e1.e1 = make_complex(0., 0.);
		sf_in[i].e1.e2 = make_complex(0., 0.);
		sf_in[i].e2.e0 = make_complex(0., 0.);
		sf_in[i].e2.e1 = make_complex(0., 0.);
		sf_in[i].e2.e2 = make_complex(0., 0.);
		sf_in[i].e3.e0 = make_complex(0., 0.);
		sf_in[i].e3.e1 = make_complex(0., 0.);
		sf_in[i].e3.e2 = make_complex(0., 0.);
	}
	return;
}

void fill_sf_with_random_eo(spinor * sf_in1, spinor * sf_in2, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.e0.re = prng_double();
		sf_in1[i].e0.e1.re = prng_double();
		sf_in1[i].e0.e2.re = prng_double();
		sf_in1[i].e1.e0.re = prng_double();
		sf_in1[i].e1.e1.re = prng_double();
		sf_in1[i].e1.e2.re = prng_double();
		sf_in1[i].e2.e0.re = prng_double();
		sf_in1[i].e2.e1.re = prng_double();
		sf_in1[i].e2.e2.re = prng_double();
		sf_in1[i].e3.e0.re = prng_double();
		sf_in1[i].e3.e1.re = prng_double();
		sf_in1[i].e3.e2.re = prng_double();

		sf_in1[i].e0.e0.im = prng_double();
		sf_in1[i].e0.e1.im = prng_double();
		sf_in1[i].e0.e2.im = prng_double();
		sf_in1[i].e1.e0.im = prng_double();
		sf_in1[i].e1.e1.im = prng_double();
		sf_in1[i].e1.e2.im = prng_double();
		sf_in1[i].e2.e0.im = prng_double();
		sf_in1[i].e2.e1.im = prng_double();
		sf_in1[i].e2.e2.im = prng_double();
		sf_in1[i].e3.e0.im = prng_double();
		sf_in1[i].e3.e1.im = prng_double();
		sf_in1[i].e3.e2.im = prng_double();

		//    sf_in2[i] = sf_in1[i];

		sf_in2[i].e0.e0.re = prng_double();
		sf_in2[i].e0.e1.re = prng_double();
		sf_in2[i].e0.e2.re = prng_double();
		sf_in2[i].e1.e0.re = prng_double();
		sf_in2[i].e1.e1.re = prng_double();
		sf_in2[i].e1.e2.re = prng_double();
		sf_in2[i].e2.e0.re = prng_double();
		sf_in2[i].e2.e1.re = prng_double();
		sf_in2[i].e2.e2.re = prng_double();
		sf_in2[i].e3.e0.re = prng_double();
		sf_in2[i].e3.e1.re = prng_double();
		sf_in2[i].e3.e2.re = prng_double();

		sf_in2[i].e0.e0.im = prng_double();
		sf_in2[i].e0.e1.im = prng_double();
		sf_in2[i].e0.e2.im = prng_double();
		sf_in2[i].e1.e0.im = prng_double();
		sf_in2[i].e1.e1.im = prng_double();
		sf_in2[i].e1.e2.im = prng_double();
		sf_in2[i].e2.e0.im = prng_double();
		sf_in2[i].e2.e1.im = prng_double();
		sf_in2[i].e2.e2.im = prng_double();
		sf_in2[i].e3.e0.im = prng_double();
		sf_in2[i].e3.e1.im = prng_double();
		sf_in2[i].e3.e2.im = prng_double();
	}
	return;
}

void fill_sf_with_random_noneo(spinor * sf_in, int size, int seed, const meta::Inputparameters& params)
{
	prng_init(seed);
	//the simple for loop through the vector does not give the right vector compared to the eo-converted one...
	//this is because if the vector is simply filled with random numbers, the nontrivial
	//even-odd structure is not contained!!
	//However, one gets the same content if one uses a loop over the eoprec spinorfieldsize and then
	//updates 2 spinors at ones
	for(int i = 0; i < size / 2; i++) {
		int n, t, global_pos;
		get_even_site(i, &n, &t, params);
		global_pos = get_global_pos(n, t, params);
		sf_in[global_pos].e0.e0.re = prng_double();
		sf_in[global_pos].e0.e1.re = prng_double();
		sf_in[global_pos].e0.e2.re = prng_double();
		sf_in[global_pos].e1.e0.re = prng_double();
		sf_in[global_pos].e1.e1.re = prng_double();
		sf_in[global_pos].e1.e2.re = prng_double();
		sf_in[global_pos].e2.e0.re = prng_double();
		sf_in[global_pos].e2.e1.re = prng_double();
		sf_in[global_pos].e2.e2.re = prng_double();
		sf_in[global_pos].e3.e0.re = prng_double();
		sf_in[global_pos].e3.e1.re = prng_double();
		sf_in[global_pos].e3.e2.re = prng_double();

		sf_in[global_pos].e0.e0.im = prng_double();
		sf_in[global_pos].e0.e1.im = prng_double();
		sf_in[global_pos].e0.e2.im = prng_double();
		sf_in[global_pos].e1.e0.im = prng_double();
		sf_in[global_pos].e1.e1.im = prng_double();
		sf_in[global_pos].e1.e2.im = prng_double();
		sf_in[global_pos].e2.e0.im = prng_double();
		sf_in[global_pos].e2.e1.im = prng_double();
		sf_in[global_pos].e2.e2.im = prng_double();
		sf_in[global_pos].e3.e0.im = prng_double();
		sf_in[global_pos].e3.e1.im = prng_double();
		sf_in[global_pos].e3.e2.im = prng_double();

		get_odd_site(i, &n, &t, params);
		global_pos = get_global_pos(n, t, params);
		sf_in[global_pos].e0.e0.re = prng_double();
		sf_in[global_pos].e0.e1.re = prng_double();
		sf_in[global_pos].e0.e2.re = prng_double();
		sf_in[global_pos].e1.e0.re = prng_double();
		sf_in[global_pos].e1.e1.re = prng_double();
		sf_in[global_pos].e1.e2.re = prng_double();
		sf_in[global_pos].e2.e0.re = prng_double();
		sf_in[global_pos].e2.e1.re = prng_double();
		sf_in[global_pos].e2.e2.re = prng_double();
		sf_in[global_pos].e3.e0.re = prng_double();
		sf_in[global_pos].e3.e1.re = prng_double();
		sf_in[global_pos].e3.e2.re = prng_double();

		sf_in[global_pos].e0.e0.im = prng_double();
		sf_in[global_pos].e0.e1.im = prng_double();
		sf_in[global_pos].e0.e2.im = prng_double();
		sf_in[global_pos].e1.e0.im = prng_double();
		sf_in[global_pos].e1.e1.im = prng_double();
		sf_in[global_pos].e1.e2.im = prng_double();
		sf_in[global_pos].e2.e0.im = prng_double();
		sf_in[global_pos].e2.e1.im = prng_double();
		sf_in[global_pos].e2.e2.im = prng_double();
		sf_in[global_pos].e3.e0.im = prng_double();
		sf_in[global_pos].e3.e1.im = prng_double();
		sf_in[global_pos].e3.e2.im = prng_double();
	}
	return;
}

void Dummyfield::convert_to_noneo(spinor * sf_noneo, spinor* eo1, spinor * eo2)
{

	int pos, t;
	for(int n = 0; n < meta::get_eoprec_spinorfieldsize(get_parameters()); n++) {
		get_even_site(n, &pos, &t, get_parameters());
		sf_noneo[get_global_pos(pos, t, get_parameters())] = eo1[n];
		get_odd_site(n, &pos, &t, get_parameters());
		sf_noneo[get_global_pos(pos, t, get_parameters())] = eo2[n];
	}
	return;
}

void Dummyfield::convert_to_eo(spinor * sf_noneo, spinor* eo1, spinor * eo2)
{

	int pos, t;
	for(int n = 0; n < meta::get_eoprec_spinorfieldsize(get_parameters()); n++) {
		get_even_site(n, &pos, &t, get_parameters());
		eo1[n] = sf_noneo[get_global_pos(pos, t, get_parameters())];
		get_odd_site(n, &pos, &t, get_parameters());
		eo2[n] = sf_noneo[get_global_pos(pos, t, get_parameters())];
	}
	return;
}

void Dummyfield::reset_outfield_eo()
{
	out_eo->load(sf_out_eo);
}

void Dummyfield::reset_outfield_eo_converted()
{
	out_eo_converted->load(sf_out_eo_converted);
}

void Dummyfield::reset_outfield_noneo()
{
	out_noneo->load(sf_out_noneo);
}

void Dummyfield::reset_outfield_noneo_converted()
{
	out_noneo_converted->load(sf_out_noneo_converted);
}

bool compare_entries(hmc_complex in1, hmc_complex in2)
{
	using namespace std;

	if(abs(in1.re - in2.re) > 1e-8) {
		cout << "\treal entries DO NOT match: in1: " << in1.re << "\tin2: " << in2.re << endl;
		return false;
	}
	if (abs(in1.im - in2.im) > 1e-8) {
		cout << "\timag entries DO NOT match: in1: " << in1.im << "\tin2: " << in2.im << endl;
		return false;
	}
	return true;
}

bool compare_entries_float(hmc_float in1, hmc_float in2)
{
	using namespace std;

	if(abs(in1 - in2) > 1e-8) {
		cout << "\tentries DO NOT match: in1: " << in1 << "\tin2: " << in2 << endl;
		return false;
	}
	return true;
}

void print_spinor(const spinor in1, const spinor in2)
{
	using namespace std;

	cout << endl;
	cout << "(" << in1.e0.e0.re << "),(" << in1.e0.e0.im << ") : (" << in2.e0.e0.re << "),(" << in2.e0.e0.im << ")"  << endl;
	cout << "(" << in1.e0.e1.re << "),(" << in1.e0.e1.im << ") : (" << in2.e0.e1.re << "),(" << in2.e0.e1.im << ")"  << endl;
	cout << "(" << in1.e0.e2.re << "),(" << in1.e0.e2.im << ") : (" << in2.e0.e2.re << "),(" << in2.e0.e2.im << ")"  << endl;
	cout << "(" << in1.e1.e0.re << "),(" << in1.e1.e0.im << ") : (" << in2.e1.e0.re << "),(" << in2.e1.e0.im << ")"  << endl;
	cout << "(" << in1.e1.e1.re << "),(" << in1.e1.e1.im << ") : (" << in2.e1.e1.re << "),(" << in2.e1.e1.im << ")"  << endl;
	cout << "(" << in1.e1.e2.re << "),(" << in1.e1.e2.im << ") : (" << in2.e1.e2.re << "),(" << in2.e1.e2.im << ")"  << endl;
	cout << "(" << in1.e2.e0.re << "),(" << in1.e2.e0.im << ") : (" << in2.e2.e0.re << "),(" << in2.e2.e0.im << ")"  << endl;
	cout << "(" << in1.e2.e1.re << "),(" << in1.e2.e1.im << ") : (" << in2.e2.e1.re << "),(" << in2.e2.e1.im << ")"  << endl;
	cout << "(" << in1.e2.e2.re << "),(" << in1.e2.e2.im << ") : (" << in2.e2.e2.re << "),(" << in2.e2.e2.im << ")"  << endl;
	cout << "(" << in1.e3.e0.re << "),(" << in1.e3.e0.im << ") : (" << in2.e3.e0.re << "),(" << in2.e3.e0.im << ")"  << endl;
	cout << "(" << in1.e3.e1.re << "),(" << in1.e3.e1.im << ") : (" << in2.e3.e1.re << "),(" << in2.e3.e1.im << ")"  << endl;
	cout << "(" << in1.e3.e2.re << "),(" << in1.e3.e2.im << ") : (" << in2.e3.e2.re << "),(" << in2.e3.e2.im << ")"  << endl;
	return;
}

bool compare_vectors(spinor * in1, spinor * in2, int size)
{
	using namespace std;

	bool check;
	bool returner = true;
	for(int i = 0; i < size; i++) {
		check = compare_entries(in1[i].e0.e0, in2[i].e0.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e0.e1, in2[i].e0.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e0.e2, in2[i].e0.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e1.e0, in2[i].e1.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e1.e1, in2[i].e1.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e1.e2, in2[i].e1.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e2.e0, in2[i].e2.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e2.e1, in2[i].e2.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e2.e2, in2[i].e2.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e3.e0, in2[i].e3.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e3.e1, in2[i].e3.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
		check = compare_entries(in1[i].e3.e2, in2[i].e3.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		if(!check) returner = false;
	}
	return returner;
}

void Dummyfield::verify_converted_vectors()
{
	//copy vectors to host
	in_noneo_converted->dump(sf_in_noneo_converted);
	in_eo1_converted->dump(sf_in_eo1_converted);
	in_eo2_converted->dump(sf_in_eo2_converted);

	bool result = true;
	logger.info() << "\tcompare in_noneo with in_noneo_converted";
	result = compare_vectors(sf_in_noneo, sf_in_noneo_converted, meta::get_spinorfieldsize(get_parameters()));
	if(result) logger.info() << "\t\t...passed!";
	else {
		logger.info() << "\t\t...did not pass!";
		BOOST_REQUIRE_EQUAL(1, 0);
	}
	logger.info() << "\tcompare in_eo1 with in_eo1_converted";
	result = compare_vectors(sf_in_eo1, sf_in_eo1_converted, meta::get_eoprec_spinorfieldsize(get_parameters()));
	if(result) logger.info() << "\t\t...passed!";
	else {
		logger.info() << "\t\t...did not pass!";
		BOOST_REQUIRE_EQUAL(1, 0);
	}

	logger.info() << "\tcompare in_eo1 with in_eo2_converted";
	result = compare_vectors(sf_in_eo2, sf_in_eo2_converted, meta::get_eoprec_spinorfieldsize(get_parameters()));
	if(result) logger.info() << "\t\t...passed!";
	else {
		logger.info() << "\t\t...did not pass!";
		BOOST_REQUIRE_EQUAL(1, 0);
	}

	return;
}

void Dummyfield::verify_output_vectors()
{
	//copy vectors to host
	out_noneo->dump(sf_out_noneo);
	out_noneo_converted->dump(sf_out_noneo_converted);
	out_eo->dump(sf_out_eo);
	out_eo_converted->dump(sf_out_eo_converted);

	bool result = true;
	logger.info() << "\tcompare out_noneo with out_noneo_converted";
	result = compare_vectors(sf_out_noneo, sf_out_noneo_converted, meta::get_spinorfieldsize(get_parameters()));
	if(result) logger.info() << "\t\t...passed!";
	else {
		logger.info() << "\t\t...did not pass!";
		BOOST_REQUIRE_EQUAL(1, 0);
	}

	logger.info() << "\tcompare out_eo with out_noneo_converted";
	result = compare_vectors(sf_out_eo, sf_out_noneo_converted, meta::get_spinorfieldsize(get_parameters()));
	if(result) logger.info() << "\t\t...passed!";
	else {
		logger.info() << "\t\t...did not pass!";
		BOOST_REQUIRE_EQUAL(1, 0);
	}

	logger.info() << "\tcompare out_noneo with out_eo_converted";
	result = compare_vectors(sf_out_noneo, sf_out_eo_converted, meta::get_spinorfieldsize(get_parameters()));
	if(result) logger.info() << "\t\t...passed!";
	else {
		logger.info() << "\t\t...did not pass!";
		BOOST_REQUIRE_EQUAL(1, 0);
	}

	logger.info() << "\tcompare out_noneo with out_noneo_converted";
	result = compare_vectors(sf_out_eo, sf_out_eo_converted, meta::get_spinorfieldsize(get_parameters()));
	if(result) logger.info() << "\t\t...passed!";
	else {
		logger.info() << "\t\t...did not pass!";
		BOOST_REQUIRE_EQUAL(1, 0);
	}

	return;
}


void Dummyfield::fill_buffers()
{
	using namespace hardware::buffers;

	// don't invoke parent function as we don't require the original buffers
	int NUM_ELEMENTS_SF_EO;
	int NUM_ELEMENTS_SF_NON_EO;
	NUM_ELEMENTS_SF_EO =  meta::get_eoprec_spinorfieldsize(get_parameters());
	NUM_ELEMENTS_SF_NON_EO =  meta::get_spinorfieldsize(get_parameters());

	sf_in_eo1 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in_eo1_converted = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in_eo2 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in_eo2_converted = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out2_eo = new spinor[NUM_ELEMENTS_SF_EO];

	sf_in_eo1_converted = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in_eo2_converted = new spinor[NUM_ELEMENTS_SF_EO];

	sf_out1_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out2_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out1_eo_converted = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out2_eo_converted = new spinor[NUM_ELEMENTS_SF_EO];

	sf_out_eo = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_out_eo_converted = new spinor[NUM_ELEMENTS_SF_NON_EO];

	sf_in_noneo = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_in_noneo_converted = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_out_noneo = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_out_noneo_converted = new spinor[NUM_ELEMENTS_SF_NON_EO];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		fill_eo_sf_with_specific_float(sf_in_eo1, sf_in_eo2, get_parameters());
		fill_noneo_sf_with_specific_float(sf_in_noneo, get_parameters());
	} else {
		fill_sf_with_random_eo(sf_in_eo1, sf_in_eo2, NUM_ELEMENTS_SF_EO, 123456);
		fill_sf_with_random_noneo(sf_in_noneo, NUM_ELEMENTS_SF_NON_EO, 123456, get_parameters());
	}

	BOOST_REQUIRE(sf_in_eo1);
	BOOST_REQUIRE(sf_in_eo2);
	BOOST_REQUIRE(sf_in_noneo);
	BOOST_REQUIRE(sf_in_eo1_converted);
	BOOST_REQUIRE(sf_in_eo2_converted);
	BOOST_REQUIRE(sf_in_noneo_converted);
	fill_with_zero(sf_out_eo, NUM_ELEMENTS_SF_NON_EO);
	fill_with_zero(sf_out_eo_converted, NUM_ELEMENTS_SF_NON_EO);
	fill_with_zero(sf_out_noneo, NUM_ELEMENTS_SF_NON_EO);
	fill_with_zero(sf_out_noneo_converted, NUM_ELEMENTS_SF_NON_EO);
	fill_with_zero(sf_out2_eo, NUM_ELEMENTS_SF_EO);
	fill_with_zero(sf_out1_eo, NUM_ELEMENTS_SF_EO);
	fill_with_zero(sf_out1_eo_converted, NUM_ELEMENTS_SF_EO);
	fill_with_zero(sf_out2_eo_converted, NUM_ELEMENTS_SF_EO);

	//create buffers
	Opencl_Module_Spinors * spinor_module = static_cast<Opencl_Module_Spinors*>(opencl_modules[0]);
	//create buffer for sf on device (and copy sf_in to both for convenience)
	in_eo1 = new Spinor(meta::get_eoprec_spinorfieldsize(get_parameters()), spinor_module->get_device());
	in_eo2 = new Spinor(meta::get_eoprec_spinorfieldsize(get_parameters()), spinor_module->get_device());
	in_eo1_converted = new Spinor(meta::get_eoprec_spinorfieldsize(get_parameters()), spinor_module->get_device());
	in_eo2_converted = new Spinor(meta::get_eoprec_spinorfieldsize(get_parameters()), spinor_module->get_device());

	in_noneo = new Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), spinor_module->get_device());
	in_noneo_converted = new Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), spinor_module->get_device());

	out_noneo = new Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), spinor_module->get_device());
	out_noneo_converted = new Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), spinor_module->get_device());
	out_eo = new Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), spinor_module->get_device());
	out_eo_converted = new Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), spinor_module->get_device());

	//copy content into buffers

	spinor_module->copy_to_eoprec_spinorfield_buffer(in_eo1, sf_in_eo1);
	spinor_module->copy_to_eoprec_spinorfield_buffer(in_eo2, sf_in_eo2);

	out_eo->load(sf_out_eo);
	in_noneo->load(sf_in_noneo);
	out_noneo->load(sf_out_noneo);

	this->reset_outfield_noneo();

	//create buffer for squarenorm on device
	sqnorm = new Plain<hmc_float>(1, spinor_module->get_device());
}

void Device::fill_kernels()
{
	Opencl_Module_Hmc::fill_kernels();
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	delete in_eo1;
	delete in_eo2;
	delete in_eo1_converted;
	delete in_eo2_converted;
	delete out_eo;
	delete sqnorm;

	delete[] sf_in_eo1;
	delete[] sf_in_eo2;
	delete[] sf_in_eo1_converted;
	delete[] sf_in_eo2_converted;
	delete[] sf_out_eo;

	delete in_noneo;
	delete in_noneo_converted;
	delete out_noneo;

	delete[] sf_in_noneo;
	delete[] sf_in_noneo_converted;
	delete[] sf_out_noneo;
	delete[] sf_out_noneo_converted;
}

void Device::clear_kernels()
{
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel2(const hardware::buffers::Plain<spinor> * out, const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2, const hardware::buffers::SU3 * gf, int gs, int ls, hmc_float kappa, hmc_float mubar)
{
	using namespace hardware::buffers;

	//suppose in1 is the even, in2 the odd input vector
	//create 3 buffers for intermediate results
	cl_int err;
	const Spinor out_tmp_eo1(meta::get_eoprec_spinorfieldsize(get_parameters()), get_device());
	const Spinor out_tmp_eo2(meta::get_eoprec_spinorfieldsize(get_parameters()), get_device());
	const Spinor tmp_eo(meta::get_eoprec_spinorfieldsize(get_parameters()), get_device());

	//create -1 on device
	hmc_complex minusone_tmp = { -1., 0.};
	hardware::buffers::Plain<hmc_complex> minusone(1, get_device());
	minusone.load(&minusone_tmp);

	//now calc out_tmp_eo1 = (R_even in1 + D_eo in2)
	this->set_zero_spinorfield_eoprec_device(&tmp_eo);
	this->set_zero_spinorfield_eoprec_device(&out_tmp_eo1);

	this->dslash_eo_device(in2, &out_tmp_eo1, gf, EO, kappa);
	this->M_tm_sitediagonal_device(in1, &tmp_eo, mubar);

	this->saxpy_eoprec_device(&out_tmp_eo1, &tmp_eo, &minusone, &out_tmp_eo1);

	//now calc out_tmp_eo2 = ( R_odd in2 + D_oe in1)
	this->set_zero_spinorfield_eoprec_device(&tmp_eo);
	this->set_zero_spinorfield_eoprec_device(&out_tmp_eo2);

	this->dslash_eo_device(in1, &out_tmp_eo2, gf, OE, kappa);
	this->M_tm_sitediagonal_device(in2, &tmp_eo, mubar);

	this->saxpy_eoprec_device(&out_tmp_eo2, &tmp_eo, &minusone, &out_tmp_eo2);

	//now, both output vectors have to be converted back to noneo
	this->convert_from_eoprec_device(&out_tmp_eo1, &out_tmp_eo2, out);
}

void Device::runTestKernel(const hardware::buffers::Plain<spinor> * out, const hardware::buffers::Plain<spinor> * in1, const hardware::buffers::SU3 * gf, int gs, int ls, hmc_float kappa, hmc_float mubar)
{
	this->M_tm_plus_device(in1, out, gf, kappa , mubar );
}

hmc_float Dummyfield::get_squarenorm_eo(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in_eo1, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in_eo2, sqnorm);
	if(which == 2) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(out_eo, sqnorm);
	if(which == 3) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in_eo1_converted, sqnorm);
	if(which == 4) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in_eo2_converted, sqnorm);
	if(which == 5) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(out_eo_converted, sqnorm);
	// get stuff from device
	hmc_float result;
	sqnorm->dump(&result);
	logger.info() << result;
	return result;
}

hmc_float Dummyfield::get_squarenorm_noneo(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in_noneo, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in_noneo_converted, sqnorm);
	if(which == 2) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(out_noneo, sqnorm);
	if(which == 3) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(out_noneo_converted, sqnorm);
	// get stuff from device
	hmc_float result;
	sqnorm->dump(&result);
	logger.info() << result;
	return result;
}

void Dummyfield::convert_input()
{
	static_cast<Device*>(opencl_modules[0])->convert_from_eoprec_device(in_eo1, in_eo2, in_noneo_converted)  ;
	static_cast<Device*>(opencl_modules[0])->convert_to_eoprec_device(in_eo1_converted, in_eo2_converted, in_noneo) ;
}

//identify cpu: eo and gpu: noneo results
void Dummyfield::verify(hmc_float cpu, hmc_float gpu)
{
	//this is too much required, since rounding errors can occur
	//  BOOST_REQUIRE_EQUAL(cpu, gpu);
	//instead, test if the two number agree within some percent
	hmc_float dev = (cpu - gpu) / cpu / 100.;
	if(abs(dev) < 1e-10) {
		logger.info() << "eo and non-eo result agree within accuary of " << 1e-10;
		logger.info() << "eo: " << cpu << "\tnoneo: " << gpu;
	} else {
		logger.info() << "eo and noneo result DO NOT agree within accuary of " << 1e-10;
		logger.info() << "eo: " << cpu << "\tnoneo: " << gpu;
		BOOST_REQUIRE_EQUAL(1, 0);
	}
}

//this calls eo dslash and eo diagonal matrix on the two eo inputs and saves it in one noneo vector
void Dummyfield::runTestKernel2()
{
	int gs = 0, ls = 0;
	if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_eoprec_spinorfieldsize(get_parameters());
		ls = 64;
	} else if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = get_device_for_task(0)->get_num_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel2(out_eo, in_eo1, in_eo2, device->get_gaugefield(), gs, ls, get_parameters().get_kappa(), meta::get_mubar(get_parameters()));
}

void Dummyfield::runTestKernel2withconvertedfields()
{
	int gs = 0, ls = 0;
	if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_eoprec_spinorfieldsize(get_parameters());
		ls = 64;
	} else if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = get_device_for_task(0)->get_num_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel2(out_eo_converted, in_eo1_converted, in_eo2_converted, device->get_gaugefield(), gs, ls, get_parameters().get_kappa(), meta::get_mubar(get_parameters()));
}

//this calls the noneo fermionmatrix M_tm_plus
void Dummyfield::runTestKernel()
{
	int gs = 0, ls = 0;
	if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_spinorfieldsize(get_parameters());
		ls = 64;
	} else if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = get_device_for_task(0)->get_num_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(out_noneo, in_noneo, device->get_gaugefield(), gs, ls, get_parameters().get_kappa(), meta::get_mubar(get_parameters()));
}

void Dummyfield::runTestKernelwithconvertedfields()
{
	int gs = 0, ls = 0;
	if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_spinorfieldsize(get_parameters());
		ls = 64;
	} else if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = get_device_for_task(0)->get_num_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(out_noneo_converted, in_noneo_converted, device->get_gaugefield(), gs, ls, get_parameters().get_kappa(), meta::get_mubar(get_parameters()));
}

BOOST_AUTO_TEST_CASE( M_noneo_eo_test )
{

	/*
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "Run kernel";
	cpu.runTestKernel();
	logger.info() << "result:";
	hmc_float cpu_res;
	cpu_res = cpu.get_squarenorm();
	logger.info() << cpu_res;

	logger.info() << "Choosing reference value and acceptance precision";
	hmc_float ref_val = params.get_test_ref_value();
	logger.info() << "reference value:\t" << ref_val;
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	logger.info() << "Compare result to reference value";
	BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
	logger.info() << "Done";
	BOOST_MESSAGE("Test done");
	*/
	logger.info() << "Test equivalence of M_tm_plus in eo- and non-eo formulation";

	int param_expect = 4;
	logger.info() << "expect parameters:";
	logger.info() << "\texec_name\tinputfile\tgpu_usage\trec12_usage";
	//get number of parameters
	int num_par = boost::unit_test::framework::master_test_suite().argc;
	if(num_par < param_expect) {
		logger.fatal() << "need more inputparameters! Got only " << num_par << ", expected " << param_expect << "! Aborting...";
		exit(-1);
	}

	//get input file that has been passed as an argument
	const char*  inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
	logger.info() << "inputfile used: " << inputfile;
	//get use_gpu = true/false that has been passed as an argument
	const char*  gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
	logger.info() << "GPU usage: " << gpu_opt;
	//get use_rec12 = true/false that has been passed as an argument
	const char* rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
	logger.info() << "rec12 usage: " << rec12_opt;

	logger.info() << "Init device";
	const char* _params_cpu[] = {"foo", inputfile, gpu_opt, rec12_opt};
	meta::Inputparameters params(param_expect, _params_cpu);
	hardware::System system(params);
	Dummyfield cpu(&system);

	if (params.get_use_gpu()) {
		logger.fatal() << "there is some problem with the input vectors if the GPU is used. To be checked!! Aborting..";
		BOOST_REQUIRE_EQUAL(1., 0.);
	}


	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);

	logger.info() << "##";
	logger.info() << "##";

	logger.info() << "convert input vectors";
	cpu.convert_input();
	logger.info() << "non-eo input:";
	logger.info() << "|phi|^2:";
	hmc_float cpu_back_noneo = cpu.get_squarenorm_noneo(0);

	logger.info() << "eo input:";
	logger.info() << "|phi_even|^2:";
	hmc_float cpu_back_eo = cpu.get_squarenorm_eo(0);
	logger.info() << "|phi_odd|^2:";
	hmc_float cpu_back2_eo = cpu.get_squarenorm_eo(1);

	logger.info() << "converted noneo input:";
	logger.info() << "|phi_conv|^2:";
	hmc_float cpu_back_noneo_converted = cpu.get_squarenorm_noneo(1);

	logger.info() << "converted eo input:";
	logger.info() << "|phi_even_converted|^2:";
	hmc_float cpu_back_eo_converted = cpu.get_squarenorm_eo(3);
	logger.info() << "|phi_odd_converted|^2:";
	hmc_float cpu_back2_eo_converted = cpu.get_squarenorm_eo(4);

	logger.info() << "##";
	logger.info() << "##";

	logger.info() << "Compare input vectors";
	cpu.verify(cpu_back_noneo_converted, cpu_back_noneo);
	cpu.verify(cpu_back_eo_converted, cpu_back_eo);
	cpu.verify(cpu_back2_eo_converted, cpu_back2_eo);
	logger.info() << "\tCompare non-eo and eo-converted input vectors entry by entry:";
	cpu.verify_converted_vectors();

	logger.info() << "##";
	logger.info() << "##";

	logger.info() << "run noneo fermionmatrix with noneo input...";
	cpu.reset_outfield_noneo();
	cpu.runTestKernel();
	logger.info() << "|M phi|^2:";
	hmc_float cpu_res_noneo;
	cpu_res_noneo = cpu.get_squarenorm_noneo(2);

	cpu.reset_outfield_noneo_converted();
	logger.info() << "run noneo fermionmatrix with converted eo input...";
	cpu.runTestKernelwithconvertedfields();
	logger.info() << "|M phi_converted|^2:";
	hmc_float cpu_res_noneo_converted;
	cpu_res_noneo_converted = cpu.get_squarenorm_noneo(3);

	logger.info() << "run eo fermionmatrix...";
	cpu.reset_outfield_eo();
	cpu.runTestKernel2();
	logger.info() << "|M_eo phi_eo|^2:";
	hmc_float cpu_res_eo;
	cpu_res_eo = cpu.get_squarenorm_eo(2);

	logger.info() << "run eo fermionmatrix with converted noneo input...";
	cpu.reset_outfield_eo_converted();
	cpu.runTestKernel2withconvertedfields();
	logger.info() << "|M_eo phi_eo_converted|^2:";
	hmc_float cpu_res_eo_converted;
	cpu_res_eo_converted = cpu.get_squarenorm_eo(5);

	logger.info() << "##";
	logger.info() << "##";

	logger.info() << "Compare eo and non-eo CPU results";
	logger.info() << "Compare output vectors:";
	logger.info() << "\teo and noneo_converted:";
	cpu.verify(cpu_res_eo, cpu_res_noneo_converted);
	logger.info() << "\tnoneo and noneo_converted:";
	cpu.verify(cpu_res_noneo_converted, cpu_res_noneo);
	logger.info() << "\teo and noneo:";
	cpu.verify(cpu_res_eo, cpu_res_noneo);
	logger.info() << "\tnoneo and eo_converted:";
	cpu.verify(cpu_res_noneo, cpu_res_eo_converted);
	logger.info() << "\tCompare output vectors entry by entry:";
	cpu.verify_output_vectors();

	logger.info() << "##";
	logger.info() << "##";
	BOOST_MESSAGE("Test done");
}
