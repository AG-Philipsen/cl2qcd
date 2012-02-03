#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fermionforce_compare_noneo_eo
#include <boost/test/unit_test.hpp>

Random rnd(15);
extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module_Hmc {

	cl_kernel testKernel;
	cl_kernel testKernel2;
	cl_kernel ae_sqn;
public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module_Hmc() {
		Opencl_Module_Hmc::init(queue, 0, params, maxcomp, double_ext); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem in1, cl_mem in2, cl_mem out, cl_mem gf, int gs, int ls, int evenodd);
	void runTestKernel2(cl_mem in1, cl_mem in2, cl_mem out, cl_mem gf, int gs, int ls);
	void fill_kernels();
	void set_float_to_gm_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out);
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type) : Gaugefield_hybrid() {
		std::stringstream tmp;
#ifdef _USEDOUBLEPREC_
		tmp << SOURCEDIR << "/tests/f_fermion_eo_input_1";
#else
		tmp << SOURCEDIR << "/tests/f_fermion_eo_input_1_single";
#endif
		params.readfile(tmp.str().c_str());

		init(1, device_type, &params);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm_eo(int which);
	hmc_float get_squarenorm_noneo(int which);
	void verify(hmc_float, hmc_float);
	void verify_converted_vectors();
	void runTestKernel(int evenodd);
	void runTestKernel2();
	void runTestKernel2withconvertedfields();
	void reset_outfield_eo();
	void reset_outfield_noneo();
	void convert_to_noneo(spinor * sf_noneo, spinor* eo1, spinor * eo2);
private:
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
	cl_mem in1_eo, in2_eo, in3_eo, in4_eo, out_eo;
	cl_mem in1_noneo, in2_noneo, out_noneo;
	cl_mem in1_noneo_converted, in2_noneo_converted;
	spinor * sf_in1_noneo;
	spinor * sf_in2_noneo;
	spinor * sf_in1_noneo_converted;
	spinor * sf_in2_noneo_converted;
	hmc_float * sf_out_noneo;
	cl_mem sqnorm;
	spinor * sf_in1_eo;
	spinor * sf_in2_eo;
	spinor * sf_in3_eo;
	spinor * sf_in4_eo;
	hmc_float * sf_out_eo;
	Matrixsu3 * gf_in;

};

BOOST_AUTO_TEST_CASE( F_FERMION )
{

	logger.info() << "Perform test only on CPU";
	logger.info() << "\tOther tests for fermion force shall test equivalence of GPU and CPU...";

	logger.info() << "Init CPU device";
	//params.print_info_inverter("m_gpu");
	Dummyfield cpu(CL_DEVICE_TYPE_CPU);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);

	logger.info() << "eo input:";
	logger.info() << "|phi_even_1|^2:";
	//hmc_float cpu_back_eo = cpu.get_squarenorm_eo(0);
	logger.info() << "|phi_even_2|^2:";
	//hmc_float cpu_back2_eo = cpu.get_squarenorm_eo(1);
	logger.info() << "|phi_odd_1|^2:";
	//hmc_float cpu_back3_eo = cpu.get_squarenorm_eo(2);
	logger.info() << "|phi_odd_2|^2:";
	//hmc_float cpu_back4_eo = cpu.get_squarenorm_eo(3);

	logger.info() << "run eo force on EVEN and ODD sites...";
	cpu.runTestKernel(EVEN);
	cpu.runTestKernel(ODD);
	logger.info() << "|force (even) + force (odd)|^2:";
	hmc_float cpu_res_eo;
	cpu_res_eo = cpu.get_squarenorm_eo(4);

	logger.info() << "non-eo input:";
	logger.info() << "|phi_1|^2:";
	hmc_float cpu_back_noneo = cpu.get_squarenorm_noneo(0);
	logger.info() << "|phi_2|^2:";
	hmc_float cpu_back2_noneo = cpu.get_squarenorm_noneo(1);

	logger.info() << "run noneo force with noneo input...";
	//cpu.runTestKernel2();
	logger.info() << "|force_noneo|^2:";
	hmc_float cpu_res_noneo;
	cpu_res_noneo = cpu.get_squarenorm_noneo(2);

	logger.info() << "converted non-eo input:";
	logger.info() << "|phi_1|^2:";
	hmc_float cpu_back_noneo_converted = cpu.get_squarenorm_noneo(3);
	logger.info() << "|phi_2|^2:";
	hmc_float cpu_back2_noneo_converted = cpu.get_squarenorm_noneo(4);

	cpu.reset_outfield_noneo();
	logger.info() << "run noneo force with converted eo input...";
	cpu.runTestKernel2withconvertedfields();
	logger.info() << "|force_noneo|^2:";
	hmc_float cpu_res_noneo_converted;
	cpu_res_noneo_converted = cpu.get_squarenorm_noneo(2);

	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Compare eo and non-eo CPU results";
	logger.info() << "Input vectors (compare only converted eo vectors to noneo):";
	cpu.verify(cpu_back_noneo_converted, cpu_back_noneo);
	cpu.verify(cpu_back2_noneo_converted, cpu_back2_noneo);
	logger.info() << "Compare non-eo and eo-converted input vectors entry by entry:";
	cpu.verify_converted_vectors();
	logger.info() << "Output vectors:";
	logger.info() << "eo and noneo_converted:";
	cpu.verify(cpu_res_eo, cpu_res_noneo_converted);
	logger.info() << "noneo and noneo_converted:";
	cpu.verify(cpu_res_noneo_converted, cpu_res_noneo);
	logger.info() << "eo and noneo:";
	cpu.verify(cpu_res_eo, cpu_res_noneo);


}

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(queue[0], get_parameters(), get_max_compute_units(0), get_double_ext(0));

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

//it is assumed that idx iterates only over half the number of sites
void get_even_site(int idx, int * out_space, int * out_t, const inputparameters * const params)
{
	const size_t NSPACE = params->get_ns();
	const size_t VOLSPACE = params->get_volspace();
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

void get_odd_site(int idx, int * out_space, int * out_t, const inputparameters * const params)
{
	const size_t NSPACE = params->get_ns();
	const size_t VOLSPACE = params->get_volspace();
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
void fill_noneo_sf_with_specific_float(spinor * sf_in, inputparameters * params)
{
	int x, y, z, t;
	int ns = params->get_ns();
	for(x = 0;  x < params->get_ns(); x++) {
		for(y = 0;  y < params->get_ns();  y++) {
			for(z = 0;  z < params->get_ns();  z++) {
				for(t = 0; t < params->get_nt(); ++t) {
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

void fill_eo_sf_with_specific_float(spinor * sf_even, spinor * sf_odd, inputparameters * params)
{
	int x, y, z, t;
	int ns = params->get_ns();
	for(x = 0;  x < params->get_ns(); x++) {
		for(y = 0;  y < params->get_ns();  y++) {
			for(z = 0;  z < params->get_ns();  z++) {
				for(t = 0; t < params->get_nt(); ++t) {
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

void fill_with_zero(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 0.;
	}
	return;
}

void fill_sf_with_random_eo(spinor * sf_in1, spinor * sf_in2, int size, int seed)
{
	Random rnd_loc(seed);
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.e0.re = rnd_loc.doub();
		sf_in1[i].e0.e1.re = rnd_loc.doub();
		sf_in1[i].e0.e2.re = rnd_loc.doub();
		sf_in1[i].e1.e0.re = rnd_loc.doub();
		sf_in1[i].e1.e1.re = rnd_loc.doub();
		sf_in1[i].e1.e2.re = rnd_loc.doub();
		sf_in1[i].e2.e0.re = rnd_loc.doub();
		sf_in1[i].e2.e1.re = rnd_loc.doub();
		sf_in1[i].e2.e2.re = rnd_loc.doub();
		sf_in1[i].e3.e0.re = rnd_loc.doub();
		sf_in1[i].e3.e1.re = rnd_loc.doub();
		sf_in1[i].e3.e2.re = rnd_loc.doub();

		sf_in1[i].e0.e0.im = rnd_loc.doub();
		sf_in1[i].e0.e1.im = rnd_loc.doub();
		sf_in1[i].e0.e2.im = rnd_loc.doub();
		sf_in1[i].e1.e0.im = rnd_loc.doub();
		sf_in1[i].e1.e1.im = rnd_loc.doub();
		sf_in1[i].e1.e2.im = rnd_loc.doub();
		sf_in1[i].e2.e0.im = rnd_loc.doub();
		sf_in1[i].e2.e1.im = rnd_loc.doub();
		sf_in1[i].e2.e2.im = rnd_loc.doub();
		sf_in1[i].e3.e0.im = rnd_loc.doub();
		sf_in1[i].e3.e1.im = rnd_loc.doub();
		sf_in1[i].e3.e2.im = rnd_loc.doub();

		//    sf_in2[i] = sf_in1[i];

		sf_in2[i].e0.e0.re = rnd_loc.doub();
		sf_in2[i].e0.e1.re = rnd_loc.doub();
		sf_in2[i].e0.e2.re = rnd_loc.doub();
		sf_in2[i].e1.e0.re = rnd_loc.doub();
		sf_in2[i].e1.e1.re = rnd_loc.doub();
		sf_in2[i].e1.e2.re = rnd_loc.doub();
		sf_in2[i].e2.e0.re = rnd_loc.doub();
		sf_in2[i].e2.e1.re = rnd_loc.doub();
		sf_in2[i].e2.e2.re = rnd_loc.doub();
		sf_in2[i].e3.e0.re = rnd_loc.doub();
		sf_in2[i].e3.e1.re = rnd_loc.doub();
		sf_in2[i].e3.e2.re = rnd_loc.doub();

		sf_in2[i].e0.e0.im = rnd_loc.doub();
		sf_in2[i].e0.e1.im = rnd_loc.doub();
		sf_in2[i].e0.e2.im = rnd_loc.doub();
		sf_in2[i].e1.e0.im = rnd_loc.doub();
		sf_in2[i].e1.e1.im = rnd_loc.doub();
		sf_in2[i].e1.e2.im = rnd_loc.doub();
		sf_in2[i].e2.e0.im = rnd_loc.doub();
		sf_in2[i].e2.e1.im = rnd_loc.doub();
		sf_in2[i].e2.e2.im = rnd_loc.doub();
		sf_in2[i].e3.e0.im = rnd_loc.doub();
		sf_in2[i].e3.e1.im = rnd_loc.doub();
		sf_in2[i].e3.e2.im = rnd_loc.doub();
	}
	return;
}

void fill_sf_with_random_noneo(spinor * sf_in, int size, int seed, inputparameters * params)
{
	Random rnd_loc(seed);
	//the simple for loop through the vector does not give the right vector compared to the eo-converted one...
	//this is because if the vector is simply filled with random numbers, the nontrivial
	//even-odd structure is not contained!!
	//However, one gets the same content if one uses a loop over the eoprec spinorfieldsize and then
	//updates 2 spinors at ones
	for(int i = 0; i < size / 2; i++) {
		int n, t, global_pos;
		get_even_site(i, &n, &t, params);
		global_pos = get_global_pos(n, t, params);
		sf_in[global_pos].e0.e0.re = rnd_loc.doub();
		sf_in[global_pos].e0.e1.re = rnd_loc.doub();
		sf_in[global_pos].e0.e2.re = rnd_loc.doub();
		sf_in[global_pos].e1.e0.re = rnd_loc.doub();
		sf_in[global_pos].e1.e1.re = rnd_loc.doub();
		sf_in[global_pos].e1.e2.re = rnd_loc.doub();
		sf_in[global_pos].e2.e0.re = rnd_loc.doub();
		sf_in[global_pos].e2.e1.re = rnd_loc.doub();
		sf_in[global_pos].e2.e2.re = rnd_loc.doub();
		sf_in[global_pos].e3.e0.re = rnd_loc.doub();
		sf_in[global_pos].e3.e1.re = rnd_loc.doub();
		sf_in[global_pos].e3.e2.re = rnd_loc.doub();

		sf_in[global_pos].e0.e0.im = rnd_loc.doub();
		sf_in[global_pos].e0.e1.im = rnd_loc.doub();
		sf_in[global_pos].e0.e2.im = rnd_loc.doub();
		sf_in[global_pos].e1.e0.im = rnd_loc.doub();
		sf_in[global_pos].e1.e1.im = rnd_loc.doub();
		sf_in[global_pos].e1.e2.im = rnd_loc.doub();
		sf_in[global_pos].e2.e0.im = rnd_loc.doub();
		sf_in[global_pos].e2.e1.im = rnd_loc.doub();
		sf_in[global_pos].e2.e2.im = rnd_loc.doub();
		sf_in[global_pos].e3.e0.im = rnd_loc.doub();
		sf_in[global_pos].e3.e1.im = rnd_loc.doub();
		sf_in[global_pos].e3.e2.im = rnd_loc.doub();

		get_odd_site(i, &n, &t, params);
		global_pos = get_global_pos(n, t, params);
		sf_in[global_pos].e0.e0.re = rnd_loc.doub();
		sf_in[global_pos].e0.e1.re = rnd_loc.doub();
		sf_in[global_pos].e0.e2.re = rnd_loc.doub();
		sf_in[global_pos].e1.e0.re = rnd_loc.doub();
		sf_in[global_pos].e1.e1.re = rnd_loc.doub();
		sf_in[global_pos].e1.e2.re = rnd_loc.doub();
		sf_in[global_pos].e2.e0.re = rnd_loc.doub();
		sf_in[global_pos].e2.e1.re = rnd_loc.doub();
		sf_in[global_pos].e2.e2.re = rnd_loc.doub();
		sf_in[global_pos].e3.e0.re = rnd_loc.doub();
		sf_in[global_pos].e3.e1.re = rnd_loc.doub();
		sf_in[global_pos].e3.e2.re = rnd_loc.doub();

		sf_in[global_pos].e0.e0.im = rnd_loc.doub();
		sf_in[global_pos].e0.e1.im = rnd_loc.doub();
		sf_in[global_pos].e0.e2.im = rnd_loc.doub();
		sf_in[global_pos].e1.e0.im = rnd_loc.doub();
		sf_in[global_pos].e1.e1.im = rnd_loc.doub();
		sf_in[global_pos].e1.e2.im = rnd_loc.doub();
		sf_in[global_pos].e2.e0.im = rnd_loc.doub();
		sf_in[global_pos].e2.e1.im = rnd_loc.doub();
		sf_in[global_pos].e2.e2.im = rnd_loc.doub();
		sf_in[global_pos].e3.e0.im = rnd_loc.doub();
		sf_in[global_pos].e3.e1.im = rnd_loc.doub();
		sf_in[global_pos].e3.e2.im = rnd_loc.doub();
	}
	return;
}

void Dummyfield::convert_to_noneo(spinor * sf_noneo, spinor* eo1, spinor * eo2)
{

	int pos, t;
	for(int n = 0; n < get_parameters()->get_eoprec_spinorfieldsize(); n++) {
		get_even_site(n, &pos, &t, get_parameters());
		sf_noneo[get_global_pos(pos, t, get_parameters())] = eo1[n];
		get_odd_site(n, &pos, &t, get_parameters());
		sf_noneo[get_global_pos(pos, t, get_parameters())] = eo2[n];
	}
	return;
}

void Dummyfield::reset_outfield_eo()
{
	size_t ae_buf_size = get_parameters()->get_gm_buf_size();
	int err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), out_eo, CL_TRUE, 0, ae_buf_size, sf_out_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	return;
}

void Dummyfield::reset_outfield_noneo()
{
	size_t ae_buf_size = get_parameters()->get_gm_buf_size();
	int err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), out_noneo, CL_TRUE, 0, ae_buf_size, sf_out_noneo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	return;
}

bool compare_entries(hmc_complex in1, hmc_complex in2)
{
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

void print_spinor(spinor in1, spinor in2)
{
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

void compare_vectors(spinor * in1, spinor * in2, int size)
{
	bool check;
	//int n, t;
	//int nodd, todd;
	for(int i = 0; i < size; i++) {
		//    get_even_site(i/2, &n, &t, params);
		//get_odd_site(i/2, &nodd, &todd, params);
		//    cout << "\npos: " << i << "\teoidx: " <<  i/2 << " even: " << "(" << n << "," << t << ") " << endl;
		//print_spinor(in1[i] , in2[i]);
		check = compare_entries(in1[i].e0.e0, in2[i].e0.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e0.e1, in2[i].e0.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e0.e2, in2[i].e0.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e1.e0, in2[i].e1.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e1.e1, in2[i].e1.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e1.e2, in2[i].e1.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e2.e0, in2[i].e2.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e2.e1, in2[i].e2.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e2.e2, in2[i].e2.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e3.e0, in2[i].e3.e0);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e3.e1, in2[i].e3.e1);
		if(!check) cout << "\terror occured at " << i << endl;
		check = compare_entries(in1[i].e3.e2, in2[i].e3.e2);
		if(!check) cout << "\terror occured at " << i << endl;
		/*
		cout << "\npos: " << i << "\teoidx: " << i/2 << " odd: (" << nodd << "," << todd << ")" << endl;
		print_spinor(in1[i+1] , in2[i+1]);
		check = compare_entries(in1[i+1].e0.e0, in2[i+1].e0.e0);

		if(!check){
		  cout << "\terror occured at " << i << endl;
		}
		*/
	}
	return;
}

void Dummyfield::verify_converted_vectors()
{
	logger.info() << "\tcompare in1_noneo with in1_noneo_converted";
	compare_vectors(sf_in1_noneo, sf_in1_noneo_converted, get_parameters()->get_spinorfieldsize());

	logger.info() << "\tcompare in2_noneo with in2_noneo_converted";
	compare_vectors(sf_in2_noneo, sf_in2_noneo_converted, get_parameters()->get_spinorfieldsize());

	return;
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	cl_context context = opencl_modules[0]->get_context();

	int NUM_ELEMENTS_SF_EO;
	int NUM_ELEMENTS_SF_NON_EO;
	NUM_ELEMENTS_SF_EO =  params.get_eoprec_spinorfieldsize();
	NUM_ELEMENTS_SF_NON_EO =  params.get_spinorfieldsize();

	int NUM_ELEMENTS_AE = params.get_gaugemomentasize() * params.get_su3algebrasize();

	sf_in1_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in2_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in3_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in4_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out_eo = new hmc_float[NUM_ELEMENTS_AE];

	sf_in1_noneo = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_in2_noneo = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_in1_noneo_converted = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_in2_noneo_converted = new spinor[NUM_ELEMENTS_SF_NON_EO];
	sf_out_noneo = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters()->get_use_cg() == true) {
		fill_sf_with_one(sf_in1_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in3_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in2_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in4_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in1_noneo, NUM_ELEMENTS_SF_NON_EO);
		fill_sf_with_one(sf_in2_noneo, NUM_ELEMENTS_SF_NON_EO);
	} else {
		fill_sf_with_random_eo(sf_in1_eo, sf_in2_eo, NUM_ELEMENTS_SF_EO, 123456);
		fill_sf_with_random_eo(sf_in3_eo, sf_in4_eo, NUM_ELEMENTS_SF_EO, 789101);
		fill_sf_with_random_noneo(sf_in1_noneo, NUM_ELEMENTS_SF_NON_EO, 123456, get_parameters());
		fill_sf_with_random_noneo(sf_in2_noneo, NUM_ELEMENTS_SF_NON_EO, 789101, get_parameters());
	}

	fill_eo_sf_with_specific_float(sf_in1_eo, sf_in2_eo, get_parameters());
	fill_eo_sf_with_specific_float(sf_in3_eo, sf_in4_eo, get_parameters());
	fill_noneo_sf_with_specific_float(sf_in1_noneo, get_parameters());
	fill_noneo_sf_with_specific_float(sf_in2_noneo, get_parameters());


	BOOST_REQUIRE(sf_in1_eo);
	BOOST_REQUIRE(sf_in2_eo);
	BOOST_REQUIRE(sf_in3_eo);
	BOOST_REQUIRE(sf_in4_eo);
	BOOST_REQUIRE(sf_in1_noneo);
	BOOST_REQUIRE(sf_in2_noneo);
	BOOST_REQUIRE(sf_in1_noneo_converted);
	BOOST_REQUIRE(sf_in2_noneo_converted);

	fill_with_zero(sf_out_eo, NUM_ELEMENTS_AE);
	fill_with_zero(sf_out_noneo, NUM_ELEMENTS_AE);

	size_t sf_buf_size_eo, sf_buf_size_noneo;
	sf_buf_size_eo = get_parameters()->get_eo_sf_buf_size();
	sf_buf_size_noneo = get_parameters()->get_sf_buf_size();
	size_t ae_buf_size = get_parameters()->get_gm_buf_size();
	//create buffer for sf on device (and copy sf_in to both for convenience)

	in1_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_eo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_eo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in3_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_eo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in4_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_eo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in1_eo, CL_TRUE, 0, sf_buf_size_eo, sf_in1_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in2_eo, CL_TRUE, 0, sf_buf_size_eo, sf_in2_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in3_eo, CL_TRUE, 0, sf_buf_size_eo, sf_in3_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in4_eo, CL_TRUE, 0, sf_buf_size_eo, sf_in4_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	out_eo = clCreateBuffer(context, CL_MEM_WRITE_ONLY, ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	this->reset_outfield_eo();


	//now convert the eo input to noneo input...
	convert_to_noneo(sf_in1_noneo_converted, sf_in1_eo, sf_in2_eo);
	convert_to_noneo(sf_in2_noneo_converted, sf_in3_eo, sf_in4_eo);

	in1_noneo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_noneo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2_noneo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_noneo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in1_noneo_converted = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_noneo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2_noneo_converted = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size_noneo, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in1_noneo, CL_TRUE, 0, sf_buf_size_noneo, sf_in1_noneo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in2_noneo, CL_TRUE, 0, sf_buf_size_noneo, sf_in2_noneo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in1_noneo_converted, CL_TRUE, 0, sf_buf_size_noneo, sf_in1_noneo_converted, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in2_noneo_converted, CL_TRUE, 0, sf_buf_size_noneo, sf_in2_noneo_converted, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	out_noneo = clCreateBuffer(context, CL_MEM_WRITE_ONLY, ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), out_noneo, CL_TRUE, 0, ae_buf_size, sf_out_noneo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	this->reset_outfield_noneo();

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	Opencl_Module_Hmc::fill_kernels();

	ae_sqn = createKernel("gaugemomentum_squarenorm") << basic_fermion_code << "types_hmc.h" << "operations_gaugemomentum.cl" << "gaugemomentum_squarenorm.cl";

	testKernel = createKernel("fermion_force_eoprec") << basic_fermion_code << "types_hmc.h"  << "operations_gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion_eo.cl";

	testKernel2 = createKernel("fermion_force") << basic_fermion_code << "types_hmc.h"  << "operations_gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in1_eo);
	clReleaseMemObject(in2_eo);
	clReleaseMemObject(in3_eo);
	clReleaseMemObject(in4_eo);
	clReleaseMemObject(out_eo);
	clReleaseMemObject(sqnorm);

	delete[] sf_in1_eo;
	delete[] sf_in2_eo;
	delete[] sf_in3_eo;
	delete[] sf_in4_eo;
	delete[] sf_out_eo;

	clReleaseMemObject(in1_noneo);
	clReleaseMemObject(in2_noneo);
	clReleaseMemObject(out_noneo);

	delete[] sf_in1_noneo;
	delete[] sf_in2_noneo;
	delete[] sf_out_noneo;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem in1, cl_mem in2, cl_mem gf, int gs, int ls, int evenodd)
{
	// convert input to SOA. TODO move to proper place, does not have to be done every time
	convertSpinorfieldToSOA_eo_device(spinorfield_soa_eo_1, in1);
	convertSpinorfieldToSOA_eo_device(spinorfield_soa_eo_2, in2);

	cl_int err;
	int eo;
	if(evenodd == ODD)
		eo = ODD;
	else
		eo = EVEN;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &spinorfield_soa_eo_1);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &spinorfield_soa_eo_2);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 3, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 4, sizeof(int), &eo);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

void Device::runTestKernel2(cl_mem out, cl_mem in1, cl_mem in2, cl_mem gf, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel2, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel2, 1, sizeof(cl_mem), &in1);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel2, 2, sizeof(cl_mem), &in2);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel2, 3, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel2, gs, ls);
}

//this is a copy of "set_float_to_gaugemomentum_squarenorm_device"
void Device::set_float_to_gm_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out)
{
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(ae_sqn, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	int clerr = clSetKernelArg(ae_sqn, 0, sizeof(cl_mem), &clmem_in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	//  /** @todo add reduction */
	clerr = clSetKernelArg(ae_sqn,  1, sizeof(cl_mem), &clmem_out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(ae_sqn  , gs2, ls2);
}


hmc_float Dummyfield::get_squarenorm_eo(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in1_eo, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in2_eo, sqnorm);
	if(which == 2) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in3_eo, sqnorm);
	if(which == 3) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in4_eo, sqnorm);
	if(which == 4) static_cast<Device*>(opencl_modules[0])->set_float_to_gm_squarenorm_device(out_eo, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << result;
	return result;
}

hmc_float Dummyfield::get_squarenorm_noneo(int which)
{

	//CP: I only use out_eo here because there is some mistake with out_noneo. However, I will not try to find it...
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in1_noneo, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in2_noneo, sqnorm);
	if(which == 2) {
		cl_mem tmp;
		cl_int err;
		cl_context context = opencl_modules[0]->get_context();
		tmp = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
		static_cast<Device*>(opencl_modules[0])->set_float_to_gm_squarenorm_device(out_eo, tmp);
		hmc_float result;
		err = clEnqueueReadBuffer(*queue, tmp, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
		logger.info() << result;
		return result;
	}
	if(which == 3) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in1_noneo_converted, sqnorm);
	if(which == 4) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in2_noneo_converted, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << result;
	return result;
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

void Dummyfield::runTestKernel(int evenodd)
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_parameters()->get_eoprec_spinorfieldsize();
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	//interprete Y = (in1_eo, in2_eo) X = (in3_eo, in4_eo)
	//Y_odd = in2_eo, Y_even = in1_eo, X_odd = in4_eo, X_even = in3_eo
	if(evenodd == ODD) {
		//this is then force(Y_odd, X_even) == force(in2, in3)
		static_cast<Device*>(opencl_modules[0])->runTestKernel(out_eo, in2_eo, in3_eo, *(get_clmem_gaugefield()), gs, ls, evenodd);
	} else {
		//this is then force(Y_even, X_odd) == force(in1, in4)
		static_cast<Device*>(opencl_modules[0])->runTestKernel(out_eo, in1_eo, in4_eo, *(get_clmem_gaugefield()), gs, ls, evenodd);
	}
}

void Dummyfield::runTestKernel2()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_parameters()->get_spinorfieldsize();
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	//CP: I only use out_eo here because there is some mistake with out_noneo. However, I will not try to find it...
	static_cast<Device*>(opencl_modules[0])->runTestKernel2(out_eo, in1_noneo, in2_noneo, *(get_clmem_gaugefield()), gs, ls);
}

void Dummyfield::runTestKernel2withconvertedfields()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_parameters()->get_spinorfieldsize();
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	static_cast<Device*>(opencl_modules[0])->runTestKernel2(out_noneo, in1_noneo_converted, in2_noneo_converted, *(get_clmem_gaugefield()), gs, ls);
}
