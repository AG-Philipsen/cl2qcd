#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Memory Access using Complex
#include <boost/test/unit_test.hpp>

class Device : public hardware::code::Opencl_Module {

	cl_kernel fillComplex;
	cl_kernel readComplex;

	void fill_kernels();
	void clear_kernels();
protected:
	virtual size_t get_read_write_size(const std::string&) const {
		return 0;
	};
	virtual uint64_t get_flop_size(const std::string&) const {
		return 0;
	};

public:
	Device(const meta::Inputparameters& params, hardware::Device * device)
		: Opencl_Module(params, device) {
		fill_kernels();
	};
	~Device() {
		clear_kernels();
	};

	void runFillKernel(const hardware::buffers::Plain<hmc_complex> * out, hmc_complex value);
	void runReadKernel(const hardware::buffers::Plain<cl_float2> * out, const hardware::buffers::Plain<hmc_complex> * in);
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, const hardware::System * system)
		: Gaugefield_hybrid(system), prng(*system) {
		init(1, device_type, prng);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	void verifyFill(const hmc_complex value);
	void runFillKernel(const hmc_complex value);
	void verifyRead();
	void runReadKernel();

private:
	void verify(cl_float2, cl_float2);
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	void clear_buffers();
	hmc_complex * h_complex;
	hardware::buffers::Plain<hmc_complex> * d_complex;
	hardware::buffers::Plain<hmc_complex> * d_readComplex;
	hardware::buffers::Plain<cl_float2> * d_float2;
	physics::PRNG prng;
};

BOOST_AUTO_TEST_CASE( CPU )
{
	const char* _params_cpu[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params_cpu(2, _params_cpu);
	hardware::System system_cpu(params_cpu);
	Dummyfield dummy(CL_DEVICE_TYPE_CPU, &system_cpu);
	dummy.runFillKernel(hmc_complex_zero);
	dummy.verifyFill(hmc_complex_zero);
	dummy.runFillKernel(hmc_complex_one);
	dummy.verifyFill(hmc_complex_one);
	BOOST_MESSAGE("Tested CPU");
}

BOOST_AUTO_TEST_CASE( GPU )
{
	const char* _params_gpu[] = {"foo", "--use_gpu=true"};
	meta::Inputparameters params_gpu(2, _params_gpu);
	hardware::System system_gpu(params_gpu);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system_gpu);
	dummy.runFillKernel(hmc_complex_zero);
	dummy.verifyFill(hmc_complex_zero);
	dummy.runFillKernel(hmc_complex_one);
	dummy.verifyFill(hmc_complex_one);
	BOOST_MESSAGE("Tested GPU");
}

void Dummyfield::init_tasks()
{
	opencl_modules = new hardware::code::Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(get_parameters(), get_device_for_task(0));

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	h_complex = new hmc_complex[meta::get_vol4d(get_parameters())];
	BOOST_REQUIRE(h_complex);

	for(int i = 0; i < meta::get_vol4d(get_parameters()); ++i) {
		hmc_complex tmp = { (hmc_float) i, (hmc_float) (meta::get_vol4d(get_parameters()) - i) };
		h_complex[i] = tmp;
	}

	d_readComplex = new hardware::buffers::Plain<hmc_complex>(meta::get_vol4d(get_parameters()), opencl_modules[0]->get_device());
	d_readComplex->load(h_complex);
	d_float2 = new hardware::buffers::Plain<cl_float2>(meta::get_vol4d(get_parameters()), opencl_modules[0]->get_device());
	d_complex = new hardware::buffers::Plain<hmc_complex>(meta::get_vol4d(get_parameters()), opencl_modules[0]->get_device());
}

void Device::fill_kernels()
{
	fillComplex = createKernel("fillComplex") << get_device()->get_gaugefield_code()->get_sources() << "tests/complex_mem_access.cl";
	readComplex = createKernel("readComplex") << get_device()->get_gaugefield_code()->get_sources() << "tests/complex_mem_access.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	delete d_readComplex;
	delete d_float2;
	delete d_complex;

	delete[] h_complex;
}

void Device::clear_kernels()
{
	clReleaseKernel(fillComplex);
}

void Device::runFillKernel(const hardware::buffers::Plain<hmc_complex> * out, hmc_complex value)
{
	cl_int err;
	err = clSetKernelArg(fillComplex, 0, sizeof(cl_mem), out->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(fillComplex, 1, sizeof(hmc_complex), &value);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	get_device()->enqueue_kernel(fillComplex, 1024);
}

void Device::runReadKernel(const hardware::buffers::Plain<cl_float2> * out, const hardware::buffers::Plain<hmc_complex> * in)
{
	cl_int err;
	err = clSetKernelArg(readComplex, 0, sizeof(cl_mem), out->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(readComplex, 1, sizeof(cl_mem), in->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	get_device()->enqueue_kernel(readComplex, 1024);
}

void Dummyfield::verify(hmc_complex left, hmc_complex right)
{
	BOOST_REQUIRE_EQUAL(left.re, right.re);
	BOOST_REQUIRE_EQUAL(left.im, right.im);
}


void Dummyfield::verifyFill(const hmc_complex value)
{
	// get stuff from device
	d_complex->dump(h_complex);

	for(int i = 0; i < meta::get_vol4d(get_parameters()); ++i) {
		verify(value, h_complex[i]);
	}
}

void Dummyfield::runFillKernel(const hmc_complex value)
{
	static_cast<Device*>(opencl_modules[0])->runFillKernel(d_complex, value);
}

// usually one shouldn't verify flaots in this way, but as we actually only store ints ...
void Dummyfield::verify(cl_float2 left, cl_float2 right)
{
	BOOST_REQUIRE_EQUAL(left.s[0], right.s[0]);
	BOOST_REQUIRE_EQUAL(left.s[1], right.s[1]);
}

void Dummyfield::verifyRead()
{
	// get stuff from device
	cl_float2 * h_float2 = new cl_float2[meta::get_vol4d(get_parameters())];
	d_float2->dump(h_float2);

	for(int i = 0; i < meta::get_vol4d(get_parameters()); ++i) {
		cl_float2 ref = {(cl_float) i, (cl_float) (meta::get_vol4d(get_parameters()) - i)};
		verify(ref, h_float2[i]);
	}

	delete[] h_float2;
}

void Dummyfield::runReadKernel()
{
	static_cast<Device*>(opencl_modules[0])->runReadKernel(d_float2, d_readComplex);
}

