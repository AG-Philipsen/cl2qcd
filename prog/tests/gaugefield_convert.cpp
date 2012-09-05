// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gaugefield_convert
#include <boost/test/unit_test.hpp>

#include "../opencl_module.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

extern std::string const version;
std::string const version = "0.1";

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, const meta::Inputparameters& params) : Gaugefield_hybrid(params) {
		init(1, device_type);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	void verify();

	void send();
	void recieve();

	Matrixsu3 * in, * out;
private:
	void fill_buffers();
	void clear_buffers();
};

BOOST_AUTO_TEST_CASE(CPU_cold)
{
	const char* _params[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params(2, _params);
	Dummyfield dummy(CL_DEVICE_TYPE_CPU, params);
	dummy.set_gaugefield_cold(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

BOOST_AUTO_TEST_CASE(CPU_hot)
{
	const char* _params[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params(2, _params);
	Dummyfield dummy(CL_DEVICE_TYPE_CPU, params);
	dummy.set_gaugefield_hot(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

BOOST_AUTO_TEST_CASE(GPU_cold)
{
	const char* _params[] = {"foo", "--use_gpu=true"};
	meta::Inputparameters params(2, _params);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, params);
	dummy.set_gaugefield_cold(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

BOOST_AUTO_TEST_CASE(GPU_hot)
{
	const char* _params[] = {"foo", "--use_gpu=true"};
	meta::Inputparameters params(2, _params);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, params);
	dummy.set_gaugefield_hot(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Opencl_Module(get_parameters());
	opencl_modules[0]->init(queue[0], get_max_compute_units(0), get_double_ext(0), 0);

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

void Dummyfield::fill_buffers()
{
	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();
	in = new Matrixsu3[NUM_ELEMENTS_AE];
	out = new Matrixsu3[NUM_ELEMENTS_AE];
}

void Dummyfield::clear_buffers()
{
	delete[] in;
	delete[] out;
}

void Dummyfield::verify()
{
	for(size_t i = 0; i < meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize(); ++i) {
		BOOST_MESSAGE("Element " << i);
		BOOST_REQUIRE_EQUAL(in[i].e00.re, out[i].e00.re);
		BOOST_REQUIRE_EQUAL(in[i].e01.re, out[i].e01.re);
		BOOST_REQUIRE_EQUAL(in[i].e02.re, out[i].e02.re);
		BOOST_REQUIRE_EQUAL(in[i].e10.re, out[i].e10.re);
		BOOST_REQUIRE_EQUAL(in[i].e11.re, out[i].e11.re);
		BOOST_REQUIRE_EQUAL(in[i].e12.re, out[i].e12.re);
		BOOST_REQUIRE_EQUAL(in[i].e20.re, out[i].e20.re);
		BOOST_REQUIRE_EQUAL(in[i].e21.re, out[i].e21.re);
		BOOST_REQUIRE_EQUAL(in[i].e22.re, out[i].e22.re);
		BOOST_REQUIRE_EQUAL(in[i].e00.im, out[i].e00.im);
		BOOST_REQUIRE_EQUAL(in[i].e01.im, out[i].e01.im);
		BOOST_REQUIRE_EQUAL(in[i].e02.im, out[i].e02.im);
		BOOST_REQUIRE_EQUAL(in[i].e10.im, out[i].e10.im);
		BOOST_REQUIRE_EQUAL(in[i].e11.im, out[i].e11.im);
		BOOST_REQUIRE_EQUAL(in[i].e12.im, out[i].e12.im);
		BOOST_REQUIRE_EQUAL(in[i].e20.im, out[i].e20.im);
		BOOST_REQUIRE_EQUAL(in[i].e21.im, out[i].e21.im);
		BOOST_REQUIRE_EQUAL(in[i].e22.im, out[i].e22.im);
	}
}

void Dummyfield::send()
{
	opencl_modules[0]->importGaugefield(in);
}

void Dummyfield::recieve()
{
	opencl_modules[0]->exportGaugefield(out);
}
