// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gaugefield_convert
#include <boost/test/unit_test.hpp>

#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"
#include "../meta/type_ops.hpp"

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, const hardware::System * system) : Gaugefield_hybrid(system) {
		init(1, device_type);
		fill_buffers();
	};

	virtual void finalize_opencl() override;

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
	hardware::System system(params);
	Dummyfield dummy(CL_DEVICE_TYPE_CPU, &system);
	dummy.set_gaugefield_cold(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

BOOST_AUTO_TEST_CASE(CPU_hot)
{
	const char* _params[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	Dummyfield dummy(CL_DEVICE_TYPE_CPU, &system);
	dummy.set_gaugefield_hot(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

BOOST_AUTO_TEST_CASE(GPU_cold)
{
	const char* _params[] = {"foo", "--use_gpu=true"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system);
	dummy.set_gaugefield_cold(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

BOOST_AUTO_TEST_CASE(GPU_hot)
{
	const char* _params[] = {"foo", "--use_gpu=true"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system);
	dummy.set_gaugefield_hot(dummy.in);

	dummy.send();
	dummy.recieve();

	dummy.verify();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

void Dummyfield::fill_buffers()
{
	size_t NUM_ELEMENTS = meta::get_vol4d(get_parameters()) * NDIM;
	in = new Matrixsu3[NUM_ELEMENTS];
	out = new Matrixsu3[NUM_ELEMENTS];
}

void Dummyfield::clear_buffers()
{
	delete[] in;
	delete[] out;
}

void Dummyfield::verify()
{
	size_t NUM_ELEMENTS = meta::get_vol4d(get_parameters()) * NDIM;
	BOOST_CHECK_EQUAL_COLLECTIONS(in, in + NUM_ELEMENTS, out, out + NUM_ELEMENTS);
}

void Dummyfield::send()
{
	static_cast<hardware::code::Gaugefield*>(opencl_modules[0])->importGaugefield(in);
}

void Dummyfield::recieve()
{
	static_cast<hardware::code::Gaugefield*>(opencl_modules[0])->exportGaugefield(out);
}
