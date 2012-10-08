/** @file
 * Testcases for the hardware::buffers::ScalarBuffer template
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "scalar_buffer.hpp"
#include "../system.hpp"
#include "../../types.h"
#include "../../meta/type_ops.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::ScalarBuffer
#include <boost/test/unit_test.hpp>

template<typename T> void test(size_t elems, hardware::Device * device)
{
	using namespace hardware::buffers;

	ScalarBuffer<T> dummy(elems, device);
	BOOST_REQUIRE_EQUAL(dummy.get_elements(), elems);
	BOOST_REQUIRE_EQUAL(dummy.get_bytes(), elems * sizeof(T));
	const cl_mem * tmp = dummy;
	BOOST_REQUIRE(tmp);
	BOOST_REQUIRE(*tmp);
	BOOST_REQUIRE_EQUAL(device->get_id(), dummy.get_device()->get_id());

	T* in = new T[elems];
	T* out = new T[elems];
	fill(in, elems, 1);
	fill(out, elems, 2);
	dummy.load(in);
	dummy.dump(out);
	BOOST_CHECK_EQUAL_COLLECTIONS(in, in + elems, out, out + elems);
}

template<typename T> void test(bool requireDouble = false)
{
	using namespace hardware;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params);
	const std::vector<Device*>& devices = system.get_devices();
for(Device * device : devices) {
		if(!requireDouble || device->is_double_supported()) {
			test<T>(1, device);
			test<T>(1024, device);
		}
	}
}

BOOST_AUTO_TEST_CASE(float_buffer)
{
	test<float>();
}

BOOST_AUTO_TEST_CASE(double_buffer)
{
	test<double>(true);
}

BOOST_AUTO_TEST_CASE(hmc_complex_buffer)
{
	test<hmc_complex>(true);
}

BOOST_AUTO_TEST_CASE(SU3_buffer)
{
	test<Matrixsu3>(true);
}
