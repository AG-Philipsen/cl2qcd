// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gaugefield_convert
#include <boost/test/unit_test.hpp>

#include "../meta/util.hpp"
#include "../meta/type_ops.hpp"
#include "../hardware/device.hpp"
#include "../hardware/system.hpp"
#include "../hardware/code/gaugefield.hpp"

void test(const hardware::System& system, const int seed)
{
	auto params = system.get_inputparameters();
	const size_t NUM_ELEMENTS = meta::get_vol4d(params) * NDIM;
for(auto device: system.get_devices()) {
		Matrixsu3 * const in = new Matrixsu3[NUM_ELEMENTS];
		Matrixsu3 * const out = new Matrixsu3[NUM_ELEMENTS];
		fill(in, NUM_ELEMENTS, seed);
		fill(out, NUM_ELEMENTS, seed + NUM_ELEMENTS);
		hardware::buffers::SU3 buffer(NUM_ELEMENTS, device);

		auto code = device->get_gaugefield_code();
		code->importGaugefield(&buffer, in);
		code->exportGaugefield(out, &buffer);

		BOOST_CHECK_EQUAL_COLLECTIONS(in, in + NUM_ELEMENTS, out, out + NUM_ELEMENTS);

		delete[] in;
		delete[] out;
	}
}

BOOST_AUTO_TEST_CASE(CPU)
{
	const char* _params[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);

	test(system, 1);
	test(system, 14);
	test(system, 21);
}

BOOST_AUTO_TEST_CASE(GPU)
{
	const char* _params[] = {"foo", "--use_cpu=false"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);

	test(system, 1);
	test(system, 14);
	test(system, 21);
}
