// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gaugefield_observables
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../opencl_module.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

extern std::string const version;
std::string const version = "0.1";

const std::string SOURCEFILE = std::string(SOURCEDIR)
#ifdef _USEDOUBLEPREC_
                               + "/tests/f_gauge_input_1";
#else
                               + "/tests/f_gauge_input_1_single";
#endif
const char * PARAMS[] = {"foo", SOURCEFILE.c_str()};
const meta::Inputparameters INPUT(2, PARAMS);

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type) : Gaugefield_hybrid(INPUT) {
		init(1, device_type);
	};

	Opencl_Module * getDevice() const;
};

//gaugeobservables(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol);

BOOST_AUTO_TEST_CASE(CPU_cold)
{
	Dummyfield dummy(CL_DEVICE_TYPE_CPU);
	dummy.set_gaugefield_cold(dummy.get_sgf());
	dummy.copy_gaugefield_to_all_tasks();

	// get reference solutions
	hmc_float ref_plaq, ref_tplaq, ref_splaq;
	ref_plaq = dummy.plaquette(&ref_tplaq, &ref_splaq);
	hmc_complex ref_pol = dummy.polyakov();

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	dummy.getDevice()->gaugeobservables(&dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	// verify
	BOOST_REQUIRE_CLOSE(ref_plaq,   dev_plaq,     0.1f);
	BOOST_REQUIRE_CLOSE(ref_tplaq,  dev_tplaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_splaq,  dev_splaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.re, dev_pol.re,   0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.im, dev_pol.im,   0.1f);
}

BOOST_AUTO_TEST_CASE(CPU_hot)
{
	Dummyfield dummy(CL_DEVICE_TYPE_CPU);
	dummy.set_gaugefield_hot(dummy.get_sgf());
	dummy.copy_gaugefield_to_all_tasks();

	// get reference solutions
	hmc_float ref_plaq, ref_tplaq, ref_splaq;
	ref_plaq = dummy.plaquette(&ref_tplaq, &ref_splaq);
	hmc_complex ref_pol = dummy.polyakov();

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	dummy.getDevice()->gaugeobservables(&dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	// verify
	BOOST_REQUIRE_CLOSE(ref_plaq,   dev_plaq,     0.1f);
	BOOST_REQUIRE_CLOSE(ref_tplaq,  dev_tplaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_splaq,  dev_splaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.re, dev_pol.re,   0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.im, dev_pol.im,   0.1f);
}

BOOST_AUTO_TEST_CASE(GPU_cold)
{
	Dummyfield dummy(CL_DEVICE_TYPE_GPU);
	dummy.set_gaugefield_cold(dummy.get_sgf());
	dummy.copy_gaugefield_to_all_tasks();

	// get reference solutions
	hmc_float ref_plaq, ref_tplaq, ref_splaq;
	ref_plaq = dummy.plaquette(&ref_tplaq, &ref_splaq);
	hmc_complex ref_pol = dummy.polyakov();

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	dummy.getDevice()->gaugeobservables(&dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	// verify
	BOOST_REQUIRE_CLOSE(ref_plaq,   dev_plaq,     0.1f);
	BOOST_REQUIRE_CLOSE(ref_tplaq,  dev_tplaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_splaq,  dev_splaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.re, dev_pol.re,   0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.im, dev_pol.im,   0.1f);
}

BOOST_AUTO_TEST_CASE(GPU_hot)
{
	Dummyfield dummy(CL_DEVICE_TYPE_GPU);
	dummy.set_gaugefield_hot(dummy.get_sgf());
	dummy.copy_gaugefield_to_all_tasks();

	// get reference solutions
	hmc_float ref_plaq, ref_tplaq, ref_splaq;
	ref_plaq = dummy.plaquette(&ref_tplaq, &ref_splaq);
	hmc_complex ref_pol = dummy.polyakov();

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	dummy.getDevice()->gaugeobservables(&dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	// verify
	BOOST_REQUIRE_CLOSE(ref_plaq,   dev_plaq,     0.1f);
	BOOST_REQUIRE_CLOSE(ref_tplaq,  dev_tplaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_splaq,  dev_splaq,    0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.re, dev_pol.re,   0.1f);
	BOOST_REQUIRE_CLOSE(ref_pol.im, dev_pol.im,   0.1f);
}

Opencl_Module * Dummyfield::getDevice() const
{
	return opencl_modules[0];
}
