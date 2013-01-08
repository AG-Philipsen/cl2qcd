/** @file
 * Tests of the flavour doublets algorithms
 */

#include "flavour_doublet.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugefield
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"
#include "../../logger.hpp"

void test_correlator(const char* params[], const std::vector<hmc_float>& ps_ref, const std::vector<hmc_float>& sc_ref, const std::vector<hmc_float>& vx_ref, const std::vector<hmc_float>& vy_ref, const std::vector<hmc_float>& vz_ref, const std::vector<hmc_float>& ax_ref, const std::vector<hmc_float>& ay_ref, const std::vector<hmc_float>& az_ref);

BOOST_AUTO_TEST_CASE(point_source_0)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=point", "--corr_dir=0"};

	hmc_float ps_tmp[] = {3.11039351013670e+38, 2.37174000000000e+05, 5.71254000000000e+05, 1.05279000000000e+06, 1.68178200000000e+06, 2.45823000000000e+06, 3.38213400000000e+06, 4.45349400000000e+06};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {2.39261039241285e+38, 2.37084000000000e+05, 5.71164000000000e+05, 1.05270000000000e+06, 1.68169200000000e+06, 2.45814000000000e+06, 3.38204400000000e+06, 4.45340400000000e+06};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = { -2.16172782113784e+17, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = { -4.32345564227568e+17, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}

BOOST_AUTO_TEST_CASE(point_source_3)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=point", "--corr_dir=3"};

	hmc_float ps_tmp[] = {1.55519675506835e+38, 1.65269400000000e+06, 1.80994200000000e+06, 1.97640600000000e+06};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {1.19630519620642e+38, 1.65260400000000e+06, 1.80985200000000e+06, 1.97631600000000e+06};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = { -1.08086391056892e+17, 0, 0, 0};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = { -2.16172782113784e+17, 0, 0, 0};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}


BOOST_AUTO_TEST_CASE(stochastic_source_0)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=volume", "--corr_dir=0"};

	hmc_float ps_tmp[] = {3.11039351013670e+38, 2.37174000000000e+05, 5.71254000000000e+05, 1.05279000000000e+06, 1.68178200000000e+06, 2.45823000000000e+06, 3.38213400000000e+06, 4.45349400000000e+06};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {2.39261039241285e+38, 2.37084000000000e+05, 5.71164000000000e+05, 1.05270000000000e+06, 1.68169200000000e+06, 2.45814000000000e+06, 3.38204400000000e+06, 4.45340400000000e+06};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = { -2.16172782113784e+17, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = { -4.32345564227568e+17, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = {0, 0, 0, 0, 0, 0, 0, 0};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}


BOOST_AUTO_TEST_CASE(stochastic_source_3)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=volume", "--corr_dir=3"};

	hmc_float ps_tmp[] = {4.12780705617089e+78, 4.38658448370506e+46, 4.80395251244701e+46, 5.24578167107861e+46};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {1.19630519620642e+38, 1.65260400000000e+06, 1.80985200000000e+06, 1.97631600000000e+06};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = { -1.08086391056892e+17, 0, 0, 0};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = { -2.16172782113784e+17, 0, 0, 0};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = {0, 0, 0, 0};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}

void check_correlator(std::string which, const std::vector<const physics::lattices::Spinorfield*>& solved, const std::vector<const physics::lattices::Spinorfield*>& sources, const meta::Inputparameters& params, const std::vector<hmc_float>& ref)
{
	using namespace std;

	auto result = physics::algorithms::calculate_correlator(which, solved, sources, params);
	logger.debug() << which;
for(auto val: result) {
		logger.debug() << scientific << setprecision(14) << val;
	}
	BOOST_REQUIRE_EQUAL(result.size(), ref.size());
	for(size_t i = 0; i < result.size(); ++i) {
		BOOST_REQUIRE_CLOSE(result[i], ref[i], 0.1);
	}
}

void test_correlator(const char* _params[], const std::vector<hmc_float>& ps_ref, const std::vector<hmc_float>& sc_ref, const std::vector<hmc_float>& vx_ref, const std::vector<hmc_float>& vy_ref, const std::vector<hmc_float>& vz_ref, const std::vector<hmc_float>& ax_ref, const std::vector<hmc_float>& ay_ref, const std::vector<hmc_float>& az_ref)
{
	using namespace physics::lattices;

	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	size_t num_sources = params.get_num_sources();
	auto sources = create_spinorfields(system, num_sources);
	for(size_t i = 0; i < num_sources; ++i) {
		pseudo_randomize<Spinorfield, spinor>(sources[i], i);
	}
	auto solved = create_spinorfields(system, num_sources);
	for(size_t i = 0; i < num_sources; ++i) {
		pseudo_randomize<Spinorfield, spinor>(solved[i], i + num_sources);
	}

	check_correlator("ps", solved, sources, params, ps_ref);
	check_correlator("sc", solved, sources, params, sc_ref);
	check_correlator("vx", solved, sources, params, vx_ref);
	check_correlator("vy", solved, sources, params, vy_ref);
	check_correlator("vz", solved, sources, params, vz_ref);
	check_correlator("ax", solved, sources, params, ax_ref);
	check_correlator("ay", solved, sources, params, ay_ref);
	check_correlator("az", solved, sources, params, az_ref);
}
