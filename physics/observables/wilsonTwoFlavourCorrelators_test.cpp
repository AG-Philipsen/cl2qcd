/** @file
 * Unit tests for the physics::observables::wilson::TwoFlavourCorrelators class
 *
 * Copyright 2014,Christopher Pinke
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "wilsonTwoFlavourCorrelators.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::observables::wilson::TwoFlavourCorrelators
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include <stdexcept>
#include "../../host_functionality/logger.hpp"
#include "../lattices/util.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

void test_correlator(const char* params[], const std::vector<hmc_float>& ps_ref, const std::vector<hmc_float>& sc_ref, const std::vector<hmc_float>& vx_ref, const std::vector<hmc_float>& vy_ref, const std::vector<hmc_float>& vz_ref, const std::vector<hmc_float>& ax_ref, const std::vector<hmc_float>& ay_ref, const std::vector<hmc_float>& az_ref);

BOOST_AUTO_TEST_CASE(point_source_0)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=point", "--corr_dir=0"};

	hmc_float ps_tmp[] = {5.9989365191736104, 5.9685966494875613, 5.9273621601759441, 5.9971283038020733, 5.9813369771401304, 5.9760321595370245, 6.0879077700116371, 6.0659567780266324};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0.0024612863844851191, -0.016374834112823547, -0.043961011270704092, -0.01523594115058316, -0.005165987427987663, 0.024049710560252593, 0.061298733087266286, -0.045272255628473722};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {4.5356982259415997, 4.4780444134285631, 4.4472419503312492, 4.5430446606989898, 4.5084478219269606, 4.4981626419407812, 4.5664936387672315, 4.5587559044286214};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = { -0.01289411252205299, -0.03919045962063928, -0.010173642800255882, 0.022418793182812542, -0.045457750951185025, 0.024744535419928094, 0.022423311850699623, 0.028991528724414503};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0.010392724423067302, -0.0092227981570463906, -0.0019124537869479354, 0.0015147185464041528, -0.039742353575948297, 0.034283909445848255, 0.015519562913546912, 0.031784190562001421};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = { 0.031118718400293559, -0.021908426433681693, -0.041120764054587693, -0.0372499095141723, 0.0091808477678250507, 0.039072907858853943, 0.06206377296194393, -0.045106698264689667};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = { -0.029066752308006504, 0.05616977588026148, -0.027229003286407247, -0.00050158467927097217, -0.0066461981729357446, 0.019173477490037887, -0.033509624999394916, -0.080189233267604784};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = { -0.047639872928628396, 0.093974382841054555, -0.043026595394032315, -0.019346309887679554, -0.00026516772944583378, 0.020659821224574604, -0.038190332787656045, -0.071688867318142765};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}

BOOST_AUTO_TEST_CASE(point_source_3)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=point", "--corr_dir=3"};

	hmc_float ps_tmp[] = {6.0056739949017723, 6.0049719331029232, 5.9919404833842922, 5.9990422472883189};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0.040695238631610849, -0.015782960265070797, -0.0039547471329859883, -0.040057681012838249};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {4.529683819089029, 4.5264627913908377, 4.5043342416499934, 4.5074637766021466};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = {0.0027656409216830317, 0.0038867393439130579, -0.019706156293010021, 0.0084848776692746999};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0.0077088724709378151, 0.00074884916511421288, -0.0098420917167491295, 0.022693120266159838};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = {0.047407708627848946, -0.0019036665046849276, 0.0049039341798998871, -0.052382751942171284};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = {0.0087981502113759005, -0.02035897271408698, -0.012418789394657199, -0.026919959774292151};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = {0.0098842213365260417, -0.021517358215967168, -0.03372254795823848, -0.0074057861522982219};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}


BOOST_AUTO_TEST_CASE(stochastic_source_0)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=volume", "--corr_dir=0"};

	hmc_float ps_tmp[] = {5.998936519, 5.9685966494875613, 5.9273621601759441, 5.9971283038020733, 5.9813369771401304, 5.9760321595370245, 6.0879077700116371, 6.0659567780266324};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0.0024612863844851191, -0.016374834112823547, -0.043961011270704092, -0.01523594115058316, -0.005165987427987663, 0.024049710560252593, 0.061298733087266286, -0.045272255628473722};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {4.5356982259415997, 4.4780444134285631, 4.4472419503312492, 4.5430446606989898, 4.5084478219269606, 4.4981626419407812, 4.5664936387672315, 4.5587559044286214};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = { -0.01289411252205299, -0.03919045962063928, -0.010173642800255882, 0.022418793182812542, -0.045457750951185025, 0.024744535419928094, 0.022423311850699623, 0.028991528724414503};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0.010392724423067302, -0.0092227981570463906, -0.0019124537869479354, 0.0015147185464041528, -0.039742353575948297, 0.034283909445848255, 0.015519562913546912, 0.031784190562001421};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = {0.031118718400293559, -0.021908426433681693, -0.041120764054587693, -0.0372499095141723, 0.0091808477678250507, 0.039072907858853943, 0.06206377296194393, -0.045106698264689667};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = { -0.029066752308006504, 0.05616977588026148, -0.027229003286407247, -0.00050158467927097217, -0.0066461981729357446, 0.019173477490037887, -0.033509624999394916, -0.080189233267604784};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = { -0.047639872928628396, 0.093974382841054555, -0.043026595394032315, -0.019346309887679554, -0.00026516772944583378, 0.020659821224574604, -0.038190332787656045, -0.071688867318142765};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}


BOOST_AUTO_TEST_CASE(stochastic_source_3)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=volume", "--corr_dir=3"};

	hmc_float ps_tmp[] = {24595.308692803472, 24592.45084, 24539.27566, 24568.2908};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));
	hmc_float sc_tmp[] = {0.040695238631610849, -0.015782960265070797, -0.0039547471329859883, -0.040057681012838249};
	std::vector<hmc_float> sc(sc_tmp, sc_tmp + sizeof(sc_tmp) / sizeof(hmc_float));
	hmc_float vx_tmp[] = {4.529683819089029, 4.5264627913908377, 4.5043342416499934, 4.5074637766021466};
	std::vector<hmc_float> vx(vx_tmp, vx_tmp + sizeof(vx_tmp) / sizeof(hmc_float));
	hmc_float vy_tmp[] = {0.0027656409216830317, 0.0038867393439130579, -0.019706156293010021, 0.0084848776692746999};
	std::vector<hmc_float> vy(vy_tmp, vy_tmp + sizeof(vy_tmp) / sizeof(hmc_float));
	hmc_float vz_tmp[] = {0.0077088724709378151, 0.00074884916511421288, -0.0098420917167491295, 0.022693120266159838};
	std::vector<hmc_float> vz(vz_tmp, vz_tmp + sizeof(vz_tmp) / sizeof(hmc_float));
	hmc_float ax_tmp[] = {0.047407708627848946, -0.0019036665046849276, 0.0049039341798998871, -0.052382751942171284};
	std::vector<hmc_float> ax(ax_tmp, ax_tmp + sizeof(ax_tmp) / sizeof(hmc_float));
	hmc_float ay_tmp[] = {0.0087981502113759005, -0.02035897271408698, -0.012418789394657199, -0.026919959774292151};
	std::vector<hmc_float> ay(ay_tmp, ay_tmp + sizeof(ay_tmp) / sizeof(hmc_float));
	hmc_float az_tmp[] = {0.0098842213365260417, -0.021517358215967168, -0.03372254795823848, -0.0074057861522982219};
	std::vector<hmc_float> az(az_tmp, az_tmp + sizeof(az_tmp) / sizeof(hmc_float));

	test_correlator(params, ps, sc, vx, vy, vz, ax, ay, az);
}

void check_correlator(std::string which, const std::vector<physics::lattices::Spinorfield*>& solved, const std::vector<physics::lattices::Spinorfield*>& sources, const hardware::System& system, const std::vector<hmc_float>& ref, physics::InterfacesHandler & interfacesHandler)
{
	using namespace std;

	auto result = physics::observables::wilson::calculate_correlator(which, solved, sources, system, interfacesHandler);
	logger.debug() << which;
for(auto val: result) {
		logger.debug() << scientific << setprecision(14) << val;
	}
	BOOST_REQUIRE_EQUAL(result.size(), ref.size());
	for(size_t i = 0; i < result.size(); ++i) {
		BOOST_CHECK_CLOSE(result[i], ref[i], 0.1);
	}
}

void test_correlator(const char* _params[], const std::vector<hmc_float>& ps_ref, const std::vector<hmc_float>& sc_ref, const std::vector<hmc_float>& vx_ref, const std::vector<hmc_float>& vy_ref, const std::vector<hmc_float>& vz_ref, const std::vector<hmc_float>& ax_ref, const std::vector<hmc_float>& ay_ref, const std::vector<hmc_float>& az_ref)
{
	using namespace physics::lattices;

	meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	const physics::PRNG prng{system, &prngParameters};

	size_t num_sources = params.get_num_sources();
	auto sources = create_spinorfields(system, num_sources, interfacesHandler);
	for(size_t i = 0; i < num_sources; ++i) {
		pseudo_randomize<Spinorfield, spinor>(sources[i], i);
	}
	auto solved = create_spinorfields(system, num_sources, interfacesHandler);
	for(size_t i = 0; i < num_sources; ++i) {
		pseudo_randomize<Spinorfield, spinor>(solved[i], i + num_sources);
	}

for(auto source: sources) {
		log_squarenorm("Source: ", *source);
	}

	check_correlator("ps", solved, sources, system, ps_ref, interfacesHandler);
	check_correlator("sc", solved, sources, system, sc_ref, interfacesHandler);
	check_correlator("vx", solved, sources, system, vx_ref, interfacesHandler);
	check_correlator("vy", solved, sources, system, vy_ref, interfacesHandler);
	check_correlator("vz", solved, sources, system, vz_ref, interfacesHandler);
	check_correlator("ax", solved, sources, system, ax_ref, interfacesHandler);
	check_correlator("ay", solved, sources, system, ay_ref, interfacesHandler);
	check_correlator("az", solved, sources, system, az_ref, interfacesHandler);

	release_spinorfields(solved);
	release_spinorfields(sources);
}
