/* StaggeredTwoFlavourCorrelators_test.cpp
 *
 *@file
 * Implementation of the staggered two flavour correlator calculation.
 *
 * Copyright 2016 Alessandro Sciarra, Tim Breitenfelder
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

#include "staggeredTwoFlavourCorrelators.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::observables::staggered::TwoFlavourCorrelators
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include "../lattices/staggeredfield_eo.hpp"
#include <stdexcept>
#include "../../host_functionality/logger.hpp"
#include "../lattices/util.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

void test_staggered_correlator(const char* params[], const std::vector<hmc_float>& ps_ref);

BOOST_AUTO_TEST_CASE(point_source)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=point", "--corr_dir=0"};

	hmc_float ps_tmp[] = {5.9989365191736104, 5.9685966494875613, 5.9273621601759441, 5.9971283038020733, 5.9813369771401304, 5.9760321595370245, 6.0879077700116371, 6.0659567780266324};
	std::vector<hmc_float> ps(ps_tmp, ps_tmp + sizeof(ps_tmp) / sizeof(hmc_float));

	test_staggered_correlator(params, ps);
}

void check_staggered_correlator(std::string which, const std::pair<std::vector<physics::lattices::Staggeredfield_eo*>, std::vector<physics::lattices::Staggeredfield_eo*> >& invertedSources, const hardware::System& system, const std::vector<hmc_float>& ref, physics::InterfacesHandler & interfacesHandler)
{
	using namespace std;

	std::vector<hmc_float> result = physics::observables::staggered::calculatePseudoscalarCorrelator(invertedSources, interfacesHandler);
	logger.debug() << which;
	for(auto val: result) {
		logger.debug() << scientific << setprecision(14) << val;
	}
	BOOST_REQUIRE_EQUAL(result.size(), ref.size());
	for(size_t i = 0; i < result.size(); ++i) {
		BOOST_CHECK_CLOSE(result[i], ref[i], 0.1);
	}
}

void test_staggered_correlator(const char* _params[], const std::vector<hmc_float>& ps_ref)
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
	std::vector<physics::lattices::Staggeredfield_eo *> sources = create_staggeredfields_eo(system, num_sources, interfacesHandler);
	for(size_t i = 0; i < num_sources; ++i) {
		pseudo_randomize<Staggeredfield_eo, su3vec>(sources[i], i);
	}
	std::vector<physics::lattices::Staggeredfield_eo *> invertedSourcesEven = create_staggeredfields_eo(system, num_sources, interfacesHandler);
	for(size_t i = 0; i < num_sources; ++i) {
		invertedSourcesEven[i]->set_cold();
	}

	std::vector<physics::lattices::Staggeredfield_eo *> invertedSourcesOdd = create_staggeredfields_eo(system, num_sources, interfacesHandler);
	for(size_t i = 0; i < num_sources; ++i) {
		invertedSourcesOdd[i]->set_cold();
	}

	auto invertedSources = std::make_pair(invertedSourcesEven, invertedSourcesOdd);


	for(auto source: sources) {
		log_squarenorm("Source: ", *source);
	}

	check_staggered_correlator("ps", invertedSources, system, ps_ref, interfacesHandler);

	release_staggeredfields_eo(invertedSourcesEven);
	release_staggeredfields_eo(invertedSourcesOdd);
	release_staggeredfields_eo(sources);
}
