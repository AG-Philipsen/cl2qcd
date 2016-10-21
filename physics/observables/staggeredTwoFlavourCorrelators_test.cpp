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


void compare_staggered_correlator(std::string which, const std::pair<std::vector<physics::lattices::Staggeredfield_eo*>, std::vector<physics::lattices::Staggeredfield_eo*> >& invertedSources, const std::vector<hmc_float>& ref, physics::InterfacesHandler & interfacesHandler)
{
	using namespace std;

	std::vector<hmc_float> correlator = physics::observables::staggered::calculatePseudoscalarCorrelator(invertedSources, interfacesHandler);
	logger.debug() << which;
	for(auto tmp_corr: correlator)
	{
		logger.debug() << scientific << setprecision(14) << tmp_corr;
	}

	BOOST_REQUIRE_EQUAL(correlator.size(), ref.size());
	for(size_t i = 0; i < correlator.size(); ++i)
	{
		BOOST_CHECK_CLOSE(correlator[i], ref[i], 0.1);
	}
}


void test_staggered_correlator(const char* _params[], const std::vector<hmc_float>& ps_ref, bool coldOrZero)
{
	using namespace physics::lattices;

	meta::Inputparameters params(6, _params);
    hardware::HardwareParametersImplementation hP(&params); //	Ns = 4, Nt = 8.
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	const physics::PRNG prng{system, &prngParameters};

	size_t num_sources = params.get_num_sources();
	std::vector<physics::lattices::Staggeredfield_eo *> invertedSourcesEven = create_staggeredfields_eo(system, num_sources, interfacesHandler);
	std::vector<physics::lattices::Staggeredfield_eo *> invertedSourcesOdd = create_staggeredfields_eo(system, num_sources, interfacesHandler);

	if(coldOrZero)
	{
		for(size_t i = 0; i < num_sources; ++i)
		{
			invertedSourcesEven[i]->set_cold();
			invertedSourcesOdd[i]->set_cold();
		}
	}
	else
	{
		for(size_t i = 0; i < num_sources; ++i)
		{
			invertedSourcesEven[i]->set_zero();
			invertedSourcesOdd[i]->set_zero();
		}

	}

	auto invertedSources = std::make_pair(invertedSourcesEven, invertedSourcesOdd);

	compare_staggered_correlator("ps", invertedSources, ps_ref, interfacesHandler);

	release_staggeredfields_eo(invertedSourcesEven);
	release_staggeredfields_eo(invertedSourcesOdd);
}

BOOST_AUTO_TEST_CASE(point_source)
{
	using namespace physics::lattices;

	const char * params[] = {"foo", "--sourcetype=point", "--num_sources=3", "--fermact=rooted_stagg", "--corr_dir=0", "--num_dev=1"};

	hmc_float referenceValuePerTimesliceZero = 0.;
	std::vector<hmc_float> referenceValueZero(8, referenceValuePerTimesliceZero);

	unsigned int Nsources = 3;
	unsigned int Nt = 8;
	unsigned int Ns = 4;
	hmc_float referenceValuePerTimesliceCold = 1./(Ns*Ns*Ns*Nt)*Nsources;
	std::vector<hmc_float> referenceValueCold(Nt, referenceValuePerTimesliceCold);

	test_staggered_correlator(params, referenceValueZero, false);
	test_staggered_correlator(params, referenceValueCold, true);
}


