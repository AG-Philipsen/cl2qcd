/** @file
 * Tests for functions working with sources.
 *
 * Copyright (c) 2012,2013 Matthias Bach
 * Copyright (c) 2014,2015 Christopher Pinke
 * Copyright (c) 2014-2016,2018 Alessandro Sciarra
 * Copyright (c) 2015 Christopher Czaban
 * Copyright (c) 2016 Tim Breitenfelder
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics
#include "sources.hpp"

#include "../interfaceImplementations/hardwareParameters.hpp"
#include "../interfaceImplementations/interfacesHandler.hpp"
#include "../interfaceImplementations/openClKernelParameters.hpp"
#include "test_util_staggered.hpp"

#include <boost/test/unit_test.hpp>
#include <sstream>

static void test_sources(std::string type, int num_sources)
{
    using namespace physics::lattices;

    std::stringstream tmp;
    tmp << "--nSources=";
    tmp << num_sources;
    std::string n_sources_string  = tmp.str();
    std::string sourcetype_string = std::string("--sourceType=") + type;
    const char* _params[]         = {"foo", n_sources_string.c_str(), sourcetype_string.c_str()};
    meta::Inputparameters params(3, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
    physics::InterfacesHandlerImplementation interfacesHandler{params};
    physics::PrngParametersImplementation prngParameters{params};
    physics::PRNG prng{system, &prngParameters};

    auto sources = create_sources(system, prng, params.get_num_sources(), interfacesHandler);

    BOOST_REQUIRE_EQUAL(params.get_num_sources(), static_cast<const int>(sources.size()));

    release_spinorfields(sources);
}

static void test_staggered_sources(std::string type, int num_sources)
{
    using namespace physics::lattices;

    std::stringstream tmp;
    tmp << "--nSources=";
    tmp << num_sources;
    std::string n_sources_string  = tmp.str();
    std::string sourcetype_string = std::string("--sourceType=") + type;
    const char* _params[]         = {"foo", n_sources_string.c_str(), sourcetype_string.c_str(),
                             "--fermionAction=rooted_stagg"};
    meta::Inputparameters params(4, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
    physics::InterfacesHandlerImplementation interfacesHandler{params};
    physics::PrngParametersImplementation prngParameters{params};
    physics::PRNG prng{system, &prngParameters};

    auto staggered_sources = create_staggered_sources(system, prng, params.get_num_sources(), interfacesHandler);

    BOOST_REQUIRE_EQUAL(params.get_num_sources(), static_cast<const int>(staggered_sources.size()));

    release_staggeredfields_eo(staggered_sources);
}

static void test_volume_source_stagg(std::string content)
{
    using namespace physics::lattices;

    std::vector<const char*> options(1, "foo");
    options.push_back("--nSpace=8");
    options.push_back("--fermionAction=rooted_stagg");
    options.push_back("--sourceType=volume");
    std::string tmp = "--sourceContent=" + content;
    options.push_back(tmp.c_str());

    meta::Inputparameters params(5, &(options[0]));
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
    physics::InterfacesHandlerImplementation interfacesHandler{params};
    physics::PrngParametersImplementation prngParameters{params};
    physics::PRNG prng{system, &prngParameters};

    Staggeredfield_eo source(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
    set_volume_source(&source, prng);

    // The following lines are to be used to produce the ref_vec file needed to get the ref_value
    //---> Comment them out when the reference values have been obtained!
    /*
    if(content == "gaussian"){
      print_staggeredfield_eo_to_textfile("ref_vec_even", &source, system);
      logger.info() << "Produced the ref_vec_even text file with the staggered field for the ref. code. Returning...";
    }
    return;
    // */

    hmc_float sqnorm = squarenorm(source);
    // Lattice is here 8^4 and we have even-odd preconditioning
    logger.info() << "source content: " << content << " and squarnorm of volume source is " << sqnorm;
    if (content == "one")
        BOOST_CHECK_CLOSE(sqnorm, 6144, 1.e-8);  // Analytic result
    else if (content == "z4")
        BOOST_CHECK_CLOSE(sqnorm, 6144, 1.e-8);  // Analytic result
    else if (content == "gaussian")
        BOOST_CHECK_CLOSE(sqnorm, 6194.3961400489661173, 3);  // Depends on random numbers -> tolerance 3%
    else if (content == "z2")
        BOOST_CHECK_CLOSE(sqnorm, 6144, 1.e-8);  // Analytic result
}

BOOST_AUTO_TEST_CASE(sources)
{
    test_sources("point", 15);
    test_sources("volume", 2);
    test_sources("timeslice", 3);
    test_sources("zslice", 1);
}

BOOST_AUTO_TEST_CASE(sources_staggered)
{
    test_staggered_sources("point", 15);
    test_staggered_sources("volume", 2);
}

BOOST_AUTO_TEST_CASE(volume_source_stagg)
{
    test_volume_source_stagg("one");
    test_volume_source_stagg("z4");
    test_volume_source_stagg("gaussian");
    test_volume_source_stagg("z2");
}
