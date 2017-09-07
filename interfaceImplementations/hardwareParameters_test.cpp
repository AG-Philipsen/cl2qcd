/**
 * Copyright 2015 Christopher Pinke
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::parameters
#include <boost/test/unit_test.hpp>

#include "hardwareParameters.hpp"

BOOST_AUTO_TEST_CASE(implementByMeansOfMetaInputparameters)
{
	const char * argv []  = {"foo"};
	meta::Inputparameters fullParameters{1, argv};
	hardware::HardwareParametersImplementation hardwareParameters( &fullParameters );

	BOOST_REQUIRE_EQUAL( hardwareParameters.useGpu(), fullParameters.get_use_gpu() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.useCpu() , fullParameters.get_use_cpu() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.splitCpu() , fullParameters.get_split_cpu() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.getMaximalNumberOfDevices() , fullParameters.get_device_count() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.getSelectedDevices().size() , fullParameters.get_selected_devices().size() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.enableProfiling() , fullParameters.get_enable_profiling());
	BOOST_REQUIRE_EQUAL( hardwareParameters.getNs() , fullParameters.get_nspace() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.getNt() , fullParameters.get_ntime() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.disableOpenCLCompilerOptimizations() , fullParameters.is_ocl_compiler_opt_disabled() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.useSameRandomNumbers() , fullParameters.get_use_same_rnd_numbers() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.useEvenOddPreconditioning() , fullParameters.get_use_eo() );
	BOOST_REQUIRE_EQUAL( hardwareParameters.getSpatialLatticeVolume(), meta::get_volspace( fullParameters ) );
	BOOST_REQUIRE_EQUAL( hardwareParameters.getLatticeVolume(), meta::get_vol4d( fullParameters ) );
}
