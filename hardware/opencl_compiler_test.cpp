/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "opencl_compiler.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OpenCL Compiler
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( PackageCreateEmpty )
{
	ClSourcePackage empty;

	BOOST_REQUIRE_EQUAL(empty.getFiles().size(), 0);
}

BOOST_AUTO_TEST_CASE( PackageAddFiles )
{
	std::string f1 = "foo", f2 = "bar"; // define pointers for later comparison
	ClSourcePackage package = ClSourcePackage() << f1;

	BOOST_REQUIRE_EQUAL(package.getFiles().size(), 1);
	BOOST_REQUIRE_EQUAL(package.getFiles()[0], f1);

	ClSourcePackage package2 = package << f2;

	auto files = package2.getFiles();
	BOOST_REQUIRE_EQUAL(files.size(), 2);
	BOOST_REQUIRE_EQUAL(files[0], f1);
	BOOST_REQUIRE_EQUAL(files[1], f2);
}

BOOST_AUTO_TEST_CASE( PackageAddPackage )
{
	ClSourcePackage package1 = ClSourcePackage() << "foo" << "bar";

	ClSourcePackage package2 = ClSourcePackage() << "bli" << package1 << "bar";

	BOOST_REQUIRE_EQUAL(package2.getFiles().size(), 4);
}

BOOST_AUTO_TEST_CASE( KernelCreation )
{
	//
	// BEFORE THE ACTUAL TEST CASE DO ALL THE OPENCL SETUP
	//

	// dump the compile results
	// the cast is safe here, we don't want to modify the value later
	setenv("GPU_DUMP_DEVICE_KERNEL", "2", 0); // can be overriden from outside

	cl_int err;

	// initialize the device
	const cl_uint MAX_PLATFORMS = 4;
	cl_platform_id platforms[ MAX_PLATFORMS ];
	cl_uint nPlatforms = 0;

	cl_device_id devices[ 32 ];
	cl_uint nDevices;

	// try to get the platform (currently it should always be one)
	err = clGetPlatformIDs( MAX_PLATFORMS, platforms, &nPlatforms );
	if( err ) {
		BOOST_TEST_MESSAGE( "Failed to get platforms: " << err );
		throw std::exception();
	}

	if( ! nPlatforms ) {
		BOOST_TEST_MESSAGE( "Failed to find an OpenCL platform." );
		throw std::exception();
	}

	if( nPlatforms > 1 ) {
		BOOST_TEST_MESSAGE( "Found more than one OpenCL platform, will use a random one." );
	}

	// currently there is always only one
	// As the AMD implementation already supports ICD we have to specify the platform
	// we want to use
	cl_context_properties context_props[3] = {
		CL_CONTEXT_PLATFORM,
		(cl_context_properties) platforms[0],
		0
	};

	cl_context context = clCreateContextFromType( context_props, CL_DEVICE_TYPE_ALL, NULL, NULL, &err );
	if( err ) {
		BOOST_TEST_MESSAGE( "Failed to create context: " << err );
		throw std::exception();
	}

	// get devices from context
	{
		size_t returned_size;
		err = clGetContextInfo( context, CL_CONTEXT_DEVICES, sizeof( devices ), devices, &returned_size);
		if( err ) {
			BOOST_TEST_MESSAGE( "Failed get Devices from Context: " << err );
			throw std::exception();
		}
		nDevices = returned_size / sizeof( cl_device_id );
	}

	if( ! nDevices ) {
		BOOST_TEST_MESSAGE( "No OpenCL devices found." );
		throw std::exception();
	}

	// just grab the first device
	cl_device_id device = devices[0];

	//
	// TESTCASE
	//

	ClSourcePackage dummyFunctions = ClSourcePackage() << "../hardware/dummyFunctions.cl";

	cl_kernel testKernel = TmpClKernel("dummyKernel", "", context, device) << dummyFunctions << "../hardware/dummyKernel.cl";

	// to verify try to get name from kernel object
	char kernel_name[64];
	err = clGetKernelInfo(testKernel, CL_KERNEL_FUNCTION_NAME, sizeof( kernel_name ), kernel_name, NULL);
	if( err ) {
		BOOST_TEST_MESSAGE( "Failed to get kernel property - CL_KERNEL_FUNCITON_NAME: " << err );
		throw std::exception();
	}
	BOOST_CHECK_EQUAL(CL_SUCCESS, err);
	BOOST_CHECK_EQUAL(std::string("dummyKernel"), std::string(kernel_name));

	//
	// BE NICE AND DO SOME CLEANUP
	//
}
