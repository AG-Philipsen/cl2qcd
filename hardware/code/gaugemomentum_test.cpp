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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_GAUGEMOMENTUM

#include "GaugemomentumTester.hpp"

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  GaugemomentumTester tester("build", "opencl_module_gaugemomentum_build_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SQUARENORM )

class SquarenormTester : public GaugemomentumTester
{
public:
  SquarenormTester(std::string inputfile) :
    GaugemomentumTester("gaugemomenta squarenorm", inputfile, 1)
  {
    hardware::buffers::Gaugemomentum in(numberOfGaugemomentumElements, device);
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentum() ));
    calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
  }
};

BOOST_AUTO_TEST_CASE(SQUARENORM_1  )
{
	SquarenormTester tester("squarenorm_input_1");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_2  )
{
	SquarenormTester tester("squarenorm_input_2");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_1  )
{
	SquarenormTester tester("squarenorm_reduction_input_1");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_2  )
{
	SquarenormTester tester("squarenorm_reduction_input_2");
}

BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_3  )
{
	SquarenormTester tester("squarenorm_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SET_ZERO )

class SetZeroTester : public GaugemomentumTester
{
public:
  SetZeroTester(std::string inputfile) :
    GaugemomentumTester("set zero", inputfile, 1)
  {
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentum() ));
    code->set_zero_gaugemomentum(gaugemomentumBuffer);
    calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
  }
};

BOOST_AUTO_TEST_CASE( SET_ZERO_1 )
{
  SetZeroTester tester("set_zero_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SAXPY )

class SaxpyTester : public GaugemomentumTester
{
public:
  SaxpyTester(std::string inputfile) :
    GaugemomentumTester("saxpy", inputfile, 1)
  {
    hardware::buffers::Gaugemomentum out(numberOfGaugemomentumElements, device);
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentum(123456) ));
    code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>( createGaugemomentum(789101) ));
    double alpha = parameters->get_tau();
    doubleBuffer->load(&alpha);

    code->saxpy_device(gaugemomentumBuffer, &out, doubleBuffer, &out);
    calcSquarenormAndStoreAsKernelResult(&out);
  }
};

BOOST_AUTO_TEST_CASE( SAXPY_1 )
{
	SaxpyTester tester("saxpy_input_1");
}

BOOST_AUTO_TEST_CASE( SAXPY_2 )
{
	SaxpyTester tester("saxpy_input_2");
}

BOOST_AUTO_TEST_CASE( SAXPY_3 )
{
	SaxpyTester tester("saxpy_input_3");
}

BOOST_AUTO_TEST_CASE( SAXPY_4 )
{
	SaxpyTester tester("saxpy_input_4");
}

BOOST_AUTO_TEST_CASE( SAXPY_5 )
{
	SaxpyTester tester("saxpy_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GENERATE_GAUSSIAN_GAUGEMOMENTA  )

class GaussianTester : public GaugemomentumTester
{
public:
  GaussianTester(std::string inputfile) :
    GaugemomentumTester("gaussian gaugemomentum", inputfile, 1, 2)
  {
	physics::PRNG prng(*system);
	auto prng_buf = prng.get_buffers().at(0);
	double result = 0.;
	double sum = 0.;
	int iterations = parameters->get_integrationsteps(0);
	ae * gm_out;

	gm_out = new ae[numberOfGaugemomentumElements * iterations];
	BOOST_REQUIRE(gm_out);

	for (int i = 0; i< iterations; i++){
	  code->generate_gaussian_gaugemomenta_device(gaugemomentumBuffer, prng_buf);
	  gaugemomentumBuffer->dump(&gm_out[i*numberOfGaugemomentumElements]);
	  sum += count_gm(&gm_out[i*numberOfGaugemomentumElements], numberOfGaugemomentumElements);
	}
	sum = sum/iterations/numberOfGaugemomentumElements/8;	
	result= sum;

	if(parameters->get_read_multiple_configs()  == false){
	  double var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_gm(&gm_out[i*numberOfGaugemomentumElements], numberOfGaugemomentumElements, sum);
	  }
	  var=var/iterations/numberOfGaugemomentumElements/8;
	  
	  result = sqrt(var);
	}

	kernelResult[0] = result;
  }
};

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_1 )
{
  GaussianTester tester("gaussian_input_1");
}

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_2 )
{
  GaussianTester tester("gaussian_input_2");
}

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_3 )
{
  GaussianTester tester("gaussian_input_3");
}

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_4 )
{
  GaussianTester tester("gaussian_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( CONVERT_TO_SOA )

class ConvertToSoaTester : public GaugemomentumTester
{
public:
  ConvertToSoaTester(std::string inputfile) :
    GaugemomentumTester("convert to soa", inputfile, 1)
  {
    BOOST_TEST_MESSAGE("NOT YET IMPLEMENTED!!");
  }
};

BOOST_AUTO_TEST_CASE( CONVERT_TO_SOA_1 )
{
  ConvertToSoaTester tester("");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( CONVERT_FROM_SOA )

class ConvertFromSoaTester : public GaugemomentumTester
{
public:
  ConvertFromSoaTester(std::string inputfile) :
    GaugemomentumTester("convert from soa", inputfile, 1)
  {
    BOOST_TEST_MESSAGE("NOT YET IMPLEMENTED!!");
  }
};

BOOST_AUTO_TEST_CASE( CONVERT_FROM_SOA_1 )
{
  ConvertFromSoaTester tester("");
}

BOOST_AUTO_TEST_SUITE_END()

