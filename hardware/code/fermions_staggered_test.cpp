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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_STAGGERED

#include "SpinorStaggeredTester.hpp"
#include "fermions_staggered.hpp"

//@todo: add parameters as tests are refactored
struct StaggeredFermionsTestParameters : public SpinorStaggeredTestParameters
{
	StaggeredFermionsTestParameters(const LatticeExtents lE):
		TestParameters(lE), SpinorStaggeredTestParameters(lE) {};
};

struct StaggeredFermionsTester : public SpinorStaggeredTester2
{
	StaggeredFermionsTester(const std::string kernelName, const ParameterCollection pC, const SpinorStaggeredTestParameters tP, const ReferenceValues rV):
	SpinorStaggeredTester2(kernelName, pC, tP, rV)
    {
		code = device->getFermionStaggeredCode();
		gaugefieldBuffer = new hardware::buffers::SU3( LatticeExtents(tP.latticeExtents).getLatticeVolume() , device);
	}
	
	virtual ~StaggeredFermionsTester(){
		code = nullptr;
	}
	
   protected:
	const hardware::code::Fermions_staggered * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

struct NonEvenOddStaggeredFermionmatrixTester : public StaggeredFermionsTester
{
	NonEvenOddStaggeredFermionmatrixTester(const std::string kernelName, const ParameterCollection pC, const SpinorStaggeredTestParameters tP, const ReferenceValues rV) :
	    StaggeredFermionsTester(kernelName, pC, tP, rV){
		in = new const hardware::buffers::Plain<su3vec>(calculateSpinorfieldSize(tP.latticeExtents), device);
		out = new const hardware::buffers::Plain<su3vec>(calculateSpinorfieldSize(tP.latticeExtents), device);
		NonEvenOddSpinorStaggeredfieldCreator ssf(tP.latticeExtents);
		in->load(ssf.createSpinorfield(SpinorFillType::one)); //@todo: make adjustable
		out->load(ssf.createSpinorfield(SpinorFillType::one));//@todo: make adjustable
	}
	virtual ~NonEvenOddStaggeredFermionmatrixTester(){
		calcSquarenormAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
protected:
	const hardware::buffers::Plain<su3vec> * in;
	const hardware::buffers::Plain<su3vec> * out;
};

struct MTester : public NonEvenOddStaggeredFermionmatrixTester
{
	MTester(const ParameterCollection pC, const SpinorStaggeredTestParameters tP) :
		NonEvenOddStaggeredFermionmatrixTester("M_staggered", pC, tP, defaultReferenceValues())
	{
			code->M_staggered_device(in, out,  gaugefieldBuffer, nonTrivialParameter); //@todo: make mass adjustable
	}
};

struct EvenOddStaggeredFermionmatrixTester : public StaggeredFermionsTester
{
	EvenOddStaggeredFermionmatrixTester(std::string kernelName, const ParameterCollection pC, const SpinorStaggeredTestParameters tP, const ReferenceValues rV) :
	   StaggeredFermionsTester(kernelName, pC, tP, rV){
		in = new const hardware::buffers::SU3vec(tP.latticeExtents, device);
		out = new const hardware::buffers::SU3vec(tP.latticeExtents, device);
		EvenOddSpinorStaggeredfieldCreator ssf(tP.latticeExtents);
		in->load(ssf.createSpinorfield(SpinorFillType::one)); //@todo: make adjustable
		out->load(ssf.createSpinorfield(SpinorFillType::one)); //@todo: make adjustable
	}
	~EvenOddStaggeredFermionmatrixTester(){
		calcSquarenormEvenOddAndStoreAsKernelResult(out);
		delete in;
		delete out;
	}
   protected:
	const hardware::buffers::SU3vec * in;
	const hardware::buffers::SU3vec * out;
};

struct DksEvenOddTester : public EvenOddStaggeredFermionmatrixTester
{
	DksEvenOddTester(const ParameterCollection pC, const SpinorStaggeredTestParameters tP, const bool evenOrOdd = ODD) : //@todo: remove the default value for evenOrOdd
		EvenOddStaggeredFermionmatrixTester("DKS_eo", pC, tP, defaultReferenceValues())
	{
			evenOrOdd ? code->D_KS_eo_device(in, out,  gaugefieldBuffer, EVEN)
			  : code->D_KS_eo_device(in, out,  gaugefieldBuffer, ODD);
	}
};

template<class TesterClass>
void callTest(const LatticeExtents lE)
{
	StaggeredFermionsTestParameters parametersForThisTest(lE);
	//todo: Work over these!
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForSpinorStaggered kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testMMatrix(const LatticeExtents lE)
{
	callTest<MTester>(lE);
}

void testDksEo(const LatticeExtents lE)
{
	callTest<DksEvenOddTester>(lE);
}

BOOST_AUTO_TEST_SUITE( M_MATRIX )

	//@todo: add tests like in "m_input_{1-12}
	BOOST_AUTO_TEST_CASE( M_MATRIX_1)
	{
		testMMatrix(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( DKS_EO )

	//@todo: add tests like in "dks_input_{1-32}
	BOOST_AUTO_TEST_CASE( DKS_EO_1)
	{
		testDksEo(LatticeExtents{ns4, nt4});
	}
	
BOOST_AUTO_TEST_SUITE_END()

