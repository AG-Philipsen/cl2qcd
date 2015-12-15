/*
 * Copyright 2012, 2013, 2014 Lars Zeidlewicz, Christopher Pinke,
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
#define BOOST_TEST_MODULE HARDWARE_CODE_MOLECULAR_DYNAMICS

#include "GaugemomentumTester.hpp"
#include "GaugefieldTester.hpp"
#include "SpinorTester.hpp"
#include "SpinorStaggeredTester.hpp"
#include "molecular_dynamics.hpp"

/**
 * NOTE: These tests were quickly refactored to use the "new" structure of tests,
 *   i.e. not to use inputfiles. In the process, no real tests were implemented.
 *   This has still to be done. That is why currently, the test results are not
 *   meaningful!
 */

struct MolecularDynamicsTester : public GaugemomentumTester
{
	MolecularDynamicsTester(std::string kernelName, const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
		GaugemomentumTester(kernelName, pC, rV, tP)
		{
			gaugefieldBuffer = new hardware::buffers::SU3( calculateGaugefieldSize(tP.latticeExtents), this->device);
			gaugefieldBuffer->load(createGaugefield(calculateGaugefieldSize(tP.latticeExtents), GaugefieldFillType::cold)); //@todo: make adjustable
			gaugefieldCode = device->getGaugefieldCode();
			molecularDynamicsCode = device->getMolecularDynamicsCode();
		}

protected:
	const hardware::code::Molecular_Dynamics * molecularDynamicsCode{nullptr};
	const hardware::code::Gaugefield * gaugefieldCode{nullptr};
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

struct GaugefieldUpdateTester : public MolecularDynamicsTester
{
	GaugefieldUpdateTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
		MolecularDynamicsTester("md_update_gaugefield", pC, rV, tP)
		{
			code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentum() ));
			double epsilon = 1.; //@todo: make adjustable
			molecularDynamicsCode->md_update_gaugefield_device(gaugemomentumBuffer, gaugefieldBuffer, epsilon);

			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);

			gaugefieldCode->plaquette_device(gaugefieldBuffer, &plaq, &tplaq, &splaq);
			plaq.dump(&kernelResult[0]);
		}
};

struct FGaugeTester : public MolecularDynamicsTester
{
	FGaugeTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
		MolecularDynamicsTester("f_gauge",pC, rV, tP)
		{
			code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
			molecularDynamicsCode->gauge_force_device( gaugefieldBuffer, gaugemomentumBuffer);
			calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FGaugeTlsymTester : public MolecularDynamicsTester
{
	FGaugeTlsymTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
		MolecularDynamicsTester("f_gauge_tlsym", pC, rV, tP)
		{
			code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
			molecularDynamicsCode->gauge_force_tlsym_device( gaugefieldBuffer, gaugemomentumBuffer);
			calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionTester : public MolecularDynamicsTester, public SpinorTester
{
	FFermionTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
		MolecularDynamicsTester("f_fermion", pC, rV, tP),
		SpinorTester("f_fermion", pC, SpinorTestParameters(tP.latticeExtents), 1, rV)
		{
			MolecularDynamicsTester::code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
			const hardware::buffers::Plain<spinor> in1(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
			const hardware::buffers::Plain<spinor> in2(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);

			in1.load(SpinorTester::createSpinorfield(SpinorFillType::one) );
			in2.load(SpinorTester::createSpinorfield(SpinorFillType::one) );

			double kappa = 1; //@todo: make adjustable
			molecularDynamicsCode->fermion_force_device( &in1, &in2, gaugefieldBuffer, gaugemomentumBuffer, kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionEvenOddTester : public MolecularDynamicsTester, public SpinorTester
{
	FFermionEvenOddTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
		MolecularDynamicsTester("f_fermion_eo", pC, rV, tP),
		SpinorTester("f_fermion_eo", pC, SpinorTestParameters(tP.latticeExtents), 1, rV)
		{
			MolecularDynamicsTester::code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
			const hardware::buffers::Spinor in1(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
			const hardware::buffers::Spinor in2(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
			fillTwoSpinorBuffers(&in1, &in2);

			double kappa = 1; //@todo: make adjustable
			bool evenOrOdd = EVEN; //@todo: make adjustable
			molecularDynamicsCode->fermion_force_eo_device( &in1, &in2, gaugefieldBuffer, gaugemomentumBuffer, evenOrOdd,  kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionEvenOddComparator : public MolecularDynamicsTester, public SpinorTester
{
	FFermionEvenOddComparator(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
		MolecularDynamicsTester("f_fermion compare even-odd and non-even-odd", pC, rV, tP),
		SpinorTester("f_fermion compare even-odd and non-even-odd", pC, SpinorTestParameters(tP.latticeExtents), 1, rV)
		{
			createBuffers(tP);
			fillBuffers();

			double kappa = 1; //@todo: make adjustable

			hmc_float cpu_res_eo;
			molecularDynamicsCode->fermion_force_eo_device(inEo1, inEo4, gaugefieldBuffer, outEo, ODD, kappa );
			MolecularDynamicsTester::code->set_float_to_gaugemomentum_squarenorm_device(outEo, MolecularDynamicsTester::doubleBuffer);
			MolecularDynamicsTester::doubleBuffer->dump(&cpu_res_eo);
			logger.info() << "|force_eo (even) + force_eo (odd)|^2:";
			logger.info() << cpu_res_eo;

			molecularDynamicsCode->fermion_force_eo_device(inEo2, inEo3, gaugefieldBuffer, outEo, EVEN, kappa );
			MolecularDynamicsTester::code->set_float_to_gaugemomentum_squarenorm_device(outEo, MolecularDynamicsTester::doubleBuffer);
			MolecularDynamicsTester::doubleBuffer->dump(&cpu_res_eo);
			logger.info() << "|force_eo (even) + force_eo (odd)|^2:";
			logger.info() << cpu_res_eo;

			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(outEo, 0);

			molecularDynamicsCode->fermion_force_device( inNonEo1, inNonEo2, gaugefieldBuffer, outNonEo, kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(outNonEo, 1);
		}
		~FFermionEvenOddComparator()
		{
			delete	outNonEo;
			delete	outEo;
			delete	inNonEo1;
			delete	inNonEo2;
			delete	inEo1;
			delete	inEo2;
			delete	inEo3;
			delete	inEo4;

			outNonEo = NULL;
			outEo = NULL;
			inNonEo1 = NULL;
			inNonEo2 = NULL;
			inEo1 = NULL;
			inEo2 = NULL;
			inEo3 = NULL;
			inEo4 = NULL;
		}

private:
	void createBuffers(const TestParameters tP)
	{
		outNonEo = new const hardware::buffers::Gaugemomentum(MolecularDynamicsTester::numberOfGaugemomentumElements, MolecularDynamicsTester::device);
		outEo = new const hardware::buffers::Gaugemomentum(MolecularDynamicsTester::numberOfGaugemomentumElements, MolecularDynamicsTester::device);
		inNonEo1 = new const hardware::buffers::Plain<spinor>(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inNonEo2 = new const hardware::buffers::Plain<spinor>(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo1 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo2 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo3 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo4 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
	}
	void fillBuffers()
	{
		MolecularDynamicsTester::code->importGaugemomentumBuffer(outNonEo, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
		MolecularDynamicsTester::code->importGaugemomentumBuffer(outEo, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));

		inNonEo1->load(SpinorTester::createSpinorfield(SpinorFillType::one));
		inNonEo2->load(SpinorTester::createSpinorfield(SpinorFillType::one));
		fillTwoSpinorBuffers(inEo1, inEo2, 123456);
		fillTwoSpinorBuffers(inEo3, inEo4, 789101);

	}

	const hardware::buffers::Gaugemomentum * outNonEo;
	const hardware::buffers::Gaugemomentum * outEo;
	const hardware::buffers::Plain<spinor> * inNonEo1;
	const hardware::buffers::Plain<spinor> * inNonEo2;
	const hardware::buffers::Spinor * inEo1;
	const hardware::buffers::Spinor * inEo2;
	const hardware::buffers::Spinor * inEo3;
	const hardware::buffers::Spinor * inEo4;
};

struct FFermionStaggeredEvenOddTester : public MolecularDynamicsTester, public SpinorStaggeredTester{
	FFermionStaggeredEvenOddTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
	   MolecularDynamicsTester("f_staggered_fermion_eo", pC, rV, tP),
	   SpinorStaggeredTester("f_staggered_fermion_eo", pC, SpinorStaggeredTestParameters(tP.latticeExtents), 1, rV )
	{

		MolecularDynamicsTester::code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( createGaugemomentumBasedOnFilltype(zero) ));
		const hardware::buffers::SU3vec in1(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		in1.load(createSpinorfield(SpinorFillType::one));
		const hardware::buffers::SU3vec in2(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		in2.load(createSpinorfield(SpinorFillType::one));

		bool evenOrOdd = true; //@todo: make a test parameter
		if(evenOrOdd){
		  molecularDynamicsCode->fermion_staggered_partial_force_device(gaugefieldBuffer, &in1, &in2, gaugemomentumBuffer, EVEN);
		}else{
		  molecularDynamicsCode->fermion_staggered_partial_force_device(gaugefieldBuffer, &in1, &in2, gaugemomentumBuffer, ODD);
		}

		MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
/*
print_staggeredfield_eo_to_textfile("ref_vec_f_stagg1_eo", createSpinorfield(spinorfieldEvenOddElements, 123));
logger.info() << "Produced the ref_vec_f_stagg1_eo text file with the staggered field for the ref. code.";
print_staggeredfield_eo_to_textfile("ref_vec_f_stagg2_eo", createSpinorfield(spinorfieldEvenOddElements, 456));
logger.info() << "Produced the ref_vec_f_stagg2_eo text file with the staggered field for the ref. code.";
*/
		}
};

template<class TesterClass>
void callTest(const LatticeExtents lE)
{
	GaugemomentumTestParameters parametersForThisTest(lE);
	//todo: Work over these!
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, defaultReferenceValues(), parametersForThisTest);
}

void testGaugefieldUpdate(const LatticeExtents lE)
{
	callTest<GaugefieldUpdateTester>(lE);
}

void testFGauge(const LatticeExtents lE)
{
	callTest<FGaugeTester>(lE);
}
void testFGaugeTlsym(const LatticeExtents lE)
{
	callTest<FGaugeTlsymTester>(lE);
}
void testNonEvenOddFermionForce(const LatticeExtents lE)
{
	callTest<FFermionTester>(lE);
}
void testEvenOddFermionForce(const LatticeExtents lE)
{
	callTest<FFermionEvenOddTester>(lE);
}
void compareEvenOddAndNonEvenOddFermionForce(const LatticeExtents lE)
{
	callTest<FFermionEvenOddComparator>(lE);
}
void testEvenOddStaggeredFermionForce(const LatticeExtents lE)
{
	callTest<FFermionStaggeredEvenOddTester>(lE);
}

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

	//@todo: add tests like in "gf_update_input_{1-5}"
	//Note: test5: "THIS TEST HAS TO BE INVESTIGATED!! IT COULD HINT TO AN ERROR IN THE FUNCTION!"
	BOOST_AUTO_TEST_CASE(GF_UPDATE_1 )
	{
		testGaugefieldUpdate(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE )

	//@todo: add tests like in "f_gauge_input_{1-2}"
	BOOST_AUTO_TEST_CASE( F_GAUGE_1 )
	{
		testFGauge(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE_TLSYM )

	//@todo: add tests like in "f_gauge_tlsym_input_{1-3}"
	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_1 )
	{
		testFGaugeTlsym(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( STOUT_SMEAR_FERMION_FORCE )

	BOOST_AUTO_TEST_CASE(STOUT_SMEAR_FERMION_FORCE_1 )
	{
		BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION )

	//@todo: add tests like in "f_fermion_input_{1-6}"
	BOOST_AUTO_TEST_CASE( F_FERMION_1 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_EO )

	//@todo: add tests like in "f_fermion_eo_input_{1-20}"
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_COMPARE_NONEO_EO )

	//@todo: add tests like in "f_fermion_compare_noneo_eo_input_{1-6}"
	BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_1 )
	{
		compareEvenOddAndNonEvenOddFermionForce(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_STAGG_FERMION_EO )

	//@todo: add tests like in "f_staggered_fermion_partial_eo_input_{1-32}"
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()


