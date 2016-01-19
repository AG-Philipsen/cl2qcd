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

const ReferenceValues calculateReferenceValues_GaugefieldUpdate(const int latticeVolume, GaugefieldFillType fillTypeIn)
{
	switch( fillTypeIn )
	{
		case GaugefieldFillType::cold :
		{
			return ReferenceValues{6. * latticeVolume};
		}
		case GaugefieldFillType::nonTrivial:
		{
			return ReferenceValues{6. * latticeVolume * 1.0000016959666707};
		}
		default:
		{
			return defaultReferenceValues();
		}
	}
}

const ReferenceValues calculateReferenceValues_FGauge(const int latticeVolume)
{
	return ReferenceValues{8. * NDIM * latticeVolume};
}

const ReferenceValues calculateReferenceValues_FGaugeTlsym(const int latticeVolume, GaugefieldFillType fillTypeIn)
{
	switch( fillTypeIn )
	{
		case GaugefieldFillType::cold :
		{
			return ReferenceValues{8. * NDIM * latticeVolume};
		}
		case GaugefieldFillType::nonTrivial:
		{
			return ReferenceValues{8.000000000020442 * NDIM * latticeVolume};
		}
		default:
		{
			return defaultReferenceValues();
		}
	}
}

const ReferenceValues calculateReferenceValues_FFermion()
{
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_FFermionEvenOdd()
{
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_FFermionEvenOddComparator()
{
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_FFermionStaggeredEvenOdd()
{
	return defaultReferenceValues();
}

struct MolecularDynamicsTestParameters : public GaugemomentumTestParameters
{
	MolecularDynamicsTestParameters(const LatticeExtents latticeExtendsIn, GaugefieldFillType fillTypeIn) :
		GaugemomentumTestParameters(latticeExtendsIn), gaugeFillType(fillTypeIn) {}
	const GaugefieldFillType gaugeFillType;
};

struct GaugeFieldUpdateTestParameters : public MolecularDynamicsTestParameters
{
	GaugeFieldUpdateTestParameters(const LatticeExtents latticeExtendsIn, GaugefieldFillType fillTypeIn, const hmc_float epsIn) :
		MolecularDynamicsTestParameters(latticeExtendsIn, fillTypeIn), eps(epsIn) {}
	const hmc_float eps;
};

struct MolecularDynamicsTester : public GaugemomentumTester
{
	MolecularDynamicsTester(std::string kernelName, const ParameterCollection pC, const ReferenceValues rV, const MolecularDynamicsTestParameters tP) :
		GaugemomentumTester(kernelName, pC, rV, tP)
		{
			GaugefieldCreator gf;
			GaugemomentumCreator gm(calculateAlgebraSize(tP.latticeExtents));
			gaugefieldBuffer = new hardware::buffers::SU3( gf.calculateGaugefieldSize(tP.latticeExtents), this->device);
			gaugefieldBuffer->load(gf.createGaugefield(gf.calculateGaugefieldSize(tP.latticeExtents), tP.gaugeFillType));
			gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(calculateGaugemomentumSize(tP.latticeExtents), this->device);
			code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(gm.Filltype::one) ));
			gaugefieldCode = device->getGaugefieldCode();
			molecularDynamicsCode = device->getMolecularDynamicsCode();
		}

protected:
	const hardware::code::Molecular_Dynamics * molecularDynamicsCode{nullptr};
	const hardware::code::Gaugefield * gaugefieldCode{nullptr};
	const hardware::buffers::SU3 * gaugefieldBuffer;
	const hardware::buffers::Gaugemomentum * gaugemomentumBuffer;
};

struct GaugefieldUpdateTester : public MolecularDynamicsTester
{
	GaugefieldUpdateTester(const ParameterCollection pC, const GaugeFieldUpdateTestParameters tP) :
		MolecularDynamicsTester("md_update_gaugefield", pC, calculateReferenceValues_GaugefieldUpdate(calculateLatticeVolume(tP.latticeExtents), tP.gaugeFillType), tP)
		{
			molecularDynamicsCode->md_update_gaugefield_device(gaugemomentumBuffer, gaugefieldBuffer, tP.eps);

			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);

			gaugefieldCode->plaquette_device(gaugefieldBuffer, &plaq, &tplaq, &splaq);
			plaq.dump(&kernelResult[0]);
		}
};

struct FGaugeTester : public MolecularDynamicsTester
{
	FGaugeTester(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_gauge",pC, calculateReferenceValues_FGauge(calculateLatticeVolume(tP.latticeExtents)), tP)
		{
			molecularDynamicsCode->gauge_force_device( gaugefieldBuffer, gaugemomentumBuffer);
			calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FGaugeTlsymTester : public MolecularDynamicsTester
{
	FGaugeTlsymTester(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_gauge_tlsym", pC, calculateReferenceValues_FGaugeTlsym(calculateLatticeVolume(tP.latticeExtents), tP.gaugeFillType), tP)
		{
			molecularDynamicsCode->gauge_force_tlsym_device( gaugefieldBuffer, gaugemomentumBuffer);
			calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionTester : public MolecularDynamicsTester
{
	FFermionTester(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_fermion", pC, calculateReferenceValues_FFermion(), tP)
		{
			NonEvenOddSpinorfieldCreator sf(tP.latticeExtents);
			const hardware::buffers::Plain<spinor> in1(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
			const hardware::buffers::Plain<spinor> in2(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);

			in1.load(sf.createSpinorfield(SpinorFillType::one) );
			in2.load(sf.createSpinorfield(SpinorFillType::one) );

			double kappa = 1; //@todo: make adjustable
			molecularDynamicsCode->fermion_force_device( &in1, &in2, gaugefieldBuffer, gaugemomentumBuffer, kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionEvenOddTester : public MolecularDynamicsTester
{
	FFermionEvenOddTester(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_fermion_eo", pC, calculateReferenceValues_FFermionEvenOdd(), tP)
		{
			EvenOddSpinorfieldCreator sf(tP.latticeExtents);
			const hardware::buffers::Spinor in1(calculateEvenOddSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
			const hardware::buffers::Spinor in2(calculateEvenOddSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
			sf.fillTwoSpinorBuffers(&in1, &in2);//to be changed to use SpinorFillType

			double kappa = 1; //@todo: make adjustable
			bool evenOrOdd = EVEN; //@todo: make adjustable
			molecularDynamicsCode->fermion_force_eo_device( &in1, &in2, gaugefieldBuffer, gaugemomentumBuffer, evenOrOdd,  kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionEvenOddComparator : public MolecularDynamicsTester
{
	FFermionEvenOddComparator(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_fermion compare even-odd and non-even-odd", pC, calculateReferenceValues_FFermionEvenOddComparator(), tP)
		{
			createBuffers(tP);
			fillBuffers(tP);

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
		outNonEo = new const hardware::buffers::Gaugemomentum(calculateGaugemomentumSize(tP.latticeExtents), MolecularDynamicsTester::device);
		outEo = new const hardware::buffers::Gaugemomentum(calculateGaugemomentumSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inNonEo1 = new const hardware::buffers::Plain<spinor>(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inNonEo2 = new const hardware::buffers::Plain<spinor>(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo1 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo2 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo3 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo4 = new const hardware::buffers::Spinor(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
	}
	void fillBuffers(const TestParameters tP)
	{
		GaugemomentumCreator gm(calculateAlgebraSize(tP.latticeExtents));
		MolecularDynamicsTester::code->importGaugemomentumBuffer(outNonEo, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(gm.Filltype::one) ));
		MolecularDynamicsTester::code->importGaugemomentumBuffer(outEo, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(gm.Filltype::one) ));

		NonEvenOddSpinorfieldCreator sf(tP.latticeExtents);
		EvenOddSpinorfieldCreator sfEo(tP.latticeExtents);
		inNonEo1->load(sf.createSpinorfield(SpinorFillType::one));
		inNonEo2->load(sf.createSpinorfield(SpinorFillType::one));
		sfEo.fillTwoSpinorBuffers(inEo1, inEo2); //to be changed to use SpinorFillType
		sfEo.fillTwoSpinorBuffers(inEo3, inEo4);

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

struct FFermionStaggeredEvenOddTester : public MolecularDynamicsTester, public SpinorStaggeredTester
{
	FFermionStaggeredEvenOddTester(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
	   MolecularDynamicsTester("f_staggered_fermion_eo", pC, calculateReferenceValues_FFermionStaggeredEvenOdd(), tP),
	   SpinorStaggeredTester("f_staggered_fermion_eo", pC, SpinorStaggeredTestParameters(tP.latticeExtents), 1, calculateReferenceValues_FFermionStaggeredEvenOdd() )
	{
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

template<typename TesterClass, typename TestParameters>
void callTest( const LatticeExtents latticeExtentsIn, GaugefieldFillType fillType, const hmc_float epsIn)
{
	TestParameters parametersForThisTest(latticeExtentsIn, fillType, epsIn);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt); //todo: could also use latticeExtendsIn here!
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents lE, GaugefieldFillType fillType)
{
	MolecularDynamicsTestParameters parametersForThisTest(lE, fillType);
	//todo: Work over these!
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForMolecularDynamics kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testGaugefieldUpdate(const LatticeExtents lE, GaugefieldFillType fillType, const hmc_float eps)
{
	callTest<GaugefieldUpdateTester, GaugeFieldUpdateTestParameters>(lE, fillType, eps);
}

void testFGauge(const LatticeExtents lE, GaugefieldFillType fillType)
{
	callTest<FGaugeTester>(lE, fillType);
}
void testFGaugeTlsym(const LatticeExtents lE, GaugefieldFillType fillType)
{
	callTest<FGaugeTlsymTester>(lE, fillType);
}
void testNonEvenOddFermionForce(const LatticeExtents lE, GaugefieldFillType fillType)
{
	callTest<FFermionTester>(lE, fillType);
}
void testEvenOddFermionForce(const LatticeExtents lE, GaugefieldFillType fillType)
{
	callTest<FFermionEvenOddTester>(lE, fillType);
}
void compareEvenOddAndNonEvenOddFermionForce(const LatticeExtents lE, GaugefieldFillType fillType)
{
	callTest<FFermionEvenOddComparator>(lE, fillType);
}
void testEvenOddStaggeredFermionForce(const LatticeExtents lE, GaugefieldFillType fillType)
{
	callTest<FFermionStaggeredEvenOddTester>(lE, fillType);
}

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

	//Note: test5: "THIS TEST HAS TO BE INVESTIGATED!! IT COULD HINT TO AN ERROR IN THE FUNCTION!"
	BOOST_AUTO_TEST_CASE(GF_UPDATE_1 )
	{
		testGaugefieldUpdate(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, .0);
	}
	BOOST_AUTO_TEST_CASE(GF_UPDATE_2 )
	{
		testGaugefieldUpdate(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, .12);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE )

	//@todo: add tests like in "f_gauge_input_{1-2}"
	BOOST_AUTO_TEST_CASE( F_GAUGE_1 )
	{
		testFGauge(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}
	BOOST_AUTO_TEST_CASE( F_GAUGE_2 )
	{
		testFGauge(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE_TLSYM )

	//@todo: add tests like in "f_gauge_tlsym_input_{1-3}"
	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_1 )
	{
		testFGaugeTlsym(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}
	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_2 )
	{
		testFGaugeTlsym(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
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
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}
//	BOOST_AUTO_TEST_CASE( F_FERMION_2 )
//	{
//		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
//	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_EO )

	//@todo: add tests like in "f_fermion_eo_input_{1-20}"
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_2 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_COMPARE_NONEO_EO )

	//@todo: add tests like in "f_fermion_compare_noneo_eo_input_{1-6}"
	BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_1 )
	{
		compareEvenOddAndNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_2 )
	{
		compareEvenOddAndNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_STAGG_FERMION_EO )

	//@todo: add tests like in "f_staggered_fermion_partial_eo_input_{1-32}"
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_2 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
	}

BOOST_AUTO_TEST_SUITE_END()


