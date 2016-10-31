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
#include "FermionTester.hpp"
#include "SpinorStaggeredTester.hpp"
#include "molecular_dynamics.hpp"

const ReferenceValues calculateReferenceValues_GaugefieldUpdate(const int latticeVolume, GaugefieldFillType fillTypeIn, const GaugeMomentumFilltype gmFillTypeIn)
{
	switch( fillTypeIn )
	{
		case GaugefieldFillType::cold :
		{
			if ( gmFillTypeIn == GaugeMomentumFilltype::One )
				return ReferenceValues{6. * latticeVolume};
			else if( gmFillTypeIn == GaugeMomentumFilltype::Ascending )
				return ReferenceValues{6. * latticeVolume};
			else return defaultReferenceValues();
		}
		case GaugefieldFillType::nonTrivial:
		{
			if ( gmFillTypeIn == GaugeMomentumFilltype::One )
				return ReferenceValues{6. * latticeVolume * 1.0000016959666707};
			else if( gmFillTypeIn == GaugeMomentumFilltype::Ascending )
				return ReferenceValues{6. * latticeVolume * 1.0000016959666707};
			else return defaultReferenceValues();
		}
		default:
			return defaultReferenceValues();
	}
}

const ReferenceValues calculateReferenceValues_FGaugeTlsym(const int latticeVolume, const GaugefieldFillType fillTypeIn, const GaugeMomentumFilltype gmFillTypeIn)
{
	switch( fillTypeIn )
	{
		case GaugefieldFillType::cold :
		{
			if ( gmFillTypeIn == GaugeMomentumFilltype::One )
				return ReferenceValues{8.000000000000002 * NDIM * latticeVolume};
			else if( gmFillTypeIn == GaugeMomentumFilltype::Ascending )
				return ReferenceValues{204. * NDIM * latticeVolume};
			else return defaultReferenceValues();
		}
		case GaugefieldFillType::nonTrivial:
		{
			if ( gmFillTypeIn == GaugeMomentumFilltype::One )
				return ReferenceValues{7.999999999997481 * NDIM * latticeVolume};
			else if( gmFillTypeIn == GaugeMomentumFilltype::Ascending )
				return ReferenceValues{203.9999999999905 * NDIM * latticeVolume};
			else return defaultReferenceValues();
		}
		default:
			return defaultReferenceValues();
	}
}

const ReferenceValues calculateReferenceValues_FGauge(const int latticeVolume, const GaugefieldFillType fillTypeIn, const GaugeMomentumFilltype gmFillTypeIn)
{
	switch( fillTypeIn )
	{
		case GaugefieldFillType::cold :
		{
			if ( gmFillTypeIn == GaugeMomentumFilltype::One )
				return ReferenceValues{8.000000000000002 * NDIM * latticeVolume};
			else if( gmFillTypeIn == GaugeMomentumFilltype::Ascending )
				return ReferenceValues{204. * NDIM * latticeVolume};
			else return defaultReferenceValues();
		}
		case GaugefieldFillType::nonTrivial:
		{
			if ( gmFillTypeIn == GaugeMomentumFilltype::One )
				return ReferenceValues{8.000000000008056 * NDIM * latticeVolume};
			else if( gmFillTypeIn == GaugeMomentumFilltype::Ascending )
				return ReferenceValues{204.0000000000302 * NDIM * latticeVolume};
			else return defaultReferenceValues();
		}
		default:
			return defaultReferenceValues();
	}
}

const ReferenceValues calculateReferenceValues_FFermion(const int latticeVolume, const GaugefieldFillType gfFillType, const GaugeMomentumFilltype gmFillType, const SpinorFillType sfFillType, const WilsonMassParameters massParameters)
{
	if ( massParameters.kappa == nonTrivialParameter )
	{
		switch ( gfFillType )
		{
			case GaugefieldFillType::cold:
			{
				if ( gmFillType == GaugeMomentumFilltype::One )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 55.117979451392 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 3649859.391025437 * latticeVolume};
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 775.9085074513919 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 3646130.956279331 * latticeVolume};
				}
			}
			break;
			case GaugefieldFillType::nonTrivial:
			{
				if ( gmFillType == GaugeMomentumFilltype::One )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 49.40336833056013 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 13927166.27677378 * latticeVolume};
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 817.7260398734077 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 14033225.67056308 * latticeVolume};
				}
			}
			break;
			default:
				return defaultReferenceValues();
		}
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_FFermionEvenOdd(const int latticeVolume, const GaugefieldFillType gfFillType, const GaugeMomentumFilltype gmFillType, const SpinorFillType sfFillType, const WilsonMassParameters massParameters)
{
	if ( massParameters.kappa == nonTrivialParameter )
	{
		switch ( gfFillType )
		{
			case GaugefieldFillType::cold:
			{
				if ( gmFillType == GaugeMomentumFilltype::One )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 31.853606862848 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 912111.817738284 * latticeVolume};
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 784.248870862848 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 910639.6003652319 * latticeVolume};
				}
			}
				break;
			case GaugefieldFillType::nonTrivial:
			{
				if ( gmFillType == GaugeMomentumFilltype::One )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 34.7169542449646 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 3490716.491209607 * latticeVolume};
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 810.8782900163883 * latticeVolume};
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 3544138.188104251 * latticeVolume};
				}
			}
				break;
			default:
				return defaultReferenceValues();
		}
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_FFermionStaggeredEvenOdd(const int latticeVolume, const GaugefieldFillType gfFillType, const GaugeMomentumFilltype gmFillType, const SpinorFillType sfFillType, const bool evenOrOdd)
{
	// in this function the ratios (5/8., 3/8.) are due to the counting of how many +1 (or -1) staggered phases are there and to the will of factorising half the lattice volume (due to EO prec.) and the number of possible directions (4)
	switch (evenOrOdd)
	{
	case EVEN:
	{
		switch( gfFillType )
		{
		case GaugefieldFillType::cold:
			{
				if ( gmFillType == GaugeMomentumFilltype::One ||  sfFillType == SpinorFillType::ascendingComplex )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 8. * latticeVolume * 4 };
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one ||  sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 204. * latticeVolume * 4 };
				}
			}
			break;
		case GaugefieldFillType::nonTrivial:
			{
				if ( gmFillType == GaugeMomentumFilltype::One )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ ( 8. + (30.81279381225022 * (5 / 8.) + 1.581192662056448 * (3 / 8.) )) * latticeVolume / 2 * 4 };
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ ( 8. + (8621.652426649009 * (5 / 8.) + 7831.388151217087 * (3 / 8.) )) * latticeVolume / 2 * 4 };
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ ( 204. + (270.9717798379283 * (5 / 8.) + 153.4222066363782 * (3 / 8.) )) * latticeVolume / 2 * 4 };
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ ( 204. + (9839.968790260735 * (5 / 8.) + 7005.071787605359 * (3 / 8.) )) * latticeVolume / 2 * 4 };
				}
			}
			break;
		default:
			return defaultReferenceValues();
		}
	}
		break;
	case ODD:
	{
		switch( gfFillType )
		{
		case GaugefieldFillType::cold:
			{
				if ( gmFillType == GaugeMomentumFilltype::One )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 8. * latticeVolume * 4 };
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 8. * latticeVolume * 4 };
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ 204. * latticeVolume * 4 };
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ 204. * latticeVolume * 4 };
				}
			}
			break;
		case GaugefieldFillType::nonTrivial:
			{
				if ( gmFillType == GaugeMomentumFilltype::One )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ ( 8. + (30.81279381225022 * (3 / 8.) + 1.581192662056448 * (5 / 8.) )) * latticeVolume / 2 * 4 };
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ ( 8. + (8621.652426649009 * (3 / 8.) + 7831.388151217087 * (5 / 8.) )) * latticeVolume / 2 * 4 };
				}
				else if ( gmFillType == GaugeMomentumFilltype::Ascending )
				{
					if ( sfFillType == SpinorFillType::one )
						return ReferenceValues{ ( 204. + (270.9717798379283 * (3 / 8.) + 153.4222066363782 * (5 / 8.) )) * latticeVolume / 2 * 4 };
					else if ( sfFillType == SpinorFillType::ascendingComplex )
						return ReferenceValues{ ( 204. + (9839.968790260735 * (3 / 8.) + 7005.071787605359 * (5 / 8.) )) * latticeVolume / 2 * 4 };
				}
			}
			break;
		default:
			return defaultReferenceValues();
		}
	}
		break;
	default:
		return defaultReferenceValues();
	}
		return defaultReferenceValues();
}

struct MolecularDynamicsTestParameters : public GaugemomentumTestParameters
{
	MolecularDynamicsTestParameters(const LatticeExtents latticeExtendsIn, GaugefieldFillType fillTypeIn, GaugeMomentumFilltype gmFillTypeIn) :
		GaugemomentumTestParameters(latticeExtendsIn), gaugeFillType(fillTypeIn), gmFillType(gmFillTypeIn), spinorFillType(SpinorFillType::one), massParameters(1.) {}
	MolecularDynamicsTestParameters(const LatticeExtents latticeExtendsIn, GaugefieldFillType gfFillTypeIn, GaugeMomentumFilltype gmFillTypeIn, SpinorFillType sfFillTypeIn, WilsonMassParameters kappaIn) :
		GaugemomentumTestParameters(latticeExtendsIn, gmFillTypeIn), gaugeFillType(gfFillTypeIn), gmFillType(gmFillTypeIn), spinorFillType(sfFillTypeIn), massParameters(kappaIn) {}
	const GaugefieldFillType gaugeFillType;
	const GaugeMomentumFilltype gmFillType;
	const SpinorFillType spinorFillType;
	WilsonMassParameters massParameters;

};

struct EvenOddMolecularDynamicsTestParameters: public MolecularDynamicsTestParameters
{
	EvenOddMolecularDynamicsTestParameters(const LatticeExtents latticeExtendsIn, GaugefieldFillType gfFillTypeIn, GaugeMomentumFilltype gmFillTypeIn, SpinorFillType sfFillTypeIn, WilsonMassParameters kappaIn, const bool evenOrOddIn ):
		MolecularDynamicsTestParameters( latticeExtendsIn, gfFillTypeIn, gmFillTypeIn, sfFillTypeIn, kappaIn), evenOrOdd(evenOrOddIn) {}
	EvenOddMolecularDynamicsTestParameters(const LatticeExtents latticeExtendsIn, GaugefieldFillType gfFillTypeIn, GaugeMomentumFilltype gmFillTypeIn, SpinorFillType sfFillTypeIn, const bool evenOrOddIn ):
		MolecularDynamicsTestParameters( latticeExtendsIn, gfFillTypeIn, gmFillTypeIn, sfFillTypeIn, 1.), evenOrOdd(evenOrOddIn) {}
	const bool evenOrOdd;
};

struct GaugeFieldUpdateTestParameters : public MolecularDynamicsTestParameters
{
	GaugeFieldUpdateTestParameters(const LatticeExtents latticeExtendsIn, GaugefieldFillType fillTypeIn, GaugeMomentumFilltype gmFillType, const hmc_float epsIn) :
		MolecularDynamicsTestParameters(latticeExtendsIn, fillTypeIn, gmFillType), eps(epsIn) {}
	const hmc_float eps;
};

struct MolecularDynamicsTester : public GaugemomentumTester
{
	MolecularDynamicsTester(std::string kernelName, const ParameterCollection pC, const ReferenceValues rV, const MolecularDynamicsTestParameters tP) :
		GaugemomentumTester(kernelName, pC, rV, tP)
		{
			GaugefieldCreator gf(tP.latticeExtents);
			GaugemomentumCreator gm(tP.latticeExtents);
			gaugefieldBuffer = new hardware::buffers::SU3( calculateGaugefieldSize(tP.latticeExtents), this->device);
			const Matrixsu3 * gf_host = gf.createGaugefield(tP.gaugeFillType);
	        device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
	        delete[] gf_host;
			gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(tP.latticeExtents, this->device);
			code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(tP.gmFillType) ));
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
		MolecularDynamicsTester("md_update_gaugefield", pC, calculateReferenceValues_GaugefieldUpdate(LatticeExtents(tP.latticeExtents).getLatticeVolume(), tP.gaugeFillType, tP.gmFillType), tP)
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
		MolecularDynamicsTester("f_gauge",pC, calculateReferenceValues_FGauge(LatticeExtents(tP.latticeExtents).getLatticeVolume(), tP.gaugeFillType, tP.gmFillType), tP)
		{
			molecularDynamicsCode->gauge_force_device( gaugefieldBuffer, gaugemomentumBuffer);
			calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FGaugeTlsymTester : public MolecularDynamicsTester
{
	FGaugeTlsymTester(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_gauge_tlsym", pC, calculateReferenceValues_FGaugeTlsym(LatticeExtents(tP.latticeExtents).getLatticeVolume(), tP.gaugeFillType, tP.gmFillType), tP)
		{
			molecularDynamicsCode->gauge_force_tlsym_device( gaugefieldBuffer, gaugemomentumBuffer);
			calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionTester : public MolecularDynamicsTester
{
	FFermionTester(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_fermion", pC, calculateReferenceValues_FFermion(tP.latticeExtents.getLatticeVolume(), tP.gaugeFillType, tP.gmFillType, tP.spinorFillType, tP.massParameters), tP)
		{
			NonEvenOddSpinorfieldCreator sf(tP.latticeExtents);
			const hardware::buffers::Plain<spinor> in1(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
			const hardware::buffers::Plain<spinor> in2(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);

			in1.load(sf.createSpinorfield(tP.spinorFillType) );
			in2.load(sf.createSpinorfield(tP.spinorFillType) );

			molecularDynamicsCode->fermion_force_device( &in1, &in2, gaugefieldBuffer, gaugemomentumBuffer, tP.massParameters.kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionEvenOddTester : public MolecularDynamicsTester
{
	FFermionEvenOddTester(const ParameterCollection pC, const EvenOddMolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_fermion_eo", pC, calculateReferenceValues_FFermionEvenOdd(tP.latticeExtents.getLatticeVolume(), tP.gaugeFillType, tP.gmFillType, tP.spinorFillType, tP.massParameters), tP)
		{
			EvenOddSpinorfieldCreator sf(tP.latticeExtents);
			const hardware::buffers::Spinor in1(tP.latticeExtents, MolecularDynamicsTester::device);
			const hardware::buffers::Spinor in2(tP.latticeExtents, MolecularDynamicsTester::device);
			sf.fillTwoSpinorBuffers(&in1, tP.spinorFillType, &in2, tP.spinorFillType);// the same SpinorFillType is used for both spinors

			molecularDynamicsCode->fermion_force_eo_device( &in1, &in2, gaugefieldBuffer, gaugemomentumBuffer, tP.evenOrOdd,  tP.massParameters.kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
		}
};

struct FFermionEvenOddComparator : public MolecularDynamicsTester
{
	FFermionEvenOddComparator(const ParameterCollection pC, const MolecularDynamicsTestParameters tP) :
		MolecularDynamicsTester("f_fermion compare even-odd and non-even-odd", pC, defaultReferenceValues(), tP)
		{
			createBuffers(tP);
			fillBuffers(tP);

			hmc_float cpu_res_eo;
			molecularDynamicsCode->fermion_force_eo_device(inEo1, inEo4, gaugefieldBuffer, outEo, ODD, tP.massParameters.kappa );
			MolecularDynamicsTester::code->set_float_to_gaugemomentum_squarenorm_device(outEo, MolecularDynamicsTester::doubleBuffer);
			MolecularDynamicsTester::doubleBuffer->dump(&cpu_res_eo);

			molecularDynamicsCode->fermion_force_eo_device(inEo2, inEo3, gaugefieldBuffer, outEo, EVEN, tP.massParameters.kappa );
			MolecularDynamicsTester::code->set_float_to_gaugemomentum_squarenorm_device(outEo, MolecularDynamicsTester::doubleBuffer);
			MolecularDynamicsTester::doubleBuffer->dump(&cpu_res_eo);
			logger.info() << "|force_eo (even) + force_eo (odd)|^2:";
			logger.info() << cpu_res_eo;

			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(outEo, 0);

			molecularDynamicsCode->fermion_force_device( inNonEo1, inNonEo2, gaugefieldBuffer, outNonEo, tP.massParameters.kappa);
			MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(outNonEo, 1);
			logger.info() << "|force|^2:";
			logger.info() << kernelResult[1];

			BOOST_REQUIRE_EQUAL(kernelResult[0], kernelResult[1]);
			kernelResult.clear();
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
	void createBuffers(const MolecularDynamicsTestParameters tP)
	{
		outNonEo = new const hardware::buffers::Gaugemomentum(tP.latticeExtents, MolecularDynamicsTester::device);
		outEo = new const hardware::buffers::Gaugemomentum(tP.latticeExtents, MolecularDynamicsTester::device);
		inNonEo1 = new const hardware::buffers::Plain<spinor>(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inNonEo2 = new const hardware::buffers::Plain<spinor>(calculateSpinorfieldSize(tP.latticeExtents), MolecularDynamicsTester::device);
		inEo1 = new const hardware::buffers::Spinor(tP.latticeExtents, MolecularDynamicsTester::device);
		inEo2 = new const hardware::buffers::Spinor(tP.latticeExtents, MolecularDynamicsTester::device);
		inEo3 = new const hardware::buffers::Spinor(tP.latticeExtents, MolecularDynamicsTester::device);
		inEo4 = new const hardware::buffers::Spinor(tP.latticeExtents, MolecularDynamicsTester::device);
	}
	void fillBuffers(const MolecularDynamicsTestParameters tP)
	{
		GaugemomentumCreator gm(tP.latticeExtents);
		MolecularDynamicsTester::code->importGaugemomentumBuffer(outNonEo, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(tP.gmFillType) ));
		MolecularDynamicsTester::code->importGaugemomentumBuffer(outEo, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(tP.gmFillType) ));

		NonEvenOddSpinorfieldCreator sf(tP.latticeExtents);
		EvenOddSpinorfieldCreator sfEo(tP.latticeExtents);
		inNonEo1->load(sf.createSpinorfield(tP.spinorFillType));
		inNonEo2->load(sf.createSpinorfield(tP.spinorFillType));
		sfEo.fillTwoSpinorBuffers(inEo1, tP.spinorFillType, inEo2, tP.spinorFillType);
		sfEo.fillTwoSpinorBuffers(inEo3, tP.spinorFillType, inEo4, tP.spinorFillType);

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

struct FFermionStaggeredEvenOddTester : public MolecularDynamicsTester
{
	FFermionStaggeredEvenOddTester(const ParameterCollection pC, const EvenOddMolecularDynamicsTestParameters tP) :
	   MolecularDynamicsTester("f_staggered_fermion_eo", pC, calculateReferenceValues_FFermionStaggeredEvenOdd(tP.latticeExtents.getLatticeVolume(), tP.gaugeFillType, tP.gmFillType, tP.spinorFillType, tP.evenOrOdd), tP)
	{
		EvenOddSpinorStaggeredfieldCreator ssf(tP.latticeExtents);
		const hardware::buffers::SU3vec in1(tP.latticeExtents, MolecularDynamicsTester::device);
		in1.load(ssf.createSpinorfield(tP.spinorFillType));
		const hardware::buffers::SU3vec in2(tP.latticeExtents, MolecularDynamicsTester::device);
		in2.load(ssf.createSpinorfield(tP.spinorFillType));

		molecularDynamicsCode->fermion_staggered_partial_force_device(gaugefieldBuffer, &in1, &in2, gaugemomentumBuffer, tP.evenOrOdd);

		MolecularDynamicsTester::calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
/*
print_staggeredfield_eo_to_textfile("ref_vec_f_stagg1_eo", createSpinorfield(spinorfieldEvenOddElements, 123));
logger.info() << "Produced the ref_vec_f_stagg1_eo text file with the staggered field for the ref. code.";
print_staggeredfield_eo_to_textfile("ref_vec_f_stagg2_eo", createSpinorfield(spinorfieldEvenOddElements, 456));
logger.info() << "Produced the ref_vec_f_stagg2_eo text file with the staggered field for the ref. code.";
*/
		}
};

template<typename TesterClass>
void callTest( const LatticeExtents lE, GaugefieldFillType fillType, GaugeMomentumFilltype gmFillType, const hmc_float epsIn)
{
	GaugeFieldUpdateTestParameters parametersForThisTest(lE, fillType, gmFillType, epsIn);
	hardware::code::OpenClKernelParametersMockup kernelParameters(lE);
	hardware::HardwareParametersMockup hardwareParameters(lE);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

template<typename TesterClass>
void callTest( const LatticeExtents lE, GaugefieldFillType fillType, GaugeMomentumFilltype gmFillType)
{
	MolecularDynamicsTestParameters parametersForThisTest(lE, fillType, gmFillType);
	hardware::code::OpenClKernelParametersMockupForMolecularDynamics kernelParameters(lE, common::action::tlsym);
	hardware::HardwareParametersMockup hardwareParameters(lE);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass tester(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents lE, GaugefieldFillType gfFillType, GaugeMomentumFilltype gmFillType, SpinorFillType sfFillType, WilsonMassParameters kappa )
{
	MolecularDynamicsTestParameters parametersForThisTest(lE, gfFillType, gmFillType, sfFillType, kappa);
	hardware::HardwareParametersMockup hardwareParameters(lE, true);
	hardware::code::OpenClKernelParametersMockupForMolecularDynamics kernelParameters(lE, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents lE, GaugefieldFillType gfFillType, GaugeMomentumFilltype gmFillType, SpinorFillType sfFillType, const bool evenOrOdd, WilsonMassParameters kappa)
{
	EvenOddMolecularDynamicsTestParameters parametersForThisTest(lE, gfFillType, gmFillType, sfFillType, kappa, evenOrOdd);
	hardware::HardwareParametersMockup hardwareParameters(lE, true);
	hardware::code::OpenClKernelParametersMockupForMolecularDynamics kernelParameters(lE, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents lE, GaugefieldFillType gfFillType, GaugeMomentumFilltype gmFillType, SpinorFillType sfFillType, const bool evenOrOdd)
{
	EvenOddMolecularDynamicsTestParameters parametersForThisTest(lE, gfFillType, gmFillType, sfFillType, evenOrOdd);
	hardware::HardwareParametersMockup hardwareParameters(lE, true);
	hardware::code::OpenClKernelParametersMockupForMolecularDynamicsStaggered kernelParameters(lE);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testGaugefieldUpdate(const LatticeExtents lE, const GaugefieldFillType fillType, const GaugeMomentumFilltype gmFillType, const hmc_float eps)
{
	callTest<GaugefieldUpdateTester>(lE, fillType, gmFillType, eps);
}
void testFGauge(const LatticeExtents lE, GaugefieldFillType gfFillType, GaugeMomentumFilltype gmFillType)
{
	callTest<FGaugeTester>(lE, gfFillType, gmFillType);
}
void testFGaugeTlsym(const LatticeExtents lE, GaugefieldFillType gfFillType, GaugeMomentumFilltype gmFillType)
{
	callTest<FGaugeTlsymTester>(lE, gfFillType, gmFillType);
}
void testNonEvenOddFermionForce(const LatticeExtents lE, GaugefieldFillType gfFillType, GaugeMomentumFilltype gmFillType, SpinorFillType sfFillType, WilsonMassParameters kappa )
{
	callTest<FFermionTester>(lE, gfFillType, gmFillType, sfFillType, kappa );
}
void testEvenOddFermionForce(const LatticeExtents lE, GaugefieldFillType fillType, GaugeMomentumFilltype gmFillType, SpinorFillType sfFillType, const bool evenOrOdd, WilsonMassParameters kappa )
{
	callTest<FFermionEvenOddTester>(lE, fillType, gmFillType, sfFillType, evenOrOdd, kappa );
}
void compareEvenOddAndNonEvenOddFermionForce(const LatticeExtents lE, GaugefieldFillType fillType, GaugeMomentumFilltype gmFillType, SpinorFillType sfFillType, WilsonMassParameters kappa)
{
	callTest<FFermionEvenOddComparator>(lE, fillType, gmFillType, sfFillType, kappa );
}
void testEvenOddStaggeredFermionForce(const LatticeExtents lE, GaugefieldFillType fillType, GaugeMomentumFilltype gmFillType, SpinorFillType sfFillType, const bool evenOrOdd )
{
	callTest<FFermionStaggeredEvenOddTester>(lE, fillType, gmFillType, sfFillType, evenOrOdd );
}

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

	//Note: test5 (but this refers to the old input files): "THIS TEST HAS TO BE INVESTIGATED!! IT COULD HINT TO AN ERROR IN THE FUNCTION!"
	BOOST_AUTO_TEST_CASE(GF_UPDATE_1 )
	{
		testGaugefieldUpdate(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, .0);
	}
	BOOST_AUTO_TEST_CASE(GF_UPDATE_2 )
	{
		testGaugefieldUpdate(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, .12);
	}
	BOOST_AUTO_TEST_CASE(GF_UPDATE_3 )
	{
		testGaugefieldUpdate(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, .0);
	}
	BOOST_AUTO_TEST_CASE(GF_UPDATE_4 )
	{
		testGaugefieldUpdate(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, .12);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE )

	BOOST_AUTO_TEST_CASE( F_GAUGE_1 )
	{
		testFGauge(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One);
	}
	BOOST_AUTO_TEST_CASE( F_GAUGE_2 )
	{
		testFGauge(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending);
	}

	BOOST_AUTO_TEST_CASE( F_GAUGE_3 )
	{
		testFGauge(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One);
	}
	BOOST_AUTO_TEST_CASE( F_GAUGE_4 )
	{
		testFGauge(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE_TLSYM )

	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_1 )
	{
		testFGaugeTlsym(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One);
	}
	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_2 )
	{
		testFGaugeTlsym(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending);
	}

	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_3 )
	{
		testFGaugeTlsym(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One);
	}
	BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_4 )
	{
		testFGaugeTlsym(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( STOUT_SMEAR_FERMION_FORCE )

	BOOST_AUTO_TEST_CASE(STOUT_SMEAR_FERMION_FORCE_1 )
	{
		BOOST_TEST_MESSAGE("NOT YET IMPLEMENTED!!");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION )

	BOOST_AUTO_TEST_CASE( F_FERMION_1 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::one, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_2 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::ascendingComplex, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_3 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending, SpinorFillType::one, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_4 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_5 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One, SpinorFillType::one, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_6 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One, SpinorFillType::ascendingComplex, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_7 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::one, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_8 )
	{
		testNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, WilsonMassParameters{nonTrivialParameter} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_EO )

	BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::one, EVEN, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_2 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::ascendingComplex, ODD, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_3 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending, SpinorFillType::one, EVEN, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_4 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, ODD, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_5 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One, SpinorFillType::one, EVEN, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_6 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One, SpinorFillType::ascendingComplex, ODD, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_7 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::one, EVEN, WilsonMassParameters{nonTrivialParameter} );
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_EO_8 )
	{
		testEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, ODD, WilsonMassParameters{nonTrivialParameter} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_COMPARE_NONEO_EO )

	BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_1 )
	{
		compareEvenOddAndNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::one, WilsonMassParameters{nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_2 )
	{
		compareEvenOddAndNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, WilsonMassParameters{nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_3 )
	{
		compareEvenOddAndNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::one, WilsonMassParameters{nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_4 )
	{
		compareEvenOddAndNonEvenOddFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, WilsonMassParameters{nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_STAGG_FERMION_EO )

	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_1 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::one, EVEN );
	}
	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_2 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::One, SpinorFillType::ascendingComplex, ODD );
	}
	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_3 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending, SpinorFillType::one, EVEN );
	}
	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_4 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, ODD );
	}
	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_5 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One, SpinorFillType::one, ODD );
	}
	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_6 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::One, SpinorFillType::ascendingComplex, EVEN );
	}
	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_7 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::one, EVEN );
	}
	BOOST_AUTO_TEST_CASE( F_STAGG_FERMION_EO_8 )
	{
		testEvenOddStaggeredFermionForce(LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, GaugeMomentumFilltype::Ascending, SpinorFillType::ascendingComplex, EVEN );
	}

BOOST_AUTO_TEST_SUITE_END()


