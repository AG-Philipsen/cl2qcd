/*
 * Copyright 2015 Christopher Pinke, Francesca Cuteri
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

# pragma once

#include "hardwareParameters.hpp"
#include "openClKernelParameters.hpp"
#include "hardwareTestUtilities.hpp"
#include "../geometry/latticeExtents.hpp"

namespace hardware
{
	class HardwareParametersMockup : public HardwareParametersInterface
	{
	public:
		HardwareParametersMockup(const int nsIn, const int ntIn) : ns(nsIn), nt(ntIn), useEvenOdd(false)
		{
			setGpuAndCpuOptions(checkBoostRuntimeArgumentsForGpuUsage());
		};
		HardwareParametersMockup(const int nsIn, const int ntIn, const bool useEvenOddIn) : ns(nsIn), nt(ntIn), useEvenOdd(useEvenOddIn)
		{
			setGpuAndCpuOptions(checkBoostRuntimeArgumentsForGpuUsage());
		};
		HardwareParametersMockup(LatticeExtents lE): ns(lE.getNs()), nt(lE.getNt()), useEvenOdd(false)
		{
			setGpuAndCpuOptions(checkBoostRuntimeArgumentsForGpuUsage());
		};
		HardwareParametersMockup(LatticeExtents lE, const bool useEvenOddIn): ns(lE.getNs()), nt(lE.getNt()), useEvenOdd(useEvenOddIn)
		{
			setGpuAndCpuOptions(checkBoostRuntimeArgumentsForGpuUsage());
		};
		~HardwareParametersMockup() {};
		virtual int getNs() const override
		{
			return ns;
		}
		virtual int getNt() const override
		{
			return nt;
		}
		virtual bool disableOpenCLCompilerOptimizations() const override
		{
			return false;
		}
		virtual bool useGpu() const override
		{
			return useGpuValue;
		}
		virtual bool useCpu() const override
		{
			return useCpuValue;
		}
		virtual int getMaximalNumberOfDevices() const override
		{
			return 1;
		}
		virtual std::vector<int> getSelectedDevices() const override
		{
			return std::vector<int>{0};
		}
		virtual bool splitCpu() const override
		{
			return false;
		}
		virtual bool enableProfiling() const override
		{
			return false;
		}
		virtual bool useSameRandomNumbers() const override
		{
			return false;
		}
		virtual bool useEvenOddPreconditioning() const override
		{
			return useEvenOdd;
		}
		virtual int getSpatialLatticeVolume() const override
		{
			return getNs() * getNs() * getNs();
		}
		virtual int getLatticeVolume() const override
		{
			return getNs() * getNs() * getNs() * getNt();
		}
	private:
		const int ns, nt;
		const bool useEvenOdd;
		bool useGpuValue, useCpuValue;
		void setGpuAndCpuOptions(const bool value)
		{
			useGpuValue = value;
			useCpuValue = !value;
		}
	};

	struct HardwareParametersMockupForDeviceSelection : public HardwareParametersMockup
	{
		HardwareParametersMockupForDeviceSelection(const int ns, const int nt, const int maximalNumberOfDevices, const std::vector<int> selectedDevices) :
			HardwareParametersMockup(ns, nt), maximalNumberOfDevices(maximalNumberOfDevices), selectedDevices(selectedDevices)
		{

		}
		virtual int getMaximalNumberOfDevices() const override
		{
			return maximalNumberOfDevices;
		}
		virtual std::vector<int> getSelectedDevices() const override
		{
			return selectedDevices;
		}
	private:
		const int maximalNumberOfDevices;
		const std::vector<int> selectedDevices;
	};

	struct HardwareParametersMockupWithoutGpus : public HardwareParametersMockup
	{
		HardwareParametersMockupWithoutGpus(const int ns, const int nt) :
			HardwareParametersMockup(ns, nt)
		{}

		virtual bool useGpu() const override
		{
			return false;
		}
		virtual bool useCpu() const override
		{
			return true;
		}
	};

	struct HardwareParametersMockupWithoutCpus : public HardwareParametersMockup
	{
		HardwareParametersMockupWithoutCpus(const int ns, const int nt) :
			HardwareParametersMockup(ns, nt)
		{}

		virtual bool useGpu() const override
		{
			return true;
		}
		virtual bool useCpu() const override
		{
			return false;
		}
	};

	struct HardwareParametersMockupWithProfiling : public HardwareParametersMockup
	{
		HardwareParametersMockupWithProfiling(const int ns, const int nt) :
			HardwareParametersMockup(ns, nt)
		{}

		virtual bool enableProfiling() const override
		{
			return true;
		}
	};

}

namespace hardware
{
	namespace code {
		class OpenClKernelParametersMockup : public OpenClKernelParametersInterface
		{
		public:
			OpenClKernelParametersMockup(int nsIn, int ntIn):
				ns(nsIn), nt(ntIn), rhoIter(0), rho(0.), useRectangles(true), useSmearing(false), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn, int rhoIterIn, double rhoIn, bool useSmearingIn):
				ns(nsIn), nt(ntIn), rhoIter(rhoIterIn), rho(rhoIn), useRectangles(false), useSmearing(useSmearingIn), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn, bool useRectanglesIn):
				ns(nsIn), nt(ntIn), rhoIter(0), rho(0.), useRectangles(useRectanglesIn), useSmearing(false), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			OpenClKernelParametersMockup(LatticeExtents lE):
				ns(lE.getNs()), nt(lE.getNt()), rhoIter(0), rho(0.), useRectangles(true), useSmearing(false), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			OpenClKernelParametersMockup(LatticeExtents lE, bool useRectanglesIn):
				ns(lE.getNs()), nt(lE.getNt()), rhoIter(0), rho(0.), useRectangles(useRectanglesIn), useSmearing(false), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			~OpenClKernelParametersMockup()	{};
			virtual int getNs() const override
			{
				return ns;
			}
			virtual int getNt() const override
			{
				return nt;
			}
			virtual size_t getPrecision() const override
			{
				return 64;
			}
			virtual bool getUseChemPotRe() const override
			{
				return false;
			}
			virtual bool getUseChemPotIm() const override
			{
				return false;
			}
			virtual bool getUseSmearing() const override
			{
				return useSmearing;
			}
			virtual double getChemPotRe() const override
			{
				return 0.;
			}
			virtual double getChemPotIm() const override
			{
				return 0.;
			}
			virtual double getRho() const override
			{
				return rho;
			}
			virtual int getRhoIter() const override
			{
				return rhoIter;
			}
			virtual bool getUseRec12() const override
			{
				return useRec12Value;
			}
			virtual bool getUseEo() const override
			{
				return false;
			}
			virtual common::action  getFermact() const override
			{
				return common::action::wilson;
			}
			virtual common::action  getGaugeact() const override
			{
				return common::action::wilson;
			}
			virtual int getMetroApproxOrd() const override
			{
				return 15.;
			}
			virtual int getMdApproxOrd() const override
			{
				return 8.;
			}
			virtual double getThetaFermionSpatial() const override
			{
				return 0.;
			}
			virtual double getThetaFermionTemporal() const override
			{
				return 0.;
			}
			virtual double getBeta() const override
			{
				return 5.69;
			}
			virtual double getKappa() const override
			{
				return 0.125;
			}
			virtual int getNumSources() const override
			{
				return 12;
			}
			virtual common::sourcecontents getSourceContent() const override
			{
				return common::sourcecontents::one;
			}
			virtual common::sourcetypes getSourceType() const override
			{
				return common::sourcetypes::point;
			}
			virtual bool getUseAniso() const override
			{
				return false;
			}
			virtual size_t getSpatialLatticeVolume() const override
			{
				return getNs() * getNs() * getNs();
			}
			virtual size_t getLatticeVolume() const override
			{
				return getNs() * getNs() * getNs() * getNt();
			}
			virtual bool getUseRectangles() const override
			{
				return useRectangles;
			}
			virtual double getC0() const override
			{
				return 1.;
			}
			virtual double getC1() const override
			{
				return 0.;
			}
			virtual double getXi0() const override
			{
				return 0.;
			}
			virtual size_t getFloatSize() const override
			{
				return getPrecision() / 8;
			}
			virtual size_t getMatSize() const override
			{
				// TODO with rec12 this becomes 6
				return 9;
			}
			virtual size_t getSpinorFieldSize() const override
			{
				return getLatticeVolume();
			}
			virtual size_t getEoprecSpinorFieldSize() const override
			{
				return getLatticeVolume() / 2;
			}
			virtual int getCorrDir() const override
			{
				return 3;
			}
			virtual bool getMeasureCorrelators() const override
			{
				return true;
			}
			virtual bool getUseMergeKernelsFermion() const override
			{
				return false;
			}
			virtual hmc_float getMuBar() const override
			{
				return 2* getKappa() * getMu();
			}
			virtual double getMass() const override
			{
				return 0.1;
			}
			virtual bool getUseSameRndNumbers() const override
			{
				return false;
			}
			virtual uint32_t getHostSeed() const override
			{
				return 4815;
			}
			virtual bool getUseMergeKernelsSpinor() const override
			{
				return false;
			}
			virtual double getMu() const override
			{
				return 0.006;
			}
			virtual double getApproxLower() const override
			{
				return 1.e-5;
			}
		protected:
			int ns, nt, rhoIter;
			double rho;
			bool useRectangles, useSmearing;
			bool useRec12Value;
		};

		class OpenClKernelParametersMockupForSpinorTests : public OpenClKernelParametersMockup
		{
		public:
			OpenClKernelParametersMockupForSpinorTests(const int nsIn, const int ntIn) :
				OpenClKernelParametersMockup(nsIn, ntIn, false), prec(64), useEvenOdd(false), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			OpenClKernelParametersMockupForSpinorTests(const int nsIn, const int ntIn, const bool useEvenOddIn) :
				OpenClKernelParametersMockup(nsIn, ntIn, false), prec(64), useEvenOdd(useEvenOddIn), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			OpenClKernelParametersMockupForSpinorTests(const LatticeExtents lE) :
				OpenClKernelParametersMockup(lE, false), prec(64), useEvenOdd(false), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			OpenClKernelParametersMockupForSpinorTests(const LatticeExtents lE, const bool useEvenOddIn) :
				OpenClKernelParametersMockup(lE, false), prec(64), useEvenOdd(useEvenOddIn), useRec12Value(checkBoostRuntimeArgumentsForRec12Usage()) {};
			~OpenClKernelParametersMockupForSpinorTests()	{};

			virtual size_t getPrecision() const override
			{
				return prec;
			}
			virtual bool getUseSmearing() const override
			{
				return false;
			}
			virtual double getRho() const override
			{
				return 0.;
			}
			virtual int getRhoIter() const override
			{
				return 0;
			}
			virtual bool getUseRec12() const override
			{
				return false;
			}
			virtual bool getUseEo() const override
			{
				return useEvenOdd;
			}
		protected:
			const size_t prec;
			const bool useEvenOdd;
			bool useRec12Value;
		};

		class OpenClKernelParametersMockupForSpinorStaggered : public OpenClKernelParametersMockupForSpinorTests
		{
		public:
			OpenClKernelParametersMockupForSpinorStaggered(const LatticeExtents lE) :
				OpenClKernelParametersMockupForSpinorTests(lE), fermact(common::action::rooted_stagg)
		{
				fermact = common::action::rooted_stagg;
		};
			OpenClKernelParametersMockupForSpinorStaggered(const LatticeExtents lE, const bool useEvenOddIn) :
				OpenClKernelParametersMockupForSpinorTests(lE, useEvenOddIn), fermact(common::action::rooted_stagg)
		{
				fermact = common::action::rooted_stagg;
		};
		virtual common::action getFermact() const override
		{
			return fermact;
		}
		common::action fermact;
		};

		class OpenClKernelParametersMockupForCorrelators : public OpenClKernelParametersMockupForSpinorTests
		{
		public:
			OpenClKernelParametersMockupForCorrelators(const int nsIn, const int ntIn, const double kappaIn, const double directionIn) :
				OpenClKernelParametersMockupForSpinorTests(nsIn, ntIn), correlatorDirection(directionIn), kappa(kappaIn) {};

			virtual int getCorrDir() const override
			{
				return correlatorDirection;
			}
			virtual double getKappa() const override
			{
				return kappa;
			}
			virtual bool getMeasureCorrelators() const override
			{
				return true;
			}
			const int correlatorDirection;
			const double kappa;
		};

		class OpenClKernelParametersMockupForStaggeredCorrelators : public OpenClKernelParametersMockupForSpinorStaggered
			{
			public:
				OpenClKernelParametersMockupForStaggeredCorrelators(const LatticeExtents lE, const double kappaIn, const double directionIn) :
					OpenClKernelParametersMockupForSpinorStaggered(lE, true), correlatorDirection(directionIn), kappa(kappaIn) {};

				virtual int getCorrDir() const override
				{
					return correlatorDirection;
				}
				virtual double getKappa() const override
				{
					return kappa;
				}
				virtual bool getMeasureCorrelators() const override
				{
					return true;
				}
				const int correlatorDirection;
				const double kappa;
			};

		struct OpenClKernelParametersMockupForSourceTests : public OpenClKernelParametersMockupForSpinorTests
		{
			OpenClKernelParametersMockupForSourceTests(const int nsIn, const int ntIn, const common::sourcecontents sC, const common::sourcetypes sT) :
				OpenClKernelParametersMockupForSpinorTests(nsIn, ntIn), sC(sC), sT(sT) {};

			virtual common::sourcecontents getSourceContent() const override
			{
				return sC;
			}
			virtual common::sourcetypes getSourceType() const override
			{
				return sT;
			}
			virtual bool getMeasureCorrelators() const override
			{
				return false;
			}
			const common::sourcecontents sC;
			const common::sourcetypes sT;
		};

		struct OpenClKernelParametersMockupForStaggeredSourceTests : public OpenClKernelParametersMockupForSpinorStaggered
		{
			OpenClKernelParametersMockupForStaggeredSourceTests(const LatticeExtents lE, const common::sourcecontents sC,  const common::sourcetypes sT) :
				OpenClKernelParametersMockupForSpinorStaggered(lE, true), sC(sC), sT(sT) {};

			virtual common::sourcecontents getSourceContent() const override
			{
				return sC;
			}
			virtual common::sourcetypes getSourceType() const override
			{
				return sT;
			}
			virtual bool getMeasureCorrelators() const override
			{
				return false;
			}
			const common::sourcecontents sC;
			const common::sourcetypes sT;
		};

		class OpenClKernelParametersMockupForTwistedMass : public OpenClKernelParametersMockupForSpinorTests
		{
		public:
			OpenClKernelParametersMockupForTwistedMass(int nsIn, int ntIn, const bool needEvenOddIn, const bool useMergedKernels= false) :
				OpenClKernelParametersMockupForSpinorTests(nsIn, ntIn, needEvenOddIn), fermact(common::action::wilson), useMergedKernels(useMergedKernels)
			{
				//NOTE: for the moment, these member are set here in order to overwrite the settings from the parent class, but this should be done nicer!
				fermact = common::action::twistedmass;
				//todo: kappa and mu should be set to 0 or so as they should not be used in the test
			}
			virtual common::action  getFermact() const override
			{
				return fermact;
			}
			virtual bool getUseMergeKernelsFermion() const override
			{
				return useMergedKernels;
			}
			common::action fermact;
			const bool useMergedKernels;
		};
		class OpenClKernelParametersMockupForDslashEvenOdd : public OpenClKernelParametersMockupForSpinorTests
		{
		public:
			OpenClKernelParametersMockupForDslashEvenOdd(int nsIn, int ntIn, const bool needEvenOddIn, const double kappaIn) :
				OpenClKernelParametersMockupForSpinorTests(nsIn, ntIn, needEvenOddIn), kappa(kappaIn), thetaFermionTemporal(0.), thetaFermionSpatial(0.),
				useChemPotIm(false), useChemPotRe(false), chemPot({0.,0.}) {}
			OpenClKernelParametersMockupForDslashEvenOdd(int nsIn, int ntIn, const bool needEvenOddIn, const double kappaIn, const double thetaTIn, const double thetaSIn, hmc_complex chemPotIn) :
				OpenClKernelParametersMockupForSpinorTests(nsIn, ntIn, needEvenOddIn), kappa(kappaIn), thetaFermionTemporal(thetaTIn), thetaFermionSpatial(thetaSIn),
				useChemPotIm(true), useChemPotRe(true), chemPot(chemPotIn) {}

			virtual double getKappa() const override
			{
				return kappa;
			}
			virtual double  getThetaFermionTemporal() const override
			{
				return thetaFermionTemporal;
			}
			virtual double getThetaFermionSpatial() const override
			{
				return thetaFermionSpatial;
			}
			virtual bool getUseChemPotIm() const override
			{
				return useChemPotIm;
			}
			virtual double getChemPotIm() const override
			{
				return chemPot.im;
			}
			virtual bool getUseChemPotRe() const override
			{
				return useChemPotRe;
			}
			virtual double getChemPotRe() const override
			{
				return chemPot.re;
			}

			double kappa;
			double thetaFermionTemporal;
			double thetaFermionSpatial;
			bool useChemPotIm, useChemPotRe;
			hmc_complex chemPot;
		};

		class OpenClKernelParametersMockupForFermionsStaggeredTests : public OpenClKernelParametersMockupForSpinorStaggered
		{
		public:
			OpenClKernelParametersMockupForFermionsStaggeredTests(LatticeExtents lE, const double massIn, const bool useEvenOddIn, const double thetaTIn , const double thetaSIn , const double ImChemPot = 0.):
				OpenClKernelParametersMockupForSpinorStaggered(lE,useEvenOddIn), mass(massIn), useEvenOdd(useEvenOddIn), thetaFermionTemporal(thetaTIn), thetaFermionSpatial(thetaSIn), useChemPotIm(true), imaginaryChemicalPotential(ImChemPot) {}
			OpenClKernelParametersMockupForFermionsStaggeredTests(LatticeExtents lE, const bool useEvenOddIn, const double thetaTIn , const double thetaSIn, const double ImChemPot ):
				OpenClKernelParametersMockupForSpinorStaggered(lE,useEvenOddIn), mass(0.), useEvenOdd(useEvenOddIn), thetaFermionTemporal(thetaTIn), thetaFermionSpatial(thetaSIn), useChemPotIm(true), imaginaryChemicalPotential(ImChemPot) {}
			OpenClKernelParametersMockupForFermionsStaggeredTests(LatticeExtents lE, const bool useEvenOddIn):
				OpenClKernelParametersMockupForSpinorStaggered(lE,useEvenOddIn), mass(0.), useEvenOdd(useEvenOddIn), thetaFermionTemporal(0.), thetaFermionSpatial(0.), useChemPotIm(false), imaginaryChemicalPotential(0.) {}

			virtual double getMass() const override
			{
				return mass;
			}
			virtual double getThetaFermionTemporal() const override
			{
				return thetaFermionTemporal;
			}
			virtual double getThetaFermionSpatial() const override
			{
				return thetaFermionSpatial;
			}
			virtual bool getUseEo() const override
			{
				return useEvenOdd;
			}
			virtual double getChemPotIm() const override
			{
				return imaginaryChemicalPotential;
			}
			virtual bool getUseChemPotIm() const override
			{
				return useChemPotIm;
			}

		protected:
			const double mass;
			const bool useEvenOdd;
			const double thetaFermionTemporal;
			const double thetaFermionSpatial;
			bool useChemPotIm;
			const double imaginaryChemicalPotential;
		};

		class OpenClKernelParametersMockupForMergedFermionKernels : public OpenClKernelParametersMockupForSpinorTests
		{
		public:
			OpenClKernelParametersMockupForMergedFermionKernels(int nsIn, int ntIn) :
				OpenClKernelParametersMockupForSpinorTests(nsIn, ntIn, true) {}
			virtual bool getUseMergeKernelsFermion() const override
			{
				return true;
			}
		};
		class OpenClKernelParametersMockupForMergedSpinorKernels : public OpenClKernelParametersMockupForSpinorTests
		{
		public:
			OpenClKernelParametersMockupForMergedSpinorKernels(int nsIn, int ntIn) :
				OpenClKernelParametersMockupForSpinorTests(nsIn, ntIn, true) {}
			virtual bool getUseMergeKernelsSpinor() const override
			{
				return true;
			}
		};
		class OpenClKernelParametersMockupForMolecularDynamics : public OpenClKernelParametersMockup
		{
		public:
			OpenClKernelParametersMockupForMolecularDynamics(LatticeExtents lE):
				OpenClKernelParametersMockup(lE), gaugeact(common::action::wilson), useEvenOdd(false) {}
			OpenClKernelParametersMockupForMolecularDynamics(LatticeExtents lE, const bool useEvenOddIn):
				OpenClKernelParametersMockup(lE), gaugeact(common::action::wilson), useEvenOdd(useEvenOddIn) {}
			OpenClKernelParametersMockupForMolecularDynamics(LatticeExtents lE, const common::action actionIn ) :
				OpenClKernelParametersMockup(lE), gaugeact(actionIn), useEvenOdd(false) {}
			virtual common::action  getGaugeact() const override
			{
				return gaugeact;
			}
			virtual double getC1() const override
			{
				return -0.083333333;
			}
			virtual bool getUseEo() const override
			{
				return useEvenOdd;
			}

		protected:
			common::action gaugeact;
			const bool useEvenOdd;
		};

		class OpenClKernelParametersMockupForMolecularDynamicsStaggered : public OpenClKernelParametersMockupForMolecularDynamics
		{
		public:
			OpenClKernelParametersMockupForMolecularDynamicsStaggered(LatticeExtents lE) :
				OpenClKernelParametersMockupForMolecularDynamics(lE), fermact(common::action::rooted_stagg)
			{
					fermact = common::action::rooted_stagg;
			};
			virtual common::action getFermact() const override
			{
				return fermact;
			}
			virtual bool getUseEo() const override
			{
				return true;
			}
		protected:
			common::action fermact;
		};
	}
}
