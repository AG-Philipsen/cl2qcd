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

#include "../hardwareParameters.hpp"

namespace hardware
{
class HardwareParametersMockup : public HardwareParametersInterface
{
public:
	HardwareParametersMockup(const int nsIn, const int ntIn) : ns(nsIn), nt(ntIn) {};
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
		return false;
	}
	virtual bool useCpu() const override
	{
		return true;
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
		return false;
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
	int ns, nt;
};
}
#include "../openClCode.hpp"

namespace hardware
{
	class OpenClCodeMockup : public OpenClCode
	{
	public:
		OpenClCodeMockup(const hardware::code::OpenClKernelParametersInterface & kernelParams) :
			kernelParameters(& kernelParams) {};
		~OpenClCodeMockup() {}
		virtual std::unique_ptr<const hardware::code::Real> getCode_real(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Real>( new hardware::code::Real{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Gaugefield> getCode_gaugefield(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Gaugefield>( new hardware::code::Gaugefield{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Prng> getCode_PRNG(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Prng>( new hardware::code::Prng{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Complex> getCode_complex(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Complex>( new hardware::code::Complex{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Spinors> getCode_Spinors(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Spinors>( new hardware::code::Spinors{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Fermions> getCode_Fermions(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Fermions>( new hardware::code::Fermions{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Gaugemomentum> getCode_Gaugemomentum(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Gaugemomentum>( new hardware::code::Gaugemomentum{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Molecular_Dynamics> getCode_Molecular_Dynamics(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Molecular_Dynamics>( new hardware::code::Molecular_Dynamics{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Correlator> getCode_Correlator(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Correlator>( new hardware::code::Correlator{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Heatbath> getCode_Heatbath(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Heatbath>( new hardware::code::Heatbath{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Kappa> getCode_Kappa(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Kappa>( new hardware::code::Kappa{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Buffer> getCode_Buffer(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Buffer>( new hardware::code::Buffer{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Spinors_staggered> getCode_Spinors_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Spinors_staggered>( new hardware::code::Spinors_staggered{*kernelParameters,deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Correlator_staggered> getCode_Correlator_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Correlator_staggered>( new hardware::code::Correlator_staggered{*kernelParameters, deviceIn} ) ;
		}
		virtual std::unique_ptr<const hardware::code::Fermions_staggered> getCode_Fermions_staggered(hardware::Device * deviceIn) const override
		{
			return std::unique_ptr<const hardware::code::Fermions_staggered>( new hardware::code::Fermions_staggered{*kernelParameters, deviceIn} ) ;
		}
	private:
		const hardware::code::OpenClKernelParametersInterface * kernelParameters;
	};
}

namespace hardware {
	namespace code {
		class OpenClKernelParametersMockup : public OpenClKernelParametersInterface
		{
		public:
			OpenClKernelParametersMockup(int nsIn, int ntIn , size_t precIn = 64) :
				ns(nsIn), nt(ntIn), rhoIter(0), rho(0.), beta(5.69), kappa(0.125), mu(0.006), mass(0.1), fermact(common::action::wilson), useRectangles(false), useSmearing(false), prec(precIn) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn, int rhoIterIn, double rhoIn, bool useSmearingIn , size_t precIn = 64):
				ns(nsIn), nt(ntIn), rhoIter(rhoIterIn), rho(rhoIn), beta(5.69), kappa(0.125), mu(0.006), mass(0.1), fermact(common::action::wilson), useRectangles(false), useSmearing(useSmearingIn), prec(precIn) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn, bool useRectanglesIn , size_t precIn = 64):
				ns(nsIn), nt(ntIn), rhoIter(0), rho(0.), beta(5.69), kappa(0.125), mu(0.006), mass(0.1), fermact(common::action::wilson), useRectangles(useRectanglesIn), useSmearing(false), prec(precIn) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn, int rhoIterIn, double rhoIn, common::action fermactIn = common::action::twistedmass, bool useRectanglesIn=false, bool useSmearingIn=false, size_t precIn = 64) :
				ns(nsIn), nt(ntIn), rhoIter(rhoIterIn), rho(rhoIn), beta(5.69), kappa(0.125), mu(0.006), mass(0.1), fermact(fermactIn), useRectangles(useRectanglesIn), useSmearing(useSmearingIn), prec(precIn) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn , double betaIn, double kappaIn, size_t precIn = 64) :
				ns(nsIn), nt(ntIn), rhoIter(0), rho(0.), beta(betaIn), kappa(kappaIn), mu(0.006), mass(0.1), fermact(common::action::wilson), useRectangles(false), useSmearing(false), prec(precIn) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn , double betaIn, double kappaIn, double rhoIn, double muIn, size_t precIn = 64) :
				ns(nsIn), nt(ntIn), rhoIter(0), rho(rhoIn), beta(betaIn), kappa(kappaIn), mu(muIn), mass(0.1), fermact(common::action::wilson), useRectangles(false), useSmearing(false), prec(precIn) {};
			OpenClKernelParametersMockup(int nsIn, int ntIn , double betaIn, double kappaIn, double rhoIn, double muIn, double massIn, size_t precIn = 64) :
				ns(nsIn), nt(ntIn), rhoIter(0), rho(rhoIn), beta(betaIn), kappa(kappaIn), mu(muIn), mass(massIn), fermact(common::action::wilson), useRectangles(false), useSmearing(false), prec(precIn) {};
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
				return prec;
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
				return false;
			}
			virtual bool getUseEo() const override
			{
				return false;
			}
			virtual common::action  getFermact() const override
			{
				return fermact;
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
				return beta;
			}
			virtual double getKappa() const override
			{
				return kappa;
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
				return prec / 8;
			}
			virtual size_t getMatSize() const override
			{
				// TODO with rec12 this becomes 6
				return 9;
			}
			virtual size_t getSpinorFieldSize() const override
			{
				return getNs() * getNs() * getNs() * getNt();
			}
			virtual size_t getEoprecSpinorFieldSize() const override
			{
				return getNs() * getNs() * getNs() * getNt() / 2;
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
				return 2* getKappa() * 0.006;
			}
			virtual double getMass() const override
			{
				return mass;
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
				return mu;
			}
			virtual double getApproxLower() const override
			{
				return 1.e-5;
			}
		protected:
			int ns, nt, rhoIter;
			double rho, beta, kappa, mu, mass;
			common::action fermact;
			bool useRectangles, useSmearing;
			size_t prec;
		};
	}
}
