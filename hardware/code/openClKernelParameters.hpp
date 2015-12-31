# pragma once

namespace hardware {
	namespace code {
		class OpenClKernelParametersInterface
		{
		public:
			virtual ~OpenClKernelParametersInterface(){};
			virtual int getNs() const = 0;
			virtual int getNt() const = 0;
			virtual size_t getPrecision() const = 0;
			virtual bool getUseChemPotRe() const = 0;
			virtual double getChemPotRe() const = 0;
			virtual bool getUseChemPotIm() const = 0;
			virtual double getChemPotIm() const = 0;
			virtual bool getUseSmearing() const = 0;
			virtual double getRho() const = 0;
			virtual int getRhoIter() const = 0;
			virtual bool getUseRec12() const = 0;
			virtual bool getUseEo() const = 0;
			virtual common::action getFermact() const = 0;
			virtual int getMetroApproxOrd() const = 0;
			virtual int getMdApproxOrd() const = 0;
			virtual double getThetaFermionSpatial() const = 0;
			virtual double getThetaFermionTemporal() const = 0;
			virtual double getBeta() const = 0;
			virtual double getKappa() const = 0;
			virtual int getNumSources() const = 0;
			virtual common::sourcecontents getSourceContent() const = 0;
			virtual common::sourcetypes getSourceType() const = 0;
			virtual bool getUseAniso() const = 0;
			virtual size_t getSpatialLatticeVolume() const = 0;
			virtual size_t getLatticeVolume() const = 0;
			virtual bool getUseRectangles() const = 0;
			virtual double getC0() const = 0;
			virtual double getC1() const = 0;
			virtual double getXi0() const = 0;
			virtual size_t getFloatSize() const = 0;
			virtual size_t getMatSize() const = 0;
			virtual size_t getSpinorFieldSize() const = 0;
			virtual size_t getEoprecSpinorFieldSize() const = 0;
			virtual int getCorrDir() const = 0;
			virtual bool getMeasureCorrelators() const = 0;
			virtual bool getUseMergeKernelsFermion() const = 0;
			virtual bool getUseMergeKernelsSpinor() const = 0;
			virtual hmc_float getMuBar() const = 0;
			virtual double getMass() const = 0;
			virtual bool getUseSameRndNumbers() const = 0;
			virtual uint32_t getHostSeed() const = 0;
			virtual double getMu() const = 0;
			virtual double getApproxLower() const = 0;
		};
	}
}

