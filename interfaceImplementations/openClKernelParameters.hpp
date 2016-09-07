# pragma once

#include "../meta/inputparameters.hpp"
#include "../meta/util.hpp"
#include "../hardware/openClKernelParameters.hpp"

namespace hardware {
	namespace code {
		class OpenClKernelParametersImplementation final : public OpenClKernelParametersInterface
		{
		public:
			OpenClKernelParametersImplementation( const meta::Inputparameters & parametersIn) : fullParameters(&parametersIn) {};
			~OpenClKernelParametersImplementation() {};
			virtual int getNs() const override
			{
				return fullParameters->get_nspace();
			}
			virtual int getNt() const override
			{
				return fullParameters->get_ntime();
			}
			virtual size_t getPrecision() const override
			{
				return fullParameters->get_precision();
			}
			virtual bool getUseChemPotRe() const override
			{
				return fullParameters->get_use_chem_pot_re();
			}
			virtual bool getUseChemPotIm() const override
			{
				return fullParameters->get_use_chem_pot_im();
			}
			virtual bool getUseSmearing() const override
			{
				return fullParameters->get_use_smearing();
			}
			virtual double getChemPotRe() const override
			{
				return fullParameters->get_chem_pot_re();
			}
			virtual double getChemPotIm() const override
			{
				return fullParameters->get_chem_pot_im();
			}
			virtual double getRho() const override
			{
				return fullParameters->get_rho();
			}
			virtual int getRhoIter() const override
			{
				return fullParameters->get_rho_iter();
			}
			virtual bool getUseRec12() const override
			{
				return fullParameters->get_use_rec12();
			}
			virtual bool getUseEo() const override
			{
				return fullParameters->get_use_eo();
			}
			virtual common::action  getFermact() const override
			{
				return fullParameters->get_fermact();
			}
			virtual common::action  getGaugeact() const override
			{
				return fullParameters->get_gaugeact();
			}
			virtual int getMetroApproxOrd() const override
			{
				return fullParameters->get_metro_approx_ord();
			}
			virtual int getMdApproxOrd() const override
			{
				return fullParameters->get_md_approx_ord();
			}
			virtual double getThetaFermionSpatial() const override
			{
				return fullParameters->get_theta_fermion_spatial();
			}
			virtual double getThetaFermionTemporal() const override
			{
				return fullParameters->get_theta_fermion_temporal();
			}
			virtual double getBeta() const override
			{
				return fullParameters->get_beta();
			}
			virtual double getKappa() const override
			{
				return fullParameters->get_kappa();
			}
			virtual int getNumSources() const override
			{
				return fullParameters->get_num_sources();
			}
			virtual common::sourcecontents getSourceContent() const override
			{
				return fullParameters->get_sourcecontent();
			}
			virtual common::sourcetypes getSourceType() const override
			{
				return fullParameters->get_sourcetype();
			}
			virtual bool getUseAniso() const override
			{
				return fullParameters->get_use_aniso();
			}
			virtual size_t getSpatialLatticeVolume() const override
			{
				return meta::get_volspace( *fullParameters );
			}
			virtual size_t getLatticeVolume() const override
			{
				return meta::get_vol4d( * fullParameters );
			}
			virtual bool getUseRectangles() const override
			{
				return meta::get_use_rectangles( * fullParameters );
			}
			virtual double getC0() const override
			{
				return meta::get_c0( * fullParameters );
			}
			virtual double getC1() const override
			{
				return meta::get_c1( * fullParameters );
			}
			virtual double getXi0() const override
			{
				return meta::get_xi_0( * fullParameters );
			}
			virtual size_t getFloatSize() const override
			{
				return meta::get_float_size( * fullParameters);
			}
			virtual size_t getMatSize() const override
			{
				return meta::get_mat_size( * fullParameters);
			}
			virtual size_t getSpinorFieldSize() const override
			{
				return meta::get_vol4d( * fullParameters); //check this!
			}
			virtual size_t getEoprecSpinorFieldSize() const override
			{
				return meta::get_vol4d( * fullParameters) / 2; //check this!
			}
			virtual int getCorrDir() const override
			{
				return fullParameters->get_corr_dir();
			}
			virtual bool getMeasureCorrelators() const override
			{
				return fullParameters->get_measure_correlators();
			}
			virtual bool getUseMergeKernelsFermion() const override
			{
				return fullParameters->get_use_merge_kernels_fermion();
			}
			virtual hmc_float getMuBar() const override
			{
				return meta::get_mubar( * fullParameters);
			}
			virtual double getMass() const override
			{
				return fullParameters->get_mass();
			}
			virtual bool getUseSameRndNumbers() const override
			{
				return fullParameters->get_use_same_rnd_numbers();
			}
			virtual uint32_t getHostSeed() const override
			{
				return fullParameters->get_host_seed();
			}
			virtual bool getUseMergeKernelsSpinor() const override
			{
				return fullParameters->get_use_merge_kernels_spinor();
			}
			virtual double getMu() const override
			{
				return fullParameters->get_mu();
			}
			virtual double getApproxLower() const override
			{
				return fullParameters->get_approx_lower();
			}
		private:
			const meta::Inputparameters * fullParameters;
		};
	}
}
