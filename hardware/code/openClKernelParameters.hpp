# pragma once

namespace hardware {
	namespace code {
		class OpenClKernelParametersInterface
		{
		public:
			virtual ~OpenClKernelParametersInterface(){};
			virtual int getNs() const = 0;
			virtual int getNt() const = 0;
			virtual int getPrecision() const = 0;
			virtual int getUseChemPotRe() const = 0;
			virtual int getChemPotRe() const = 0;
			virtual int getUseChemPotIm() const = 0;
			virtual int getChemPotIm() const = 0;
			virtual int getUseSmearing() const = 0;
			virtual int getRho() const = 0;
			virtual int getRhoIter() const = 0;
			virtual int getUseRec12() const = 0;
			virtual int getUseEo() const = 0;
			virtual int getFermact() const = 0;
			virtual int getMetroApproxOrd() const = 0;
			virtual int getMdApproxOrd() const = 0;
			virtual int getThetaFermionSpatial() const = 0;
			virtual int getThetaFermionTemporal() const = 0;
			virtual double getBeta() const = 0;
			virtual double getKappa() const = 0;
			virtual int getNumSources() const = 0;
			virtual int getSourceContent() const = 0;
			virtual int getUseAniso() const = 0;
		};
	}
}

#include "../meta/inputparameters.hpp"
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
			virtual int getPrecision() const override
			{
				return fullParameters->get_precision();
			}
			virtual int getUseChemPotRe() const override
			{
				return fullParameters->get_use_chem_pot_re();
			}
			virtual int getUseChemPotIm() const override
			{
				return fullParameters->get_use_chem_pot_im();
			}
			virtual int getUseSmearing() const override
			{
				return fullParameters->get_use_smearing();
			}
			virtual int getChemPotRe() const override
			{
				return fullParameters->get_chem_pot_re();
			}
			virtual int getChemPotIm() const override
			{
				return fullParameters->get_chem_pot_im();
			}
			virtual int getRho() const override
			{
				return fullParameters->get_rho();
			}
			virtual int getRhoIter() const override
			{
				return fullParameters->get_rho_iter();
			}
			virtual int getUseRec12() const override
			{
				return fullParameters->get_use_rec12();
			}
			virtual int getUseEo() const override
			{
				return fullParameters->get_use_eo();
			}
			virtual int  getFermact() const override
			{
				return fullParameters->get_fermact();
			}
			virtual int getMetroApproxOrd() const override
			{
				return fullParameters->get_metro_approx_ord();
			}
			virtual int getMdApproxOrd() const override
			{
				return fullParameters->get_md_approx_ord();
			}
			virtual int getThetaFermionSpatial() const override
			{
				return fullParameters->get_theta_fermion_spatial();
			}
			virtual int getThetaFermionTemporal() const override
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
			virtual int getSourceContent() const override
			{
				return fullParameters->get_sourcecontent();
			}
			virtual int getUseAniso() const override
			{
				return fullParameters->get_use_aniso();
			}
		private:
			const meta::Inputparameters * fullParameters;
		};
	}
}
