/** @file
 * physics::observables::wilson::TwoFlavourChiralCondensate class
 *
 * Copyright 2014 Christopher Pinke
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

#ifndef WILSONTWOFLAVOURCHIRALCONDENSATE_HPP_
#define WILSONTWOFLAVOURCHIRALCONDENSATE_HPP_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include <fstream>
#include <cmath>
#include "../../meta/inputparameters.hpp"

#include "../algorithms/inversion.hpp"
#include "../algorithms/flavour_doublet.hpp"
#include "../sources.hpp"

//into .cpp file
#include <cassert>

namespace physics{

	namespace observables{

		namespace wilson{

			class TwoFlavourChiralCondensate
			{
			public:
				TwoFlavourChiralCondensate(const meta::Inputparameters * parametersIn):
					chiralCondensate(0.)
				{
	  			parameters = parametersIn;
					checkInputparameters();
				}
				TwoFlavourChiralCondensate() = delete;
				
				double getChiralCondensate()
				{
					return chiralCondensate;
				}
				
				void flavour_doublet_chiral_condensate(const std::vector<physics::lattices::Spinorfield*>& inverted, const std::vector<physics::lattices::Spinorfield*>& sources, std::string pbp_fn, int number, const hardware::System& system)
				{
					if( parameters->get_pbp_version() == meta::Inputparameters::tm_one_end_trick ) 
					{
// 						flavour_doublet_chiral_condensate_tm(inverted, pbp_fn, number, system);
					} 
					else
					{
						flavour_doublet_chiral_condensate_std(inverted, sources, pbp_fn, number, system);
					}
				}
				
				void measureChiralCondensate(const physics::lattices::Gaugefield * gaugefield)
				{
					auto system = gaugefield->getSystem();
					auto prng = gaugefield->getPrng();
				
					std::string currentConfigurationName = "replace";
					filenameForChiralCondensateData = meta::get_ferm_obs_pbp_file_name(*parameters, currentConfigurationName);
					int sourceNumber = 0;

					for (; sourceNumber < parameters->get_num_sources(); sourceNumber++) {
						auto sources = physics::create_sources(*system, *prng, 1);
						auto result = physics::lattices::create_spinorfields(*system, sources.size());
						physics::algorithms::perform_inversion(&result, gaugefield, sources, *system);
						flavour_doublet_chiral_condensate(result, sources, filenameForChiralCondensateData, gaugefield->get_parameters_source().trajectorynr_source, *system);
						physics::lattices::release_spinorfields(result);
						physics::lattices::release_spinorfields(sources);
					}
				}
				

			private:
				const meta::Inputparameters * parameters;
				double chiralCondensate;
				std::ofstream outputToFile;
				std::string filenameForChiralCondensateData;
				
				void writeChiralCondensateToFile(int iter,  const std::string& filename);
				
				void checkInputparameters()
				{
					if(! parameters->get_measure_pbp() ) 
					{
						throw std::logic_error("Chiral condensate calculation disabled in parameter setting. Aborting...");
					}
					
					if( ! (  ( parameters->get_pbp_version() == meta::Inputparameters::std) ||
						       ((parameters->get_fermact() == meta::Inputparameters::twistedmass && parameters->get_pbp_version() == meta::Inputparameters::tm_one_end_trick ) )
								) )
						throw std::logic_error("No valid chiral condensate version has ben selected.");
				}
				
				void flavour_doublet_chiral_condensate_std(const std::vector<physics::lattices::Spinorfield*>& solved_fields, const std::vector<physics::lattices::Spinorfield*>& sources, std::string pbp_fn, int number, const hardware::System& system)
				{
					using namespace physics::lattices;

					auto params = system.get_inputparameters();

					// Output
					using namespace std;
					ofstream of;
					of.open(pbp_fn.c_str(), std::ios_base::app);
					if(!of.is_open()) {
						throw File_Exception(pbp_fn);
					}

					hmc_float result = 0.;
					/**
					* In the pure Wilson case one can evaluate <pbp> with stochastic estimators according to:
					* <pbp> = <ubu + dbd> = 2<ubu>  
					*       = 2 Tr_(space, colour, dirac) ( D^-1 )
					*       = lim_r->inf 2/r (Xi_r, Phi_r)
					* where the estimators satisfy
					* D^-1(x,y)_(a,b, A,B) = lim_r->inf Phi_r(x)_a,A (Xi_r(y)_b,B)^dagger
					* and Phi fulfills
					* D Phi = Xi
					* (X,Y) denotes the normal scalar product
					* In the twisted-mass case one can evaluate <pbp> with stochastic estimators similar to the pure wilson case.
					* However, one first has to switch to the twisted basis:
					* <pbp> -> <chibar i gamma_5 tau_3 chi>
					*       = <ub i gamma_5 u> - <db i gamma_5 d>
					*       = Tr( i gamma_5 (D^-1_u - D^-1_d ) )
					*       = Tr( i gamma_5 (D^-1_u - gamma_5 D^-1_u^dagger gamma_5) )
					*       = Tr( i gamma_5 (D^-1_u -  D^-1_u^dagger ) )
					*       = - Nf Im Tr ( gamma_5 D^-1_u)
					*       = - lim_r->inf Nf/r  (gamma_5 Xi_r, Phi_r)
					* NOTE: The basic difference compared to the pure Wilson case is the gamma_5 and that one takes the negative imaginary part!
					*/
					// Need 2 spinors at once..
					logger.debug() << "init buffers for chiral condensate calculation...";
					Spinorfield phi(system);
					Spinorfield xi(system);
					assert(solved_fields.size() == sources.size());
					/*
					* Normalize for VOL4D, NF and spinor entries (NC * ND = 12)
					* In addition, a factor of 2 kappa should be inserted to convert to the physical basis.
					* The additional factor of 2 is inserted to fit the reference values.
					*/
					hmc_float norm = 4. * params.get_kappa() * 2. / meta::get_vol4d(params) / 2. / 12.;
					/**
					* Currently, there is a problem with the sign if even-odd is used (issue #389).
					* This seems to cause the a wrong sign in the chiral condensate, which will be compensated for now.
					*/
					if(params.get_use_eo() ){
						norm *= -1.;
					}
					logger.info() << "chiral condensate:" ;
					for(size_t i = 0; i < solved_fields.size(); ++i) {
						copyData(&phi, solved_fields[i]);
						copyData(&xi, sources[i]);

						if(params.get_fermact() == meta::Inputparameters::twistedmass) {
							xi.gamma5();
						}
						
						hmc_complex tmp = scalar_product(xi, phi);
						tmp.re *= norm;
						tmp.im *= norm;
						switch(params.get_fermact()) {
							case  meta::Inputparameters::wilson:
								result = tmp.re;
								break;
							case meta::Inputparameters::twistedmass:
								result = (-1.) * tmp.im;
								break;
							default:
								throw std::invalid_argument("chiral condensate not implemented for given fermion action");
						}
						logger.info() << number << "\t" << scientific << setprecision(14) << result;
						of << number << "\t" << scientific << setprecision(14) << result << endl;
					}
					chiralCondensate = result;
				}

      };

			double measureChiralCondensate(const physics::lattices::Gaugefield * gaugefield, int iteration)
			{
				TwoFlavourChiralCondensate condensate(gaugefield->getParameters() );
				condensate.measureChiralCondensate(gaugefield);
				return condensate.getChiralCondensate();
			}
    }
  }
}
#endif
