/** @file
 * Implementation of the staggered chiral condensate calculation.
 *
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

#include "staggeredChiralCondensate.hpp"

#include "../lattices/staggeredfield_eo.hpp"
#include "../lattices/scalar_complex.hpp"
#include "../fermionmatrix/fermionmatrix_stagg.hpp"
#include "../sources.hpp"
#include "../algorithms/solver_shifted.hpp"

hmc_complex physics::observables::staggered::measureChiralCondensate(const physics::lattices::Gaugefield& gf, const physics::PRNG& prng, const hardware::System& system,
                                                                     physics::InterfacesHandler& interfacesHandler)
{
	/**
	 * The chiral condensate in the RHMC algorithm turns out to be
	 * @code
	 * <pbp> = 1/VOL4D*N_flavour/4 * <Tr(M^{-1})>
	 * @endcode
	 * where the <...> is the average on the RHMC ensembles. Then, since to
	 * calculate directly the inverse staggered Dirac matrix M^{-1} is too costly,
	 * noise estimators are used. This means that
	 * @code
	 * Tr(M^{-1}) = < eta^dag * M^{-1} * eta > 
	 * @endcode
	 * where eta is the noise estimator whose each number (each component per each site)
	 * is drawn independently from the others. The last equation above is only true if
	 * the draw procedure is really independent and in general it is approximated via
	 * @code
	 * Tr(M^{-1}) \approx 1/L \sum_{i=1}^L eta^dag_i * M^{-1} * eta_i
	 * @endcode
	 * where L is the number of noise sources we use. A bigger number of sources will
	 * give a better estimate, but at the same time will be more time-consuming. Let us
	 * rewrite the rhs of the last expression above in a more convenient way to be implemented.
	 * @code
	 * eta^dag_i * M^{-1} * eta_i = eta^dag_i * (M^dag*M)^{-1} M^dag * eta_i
	 * @endcode
	 * Using the even-odd notation, we have @code eta = (eta_e, eta_o) @endcode and it can
	 * be shown that
	 * @code
	 * (M^dag * eta) = (m*eta_e - Deo*eta_o, m*eta_o - Doe*eta_e)
	 * 
	 * (M^dag*M)^{-1} M^dag * eta_i = ( [(M^dag*M)^{-1}]ee * (m*eta_e - Deo*eta_o),
	 *                                  [(M^dag*M)^{-1}]oo * (m*eta_o - Doe*eta_e) ) \equiv (chi_e, chi_o)
	 * @endcode
	 * At this point one could easily implement the above expression, but this require two
	 * inversion per source and we know that the inversion is the most time-consuming part.
	 * One of these two can be avoided. In fact, chi_o can be obtained from chi_e as
	 * @code
	 * chi_o = 1/m * (eta_o - Doe * chi_e)
	 * @endcode
	 * where m is the mass of the particles. In the end we will have
	 * @code
	 * Tr(M^{-1}) \approx 1/L \sum_{i=1}^L (eta_e^dag_i * chi_e + eta_o^dag_i * chi_o)
	 * @endcode
	 * where L in the following is obtained from the parameterInterface.
	 */
	
	using namespace physics::lattices;
	using namespace physics::algorithms::solvers;
	
	const physics::observables::StaggeredChiralCondensateParametersInterface& parametersInterface = interfacesHandler.getStaggeredChiralCondensateParametersInterface();
	const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<physics::lattices::Staggeredfield_eo>();
	const int number_sources = parametersInterface.getNumberOfSources();
	const hmc_float mass = additionalParameters.getMass();
	//Result
	hmc_complex pbp = {0.0, 0.0}; 
	
	for(int i=0; i<number_sources; i++){
	  //Noise sources
	  Staggeredfield_eo eta_e(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	  Staggeredfield_eo eta_o(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	  //Auxiliary fields
	  std::vector<std::shared_ptr<Staggeredfield_eo> > chi_e; //This is the type to be used in the inverter
	  chi_e.emplace_back(std::make_shared<Staggeredfield_eo>(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>()));
	  Staggeredfield_eo chi_o(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
	  //Fermionmatrix objects
	  physics::fermionmatrix::D_KS_eo Deo(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), EVEN);
	  physics::fermionmatrix::D_KS_eo Doe(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), ODD);
	  physics::fermionmatrix::MdagM_eo MdagM(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
	  
	  /***********************************************************************************************/
	  set_volume_source(&eta_e, prng);
	  set_volume_source(&eta_o, prng);
	  
	  //Calculate chi_e = [(M^dag*M)^{-1}]ee * (m*eta_e - Deo*eta_o)   using chi_o as temporary field
	  Deo(&chi_o, gf, eta_o);
	  saxpby(&chi_o, mass, eta_e, -1.0, chi_o);
	  //Here the CGM as standard CG is used
	  std::vector<hmc_float> sigma(1, 0.0); //only one shift set to 0.0
	  cg_m(chi_e, MdagM, gf, sigma, chi_o, system, interfacesHandler, parametersInterface.getSolverPrecision(), additionalParameters);
	  
	  //Calculate chi_o = 1/m * (eta_o - Doe * chi_e)
	  Doe(&chi_o, gf, *(chi_e[0]));
	  saxpby(&chi_o, 1.0/mass, eta_o, -1.0/mass, chi_o);
	  
	  //Build up the trace of M^{-1}, namely pbp += (eta_e^dag_i * chi_e + eta_o^dag_i * chi_o)
	  Scalar<hmc_complex> tmp1(system);
	  Scalar<hmc_complex> tmp2(system);
	  scalar_product(&tmp1, eta_e, *(chi_e[0]));
	  scalar_product(&tmp2, eta_o, chi_o);
	  add(&tmp1, tmp1, tmp2);
	  pbp += tmp1.get();
	}
	
	//Multiply by the overall factor, namely pbp = 1/VOL4D*N_flavour/4 * <Tr(M^{-1})>
	//Here we also divide by number_sources that has not yet been done.
	hmc_float factor = 1.0/parametersInterface.get4dVolume() * parametersInterface.getNumberOfTastes() * 0.25 / number_sources;
	pbp.re *= factor;
	pbp.im *= factor;
	
	return pbp;
}


void physics::observables::staggered::measureChiralCondensateAndWriteToFile(const physics::lattices::Gaugefield& gf, int iteration,
                                                                            physics::InterfacesHandler& interfacesHandler)
{
    const physics::observables::StaggeredChiralCondensateParametersInterface& parametersInterface = interfacesHandler.getStaggeredChiralCondensateParametersInterface();
    if( !parametersInterface.measurePbp() )
    {
        throw std::logic_error("Chiral condensate calculation disabled in parameter setting. Aborting...");
    }
    std::ofstream outputToFile;
    std::string configurationName = gf.getName(iteration);
    //TODO: improve this!
    std::string filenameForChiralCondensateData = parametersInterface.getPbpFilename(configurationName);
    logger.info () << "  Write chiral condensate data to file \"" << filenameForChiralCondensateData << "\" ...";
    outputToFile.open(filenameForChiralCondensateData.c_str(), std::ios_base::app);
    if(!outputToFile.is_open()) {
        throw File_Exception(filenameForChiralCondensateData);
    }
    outputToFile << iteration << "\t";
    outputToFile.precision(15);
    outputToFile.setf( std::ios::scientific, std::ios::floatfield );
    std::vector<hmc_complex> pbp(parametersInterface.getPbpNumberOfMeasurements());
    for(size_t i=0; i<pbp.size(); i++){
        pbp[i] = physics::observables::staggered::measureChiralCondensate(gf, *(gf.getPrng()), *(gf.getSystem()), interfacesHandler);
        outputToFile << pbp[i].re << "   ";
    }
    outputToFile << std::endl;
    outputToFile.close();  
    logger.info () << "  ...done!";
}

