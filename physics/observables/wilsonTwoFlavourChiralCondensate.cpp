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

#include "wilsonTwoFlavourChiralCondensate.hpp"

#include "../algorithms/inversion.hpp"
#include "../sources.hpp"

class TwoFlavourChiralCondensate
{
public:
  TwoFlavourChiralCondensate(const physics::lattices::Gaugefield * gaugefield,
                             const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& parametersInterface,
                             std::string configurationName = "conf.default", int iteration = 0);
  TwoFlavourChiralCondensate() = delete;
  ~TwoFlavourChiralCondensate();
  
  std::vector<double> getChiralCondensate();
  void measureChiralCondensate(const physics::lattices::Gaugefield * gaugefield, physics::InterfacesHandler & interfacesHandler);
  void writeChiralCondensateToFile();
  
private:
  const physics::lattices::Gaugefield * gaugefield;
  const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& parametersInterface;
  const hardware::System * system;
  const physics::PRNG * prng;
  int trajectoryNumber;
  std::vector<double> chiralCondensate;
  std::ofstream outputToFile;
  std::string filenameForChiralCondensateData;
  std::string configurationName;
  
  void checkInputparameters();
  double norm_std() const ;
  double norm_tm() const;
  double flavourChiralCondensate_std(const physics::lattices::Spinorfield* phi, const physics::lattices::Spinorfield* xi);
  double flavour_doublet_chiral_condensate_tm(const physics::lattices::Spinorfield* phi);
  void openFileForWriting();
  void flavour_doublet_chiral_condensate(const physics::lattices::Spinorfield* inverted, const physics::lattices::Spinorfield* sources);
};

TwoFlavourChiralCondensate::TwoFlavourChiralCondensate(const physics::lattices::Gaugefield * gaugefieldIn,
                                                       const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& parametersInterface,
                                                       std::string configurationNameIn, int trajectoryNumberIn)
      : gaugefield(gaugefieldIn), parametersInterface(parametersInterface), system(gaugefield->getSystem() ),
        prng(gaugefield->getPrng()), chiralCondensate(), configurationName(configurationNameIn)
{
	checkInputparameters();
	openFileForWriting();
	trajectoryNumber = trajectoryNumberIn;
}

void checkFermionAction(const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& parametersInterface)
{
	if( !( (parametersInterface.getFermionicActionType() == common::action::twistedmass) || ( parametersInterface.getFermionicActionType() == common::action::wilson) )	)
		throw std::logic_error("Chiral condensate not implemented for chosen fermion action.");
}

void checkChiralCondensateVersion(const physics::observables::WilsonTwoFlavourChiralCondensateParametersInterface& parametersInterface)
{
	if( !( ( parametersInterface.getPbpVersion() == common::pbp_version::std) ||
		 ( (parametersInterface.getFermionicActionType() == common::action::twistedmass && parametersInterface.getPbpVersion() == common::pbp_version::tm_one_end_trick ) ) ) )
		throw std::logic_error("No valid chiral condensate version has been selected.");
}

void TwoFlavourChiralCondensate::checkInputparameters()
{
	if(! parametersInterface.measurePbp() )
	{
		throw std::logic_error("Chiral condensate calculation disabled in parameter setting. Aborting...");
	}
	
	checkFermionAction(parametersInterface);
	checkChiralCondensateVersion(parametersInterface);
}

TwoFlavourChiralCondensate::~TwoFlavourChiralCondensate()
{
	outputToFile.close();
}

std::vector<double> TwoFlavourChiralCondensate::getChiralCondensate()
{
	return chiralCondensate;
}
				
void TwoFlavourChiralCondensate::measureChiralCondensate(const physics::lattices::Gaugefield * gaugefield, physics::InterfacesHandler & interfacesHandler)
{
	logger.info() << "chiral condensate:" ;
	for (int sourceNumber = 0; sourceNumber < parametersInterface.getNumberOfSources(); sourceNumber++) {
		auto sources = physics::create_sources(*system, *prng, 1, interfacesHandler);
		auto result = physics::lattices::create_spinorfields(*system, sources.size(), interfacesHandler);
		physics::algorithms::perform_inversion(&result, gaugefield, sources, *system, interfacesHandler);
		flavour_doublet_chiral_condensate(result[0], sources[0]);
		physics::lattices::release_spinorfields(result);
		physics::lattices::release_spinorfields(sources);
	}
}

void printChiralCondensate(int trajectoryNumber, double value)
{
	logger.info() << trajectoryNumber << "\t" << std::scientific << std::setprecision(14) << value;
}

void TwoFlavourChiralCondensate::writeChiralCondensateToFile()
{
    logger.info () << "Write chiral condensate data to file \"" << filenameForChiralCondensateData << "\" ...";
    outputToFile << trajectoryNumber << "\t";
    outputToFile.precision(15);
    outputToFile.setf( std::ios::scientific, std::ios::floatfield );
    for (int i = 0; i < (int) chiralCondensate.size(); i++)
    {
        outputToFile << chiralCondensate[i] << "   ";
    }
    outputToFile << std::endl;
}

double TwoFlavourChiralCondensate::norm_std() const 
{
	/**
	 * Normalize for VOL4D, NF and spinor entries (NC * ND = 12)
	 * In addition, a factor of 2 kappa should be inserted to convert to the physical basis.
	 * The additional factor of 2 is inserted to fit the reference values.
	 * The additional factor of 2 is Nf.
	 */
	double norm =  2. * 2. * parametersInterface.getKappa() * 2. / parametersInterface.get4dVolume() / 2. / 12.;
	/**
	 * Currently, there is a problem with the sign if even-odd is used (issue #389).
	 * This seems to cause the a wrong sign in the chiral condensate, which will be compensated for now.
	 */
	if(parametersInterface.useEvenOdd()){
		norm *= -1.;
	}
	return norm;
}

double TwoFlavourChiralCondensate::flavourChiralCondensate_std(const physics::lattices::Spinorfield* phi, const physics::lattices::Spinorfield* xi)
{
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
	double result;

	if(parametersInterface.getFermionicActionType() == common::action::twistedmass) {
		xi->gamma5();
	}
	hmc_complex tmp = scalar_product(*xi, *phi);

	switch(parametersInterface.getFermionicActionType()) {
		case  common::action::wilson:
			result = tmp.re * norm_std();
			break;
		case common::action::twistedmass:
			result = (-1.) * tmp.im * norm_std();
			break;
		default:
			throw std::invalid_argument("chiral condensate not implemented for given fermion action");
	}

	return result;
}

void TwoFlavourChiralCondensate::openFileForWriting()
{
	filenameForChiralCondensateData = parametersInterface.getPbpFilename(configurationName);
	outputToFile.open(filenameForChiralCondensateData.c_str(), std::ios_base::app);
	if(!outputToFile.is_open()) {
		throw File_Exception(filenameForChiralCondensateData);
	}
}

void TwoFlavourChiralCondensate::flavour_doublet_chiral_condensate(const physics::lattices::Spinorfield* inverted, const physics::lattices::Spinorfield* sources)
{
    double result = 0.;
    if( parametersInterface.getPbpVersion() == common::pbp_version::tm_one_end_trick )
	{
	  result = flavour_doublet_chiral_condensate_tm(inverted);
	} 
	else
	{
	  result = flavourChiralCondensate_std(inverted, sources);
	}
 	printChiralCondensate(trajectoryNumber, result );
	chiralCondensate.push_back(result);
}

std::vector<double> physics::observables::wilson::measureTwoFlavourChiralCondensateAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, std::string currentConfigurationName, physics::InterfacesHandler & interfacesHandler)
{
    TwoFlavourChiralCondensate condensate(gaugefield, interfacesHandler.getWilsonTwoFlavourChiralCondensateParametersInterface(),
                                          currentConfigurationName, gaugefield->get_trajectoryNumberAtInit());
	condensate.measureChiralCondensate(gaugefield, interfacesHandler);
	condensate.writeChiralCondensateToFile();
	return condensate.getChiralCondensate();
}

std::vector<double> physics::observables::wilson::measureTwoFlavourChiralCondensateAndWriteToFile(const physics::lattices::Gaugefield * gaugefield, int iteration, physics::InterfacesHandler & interfacesHandler)
{
    TwoFlavourChiralCondensate condensate(gaugefield, interfacesHandler.getWilsonTwoFlavourChiralCondensateParametersInterface(),
	                                      gaugefield->getName(iteration), iteration);
	condensate.measureChiralCondensate(gaugefield, interfacesHandler);
	condensate.writeChiralCondensateToFile();
	return condensate.getChiralCondensate();
}

double TwoFlavourChiralCondensate::norm_tm() const 
{
	/*
	 * Normalize for VOL4D, Nf and spinor entries (NC * ND = 12)
	 * In addition, a factor of 2 kappa should be inserted to convert to the physical basis.
	 * The additional factor of 2 is inserted to fit the reference values.
	 */
	double norm = 4. * parametersInterface.getKappa()  / parametersInterface.get4dVolume()  * parametersInterface.getMubar() * 2. / 2. / 12.;
	return norm;
}

double TwoFlavourChiralCondensate::flavour_doublet_chiral_condensate_tm(const physics::lattices::Spinorfield* phi)
{
	/**
	 * For twisted-mass fermions one can also employ the one-end trick, which origins from
	 * D_d - D_u = - 4 i kappa amu gamma_5 <-> D^-1_u - D^-1_d = - 4 i kappa amu gamma_5 (D^-1_u)^dagger D^-1_u
	 * With this, the chiral condensate is:
	 * <pbp> = ... 
	 *       = Tr( i gamma_5 (D^-1_u - D^-1_d ) )
	 *       = - 4 i kappa amu Tr ( i gamma_5 gamma_5 (D^-1_u)^dagger D^-1_u )
	 *       = 4 kappa amu Tr ((D^-1_u)^dagger D^-1_u) 
	 *       = 4 kappa amu lim_r->inf 1/R (Phi_r, Phi_r)
	 *       = 4 kappa amu lim_r->inf 1/R |Phi_r|^2
	 * NOTE: Here one only needs Phi...
	 */
	double result = 0.;
	logger.info() << "chiral condensate:" ;
	hmc_float tmp = squarenorm(*phi);
	result = tmp * norm_tm();

	return result;
}
