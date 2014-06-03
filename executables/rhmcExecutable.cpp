/*
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

#include "rhmcExecutable.h"

rhmcExecutable::rhmcExecutable(int argc, const char* argv[]) :  generationExecutable(argc, argv, "rhmc")
{
	using namespace physics::algorithms;
	
	checkRhmcParameters(parameters);
	initializationTimer.reset();
	printParametersToScreenAndFile();
	setIterationParameters();
	initializationTimer.add();
	logger.info() << "Generation of Rational Approximations...";
	if(parameters.get_read_rational_approximations_from_file()){
		approx_hb  = new Rational_Approximation(parameters.get_approx_heatbath_file());
		approx_md  = new Rational_Approximation(parameters.get_approx_md_file());
		approx_met = new Rational_Approximation(parameters.get_approx_metropolis_file());
	}else{
		//This is the approx. to be used to generate the initial (pseudo)fermionic field
		approx_hb = new Rational_Approximation(parameters.get_metro_approx_ord(),
				parameters.get_num_tastes(), 8, parameters.get_approx_lower(),
				parameters.get_approx_upper(), false);
		//This is the approx. to be used to generate the initial (pseudo)fermionic field
		approx_md = new Rational_Approximation(parameters.get_md_approx_ord(),
				parameters.get_num_tastes(), 4, parameters.get_approx_lower(),
				parameters.get_approx_upper(), true);
		//This is the approx. to be used to generate the initial (pseudo)fermionic field
		approx_met = new Rational_Approximation(parameters.get_metro_approx_ord(),
				  parameters.get_num_tastes(), 4, parameters.get_approx_lower(),
				  parameters.get_approx_upper(), true);
	}
}

rhmcExecutable::~rhmcExecutable()
{
  using namespace std;
  logger.info() << "Acceptance rate: " << fixed <<  setprecision(1) << percent(acceptanceRate, parameters.get_hmcsteps()) << "%";
}

void rhmcExecutable::printParametersToScreenAndFile()
{
  meta::print_info_rhmc(parameters);
  writeRhmcLogfile();
}

void rhmcExecutable::writeRhmcLogfile()
{
	outputToFile.open(filenameForLogfile,
			std::ios::out | std::ios::app);
	if (outputToFile.is_open()) {
		meta::print_info_rhmc(&outputToFile, parameters);
		outputToFile.close();
	} else {
		throw File_Exception(filenameForLogfile);
	}
}

void rhmcExecutable::setIterationParameters()
{
	generationExecutable::setIterationParameters();
	generationSteps += parameters.get_rhmcsteps();
}

void rhmcExecutable::thermalizeAccordingToSpecificAlgorithm()
{
	logger.warn() << "Thermalization is not yet implemented for RHMC algorithm: --> NOP!";
}

void rhmcExecutable::generateAccordingToSpecificAlgorithm()
{
	const double randomNumber = prng->get_double();
	observables = physics::algorithms::perform_rhmc_step(*approx_hb, *approx_md, *approx_met, gaugefield, iteration, randomNumber, *prng, *system);
	acceptanceRate += observables.accept;
}

void rhmcExecutable::performOnlineMeasurements()
{
	if( ( (iteration + 1) % writeFrequency ) == 0 ) {
		std::string gaugeout_name = meta::get_rhmc_obs_file_name(parameters, "");
		printRhmcObservables(gaugeout_name);
	}
}

void rhmcExecutable::printRhmcObservables(const std::string& filename)
{
  printRhmcObservablesToFile(filename);
  printRhmcObservablesToScreen();
}

void rhmcExecutable::printRhmcObservablesToFile(const std::string& filename)
{
	const hmc_float exp_deltaH = std::exp(observables.deltaH);
	outputToFile.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!outputToFile.is_open()) throw File_Exception(filename);
	outputToFile << iteration << "\t";
	outputToFile.width(8);
	outputToFile.precision(15);
	outputToFile << observables.plaq << "\t" << observables.tplaq << "\t" << observables.splaq;
	outputToFile << "\t" << observables.poly.re << "\t" << observables.poly.im << "\t" << sqrt(observables.poly.re * observables.poly.re + observables.poly.im * observables.poly.im);
	if(parameters.get_measure_pbp()) {
	     std::vector<hmc_complex> pbp(4);
	     for(size_t i=0; i<pbp.size(); i++){
	         pbp[i] = physics::observables::staggered::measureChiralCondensate(*gaugefield, *prng, *system);
		 outputToFile <<  "\t" << pbp[i].re << "\t" << pbp[i].im;
	     }
	}
	outputToFile <<  "\t" << observables.deltaH << "\t" << exp_deltaH << "\t" << observables.prob << "\t" << observables.accept;
	//print number of iterations used in inversions with full and force precision
	/**
	 * @todo: The counters should be implemented once the solver class is used!"
	 * until then, only write "0"!
	 */
	int iter0 = 0;
	int iter1 = 0;
	outputToFile << "\t" << iter0 << "\t" << iter1;
	if(parameters.get_use_mp() ) {
		outputToFile << "\t" << iter0 << "\t" << iter1;
	}
	if(meta::get_use_rectangles(parameters) ) {
		outputToFile << "\t" << observables.rectangles;
	}
	outputToFile << std::endl;
	outputToFile.close();
}

void rhmcExecutable::printRhmcObservablesToScreen()
{
	logger.info() << "\tRHMC [OBS]:\t" << iteration << std::setw(8) << std::setfill(' ') << "\t" << std::setprecision(15) 
		      << observables.plaq << "\t" << observables.poly.re << "\t" << observables.poly.im;
}

void rhmcExecutable::checkRhmcParameters(const meta::Inputparameters& p)
{
	if(p.get_fermact() != meta::action::rooted_stagg)
	  throw Invalid_Parameters("Fermion action not suitable for RHMC!", meta::action::rooted_stagg, p.get_fermact());
	if(!p.get_use_eo())
	  throw Invalid_Parameters("RHMC available only WITH eo-prec!", "use_eo=1", p.get_use_eo());
	if(p.get_use_mp())
	  throw Invalid_Parameters("RHMC available only WITHOUT mass preconditionig!", "use_mp=0", p.get_use_mp());
	if(p.get_use_chem_pot_re())
	  throw Invalid_Parameters("RHMC available only WITHOUT real chemical potential!", "use_chem_pot_re=0", p.get_use_chem_pot_re());
	if(p.get_use_chem_pot_im())
	  throw Invalid_Parameters("RHMC available only WITHOUT imaginary chemical potential!", "use_chem_pot_im=0", p.get_use_chem_pot_im());
	if(p.get_num_tastes()%4 == 0)
	  throw Invalid_Parameters("RHMC not working with multiple of 4 tastes (there is no need of the rooting trick)!", "num_tastes%4 !=0", "num_tastes=" + std::to_string(p.get_num_tastes()));
	if(p.get_cg_iteration_block_size() == 0 || p.get_findminmax_iteration_block_size() == 0)
	  throw Invalid_Parameters("Iteration block sizes CANNOT be zero!", "cg_iteration_block_size!=0 && findminmax_iteration_block_size!=0", p.get_cg_iteration_block_size()==0 ? "cg_iteration_block_size=0" : "findminmax_iteration_block_size=0");
	if(p.get_approx_upper() != 1)
	  throw Invalid_Parameters("RHMC not available if the Rational expansion is not calculated in [..,1]!", "approx_upper=1", p.get_approx_upper());
	if(p.get_approx_lower() >= 1)
	   throw Invalid_Parameters("The lower bound of the Rational expansion >= than the upper one does not make sense!", "approx_lower < approx_upper", std::to_string(p.get_approx_lower())+" >= 1");
	if(p.get_writefrequency() == 0)
	  throw Invalid_Parameters("Write frequency CANNOT be zero!", "writefrequency!=0", p.get_writefrequency());
	if(p.get_savefrequency() == 0){
	  logger.warn() << "savefrequency==0 and hence NEITHER gaugefield NOR prng will be saved!!";
	  logger.warn() << "The simulation will start in 5 seconds..";
	  sleep(2); logger.warn() << "..3.."; sleep(1); logger.warn() << "..2.."; sleep(1);
	  logger.warn() << "..1, go!"; sleep(1);
	}
}
