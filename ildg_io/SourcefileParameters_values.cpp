/** @file
 *
 * Copyright 2014, Christopher Pinke
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

#include "SourcefileParameters_values.hpp"

sourcefileparameters_values::sourcefileparameters_values()
{
	set_defaults();
}

sourcefileparameters_values::sourcefileparameters_values(const meta::Inputparameters * parameters, int trajectoryNumber, double plaquette, Checksum checksumIn, std::string hmcVersion)
{
	set_defaults();
	
	lx_source = parameters->get_nspace();
	ly_source = parameters->get_nspace();
	lz_source = parameters->get_nspace();
	lt_source = parameters->get_ntime();
	prec_source = parameters->get_precision();
	trajectorynr_source = trajectoryNumber;
	plaquettevalue_source = plaquette;
	beta_source = parameters->get_beta();
	kappa_source = parameters->get_kappa();
	mu_source = parameters->get_mu();
	
	hmcversion_source = hmcVersion;
	
	checksum = checksumIn;
}

void sourcefileparameters_values::printMetaDataToScreen(std::string sourceFilename)
{
  logger.info() << "*************************************************************" ;
  logger.info() << "*************************************************************" ;
  logger.info() << "Metadata from file " << sourceFilename << ":";
  logger.info() << "\treading XML-data gave:";
  logger.info() << "\t\tfield type:\t" << field_source ;
  logger.info() << "\t\tprecision:\t" << prec_source ;
  logger.info() << "\t\tlx:\t\t" << lx_source ;
  logger.info() << "\t\tly:\t\t" << ly_source ;
  logger.info() << "\t\tlz:\t\t" << lz_source ;
  logger.info() << "\t\tlt:\t\t" << lt_source ;
  logger.info() << "\t\tflavours:\t" << flavours_source ;
  logger.info() << "\treading XLF-data gave:";
  logger.info() << "\t\tplaquette:\t" << plaquettevalue_source;
  logger.info() << "\t\ttrajectorynr:\t" << trajectorynr_source;
  logger.info() << "\t\tbeta:\t\t" << beta_source;
  logger.info() << "\t\tkappa:\t\t" << kappa_source;
  logger.info() << "\t\tmu:\t\t" << mu_source;
  logger.info() << "\t\tc2_rec:\t\t" << c2_rec_source;
  logger.info() << "\t\ttime:\t\t" << time_source;
  logger.info() << "\t\thmc-version:\t" << hmcversion_source;
  logger.info() << "\t\tmubar:\t\t" << mubar_source;
  logger.info() << "\t\tepsilonbar:\t" << epsilonbar_source;
  logger.info() << "\t\tdate:\t\t" << date_source;
  if(numberOfFermionFieldsRead != 0) {
    logger.info() << "\treading inverter-data gave:";
    logger.info() << "\t\tsolvertype:\t" << solvertype_source;
    logger.info() << "\t\tepssq:\t\t" << std::setprecision(30) << epssq_source;
    logger.info() << "\t\tnoiter:\t\t" << noiter_source;
    logger.info() << "\t\tkappa_solver:\t" << kappa_solver_source;
    logger.info() << "\t\tmu_solver:\t" << mu_solver_source;
    logger.info() << "\t\ttime_solver:\t" << time_solver_source;
    logger.info() << "\t\thmc-ver_solver:\t" << hmcversion_solver_source;
    logger.info() << "\t\tdate_solver:\t" << date_solver_source;
  }
  logger.info() << "\tfile-checksum:\t" << checksum;
  logger.info() << "*************************************************************" ;
}

void sourcefileparameters_values::set_defaults()
{
	lx_source = 0;
	ly_source = 0;
	lz_source = 0;
	lt_source = 0;
	prec_source = 0;
	num_entries_source = 0;
	flavours_source = 0;
	trajectorynr_source = 0;
	time_source = 0;
	time_solver_source = 0;
	noiter_source = 0;
	plaquettevalue_source = 0;
	beta_source = 0;
	kappa_source = 0;
	mu_source = 0;
	c2_rec_source = 0;
	mubar_source = 0;
	epsilonbar_source = 0;
	epssq_source = 0;
	kappa_solver_source = 0;
	mu_solver_source = 0;
	
	numberOfFermionFieldsRead = 0;
	
	field_source = "";
	date_source = "";
	hmcversion_source = "";
	solvertype_source = "";
	hmcversion_solver_source = "";
	date_solver_source = "";
}