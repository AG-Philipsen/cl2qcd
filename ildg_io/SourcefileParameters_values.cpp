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
	
	lx = parameters->get_nspace();
	ly = parameters->get_nspace();
	lz = parameters->get_nspace();
	lt = parameters->get_ntime();
	prec = parameters->get_precision();
	trajectorynr = trajectoryNumber;
	plaquettevalue = plaquette;
	beta = parameters->get_beta();
	kappa = parameters->get_kappa();
	mu = parameters->get_mu();
	
	hmcversion = hmcVersion;
	
	checksum = checksumIn;
}

void sourcefileparameters_values::printMetaDataToScreen(std::string sourceFilename)
{
  logger.info() << "*************************************************************" ;
  logger.info() << "*************************************************************" ;
  logger.info() << "Metadata from file " << sourceFilename << ":";
  logger.info() << "\treading XML-data gave:";
  logger.info() << "\t\tfield type:\t" << field ;
  logger.info() << "\t\tprecision:\t" << prec ;
  logger.info() << "\t\tlx:\t\t" << lx ;
  logger.info() << "\t\tly:\t\t" << ly ;
  logger.info() << "\t\tlz:\t\t" << lz ;
  logger.info() << "\t\tlt:\t\t" << lt ;
  logger.info() << "\t\tflavours:\t" << flavours ;
  logger.info() << "\treading XLF-data gave:";
  logger.info() << "\t\tplaquette:\t" << plaquettevalue;
  logger.info() << "\t\ttrajectorynr:\t" << trajectorynr;
  logger.info() << "\t\tbeta:\t\t" << beta;
  logger.info() << "\t\tkappa:\t\t" << kappa;
  logger.info() << "\t\tmu:\t\t" << mu;
  logger.info() << "\t\tc2_rec:\t\t" << c2_rec;
  logger.info() << "\t\ttime:\t\t" << time;
  logger.info() << "\t\thmc-version:\t" << hmcversion;
  logger.info() << "\t\tmubar:\t\t" << mubar;
  logger.info() << "\t\tepsilonbar:\t" << epsilonbar;
  logger.info() << "\t\tdate:\t\t" << date;
  if(numberOfFermionFieldsRead != 0) {
    logger.info() << "\treading inverter-data gave:";
    logger.info() << "\t\tsolvertype:\t" << solvertype;
    logger.info() << "\t\tepssq:\t\t" << std::setprecision(30) << epssq;
    logger.info() << "\t\tnoiter:\t\t" << noiter;
    logger.info() << "\t\tkappa_solver:\t" << kappa_solver;
    logger.info() << "\t\tmu_solver:\t" << mu_solver;
    logger.info() << "\t\ttime_solver:\t" << time_solver;
    logger.info() << "\t\thmc-ver_solver:\t" << hmcversion_solver;
    logger.info() << "\t\tdate_solver:\t" << date_solver;
  }
  logger.info() << "\tfile-checksum:\t" << checksum;
  logger.info() << "*************************************************************" ;
}

void sourcefileparameters_values::set_defaults()
{
	lx = 0;
	ly = 0;
	lz = 0;
	lt = 0;
	prec = 0;
	num_entries = 0;
	flavours = 0;
	trajectorynr = 0;
	time = 0;
	time_solver = 0;
	noiter = 0;
	plaquettevalue = 0;
	beta = 0;
	kappa = 0;
	mu = 0;
	c2_rec = 0;
	mubar = 0;
	epsilonbar = 0;
	epssq = 0;
	kappa_solver = 0;
	mu_solver = 0;
	
	numberOfFermionFieldsRead = 0;
	
	field = "";
	date = "";
	hmcversion = "";
	solvertype = "";
	hmcversion_solver = "";
	date_solver = "";
}

#include <time.h>

time_t getCurrentTime()
{
	time_t current_time;
	return time ( &current_time );
}

const char * getDateFromTime(time_t currentTime)
{
	return ctime (&currentTime);
}

std::string sourcefileparameters_values::getInfo_xlfInfo()
{
	time_t current_time = getCurrentTime();
	const char * date = getDateFromTime(current_time);
	
	std::string xlfInfo = "";
	xlfInfo += "plaquette = " + boost::lexical_cast<std::string>(this->plaquettevalue) + "\n";
	xlfInfo += "trajectory nr = " + boost::lexical_cast<std::string>(this->trajectorynr) + "\n";
	xlfInfo += "beta = " + boost::lexical_cast<std::string>(this->beta) + "\n";
	xlfInfo += "kappa = " + boost::lexical_cast<std::string>(this->kappa) + "\n";
	xlfInfo += "mu = " + boost::lexical_cast<std::string>(this->mu) + "\n";
	xlfInfo += "c2_rec = " + boost::lexical_cast<std::string>(this->c2_rec) + "\n";
	xlfInfo += "time = " + boost::lexical_cast<std::string>(current_time) + "\n";
	xlfInfo += "hmcversion = " + boost::lexical_cast<std::string>(this->hmcversion) + "\n";
	xlfInfo += "mubar = " + boost::lexical_cast<std::string>(this->mubar) + "\n";
	xlfInfo += "epsilonbar = " + boost::lexical_cast<std::string>(this->epsilonbar) + "\n";
	xlfInfo += "date = " + boost::lexical_cast<std::string>(this->date) + "\n";
	
	return xlfInfo;
}

std::string sourcefileparameters_values::getInfo_ildgFormat_gaugefield()
{
	std::string description = "su3gauge";
	
	std::string ildgFormat = "";
	ildgFormat += "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version>1.0</version>\n";
	ildgFormat += "  <field>" + description + "</field>\n";
	ildgFormat += "  <precision>" + boost::lexical_cast<std::string>(this->prec) + "</precision>\n";
	ildgFormat += "  <lx>" + boost::lexical_cast<std::string>(this->lx) + "</lx>\n";
	ildgFormat += "  <ly>" + boost::lexical_cast<std::string>(this->ly) + "</ly>\n";
	ildgFormat += "  <lz>" + boost::lexical_cast<std::string>(this->lz) + "</lz>\n";
	ildgFormat += "  <lt>" + boost::lexical_cast<std::string>(this->lt) + "</lt>\n";
	ildgFormat += "</ildgFormat>";
	
	return ildgFormat;
}

std::string sourcefileparameters_values::getInfo_scidacChecksum()
{
	std::string scidac_checksum;
	{
		std::ostringstream tmp;
		tmp << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<scidacChecksum>\n<version>1.0</version>\n";
		tmp << "<suma>" << std::hex << this->checksum.get_suma() << "</suma>\n";
		tmp << "<sumb>" << std::hex << this->checksum.get_sumb() << "</sumb>\n";
		tmp << "</scidacChecksum>";
		scidac_checksum = tmp.str();
	}
	return scidac_checksum;
}
	

