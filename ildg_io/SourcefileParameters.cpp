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

#include "SourcefileParameters.hpp"

#include "../executables/exceptions.h"
#include "../host_functionality/logger.hpp"

Sourcefileparameters::Sourcefileparameters()
{
	set_defaults();
}

Sourcefileparameters::Sourcefileparameters(const meta::Inputparameters * parameters, int trajectoryNumber, double plaquette, Checksum checksumIn, std::string hmcVersion)
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

void Sourcefileparameters::printMetaDataToScreen(std::string sourceFilename)
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

void Sourcefileparameters::set_defaults()
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

std::string Sourcefileparameters::getInfo_xlfInfo()
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
	xlfInfo += "date = " + boost::lexical_cast<std::string>(date) + "\n";
	
	return xlfInfo;
}

std::string Sourcefileparameters::getInfo_ildgFormat_gaugefield()
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

std::string Sourcefileparameters::getInfo_scidacChecksum()
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
	
static void checkMajorParameter_int(int int1, int int2, std::string message)
{
	if(int1 != int2) {
		throw Invalid_Parameters("Major parameter \"" + message + "\" does not match! ", int1, int2);
	}
}

static void checkMinorParameter_double(double value1, double value2, std::string message)
{
	if(value1 != value2) {
		logger.warn() << "Minor parameter \"" + message + "\" does not match! ";
		logger.warn() << "\tExpected: " << value1 << "\tFound: " << value2;
	}
}
	
//todo: add check on plaquette
void Sourcefileparameters::checkAgainstInputparameters(const meta::Inputparameters * parameters)
{
	logger.info() << "Checking sourcefile parameters against inputparameters...";
	
	checkMajorParameter_int(this->lt, parameters->get_ntime(), "lt");
	checkMajorParameter_int(this->lx, parameters->get_nspace(), "lx");
	checkMajorParameter_int(this->ly, parameters->get_nspace(), "ly");
	checkMajorParameter_int(this->lz, parameters->get_nspace(), "lz");
	checkMajorParameter_int(this->prec, parameters->get_precision(), "precision");	
	
	checkMinorParameter_double(this->beta, parameters->get_beta(), "beta" );
	checkMinorParameter_double(this->kappa, parameters->get_kappa(), "kappa" );
	checkMinorParameter_double(this->mu, parameters->get_mu(), "mu" );
}

void Sourcefileparameters::checkAgainstChecksum(Checksum checksum, bool ignoreChecksumErrors, std::string filename)
{
	if(checksum != this->checksum) {
		logger.error() << "Checksum of data does not match checksum given in file.";
		logger.error() << "Calculated Checksum: " << checksum;
		logger.error() << "Embedded Checksum:   " << this->checksum;
		if(!ignoreChecksumErrors) 
		{
			throw File_Exception(filename);
		}
	}
}

size_t Sourcefileparameters::getSizeInBytes() noexcept
{
	return (size_t) num_entries * prec / 8;
}