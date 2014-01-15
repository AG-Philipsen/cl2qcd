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

#include "ildg_write_gaugefield.h"

#include <assert.h>

#include "../logger.hpp"
#include <sstream>

void write_gaugefield (
  char * binary_data, n_uint64_t num_bytes, Checksum checksum,
  int lx, int ly, int lz, int lt, int prec, int trajectorynr, hmc_float plaquettevalue, hmc_float beta, hmc_float kappa, hmc_float mu, hmc_float c2_rec, hmc_float epsilonbar, hmc_float mubar,
  const char * hmc_version, const char * filename)
{

	logger.info() << "writing gaugefield to lime-file...";

	time_t current_time;
	FILE *outputfile;
	outputfile = fopen(filename, "w");
	int MB_flag;
	int ME_flag;
	n_uint64_t length_xlf_info = 0, length_ildg_format = 0, length_scidac_checksum = 0;
	LimeRecordHeader * header_ildg_format, *header_scidac_checksum, * header_ildg_binary_data, * header_xlf_info;

	//set values
	const char * field_out = "su3gauge";
	time_t rawtime;
	//here is the problem with current_time: rawtime has to be converted into a meaningful format
	current_time = rawtime;
	time ( &rawtime );
	const char * date = ctime (&rawtime);

	// TODO replace this whole block by something templated
	//get binary data
	//here it must not be assumed that the argument prec and sizeof(hmc_float) are the same!!
	logger.debug() << "  num_bytes = " << num_bytes;

	if(sizeof(hmc_float) * 8 != prec) {
		throw Invalid_Parameters("Precision does not match executables.", sizeof(hmc_float) * 8, prec);
	}

	//write xlf-info to string, should look like this
	/*
	char xlf_info [] = "plaquette = 6.225960e-01\n trajectory nr = 67336\n beta = 6.000000, kappa = 0.177000, mu = 0.500000, c2_rec = 0.000000\n time = 1278542490\n hmcversion = 5.1.5\n mubar = 0.000000\n epsilonbar = 0.000000\n date = Thu Jul  8 00:41:30 2010\n";
	*/

	//note: it is not clear what "time" is supposed to be
	char dummystring[1000];
	char xlf_info[1000];
	sprintf(xlf_info, "%s", "plaquette = ");
	sprintf(dummystring, "%f", plaquettevalue);
	strcat(xlf_info, dummystring);
	sprintf(dummystring, "%s", "\n trajectory nr = ");
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n trajectory nr = ");
	sprintf(dummystring, "%i", trajectorynr);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n beta = ");
	sprintf(dummystring, "%f", beta);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, ", kappa = ");
	sprintf(dummystring, "%f", kappa);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, ", mu = ");
	sprintf(dummystring, "%f", mu);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, ", c2_rec = ");
	sprintf(dummystring, "%f", c2_rec);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n time = ");
	sprintf(dummystring, "%i", (int) current_time);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n hmcversion = ");
	sprintf(dummystring, "%s", hmc_version);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n mubar = ");
	sprintf(dummystring, "%f", mubar);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n epsilonbar = ");
	sprintf(dummystring, "%f", epsilonbar);
	strcat(xlf_info, dummystring);
	strcat(xlf_info, "\n date = ");
	sprintf(dummystring, "%s", date);
	strcat(xlf_info, dummystring);

	length_xlf_info = strlen(xlf_info);

	//write scidac checksum, this is stubb
	std::string scidac_checksum;
	{
		std::ostringstream tmp;
		tmp << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<scidacChecksum>\n<version>1.0</version>\n";
		tmp << "<suma>" << std::hex << checksum.get_suma() << "</suma>\n";
		tmp << "<sumb>" << std::hex << checksum.get_sumb() << "</sumb>\n";
		tmp << "</scidacChecksum>";
		scidac_checksum = tmp.str();
	}
	length_scidac_checksum = scidac_checksum.length();

	//write ildg_format to string, should look like this:
	/*
	char ildg_format [] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version>1.0</version>\n  <field>su3gauge</field>\n  <precision>64</precision>\n  <lx>4</lx>\n  <ly>4</ly>\n  <lz>4</lz>\n  <lt>4</lt>\n</ildgFormat>";
	*/
	char dummystring2[1000];
	char ildg_format[1000];
	sprintf(ildg_format, "%s", "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<ildgFormat xmlns=\"http://www.lqcd.org/ildg\"\n            xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n            xsi:schemaLocation=\"http://www.lqcd.org/ildg filefmt.xsd\">\n  <version>1.0</version>\n  <field>");
	strcat(ildg_format, field_out);
	strcat(ildg_format, "</field>\n  <precision>");
	sprintf(dummystring2, "%i", prec);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</precision>\n  <lx>");
	sprintf(dummystring2, "%i", lx);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</lx>\n  <ly>");
	sprintf(dummystring2, "%i", ly);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</ly>\n  <lz>");
	sprintf(dummystring2, "%i", lz);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</lz>\n  <lt>");
	sprintf(dummystring2, "%i", lt);
	strcat(ildg_format, dummystring2);
	strcat(ildg_format, "</lt>\n</ildgFormat>");

	length_ildg_format = strlen(ildg_format);


	//write the lime file
	MB_flag = 1;

	LimeWriter *writer;
	writer = limeCreateWriter(outputfile);

	const char * types[] = {"xlf-info", "ildg-format", "ildg-binary-data", "scidac-checksum"};

	//xlf-info
	ME_flag = 1;
	header_xlf_info = limeCreateHeader(MB_flag, ME_flag, (char*) types[0], length_xlf_info);
	limeWriteRecordHeader(header_xlf_info, writer);
	limeDestroyHeader(header_xlf_info);
	limeWriteRecordData( xlf_info, &length_xlf_info, writer);
	logger.debug() << "  xlf-info written";

	//ildg-format
	ME_flag = 2;
	header_ildg_format = limeCreateHeader(MB_flag, ME_flag, (char*) types[1], length_ildg_format);
	limeWriteRecordHeader(header_ildg_format, writer);
	limeDestroyHeader(header_ildg_format);
	limeWriteRecordData( ildg_format, &length_ildg_format, writer);
	logger.debug() << "  ildg-format written";

	//binary data
	ME_flag = 3;
	header_ildg_binary_data = limeCreateHeader(MB_flag, ME_flag, (char*) types[2], num_bytes);
	limeWriteRecordHeader(header_ildg_binary_data, writer);
	limeDestroyHeader(header_ildg_binary_data);
	limeWriteRecordData(binary_data, &num_bytes, writer);
	logger.debug() << "  ildg_binary_data written";

	//scidac-checksum
	ME_flag = 4;
	header_scidac_checksum = limeCreateHeader(MB_flag, ME_flag, (char*) types[3], length_scidac_checksum);
	limeWriteRecordHeader(header_scidac_checksum, writer);
	limeDestroyHeader(header_scidac_checksum);
	limeWriteRecordData(const_cast<char*>(scidac_checksum.c_str()), &length_scidac_checksum, writer);
	logger.debug() << "  scidac-checksum written";

	//closing
	fclose(outputfile);
	limeDestroyWriter(writer);
	logger.info() << "  " << (float) ( (float) (length_xlf_info + length_ildg_format + num_bytes + length_scidac_checksum) / 1024 / 1024 ) << " MBytes were written to the lime file " << filename;
}
