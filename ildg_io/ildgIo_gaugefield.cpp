/*
 * Copyright 2012, 2013, 2014 Lars Zeidlewicz, Christopher Pinke,
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

#include "ildgIo_gaugefield.hpp"

#include "limeUtilities.hpp"
#include "../host_functionality/logger.hpp"
#include "../host_functionality/host_geometry.h"
#include "../meta/util.hpp"
#include "../meta/version.hpp"

using namespace ildgIo;

void checkLimeFileForFieldType(std::string fieldTypeIn) throw(std::logic_error)
{
	if (fieldTypeIn != "su3gauge")
	{
		throw std::logic_error("LIME file does not seem to include gaugefield data. Aborting...");
	}
}

size_t sizeOfGaugefieldBuffer(size_t num_entries)
{
	return num_entries * sizeof(hmc_float);
}

//todo: make char ** std::vector<char*>
IldgIoReader_gaugefield::IldgIoReader_gaugefield(std::string sourceFilenameIn, const meta::Inputparameters * parametersIn, Matrixsu3 ** destination) : LimeFileReader(sourceFilenameIn, parametersIn->get_precision())
{
	if ( limeFileProp.numberOfBinaryDataEntries >= 1 )
	{
		*destination = new Matrixsu3[getNumberOfElements_gaugefield(parametersIn)];
	
		char * gf_ildg;
		
		checkLimeFileForFieldType(parameters.field);
		size_t numberOfBytes = sizeOfGaugefieldBuffer(parameters.num_entries);
		
		extractDataFromLimeFile(&gf_ildg, numberOfBytes);
		
		Checksum checksum = ildgIo::calculate_ildg_checksum(gf_ildg, parameters.getSizeInBytes(), *parametersIn);
	
		ildgIo::copy_gaugefield_from_ildg_format(*destination, gf_ildg, parameters.num_entries, *parametersIn);
		
		delete[] gf_ildg;
		
		parameters.checkAgainstChecksum(checksum, parametersIn->get_ignore_checksum_errors(), sourceFilenameIn);
		parameters.checkAgainstInputparameters(parametersIn);
	}
	else 
	{
		throw std::logic_error("LIME file does not seem to include binary data. Aborting...");
	}
}

void* createVoidPointerFromString(std::string stringIn) noexcept
{
	return (void*) stringIn.c_str();
}

IldgIoWriter_gaugefield::IldgIoWriter_gaugefield(Matrixsu3 * data, const meta::Inputparameters * parameters, std::string filenameIn, int trajectoryNumber, double plaquetteValue): LimeFileWriter(filenameIn)
{
	logger.info() << "writing gaugefield to lime-file \""  + filenameIn + "\"";
	
	size_t numberOfElements = getNumberOfElements_gaugefield(parameters);
	n_uint64_t num_bytes = getSizeInBytes_gaugefield(numberOfElements);
	char * binary_data = new char[num_bytes];
	copy_gaugefield_to_ildg_format(binary_data, data, *parameters);
	
	const Checksum checksum = calculate_ildg_checksum(binary_data, num_bytes, *parameters);

	Sourcefileparameters srcFileParameters(parameters, trajectoryNumber, plaquetteValue, checksum, version);	
	
	std::string xlfInfo = srcFileParameters.getInfo_xlfInfo();
	std::string scidac_checksum = srcFileParameters.getInfo_scidacChecksum();
	std::string ildgFormat = srcFileParameters.getInfo_ildgFormat_gaugefield();
	
	writeMemoryToLimeFile( createVoidPointerFromString(xlfInfo), xlfInfo.size(), limeEntryTypes["xlf"]);
	writeMemoryToLimeFile( createVoidPointerFromString(ildgFormat), ildgFormat.size(), limeEntryTypes["ildg"]);
	writeMemoryToLimeFile( binary_data, num_bytes, limeEntryTypes["ildg binary data"]);
	writeMemoryToLimeFile( createVoidPointerFromString(scidac_checksum), scidac_checksum.size(), limeEntryTypes["scidac checksum"]);
	
	delete[] binary_data;		
}

Checksum ildgIo::calculate_ildg_checksum(const char * buf, size_t nbytes, const meta::Inputparameters& inputparameters)
{
	const size_t elem_size = 4 * sizeof(Matrixsu3);

	const size_t NT = inputparameters.get_ntime();
	const size_t NS = inputparameters.get_nspace();

	if(nbytes != (NT * NS * NS * NS * elem_size)) {
		logger.error() << "Buffer does not contain a gaugefield!";
		throw Invalid_Parameters("Buffer size not match possible gaugefield size", (NT * NS * NS * NS * elem_size), nbytes);
	}

	Checksum checksum;

	size_t offset = 0;
	for(uint32_t t = 0; t < NT; ++t) {
		for(uint32_t z = 0; z < NS; ++z) {
			for(uint32_t y = 0; y < NS; ++y) {
				for(uint32_t x = 0; x < NS; ++x) {
					assert(offset < nbytes);
					uint32_t rank = ((t * NS + z) * NS + y) * NS + x;
					checksum.accumulate(&buf[offset], elem_size, rank);
					offset += elem_size;
				}
			}
		}
	}

	logger.debug() << "Calculated Checksum: " << checksum;
	return checksum;
}

static hmc_float make_float_from_big_endian(const char* in)
{
	hmc_float result;
	char * const raw_out = reinterpret_cast<char *>(&result);

	for(size_t i = 0; i < sizeof(hmc_float); ++i) {
#if BIG_ENDIAN_ARCH
		raw_out[i] = in[i];
#else
		raw_out[i] = in[sizeof(hmc_float) - 1 - i];
#endif
	}

	return result;
}

void ildgIo::copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, char * gaugefield_tmp, int check, const meta::Inputparameters& parameters)
{
	//little check if arrays are big enough
	if ((int) (meta::get_vol4d(parameters) *NDIM * NC * NC * 2) != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	const size_t NSPACE = parameters.get_nspace();
	int cter = 0;
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//save current link in a complex array
						hmc_complex tmp [NC][NC];
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								size_t pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								tmp[m][n].re = make_float_from_big_endian(&gaugefield_tmp[pos * sizeof(hmc_float)]);
								tmp[m][n].im = make_float_from_big_endian(&gaugefield_tmp[(pos + 1) * sizeof(hmc_float)]);
								cter++;
							}
						}
						//store su3matrix tmp in our format
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);

						//copy hmc_su3matrix to Matrixsu3 format
						Matrixsu3 destElem;
						destElem.e00 = tmp[0][0];
						destElem.e01 = tmp[0][1];
						destElem.e02 = tmp[0][2];
						destElem.e10 = tmp[1][0];
						destElem.e11 = tmp[1][1];
						destElem.e12 = tmp[1][2];
						destElem.e20 = tmp[2][0];
						destElem.e21 = tmp[2][1];
						destElem.e22 = tmp[2][2];

						gaugefield[get_global_link_pos((l + 1) % NDIM, spacepos, t, parameters)] = destElem;
					}
				}
			}
		}
	}

	if(cter * 2 != check) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values! there were " << cter * 2 << " vals set and not " << check << ".";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}
}

static void make_big_endian_from_float(char* out, const hmc_float in)
{
	char const * const raw_in = reinterpret_cast<char const *>(&in);

	for(size_t i = 0; i < sizeof(hmc_float); ++i) {
		out[i] = raw_in[sizeof(hmc_float) - 1 - i];
	}
}

void ildgIo::copy_gaugefield_to_ildg_format(char * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters)
{
	const size_t NSPACE = parameters.get_nspace();
	for (int t = 0; t < parameters.get_ntime(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						int coord[4];
						coord[0] = t;
						coord[1] = z;
						coord[2] = y;
						coord[3] = x;
						int spacepos = get_nspace(coord, parameters);
						hmc_complex destElem [NC][NC];

						Matrixsu3 srcElem = source_in[get_global_link_pos((l + 1) % NDIM, spacepos, t, parameters)];
						destElem[0][0] = srcElem.e00;
						destElem[0][1] = srcElem.e01;
						destElem[0][2] = srcElem.e02;
						destElem[1][0] = srcElem.e10;
						destElem[1][1] = srcElem.e11;
						destElem[1][2] = srcElem.e12;
						destElem[2][0] = srcElem.e20;
						destElem[2][1] = srcElem.e21;
						destElem[2][2] = srcElem.e22;

						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
								size_t pos = get_su3_idx_ildg_format(n, m, x, y, z, t, l, parameters);
								make_big_endian_from_float(&dest[pos * sizeof(hmc_float)], destElem[m][n].re);
								make_big_endian_from_float(&dest[(pos + 1) * sizeof(hmc_float)], destElem[m][n].im);
							}
						}
					}
				}
			}
		}
	}
}

size_t ildgIo::getNumberOfElements_gaugefield(const meta::Inputparameters * parameters)
{
	return meta::get_vol4d(*parameters) * NDIM;
}

n_uint64_t ildgIo::getSizeInBytes_gaugefield(size_t numberOfElements)
{
	return 2 * NC * NC * numberOfElements * sizeof(hmc_float);
}
