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

#include "../host_functionality/logger.hpp"
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
IldgIoReader_gaugefield::IldgIoReader_gaugefield(std::string sourceFilenameIn, const IldgIoParameters * parametersIn, Matrixsu3 ** destination) : LimeFileReader(sourceFilenameIn, parametersIn->getPrecision())
{
	if ( limeFileProp.numberOfBinaryDataEntries >= 1 )
	{
		*destination = new Matrixsu3[parametersIn->getNumberOfElements()];
	
		char * gf_ildg;
		
		checkLimeFileForFieldType(parameters.field);
		size_t numberOfBytes = sizeOfGaugefieldBuffer(parameters.num_entries);
		
		extractDataFromLimeFile(&gf_ildg, numberOfBytes);
		
		Checksum checksum = ildgIo::calculate_ildg_checksum(gf_ildg, parameters.getSizeInBytes(), parametersIn->getNt(), parametersIn->getNs() );
	
		copy_gaugefield_from_ildg_format(*destination, gf_ildg, parameters.num_entries, *parametersIn);
		
		delete[] gf_ildg;
		
		parameters.checkAgainstChecksum(checksum, parametersIn->ignoreChecksumErrors(), sourceFilenameIn);
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

IldgIoWriter_gaugefield::IldgIoWriter_gaugefield(const std::vector<Matrixsu3> & data, const IldgIoParameters * parameters, std::string filenameIn, int trajectoryNumber, double plaquetteValue): LimeFileWriter(filenameIn)
{
	logger.info() << "writing gaugefield at tr. " << trajectoryNumber << " to lime-file \""  + filenameIn + "\"";

	n_uint64_t num_bytes = getSizeInBytes_gaugefield(parameters->getNumberOfElements());
	std::vector<char> binary_data(num_bytes);
	char * binary_data_ptr = &binary_data[0];
	
	copy_gaugefield_to_ildg_format(binary_data, data, *parameters);
	
	const Checksum checksum = calculate_ildg_checksum(binary_data_ptr, num_bytes, parameters->getNt(), parameters->getNs() );

	Sourcefileparameters srcFileParameters(parameters, trajectoryNumber, plaquetteValue, checksum, version);
	
	std::string xlfInfo = srcFileParameters.getInfo_xlfInfo();
	std::string scidac_checksum = srcFileParameters.getInfo_scidacChecksum();
	std::string ildgFormat = srcFileParameters.getInfo_ildgFormat_gaugefield();
	
	writeMemoryToLimeFile( createVoidPointerFromString(xlfInfo), xlfInfo.size(), limeEntryTypes["xlf"]);
	writeMemoryToLimeFile( createVoidPointerFromString(ildgFormat), ildgFormat.size(), limeEntryTypes["ildg"]);
	writeMemoryToLimeFile( binary_data_ptr, num_bytes, limeEntryTypes["ildg binary data"]);
	writeMemoryToLimeFile( createVoidPointerFromString(scidac_checksum), scidac_checksum.size(), limeEntryTypes["scidac checksum"]);

	closeLimeFile();
	logger.info() << "Checking that gaugefield was successfully written by re-reading it...";
	Matrixsu3 * gaugefieldTmp = NULL;
	IldgIoReader_gaugefield reader(filenameIn, parameters, &gaugefieldTmp);
	delete[] gaugefieldTmp;
	logger.info() << "...done";
}

Checksum ildgIo::calculate_ildg_checksum(const char * buf, size_t nbytes, const size_t NT, const size_t NS)
{
	const size_t elem_size = 4 * sizeof(Matrixsu3);

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

void ildgIo::copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, char * gaugefield_tmp, int expectedNumberOfEntries, const IldgIoParameters& parameters)
{
	//little check if arrays are big enough
	if ((int) ( parameters.getNumberOfElements() * NC * NC * 2) != expectedNumberOfEntries) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values!!\nCheck global settings!!";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}

	const size_t NSPACE = parameters.getNs();
	int cter = 0;
	for (int t = 0; t < parameters.getNt(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						//save current link in a complex array
						hmc_complex tmp [NC][NC];
						for (int m = 0; m < NC; m++) {
							for (int n = 0; n < NC; n++) {
	 							uint pos = LinkIndex(Index(z,y,x,t,LatticeExtents(parameters.getNs(),parameters.getNt())),static_cast<Direction>(l)).get_su3_idx_ildg_format(n,m);
								tmp[m][n].re = make_float_from_big_endian(&gaugefield_tmp[pos * sizeof(hmc_float)]);
								tmp[m][n].im = make_float_from_big_endian(&gaugefield_tmp[(pos + 1) * sizeof(hmc_float)]);
								cter++;
							}
						}

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

						//store su3matrix tmp in our format
						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						gaugefield[uint(LinkIndex(Index(z,y,x,t,LatticeExtents(parameters.getNs(),parameters.getNt())),static_cast<Direction>((l + 1) % NDIM)))] = destElem;
					}
				}
			}
		}
	}

	if(cter * 2 != expectedNumberOfEntries) {
		std::stringstream errstr;
		errstr << "Error in setting gaugefield to source values! there were " << cter * 2 << " values set and not " << expectedNumberOfEntries << ".";
		throw Print_Error_Message(errstr.str(), __FILE__, __LINE__);
	}
}

static void make_big_endian_from_float(char* out, const hmc_float in)
{
	char const * const raw_in = reinterpret_cast<char const *>(&in);

	for(size_t i = 0; i < sizeof(hmc_float); ++i) {
#if BIG_ENDIAN_ARCH
		out[i] = raw_in[i];
#else
		out[i] = raw_in[sizeof(hmc_float) - 1 - i];
#endif
	}
}

void ildgIo::copy_gaugefield_to_ildg_format(std::vector<char> & dest, const std::vector<Matrixsu3> & source_in, const IldgIoParameters& parameters)
{
	const size_t NSPACE = parameters.getNs();
	for (int t = 0; t < parameters.getNt(); t++) {
		for (size_t x = 0; x < NSPACE; x++) {
			for (size_t y = 0; y < NSPACE; y++) {
				for (size_t z = 0; z < NSPACE; z++) {
					for (int l = 0; l < NDIM; l++) {
						hmc_complex destElem [NC][NC];

						//our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
						//CP: interchange x<->z temporarily because spacepos has to be z + y * NSPACE + x * NSPACE * NSPACE!!
						Matrixsu3 srcElem = source_in[uint(LinkIndex(Index(z,y,x,t,LatticeExtents(parameters.getNs(),parameters.getNt())),static_cast<Direction>((l + 1) % NDIM)))];
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
								uint pos = LinkIndex(Index(z,y,x,t,LatticeExtents(parameters.getNs(),parameters.getNt())),static_cast<Direction>(l)).get_su3_idx_ildg_format(n,m);
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
