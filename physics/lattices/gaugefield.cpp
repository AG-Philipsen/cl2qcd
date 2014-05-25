/** @file
 * Implementation of the physics::lattices::Gaugefield class
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 *
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

#include "gaugefield.hpp"
#include "../../meta/version.hpp"
#include "../../meta/util.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../host_functionality/host_operations_gaugefield.h"
#include "../../ildg_io/ildg_write_gaugefield.h"
#include <cassert>
#include "../../ildg_io/checksum.h"
#include <fstream>
#include "../../hardware/device.hpp"
#include "../../hardware/buffers/halo_update.hpp"
#include "../observables/gaugeObservables.h"


/**
 * Version number.
 *
 * @deprecated move this into some specific header or so
 */
extern std::string const version;

static std::vector<const hardware::buffers::SU3 *> allocate_buffers(const hardware::System& system);
static void release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers);

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, const physics::PRNG& prng);
static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers);
static void set_cold(Matrixsu3 * field, size_t elems);
static void set_hot(Matrixsu3 * field, const physics::PRNG& prng, size_t elems);
static void copy_gaugefield_to_ildg_format(char * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters);
static void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, char * gaugefield_tmp, int check, const meta::Inputparameters& parameters);
static void check_sourcefileparameters(const meta::Inputparameters& parameters, const hmc_float, sourcefileparameters& parameters_source);
static Checksum calculate_ildg_checksum(const char * buf, size_t nbytes, const meta::Inputparameters& inputparameters);
static hmc_float make_float_from_big_endian(const char* in);
static void make_big_endian_from_float(char* out, const hmc_float in);
static void send_gaugefield_to_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, const Matrixsu3 * const gf_host, const meta::Inputparameters& params);
static void fetch_gaugefield_from_buffers(Matrixsu3 * const gf_host, const std::vector<const hardware::buffers::SU3 *> buffers, const meta::Inputparameters& params);
static void update_halo_soa(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system);
static void update_halo_aos(const std::vector<const hardware::buffers::SU3 *> buffers, const meta::Inputparameters& params);

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const physics::PRNG& prng)
  : system(system), prng(prng), buffers(allocate_buffers(system)), unsmeared_buffers(),  parameters(&system.get_inputparameters()) , parameters_source()
{
	initializeBasedOnParameters();
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const physics::PRNG& prng, bool hot)
  : system(system), prng(prng), buffers(allocate_buffers(system)), unsmeared_buffers(), parameters(&system.get_inputparameters()), parameters_source() 
{
	initializeHotOrCold(hot);
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const physics::PRNG& prng, std::string ildgfile)
  : system(system), prng(prng), buffers(allocate_buffers(system)), unsmeared_buffers(),  parameters(&system.get_inputparameters()), parameters_source() 
{
	initializeFromILDGSourcefile(ildgfile);
}

void physics::lattices::Gaugefield::initializeBasedOnParameters()
{
	switch(parameters->get_startcondition()) {
		case meta::Inputparameters::start_from_source:
			initializeFromILDGSourcefile(parameters->get_sourcefile());
			break;
		case meta::Inputparameters::cold_start:
			set_cold(buffers);
			break;
		case meta::Inputparameters::hot_start:
			set_hot(buffers, prng);
			update_halo();
			break;
	}
}

void physics::lattices::Gaugefield::initializeHotOrCold(bool hot)
{
	if(hot) {
		set_hot(buffers, prng);
		update_halo();
	} else {
		set_cold(buffers);
	}
}

void physics::lattices::Gaugefield::initializeFromILDGSourcefile(std::string ildgfile)
{
	auto parameters = system.get_inputparameters();
	Matrixsu3 * gf_host = new Matrixsu3[meta::get_vol4d(parameters) * 4];

	char * gf_ildg; // filled by readsourcefile
	parameters_source.readsourcefile(ildgfile.c_str(), parameters.get_precision(), &gf_ildg);

	Checksum checksum = calculate_ildg_checksum(gf_ildg, parameters_source.num_entries_source * sizeof(hmc_float), parameters);
	logger.debug() << "Calculated Checksum: " << checksum;

	if(checksum != parameters_source.checksum) {
		logger.error() << "Checksum of data does not match checksum given in file.";
		logger.error() << "Calculated Checksum: " << checksum;
		logger.error() << "Embedded Checksum:   " << parameters_source.checksum;
		if(!parameters.get_ignore_checksum_errors()) {
			throw File_Exception(ildgfile);
		}
	}

	copy_gaugefield_from_ildg_format(gf_host, gf_ildg, parameters_source.num_entries_source, parameters);

	send_gaugefield_to_buffers(buffers, gf_host, parameters);

	delete[] gf_ildg;
	delete[] gf_host;

	physics::observables::gaugeObservables obs(&parameters);
	hmc_float plaq = obs.measurePlaquette(this);
	check_sourcefileparameters(parameters, plaq, parameters_source);
}

static std::vector<const hardware::buffers::SU3 *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::SU3;

	std::vector<const SU3 *> buffers;

	auto const devices = system.get_devices();
	for(auto device: devices) {
		buffers.push_back(new SU3(get_vol4d(device->get_mem_lattice_size()) * 4, device));
	}
	return buffers;
}

static void release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers)
{
for(auto buffer: *buffers) {
		delete buffer;
	}
	buffers->clear();
}

physics::lattices::Gaugefield::~Gaugefield()
{
	release_buffers(&buffers);
	release_buffers(&unsmeared_buffers);
}

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, const physics::PRNG& prng)
{
	using hardware::Device;

for(auto buffer: buffers) {
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_hot(tmp, prng, elems);
		Device * device = buffer->get_device();
		device->get_gaugefield_code()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers)
{
	using hardware::Device;

for(auto buffer: buffers) {
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_cold(tmp, elems);
		Device * device = buffer->get_device();
		device->get_gaugefield_code()->importGaugefield(buffer, tmp);
		device->synchronize();
		delete[] tmp;
	}
}

void set_cold(Matrixsu3 * field, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = unit_matrixsu3();
	}
}

void set_hot(Matrixsu3 * field, const physics::PRNG& prng, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = random_matrixsu3(prng);
	}
}

void physics::lattices::Gaugefield::save(int number)
{
	std::string outputfile = meta::create_configuration_name(system.get_inputparameters());
	save(outputfile, number);
}

void physics::lattices::Gaugefield::saveToSpecificFile(int number)
{
	std::string outputfile = meta::create_configuration_name(system.get_inputparameters(), number);
	save(outputfile, number);
}

void physics::lattices::Gaugefield::save(std::string outputfile, int number)
{
	logger.info() << "saving current gauge configuration to file \"" << outputfile << "\"";
	auto parameters = system.get_inputparameters();
	const size_t NTIME = parameters.get_ntime();
	const size_t gaugefield_buf_size = 2 * NC * NC * NDIM * meta::get_volspace(parameters) * NTIME * sizeof(hmc_float);
	char * gaugefield_buf = new char[gaugefield_buf_size];

	//these are not yet used...
	hmc_float c2_rec = 0, epsilonbar = 0, mubar = 0;

	{
		Matrixsu3 * host_buf = new Matrixsu3[meta::get_vol4d(parameters) * NDIM];
		fetch_gaugefield_from_buffers(host_buf, buffers, parameters);
		copy_gaugefield_to_ildg_format(gaugefield_buf, host_buf, parameters);
		delete host_buf;
	}

	physics::observables::gaugeObservables obs(&parameters);
	hmc_float plaq = obs.measurePlaquette(this);

	const size_t NSPACE = parameters.get_nspace();

	const Checksum checksum = calculate_ildg_checksum(gaugefield_buf, gaugefield_buf_size, parameters);

	write_gaugefield(gaugefield_buf, gaugefield_buf_size, checksum, NSPACE, NSPACE, NSPACE, NTIME, parameters.get_precision(), number, plaq, parameters.get_beta(), parameters.get_kappa(), parameters.get_mu(), c2_rec, epsilonbar, mubar, version.c_str(), outputfile.c_str());

	delete[] gaugefield_buf;
}

static void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, char * gaugefield_tmp, int check, const meta::Inputparameters& parameters)
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

static void copy_gaugefield_to_ildg_format(char * dest, Matrixsu3 * source_in, const meta::Inputparameters& parameters)
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

static void check_sourcefileparameters(const meta::Inputparameters& parameters, const hmc_float plaquette, sourcefileparameters& parameters_source)
{
	logger.info() << "Checking sourcefile parameters against inputparameters...";
	//checking major parameters
	int tmp1, tmp2;
	std::string testobj;
	std::string msg = "Major parameters do not match: ";

	testobj = msg + "NT";
	tmp1 = parameters.get_ntime();
	tmp2 = parameters_source.lt_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "NX";
	tmp1 = parameters.get_nspace();
	tmp2 = parameters_source.lx_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "NY";
	tmp1 = parameters.get_nspace();
	tmp2 = parameters_source.ly_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "NZ";
	tmp1 = parameters.get_nspace();
	tmp2 = parameters_source.lz_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}
	testobj = msg + "PRECISION";
	tmp1 = parameters.get_precision();
	tmp2 = parameters_source.prec_source;
	if(tmp1 != tmp2) {
		throw Invalid_Parameters(testobj , tmp1, tmp2);
	}

	//checking minor parameters
	msg = "Minor parameters do not match: ";
	hmc_float float1, float2;
	testobj = msg + "plaquette";
	float1 = plaquette;
	float2 = parameters_source.plaquettevalue_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}
	testobj = msg + "beta";
	float1 = parameters.get_beta();
	float2 = parameters_source.beta_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}
	testobj = msg + "kappa";
	float1 = parameters.get_kappa();
	float2 = parameters_source.kappa_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}
	testobj = msg + "mu";
	float1 = parameters.get_mu();
	float2 = parameters_source.mu_source;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}

	logger.info() << "...done";
	return;
}

const std::vector<const hardware::buffers::SU3 *> physics::lattices::Gaugefield::get_buffers() const noexcept
{
	return buffers;
}

static Checksum calculate_ildg_checksum(const char * buf, size_t nbytes, const meta::Inputparameters& inputparameters)
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

	return checksum;
}

static hmc_float make_float_from_big_endian(const char* in)
{
	union {
		char b[sizeof(hmc_float)];
		hmc_float f;
	} val;

	for(size_t i = 0; i < sizeof(hmc_float); ++i) {
		val.b[i] = in[sizeof(hmc_float) - 1 - i];
	}
	return val.f;
}

static void make_big_endian_from_float(char* out, const hmc_float in)
{
	union {
		char b[sizeof(hmc_float)];
		hmc_float f;
	} val;

	val.f = in;

	for(size_t i = 0; i < sizeof(hmc_float); ++i) {
		out[i] = val.b[sizeof(hmc_float) - 1 - i];
	}
}

void physics::lattices::Gaugefield::smear()
{
	auto parameters = system.get_inputparameters();

	unsmeared_buffers = allocate_buffers(system);

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		auto device = buf->get_device();
		auto gf_code = device->get_gaugefield_code();

		hardware::buffers::copyData(unsmeared_buffers[i], buf);

		int rho_iter = parameters.get_rho_iter();
		logger.debug() << "\t\tperform " << rho_iter << " steps of stout-smearing to the gaugefield...";

		//one needs a temporary gf to apply the smearing to
		const hardware::buffers::SU3 gf_tmp(buf->get_elements(), device);
		for(int i = 0; i < rho_iter - 1; i += 2) {
			gf_code->stout_smear_device(buf, &gf_tmp);
			gf_code->stout_smear_device(&gf_tmp, buf);
		}
		//if rho_iter is odd one has to copy ones more
		if(rho_iter % 2 == 1) {
			gf_code->stout_smear_device(buf, &gf_tmp);
			hardware::buffers::copyData(buf, &gf_tmp);
		}
	}
}

void physics::lattices::Gaugefield::smear() const
{
  this->smear();
}

void physics::lattices::Gaugefield::unsmear() const
{
  this->unsmear();
}

void physics::lattices::Gaugefield::unsmear()
{
	if(unsmeared_buffers.size() == 0) {
		logger.warn() << "Tried to unsmear gaugefield that is not smeared.";
		return;
	}

	auto parameters = system.get_inputparameters();

	unsmeared_buffers = allocate_buffers(system);

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		hardware::buffers::copyData(buf, unsmeared_buffers[i]);
	}

	release_buffers(&unsmeared_buffers);
}

sourcefileparameters physics::lattices::Gaugefield::get_parameters_source()
{
  return parameters_source;
}

static void send_gaugefield_to_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, const Matrixsu3 * const gf_host, const meta::Inputparameters& params) {
	if(buffers.size() == 1) {
		auto device = buffers[0]->get_device();
		device->get_gaugefield_code()->importGaugefield(buffers[0], gf_host);
		device->synchronize();
	} else {
		auto const _device = buffers.at(0)->get_device();
		auto const local_size = _device->get_local_lattice_size();
		size_4 const halo_size(local_size.x, local_size.y, local_size.z, _device->get_halo_size());
		auto const grid_size = _device->get_grid_size();
		if(grid_size.x != 1 || grid_size.y != 1 || grid_size.z != 1) {
			throw Print_Error_Message("Not implemented!", __FILE__, __LINE__);
		}
		for(auto const buffer: buffers) {
			auto device = buffer->get_device();
			Matrixsu3 * mem_host = new Matrixsu3[buffer->get_elements()];

			size_4 offset(0, 0, 0, device->get_grid_pos().t * local_size.t);
			logger.debug() << offset;
			const size_t local_volume = get_vol4d(local_size) * NDIM;
			memcpy(mem_host, &gf_host[get_global_link_pos(0, offset, params)], local_volume * sizeof(Matrixsu3));

			const size_t halo_volume = get_vol4d(halo_size) * NDIM;
			size_4 halo_offset(0, 0, 0, (offset.t + local_size.t) % params.get_ntime());
			logger.debug() << halo_offset;
			memcpy(&mem_host[local_volume], &gf_host[get_global_link_pos(0, halo_offset, params)], halo_volume * sizeof(Matrixsu3));

			halo_offset = size_4(0, 0, 0, (offset.t + params.get_ntime() - halo_size.t) % params.get_ntime());
			logger.debug() << halo_offset;
			memcpy(&mem_host[local_volume + halo_volume], &gf_host[get_global_link_pos(0, halo_offset, params)], halo_volume * sizeof(Matrixsu3));

			device->get_gaugefield_code()->importGaugefield(buffer, mem_host);

			delete[] mem_host;
		}
	}
}

static void fetch_gaugefield_from_buffers(Matrixsu3 * const gf_host, const std::vector<const hardware::buffers::SU3 *> buffers, const meta::Inputparameters& params)
{
	if(buffers.size() == 1) {
		auto device = buffers[0]->get_device();
		device->get_gaugefield_code()->exportGaugefield(gf_host, buffers[0]);
		device->synchronize();
	} else {
		auto const _device = buffers.at(0)->get_device();
		auto const local_size = _device->get_local_lattice_size();
		auto const grid_size = _device->get_grid_size();
		if(grid_size.x != 1 || grid_size.y != 1 || grid_size.z != 1) {
			throw Print_Error_Message("Not implemented!", __FILE__, __LINE__);
		}
		for(auto const buffer: buffers) {
			// fetch local part for each device
			auto device = buffer->get_device();
			Matrixsu3 * mem_host = new Matrixsu3[buffer->get_elements()];

			device->get_gaugefield_code()->exportGaugefield(mem_host, buffer);
			size_4 offset(0, 0, 0, device->get_grid_pos().t * local_size.t);
			const size_t local_volume = get_vol4d(local_size) * NDIM;
			memcpy(&gf_host[get_global_link_pos(0, offset, params)], mem_host, local_volume * sizeof(Matrixsu3));

			delete[] mem_host;
		}
	}
}

void physics::lattices::Gaugefield::update_halo() const
{
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			update_halo_soa(buffers, system);
		} else {
			update_halo_aos(buffers, system.get_inputparameters());
		}
	}
}

static void update_halo_aos(const std::vector<const hardware::buffers::SU3 *> buffers, const meta::Inputparameters& params)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<Matrixsu3>(buffers, params, NDIM);
}

static void update_halo_soa(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(!buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo_soa<Matrixsu3>(buffers, system, .5, 2 * NDIM);
}

const physics::PRNG * physics::lattices::Gaugefield::getPrng() const 
{
  return &prng;
}

const hardware::System * physics::lattices::Gaugefield::getSystem() const
{
  return &system;
}

const meta::Inputparameters * physics::lattices::Gaugefield::getParameters() const
{
  return parameters;
}

