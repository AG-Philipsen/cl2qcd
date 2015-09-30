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
#include "../../host_functionality/logger.hpp"
#include "../../host_functionality/host_operations_gaugefield.h"
#include "../../hardware/device.hpp"
#include "../../hardware/buffers/halo_update.hpp"
#include "../observables/gaugeObservables.h"
#include "../../ildg_io/ildgIo.hpp"
#include "../../hardware/code/gaugefield.hpp"

static std::vector<const hardware::buffers::SU3 *> allocate_buffers(const hardware::System& system);
static void release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers);

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, const physics::PRNG& prng);
static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers);
static void set_cold(Matrixsu3 * field, size_t elems);
static void set_hot(Matrixsu3 * field, const physics::PRNG& prng, size_t elems);
static void send_gaugefield_to_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, const Matrixsu3 * const gf_host, const LatticeObjectParametersInterface * params);
static void fetch_gaugefield_from_buffers(Matrixsu3 * const gf_host, const std::vector<const hardware::buffers::SU3 *> buffers, const LatticeObjectParametersInterface * params);
static void update_halo_soa(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system);
static void update_halo_aos(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system);

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const LatticeObjectParametersInterface * parameters, const physics::PRNG& prng)
  : system(system), prng(prng), buffers(allocate_buffers(system)), unsmeared_buffers(),  latticeObjectParameters(parameters)
{
	initializeBasedOnParameters();
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const LatticeObjectParametersInterface * parameters, const physics::PRNG& prng, bool hot)
  : system(system), prng(prng), buffers(allocate_buffers(system)), unsmeared_buffers(), latticeObjectParameters(parameters)
{
	initializeHotOrCold(hot);
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const LatticeObjectParametersInterface * parameters, const physics::PRNG& prng, std::string ildgfile)
  : system(system), prng(prng), buffers(allocate_buffers(system)), unsmeared_buffers(),  latticeObjectParameters(parameters)
{
	initializeFromILDGSourcefile(ildgfile);
}

void physics::lattices::Gaugefield::initializeBasedOnParameters()
{
	switch(latticeObjectParameters->getStartcondition()) {
		case common::startcondition::start_from_source:
			initializeFromILDGSourcefile(latticeObjectParameters->getSourcefileName());
			break;
		case common::startcondition::cold_start:
			initializeHotOrCold(false);
			break;
		case common::startcondition::hot_start:
			initializeHotOrCold(true);
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
	trajectoryNumberAtInit = 0;
}

//move to namespace ildg_io
static void check_plaq(const hmc_float plaquette, double plaqSourcefile)
{
	logger.info() << "Checking plaquette against sourcefile value...";
	std::string msg = "Minor parameters do not match: ";
	hmc_float float1, float2;
	std::string testobj = msg + "plaquette";
	float1 = plaquette;
	float2 = plaqSourcefile;
	if(float1 != float2) {
		logger.warn() << testobj;
		logger.warn() << "\tExpected: " << float1 << "\tFound: " << float2;
	}

	logger.info() << "...done";
	return;
}

void physics::lattices::Gaugefield::initializeFromILDGSourcefile(std::string ildgfile)
{
	double plaqSourcefile;

	Matrixsu3 * gf_host = ildgIo::readGaugefieldFromSourcefile(ildgfile, latticeObjectParameters, trajectoryNumberAtInit, plaqSourcefile);

	send_gaugefield_to_buffers(buffers, gf_host, latticeObjectParameters);

	delete[] gf_host;

	//todo: move this to ildgIo
	hmc_float plaq = physics::observables::measurePlaquette(this);
	check_plaq(plaq, plaqSourcefile);
}

static std::vector<const hardware::buffers::SU3 *> allocate_buffers(const hardware::System& system)
{
	using hardware::buffers::SU3;

	std::vector<const SU3 *> buffers;

	auto const devices = system.get_devices();
	for(auto device: devices) 
	{
		buffers.push_back(new SU3(get_vol4d(device->get_mem_lattice_size()) * 4, device));
	}
	return buffers;
}

static void release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers)
{
	for(auto buffer: *buffers) 
	{
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

	for(auto buffer: buffers) 
	{
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

	for(auto buffer: buffers) 
	{
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

std::string physics::lattices::getConfigurationName( std::string prefix, std::string postfix, int numberOfDigitsInName, int number)
{
	std::stringstream middle;
	if (number == -1)
	{
		middle << "save";
	}
	else
	{
		middle.fill('0');
		middle.width(numberOfDigitsInName);
		middle << std::right << number;
	}
	std::stringstream outfilename;
	outfilename << prefix << middle << postfix;
	return outfilename.str();
}

std::string physics::lattices::Gaugefield::getName(int number) const noexcept
{
	return getConfigurationName( latticeObjectParameters->getNamePrefix(), latticeObjectParameters->getNamePostfix(), latticeObjectParameters->getNumberOfDigitsInName(), number);
}

void physics::lattices::Gaugefield::save(int number)
{
	save(getName(), number);
}

void physics::lattices::Gaugefield::saveToSpecificFile(int number)
{
	save(getName(number), number);
}

void physics::lattices::Gaugefield::save(std::string outputfile, int number)
{
	logger.info() << "saving current gauge configuration to file \"" << outputfile << "\"";
	size_t numberOfElements = latticeObjectParameters->getNumberOfElements();
	Matrixsu3 * host_buf = new Matrixsu3[numberOfElements];
	fetch_gaugefield_from_buffers(host_buf, buffers, latticeObjectParameters);
	double plaq = physics::observables::measurePlaquette(this);

	//http://stackoverflow.com/questions/2434196/how-to-initialize-stdvector-from-c-style-array
	std::vector<Matrixsu3> tmp(numberOfElements);
	tmp.assign(host_buf, host_buf + numberOfElements);
	
	ildgIo::writeGaugefieldToFile(outputfile, tmp, latticeObjectParameters, number, plaq);

	delete host_buf;
}

const std::vector<const hardware::buffers::SU3 *> physics::lattices::Gaugefield::get_buffers() const noexcept
{
	return buffers;
}

void physics::lattices::Gaugefield::smear()
{
	unsmeared_buffers = allocate_buffers(system);

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		auto device = buf->get_device();
		auto gf_code = device->get_gaugefield_code();

		hardware::buffers::copyData(unsmeared_buffers[i], buf);

		int rho_iter = latticeObjectParameters->getSmearingSteps();
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

	unsmeared_buffers = allocate_buffers(system);

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		hardware::buffers::copyData(buf, unsmeared_buffers[i]);
	}

	release_buffers(&unsmeared_buffers);
}

static void send_gaugefield_to_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, const Matrixsu3 * const gf_host, const LatticeObjectParametersInterface * params) {
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
			memcpy(mem_host, &gf_host[get_global_link_pos(0, offset, params->getNt(), params->getNs())], local_volume * sizeof(Matrixsu3));

			const size_t halo_volume = get_vol4d(halo_size) * NDIM;
			size_4 halo_offset(0, 0, 0, (offset.t + local_size.t) % params->getNt());
			logger.debug() << halo_offset;
			memcpy(&mem_host[local_volume], &gf_host[get_global_link_pos(0, halo_offset,  params->getNt(), params->getNs())], halo_volume * sizeof(Matrixsu3));

			halo_offset = size_4(0, 0, 0, (offset.t + params->getNt() - halo_size.t) % params->getNt());
			logger.debug() << halo_offset;
			memcpy(&mem_host[local_volume + halo_volume], &gf_host[get_global_link_pos(0, halo_offset,  params->getNt(), params->getNs())], halo_volume * sizeof(Matrixsu3));

			device->get_gaugefield_code()->importGaugefield(buffer, mem_host);

			delete[] mem_host;
		}
	}
}

static void fetch_gaugefield_from_buffers(Matrixsu3 * const gf_host, const std::vector<const hardware::buffers::SU3 *> buffers, const LatticeObjectParametersInterface * params)
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
			memcpy(&gf_host[get_global_link_pos(0, offset, params->getNt(), params->getNs())], mem_host, local_volume * sizeof(Matrixsu3));

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
			update_halo_aos(buffers, system);
		}
	}
}

static void update_halo_aos(const std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system)
{
	// check all buffers are non-soa
	for(auto const buffer: buffers) {
		if(buffer->is_soa()) {
			throw Print_Error_Message("Mixed SoA-AoS configuration halo update is not implemented, yet.", __FILE__, __LINE__);
		}
	}

	hardware::buffers::update_halo<Matrixsu3>(buffers, system, NDIM);
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

const LatticeObjectParametersInterface * physics::lattices::Gaugefield::getParameters() const
{
  return latticeObjectParameters;
}

int physics::lattices::Gaugefield::get_trajectoryNumberAtInit() const
{
	return trajectoryNumberAtInit;
}

