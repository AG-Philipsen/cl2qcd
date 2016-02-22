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
#include "../../ildg_io/ildgIo.hpp"
#include "../../hardware/code/gaugefield.hpp"

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, const physics::PRNG& prng);
static void set_cold(std::vector<const hardware::buffers::SU3 *> buffers);
static void set_cold(Matrixsu3 * field, size_t elems);
static void set_hot(Matrixsu3 * field, const physics::PRNG& prng, size_t elems);

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const GaugefieldParametersInterface * parameters, const physics::PRNG& prng)
  : system(system), prng(prng),  latticeObjectParameters(parameters), gaugefield(system), buffers(gaugefield.allocate_buffers()), unsmeared_buffers()
{
	initializeBasedOnParameters();
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const GaugefieldParametersInterface * parameters, const physics::PRNG& prng, bool hot)
  : system(system), prng(prng), latticeObjectParameters(parameters), gaugefield(system), buffers(gaugefield.allocate_buffers()), unsmeared_buffers()
{
	initializeHotOrCold(hot);
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const GaugefieldParametersInterface * parameters, const physics::PRNG& prng, std::string ildgfile)
  : system(system), prng(prng),  latticeObjectParameters(parameters), gaugefield(system), buffers(gaugefield.allocate_buffers()), unsmeared_buffers()
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



void physics::lattices::Gaugefield::initializeFromILDGSourcefile(std::string ildgfile)
{
	Matrixsu3 * gf_host = ildgIo::readGaugefieldFromSourcefile(ildgfile, latticeObjectParameters, trajectoryNumberAtInit);

	gaugefield.send_gaugefield_to_buffers(buffers, gf_host);

	delete[] gf_host;
}

physics::lattices::Gaugefield::~Gaugefield()
{
	gaugefield.release_buffers(&buffers);
	gaugefield.release_buffers(&unsmeared_buffers);
}

static void set_hot(std::vector<const hardware::buffers::SU3 *> buffers, const physics::PRNG& prng)
{
	using hardware::Device;

	for(auto buffer: buffers) 
	{
		size_t elems = buffer->get_elements();
		Matrixsu3 * tmp = new Matrixsu3[elems];
		set_hot(tmp, prng, elems);
		const Device * device = buffer->get_device();
		device->getGaugefieldCode()->importGaugefield(buffer, tmp);
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
		const Device * device = buffer->get_device();
		device->getGaugefieldCode()->importGaugefield(buffer, tmp);
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
	gaugefield.fetch_gaugefield_from_buffers(buffers, host_buf);

	//http://stackoverflow.com/questions/2434196/how-to-initialize-stdvector-from-c-style-array
	std::vector<Matrixsu3> tmp(numberOfElements);
	tmp.assign(host_buf, host_buf + numberOfElements);
	
	ildgIo::writeGaugefieldToFile(outputfile, tmp, latticeObjectParameters, number);

	delete host_buf;
}

const std::vector<const hardware::buffers::SU3 *> physics::lattices::Gaugefield::get_buffers() const noexcept
{
	return buffers;
}

void physics::lattices::Gaugefield::smear()
{
	unsmeared_buffers = gaugefield.allocate_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		auto device = buf->get_device();
		auto gf_code = device->getGaugefieldCode(); // should be like: get_gaugefield_code( HardwareParameters_gaugefield )

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

	unsmeared_buffers = gaugefield.allocate_buffers();

	for(size_t i = 0; i < buffers.size(); ++i) {
		auto buf = buffers[i];
		hardware::buffers::copyData(buf, unsmeared_buffers[i]);
	}

	gaugefield.release_buffers(&unsmeared_buffers);
}

void physics::lattices::Gaugefield::update_halo() const
{
	if(buffers.size() > 1) { // for a single device this will be a noop
		// currently either all or none of the buffers must be SOA
		if(buffers[0]->is_soa()) {
			gaugefield.update_halo_soa(buffers, system);
		} else {
			gaugefield.update_halo_aos(buffers, system);
		}
	}
}

const physics::PRNG * physics::lattices::Gaugefield::getPrng() const 
{
  return &prng;
}

const hardware::System * physics::lattices::Gaugefield::getSystem() const
{
  return &system;
}

const physics::lattices::GaugefieldParametersInterface * physics::lattices::Gaugefield::getParameters() const
{
  return latticeObjectParameters;
}

int physics::lattices::Gaugefield::get_trajectoryNumberAtInit() const
{
	return trajectoryNumberAtInit;
}

