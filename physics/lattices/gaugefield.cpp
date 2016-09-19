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
#include "../utilities.hpp"

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const GaugefieldParametersInterface * parameters, const physics::PRNG& prng)
  : system(system), prng(prng),  latticeObjectParameters(parameters), gaugefield(system)
{
	initializeBasedOnParameters();
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const GaugefieldParametersInterface * parameters, const physics::PRNG& prng, bool hot)
  : system(system), prng(prng), latticeObjectParameters(parameters), gaugefield(system)
{
	initializeHotOrCold(hot);
}

physics::lattices::Gaugefield::Gaugefield(const hardware::System& system, const GaugefieldParametersInterface * parameters, const physics::PRNG& prng, std::string ildgfile)
  : system(system), prng(prng),  latticeObjectParameters(parameters), gaugefield(system)
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
		gaugefield.set_hot();
		update_halo();
	} else {
		gaugefield.set_cold();
	}
	trajectoryNumberAtInit = 0;
}



void physics::lattices::Gaugefield::initializeFromILDGSourcefile(std::string ildgfile)
{
	Matrixsu3 * gf_host = ildgIo::readGaugefieldFromSourcefile(ildgfile, latticeObjectParameters, trajectoryNumberAtInit);

	gaugefield.send_gaugefield_to_buffers(gf_host);

	delete[] gf_host;
}

physics::lattices::Gaugefield::~Gaugefield()
{}

std::string physics::lattices::Gaugefield::getName(int number) const noexcept
{
	return physics::buildCheckpointName( latticeObjectParameters->getNamePrefix(), latticeObjectParameters->getNamePostfix(), latticeObjectParameters->getNumberOfDigitsInName(), number);
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
	gaugefield.fetch_gaugefield_from_buffers(host_buf);

	//http://stackoverflow.com/questions/2434196/how-to-initialize-stdvector-from-c-style-array
	std::vector<Matrixsu3> tmp(numberOfElements);
	tmp.assign(host_buf, host_buf + numberOfElements);
	
	ildgIo::writeGaugefieldToFile(outputfile, tmp, latticeObjectParameters, number);

	delete host_buf;
}

const std::vector<const hardware::buffers::SU3 *> physics::lattices::Gaugefield::get_buffers() const noexcept
{
	return gaugefield.get_buffers();
}

void physics::lattices::Gaugefield::smear()
{
	gaugefield.smear(latticeObjectParameters->getSmearingSteps());
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
	gaugefield.unsmear();
}

void physics::lattices::Gaugefield::update_halo() const
{
	gaugefield.update_halo();
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

