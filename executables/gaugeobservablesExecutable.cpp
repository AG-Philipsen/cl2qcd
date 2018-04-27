/*
 * Copyright (c) 2013-2015 Christopher Pinke
 * Copyright (c) 2016,2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#include "gaugeobservablesExecutable.hpp"

#include "../meta/parametersConfig.hpp"
#include "exceptions.hpp"

gaugeobservablesExecutable::gaugeobservablesExecutable(int argc, const char* argv[])
    : measurementExecutable(argc, argv, "gaugeobservables")
{
    if (parameters.get_startcondition() != common::startcondition::start_from_source) {
        throw std::logic_error("Invalid startcondition specified for gaugeobservables executable! Must be continue!");
    }
    printParametersToScreenAndFile();
}

void gaugeobservablesExecutable::writeGaugeobservablesLogfile()
{
    outputToFile.open(filenameForLogfile, std::ios::out | std::ios::app);
    if (outputToFile.is_open()) {
        meta::print_info_observables_gauge_io(&outputToFile, parameters);
        outputToFile.close();
    } else {
        throw File_Exception(filenameForLogfile);
    }
}

void gaugeobservablesExecutable::printParametersToScreenAndFile()
{
    meta::print_info_observables_gauge_io(parameters);
    writeGaugeobservablesLogfile();
}

void gaugeobservablesExecutable::performApplicationSpecificMeasurements()
{
    physics::observables::measureGaugeObservablesAndWriteToFile(gaugefield, gaugefield->get_trajectoryNumberAtInit(),
                                                                interfacesHandler
                                                                    ->getGaugeObservablesParametersInterface());
}
