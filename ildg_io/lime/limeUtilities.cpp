/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2018 Alessandro Sciarra
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

#include "limeUtilities.hpp"

#include "../executables/exceptions.hpp"
#include "../host_functionality/logger.hpp"

#include <boost/lexical_cast.hpp>

LimeHeaderData::LimeHeaderData(LimeReader* r)
{
    numberOfBytes = limeReaderBytes(r);
    limeEntryType = (limeReaderType(r));
    bytes_pad     = limeReaderPadBytes(r);
    MB_flag       = limeReaderMBFlag(r);
    ME_flag       = limeReaderMEFlag(r);
}

LimeFileProperties::LimeFileProperties()
    : numberOfEntries(0), numberOfBinaryDataEntries(0), numberOfFermionicEntries(0), readMetaData(false)
{
}

LimeFileProperties::LimeFileProperties(int numberOfEntries, int numberOfBinaryDataEntries)
    : numberOfEntries(numberOfEntries), numberOfBinaryDataEntries(numberOfBinaryDataEntries)
{
}

void LimeFileProperties::operator+=(LimeFileProperties other)
{
    this->numberOfEntries += other.numberOfEntries;
    this->numberOfBinaryDataEntries += other.numberOfBinaryDataEntries;
    this->numberOfFermionicEntries += other.numberOfFermionicEntries;
}

LimeFilePropertiesCollector::~LimeFilePropertiesCollector()
{
    logger.trace() << "Found " << numberOfEntries << " LIME records.";
    logger.trace() << "Found " << numberOfBinaryDataEntries << " binary entries in LIME file";
    if (numberOfFermionicEntries > 0) {
        logger.trace() << "\tfile contains " << numberOfFermionicEntries << " fermion entries.";
    } else {
        logger.trace() << "\tfile does not contain informations about fermions";
    }
}

LimeEntryTypes::Mapper LimeEntryTypes::mapper = {{"propagator", "propagator-info"},
                                                 {"xlf", "xlf-info"},
                                                 {"inverter", "inverter-info"},
                                                 {"gauge-checksum-copy", "gauge-scidac-checksum-copy"},
                                                 {"etmc-propagator", "etmc-propagator-format"},
                                                 {"scidac binary data", "scidac-binary-data"},
                                                 {"scidac checksum", "scidac-checksum"},
                                                 {"ildg", "ildg-format"},
                                                 {"ildg binary data", "ildg-binary-data"}};

void LimeFileReader_basic::openFile()
{
    outputfile = fopen(filename.c_str(), "r");
    limeReader = limeCreateReader(outputfile);
}

void LimeFileReader_basic::closeFile()
{
    limeDestroyReader(limeReader);
    fclose(outputfile);
}
