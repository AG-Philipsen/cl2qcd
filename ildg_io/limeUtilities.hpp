/** @file
 * Utilities for handling LIME files
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

#ifndef _LIMEUTILITIES_HPP_
#define _LIMEUTILITIES_HPP_

#include <string>

extern "C" {
#include <lime.h>
}

class LimeHeaderData
{
public:
  LimeHeaderData(LimeReader *r);

  n_uint64_t numberOfBytes;
  size_t bytes_pad;
  int MB_flag, ME_flag;
  std::string limeEntryType;
};

class LimeFileProperties
{
public:
	LimeFileProperties();
	LimeFileProperties(int numberOfEntries,  int numberOfBinaryDataEntries);
	
  void operator+=(LimeFileProperties other);
  int numberOfEntries;
  int numberOfBinaryDataEntries;
	int numberOfFermionicEntries;
	bool readMetaData;
};

class LimeFilePropertiesCollector: public LimeFileProperties
{
 public:
  ~LimeFilePropertiesCollector();  
};

#endif