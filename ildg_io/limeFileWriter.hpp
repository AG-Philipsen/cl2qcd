/** @file
 * Reader for LIME files
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

#ifndef _LIMEFILEWRITER_HPP_
#define _LIMEFILEWRITER_HPP_

#include <string>
#include "limeUtilities.hpp"

class LimeFileWriter
{
	//todo: make this protected
public:
	LimeFileWriter(std::string filenameIn);
	~LimeFileWriter();
	void writeMemoryToLimeFile(void * memoryPointer, n_uint64_t bytes, std::string description);
protected:
	LimeEntryTypes limeEntryTypes;
private:
	std::string filename;
	FILE *outputfile;
	int MB_flag;
	int ME_flag;
	n_uint64_t writtenBytes;
	LimeWriter *writer;
	void writeLimeHeaderToLimeFile(LimeRecordHeader * header);
};

#endif