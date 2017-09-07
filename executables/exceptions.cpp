/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
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

#include "exceptions.h"

#include "../meta/inputparameters.hpp"
#include <sstream>

File_Exception::File_Exception(const char* name) : std::runtime_error(name)
{
	filename = name;
	return;
}

File_Exception::File_Exception(std::string name) : std::runtime_error(name)
{
	filename = name.c_str();
	return;
}

std::string File_Exception::get_filename()
{
	return filename;
}
const char* File_Exception::what() const noexcept
{
	return filename.c_str();
}

Invalid_Fermact::Invalid_Fermact(int fermact, bool muset, bool cswset) : std::invalid_argument("Invalid fermion action.")
{

	error_message = "Invalid setting for fermionic sector: ";

	if( muset == true && fermact != common::action::twistedmass )
		error_message.append("Value of twisted mass parameter mu is set but fermion action is not twisted mass.");

	if( cswset == true && fermact != common::action::clover )
		error_message.append("Value of clover parameter csw is set but fermion action is not clover. ");

	if( cswset == true && muset == true )
		error_message.append("Values for both csw and mu are set. This is not possible in this program. ");

	return;
}

const char* Invalid_Fermact::what() const noexcept
{
	return error_message.c_str();
}


Invalid_Gaugeact::Invalid_Gaugeact() : std::invalid_argument("Invalid gauge action.")
{

	error_message = "Invalid setting for gauge sector: ";

	return;
}

const char* Invalid_Gaugeact::what() const noexcept
{
	return error_message.c_str();
}

Invalid_Parameters::Invalid_Parameters(std::string descr, std::string expected, std::string found) : std::invalid_argument(descr)
{
	std::stringstream msg;
	msg << descr << " Expected: " << expected << " But found: " << found;
	error_message = msg.str();
	return;
}

Invalid_Parameters::Invalid_Parameters(std::string descr, std::string expected, int found) : std::invalid_argument(descr)
{
	std::stringstream msg;
	msg << descr << " Expected: " << expected << " But found: " << found;
	error_message = msg.str();
	return;
}

Invalid_Parameters::Invalid_Parameters(std::string descr, std::string expected, double found) : std::invalid_argument(descr)
{
	std::stringstream msg;
	msg << descr << " Expected: " << expected << " But found: " << found;
	error_message = msg.str();
	return;
}

Invalid_Parameters::Invalid_Parameters(std::string descr, int expected, int found) : std::invalid_argument(descr)
{
	std::stringstream msg;
	msg << descr << " Expected: " << expected << " But found: " << found;
	error_message = msg.str();
	return;
}

const char* Invalid_Parameters::what() const noexcept
{
	return error_message.c_str();
}


Opencl_Error::Opencl_Error(int clerr) : std::runtime_error("OpenCL Error")
{
	std::stringstream msg;
	msg << "OpenCL reported an error, error code: " << clerr;
	error_message = msg.str();
	return;
}

Opencl_Error::Opencl_Error(int clerr, std::string clname) : std::runtime_error("OpenCL Error")
{
	std::stringstream msg;
	msg << "OpenCL reported an error in " << clname << ", error code: " << clerr;
	error_message = msg.str();
	return;
}

Opencl_Error::Opencl_Error(int clerr, std::string clname, std::string filename, int linenumber) : std::runtime_error("OpenCL Error")
{
	std::stringstream msg;
	msg << "OpenCL failed. Error code " << clerr << " in " << clname << " at " << filename << ":" << linenumber;
	error_message = msg.str();
	return;
}

const char* Opencl_Error::what() const noexcept
{
	return error_message.c_str();
}


Print_Error_Message::Print_Error_Message(std::string msg)
{
	error_message = msg;
	return;
}

Print_Error_Message::Print_Error_Message(std::string msg, std::string filename, int linenumber)
{
	std::stringstream streammsg;
	streammsg << msg << " Occurred at " << filename << ":" << linenumber;
	error_message = streammsg.str();
	return;
}

const char* Print_Error_Message::what() const noexcept
{
	return error_message.c_str();
}


