/** @file
 * Exception handling
 *
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

#ifndef _FILE_EXCEPTIONH_
#define _FILE_EXCEPTIONH_

#include "common_header_files/globaldefs.h"
#include <stdexcept>

#include <string>

class File_Exception : public std::runtime_error {
public:
	File_Exception(const char* name);
	File_Exception(std::string name);
	std::string get_filename();
	virtual const char* what() const noexcept override;
private:
	std::string filename;
};

class Invalid_Fermact : public std::invalid_argument {
public:
	Invalid_Fermact(int fermact, bool muset, bool cswset);
	virtual const char* what() const noexcept override;
private:
	std::string error_message;
};

class Invalid_Gaugeact : public std::invalid_argument {
public:
	Invalid_Gaugeact();
	virtual const char* what() const noexcept override;
private:
	std::string error_message;
};

class Invalid_Parameters : public std::invalid_argument {
public:
	Invalid_Parameters(std::string descr, std::string expected, std::string found);
	Invalid_Parameters(std::string descr, std::string expected, int found);
	Invalid_Parameters(std::string descr, std::string expected, double found);
	Invalid_Parameters(std::string descr, int expected, int found);
	virtual const char* what() const noexcept override;
private:
	std::string error_message;
};

class Opencl_Error : public std::runtime_error {
public:
	Opencl_Error(int clerr);
	Opencl_Error(int clerr, std::string clname);
	Opencl_Error(int clerr, std::string clname, std::string filename, int linenumber);
	virtual const char* what() const noexcept override;
private:
	std::string error_message;
};

class Print_Error_Message : public std::exception {
public:
	Print_Error_Message(std::string msg);
	Print_Error_Message(std::string msg, std::string filename, int linenumber);
	virtual const char* what() const noexcept override;
private:
	std::string error_message;
};

#endif
