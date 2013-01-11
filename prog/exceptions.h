/** @file
 * Exception handling
 */

#ifndef _FILE_EXCEPTIONH_
#define _FILE_EXCEPTIONH_

#include "globaldefs.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>

class File_Exception {
public:
	File_Exception(const char* name);
	File_Exception(std::string name);
	std::string get_filename();
private:
	std::string filename;
};

class Invalid_Fermact {
public:
	Invalid_Fermact(int fermact, bool muset, bool cswset);
	std::string what();
private:
	std::string error_message;
};

class Invalid_Gaugeact {
public:
	Invalid_Gaugeact();
	std::string what();
private:
	std::string error_message;
};

class Invalid_Parameters  {
public:
	Invalid_Parameters(std::string descr, std::string expected, std::string found);
	Invalid_Parameters(std::string descr, std::string expected, int found);
	Invalid_Parameters(std::string descr, int expected, int found);
	std::string what();
private:
	std::string error_message;
};

class Opencl_Error {
public:
	Opencl_Error(int clerr);
	Opencl_Error(int clerr, std::string clname);
	Opencl_Error(int clerr, std::string clname, std::string filename, int linenumber);
	std::string what();
private:
	std::string error_message;
};

class Print_Error_Message {
public:
	Print_Error_Message(std::string msg);
	Print_Error_Message(std::string msg, std::string filename, int linenumber);
	std::string what();
private:
	std::string error_message;
};

#endif
