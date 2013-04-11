/** @file
 * Exception handling
 */

#ifndef _FILE_EXCEPTIONH_
#define _FILE_EXCEPTIONH_

#include "globaldefs.h"
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
