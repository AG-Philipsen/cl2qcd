/** @file
 * Exception handling
 */

#ifndef _FILE_EXCEPTIONH_
#define _FILE_EXCEPTIONH_

#include "globaldefs.h"

#include <cstdlib>
#include <iostream>
#include <cstring>
#include <sstream>

class File_Exception {
 public:
  File_Exception(char* name);
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

class Invalid_Parameters  {
 public:
  Invalid_Parameters(std::string descr, std::string expected, std::string found);
  Invalid_Parameters(std::string descr, std::string expected, int found);
  Invalid_Parameters(std::string descr, int expected, int found);
  std::string what();
 private:
  std::string error_message;
};

#endif
