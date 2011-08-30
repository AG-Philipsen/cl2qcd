#include "exceptions.h"

File_Exception::File_Exception(char* name) {
  filename = name;
  return;
}

File_Exception::File_Exception(std::string name) {
  filename = name.c_str();
  return;
}

std::string File_Exception::get_filename(){
  return filename;
}


Invalid_Fermact::Invalid_Fermact(int fermact, bool muset, bool cswset){

  error_message = "Invalid setting for fermionic sector: ";

  if( muset == true && fermact != TWISTEDMASS ) 
    error_message.append("Value of twisted mass parameter mu is set but fermion action is not twisted mass.");

  if( cswset == true && fermact != CLOVER ) 
    error_message.append("Value of clover parameter csw is set but fermion action is not clover. ");

  if( cswset == true && muset == true ) 
    error_message.append("Values for both csw and mu are set. This is not possible in this program. ");

  return;
}

std::string Invalid_Fermact::what(){
  return error_message;
}

Invalid_Parameters::Invalid_Parameters(std::string descr, std::string expected, std::string found){
  std::stringstream msg;
  msg<< descr<<" Expected: "<<expected<<" But found: "<<found;
  error_message = msg.str();
  return;
}

Invalid_Parameters::Invalid_Parameters(std::string descr, std::string expected, int found){
  std::stringstream msg;
  msg<< descr<<" Expeceted: "<<expected<<" But found: "<<found;
  error_message = msg.str();
  return;
}

Invalid_Parameters::Invalid_Parameters(std::string descr, int expected, int found){
  std::stringstream msg;
  msg<< descr<<" Expeceted: "<<expected<<" But found: "<<found;
  error_message = msg.str();
  return;
}

std::string Invalid_Parameters::what(){
  return error_message;
}


Opencl_Error::Opencl_Error(int clerr){
  std::stringstream msg;
  msg<<"OpenCL reported an error, error code: "<<clerr;
  error_message = msg.str();
  return;
}

Opencl_Error::Opencl_Error(int clerr, std::string clname){
  std::stringstream msg;
  msg<<"OpenCL reported an error in "<<clname<<", error code: "<<clerr;
  error_message = msg.str();
  return;
}

Opencl_Error::Opencl_Error(int clerr, std::string clname, std::string filename, int linenumber){
  std::stringstream msg;
  msg<<"OpenCL failed. Error code "<<clerr<<" in "<<clname<<" at "<<filename<<":"<<linenumber;
  error_message = msg.str();
  return;
}

std::string Opencl_Error::what(){
  return error_message;
}


Print_Error_Message::Print_Error_Message(std::string msg){
  error_message = msg;
  return;
}

Print_Error_Message::Print_Error_Message(std::string msg, std::string filename, int linenumber){
  std::stringstream streammsg;
  streammsg<<msg<<" Occurred at "<<filename<<":"<<linenumber;
  error_message = streammsg.str();
  return;
}

std::string Print_Error_Message::what(){
  return error_message;
}


