#include "exceptions.h"

File_Exception::File_Exception(char* name) {
  filename = name;
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
  msg<< descr<<" Expeceted: "<<expected<<" But found: "<<found;
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
