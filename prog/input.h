#ifndef _INPUTH_
#define _INPUTH_

#include "hmcerrs.h"
#include "types.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

class inputparameters {
 public:
  inputparameters(){set_defaults();};
  hmc_error readfile(char* ifn);
  hmc_error set_defaults();
  hmc_float get_kappa();
  hmc_float get_beta();
  hmc_float get_mu();
  int get_cgmax();
  int get_prec();
  int get_readsource();
  int get_thermalizationsteps();
  int get_heatbathsteps();
  void display_sourcefile();
  //CP
  //this is out of laziness
  std::string sourcefile;
 private:
  hmc_float kappa;
  hmc_float beta;
  hmc_float mu;
  int cgmax;
  int prec;
  int readsource;
  int thermalizationsteps;
  int heatbathsteps;
  void val_assign(hmc_float* out, std::string line);
  void val_assign(int * out, std::string line);
  void val_assign(std::string * out, std::string line);
};

#endif
