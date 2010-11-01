#ifndef _INPUTH_
#define _INPUTH_

#include "hmcerrs.h"
#include "types.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

class inputparameters {
 public:
  inputparameters(){set_defaults();};
  hmc_error readfile(char* ifn);
  hmc_error set_defaults();
  hmc_float get_kappa();
  hmc_float get_beta();
  hmc_float get_mu();
  int get_cgmax();
 private:
  hmc_float kappa;
  hmc_float beta;
  hmc_float mu;
  int cgmax;
  void val_assign(hmc_float* out, std::string line);
  void val_assign(int * out, std::string line);
};

#endif
