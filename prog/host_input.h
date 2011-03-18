#ifndef _INPUTH_
#define _INPUTH_

#include "hmcerrs.h"
#include "types.h"
#include "globaldefs.h"

#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <cstring>
using namespace std;

class inputparameters {
 public:
  inputparameters(){set_defaults();};
  hmc_error readfile(char* ifn);
  hmc_error set_defaults();
  hmc_float get_kappa();
  hmc_float get_beta();
  hmc_float get_theta_fermion();
  hmc_float get_theta_gaugefield();
  hmc_float get_mu();
  hmc_float get_chem_pot_re();
  hmc_float get_chem_pot_im();
  int get_cgmax();
  int get_prec();
  int get_startcondition();
  int get_thermalizationsteps();
  int get_heatbathsteps();
	int get_overrelaxsteps();
  int get_saveconfigs();
  int get_savefrequency();
  int get_writefrequency();
  void display_sourcefile();
  void display_sourcefilenumber();
  //CP
  //this is out of laziness
  std::string sourcefile;
  std::string sourcefilenumber;
 private:
  hmc_float kappa;
  hmc_float beta;
  hmc_float mu;
  hmc_float theta_fermion;
  hmc_float theta_gaugefield;
  hmc_float chem_pot_re;
  hmc_float chem_pot_im;
  int cgmax;
  int prec;
  int startcondition;
  int thermalizationsteps;
  int heatbathsteps;
	int overrelaxsteps;
  int savefrequency;
  int saveconfigs;
  int writefrequency;
  void val_assign(hmc_float* out, std::string line);
  void val_assign(int * out, std::string line);
  void sourcefilenumber_assign(std::string * out);
  void cond_assign(int * out, std::string line);
  void savecond_assign(int * out, std::string line);
  void val_assign(std::string * out, std::string line);
};

#endif
