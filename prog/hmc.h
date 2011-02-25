#ifndef _HMCH_
#define _HMCH_
//should only be included in main prog

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_geometry.h"
#include "host_testing.h"
#include "host_gaugeobservables.h"
#include "host_gaugefieldoperations.h"
#include "host_input.h"
#include "host_readgauge.h"
#include "host_random.h"
#include "host_update_heatbath.h"
#include "host_use_timer.h"
#include "host_writegaugefield.h"
#include "opencl.h"
#include <CL/cl.h>

#ifdef _OPENMP
# include <omp.h>
#endif

string const version = "0.1";

using namespace std;

//global random number generator
Random rnd (seed);

//couple of timers
usetimer totaltime;
usetimer inittime;
usetimer polytime;
usetimer plaqtime;
usetimer updatetime;
usetimer overrelaxtime;
usetimer copytime;

//to save gaugeobservables
hmc_float plaq, splaq, tplaq;
hmc_complex pol;

void print_hello(char* name){
  std::cout<<"This is hmc program "<<name<<", version "<<version<<"."<<endl;
  return;
}

void print_info(inputparameters* params){
  printf("**********************************************************\n");
  printf("Compile time parameters:\n");
  printf("NSPACE:  %d\n",NSPACE);
  printf("NTIME:   %d\n",NTIME);
  printf("NDIM:    %d\n",NDIM);
  printf("NCOLOR:  %d\n",NC);
  printf("NSPIN:   %d\n",NSPIN);
  printf("\n");
  printf("Run time parameters:\n");
  printf("kappa = %f\n",(*params).get_kappa());
  printf("mu    = %f\n",(*params).get_mu());
  printf("beta  = %f\n",(*params).get_beta());
  printf("CGmax = %d\n",(*params).get_cgmax());
  printf("prec = \t%d\n",(*params).get_prec());
  printf("thermsteps = \t%d\n",(*params).get_thermalizationsteps());
  printf("heatbathsteps = %d\n",(*params).get_heatbathsteps());
  printf("\n");
  if ((*params).get_startcondition()==START_FROM_SOURCE) {
    printf("sourcefile = ");
    (*params).display_sourcefile();
    printf("\n");
  }
  if ((*params).get_startcondition()==COLD_START) {
    printf("cold start\n");
  }
  if ((*params).get_startcondition()==HOT_START) {
    printf("hot start\n");
  }
  printf("**********************************************************\n");
  printf("\n");
  return;
}

#endif
