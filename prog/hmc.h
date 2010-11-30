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
#include "operations.h"
#include "geometry.h"
#include "testing.h"
#include "gaugeobservables.h"
#include "gaugefieldoperations.h"
#include "input.h"
#include "readgauge.h"
#include "random.h"
#include "update.h"
#include "use_timer.h"
#include "opencl.h"
#include <CL/cl.hpp>

#ifdef _OPENMP
# include <omp.h>
#endif

string const version = "0.1";

using namespace std;

// global random number thing
Ran myran;


void print_hello(char* name){
  //  printf("\n%s says: \"when I'm grown up, I will be a complete HMC simulation...\"\n\n",name);
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
  if ((*params).get_startcondition()) {
    printf("sourcefile = ");
    (*params).display_sourcefile();
    printf("\n");
  }
  else {
    printf("no sourcefile\n");
  }
  printf("**********************************************************\n");
  printf("\n");
  return;
}

#endif
