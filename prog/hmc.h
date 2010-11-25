#ifndef _HMCH_
#define _HMCH_
//should only be included in main prog

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <boost/lexical_cast.hpp>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "operations.h"
#include "geometry.h"
#include "testing.h"
#include "gaugeobservables.h"
#include "input.h"
#include "readgauge.h"
#include "random.h"
#include "update.h"
#include "timer.h"

#include <CL/cl.hpp>
#include "opencl.h"

#ifdef _OPENMP
# include <omp.h>
#endif

using namespace std;

// global random number thing
Ran myran;

// opencl global vars
cl_device_type wanted_device=CL_DEVICE_TYPE_GPU; 
cl_context context;
cl_command_queue cmdqueue;
cl_program clprogram;
string cl_kernels_file = "opencl_kernels.cl";

void print_hello(char* name){
  printf("\n%s says: \"when I'm grown up, I will be a complete HMC simulation...\"\n\n",name);
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
  if ((*params).get_readsource()) {
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


void print_info_source(sourcefileparameters* params){
  printf("**********************************************************\n");
  printf("Sourcefile parameters: (list not complete)\n");
  printf("field:  %s\n",(*params).field_source);
  printf("LX:  \t%d\n",(*params).lx_source);
  printf("LY:  \t%d\n",(*params).ly_source);
  printf("LZ:  \t%d\n",(*params).lz_source);
  printf("LT:  \t%d\n",(*params).lt_source);
  printf("entries: %d\n", (*params).num_entries_source);
  printf("beta:  \t%f\n",(*params).beta_source);
  printf("mu:  \t%f\n",(*params).mu_source);
  printf("kappa:  %f\n",(*params).kappa_source);
  printf("mubar:  %f\n",(*params).mubar_source);
  printf("plaq: \t%f\n",(*params).plaquettevalue_source);
  printf("**********************************************************\n");
  printf("\n");
  return;
}


#endif
