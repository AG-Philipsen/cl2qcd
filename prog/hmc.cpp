#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "globals.h"
#include "hmcerrs.h"
#include "types.h"
#include "operations.h"
#include "geometry.h"
#include "testing.h"
#include "gaugeobservables.h"
#include "input.h"

using namespace std;

void print_hello(char* name){
  printf("\n%s says: \"when I'm grown up, I will be a complete HMC simulation...\"\n\n",name);
  return;
}

void print_info(inputparameters* params){
  printf("*********************************************\n");
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
  printf("*********************************************\n");
  printf("\n");

  return;
}

int main(int argc, char* argv[]) {

  char* progname = argv[0];

  print_hello(progname);

  char* inputfile = argv[1];

  //read input to variables

  inputparameters parameters;
  parameters.readfile(inputfile);

  print_info(&parameters);

  //do stuff with the gauge field
  hmc_gaugefield gaugefield;
  set_gaugefield_cold(&gaugefield);
  printf("plaquette: %f\n",plaquette(&gaugefield));
  printf("Polyakov loop: (%f,%f)\n",polyakov(&gaugefield).re,polyakov(&gaugefield).im);

  //init everything

  //perform updates

  //output stuff

  return HMC_SUCCESS;
}
