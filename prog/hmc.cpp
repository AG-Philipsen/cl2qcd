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

using namespace std;

void print_hello(char* name){
  printf("\n%s says: \"when I'm grown up, I will be a complete HMC simulation...\"\n\n",name);
  return;
}

void print_info(){
  printf("NSPACE:  %d\n",NSPACE);
  printf("NTIME:   %d\n",NTIME);
  printf("NDIM:    %d\n",NDIM);
  printf("NCOLOR:  %d\n",NC);
  printf("NSPIN:   %d\n",NSPIN);
  printf("+ input parameters read...\n");
  return;
}

int main(int argc, char* argv[]) {

  char* progname = argv[0];

  print_hello(progname);

  //read input to variables

  print_info();

  hmc_gaugefield gaugefield;
  set_gaugefield_cold(&gaugefield);
  printf("plaquette: %f\n",plaquette(&gaugefield));
  printf("Polyakov loop: (%f,%f)\n",polyakov(&gaugefield).re,polyakov(&gaugefield).im);

  //init everything

  //perform updates

  //output stuff

  return HMC_SUCCESS;
}
