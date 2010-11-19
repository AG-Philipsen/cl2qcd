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
#include "readgauge.h"
#include "random.h"
#include "update.h"
#include "timer.h"

// global random number thing
Ran myran;

using namespace std;

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
    printf("no sourcefile");
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

int main(int argc, char* argv[]) {

  Timer timer;
  int time_measurements[10];
  char* progname = argv[0];
  print_hello(progname);

  char* inputfile = argv[1];
  inputparameters parameters;
  parameters.readfile(inputfile);

  sourcefileparameters parameters_source;
  hmc_gaugefield * gaugefield;
  gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialitation
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int err;
  if(parameters.get_readsource()){
    timer.reset();
    //tmp gauge field
    hmc_float * gaugefield_tmp;
    gaugefield_tmp = (hmc_float*) malloc(sizeof(hmc_float)*NDIM*NC*NC*NTIME*VOLSPACE);
    err = parameters_source.readsourcefile(&(parameters.sourcefile)[0], parameters.get_prec(), &gaugefield_tmp);
    err = set_gaugefield_source(gaugefield, gaugefield_tmp, parameters_source.num_entries_source);
    free(gaugefield_tmp);
    time_measurements[0] = (int) timer.getTime();
  }
  else{
    timer.reset();
    set_gaugefield_cold(gaugefield);
    time_measurements[0] = (int) timer.getTime();
  }
  
  print_info(&parameters);
  if (parameters.get_readsource() && err == 0){
    print_info_source(&parameters_source);
  }
  else if (parameters.get_readsource() && err != 0){
    printf("error in setting vals from source!!! check global settings!!!\n\n");
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Thermalization
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(!parameters.get_readsource()){
    timer.reset();
    cout << endl << "perform thermalization" << endl;
    for (int i = 0; i<parameters.get_thermalizationsteps(); i++){
      heatbath_update (gaugefield, parameters.get_beta());
    }
    time_measurements[1] = (int) timer.getTime();
  }
  else
    time_measurements[1] = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Test-Measurements
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  printf("stuff..\n");
  
  timer.reset();
  printf("plaquette: %f\n",plaquette(gaugefield));
  time_measurements[2] = (int) timer.getTimeAndReset();
  
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",spatial_polyakov(gaugefield,1).re,spatial_polyakov(gaugefield,1).im);
  printf("Polyakov loop in y: (%f,%f)\n",spatial_polyakov(gaugefield,2).re,spatial_polyakov(gaugefield,2).im);
  printf("Polyakov loop in z: (%f,%f)\n",spatial_polyakov(gaugefield,3).re,spatial_polyakov(gaugefield,3).im);
  
  time_measurements[3] = (int) timer.getTimeAndReset();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Gaugefield-updates
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  printf("\nheatbath...\n");
  timer.reset();
  for (int i = 0; i<parameters.get_heatbathsteps(); i++){
      heatbath_update (gaugefield, parameters.get_beta());
      heatbath_overrelax (gaugefield, parameters.get_beta());
  }
  time_measurements[4] = (int) timer.getTimeAndReset();
  printf("now the stuff is..\n");
  printf("plaquette: %f\n",plaquette(gaugefield));
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Output
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  cout << endl << "Time (mus): " << endl << endl;
  if (parameters.get_readsource()){
    cout  << "reading sourcefile:\t" << time_measurements[0] << endl;
  }
  else cout  << "Initialitation:\t" << time_measurements[0] << endl;
  if(!parameters.get_readsource())
    cout  << "thermalization:\t" << time_measurements[1] << endl;
  cout  << "average plaquette:\t" << time_measurements[2] << endl;
  cout  << "average polyakovloop:\t" << time_measurements[3]/7 << endl;
  cout  << "average updating took:\t" << time_measurements[4]/parameters.get_heatbathsteps()/2 + time_measurements[1]/parameters.get_thermalizationsteps() << endl << endl;
  

  return HMC_SUCCESS;
}
