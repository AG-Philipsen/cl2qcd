#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <string>
#include <boost/lexical_cast.hpp>

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

int main(int argc, char* argv[]) {

  Timer timer;
  int time_measurements[10];
  for (int i = 0; i<10; i++){
    time_measurements[i] = 0;
  }

  char* progname = argv[0];
  print_hello(progname);

  char* inputfile = argv[1];
  inputparameters parameters;
  parameters.readfile(inputfile);

  sourcefileparameters parameters_source;
  hmc_gaugefield * gaugefield;
  gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialisation
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
    return -1;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Thermalization
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  printf("stuff..\n");
  
  timer.reset();
  printf("plaquette: %f\n",plaquette(gaugefield));
  time_measurements[2] = (int) timer.getTimeAndReset();
  
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  
  time_measurements[3] = (int) timer.getTimeAndReset();
  
  //Testing
  for (int i = 0; i<parameters.get_heatbathsteps()/2-1; i++){
    plaquette(gaugefield);
    time_measurements[2] += (int) timer.getTimeAndReset();
  
    polyakov(gaugefield);
    polyakov_x(gaugefield);
    polyakov_y(gaugefield);
    polyakov_z(gaugefield);
    
    time_measurements[3] += (int) timer.getTimeAndReset();
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Gaugefield-updates in thermalized system
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  printf("\nheatbath...\n");
  time_measurements[4] = 0;
  time_measurements[5] = 0;
  timer.reset();
  for (int i = 0; i<parameters.get_heatbathsteps(); i++){
      heatbath_update (gaugefield, parameters.get_beta());
      time_measurements[4] += (int) timer.getTimeAndReset();
//       heatbath_overrelax (gaugefield, parameters.get_beta());
  }
  
  for (int i = 0; i<parameters.get_heatbathsteps(); i++){
      heatbath_overrelax (gaugefield, parameters.get_beta());
      time_measurements[5] += (int) timer.getTimeAndReset();
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // More measurements
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  printf("now the stuff is..\n");  
  
  timer.reset();
  printf("plaquette: %f\n",plaquette(gaugefield));
  time_measurements[2] += (int) timer.getTimeAndReset();
  
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  
  time_measurements[3] += (int) timer.getTimeAndReset();
  
    //Testing
  for (int i = 0; i<parameters.get_heatbathsteps()/2-1; i++){
    plaquette(gaugefield);
    time_measurements[2] += (int) timer.getTimeAndReset();
  
    polyakov(gaugefield);
    polyakov_x(gaugefield);
    polyakov_y(gaugefield);
    polyakov_z(gaugefield);
    
    time_measurements[3] += (int) timer.getTimeAndReset();
  }  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  int totaltime;
  for (int i = 0; i<10; i++){
    totaltime += time_measurements[i];
  }
  //format is valid for measurements up to ca. 3 hours
  printf("\n");
  printf("**************************************************************\n");
  printf("total runtime:\t\t%i\n\n",  totaltime );
  printf("Times:\t\ttot\t\tavg\t\tsite\tperc\n");
  
  if (parameters.get_readsource()){
    printf("sourc:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[0], time_measurements[0], time_measurements[0]/VOL4D, (float)time_measurements[0]/totaltime*100 );
  }
  else printf("Init.:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[0], time_measurements[0], time_measurements[0]/VOL4D,(float)time_measurements[0]/totaltime*100 );
  if(!parameters.get_readsource())
    printf("therm:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[1], time_measurements[1]/parameters.get_thermalizationsteps() , time_measurements[1]/parameters.get_thermalizationsteps()/VOL4D , (float)time_measurements[1]/totaltime*100 );
  
  printf("plaq.:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[2], time_measurements[2]/parameters.get_heatbathsteps(), time_measurements[2]/parameters.get_heatbathsteps()/VOL4D, (float) time_measurements[2]*100/totaltime);
  printf("poly.:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[3], time_measurements[3]/parameters.get_heatbathsteps()/4 , time_measurements[3]/parameters.get_heatbathsteps()/4/VOL4D , (float) time_measurements[3]*100/totaltime);
  printf("updte:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[4] ,(time_measurements[4]/parameters.get_heatbathsteps() ) , (time_measurements[4]/parameters.get_heatbathsteps() )/VOL4D 
  	, (float) (time_measurements[4])*100/totaltime);
  printf("overr:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[5] ,(time_measurements[5]/parameters.get_heatbathsteps()) , (time_measurements[5]/parameters.get_heatbathsteps()) /VOL4D
  	, (float) (time_measurements[5])*100/totaltime);
  
  printf("\nuptot:\t%12i\t%12i\t%12i\t%.3f\n", time_measurements[4] + time_measurements[1],(time_measurements[4]/parameters.get_heatbathsteps() + time_measurements[1]/parameters.get_thermalizationsteps()) , (time_measurements[4]/parameters.get_heatbathsteps() + time_measurements[1]/parameters.get_thermalizationsteps())/VOL4D
  	, (float) (time_measurements[4] + time_measurements[1])*100/totaltime);
  printf("**************************************************************\n");
 
  //same some data to file
  ofstream out;
  string filename = ("time_");
  string space = ("_");
  filename += boost::lexical_cast<std::string>(NTIME);
  filename += space;
  filename += boost::lexical_cast<std::string>(NSPACE);
  filename += space;
#ifdef _RECONSTRUCT_TWELVE_
  filename += boost::lexical_cast<std::string>(1);
#else
  filename += boost::lexical_cast<std::string>(0);
#endif
  out.open(filename); 
  if (out.is_open())
  {
    out << NTIME << "\t" << NSPACE << "\t" << totaltime << "\t" << time_measurements[0] << "\t" <<  time_measurements[1]/parameters.get_thermalizationsteps() << "\t" << time_measurements[2]/parameters.get_heatbathsteps()<< "\t" << time_measurements[3]/parameters.get_heatbathsteps()/4 << "\t" << (time_measurements[4]/parameters.get_heatbathsteps() ) << "\t" << (time_measurements[5]/parameters.get_heatbathsteps())<< endl;
    out.close();
  }
  else cout << "Unable to open file for output" << endl;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // free variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  free (gaugefield);
  
  
  return HMC_SUCCESS;
}
