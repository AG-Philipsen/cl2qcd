#include "hmc.h"

int main(int argc, char* argv[]) {

  //   init_opencl();
  // return 0;

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
  // Initialization
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
  timer.reset();
  if(!parameters.get_readsource()){
    cout << endl << "perform thermalization" << endl;
    for (int i = 0; i<parameters.get_thermalizationsteps(); i++){
      heatbath_update (gaugefield, parameters.get_beta());cout << i << " " << endl;
    }
    
  }
  time_measurements[1] = (int) timer.getTime();

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
//        heatbath_overrelax (gaugefield, parameters.get_beta());
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
 
  //save some data to file
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
