#include "hmc.h"

int main(int argc, char* argv[]) {

  char* progname = argv[0];
  print_hello(progname);

  opencl gpu(CL_DEVICE_TYPE_GPU);

  char* inputfile = argv[1];
  inputparameters parameters;
  parameters.readfile(inputfile);
  print_info(&parameters);

  usetimer timer(&parameters);
  timer.start();

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialization
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  sourcefileparameters parameters_source;
  hmc_gaugefield * gaugefield;
  gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));

  init_gaugefield(gaugefield,&parameters,&timer);


  //LZ: test opencl implementation -- in future, we will need timer measurements here...
  print_gaugeobservables(gaugefield);

  gpu.copy_gaugefield_to_device(gaugefield);
  gpu.run_heatbath(10,parameters.get_beta());
  gpu.get_gaugefield_from_device(gaugefield);

  print_gaugeobservables(gaugefield);

  return 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Thermalization
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  timer.reset();
  if(!parameters.get_startcondition()==START_FROM_SOURCE){
    cout << endl << "perform thermalization" << endl;
    for (int i = 0; i<parameters.get_thermalizationsteps(); i++){
      heatbath_update (gaugefield, parameters.get_beta());cout << i << " " << endl;
    }
    
  }
  timer.measurements(1);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Test-Measurements
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  printf("Initial measurements...\n");
  
  timer.reset();
  printf("plaquette: %f\n",plaquette(gaugefield));
  timer.getTimeAndReset(2);
  
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  
  timer.getTimeAndReset(3);
  
  //Testing
  for (int i = 0; i<parameters.get_heatbathsteps()/2-1; i++){
    plaquette(gaugefield);
    timer.add(2);
  
    polyakov(gaugefield);
    polyakov_x(gaugefield);
    polyakov_y(gaugefield);
    polyakov_z(gaugefield);
    
    timer.add(3);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Gaugefield-updates in thermalized system on CPU
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  printf("\nheatbath...\n");
  timer.zero(4);
  timer.zero(5);
  timer.reset();
  for (int i = 0; i<parameters.get_heatbathsteps(); i++){
      heatbath_update (gaugefield, parameters.get_beta());
      //      printf("%d %f %f\n",i,plaquette(gaugefield),polyakov(gaugefield).re);
      timer.add(4);
//        heatbath_overrelax (gaugefield, parameters.get_beta());
  }
  
  for (int i = 0; i<parameters.get_heatbathsteps(); i++){
      heatbath_overrelax (gaugefield, parameters.get_beta());
      timer.add(5);
  }
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // More measurements
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  printf("now the stuff is..\n");  
  
  timer.reset();
  printf("plaquette: %f\n",plaquette(gaugefield));
  timer.add(2);
  
  printf("Polyakov loop in t: (%f,%f)\n",polyakov(gaugefield).re,polyakov(gaugefield).im);
  printf("Polyakov loop in x: (%f,%f)\n",polyakov_x(gaugefield).re,polyakov_x(gaugefield).im);
  printf("Polyakov loop in y: (%f,%f)\n",polyakov_y(gaugefield).re,polyakov_y(gaugefield).im);
  printf("Polyakov loop in z: (%f,%f)\n",polyakov_z(gaugefield).re,polyakov_z(gaugefield).im);
  
  timer.add(3);
  
    //Testing
  for (int i = 0; i<parameters.get_heatbathsteps()/2-1; i++){
    plaquette(gaugefield);
    timer.add(2);
  
    polyakov(gaugefield);
    polyakov_x(gaugefield);
    polyakov_y(gaugefield);
    polyakov_z(gaugefield);
    
    timer.add(3);
  }  
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  timer.output();
 

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // free variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  free(gaugefield);

  return HMC_SUCCESS;
}
