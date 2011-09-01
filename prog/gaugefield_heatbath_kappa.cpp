#include "gaugefield_heatbath_kappa.h"


void Gaugefield_heatbath_kappa::init_devices(){

  init_random_arrays();

  opencl_modules = new Opencl_Module [get_num_tasks()];
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    opencl_modules[ntask].init(queue[ntask], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(ntask), get_double_ext(ntask));
  }

  return;
}

void Gaugefield_heatbath_kappa::init_random_arrays(){
  // Prepare random number arrays, for each task and device separately
  numrndstates    = new int          [get_num_tasks()];
  sizeof_rndarray = new size_t       [get_num_tasks()];
  rndarray        = new hmc_ocl_ran* [get_num_tasks()];
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    if(get_device_type(ntask) == CL_DEVICE_TYPE_GPU)
      numrndstates[ntask] = 5120;
    else
      numrndstates[ntask] = 64;
    rndarray[ntask] = new hmc_ocl_ran [numrndstates[ntask]];
    sizeof_rndarray[ntask] = sizeof(hmc_ocl_ran)*numrndstates[ntask];
    init_random_seeds(rndarray[ntask], "rand_seeds", numrndstates[ntask]);
  }

  return;
}

void Gaugefield_heatbath_kappa::delete_variables(){

  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
      delete [] rndarray[ntask];
  }
  delete [] numrndstates;
  delete [] sizeof_rndarray;
  delete [] rndarray;

  Gaugefield_hybrid::delete_variables();
  
  return;
}

void Gaugefield_heatbath_kappa::finalize_opencl(){

  Gaugefield_hybrid::finalize_opencl();

  return;
}

hmc_ocl_ran* Gaugefield_heatbath_kappa::get_rndarray(int ntask){
  if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("rndarray index out of range",__FILE__,__LINE__); 
  return rndarray[ntask];
}
