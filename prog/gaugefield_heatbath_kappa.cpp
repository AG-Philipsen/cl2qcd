#include "gaugefield_heatbath_kappa.h"


void Gaugefield_heatbath_kappa::init_devices(){

  opencl_modules = new Opencl_Module [get_num_tasks()];
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    opencl_modules[ntask].init(queue[ntask], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(ntask), get_double_ext(ntask));
  }

  return;
}

void Gaugefield_heatbath_kappa::delete_variables(){
  Gaugefield_hybrid::delete_variables();
  
  return;
}

void Gaugefield_heatbath_kappa::finalize_opencl(){

  Gaugefield_hybrid::finalize_opencl();

  return;
}
