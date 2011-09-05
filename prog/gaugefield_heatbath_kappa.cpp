#include "gaugefield_heatbath_kappa.h"


Opencl_Module_Heatbath* Gaugefield_heatbath_kappa::get_devices_heatbath(int task) {
  return (Opencl_Module_Heatbath*)opencl_modules[task];
}

Opencl_Module_Kappa* Gaugefield_heatbath_kappa::get_devices_kappa(int task) {
  return (Opencl_Module_Kappa*)opencl_modules[task];
}

void Gaugefield_heatbath_kappa::init_devices(){

  opencl_modules = new Opencl_Module* [get_num_tasks()];

  opencl_modules[0] = new Opencl_Module_Heatbath[1];
  get_devices_heatbath(0)->init(queue[0], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(0), get_double_ext(0));
    
  opencl_modules[1] = new Opencl_Module_Kappa[1];
  get_devices_kappa(1)->init(queue[1], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(1), get_double_ext(1));
  


  return;
}

void Gaugefield_heatbath_kappa::perform_tasks(int nheat, int nover){

  for(int iter = 0; iter < nheat; iter++) {
    get_devices_heatbath(0)->run_heatbath();
      for(int iter_over = 0; iter_over < nover; iter_over++) 
	get_devices_heatbath(0)->run_overrelax();
  }

  get_devices_kappa(1)->run_kappa_clover(get_parameters()->get_beta());

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
