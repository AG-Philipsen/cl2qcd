#include "gaugefield_inverter.h"


Opencl_Module_Fermions* Gaugefield_inverter::get_task_solver() {
  return (Opencl_Module_Fermions*)opencl_modules[task_solver];
}

Opencl_Module_Correlator* Gaugefield_inverter::get_task_correlator() {
  return (Opencl_Module_Correlator*)opencl_modules[task_correlator];
}

void Gaugefield_inverter::init_tasks(){

  switch (get_num_tasks()) {
  case 2 :
    task_solver = 0;
    task_correlator = 1;
    break;
  case 1:
    task_solver = 0;
    task_correlator = 0;
    break;
  default:
    throw Print_Error_Message("We need exactly 2 tasks");
  }

  opencl_modules = new Opencl_Module* [get_num_tasks()];

  //LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
  opencl_modules[task_solver] = new Opencl_Module_Fermions[1];
  get_task_solver()->init(queue[task_solver], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_solver), get_double_ext(task_solver));

  opencl_modules[task_correlator] = new Opencl_Module_Correlator[1];
  get_task_correlator()->init(queue[task_correlator], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_correlator), get_double_ext(task_correlator));

  return;
}

void Gaugefield_inverter::delete_variables(){
  Gaugefield_hybrid::delete_variables();  
  return;
}

void Gaugefield_inverter::finalize_opencl(){

  Gaugefield_hybrid::finalize_opencl();
  return;
}


void Gaugefield_inverter::perform_inversion(usetimer* solver_timer){
  (*solver_timer).reset();

  cout<<"buh"<<endl;

  (*solver_timer).add();
  return;
}

void Gaugefield_inverter::flavour_doublet_correlators(string corr_fn){

  cout<<"buh"<<endl;

  return;
}
