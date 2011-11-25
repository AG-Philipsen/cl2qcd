#include "../gaugefield_inverter.h"


Opencl_Module_Fermions* Gaugefield_inverter::get_task_solver()
{
	return (Opencl_Module_Fermions*)opencl_modules[task_solver];
}

Opencl_Module_Correlator* Gaugefield_inverter::get_task_correlator()
{
	return (Opencl_Module_Correlator*)opencl_modules[task_correlator];
}

void Gaugefield_inverter::init_tasks()
{
  //this is changed compared to "gaugefield_inverter.cpp" because not all kernels are needed

	task_solver = 0;
	task_correlator = 1;

	opencl_modules = new Opencl_Module* [get_num_tasks()];

	//LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
	opencl_modules[task_solver] = new Opencl_Module_Fermions[1];
	get_task_solver()->init(queue[task_solver], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_solver), get_double_ext(task_solver));
	/*
	opencl_modules[task_correlator] = new Opencl_Module_Correlator[1];
	get_task_correlator()->init(queue[task_correlator], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_correlator), get_double_ext(task_correlator));
	*/


}

void Gaugefield_inverter::delete_variables()
{
  //Gaugefield_hybrid::delete_variables();
}

void Gaugefield_inverter::finalize_opencl()
{
  // task_correlator is not finalized because it is not initialized!!
  for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
    if (ntask == task_correlator) continue;
    opencl_modules[ntask]->finalize();
  }

  Gaugefield_hybrid::finalize_opencl();
}

