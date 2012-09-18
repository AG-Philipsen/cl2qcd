#include "gaugefield_heatbath.h"


Opencl_Module_Heatbath* Gaugefield_heatbath::get_task_heatbath()
{
	return (Opencl_Module_Heatbath*)opencl_modules[task_heatbath];
}

void Gaugefield_heatbath::init_tasks()
{

	//CP: At the moment, this is meant to be only for 1 task
	if(get_num_tasks() != 1) throw Print_Error_Message("We need exactly 1 task");
	task_heatbath = 0;

	//this must be a pointer in general, since most likely not every entry in the array is of the same size
	opencl_modules = new Opencl_Module* [get_num_tasks()];

	opencl_modules[task_heatbath] = new Opencl_Module_Heatbath(get_parameters());
	get_task_heatbath()->init(queue[task_heatbath], get_max_compute_units(task_heatbath), get_double_ext(task_heatbath), task_heatbath);

	return;
}

void Gaugefield_heatbath::perform_tasks(int nover)
{
	get_task_heatbath()->run_heatbath();
	for(int iter_over = 0; iter_over < nover; iter_over++)
		get_task_heatbath()->run_overrelax();
	return;
}

void Gaugefield_heatbath::delete_variables()
{
	Gaugefield_hybrid::delete_variables();
	return;
}

void Gaugefield_heatbath::finalize_opencl()
{
	/// @todo this must be generalized if more than one device is used for one task
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		opencl_modules[ntask]->finalize();
	}
	Gaugefield_hybrid::finalize_opencl();
	return;
}
