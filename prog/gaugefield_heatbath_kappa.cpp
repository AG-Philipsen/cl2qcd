#include "gaugefield_heatbath_kappa.h"


Opencl_Module_Heatbath* Gaugefield_heatbath_kappa::get_task_heatbath() {
  return (Opencl_Module_Heatbath*)opencl_modules[task_heatbath];
}

Opencl_Module_Kappa* Gaugefield_heatbath_kappa::get_task_kappa() {
  return (Opencl_Module_Kappa*)opencl_modules[task_kappa];
}

void Gaugefield_heatbath_kappa::init_tasks(){

  switch (get_num_tasks()) {
  case 2 :
    task_heatbath = 0;
    task_kappa = 1;
    break;
  case 1:
    task_heatbath = 0;
    task_kappa = 0;
    break;
  default:
    throw Print_Error_Message("We need 1 or 2 tasks");
  }

  opencl_modules = new Opencl_Module* [get_num_tasks()];

  //LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
  opencl_modules[task_kappa] = new Opencl_Module_Kappa[1];
  get_task_kappa()->init(queue[task_kappa], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_kappa), get_double_ext(task_kappa));

  opencl_modules[task_heatbath] = new Opencl_Module_Heatbath[1];
  get_task_heatbath()->init(queue[task_heatbath], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_heatbath), get_double_ext(task_heatbath));

  return;
}

void Gaugefield_heatbath_kappa::perform_heatbath(int nheat, int nover){

  for(int iter = 0; iter < nheat; iter++) {
    get_task_heatbath()->run_heatbath();
    for(int iter_over = 0; iter_over < nover; iter_over++) 
			get_task_heatbath()->run_overrelax();
  }

  return;
}


void Gaugefield_heatbath_kappa::perform_tasks(int nheat, int nover){

  for(int iter = 0; iter < nheat; iter++) {
    get_task_heatbath()->run_heatbath();
    for(int iter_over = 0; iter_over < nover; iter_over++) 
			get_task_heatbath()->run_overrelax();
  }

  get_task_kappa()->run_kappa_clover(get_parameters()->get_beta());

  return;
}

void Gaugefield_heatbath_kappa::perform_tasks(int nheat, int nover, int* nheat_optimal){

  uint64_t time_for_heatbath;
  uint64_t time_for_kappa;

  usetimer timer;
  timer.reset();
  for(int iter = 0; iter < nheat; iter++) {
    get_task_heatbath()->run_heatbath();
      for(int iter_over = 0; iter_over < nover; iter_over++) 
	get_task_heatbath()->run_overrelax();
  }
  timer.add();
  time_for_heatbath = timer.getTime() / nheat;

  timer.reset();
  get_task_kappa()->run_kappa_clover(get_parameters()->get_beta());
  timer.add();
  time_for_kappa = timer.getTime();

  logger.info()<<"Autotuning:";
  logger.info()<<"\tAverage time for one heatbath iteration (including overrelaxation):  " << time_for_heatbath << " microseconds";
  logger.info()<<"\tTime for calculation of transport coefficient:                       " << time_for_kappa    << " microseconds";

  *nheat_optimal = time_for_kappa / time_for_heatbath;
  if(*nheat_optimal <= 0) *nheat_optimal = 1;
  logger.info()<<"\tThus set the optimal number of heatbath steps per tk calculation to: " << *nheat_optimal;

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


void Gaugefield_heatbath_kappa::print_kappa(int iter)
{
  hmc_float kappa_clover =   get_task_kappa()->get_kappa_clover();
  logger.info() << iter << '\t' << kappa_clover;
  return;
}

void Gaugefield_heatbath_kappa::print_kappa(int iter, std::string filename)
{

  hmc_float kappa_clover =   get_task_kappa()->get_kappa_clover();

	std::fstream outfile;
	outfile.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!outfile.is_open()) throw File_Exception(filename);
	outfile.width(8);
	outfile << iter;
	outfile << "\t";
	outfile.precision(15);
	outfile << kappa_clover << std::endl;
	outfile.close();
	return;
}
