#include "gaugefield_hmc_tmp.h"


Opencl_Module_Hmc* Gaugefield_hmc::get_task_hmc(int dev)
{
	//if more than one device is used, the element dev from the array must be called here!!
	return (Opencl_Module_Hmc*)opencl_modules[task_hmc];
}

void Gaugefield_hmc::init_tasks()
{
	task_hmc = 0;

	opencl_modules = new Opencl_Module* [get_num_tasks()];

	//LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
	opencl_modules[task_hmc] = new Opencl_Module_Hmc[1];
	get_task_hmc(0)->init(queue[task_hmc], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_hmc), get_double_ext(task_hmc));

	return;
}

void Gaugefield_hmc::delete_variables()
{
	Gaugefield_hybrid::delete_variables();
	return;
}

void Gaugefield_hmc::finalize_opencl()
{

	Gaugefield_hybrid::finalize_opencl();
	return;
}

void Gaugefield_hmc::perform_hmc_step(hmc_observables *obs, int iter, hmc_float rnd_number, usetimer* solver_timer)
{

	return;
}

void Gaugefield_hmc::print_hmcobservables(hmc_observables obs, int iter, std::string filename)
{
	hmc_float exp_deltaH = exp(obs.deltaH);
	logger.trace() << "Observables: " << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq << "\t" << obs.poly.re << "\t" << obs.poly.im <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept ;
//  printf("Observables:%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%d\n",iter,obs.plaq,obs.tplaq,obs.splaq,obs.poly.re,obs.poly.im,obs.deltaH, exp_deltaH, obs.prob, obs.accept );
	std::fstream hmcout;
	hmcout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!hmcout.is_open()) throw File_Exception(filename);
	hmcout.width(8);
	hmcout << iter;
	hmcout << "\t";
	hmcout.precision(15);
	hmcout << obs.plaq << "\t" << obs.tplaq << "\t" << obs.splaq << "\t" << obs.poly.re << "\t" << obs.poly.im << "\t" << sqrt(obs.poly.re * obs.poly.re + obs.poly.im * obs.poly.im) <<  "\t" << obs.deltaH << "\t" << exp_deltaH << "\t" << obs.prob << "\t" << obs.accept << std::endl;
	hmcout.close();
	return;
}