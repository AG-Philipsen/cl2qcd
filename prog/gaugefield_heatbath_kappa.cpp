#include "gaugefield_heatbath_kappa.h"


hardware::code::Heatbath* Gaugefield_heatbath_kappa::get_task_heatbath()
{
	return (hardware::code::Heatbath*)opencl_modules[task_heatbath];
}

hardware::code::Kappa* Gaugefield_heatbath_kappa::get_task_kappa()
{
	return (hardware::code::Kappa*)opencl_modules[task_kappa];
}

void Gaugefield_heatbath_kappa::init_tasks()
{

	if (get_num_tasks() != 2)
		throw Print_Error_Message("We need 2 tasks");

	if (get_num_devices() != 2)
		logger.warn() << "Calculation of transport coefficients has not been tested for " << get_num_devices() ;

	opencl_modules = new hardware::code::Opencl_Module* [get_num_tasks()];

	//LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
	opencl_modules[task_kappa] = get_device_for_task(task_kappa)->get_kappa_code();

	opencl_modules[task_heatbath] = get_device_for_task(task_heatbath)->get_heatbath_code();
}

void Gaugefield_heatbath_kappa::perform_heatbath(int nheat, int nover)
{

	auto gf = get_device_for_task(task_heatbath)->get_gaugefield_code()->get_gaugefield();
	auto prng = &get_device_for_task(task_heatbath)->get_prng_code()->get_prng_buffer();

	for(int iter = 0; iter < nheat; iter++) {
		get_task_heatbath()->run_heatbath(gf, prng);
		for(int iter_over = 0; iter_over < nover; iter_over++)
			get_task_heatbath()->run_overrelax(gf, prng);
	}

	return;
}


void Gaugefield_heatbath_kappa::perform_tasks(int nheat, int nover)
{
	auto gf = get_device_for_task(task_heatbath)->get_gaugefield_code()->get_gaugefield();
	auto prng = &get_device_for_task(task_heatbath)->get_prng_code()->get_prng_buffer();

	for(int iter = 0; iter < nheat; iter++) {
		get_task_heatbath()->run_heatbath(gf, prng);
		for(int iter_over = 0; iter_over < nover; iter_over++)
			get_task_heatbath()->run_overrelax(gf, prng);
	}

	get_task_kappa()->run_kappa_clover(get_device_for_task(task_kappa)->get_gaugefield_code()->get_gaugefield(), get_parameters().get_beta());
}

void Gaugefield_heatbath_kappa::perform_tasks(int nheat, int nover, int* nheat_optimal)
{
	auto gf = get_device_for_task(task_heatbath)->get_gaugefield_code()->get_gaugefield();
	auto prng = &get_device_for_task(task_heatbath)->get_prng_code()->get_prng_buffer();

	uint64_t time_for_heatbath;
	uint64_t time_for_kappa;

	usetimer timer;
	timer.reset();
	for(int iter = 0; iter < nheat; iter++) {
		get_task_heatbath()->run_heatbath(gf, prng);
		for(int iter_over = 0; iter_over < nover; iter_over++)
			get_task_heatbath()->run_overrelax(gf, prng);
	}
	timer.add();
	time_for_heatbath = timer.getTime() / nheat;

	timer.reset();
	get_task_kappa()->run_kappa_clover(get_device_for_task(task_kappa)->get_gaugefield_code()->get_gaugefield(), get_parameters().get_beta());
	timer.add();
	time_for_kappa = timer.getTime();

	logger.info() << "Autotuning:";
	logger.info() << "\tAverage time for one heatbath iteration (including overrelaxation):  " << time_for_heatbath << " microseconds";
	logger.info() << "\tTime for calculation of transport coefficient:                       " << time_for_kappa    << " microseconds";

	*nheat_optimal = time_for_kappa / time_for_heatbath;
	if(*nheat_optimal <= 0) *nheat_optimal = 1;
	logger.info() << "\tThus set the optimal number of heatbath steps per tk calculation to: " << *nheat_optimal;

	return;
}

void Gaugefield_heatbath_kappa::delete_variables()
{
	Gaugefield_hybrid::delete_variables();

	return;
}

void Gaugefield_heatbath_kappa::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
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
