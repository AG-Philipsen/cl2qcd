#include "gaugefield_inverter.h"


Opencl_Module_Fermions* Gaugefield_inverter::get_task_solver() {
  return (Opencl_Module_Fermions*)opencl_modules[task_solver];
}

Opencl_Module_Correlator* Gaugefield_inverter::get_task_correlator() {
  return (Opencl_Module_Correlator*)opencl_modules[task_correlator];
}

void Gaugefield_inverter::init_tasks(){

  solution_buffer = 0;

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
	logger.debug() << "free solution buffer";
	delete [] solution_buffer;
  return;
}

void Gaugefield_inverter::sync_solution_buffer(){
	size_t sfsize = 12*get_parameters()->get_spinorfieldsize()*sizeof(spinor);
	get_task_correlator()->copy_buffer_to_device(solution_buffer, get_task_correlator()->get_clmem_corr(), sfsize);
	return;
}

void Gaugefield_inverter::perform_inversion(usetimer* solver_timer){
	
	int use_eo = get_parameters()->get_use_eo();
	
	//decide on type of sources
	if(get_parameters()->get_use_pointsource() == true){
		//allocate host-memory for solution-buffer
		size_t sfsize = get_parameters()->get_spinorfieldsize()*sizeof(spinor);
		size_t bufsize = 12*sfsize;
		logger.debug() << "allocate memory for solution-buffer on host of size " << bufsize/1024./1024./1024. << " GByte";
		if(solution_buffer == 0) solution_buffer = new spinor [12*get_parameters()->get_spinorfieldsize()];
		spinor* sftmp = new spinor [get_parameters()->get_spinorfieldsize()];
		
		for(int k=0; k<12; k++) {
			//create source
			logger.debug() << "create pointsource";
			get_task_correlator()->create_point_source_device(get_task_correlator()->get_clmem_source(), k,get_parameters()->get_source_pos_spatial(),get_parameters()->get_source_pos_temporal());
			//copy source from one device to the other
			logger.debug() << "copy pointsource between devices";
			///@todo is this possible without the host in between?
			get_task_correlator()->get_buffer_from_device(get_task_correlator()->get_clmem_source(), sftmp, sfsize);
			get_task_solver()->copy_buffer_to_device(sftmp, get_task_solver()->get_clmem_source(), sfsize);
		
			logger.debug() << "calling solver..";
			if(use_eo == false)
				get_task_solver()->solver(M_call, get_task_solver()->get_clmem_inout(), get_task_solver()->get_clmem_source(), *get_clmem_gaugefield(), solver_timer, get_parameters()->get_cgmax());
			else
				get_task_solver()->solver(Aee_call, get_task_solver()->get_clmem_inout(), get_task_solver()->get_clmem_source(), *get_clmem_gaugefield(), solver_timer, get_parameters()->get_cgmax());

			//add solution to solution-buffer
			//NOTE: this is a blocking call!
			logger.debug() << "add solution...";
			get_task_solver()->get_buffer_from_device(get_task_solver()->get_clmem_inout(), &solution_buffer[k*VOL4D], sfsize);
		}
		delete [] sftmp;
	}
	else{
		//use stochastic sources
		//allocate host-memory for solution-buffer
		int sfsize = get_parameters()->get_spinorfieldsize()*sizeof(spinor);
		size_t bufsize = sfsize*get_parameters()->get_num_sources();
		logger.debug() << "allocate memory for solution-buffer on host of size " << bufsize/1024./1024./1024. << " GByte";
		solution_buffer = (spinorfield*) malloc(bufsize);
		spinorfield* sftmp = (spinorfield*) malloc(bufsize);
		
		throw Print_Error_Message("Stochastic Sources not yet implemented.");
		for(int k=0; k<get_parameters()->get_num_sources(); k++) {
			///@todo this is the same as above, just with a different source-generating-function
		}
	}
	
  return;
}

void Gaugefield_inverter::flavour_doublet_correlators(string corr_fn){
	//suppose that the buffer on the device has been filled with the prior calculated solutions of the solver
	logger.debug() << "start calculating correlators...";
	get_task_correlator()->ps_correlator_device(get_task_correlator()->get_clmem_corr());

  return;
}
