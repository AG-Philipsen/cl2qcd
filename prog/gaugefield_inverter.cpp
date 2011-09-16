#include "gaugefield_inverter.h"


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

	//allocate host-memory for solution- and source-buffer
	int num_sources;
	if(get_parameters()->get_use_pointsource() == true)
		num_sources = 12;
	else
		num_sources = get_parameters()->get_num_sources();

	size_t bufsize = num_sources * get_parameters()->get_spinorfieldsize() * sizeof(spinor);
	logger.debug() << "allocate memory for solution-buffer on host of size " << bufsize / 1024. / 1024. / 1024. << " GByte";
	solution_buffer = new spinor [num_sources*get_parameters()->get_spinorfieldsize()];
	source_buffer = new spinor [num_sources*get_parameters()->get_spinorfieldsize()];

	task_solver = 0;
	task_correlator = 1;

	opencl_modules = new Opencl_Module* [get_num_tasks()];


	//LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
	opencl_modules[task_solver] = new Opencl_Module_Fermions[1];
	get_task_solver()->init(queue[task_solver], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_solver), get_double_ext(task_solver));

	opencl_modules[task_correlator] = new Opencl_Module_Correlator[1];
	get_task_correlator()->init(queue[task_correlator], get_clmem_gaugefield(), get_parameters(), get_max_compute_units(task_correlator), get_double_ext(task_correlator));

	return;
}

void Gaugefield_inverter::delete_variables()
{
	Gaugefield_hybrid::delete_variables();
	return;
}

void Gaugefield_inverter::finalize_opencl()
{

	Gaugefield_hybrid::finalize_opencl();
	logger.debug() << "free solution buffer";
	delete [] solution_buffer;
	logger.debug() << "free source buffer";
	delete [] source_buffer;
	return;
}

void Gaugefield_inverter::sync_solution_buffer()
{
	size_t sfsize = 12 * get_parameters()->get_spinorfieldsize() * sizeof(spinor);
	get_task_correlator()->copy_buffer_to_device(solution_buffer, get_task_correlator()->get_clmem_corr(), sfsize);
	return;
}

void Gaugefield_inverter::perform_inversion(usetimer* solver_timer)
{

	int use_eo = get_parameters()->get_use_eo();

	//decide on type of sources
	int num_sources;
	if(get_parameters()->get_use_pointsource() == true)
		num_sources = 12;
	else
		num_sources = get_parameters()->get_num_sources();

	//allocate host-memory for tmp-buffer
	size_t sfsize = get_parameters()->get_spinorfieldsize() * sizeof(spinor);
	spinor* sftmp = new spinor [get_parameters()->get_spinorfieldsize()];

	for(int k = 0; k < num_sources; k++) {
		//copy source from to device
		//NOTE: this is a blocking call!
		logger.debug() << "copy pointsource between devices";
		get_task_solver()->copy_buffer_to_device(&source_buffer[k*VOL4D], get_task_solver()->get_clmem_source(), sfsize);

		logger.debug() << "calling solver..";
		if(use_eo == false)
			get_task_solver()->solver(M_call, get_task_solver()->get_clmem_inout(), get_task_solver()->get_clmem_source(), *get_clmem_gaugefield(), solver_timer);
		else
			get_task_solver()->solver(Aee_call, get_task_solver()->get_clmem_inout(), get_task_solver()->get_clmem_source(), *get_clmem_gaugefield(), solver_timer);

		//add solution to solution-buffer
		//NOTE: this is a blocking call!
		logger.debug() << "add solution...";
		get_task_solver()->get_buffer_from_device(get_task_solver()->get_clmem_inout(), &solution_buffer[k*VOL4D], sfsize);
	}
	delete [] sftmp;

	return;
}

void Gaugefield_inverter::flavour_doublet_correlators(string corr_fn)
{
	//suppose that the buffer on the device has been filled with the prior calculated solutions of the solver
	logger.debug() << "start calculating correlators...";

	ofstream of;
	of.open(corr_fn.c_str(), ios_base::app);
	if( !of.is_open() ) throw File_Exception(corr_fn);
	of << "# flavour doublet correlators" << endl;
	if(get_parameters()->get_corr_dir() == 3) {
		of << "# format: J P z real complex"  << endl;
		of << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), z spatial distance, value (aggregate x y z)" << endl;
	} else {
		of << "# format: J P t real complex"  << endl;
		of << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), t timelike distance, value (aggregate x y z)" << endl;
	}


	int num_corr_entries =  0;
	switch (get_parameters()->get_corr_dir()) {
		case 0 :
			num_corr_entries = get_parameters()->get_nt();
			break;
		case 3 :
			num_corr_entries = get_parameters()->get_ns();
			break;
		default :
			stringstream errmsg;
			errmsg << "Correlator direction " << get_parameters()->get_corr_dir() << " has not been implemented.";
			throw Print_Error_Message(errmsg.str());
	}

	size_t buffersize = sizeof(hmc_float) * num_corr_entries;

	cl_mem result_ps = get_task_correlator()->create_rw_buffer(buffersize);
	cl_mem result_sc = get_task_correlator()->create_rw_buffer(buffersize);
	cl_mem result_vx = get_task_correlator()->create_rw_buffer(buffersize);
	cl_mem result_vy = get_task_correlator()->create_rw_buffer(buffersize);
	cl_mem result_vz = get_task_correlator()->create_rw_buffer(buffersize);
	cl_mem result_ax = get_task_correlator()->create_rw_buffer(buffersize);
	cl_mem result_ay = get_task_correlator()->create_rw_buffer(buffersize);
	cl_mem result_az = get_task_correlator()->create_rw_buffer(buffersize);

	hmc_float* host_result = new hmc_float [num_corr_entries];
	hmc_float* host_result_y = new hmc_float [num_corr_entries];
	hmc_float* host_result_z = new hmc_float [num_corr_entries];

	logger.info() << "calculate correlators..." ;
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ps"), get_task_correlator()->get_clmem_corr(), result_ps);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("sc"), get_task_correlator()->get_clmem_corr(), result_sc);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vx"), get_task_correlator()->get_clmem_corr(), result_vx);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vy"), get_task_correlator()->get_clmem_corr(), result_vy);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vz"), get_task_correlator()->get_clmem_corr(), result_vz);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ax"), get_task_correlator()->get_clmem_corr(), result_ax);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ay"), get_task_correlator()->get_clmem_corr(), result_ay);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("az"), get_task_correlator()->get_clmem_corr(), result_az);

	//the pseudo-scalar (J=0, P=1)
	logger.info() << "pseudo scalar correlator:" ;
	get_task_correlator()->get_buffer_from_device(result_ps, host_result, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		logger.info() << j << "\t" << scientific << setprecision(14) << host_result[j];
		of << scientific << setprecision(14) << "0 1\t" << j << "\t" << host_result[j] << endl;
	}


	//the scalar (J=0, P=0)
	get_task_correlator()->get_buffer_from_device(result_sc, host_result, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "0 0\t" << j << "\t" << host_result[j] << endl;
	}


	//the vector (J=1, P=1)
	get_task_correlator()->get_buffer_from_device(result_vx, host_result, buffersize);
	get_task_correlator()->get_buffer_from_device(result_vy, host_result_y, buffersize);
	get_task_correlator()->get_buffer_from_device(result_vz, host_result_z, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 1\t" << j << "\t" << (host_result[j] + host_result_y[j] + host_result_z[j]) / 3. << "\t" << host_result[j] << "\t" << host_result_y[j] << "\t" << host_result_z[j] << endl;
	}


	//the axial vector (J=1, P=0)
	get_task_correlator()->get_buffer_from_device(result_ax, host_result, buffersize);
	get_task_correlator()->get_buffer_from_device(result_ay, host_result_y, buffersize);
	get_task_correlator()->get_buffer_from_device(result_az, host_result_z, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 0\t" << j << "\t" << (host_result[j] + host_result_y[j] + host_result_z[j]) / 3. << "\t" << host_result[j] << "\t" << host_result_y[j] << "\t" << host_result_z[j] << endl;
	}


	of << endl;
	of.close();
	delete [] host_result;
	delete [] host_result_y;
	delete [] host_result_z;
	cl_int clerr = CL_SUCCESS;
	clerr = clReleaseMemObject(result_ps);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_sc);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_vx);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_vy);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_vz);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_ax);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_ay);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	clerr = clReleaseMemObject(result_az);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseMemObject", __FILE__, __LINE__);

	return;
}

void Gaugefield_inverter::create_sources()
{
	//create sources on the correlator-device and save them on the host
	size_t sfsize = get_parameters()->get_spinorfieldsize() * sizeof(spinor);
	if(get_parameters()->get_use_pointsource() == true) {
		logger.debug() << "start creating point-sources...";
		for(int k = 0; k < 12; k++) {
			get_task_correlator()->create_point_source_device(get_task_correlator()->get_clmem_source(), k, get_parameters()->get_source_pos_spatial(), get_parameters()->get_source_pos_temporal());
			logger.debug() << "copy pointsource to host";
			get_task_correlator()->get_buffer_from_device(get_task_correlator()->get_clmem_source(), &source_buffer[k*VOL4D], sfsize);
		}
	} else {
		logger.debug() << "start creating stochastic-sources...";
		int num_sources = get_parameters()->get_num_sources();
		for(int k = 0; k < num_sources; k++) {
			get_task_correlator()->create_stochastic_source_device(get_task_correlator()->get_clmem_source());
			logger.debug() << "copy stochastic-source to host";
			get_task_correlator()->get_buffer_from_device(get_task_correlator()->get_clmem_source(), &source_buffer[k*VOL4D], sfsize);
		}
	}
	return;
}
