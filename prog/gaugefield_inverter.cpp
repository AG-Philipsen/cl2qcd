#include "gaugefield_inverter.h"

#include "meta/util.hpp"


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
	if(get_parameters().get_use_pointsource() == true)
		num_sources = 12;
	else
		num_sources = get_parameters().get_num_sources();

	size_t bufsize = num_sources * meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	logger.debug() << "allocate memory for solution-buffer on host of size " << bufsize / 1024. / 1024. / 1024. << " GByte";
	solution_buffer = new spinor [num_sources * meta::get_spinorfieldsize(get_parameters())];
	source_buffer = new spinor [num_sources * meta::get_spinorfieldsize(get_parameters())];

	task_solver = 0;
	task_correlator = 1;

	opencl_modules = new Opencl_Module* [get_num_tasks()];


	//LZ: right now, each task carries exactly one opencl device -> thus the below allocation with [1]. Could be generalized in future
	opencl_modules[task_solver] = new Opencl_Module_Fermions(get_parameters(), get_device_for_task(task_solver));
	get_task_solver()->init();

	opencl_modules[task_correlator] = new Opencl_Module_Correlator(get_parameters(), get_device_for_task(task_correlator));
	get_task_correlator()->init();

	clmem_corr = new hardware::buffers::Plain<spinor>(num_sources * meta::get_spinorfieldsize(get_parameters()), get_task_correlator()->get_device());
	clmem_source_solver = new hardware::buffers::Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), get_task_solver()->get_device());
	clmem_source_corr = new hardware::buffers::Plain<spinor>(meta::get_spinorfieldsize(get_parameters()), get_task_correlator()->get_device());
}

void Gaugefield_inverter::delete_variables()
{
	Gaugefield_hybrid::delete_variables();
}

void Gaugefield_inverter::finalize_opencl()
{
	delete clmem_source_solver;
	delete clmem_source_corr;
	delete clmem_corr;

	/// @todo this must be generalized if more than one device is used for one task
	for(int ntask = 0; ntask < get_num_tasks(); ntask++) {
		opencl_modules[ntask]->finalize();
	}
	Gaugefield_hybrid::finalize_opencl();

	logger.debug() << "free solution buffer";
	delete [] solution_buffer;
	logger.debug() << "free source buffer";
	delete [] source_buffer;
}

void Gaugefield_inverter::perform_inversion(usetimer* solver_timer)
{
	int use_eo = get_parameters().get_use_eo();

	//decide on type of sources
	int num_sources;
	if(get_parameters().get_use_pointsource() == true)
		num_sources = 12;
	else
		num_sources = get_parameters().get_num_sources();

	Opencl_Module_Fermions * solver = get_task_solver();

	hardware::buffers::Plain<spinor> clmem_res(meta::get_spinorfieldsize(get_parameters()), solver->get_device());

	//apply stout smearing if wanted
	if(get_parameters().get_use_smearing() == true) {
		solver->smear_gaugefield(solver->get_gaugefield(), std::vector<const hardware::buffers::SU3 *>());
	}

	//for CG, one needs a hermitian matrix...
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		logger.fatal() << "CG usage requires a hermitian matrix. This is not implemented yet...";
		//the call shoul be like this
		//::QplusQminus_eo f_eo(solver);
	}

	for(int k = 0; k < num_sources; k++) {
		//copy source from to device
		//NOTE: this is a blocking call!
		logger.debug() << "copy pointsource between devices";
		get_clmem_source_solver()->load(&source_buffer[k * meta::get_vol4d(get_parameters())]);
		logger.debug() << "calling solver..";
		if(use_eo) {
			::Aee f_eo(solver);
			solver->solver(f_eo, &clmem_res, get_clmem_source_solver(), solver->get_gaugefield(), solver_timer);
		} else {
			::M f_neo(solver);
			solver->solver(f_neo, &clmem_res, get_clmem_source_solver(), solver->get_gaugefield(), solver_timer);
		}

		//add solution to solution-buffer
		//NOTE: this is a blocking call!
		logger.debug() << "add solution...";
		clmem_res.dump(&solution_buffer[k * meta::get_vol4d(get_parameters())]);
	}

	if(get_parameters().get_use_smearing() == true) {
		solver->unsmear_gaugefield(solver->get_gaugefield());
	}
}

void Gaugefield_inverter::flavour_doublet_correlators(std::string corr_fn)
{
	using namespace std;

	//for now, make sure clmem_corr is properly filled; maybe later we can increase performance a bit by playing with this...
	get_clmem_corr()->load(solution_buffer);

	//suppose that the buffer on the device has been filled with the prior calculated solutions of the solver
	logger.debug() << "start calculating correlators...";

	ofstream of;
	of.open(corr_fn.c_str(), ios_base::app);
	if( !of.is_open() ) throw File_Exception(corr_fn);
	of << "# flavour doublet correlators" << endl;
	if(get_parameters().get_corr_dir() == 3) {
		of << "# format: J P z real complex"  << endl;
		of << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), z spatial distance, value (aggregate x y z)" << endl;
	} else {
		of << "# format: J P t real complex"  << endl;
		of << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), t timelike distance, value (aggregate x y z)" << endl;
	}


	int num_corr_entries =  0;
	switch (get_parameters().get_corr_dir()) {
		case 0 :
			num_corr_entries = get_parameters().get_ntime();
			break;
		case 3 :
			num_corr_entries = get_parameters().get_nspace();
			break;
		default :
			stringstream errmsg;
			errmsg << "Correlator direction " << get_parameters().get_corr_dir() << " has not been implemented.";
			throw Print_Error_Message(errmsg.str());
	}

	size_t buffersize = sizeof(hmc_float) * num_corr_entries;

	hmc_float* host_result_ps = new hmc_float [num_corr_entries];
	hmc_float* host_result_sc = new hmc_float [num_corr_entries];
	hmc_float* host_result_vx = new hmc_float [num_corr_entries];
	hmc_float* host_result_vy = new hmc_float [num_corr_entries];
	hmc_float* host_result_vz = new hmc_float [num_corr_entries];
	hmc_float* host_result_ax = new hmc_float [num_corr_entries];
	hmc_float* host_result_ay = new hmc_float [num_corr_entries];
	hmc_float* host_result_az = new hmc_float [num_corr_entries];

	cl_mem result_ps;
	cl_mem result_sc;
	cl_mem result_vx;
	cl_mem result_vy;
	cl_mem result_vz;
	cl_mem result_ax;
	cl_mem result_ay;
	cl_mem result_az;
	//LZ usually the correlator is calculated on CPU, then we don't need to copy all those buffers...
	bool needcopy = true;
	if( get_task_correlator()->get_device()->get_device_type() == CL_DEVICE_TYPE_CPU ) {
		needcopy = false;
		result_ps = get_task_correlator()->create_uhp_buffer(buffersize, host_result_ps);
		result_sc = get_task_correlator()->create_uhp_buffer(buffersize, host_result_sc);
		result_vx = get_task_correlator()->create_uhp_buffer(buffersize, host_result_vx);
		result_vy = get_task_correlator()->create_uhp_buffer(buffersize, host_result_vy);
		result_vz = get_task_correlator()->create_uhp_buffer(buffersize, host_result_vz);
		result_ax = get_task_correlator()->create_uhp_buffer(buffersize, host_result_ax);
		result_ay = get_task_correlator()->create_uhp_buffer(buffersize, host_result_ay);
		result_az = get_task_correlator()->create_uhp_buffer(buffersize, host_result_az);
	} else {
		result_ps = get_task_correlator()->create_rw_buffer(buffersize);
		result_sc = get_task_correlator()->create_rw_buffer(buffersize);
		result_vx = get_task_correlator()->create_rw_buffer(buffersize);
		result_vy = get_task_correlator()->create_rw_buffer(buffersize);
		result_vz = get_task_correlator()->create_rw_buffer(buffersize);
		result_ax = get_task_correlator()->create_rw_buffer(buffersize);
		result_ay = get_task_correlator()->create_rw_buffer(buffersize);
		result_az = get_task_correlator()->create_rw_buffer(buffersize);
	}

	logger.info() << "calculate correlators..." ;
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ps"), get_clmem_corr(), result_ps);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("sc"), get_clmem_corr(), result_sc);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vx"), get_clmem_corr(), result_vx);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vy"), get_clmem_corr(), result_vy);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vz"), get_clmem_corr(), result_vz);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ax"), get_clmem_corr(), result_ax);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ay"), get_clmem_corr(), result_ay);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("az"), get_clmem_corr(), result_az);

	//the correlator_device calcultions are all non-blocking, hence we need a barrier here: (why? -- MB)
	if( needcopy == false )
		get_device_for_task(task_correlator)->synchronize();

	//the pseudo-scalar (J=0, P=1)
	logger.info() << "pseudo scalar correlator:" ;
	if( needcopy == true )
		get_task_correlator()->get_buffer_from_device(result_ps, host_result_ps, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		logger.info() << j << "\t" << scientific << setprecision(14) << host_result_ps[j];
		of << scientific << setprecision(14) << "0 1\t" << j << "\t" << host_result_ps[j] << endl;
	}


	//the scalar (J=0, P=0)
	if( needcopy == true )
		get_task_correlator()->get_buffer_from_device(result_sc, host_result_sc, buffersize);
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "0 0\t" << j << "\t" << host_result_sc[j] << endl;
	}


	//the vector (J=1, P=1)
	if( needcopy == true ) {
		get_task_correlator()->get_buffer_from_device(result_vx, host_result_vx, buffersize);
		get_task_correlator()->get_buffer_from_device(result_vy, host_result_vy, buffersize);
		get_task_correlator()->get_buffer_from_device(result_vz, host_result_vz, buffersize);
	}
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 1\t" << j << "\t" << (host_result_vx[j] + host_result_vy[j] + host_result_vz[j]) / 3. << "\t" << host_result_vx[j] << "\t" << host_result_vy[j] << "\t" << host_result_vz[j] << endl;
	}


	//the axial vector (J=1, P=0)
	if( needcopy == true ) {
		get_task_correlator()->get_buffer_from_device(result_ax, host_result_ax, buffersize);
		get_task_correlator()->get_buffer_from_device(result_ay, host_result_ay, buffersize);
		get_task_correlator()->get_buffer_from_device(result_az, host_result_az, buffersize);
	}
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 0\t" << j << "\t" << (host_result_ax[j] + host_result_ay[j] + host_result_az[j]) / 3. << "\t" << host_result_ax[j] << "\t" << host_result_ay[j] << "\t" << host_result_az[j] << endl;
	}


	of << endl;
	of.close();
	delete [] host_result_ps;
	delete [] host_result_sc;
	delete [] host_result_vx;
	delete [] host_result_vy;
	delete [] host_result_vz;
	delete [] host_result_ax;
	delete [] host_result_ay;
	delete [] host_result_az;

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
}

void Gaugefield_inverter::create_sources()
{
	//create sources on the correlator-device and save them on the host
	if(get_parameters().get_use_pointsource() == true) {
		logger.debug() << "start creating point-sources...";
		for(int k = 0; k < 12; k++) {
			get_task_correlator()->create_point_source_device(get_clmem_source_corr(), k, get_source_pos_spatial(get_parameters()), get_parameters().get_pointsource_t());
			logger.debug() << "copy pointsource to host";
			get_clmem_source_corr()->dump(&source_buffer[k * meta::get_vol4d(get_parameters())]);
		}
	} else {
		logger.debug() << "start creating stochastic-sources...";
		int num_sources = get_parameters().get_num_sources();
		for(int k = 0; k < num_sources; k++) {
			get_task_correlator()->create_stochastic_source_device(get_clmem_source_corr());
			logger.debug() << "copy stochastic-source to host";
			get_clmem_source_corr()->dump(&source_buffer[k * meta::get_vol4d(get_parameters())]);
		}
	}
}

const hardware::buffers::Plain<spinor> * Gaugefield_inverter::get_clmem_corr()
{
	return clmem_corr;
}

const hardware::buffers::Plain<spinor> * Gaugefield_inverter::get_clmem_source_corr()
{
	return clmem_source_corr;
}

const hardware::buffers::Plain<spinor> * Gaugefield_inverter::get_clmem_source_solver()
{
	return clmem_source_solver;
}
