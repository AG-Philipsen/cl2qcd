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
  int num_sources = get_parameters().get_num_sources();

	size_t bufsize = num_sources * meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	logger.debug() << "allocate memory for solution-buffer on host of size " << bufsize / 1024. / 1024. / 1024. << " GByte";
	solution_buffer = new spinor [num_sources * meta::get_spinorfieldsize(get_parameters())];
	source_buffer = new spinor [num_sources * meta::get_spinorfieldsize(get_parameters())];

	task_solver = 0;
	task_correlator = 1;

	opencl_modules = new Opencl_Module* [get_num_tasks()];

	opencl_modules[task_solver] = new Opencl_Module_Fermions(get_parameters(), get_device_for_task(task_solver));
	opencl_modules[task_correlator] = new Opencl_Module_Correlator(get_parameters(), get_device_for_task(task_correlator));

	clmem_corr = new hardware::buffers::Plain<spinor>(get_parameters().get_num_sources() * meta::get_spinorfieldsize(get_parameters()), get_task_correlator()->get_device());
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

	Gaugefield_hybrid::finalize_opencl();

	logger.debug() << "free solution buffer";
	delete [] solution_buffer;
	logger.debug() << "free source buffer";
	delete [] source_buffer;
}

void Gaugefield_inverter::invert_M_nf2_upperflavour(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::Plain<spinor> * source, const hardware::buffers::SU3 * gf, usetimer * solvertimer)
{
	/** This solves the sparse-matrix system
	 *  A x = b
	 *  with  x == inout
	 *        A == f
	 *        b == source
	 * using a Krylov-Solver (BiCGStab or CG)
	 */
	int converged = -1;
	Opencl_Module_Fermions * solver = get_task_solver();
	if(get_parameters().get_profile_solver() ) (*solvertimer).reset();

	if( !get_parameters().get_use_eo() ){
	  //noneo case
	  //Trial solution
	  ///@todo this should go into a more general function
	  solver->set_spinorfield_cold_device(inout);
	  if(get_parameters().get_solver() == meta::Inputparameters::cg) {
	    //to use cg, one needs an hermitian matrix, which is QplusQminus
	    //the source must now be gamma5 b, to obtain the desired solution in the end
	    solver->gamma5_device(source);
	    ::QplusQminus f_neo(solver);
	    converged = solver->cg(f_neo, inout, source, gf, get_parameters().get_solver_prec());
	    //now, calc Qminus inout to obtain x = A^⁻1 b
	    //therefore, use source as an intermediate buffer
	    solver->Qminus(inout, source, gf, get_parameters().get_kappa(), meta::get_mubar(get_parameters() ));
	    //save the result to inout
	    hardware::buffers::copyData(inout, source);
	  } else {
	    ::M f_neo(solver);
	    converged = solver->bicgstab(f_neo, inout, source, gf, get_parameters().get_solver_prec());
	  }
	}
	else{
	  /**
	   * If even-odd-preconditioning is used, the inversion is split up
	   * into even and odd parts using Schur decomposition, assigning the
	   * non-trivial inversion to the even sites (see DeGran/DeTar p 174ff).
	   */
	  //init some helping buffers
	  const hardware::buffers::Spinor clmem_source_even  (meta::get_eoprec_spinorfieldsize(get_parameters()), solver->get_device());
	  const hardware::buffers::Spinor clmem_source_odd  (meta::get_eoprec_spinorfieldsize(get_parameters()), solver->get_device());
	  const hardware::buffers::Spinor clmem_tmp_eo_1  (meta::get_eoprec_spinorfieldsize(get_parameters()), solver->get_device());
	  const hardware::buffers::Spinor clmem_tmp_eo_2  (meta::get_eoprec_spinorfieldsize(get_parameters()), solver->get_device());
	  const hardware::buffers::Plain<hmc_complex> clmem_one (1, solver->get_device());
	  hmc_complex one = hmc_complex_one;
	  clmem_one.load(&one);
	  
	  //convert source and input-vector to eoprec-format
	  solver->convert_to_eoprec_device(&clmem_source_even, &clmem_source_odd, source);
	  //prepare sources
	  /**
	   * This changes the even source according to (with A = M + D):
	   *  b_e = b_e - D_eo M_inv b_o
	   */
	  
	  if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
	    //in this case, the diagonal matrix is just 1 and falls away.
	    solver->dslash_eo_device(&clmem_source_odd, &clmem_tmp_eo_1, gf, EVEN);
	    solver->saxpy_eoprec_device(&clmem_source_even, &clmem_tmp_eo_1, &clmem_one, &clmem_source_even);
	  } else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
	    solver->M_tm_inverse_sitediagonal_device(&clmem_source_odd, &clmem_tmp_eo_1);
	    solver->dslash_eo_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, gf, EVEN);
	    solver->saxpy_eoprec_device(&clmem_source_even, &clmem_tmp_eo_2, &clmem_one, &clmem_source_even);
	  }
	  
	  //Trial solution
	  ///@todo this should go into a more general function
	  solver->set_eoprec_spinorfield_cold_device(solver->get_inout_eo());
	  logger.debug() << "start eoprec-inversion";
	  //even solution
	  if(get_parameters().get_solver() == meta::Inputparameters::cg){
	    //to use cg, one needs an hermitian matrix, which is QplusQminus
	    //the source must now be gamma5 b, to obtain the desired solution in the end
	    solver->gamma5_eo_device(&clmem_source_even);
	    ::QplusQminus_eo f_eo(solver);
	    converged = solver->cg_eo(f_eo, solver->get_inout_eo(), &clmem_source_even, gf, get_parameters().get_solver_prec());
	    //now, calc Qminus inout to obtain x = A^⁻1 b
	    //therefore, use source as an intermediate buffer
	    solver->Qminus_eo(solver->get_inout_eo(), &clmem_source_even, gf, get_parameters().get_kappa(), meta::get_mubar(get_parameters() ));
	    //save the result to inout
	    hardware::buffers::copyData(solver->get_inout_eo(), &clmem_source_even);
	  } else{
	    ::Aee f_eo(solver);
	    converged = solver->bicgstab_eo(f_eo, solver->get_inout_eo(), &clmem_source_even, gf, get_parameters().get_solver_prec());
	  }
	  
	  //odd solution
	  /** The odd solution is obtained from the even one according to:
	   *  x_o = M_inv D x_e - M_inv b_o
	   */
	  if(get_parameters().get_fermact() == meta::Inputparameters::wilson) {
	    //in this case, the diagonal matrix is just 1 and falls away.
	    solver->dslash_eo_device(solver->get_inout_eo(), &clmem_tmp_eo_1, gf, ODD);
	    solver->saxpy_eoprec_device(&clmem_tmp_eo_1, &clmem_source_odd, &clmem_one, &clmem_tmp_eo_1);
	  } else if(get_parameters().get_fermact() == meta::Inputparameters::twistedmass) {
	    solver->dslash_eo_device(solver->get_inout_eo(), &clmem_tmp_eo_2, gf, ODD);
	    solver->M_tm_inverse_sitediagonal_device(&clmem_tmp_eo_2, &clmem_tmp_eo_1);
	    solver->M_tm_inverse_sitediagonal_device(&clmem_source_odd, &clmem_tmp_eo_2);
	    solver->saxpy_eoprec_device(&clmem_tmp_eo_1, &clmem_tmp_eo_2, &clmem_one, &clmem_tmp_eo_1);
	  }
	  //CP: whole solution
	  //CP: suppose the even sol is saved in inout_eoprec, the odd one in clmem_tmp_eo_1
	  solver->convert_from_eoprec_device(solver->get_inout_eo(), &clmem_tmp_eo_1, inout);
	}

	if(get_parameters().get_profile_solver() ) {
		solver->get_device()->synchronize();
		(*solvertimer).add();
	}

	if (converged < 0) {
		if(converged == -1) logger.fatal() << "\t\t\tsolver did not solve!!";
		else logger.fatal() << "\t\t\tsolver got stuck after " << abs(converged) << " iterations!!";
	} else logger.debug() << "\t\t\tsolver solved in " << converged << " iterations!";
}

void Gaugefield_inverter::perform_inversion(usetimer* solver_timer)
{
  int num_sources = get_parameters().get_num_sources();

	Opencl_Module_Fermions * solver = get_task_solver();
	const hardware::buffers::Plain<spinor> clmem_res(meta::get_spinorfieldsize(get_parameters()), solver->get_device());
	const hardware::buffers::Plain<spinor> clmem_source(meta::get_spinorfieldsize(get_parameters()), solver->get_device());

	//apply stout smearing if wanted
	if(get_parameters().get_use_smearing() == true)
		solver->smear_gaugefield(solver->get_gaugefield(), std::vector<const hardware::buffers::SU3 *>());

	for(int k = 0; k < num_sources; k++) {
		//copy source from to device
		//NOTE: this is a blocking call!
		logger.debug() << "copy pointsource between devices";
		clmem_source.load(&source_buffer[k * meta::get_vol4d(get_parameters())]);
		logger.debug() << "calling solver..";
		invert_M_nf2_upperflavour( &clmem_res, &clmem_source, solver->get_gaugefield(), solver_timer);
		//add solution to solution-buffer
		//NOTE: this is a blocking call!
		logger.debug() << "add solution...";
		clmem_res.dump(&solution_buffer[k * meta::get_vol4d(get_parameters())]);
	}

	if(get_parameters().get_use_smearing() == true) 
		solver->unsmear_gaugefield(solver->get_gaugefield());
}

void Gaugefield_inverter::flavour_doublet_correlators(std::string corr_fn)
{
	using namespace std;
	using namespace hardware::buffers;

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

	hmc_float* host_result_ps = new hmc_float [num_corr_entries];
	hmc_float* host_result_sc = new hmc_float [num_corr_entries];
	hmc_float* host_result_vx = new hmc_float [num_corr_entries];
	hmc_float* host_result_vy = new hmc_float [num_corr_entries];
	hmc_float* host_result_vz = new hmc_float [num_corr_entries];
	hmc_float* host_result_ax = new hmc_float [num_corr_entries];
	hmc_float* host_result_ay = new hmc_float [num_corr_entries];
	hmc_float* host_result_az = new hmc_float [num_corr_entries];

	const Plain<hmc_float> result_ps(num_corr_entries, get_task_correlator()->get_device());
	const Plain<hmc_float> result_sc(num_corr_entries, get_task_correlator()->get_device());
	const Plain<hmc_float> result_vx(num_corr_entries, get_task_correlator()->get_device());
	const Plain<hmc_float> result_vy(num_corr_entries, get_task_correlator()->get_device());
	const Plain<hmc_float> result_vz(num_corr_entries, get_task_correlator()->get_device());
	const Plain<hmc_float> result_ax(num_corr_entries, get_task_correlator()->get_device());
	const Plain<hmc_float> result_ay(num_corr_entries, get_task_correlator()->get_device());
	const Plain<hmc_float> result_az(num_corr_entries, get_task_correlator()->get_device());

	logger.info() << "calculate correlators..." ;
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ps"), get_clmem_corr(), &result_ps);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("sc"), get_clmem_corr(), &result_sc);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vx"), get_clmem_corr(), &result_vx);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vy"), get_clmem_corr(), &result_vy);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("vz"), get_clmem_corr(), &result_vz);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ax"), get_clmem_corr(), &result_ax);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("ay"), get_clmem_corr(), &result_ay);
	get_task_correlator()->correlator_device(get_task_correlator()->get_correlator_kernel("az"), get_clmem_corr(), &result_az);

	//the pseudo-scalar (J=0, P=1)
	logger.info() << "pseudo scalar correlator:" ;
	result_ps.dump(host_result_ps);
	for(int j = 0; j < num_corr_entries; j++) {
		logger.info() << j << "\t" << scientific << setprecision(14) << host_result_ps[j];
		of << scientific << setprecision(14) << "0 1\t" << j << "\t" << host_result_ps[j] << endl;
	}


	//the scalar (J=0, P=0)
	result_sc.dump(host_result_sc);
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "0 0\t" << j << "\t" << host_result_sc[j] << endl;
	}


	//the vector (J=1, P=1)
	result_vx.dump(host_result_vx);
	result_vy.dump(host_result_vy);
	result_vz.dump(host_result_vz);
	for(int j = 0; j < num_corr_entries; j++) {
		of << scientific << setprecision(14) << "1 1\t" << j << "\t" << (host_result_vx[j] + host_result_vy[j] + host_result_vz[j]) / 3. << "\t" << host_result_vx[j] << "\t" << host_result_vy[j] << "\t" << host_result_vz[j] << endl;
	}


	//the axial vector (J=1, P=0)
	result_ax.dump(host_result_ax);
	result_ay.dump(host_result_ay);
	result_az.dump(host_result_az);
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
}

void Gaugefield_inverter::create_sources()
{
	//create sources on the correlator-device and save them on the host
	if(get_parameters().get_use_pointsource() == true) {
		logger.debug() << "start creating point-sources...";
		for(int k = 0; k < get_parameters().get_num_sources(); k++) {
		  get_task_correlator()->create_point_source_device(get_clmem_source_corr(), k, get_source_pos_spatial(get_parameters()), get_parameters().get_source_t());
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
