#include "inverter.h"

void print_hello_bench(char* name)
{
	std::cout << "This is the Inverter-Benchmarking Programme " << name << endl;
	std::cout << "NOTE: All output times are only taken in benchmarkig-loop (except OpenCL-init) " << name << endl;
	return;
}

int main(int argc, char* argv[])
{

	if(argc != 2) {
		cerr << "need file name for input parameters" << endl;
		return HMC_FILEERROR;
	}

	char* progname = argv[0];
	print_hello_bench(progname);

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	print_info(&parameters,&cout);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sourcefileparameters parameters_source;

	Gaugefield_inversion gaugefield;
	cl_device_type devicetypes[1];

	//CP: this is the timer used outside the main loop
	usetimer noop;
	
#ifdef _USEGPU_
	devicetypes[0] = CL_DEVICE_TYPE_GPU;
#else
	devicetypes[0] = CL_DEVICE_TYPE_CPU;
#endif
	cerr << "init gaugefield" << endl;
	gaugefield.init(1, devicetypes, &parameters, &inittime);
	//cerr << "print initial gaugeobservables..." << endl;
	//	gaugefield.print_gaugeobservables(&polytime, &plaqtime);
	cerr << "copy gaugefield" << endl;
	gaugefield.copy_gaugefield_to_devices(&noop);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// inverter-benchmark
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

//kernels existing in Opencl_fermions, sorted roughly by groups
/*
	//fermionmatrix
	cl_kernel M;
	cl_kernel gamma5;
	cl_kernel Qplus;
	cl_kernel Qminus;
	cl_kernel ps_correlator;
	cl_kernel set_spinorfield_cold;
	cl_kernel gamma5_eoprec;
	cl_kernel M_diag;
	cl_kernel dslash;
	cl_kernel M_sitediagonal;
	cl_kernel M_inverse_sitediagonal;
	cl_kernel dslash_eoprec;
		
	//BLAS
	cl_kernel set_eoprec_spinorfield_cold;
	cl_kernel convert_from_eoprec;
	cl_kernel saxpy;
	cl_kernel saxsbypz;
	cl_kernel set_zero_spinorfield;
	cl_kernel convert_to_kappa_format;
	cl_kernel convert_from_kappa_format;
	cl_kernel convert_to_kappa_format_eoprec;
	cl_kernel convert_from_kappa_format_eoprec;
	cl_kernel create_point_source;
	cl_kernel saxpy_eoprec;
	cl_kernel saxsbypz_eoprec;
	cl_kernel set_zero_spinorfield_eoprec;
	cl_kernel create_point_source_eoprec;
			
	//Scalar Product
	cl_kernel scalar_product;
	cl_kernel scalar_product_reduction;
	cl_kernel global_squarenorm;
	cl_kernel global_squarenorm_reduction;
	cl_kernel scalar_product_eoprec;
	cl_kernel global_squarenorm_eoprec;
	
	//Single
	cl_kernel ratio;
	cl_kernel product;
*/

	//This is just copied from the inverter-function...
	
		//global and local work sizes;
	//LZ: should eventually be moved inside opencl_fermions class
#ifdef _USEGPU_
	const size_t ls = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
#else
	const size_t ls = 1; // nothing else makes sens on CPU
#endif

#ifdef _USEGPU_
	size_t gs = 4 * NUMTHREADS * gaugefield.get_devices()[0].max_compute_units; /// @todo autotune
#else
	size_t gs = gaugefield.get_devices()[0].max_compute_units;
#endif

	const cl_uint num_groups = (gs + ls - 1) / ls;
	gs = ls * num_groups;



	int use_eo = gaugefield.get_parameters()->get_use_eo();

	logger.trace() << "prepare fields for benchmarking";
  if(use_eo==FALSE){
		gaugefield.get_devices_fermions()[0].set_spinorfield_cold_device(ls, gs, &noop);
		gaugefield.get_devices_fermions()[0].create_point_source_device(0,0, 0,ls, gs, &noop);
      
  }
  else{
		gaugefield.get_devices_fermions()[0].set_eoprec_spinorfield_cold_device(ls, gs, &noop);		
		gaugefield.get_devices_fermions()[0].create_point_source_eoprec_device(0,0,0,ls, gs, &noop, &noop, &noop);
	}
	int steps = parameters.get_hmcsteps();
	int cgmax = parameters.get_cgmax();
	logger.trace() << "perform " << steps << " steps of Inverter-Benchmarking on device..";

	totaltime.reset();
	for(int i = 0; i<steps; i++){
		if(use_eo == FALSE) {
// 			logger.trace() << "benchmark fermionmatrix kernels";
			Mtimer.reset();
// 			gaugefield.get_devices_fermions()[0].M_device(gaugefield.get_devices_fermions()[0].get_clmem_inout(), gaugefield.get_devices_fermions()[0].get_clmem_tmp(), ls, gs,  &noop, &noop, &noop);
			Mtimer.add();
// 			logger.trace() << "benchmark BLAS kernels";
			latimer.reset();
			
			latimer.add();
// 			logger.trace() << "benchmark ScalarProduct kernels";
			scalarprodtimer.reset();
			
			scalarprodtimer.add();
// 			logger.trace() << "benchmark Single kernels";
			singletimer.reset();
			
			singletimer.add();
			
			//the input spinorfield should always be the same
			gaugefield.get_devices_fermions()[0].set_spinorfield_cold_device(ls, gs, &noop);
			solvertimer.reset();
			gaugefield.get_devices_fermions()[0].bicgstab_device(&copytimer, &singletimer, &Mtimer, &scalarprodtimer, &latimer, &dslashtimer, &Mdiagtimer, ls, gs, cgmax);
// 			gaugefield.get_devices_fermions()[0].cg_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
			solvertimer.add();
		}
		else{
// 			logger.trace() << "benchmark fermionmatrix kernels";
			Mtimer.reset();
			
			
			
			Mtimer.add();
// 			logger.trace() << "benchmark BLAS kernels";
			latimer.reset();
			
			latimer.add();
// 			logger.trace() << "benchmark ScalarProduct kernels";
			scalarprodtimer.reset();
			
			scalarprodtimer.add();
// 			logger.trace() << "benchmark Single kernels";
			singletimer.reset();
			
			singletimer.add();
		}
		
	}
	totaltime.add();
	cout << "inverter-Benchmarking done" << endl;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	time_output_inverter(
   &totaltime,  &inittime,  &polytime,  &plaqtime,  &updatetime,  &overrelaxtime,  &copytime
  ,  &ferm_inittime, &singletimer, &Mtimer, &copytimer, &scalarprodtimer, &latimer,  &solvertimer,  &dslashtimer,  &Mdiagtimer);

	int err = gaugefield.finalize();
	if(err!=HMC_SUCCESS)
		logger.fatal() << "error in finalizing gaugefield";
	
	
	
	return HMC_SUCCESS;
}
