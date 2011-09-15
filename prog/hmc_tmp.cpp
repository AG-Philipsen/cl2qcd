#include "hmc_tmp.h"

int main(int argc, char* argv[])
{
	try {

		if(argc != 2) throw Print_Error_Message("Need file name for input parameters", __FILE__, __LINE__);

		char* inputfile = argv[1];
		inputparameters parameters;
		parameters.readfile(inputfile);
		parameters.print_info_inverter(argv[0]);

		ofstream ofile;
		ofile.open("hmc.log");
		if(ofile.is_open()) {
			parameters.print_info_hmc(argv[0], &ofile);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for hmc.";
		}

		//name of file to store gauge observables, print initial information
		/** @todo think about what is a senseful filename*/
		stringstream gaugeout_name;
		gaugeout_name << "hmc_output";
		
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		hmc_observables obs;
		Gaugefield_hmc gaugefield;

		//use 1 task: the hmc-algorithm
		int numtasks = 1;
		if(parameters.get_num_dev() == 2 )
			logger.warn() << "Only 1 device demanded by input file. All calculations performed on primary device.";

		cl_device_type primary_device;
		switch ( parameters.get_use_gpu() ) {
			case true :
				primary_device = CL_DEVICE_TYPE_GPU;
				break;
			case false :
				primary_device = CL_DEVICE_TYPE_CPU;
				break;
		}

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device, &parameters);


		logger.info() << "Gaugeobservables:";
		gaugefield.print_gaugeobservables(0);
		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// hmc
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();
		/** @todo usage of solver_timer has to be checked. No output yet */
		usetimer solver_timer;
		
		int hmc_iter = parameters.get_hmcsteps();
		int iter;	
		hmc_float acc_rate = 0.;
		int writefreq = parameters.get_writefrequency();
		int savefreq = parameters.get_savefrequency();
		//This is the random-number generator for the metropolis-step
		Random hmc_rnd_gen (parameters.get_host_seed());

		logger.trace() << "perform HMC on device(s)... ";

		//main hmc-loop
		for(iter = 0; iter < hmc_iter; iter ++) {
			//generate new random-number for Metropolis step
			hmc_float rnd_number = hmc_rnd_gen.doub();
			gaugefield.perform_hmc_step(&obs, iter, rnd_number, &solver_timer);
			acc_rate += obs.accept;
			if( ( (iter + 1) % writefreq ) == 0 ) {
				gaugefield.print_hmcobservables(obs, iter, gaugeout_name.str());
			}
			if( parameters.get_saveconfigs() == true && ( (iter + 1) % savefreq ) == 0 ) {
				//CP: I dont think this is necessary at the moment...
// 				gaugefield.synchronize(0);
				gaugefield.save(iter);
			}
		}
		logger.trace() << "HMC done";
		logger.trace() << "Acceptance rate: " << fixed <<  setprecision(1) << percent(acc_rate, hmc_iter) << "%";
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		uint64_t totaltime = total_timer.getTime();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);
		//print times from the devices...
		logger.trace() << "## Device: HMC [0]";
		(gaugefield.get_task_hmc(0))->print_copy_times(totaltime);
		/// @todo add more than one device if used

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// free variables
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		gaugefield.finalize();

	} //try
	//exceptions from Opencl classes
	catch (Opencl_Error& e) {
		logger.fatal() << e.what();
		exit(1);
	} catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(1);
	} catch (Print_Error_Message& em) {
		logger.fatal() << em.what();
		exit(1);
	}

	return 0;

}
