#include "gaugeobservables.h"

#include "meta/util.hpp"

int main(int argc, const char* argv[])
{
	try {
	  //@todo Debug this!
	  /*
	  logger.info() << "This executable requires the following parameter value(s) to work properly:";
	  logger.info() << "startcondition:\tcontinue";
	  logger.info() << "These will be parsed and possibly overwrite inputfile values!";
	  std::string startcond_opt = "--startcondition=continue";
	  const char* _params[] = {(const char*) argv,   startcond_opt.c_str()};
	  meta::Inputparameters parameters( (int) argc+1, _params);
	  */
	  meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		meta::print_info_gaugeobservables(argv[0], parameters);

		ofstream ofile;
		ofile.open("gaugeobservables.log");
		if(ofile.is_open()) {
			meta::print_info_gaugeobservables(argv[0], &ofile, parameters);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for gaugeobservables.";
		}

		//name of file to store gauge observables, print initial information
		/** @todo think about what is a senseful filename*/
		stringstream gaugeout_name;
		gaugeout_name << "gaugeobservables.data";

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();

		hardware::System system(parameters);
		Gaugefield_hybrid gaugefield(&system);

		//use 1 task:
		int numtasks = 1;
		if(parameters.get_device_count() == 2 )
			logger.warn() << "Only 1 device demanded by input file. All calculations performed on primary device.";

		cl_device_type primary_device = parameters.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// gaugeobservables
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		int iter_end = parameters.get_config_read_end();
		int iter_start = parameters.get_config_read_start();
		int iter_incr = parameters.get_config_read_incr();
		int iter = 0;

		logger.info() << "Measure gaugeobservables on device(s)... ";

		if(parameters.get_read_multiple_configs()){
		  //main loop
		  for(iter = iter_start; iter < iter_end; iter+=iter_incr) {
		    std::string config_name = gaugefield.create_configuration_name(iter);
		    logger.info() << "Measure gaugeobservables of configuration: " << config_name;
		    gaugefield.init_gaugefield(config_name.c_str());
		    gaugefield.synchronize(0);
		    if(parameters.get_print_to_screen() ){
		      gaugefield.print_gaugeobservables(iter);
		    }
		    gaugefield.print_gaugeobservables_from_task(iter, 0, gaugeout_name.str());
		  }
		}
		else{
		  //in this case only the config from the initialization is taken into account
		  logger.info() << "Measure gaugeobservables of configuration: " << parameters.get_sourcefile();
		  //@todo: adjust the "iter" here to be the number from the sourcefile!!
		  if(parameters.get_print_to_screen() ){
		    gaugefield.print_gaugeobservables(iter);
		  }
		  gaugefield.print_gaugeobservables_from_task(iter, 0, gaugeout_name.str());
		}
		  logger.info() << "... done";
		  perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		uint64_t totaltime = total_timer.getTime();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

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
