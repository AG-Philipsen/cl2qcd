#include <fstream>
#include "../meta/util.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../physics/lattices/staggeredfield_eo.hpp"
#include "../physics/fermionmatrix/fermionmatrix_stagg.hpp"
#include "../klepsydra/klepsydra.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/fermions_staggered.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, const char* argv[])
{
	using physics::fermionmatrix::DKS_eo;

	try {
		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		//meta::print_info_inverter(argv[0], parameters);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		hardware::System system(parameters);

		physics::PRNG prng(system);
		physics::lattices::Gaugefield gf(system, prng);
		physics::lattices::Staggeredfield_eo sf1(system);
		physics::lattices::Staggeredfield_eo sf2(system);

		// update gaugefield buffers once to have update links fully initialized
		gf.update_halo();

		logger.info() << "Gaugeobservables:";
		print_gaugeobservables(gf, 0);

		DKS_eo(&sf2, gf, sf1, EVEN);
		DKS_eo(&sf1, gf, sf2, ODD);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// DKS_eo-benchmark
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for(auto dev: system.get_devices()) {
			dev->synchronize();
		}

		int hmc_iter = parameters.get_hmcsteps();

		logger.info() << "Perform DKS_eo (EVEN + ODD) " << hmc_iter << " times.";
		klepsydra::Monotonic timer;
		for(int iter = 0; iter < hmc_iter; ++iter) {
			DKS_eo(&sf2, gf, sf1, EVEN);
			DKS_eo(&sf1, gf, sf2, ODD);
		}
		for(auto dev: system.get_devices()) {
			dev->synchronize();
		}
		auto elapsed_mus = timer.getTime();
		logger.trace() << "DKS_eo benchmarking done" ;

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		auto fermion_code = system.get_devices()[0]->get_fermion_staggered_code();
		size_t flop_count = fermion_code->get_flop_size("D_KS_eo");
		size_t byte_count = fermion_code->get_read_write_size("D_KS_eo");
		double gflops = static_cast<double>(flop_count) * 2 * hmc_iter / elapsed_mus / 1e3;
		double gbytes = static_cast<double>(byte_count) * 2 * hmc_iter / elapsed_mus / 1e3;
		logger.info() << "D_KS performance: " << gflops << " GFLOPS";
		logger.info() << "D_KS memory: " << gbytes << " GB/S";
		logger.info() << "Measured TIME: " << elapsed_mus / 1.e3 << "msec";

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
}

