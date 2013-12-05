#include "heatbathExecutable.h"

class heatbathExecutable: public generalExecutable
{
public:
	heatbathExecutable(int argc, const char* argv[]) : generalExecutable(argc, argv)
	{
		initializationTimer.reset();
		prng = new physics::PRNG(*system);
		gaugefield = new physics::lattices::Gaugefield(*system, *prng);
		meta::print_info_heatbath(ownName, parameters);
		writeHeatbathLogfile();
		setIterationParameters();
		initializationTimer.add();
	}

	void performHeatbathAndMeasureGaugeObservables()
	{
		performanceTimer.reset();
		performThermalization();
		performHeatbathAndMeasurements();
		performanceTimer.add();
	}

private:
	physics::PRNG * prng;
	physics::lattices::Gaugefield * gaugefield;
	const std::string 	filenameForHeatbathLogfile 		= "heatbath.log";
	std::string filenameForGaugeobservables;
	int thermalizationSteps;
	int heatbathSteps;
	int overrelaxSteps;
	int writeFrequency;
	int saveFrequency;

	void writeHeatbathLogfile()
	{
		outputToFile.open(filenameForHeatbathLogfile, std::ios::out | std::ios::app);
		if (outputToFile.is_open()) {
			meta::print_info_heatbath(ownName, &outputToFile, parameters);
			outputToFile.close();
		} else {
			throw File_Exception(filenameForHeatbathLogfile);
		}
	}

	void setIterationParameters()
	{
		thermalizationSteps = 	parameters.get_thermalizationsteps();
		heatbathSteps = 		parameters.get_heatbathsteps();
		overrelaxSteps = 		parameters.get_overrelaxsteps();
		writeFrequency = 		parameters.get_writefrequency();
		saveFrequency = 		parameters.get_savefrequency();
	}

	void performThermalization()
	{
		logger.info() << "Start thermalization";
		int iteration = 0;
		writeGaugeObservablesToScreen(iteration);
		for (; iteration < thermalizationSteps; iteration++)
		{
			physics::algorithms::heatbath(*gaugefield, *prng);
		}
		logger.info() << "thermalization done";
	}

	void writeGaugeObservablesToFile(int& iteration)
	{
		if (((iteration + 1) % writeFrequency) == 0) {
			filenameForGaugeobservables = meta::get_gauge_obs_file_name(parameters, "");
			print_gaugeobservables(*gaugefield, iteration, filenameForGaugeobservables);
		}
	}

	void writeGaugeObservablesToScreen(int& iteration)
	{
		if ( parameters.get_print_to_screen() || (iteration == 0) ) {
				print_gaugeobservables(*gaugefield, iteration);
		}
	}

	void writeGaugeObservablesToScreenAndFile(int iteration)
	{
		writeGaugeObservablesToScreen(iteration);
		writeGaugeObservablesToFile(iteration);
	}

	void saveGaugefield(int iteration)
	{
		if (((saveFrequency != 0) && ((iteration + 1) % saveFrequency) == 0)
				|| (iteration == heatbathSteps - 1)) {
			gaugefield->save(iteration + 1);
		}
	}

	void measureGaugeObservables(int& iteration)
	{
		writeGaugeObservablesToScreenAndFile(iteration);
	}

	void performHeatbathAndMeasurements() {
		logger.info() << "Start heatbath";
		for (int iteration = 0; iteration < heatbathSteps; iteration++) {
			physics::algorithms::heatbath(*gaugefield, *prng, overrelaxSteps);
			measureGaugeObservables(iteration);
			saveGaugefield(iteration);
		}
		logger.info() << "heatbath done";
	}
};

int main(int argc, const char* argv[])
{
	try {
		heatbathExecutable heatbathInstance(argc, argv);
		heatbathInstance.performHeatbathAndMeasureGaugeObservables();
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
	} catch (Invalid_Parameters& es) {
		logger.fatal() << es.what();
		exit(1);
	}

	return 0;

}
