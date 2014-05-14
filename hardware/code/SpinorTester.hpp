#include "../hardware/code/kernelTester.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/prng.hpp"
#include "../device.hpp"
#include "spinors.hpp"
#include "complex.hpp"

class SpinorTester : public KernelTester {
public:
  //todo: move to .cpp eventually
	SpinorTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
		inputfile(getSpecificInputfile(inputfileIn)), KernelTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues)
		{
		//todo: this object should be a member of KernelTester!
		meta::Inputparameters parameters = createParameters(inputfile);

		system = new hardware::System(parameters);
		device = system->get_devices()[0];

		code = device->get_spinor_code();
		
		prng = new physics::PRNG(*system);

		doubleBuffer = new hardware::buffers::Plain<double> (1, device);
		
		NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(parameters);
		NUM_ELEMENTS_EO = hardware::code::get_eoprec_spinorfieldsize(parameters);
		(parameters.get_solver() == meta::Inputparameters::cg) ? useRandom = false : useRandom =true;
		(parameters.get_read_multiple_configs() ) ? evenOrOdd = true : evenOrOdd = false;
		alpha_host = {parameters.get_beta(), parameters.get_rho()};
		beta_host = {parameters.get_kappa(), parameters.get_mu()};
		ns = parameters.get_nspace();
		nt = parameters.get_ntime();
		iterations = parameters.get_integrationsteps(0);
		parameters.get_read_multiple_configs() ? calcVariance=false : calcVariance = true;
	}
protected:
	std::string inputfile;
	
        std::string getSpecificInputfile(std::string inputfileIn);	
        spinor * createSpinorfield(size_t numberOfElements, int seed = 123456);
        void fill_with_one(spinor * in, int size);
        void fill_with_random(spinor * in, int size, int seed);
        spinor * createSpinorfieldWithOnesAndZerosDependingOnSiteParity();	
        void fill_with_one_eo(spinor * in, int size, bool eo);
        hmc_float count_sf(spinor * in, int size);
        hmc_float calc_var(hmc_float in, hmc_float mean);
        hmc_float calc_var_sf(spinor * in, int size, hmc_float sum);
        void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<spinor> * in);
        void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::Spinor * in);

	//todo: use the fct. from geometry.h here!
	int get_nspace(int* coord)
	{
		int n = 0;
		n = ns * ns * coord[3] + ns * coord[2] + coord[1];
		return n;
	}
	
	int get_global_pos(int spacepos, int t)
	{
		return spacepos + ns*ns*ns * t;
	}
	
	const hardware::System * system;
	hardware::Device * device;
	const hardware::code::Spinors * code;
	physics::PRNG * prng;

        hardware::buffers::Plain<double> * doubleBuffer;
	
	size_t NUM_ELEMENTS_SF;
	size_t NUM_ELEMENTS_EO;
	bool useRandom;
	bool evenOrOdd;
	bool calcVariance;
	hmc_complex alpha_host;
	hmc_complex beta_host;
	int ns;
	int nt;
	int iterations;
};
