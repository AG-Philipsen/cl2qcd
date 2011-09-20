#include "inputparameters.h"

#include "host_geometry.h"

void inputparameters::set_defaults()
{
	//global parameters
	prec = 64;
	use_rec12 = false;
	use_gpu = false;
	use_chem_pot_re = false;
	use_chem_pot_im = false;
	use_smearing = false;
	nspace = 4;
	ntime = 8;
	volspace = nspace * nspace * nspace;
	vol4d = volspace * ntime;
	spinorsize = 12;
	halfspinorsize = 6;
	spinorfieldsize = vol4d;
	eoprec_spinorfieldsize = spinorfieldsize / 2;
	su3algebrasize = NC * NC - 1;
	gaugemomentasize = NDIM * vol4d;
	gaugefieldsize = NC * NC * NDIM * vol4d;

	startcondition = COLD_START;
	saveconfigs = false;
	writefrequency = 1;
	savefrequency = 100;
	num_dev = 1;
	sourcefilenumber = "00000";
	print_to_screen = false;
	//This is obvious!!!
	host_seed = 4815162342;
	use_autotuning = false;

#ifdef _PROFILING_
	mat_size = 9;
	float_size = 8;
#endif
	//gaugefield parameters
	beta = 4.0;
	theta_gaugefield = 0.;
	rho = 0.;
	rho_iter = 0;

	//heatbath parameters
	thermalizationsteps = 0;
	heatbathsteps = 1000;
	overrelaxsteps = 1;

	//fermionic parameters
	fermact = WILSON;
	kappa = 0.125;
	mu = 0.006;
	csw = 0.;
	cgmax = 1000;
	theta_fermion_spatial = 0.;
	theta_fermion_temporal = 0.;
	chem_pot_re = 0.;
	chem_pot_im = 0.;
	use_eo = true;
	//at the moment, only 2 solvers are implemented..
	use_cg = false;
	use_bicgstab_save = false;
	use_pointsource = true;
	num_sources = 12;
	pointsource_x = 0;
	pointsource_y = 0;
	pointsource_z = 0;
	pointsource_t = 0;
#ifdef _USEDOUBLEPREC_
	solver_prec = 1e-23;
#else
	solver_prec = 1e-16;
#endif
	iter_refresh = 100;

	//HMC specific parameters
	tau = 0.5;
	reversibility_check = false;
	integrationsteps1 = 10;
	integrationsteps2 = integrationsteps1;
	hmcsteps = 10;
	num_timescales = 1;
	integrator = LEAPFROG;
	lambda1 = 0.1931833275037836;
	lambda2 = 0.1931833275037836;

	//direction for the correlator
	corr_dir = 3;

	return;
}

void inputparameters::readfile(const char* ifn)
{
	try {
		std::ifstream infile;
		infile.open(ifn);
		if(!infile.is_open()) throw File_Exception(ifn);

		bool muset  = false;
		bool cswset = false;

		while (infile.good()) {
			char linebuf[256];
			infile.getline(linebuf, 256);
			std::string line(linebuf);
			if(line.length() == 255) {
				/// @todo Handle such lines properly
				logger.fatal() << "The file contains a line longer than 255 characters - bailing out";
			}
			if(line.find("#") != std::string::npos) continue; //allow comments
			if(line.find("kappa") != std::string::npos) val_assign(&kappa, line);
			if(line.find("Kappa") != std::string::npos) val_assign(&kappa, line);
			if(line.find("mu") != std::string::npos) {
				val_assign(&mu, line);
				muset = true;
			}
			if(line.find("Mu") != std::string::npos) {
				val_assign(&mu, line);
				muset = true;
			}
			if(line.find("csw") != std::string::npos) {
				val_assign(&csw, line);
				cswset = true;
			}
			if(line.find("Csw") != std::string::npos) {
				val_assign(&csw, line);
				cswset = true;
			}
			if(line.find("beta") != std::string::npos) val_assign(&beta, line);
			if(line.find("tau") != std::string::npos) val_assign(&tau, line);
			if(line.find("Beta") != std::string::npos) val_assign(&beta, line);

			if(line.find("theta_fermion_spatial") != std::string::npos) val_assign(&theta_fermion_spatial, line);
			if(line.find("theta_fermion_temporal") != std::string::npos) val_assign(&theta_fermion_temporal, line);
			if(line.find("theta_spatial") != std::string::npos) val_assign(&theta_fermion_spatial, line);
			if(line.find("theta_temporal") != std::string::npos) val_assign(&theta_fermion_temporal, line);
			if(line.find("theta_s") != std::string::npos) val_assign(&theta_fermion_spatial, line);
			if(line.find("theta_t") != std::string::npos) val_assign(&theta_fermion_temporal, line);
			if(line.find("thetas") != std::string::npos) val_assign(&theta_fermion_spatial, line);
			if(line.find("thetat") != std::string::npos) val_assign(&theta_fermion_temporal, line);
			if(line.find("ThetaS") != std::string::npos) val_assign(&theta_fermion_spatial, line);
			if(line.find("ThetaT") != std::string::npos) val_assign(&theta_fermion_temporal, line);

			if(line.find("cgmax") != std::string::npos) val_assign(&cgmax, line);
			if(line.find("CGmax") != std::string::npos) val_assign(&cgmax, line);
			if(line.find("Cgmax") != std::string::npos) val_assign(&cgmax, line);

			if(line.find("Solver") != std::string::npos) solver_assign(&use_cg, line);
			if(line.find("solver") != std::string::npos) solver_assign(&use_cg, line);
			if(line.find("SOLVER") != std::string::npos) solver_assign(&use_cg, line);

			if(line.find("writefrequency") != std::string::npos) val_assign(&writefrequency, line);
			if(line.find("savefrequency") != std::string::npos) val_assign(&savefrequency, line);
			if(line.find("saveconfigs") != std::string::npos) bool_assign(&saveconfigs, line);
			if(line.find("prec") != std::string::npos) val_assign(&prec, line);
			if(line.find("Prec") != std::string::npos) val_assign(&prec, line);
			if(line.find("readsource") != std::string::npos) startcond_assign(&startcondition, line);
			if(line.find("startcondition") != std::string::npos) startcond_assign(&startcondition, line);
			if(line.find("sourcefile") != std::string::npos) {
				val_assign(&sourcefile, line);
				sourcefilenumber_assign(&sourcefilenumber);
			}

			if(line.find("thermalizationsteps") != std::string::npos) val_assign(&thermalizationsteps, line);
			if(line.find("heatbathsteps") != std::string::npos) val_assign(&heatbathsteps, line);
			if(line.find("thermsteps") != std::string::npos) val_assign(&thermalizationsteps, line);
			if(line.find("thermalization") != std::string::npos) val_assign(&thermalizationsteps, line);
			if(line.find("overrelaxsteps") != std::string::npos) val_assign(&overrelaxsteps, line);
			if(line.find("overrelax") != std::string::npos) val_assign(&overrelaxsteps, line);
			if(line.find("oversteps") != std::string::npos) val_assign(&overrelaxsteps, line);
			if(line.find("hmcsteps") != std::string::npos) val_assign(&hmcsteps, line);
			if(line.find("integrationsteps1") != std::string::npos) val_assign(&integrationsteps1, line);
			if(line.find("integrationsteps2") != std::string::npos) val_assign(&integrationsteps2, line);
			if(line.find("num_timescales") != std::string::npos) val_assign(&num_timescales, line);
			if(line.find("num_dev") != std::string::npos) val_assign(&num_dev, line);

			if(line.find("fermaction") != std::string::npos) fermact_assign(&fermact, line);
			if(line.find("fermionaction") != std::string::npos) fermact_assign(&fermact, line);
			if(line.find("fermact") != std::string::npos) fermact_assign(&fermact, line);

			if(line.find("integrator") != std::string::npos) integrator_assign(&integrator, line);

			if(line.find("evenodd") != std::string::npos) bool_assign(&use_eo, line);
			if(line.find("even_odd") != std::string::npos) bool_assign(&use_eo, line);
			if(line.find("even-odd") != std::string::npos) bool_assign(&use_eo, line);
			if(line.find("even-odd-preconditioning") != std::string::npos) bool_assign(&use_eo, line);
			if(line.find("use_eo") != std::string::npos) bool_assign(&use_eo, line);
			if(line.find("use_evenodd") != std::string::npos) bool_assign(&use_eo, line);
			if(line.find("use_rec12") != std::string::npos) bool_assign(&use_rec12, line);
			if(line.find("REC12") != std::string::npos) bool_assign(&use_rec12, line);
			if(line.find("use_gpu") != std::string::npos) bool_assign(&use_gpu, line);
			if(line.find("GPU") != std::string::npos) bool_assign(&use_gpu, line);

			if(line.find("NS") != std::string::npos) val_assign(&nspace, line);
			if(line.find("NSPACE") != std::string::npos) val_assign(&nspace, line);
			if(line.find("NT") != std::string::npos) val_assign(&ntime, line);
			if(line.find("NTIME") != std::string::npos) val_assign(&ntime, line);
			if(line.find("print_to_screen") != std::string::npos) bool_assign(&print_to_screen, line);

			if(line.find("use_pointsource") != std::string::npos) bool_assign(&use_pointsource, line);
			if(line.find("num_sources") != std::string::npos) val_assign(&num_sources, line);
			if(line.find("pointsource_x") != std::string::npos) val_assign(&pointsource_x, line);
			if(line.find("pointsource_y") != std::string::npos) val_assign(&pointsource_y, line);
			if(line.find("pointsource_z") != std::string::npos) val_assign(&pointsource_z, line);
			if(line.find("pointsource_t") != std::string::npos) val_assign(&pointsource_t, line);

			if(line.find("use_smearing") != std::string::npos) bool_assign(&use_smearing, line);
			if(line.find("rho") != std::string::npos) val_assign(&rho, line);
			if(line.find("rho_iter") != std::string::npos) val_assign(&rho_iter, line);
			if(line.find("reversibility_check") != std::string::npos) bool_assign(&reversibility_check, line);
			if(line.find("rev_check") != std::string::npos) bool_assign(&reversibility_check, line);

			if(line.find("autotuning") != std::string::npos) bool_assign(&use_autotuning, line);

			if(line.find("corr_dir") != std::string::npos) val_assign(&corr_dir, line);
			if(line.find("correlator_direction") != std::string::npos) val_assign(&corr_dir, line);

		}

		if( ( (muset  == true) && (fermact != TWISTEDMASS) ) ||
		    ( (cswset == true) && (fermact != CLOVER     ) ) ||
		    ( (cswset == true) && (muset   == true       ) )  )
			throw Invalid_Fermact(fermact, muset, cswset);

		//check the read-in values against the compile time values
		this->set_settings_global();
		this->check_settings_global();

#ifdef _PROFILING_
		//set variables needed for Profiling according to the input-parameters
		if(this->get_use_rec12() != false) this->mat_size = 6;
		if(this->get_prec() == 32 ) this->float_size = 4;
#endif

	} //try
	catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(EXIT_INPUTPARAMETERS);
	} catch (std::string line) {
		logger.fatal() << "Read invalid option in parameter file. Critical line was:";
		logger.fatal() << line;
		logger.fatal() << "Aborting.";
		exit(EXIT_INPUTPARAMETERS);
	} catch (Invalid_Fermact& e) {
		logger.fatal() << e.what();
		logger.fatal() << "Aborting.";
		exit(EXIT_INPUTPARAMETERS);
	}
	return;
}

void inputparameters::val_assign(hmc_float * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);
	(*out) = atof(value.c_str());
	return;
}

void inputparameters::val_assign(int * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);
	(*out) = atoi(value.c_str());
	return;
}

void inputparameters::sourcefilenumber_assign(std::string * out)
{
	//it is supposed that the file is called conf.xxxxx
	size_t length;
	char buffer[20];
	string str ("Test string...");
	length = sourcefile.copy(buffer, 5, 5);
	buffer[length] = '\0';
	//atoi should neglect any letters left in the string
	int tmp = atoi(buffer);
	char buffer2[20];
	//there is no check if the number is bigger than 99999!!
	sprintf(buffer2, "%.5i", tmp + 1);
	(*out) = buffer2;
}

void inputparameters::startcond_assign(int * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);

	if(value.find("cold") != std::string::npos) {
		(*out) = COLD_START;
		return;
	}
	if(value.find("hot") != std::string::npos) {
		(*out) = HOT_START;
		return;
	}
	if(value.find("source") != std::string::npos) {
		(*out) = START_FROM_SOURCE;
		return;
	}
	if(value.find("continue") != std::string::npos) {
		(*out) = START_FROM_SOURCE;
		return;
	}
	throw line;
	return;
}

void inputparameters::fermact_assign(int * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);

	if(value.find("TWISTEDMASS") != std::string::npos) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(value.find("twistedmass") != std::string::npos) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(value.find("Twistedmass") != std::string::npos) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(value.find("TwistedMass") != std::string::npos) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(value.find("clover") != std::string::npos) {
		(*out) = CLOVER;
		return;
	}
	if(value.find("CLOVER") != std::string::npos) {
		(*out) = CLOVER;
		return;
	}
	if(value.find("Clover") != std::string::npos) {
		(*out) = CLOVER;
		return;
	}
	if(value.find("WILSON") != std::string::npos) {
		(*out) = WILSON;
		return;
	}
	if(value.find("Wilson") != std::string::npos) {
		(*out) = WILSON;
		return;
	}
	if(value.find("wilson") != std::string::npos) {
		(*out) = WILSON;
		return;
	}
	if(value.find("unimproved") != std::string::npos) {
		(*out) = WILSON;
		return;
	}
	throw line;
	return;
}

void inputparameters::solver_assign(bool * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);

	if(value.find("BiCGStab") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("bicgstab") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("BICGSTAB") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("BiCGStab_save") != std::string::npos) {
		(*out) = false;
		this->use_bicgstab_save = true;
		return;
	}
	if(value.find("bicgstab_save") != std::string::npos) {
		(*out) = false;
		this->use_bicgstab_save = true;
		return;
	}
	if(value.find("BICGSTAB_save") != std::string::npos) {
		(*out) = false;
		this->use_bicgstab_save = true;
		return;
	}
	if(value.find("BICGSTAB_SAVE") != std::string::npos) {
		(*out) = false;
		this->use_bicgstab_save = true;
		return;
	}
	//LZ: note that the ordering of false/true is crucial here
	//    as any "bicgstab" hit will also be a "cg" hit...
	if(value.find("Cg") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("cg") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("CG") != std::string::npos) {
		(*out) = true;
		return;
	}

	throw line;
	return;
}

void inputparameters::integrator_assign(int * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);

	if(value.find("LEAPFROG") != std::string::npos) {
		(*out) = LEAPFROG;
		return;
	}
	if(value.find("leapfrog") != std::string::npos) {
		(*out) = LEAPFROG;
		return;
	}
	if(value.find("2MN") != std::string::npos) {
		(*out) = TWOMN;
		return;
	}
	if(value.find("2mn") != std::string::npos) {
		(*out) = TWOMN;
		return;
	}
	throw line;
	return;
}

void inputparameters::bool_assign(bool * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);

	if(value.find("1") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("on") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("ON") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("yes") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("TRUE") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("true") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("True") != std::string::npos) {
		(*out) = true;
		return;
	}
	if(value.find("0") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("off") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("OFF") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("no") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("false") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("FALSE") != std::string::npos) {
		(*out) = false;
		return;
	}
	if(value.find("False") != std::string::npos) {
		(*out) = false;
		return;
	}
	throw line;
	return;
}


void inputparameters::val_assign(std::string * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);
	(*out) = value.c_str();
	return;
}

hmc_float inputparameters::get_kappa() const
{
	return kappa;
}

void inputparameters::set_mubar_negative()
{
	mubar *= -1.;
}

void inputparameters::calc_mubar()
{
	mubar = 2.*kappa * mu;
}

hmc_float inputparameters::get_mubar() const
{
	return mubar;
}

hmc_float inputparameters::get_theta_fermion_spatial() const
{
	return theta_fermion_spatial;
}

hmc_float inputparameters::get_theta_fermion_temporal() const
{
	return theta_fermion_temporal;
}

hmc_float inputparameters::get_theta_gaugefield() const
{
	return theta_gaugefield;
}

hmc_float inputparameters::get_beta() const
{
	return beta;
}

hmc_float inputparameters::get_mu() const
{
	return mu;
}

hmc_float inputparameters::get_tau() const
{
	return tau;
}

hmc_float inputparameters::get_csw() const
{
	return csw;
}

hmc_float inputparameters::get_chem_pot_re() const
{
	return chem_pot_re;
}

hmc_float inputparameters::get_chem_pot_im() const
{
	return chem_pot_im;
}

hmc_float inputparameters::get_rho() const
{
	return rho;
}

hmc_float inputparameters::get_lambda1() const
{
	return lambda1;
}

hmc_float inputparameters::get_lambda2() const
{
	return lambda2;
}

int inputparameters::get_cgmax() const
{
	return cgmax;
}

int inputparameters::get_ns() const
{
	return nspace;
}

int inputparameters::get_nt() const
{
	return ntime;
}

int inputparameters::get_prec() const
{
	return prec;
}

int inputparameters::get_num_dev() const
{
	return num_dev;
}

int inputparameters::get_thermalizationsteps() const
{
	return thermalizationsteps;
}

int inputparameters::get_heatbathsteps() const
{
	return heatbathsteps;
}

int inputparameters::get_overrelaxsteps() const
{
	return overrelaxsteps;
}

int inputparameters::get_hmcsteps() const
{
	return hmcsteps;
}

int inputparameters::get_integrationsteps1() const
{
	return integrationsteps1;
}

int inputparameters::get_integrationsteps2() const
{
	return integrationsteps2;
}


int inputparameters::get_writefrequency() const
{
	return writefrequency;
}

int inputparameters::get_savefrequency() const
{
	return savefrequency;
}

int inputparameters::get_startcondition() const
{
	return startcondition;
}

int inputparameters::get_fermact() const
{
	return fermact;
}

int inputparameters::get_integrator() const
{
	return integrator;
}

int inputparameters::get_num_timescales() const
{
	return num_timescales;
}

bool inputparameters::get_saveconfigs() const
{
	return saveconfigs;
}

void inputparameters::display_sourcefile() const
{
	cout << sourcefile;
}

void inputparameters::display_sourcefilenumber() const
{
	cout << sourcefilenumber;
}

bool inputparameters::get_use_eo() const
{
	return use_eo;
}

bool inputparameters::get_use_rec12() const
{
	return use_rec12;
}

bool inputparameters::get_use_gpu() const
{
	return use_gpu;
}

int inputparameters::get_volspace() const
{
	return volspace;
}

int inputparameters::get_vol4d() const
{
	return vol4d;
}

int inputparameters::get_spinorfieldsize() const
{
	return spinorfieldsize;
}

int inputparameters::get_sf_buf_size() const
{
	return sf_buf_size;
}

int inputparameters::get_eo_sf_buf_size() const
{
	return eo_sf_buf_size;
}

int inputparameters::get_gf_buf_size() const
{
	return gf_buf_size;
}

int inputparameters::get_gm_buf_size() const
{
	return gm_buf_size;
}
int inputparameters::get_spinorsize() const
{
	return spinorsize;
}

int inputparameters::get_halfspinorsize() const
{
	return halfspinorsize;
}

int inputparameters::get_gaugemomentasize() const
{
	return gaugemomentasize;
}

int inputparameters::get_su3algebrasize() const
{
	return su3algebrasize;
}

int inputparameters::get_gaugefieldsize() const
{
	return gaugefieldsize;
}

int inputparameters::get_eoprec_spinorfieldsize() const
{
	return eoprec_spinorfieldsize;
}

bool inputparameters::get_use_chem_pot_re() const
{
	return use_chem_pot_re;
}

bool inputparameters::get_use_chem_pot_im() const
{
	return use_chem_pot_im;
}

bool inputparameters::get_use_smearing() const
{
	return use_smearing;
}

int inputparameters::get_host_seed() const
{
	return host_seed;
}

bool inputparameters::get_print_to_screen() const
{
	return print_to_screen;
}

int inputparameters::get_rho_iter() const
{
	return rho_iter;
}

bool inputparameters::get_use_cg() const
{
	return use_cg;
}

bool inputparameters::get_use_bicgstab_save() const
{
	return use_bicgstab_save;
}

bool inputparameters::get_use_pointsource() const
{
	return use_pointsource;
}

int inputparameters::get_num_sources() const
{
	return num_sources;
}

int inputparameters::get_source_pos_spatial() const
{
	int coord [4];
	coord[1] = pointsource_x;
	coord[2] = pointsource_y;
	coord[3] = pointsource_z;

	return get_nspace(coord, this);
}

int inputparameters::get_source_pos_temporal() const
{
	return pointsource_t;
}

#ifdef _PROFILING_
int inputparameters::get_mat_size() const
{
	return mat_size;
}
int inputparameters::get_float_size() const
{
	return float_size;
}
#endif

bool inputparameters::get_use_autotuning() const
{
	return use_autotuning;
}

int inputparameters::get_corr_dir() const
{
	return corr_dir;
}

hmc_float inputparameters::get_solver_prec() const
{
	return solver_prec;
}

uint inputparameters::get_iter_refresh() const
{
	return iter_refresh;
}

bool inputparameters::get_reversibility_check() const
{
	return reversibility_check;
}

void inputparameters::print_info_global() const
{
	logger.info() << "## **********************************************************";
	logger.info() << "## Global parameters:";
	logger.info() << "## NSPACE:  " << this->get_ns();
	logger.info() << "## NTIME:   " << this->get_nt();
	logger.info() << "## NDIM:    " << NDIM;
	logger.info() << "## NCOLOR:  " << NC;
	logger.info() << "## NSPIN:   " << NSPIN;
	logger.info() << "## **********************************************************";
	logger.info() << "## Computational parameters:";
	logger.info() << "## PREC:    " << this->get_prec();
	if(this->get_use_rec12() == true) {
		logger.info() << "## REC12:   ON";
	} else {
		logger.info() << "## REC12:   OFF";
	}
	if(this->get_use_gpu() == true) {
		logger.info() << "## USE GPU: ON";
	} else {
		logger.info() << "## USE GPU: OFF";
	}
	logger.info() << "## Number of devices demanded for calculations: " << this->get_num_dev()  ;
	logger.info() << "## **********************************************************";
	logger.info() << "## I/O parameters:";
	logger.info() << "## SvConf:  " << this->get_saveconfigs();
	logger.info() << "## WrFreq:  " << this->get_writefrequency();
	logger.info() << "## SvFreq:  " << this->get_savefrequency();
	if (this->get_startcondition() == START_FROM_SOURCE) {
		string sf = this->sourcefile;
		logger.info() << "## sourcefile = " << sf;
	}
	if (this->get_startcondition() == COLD_START) {
		logger.info() << "## cold start";
	}
	if (this->get_startcondition() == HOT_START) {
		logger.info() << "## hot start";
	}
	if(this->get_use_smearing() == 1) {
		logger.info() << "## **********************************************************";
		logger.info() << "## Apply Smearing with:";
		logger.info() << "## rho:      " << this->get_rho();
		logger.info() << "## rho_iter: " << this->get_rho_iter();
	}
}

void inputparameters::print_info_global(ostream* os) const
{
	*os  << "## **********************************************************" << endl;
	*os  << "## Global parameters:" << endl;
	*os  << "## NSPACE:  " << this->get_ns() << endl;
	*os  << "## NTIME:   " << this->get_nt() << endl;
	*os  << "## NDIM:    " << NDIM << endl;
	*os  << "## NCOLOR:  " << NC << endl;
	*os  << "## NSPIN:   " << NSPIN << endl;
	*os  << "## **********************************************************" << endl;
	*os  << "## Computational parameters:" << endl;
	*os  << "## PREC:    " << this->get_prec() << endl;
	if(this->get_use_rec12() == true) {
		*os << "## REC12:   ON"  << endl;
	} else {
		*os << "## REC12:   OFF"  << endl;
	}
	if(this->get_use_gpu() == true) {
		*os << "## USE GPU: ON"  << endl;
	} else {
		*os << "## USE GPU: OFF"  << endl;
	}
	*os  << "## Number of devices demanded for calculations: " << this->get_num_dev()  << endl;
	*os  << "## **********************************************************" << endl;
	*os  << "## I/O parameters:" << endl;
	*os  << "## SvConf:  " << this->get_saveconfigs() << endl;
	*os  << "## WrFreq:  " << this->get_writefrequency() << endl;
	*os  << "## SvFreq:  " << this->get_savefrequency() << endl;
	if (this->get_startcondition() == START_FROM_SOURCE) {
		string sf = this->sourcefile;
		*os  << "## sourcefile = " << sf << endl;
	}
	if (this->get_startcondition() == COLD_START) {
		*os  << "## cold start" << endl;
	}
	if (this->get_startcondition() == HOT_START) {
		*os  << "## hot start" << endl;
	}
	if(this->get_use_smearing() == true) {
		*os  << "## **********************************************************" << endl;
		*os  << "## Apply Smearing with:" << endl;
		*os  << "## rho:      " << this->get_rho() << endl;
		*os  << "## rho_iter: " << this->get_rho_iter() << endl;
	}
}


void inputparameters::print_info_heatbath(char* progname) const
{
	logger.info() << "## Starting heatbath program, executable name: " << progname;
	this->print_info_global();
	logger.info() << "## **********************************************************";
	logger.info() << "## Simulation parameters:";
	logger.info() << "## beta           = " << this->get_beta();
	logger.info() << "## thermsteps     = " << this->get_thermalizationsteps() ;
	logger.info() << "## heatbathsteps  = " << this->get_heatbathsteps();
	logger.info() << "## overrelaxsteps = " << this->get_overrelaxsteps();
	logger.info() << "## **********************************************************";
	return;
}


void inputparameters::print_info_heatbath(char* progname, ostream* os) const
{
	*os  << "## Starting heatbath program, executable name: " << progname << endl;
	this->print_info_global(os);
	*os  << "## **********************************************************" << endl;
	*os  << "## Simulation parameters:" << endl;
	*os  << "## beta           = " << this->get_beta() << endl;
	*os  << "## thermsteps     = " << this->get_thermalizationsteps() << endl;
	*os  << "## heatbathsteps  = " << this->get_heatbathsteps() << endl;
	*os  << "## overrelaxsteps = " << this->get_overrelaxsteps() << endl;
	*os  << "## **********************************************************" << endl;
	return;
}


void inputparameters::print_info_tkkappa(char* progname, ostream* os) const
{
	*os << "## Starting tk_kappa program, " << progname << endl;
	this->print_info_global(os);
	*os  << "## **********************************************************" << endl;
	*os  << "## Simulation parameters:" << endl;
	*os  << "## beta           = " << this->get_beta() << endl;
	*os  << "## thermsteps     = " << this->get_thermalizationsteps() << endl;
	*os  << "## heatbathsteps  = " << this->get_heatbathsteps() << endl;
	*os  << "## overrelaxsteps = " << this->get_overrelaxsteps() << endl;
	if(use_autotuning == true)
		*os << "## autotuning for hybrid program is switched on" << endl;
	*os  << "## TODO: INSERT SPECIFIC PARAMETERS!!!!!" << endl;
	*os  << "## **********************************************************" << endl;
	return;
}

void inputparameters::print_info_tkkappa(char* progname) const
{
	logger.info() << "## Starting tk_kappa program, " << progname ;
	this->print_info_global();
	logger.info() << "## **********************************************************";
	logger.info() << "## Simulation parameters:";
	logger.info() << "## beta           = " << this->get_beta();
	logger.info() << "## thermsteps     = " << this->get_thermalizationsteps() ;
	logger.info() << "## heatbathsteps  = " << this->get_heatbathsteps();
	logger.info() << "## overrelaxsteps = " << this->get_overrelaxsteps();
	if(use_autotuning == true)
		logger.info() << "## autotuning for hybrid program is switched on";
	logger.info() << "## TODO: INSERT SPECIFIC PARAMETERS!!!!!";
	logger.info() << "## **********************************************************";
	return;
}


void inputparameters::print_info_fermion() const
{
	logger.info() << "## **********************************************************";
	logger.info() << "## Fermionic parameters:";
	logger.info() << "##" ;
	logger.info() << "## Boundary Conditions:";
	logger.info() << "## theta_fermion_spatial  = " << this->get_theta_fermion_spatial();
	logger.info() << "## theta_fermion_temporal = " << this->get_theta_fermion_temporal();
	logger.info() << "##" ;
	logger.info() << "## Chemical Potential:" ;
	if(this->get_use_chem_pot_re() == true)
		logger.info() << "## chem_pot_re  = " << this->get_chem_pot_re();
	else
		logger.info() << "## do not use real chem. pot.";
	if(this->get_use_chem_pot_im() == true)
		logger.info() << "## chem_pot_im = " << this->get_chem_pot_im();
	else
		logger.info() << "## do not use imag. chem. pot.";
	logger.info() << "##" ;
	if(this->get_fermact() == WILSON) {
		logger.info() <<  "## fermion action: unimproved Wilson";
		logger.info() << "## kappa  = " << this->get_kappa();
	}
	if(this->get_fermact() == TWISTEDMASS) {
		logger.info() <<  "## fermion action: twisted mass Wilson";
		logger.info() << "## kappa  = " << this->get_kappa();
		logger.info() << "## mu     = " << this->get_mu();
	}
	if(this->get_fermact() == CLOVER) {
		logger.info() <<  "## fermion action: clover Wilson";
		logger.info() << "## kappa  = " << this->get_kappa();
		logger.info() << "## csw    = " << this->get_csw();
	}
	logger.info() << "##" ;
	logger.info() << "## Inverter parameters:";
	if(this->get_use_pointsource() == true) {
		logger.info() << "## Use pointsource for inversion" ;
		logger.info() << "## Position (x,y,z,t): " << pointsource_x << " " <<  pointsource_y << " " <<  pointsource_z << " " <<  pointsource_t;
	}
	if(this->get_use_pointsource() == false) {
		logger.info() << "## Use stochastic sources for inversion" ;
		logger.info() << "## Number of sources: " << this->get_num_sources();
	}
	if(this->get_use_cg() == true)
		logger.info() << "## Use CG-solver for inversions" ;
	if(this->get_use_cg() == false)
		logger.info() << "## Use BiCGStab for inversions";
	if(this->get_use_eo() == true)
		logger.info() << "## Use even-odd preconditioning" ;
	if(this->get_use_eo() == false)
		logger.info() << "## Do NOT use even-odd preconditioning";
	logger.info() << "## cgmax  = " << this->get_cgmax();

	//print extra warning if BC are set to default since this is a serious source of errors...
	if ( this->get_theta_fermion_spatial() == 0. && this->get_theta_fermion_temporal() == 0.) {
		logger.warn() << "\nNOTE: BCs have been set to periodic values by default!!\nTo change this use e.g. ThetaT/ThetaS in the input-file.\n";
	}
}

void inputparameters::print_info_fermion(ostream * os) const
{
	*os  << "## **********************************************************" << endl;
	*os  << "## Fermionic parameters:" << endl;
	*os  << "##" << endl;
	*os  << "## Boundary Conditions:" << endl;
	*os  << "## theta_fermion_spatial  = " << this->get_theta_fermion_spatial() << endl;
	*os  << "## theta_fermion_temporal = " << this->get_theta_fermion_temporal() << endl;

	*os  << "##" << endl;
	*os  << "## Chemical Potential:" << endl;
	if(this->get_use_chem_pot_re() == 1)
		*os  << "## chem_pot_re  = " << this->get_chem_pot_re() << endl;
	else
		*os  << "## do not use real chem. pot." << endl;
	if(this->get_use_chem_pot_im() == 1)
		*os  << "## chem_pot_im = " << this->get_chem_pot_im() << endl;
	else
		*os  << "## do not use imag. chem. pot." << endl;
	*os  << "##" << endl;
	if(this->get_fermact() == WILSON) {
		*os <<  "## fermion action: unimproved Wilson" << endl;
		*os  << "## kappa  = " << this->get_kappa() << endl;
	}
	if(this->get_fermact() == TWISTEDMASS) {
		*os <<  "## fermion action: twisted mass Wilson" << endl;
		*os  << "## kappa  = " << this->get_kappa() << endl;
		*os  << "## mu     = " << this->get_mu() << endl;
	}
	if(this->get_fermact() == CLOVER) {
		*os <<  "## fermion action: clover Wilson" << endl;
		*os  << "## kappa  = " << this->get_kappa() << endl;
		*os  << "## csw    = " << this->get_csw() << endl;
	}
	*os  << "##" << endl;
	*os  << "## Inverter parameters:" << endl;
	if(this->get_use_pointsource() == true) {
		*os  << "## Use pointsource for inversion"  << endl;
		*os  << "## Position (x,y,z,t): " << pointsource_x << " " <<  pointsource_y << " " <<  pointsource_z << " " <<  pointsource_t  << endl;
	}
	if(this->get_use_pointsource() == false) {
		*os  << "## Use stochastic sources for inversion"  << endl;
		*os  << "## Number of sources: " << this->get_num_sources()  << endl;
	}
	if(this->get_use_cg() == true)
		*os << "## Use CG-solver for inversions"  << endl;
	if(this->get_use_cg() == false)
		*os << "## Use BiCGStab for inversions" << endl;

	if(this->get_use_eo() == true)
		*os  << "## Use even-odd preconditioning" << endl;
	if(this->get_use_eo() == false)
		*os  << "## Do NOT use even-odd preconditioning" << endl;
	*os << "## cgmax  = " << this->get_cgmax() << endl;

	//print extra warning if BC are set to default since this is a serious source of errors...
	if ( this->get_theta_fermion_spatial() == 0. && this->get_theta_fermion_temporal() == 0.) {
		*os << "\nNOTE: BCs have been set to periodic values by default!!\nTo change this use e.g. ThetaT/ThetaS in the input-file.\n";
	}
}


void inputparameters::print_info_inverter(char* progname) const
{
	logger.info() << "## Starting inverter program, executable name: " << progname;
	this->print_info_global();
	this->print_info_fermion();
	logger.info() << "## **********************************************************";
	return;
}

void inputparameters::print_info_inverter(char* progname, ostream* os) const
{
	*os << "## Starting inverter program, executable name: " << progname << endl;
	this->print_info_global(os);
	this->print_info_fermion(os);
	*os << "## **********************************************************" << endl;
	return;
}

void inputparameters::print_info_hmc(char* progname) const
{

	logger.info() << "## Starting hmc program, executable name: " << progname ;
	this->print_info_global();
	this->print_info_fermion();
	logger.info() << "##  ";
	logger.info() << "## HMC parameters: " ;
	logger.info() << "## tau  = " << this->get_tau();
	logger.info() << "## HMC steps  = " << this->get_hmcsteps();
	logger.info() << "## # Timescales  = " << this->get_num_timescales();
	logger.info() << "## integrationsteps1  = " << this->get_integrationsteps1();
	if(this->get_num_timescales() == 2)
		logger.info() << "## integrationsteps2  = " << this->get_integrationsteps2();
	logger.info() << "## **********************************************************";
	return;
}

void inputparameters::print_info_hmc(char* progname, ostream* os) const
{
	*os << "## Starting hmc program, executable name: " << progname << endl;
	this->print_info_global(os);
	this->print_info_fermion(os);
	*os << "##  " << '\n';
	*os << "## HMC parameters: "  << '\n';
	*os << "## tau  = " << this->get_tau() << '\n';
	*os << "## HMC steps  = " << this->get_hmcsteps() << '\n';
	*os << "## # Timescales  = " << this->get_num_timescales() << '\n';
	*os << "## integrationsteps1  = " << this->get_integrationsteps1() << '\n';
	if(this->get_num_timescales() == 2)
		*os << "## integrationsteps2  = " << this->get_integrationsteps2() << '\n';
	*os << "## **********************************************************\n";
	*os << std::endl;
	return;
}

void inputparameters::set_settings_global()
{
	logger.trace() << "setting global parameters...";
	this->volspace = this->nspace * this->nspace * this->nspace;
	this->vol4d = this->volspace * this->ntime;
	this->spinorfieldsize = this->vol4d;
	this->eoprec_spinorfieldsize = this->spinorfieldsize / 2;
	this->gaugemomentasize = NDIM * this->vol4d;
	if(get_use_rec12() == true)
		this->gaugefieldsize = NC * (NC - 1) * NDIM * this->vol4d;
	else
		this->gaugefieldsize = NC * NC * NDIM * this->vol4d;

	//set sizes of buffers
	this->sf_buf_size = this->spinorfieldsize * sizeof(spinor);
	this->eo_sf_buf_size = this->eoprec_spinorfieldsize * sizeof(spinor);
	this->gf_buf_size = this->gaugefieldsize * sizeof(hmc_complex);
	this->gm_buf_size = this->gaugemomentasize * sizeof(ae);
}

void inputparameters::check_settings_global() const
{
	logger.trace() << "checking compile-time-parameters against input-parameters...";

	try {

		//compile time parameters
		//numerical precision
#ifdef _USEDOUBLEPREC_
		if( this->get_prec() != 64) throw Invalid_Parameters("Numerical precision", "64", this->get_prec());
#else
		if( this->get_prec() != 32) throw Invalid_Parameters("Numerical precision", "32", this->get_prec());
#endif

		if( this->get_use_rec12() == true) throw Invalid_Parameters("Reconstruct12.", "OFF", "ON");

	}//try
	catch (Invalid_Parameters& e) {
		logger.fatal() << "Error in check_setting_global():";
		logger.fatal() << e.what();
		logger.fatal() << "Aborting.";
		exit(EXIT_INPUTPARAMETERS);
	}
}
