#include "inputparameters.h"

#include "logger.hpp"

hmc_error inputparameters::set_defaults()
{
	kappa = 0.125;
	beta = 4.0;
	mu = 0.006;
	csw = 0.;
	cgmax = 1000;
	prec = 64;
	theta_fermion_spatial = 0.;
	theta_fermion_temporal = 0.;	
	theta_gaugefield = 0.;
	chem_pot_re = 0.;
	chem_pot_im = 0.;
	tau = 0.5;
	integrationsteps1 = 10;
	integrationsteps2 = 10;
	startcondition = COLD_START;
	thermalizationsteps = 0;
	heatbathsteps = 1000;
	hmcsteps = 10;
	overrelaxsteps = 1;
	writefrequency = 1;
	savefrequency = 100;
	saveconfigs = FALSE;
	use_eo = TRUE;
	//sourcefile = "\0";
	sourcefilenumber = "00000";
	fermact = WILSON;
	num_dev = 1;
	perform_heatbath = 1;
	return HMC_SUCCESS;
}

hmc_error inputparameters::readfile(char* ifn)
{
	std::ifstream infile;
	infile.open(ifn);
	if(!infile.is_open()) {
		printf("Could not open input file: %s\n", ifn);
		exit(HMC_FILEERROR);
	}

	int muset = FALSE;
	int cswset = FALSE;

	while (infile.good()) {
		std::string line;
		infile >> line;
		if(line.find("#") != std::string::npos) continue; //allow comments
		if(line.find("kappa") != std::string::npos) val_assign(&kappa, line);
		if(line.find("Kappa") != std::string::npos) val_assign(&kappa, line);
		if(line.find("mu") != std::string::npos) {
			val_assign(&mu, line);
			muset = TRUE;
		}
		if(line.find("Mu") != std::string::npos) {
			val_assign(&mu, line);
			muset = TRUE;
		}
		if(line.find("csw") != std::string::npos) {
			val_assign(&csw, line);
			cswset = TRUE;
		}
		if(line.find("Csw") != std::string::npos) {
			val_assign(&csw, line);
			cswset = TRUE;
		}
		if(line.find("beta") != std::string::npos) val_assign(&beta, line);
		if(line.find("tau") != std::string::npos) val_assign(&tau, line);
		if(line.find("Beta") != std::string::npos) val_assign(&beta, line);
		
		if(line.find("theta_fermion_spatial") != std::string::npos) val_assign(&theta_fermion_spatial, line);
		if(line.find("theta_fermion_temporal") != std::string::npos) val_assign(&theta_fermion_temporal, line);
		
		
		if(line.find("cgmax") != std::string::npos) val_assign(&cgmax, line);
		if(line.find("CGmax") != std::string::npos) val_assign(&cgmax, line);
		if(line.find("Cgmax") != std::string::npos) val_assign(&cgmax, line);
		if(line.find("writefrequency") != std::string::npos) val_assign(&writefrequency, line);
		if(line.find("savefrequency") != std::string::npos) val_assign(&savefrequency, line);
		if(line.find("saveconfigs") != std::string::npos) savecond_assign(&saveconfigs, line);
		if(line.find("prec") != std::string::npos) val_assign(&prec, line);
		if(line.find("Prec") != std::string::npos) val_assign(&prec, line);
		if(line.find("readsource") != std::string::npos) cond_assign(&startcondition, line);
		if(line.find("startcondition") != std::string::npos) cond_assign(&startcondition, line);
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
		if(line.find("num_dev") != std::string::npos) val_assign(&num_dev, line);
		if(line.find("perform_heatbath") != std::string::npos) val_assign(&perform_heatbath, line);

		if(line.find("fermaction") != std::string::npos) fermact_assign(&fermact, line);
		if(line.find("fermionaction") != std::string::npos) fermact_assign(&fermact, line);
		if(line.find("fermact") != std::string::npos) fermact_assign(&fermact, line);

		if(line.find("evenodd") != std::string::npos) eocond_assign(&use_eo, line);
		if(line.find("even_odd") != std::string::npos) eocond_assign(&use_eo, line);
		if(line.find("even-odd") != std::string::npos) eocond_assign(&use_eo, line);
		if(line.find("even-odd-preconditioning") != std::string::npos) eocond_assign(&use_eo, line);
		if(line.find("use_eo") != std::string::npos) eocond_assign(&use_eo, line);
		if(line.find("use_evenodd") != std::string::npos) eocond_assign(&use_eo, line);
		
		

	}

	if(muset == TRUE && fermact != TWISTEDMASS) {
		logger.fatal() << "Setting a value for mu is not allowed for fermion action other than twisted mass. Aborting...";
		exit(HMC_STDERR);
	}
	if(cswset == TRUE && fermact != CLOVER) {
		logger.fatal() << "Setting a value for csw is not allowed for fermion action other than clover. Aborting...";
		exit(HMC_STDERR);
	}
	if(cswset == TRUE && muset == TRUE) {
		logger.fatal() << "Setting values for both csw and mu is currently not allowed. Aborting...";
		exit(HMC_STDERR);
	}

	return HMC_SUCCESS;
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

void inputparameters::cond_assign(int * out, std::string line)
{
	if(std::strstr(line.c_str(), "cold") != NULL) {
		(*out) = COLD_START;
		return;
	}
	if(std::strstr(line.c_str(), "hot") != NULL) {
		(*out) = HOT_START;
		return;
	}
	if(std::strstr(line.c_str(), "source") != NULL) {
		(*out) = START_FROM_SOURCE;
		return;
	}
	if(std::strstr(line.c_str(), "continue") != NULL) {
		(*out) = START_FROM_SOURCE;
		return;
	}
	printf("invalid startcondition\n");
	exit(HMC_STDERR);
	return;
}

void inputparameters::fermact_assign(int * out, std::string line)
{
	if(std::strstr(line.c_str(), "TWISTEDMASS") != NULL) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(), "twistedmass") != NULL) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(), "Twistedmass") != NULL) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(), "TwistedMass") != NULL) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(), "clover") != NULL) {
		(*out) = CLOVER;
		return;
	}
	if(std::strstr(line.c_str(), "CLOVER") != NULL) {
		(*out) = CLOVER;
		return;
	}
	if(std::strstr(line.c_str(), "Clover") != NULL) {
		(*out) = CLOVER;
		return;
	}
	if(std::strstr(line.c_str(), "WILSON") != NULL) {
		(*out) = WILSON;
		return;
	}
	if(std::strstr(line.c_str(), "Wilson") != NULL) {
		(*out) = WILSON;
		return;
	}
	if(std::strstr(line.c_str(), "wilson") != NULL) {
		(*out) = WILSON;
		return;
	}
	if(std::strstr(line.c_str(), "unimproved") != NULL) {
		(*out) = WILSON;
		return;
	}
	printf("invalid fermion action\n");
	exit(HMC_STDERR);
	return;
}


void inputparameters::savecond_assign(int * out, std::string line)
{
	if(std::strstr(line.c_str(), "yes") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "true") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "TRUE") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "True") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "no") != NULL) {
		(*out) = FALSE;
		return;
	}
	if(std::strstr(line.c_str(), "false") != NULL) {
		(*out) = FALSE;
		return;
	}
	if(std::strstr(line.c_str(), "FALSE") != NULL) {
		(*out) = FALSE;
		return;
	}
	if(std::strstr(line.c_str(), "False") != NULL) {
		(*out) = FALSE;
		return;
	}
	printf("invalid save condition\n");
	exit(HMC_STDERR);
	return;
}

void inputparameters::eocond_assign(int * out, std::string line)
{
	if(std::strstr(line.c_str(), "yes") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "true") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "TRUE") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "True") != NULL) {
		(*out) = TRUE;
		return;
	}
	if(std::strstr(line.c_str(), "no") != NULL) {
		(*out) = FALSE;
		return;
	}
	if(std::strstr(line.c_str(), "false") != NULL) {
		(*out) = FALSE;
		return;
	}
	if(std::strstr(line.c_str(), "FALSE") != NULL) {
		(*out) = FALSE;
		return;
	}
	if(std::strstr(line.c_str(), "False") != NULL) {
		(*out) = FALSE;
		return;
	}
	printf("invalid even-odd condition\n");
	exit(HMC_STDERR);
	return;
}


void inputparameters::val_assign(std::string * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);
	(*out) = value.c_str();
	return;
}

hmc_float inputparameters::get_kappa()
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

hmc_float inputparameters::get_mubar()
{
	return mubar;
}

hmc_float inputparameters::get_theta_fermion_spatial()
{
	return theta_fermion_spatial;
}

hmc_float inputparameters::get_theta_fermion_temporal()
{
	return theta_fermion_temporal;
}

hmc_float inputparameters::get_theta_gaugefield()
{
	return theta_gaugefield;
}

hmc_float inputparameters::get_beta()
{
	return beta;
}

hmc_float inputparameters::get_mu()
{
	return mu;
}

hmc_float inputparameters::get_tau()
{
	return tau;
}

hmc_float inputparameters::get_csw()
{
	return csw;
}

hmc_float inputparameters::get_chem_pot_re()
{
	return chem_pot_re;
}

hmc_float inputparameters::get_chem_pot_im()
{
	return chem_pot_im;
}

int inputparameters::get_cgmax()
{
	return cgmax;
}

int inputparameters::get_prec()
{
	return prec;
}

int inputparameters::get_num_dev()
{
	return num_dev;
}

int inputparameters::get_thermalizationsteps()
{
	return thermalizationsteps;
}

int inputparameters::get_heatbathsteps()
{
	return heatbathsteps;
}

int inputparameters::get_overrelaxsteps()
{
	return overrelaxsteps;
}

int inputparameters::get_hmcsteps()
{
	return hmcsteps;
}

int inputparameters::get_integrationsteps1()
{
	return integrationsteps1;
}

int inputparameters::get_integrationsteps2()
{
	return integrationsteps2;
}


int inputparameters::get_writefrequency()
{
	return writefrequency;
}

int inputparameters::get_savefrequency()
{
	return savefrequency;
}

int inputparameters::get_startcondition()
{
	return startcondition;
}

int inputparameters::get_fermact()
{
	return fermact;
}

int inputparameters::get_perform_heatbath()
{
	return perform_heatbath;
}

int inputparameters::get_saveconfigs()
{
	return saveconfigs;
}

void inputparameters::display_sourcefile()
{
	cout << sourcefile;
}

void inputparameters::display_sourcefilenumber()
{
	cout << sourcefilenumber;
}

int inputparameters::get_use_eo()
{
	return use_eo;
}

void inputparameters::print_info_heatbath(char* progname){
  logger.info() << "## Starting heatbath program, executable name: " << progname;
  logger.info() << "## **********************************************************";
  logger.info() << "## Compile time parameters:";
  logger.info() << "## NSPACE:  " << NSPACE;
  logger.info() << "## NTIME:   " << NTIME;
  logger.info() << "## NDIM:    " << NDIM;
  logger.info() << "## NCOLOR:  " << NC;
  logger.info() << "## NSPIN:   " << NSPIN;
  logger.info() << "##";
  logger.info() << "## Run time parameters:";
  logger.info() << "## beta  = " << this->get_beta();
  logger.info() << "## prec  = " << this->get_prec();
  logger.info() << "## thermsteps     = " << this->get_thermalizationsteps() ;
  logger.info() << "## heatbathsteps  = " << this->get_heatbathsteps();
  logger.info() << "## overrelaxsteps = " << this->get_overrelaxsteps();
  logger.info() << "##";
  
  logger.info() << "## number of devices demanded for calculations: " << this->get_num_dev()  ;
  logger.info() << "##" ;	
	
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
  logger.info() << "## **********************************************************";
  return;
}

void inputparameters::print_info_heatbath(char* progname, ostream* os){
  *os << "## Starting heatbath program, executable name: " << progname <<endl;
  *os << "## **********************************************************"<<endl;
  *os << "## Compile time parameters:"<<endl;
  *os << "## NSPACE:  " << NSPACE<<endl;
  *os << "## NTIME:   " << NTIME<<endl;
  *os << "## NDIM:    " << NDIM<<endl;
  *os << "## NCOLOR:  " << NC<<endl;
  *os << "## NSPIN:   " << NSPIN<<endl;
  *os << "##"<<endl;
  *os << "## Run time parameters:"<<endl;
  *os << "## beta  = " << this->get_beta()<<endl;
  *os << "## prec  = " << this->get_prec()<<endl;
  *os << "## thermsteps     = " << this->get_thermalizationsteps()<<endl ;
  *os << "## heatbathsteps  = " << this->get_heatbathsteps()<<endl;
  *os << "## overrelaxsteps = " << this->get_overrelaxsteps()<<endl;
  *os << "##"<<endl;
  
  *os << "## number of devices demanded for calculations: " << this->get_num_dev()<<endl  ;
  *os << "##" <<endl;	
	
  if (this->get_startcondition() == START_FROM_SOURCE) {
    string sf = this->sourcefile;
    *os << "## sourcefile = " << sf<<endl;
  }
  if (this->get_startcondition() == COLD_START) {
    *os << "## cold start"<<endl;
  }
  if (this->get_startcondition() == HOT_START) {
    *os << "## hot start"<<endl;
  }
  *os << "## **********************************************************"<<endl;
  return;
}

void inputparameters::print_info_inverter(char* progname){
  logger.info() << "## Starting inverter program, executable name: " << progname;
  logger.info() << "## **********************************************************";
  logger.info() << "## Compile time parameters:";
  logger.info() << "## NSPACE:  " << NSPACE;
  logger.info() << "## NTIME:   " << NTIME;
  logger.info() << "## NDIM:    " << NDIM;
  logger.info() << "## NCOLOR:  " << NC;
  logger.info() << "## NSPIN:   " << NSPIN;
  logger.info() << "##" ;
  logger.info() << "## Run time parameters:";

	if(this->get_fermact()==WILSON) {
	  logger.info()<<  "## fermion action: unimproved Wilson";
	  logger.info() << "## kappa  = "<<this->get_kappa();
	}
	if(this->get_fermact()==TWISTEDMASS) {
	  logger.info()<<  "## fermion action: twisted mass Wilson";
	  logger.info() << "## kappa  = "<<this->get_kappa();
	  logger.info() << "## mu     = "<<this->get_mu();
	}
	if(this->get_fermact()==CLOVER) {
	  logger.info()<<  "## fermion action: clover Wilson";
	  logger.info() << "## kappa  = "<<this->get_kappa();
	  logger.info() << "## csw    = "<<this->get_csw();
	}
	logger.info() << "## prec  = " << this->get_prec() ;
	logger.info() << "##" ;
	if(this->get_use_eo()==TRUE)
	  logger.info() << "## use even-odd preconditioning" ;
	if(this->get_use_eo()==FALSE) 
	  logger.info() << "## do not use even-odd preconditioning";
	logger.info() << "##" ;

	logger.info() << "## Boundary Conditions:";
	logger.info() << "## theta_fermion_spatial  = "<<this->get_theta_fermion_spatial();
	logger.info() << "## theta_fermion_temporal = "<<this->get_theta_fermion_temporal();
	logger.info() << "##" ;

	logger.info() << "## Chemical Potential:" ;
	#ifdef _CP_REAL_
	logger.info() << "## chem_pot_re  = "<<this->get_chem_pot_re();
	#else
	logger.info() << "## do not use real chem. pot.";
	#endif
	#ifdef _CP_IMAG_
	logger.info() << "## chem_pot_im = "<<this->get_chem_pot_im();
	#else
	logger.info() << "## do not use imag. chem. pot.";
	#endif
	logger.info() << "##" ;
	
	logger.info()<<"## number of devices desired for calculations: " << this->get_num_dev() ;
	logger.info() << "##" ;
	
	if (this->get_startcondition() == START_FROM_SOURCE) {
		string sf = this->sourcefile;
		logger.info() << "## sourcefile = " << sf ;
	}
	if (this->get_startcondition() == COLD_START) {
		logger.info() << "## WARNING: cold start - no configuration read";
	}
	if (this->get_startcondition() == HOT_START) {
		logger.info() << "## WARNING: hot start - no configuration read";
	}
	logger.info() << "## **********************************************************";
  return;
}

void inputparameters::print_info_inverter(char* progname, ostream* os){
  *os << "## Starting inverter program, executable name: " << progname << endl;
  *os << "## **********************************************************"<< endl;
  *os << "## Compile time parameters:"<< endl;
  *os << "## NSPACE:  " << NSPACE<< endl;
  *os << "## NTIME:   " << NTIME<< endl;
  *os << "## NDIM:    " << NDIM<< endl;
  *os << "## NCOLOR:  " << NC<< endl;
  *os << "## NSPIN:   " << NSPIN<< endl;
  *os << "##" << endl;
  *os << "## Run time parameters:"<< endl;

	if(this->get_fermact()==WILSON) {
	  *os<<  "## fermion action: unimproved Wilson"<< endl;
	  *os << "## kappa  = "<<this->get_kappa()<< endl;
	}
	if(this->get_fermact()==TWISTEDMASS) {
	  *os<<  "## fermion action: twisted mass Wilson"<< endl;
	  *os << "## kappa  = "<<this->get_kappa()<< endl;
	  *os << "## mu     = "<<this->get_mu()<< endl;
	}
	if(this->get_fermact()==CLOVER) {
	  *os<<  "## fermion action: clover Wilson"<< endl;
	  *os << "## kappa  = "<<this->get_kappa()<< endl;
	  *os << "## csw    = "<<this->get_csw()<< endl;
	}
	*os << "## prec  = " << this->get_prec() << endl;
	*os << "##" << endl;
	if(this->get_use_eo()==TRUE)
	  *os << "## use even-odd preconditioning"<< endl ;
	if(this->get_use_eo()==FALSE) 
	  *os << "## do not use even-odd preconditioning"<< endl;
	*os << "##" << endl;

	*os << "## Boundary Conditions:"<< endl;
	*os << "## theta_fermion_spatial  = "<<this->get_theta_fermion_spatial()<< endl;
	*os << "## theta_fermion_temporal = "<<this->get_theta_fermion_temporal()<< endl;
	*os << "##" << endl;

	*os << "## Chemical Potential:" << endl;
	#ifdef _CP_REAL_
	*os << "## chem_pot_re  = "<<this->get_chem_pot_re()<< endl;
	#else
	*os << "## do not use real chem. pot."<< endl;
	#endif
	#ifdef _CP_IMAG_
	*os << "## chem_pot_im = "<<this->get_chem_pot_im()<< endl;
	#else
	*os << "## do not use imag. chem. pot."<< endl;
	#endif
	*os << "##" << endl;
	
	*os<<"## number of devices desired for calculations: " << this->get_num_dev() << endl;
	*os << "##" << endl;
	
	if (this->get_startcondition() == START_FROM_SOURCE) {
		string sf = this->sourcefile;
		*os << "## sourcefile = " << sf << endl;
	}
	if (this->get_startcondition() == COLD_START) {
		*os << "## WARNING: cold start - no configuration read"<< endl;
	}
	if (this->get_startcondition() == HOT_START) {
		*os << "## WARNING: hot start - no configuration read"<< endl;
	}
	*os << "## **********************************************************"<< endl;
  return;
}

void inputparameters::print_info_tkkappa(char* progname, ostream* os){
	*os << "## Starting tk_kappa program, " << progname << endl;
	*os << "## **********************************************************\n";
	*os << "## Compile time parameters:\n";
	*os << "## NSPACE:  " << NSPACE << '\n';
	*os << "## NTIME:   " << NTIME << '\n';
	*os << "## NDIM:    " << NDIM << '\n';
	*os << "## NCOLOR:  " << NC << '\n';
	*os << "## NSPIN:   " << NSPIN << '\n';
	*os << "##" << '\n';
	*os << "## Run time parameters:\n";
	*os << "## beta  = " << this->get_beta() << '\n';
	*os << "## prec  = " << this->get_prec() << '\n';
	*os << "## thermsteps     = " << this->get_thermalizationsteps() << '\n';
	*os << "## heatbathsteps  = " << this->get_heatbathsteps() << '\n';
	*os << "## overrelaxsteps = " << this->get_overrelaxsteps() << '\n';
	*os << "##" << '\n';
	if (this->get_startcondition() == START_FROM_SOURCE) {
		*os << "## sourcefile = ";
		string sf=this->sourcefile;
		*os << sf << '\n';
	}
	if (this->get_startcondition() == COLD_START) {
		*os << "## cold start\n";
	}
	if (this->get_startcondition() == HOT_START) {
		*os << "## hot start\n";
	}
	*os << "## **********************************************************\n";
	*os << std::endl;
  return;
}

void inputparameters::print_info_tkkappa(char* progname){
	logger.info() << "## Starting tk_kappa program, " << progname ;
	logger.info() << "## **********************************************************";
	logger.info() << "## Compile time parameters:";
	logger.info() << "## NSPACE:  " << NSPACE;
	logger.info() << "## NTIME:   " << NTIME;
	logger.info() << "## NDIM:    " << NDIM;
	logger.info() << "## NCOLOR:  " << NC;
	logger.info() << "## NSPIN:   " << NSPIN;
	logger.info() << "##";
	logger.info() << "## Run time parameters:";
	logger.info() << "## beta  = " << this->get_beta();
	logger.info() << "## prec  = " << this->get_prec();
	logger.info() << "## thermsteps     = " << this->get_thermalizationsteps();
	logger.info() << "## heatbathsteps  = " << this->get_heatbathsteps();
	logger.info() << "## overrelaxsteps = " << this->get_overrelaxsteps();
	logger.info() << "##";
	if (this->get_startcondition() == START_FROM_SOURCE) {
	  string sf=this->sourcefile;
	  logger.info() << "## sourcefile = "<< sf;
	}
	if (this->get_startcondition() == COLD_START) {
		logger.info() << "## cold start";
	}
	if (this->get_startcondition() == HOT_START) {
		logger.info() << "## hot start";
	}
	logger.info() << "## **********************************************************";
  return;
}

void inputparameters::print_info_hmc(char* progname){

	logger.info() << "## Starting hmc program, executable name: " << progname ;
	logger.info() << "## **********************************************************";
	logger.info() << "## Compile time parameters:";
	logger.info() << "## NSPACE:  " << NSPACE ;
	logger.info() << "## NTIME:   " << NTIME ;
	logger.info() << "## NDIM:    " << NDIM ;
	logger.info() << "## NCOLOR:  " << NC ;
	logger.info() << "## NSPIN:   " << NSPIN ;
	logger.info() << "##";
	logger.info() << "## Run time parameters:";

	if(this->get_fermact()==WILSON) {
	  logger.info()<<  "## fermion action: unimproved Wilson";
	  logger.info() << "## kappa  = "<<this->get_kappa();
	}
	if(this->get_fermact()==TWISTEDMASS) {
	  logger.info()<<  "## fermion action: twisted mass Wilson";
	  logger.info() << "## kappa  = "<<this->get_kappa();
	  logger.info() << "## mu     = "<<this->get_mu();
	}
	if(this->get_fermact()==CLOVER) {
	  logger.info()<<  "## fermion action: clover Wilson";
	  logger.info() << "## kappa  = "<<this->get_kappa();
	  logger.info() << "## csw    = "<<this->get_csw();
	}
	logger.info() << "## prec  = " << this->get_prec();
	logger.info() << "##";
	if(this->get_use_eo()==TRUE)
	  logger.info() << "## use even-odd preconditioning";
	if(this->get_use_eo()==FALSE) 
	  logger.info() << "## do not use even-odd preconditioning";
	logger.info() << "##";

	logger.info() << "## Boundary Conditions:" ;
	logger.info() << "## theta_fermion_spatial  = "<<this->get_theta_fermion_spatial();
	logger.info() << "## theta_fermion_temporal = "<<this->get_theta_fermion_temporal();
	logger.info() << "##";

	logger.info() << "## Chemical Potential:";
	#ifdef _CP_REAL_
	logger.info() << "## chem_pot_re  = "<<this->get_chem_pot_re();
	#else
	logger.info() << "## do not use real chem. pot.";
	#endif
	#ifdef _CP_IMAG_
	logger.info() << "## chem_pot_im = "<<this->get_chem_pot_im();
	#else
	logger.info() << "## do not use imag. chem. pot.";
	#endif
	logger.info() << "##" ;	
	
	logger.info()<<"## number of devices desired for calculations: " << this->get_num_dev() ;
	logger.info() << "##";
	
	if (this->get_startcondition() == START_FROM_SOURCE) {
		string sf = this->sourcefile;
		logger.info() << "## sourcefile = "<< sf;
	}
	if (this->get_startcondition() == COLD_START) {
		logger.info() << "## WARNING: cold start - no configuration read";
	}
	if (this->get_startcondition() == HOT_START) {
		logger.info() << "## WARNING: hot start - no configuration read";
	}
	logger.info() << "##  ";
	logger.info() << "## HMC parameters: " ;
	logger.info() << "## tau  = " << this->get_tau();
	logger.info() << "## HMC steps  = " << this->get_hmcsteps();
	logger.info() << "## integrationsteps1  = " << this->get_integrationsteps1();
	logger.info() << "## integrationsteps2  = " << this->get_integrationsteps2();
	logger.info() << "## **********************************************************";
    return;
 }

 void inputparameters::print_info_hmc(char* progname, ostream* os){
   *os << "## Starting hmc program, executable name: " << progname << endl;
	*os << "## **********************************************************\n";
	*os << "## Compile time parameters:\n";
	*os << "## NSPACE:  " << NSPACE << '\n';
	*os << "## NTIME:   " << NTIME << '\n';
	*os << "## NDIM:    " << NDIM << '\n';
	*os << "## NCOLOR:  " << NC << '\n';
	*os << "## NSPIN:   " << NSPIN << '\n';
	*os << "##" << '\n';
	*os << "## Run time parameters:\n";

	if(this->get_fermact()==WILSON) {
	  *os<<  "## fermion action: unimproved Wilson"<<'\n';
	  *os << "## kappa  = "<<this->get_kappa()<< '\n';
	}
	if(this->get_fermact()==TWISTEDMASS) {
	  *os<<  "## fermion action: twisted mass Wilson"<<'\n';
	  *os << "## kappa  = "<<this->get_kappa()<< '\n';
	  *os << "## mu     = "<<this->get_mu()<< '\n';
	}
	if(this->get_fermact()==CLOVER) {
	  *os<<  "## fermion action: clover Wilson"<<'\n';
	  *os << "## kappa  = "<<this->get_kappa()<< '\n';
	  *os << "## csw    = "<<this->get_csw()<< '\n';
	}
	*os << "## prec  = " << this->get_prec() << '\n';
	*os << "##" << '\n';
	if(this->get_use_eo()==TRUE)
	  *os << "## use even-odd preconditioning" << '\n';
	if(this->get_use_eo()==FALSE) 
	  *os << "## do not use even-odd preconditioning" << '\n';
	*os << "##" << '\n';

	*os << "## Boundary Conditions:" << '\n';
	*os << "## theta_fermion_spatial  = "<<this->get_theta_fermion_spatial()<< '\n';
	*os << "## theta_fermion_temporal = "<<this->get_theta_fermion_temporal()<< '\n';
	*os << "##" << '\n';

	*os << "## Chemical Potential:" << '\n';
	#ifdef _CP_REAL_
	*os << "## chem_pot_re  = "<<this->get_chem_pot_re()<< '\n';
	#else
	*os << "## do not use real chem. pot."<< '\n';
	#endif
	#ifdef _CP_IMAG_
	*os << "## chem_pot_im = "<<this->get_chem_pot_im()<< '\n';
	#else
	*os << "## do not use imag. chem. pot."<< '\n';
	#endif
	*os << "##" << '\n';	
	
	*os<<"## number of devices desired for calculations: " << this->get_num_dev() << "\n" ;
	*os << "##" << '\n';
	
	if (this->get_startcondition() == START_FROM_SOURCE) {
		*os << "## sourcefile = ";
		string sf = this->sourcefile;
		*os << sf << '\n';
	}
	if (this->get_startcondition() == COLD_START) {
		*os << "## WARNING: cold start - no configuration read\n";
	}
	if (this->get_startcondition() == HOT_START) {
		*os << "## WARNING: hot start - no configuration read\n";
	}
	*os << "##  " << '\n';
	*os << "## HMC parameters: "  << '\n';
	*os << "## tau  = " << this->get_tau() << '\n';
	*os << "## HMC steps  = " << this->get_hmcsteps() << '\n';
	*os << "## integrationsteps1  = " << this->get_integrationsteps1() << '\n';
	*os << "## integrationsteps2  = " << this->get_integrationsteps2() << '\n';
	*os << "## **********************************************************\n";
	*os << std::endl;
    return;
}
