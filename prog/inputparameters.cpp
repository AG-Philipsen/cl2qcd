#include "inputparameters.h"

hmc_error inputparameters::set_defaults()
{
	kappa = 0.125;
	beta = 4.0;
	mu = 0.006;
	csw = 0.;
	cgmax = 1000;
	prec = 64;
	theta_fermion = 0.;
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
	//sourcefile = "\0";
	sourcefilenumber = "00000";
	fermact = WILSON;
	return HMC_SUCCESS;
}

hmc_error inputparameters::readfile(char* ifn)
{
	std::ifstream infile;
	infile.open(ifn);
	if(!infile.is_open()) {
		printf("Could not open input file: %s\n",ifn);
		exit(HMC_FILEERROR);
	}

	int muset = FALSE;
	int cswset = FALSE;

	while (infile.good()) {
		std::string line;
		infile>>line;
		if(line.find("#")!=std::string::npos) continue; //allow comments
		if(line.find("kappa")!=std::string::npos) val_assign(&kappa,line);
		if(line.find("Kappa")!=std::string::npos) val_assign(&kappa,line);
		if(line.find("mu")!=std::string::npos) {
		  val_assign(&mu,line);
		  muset = TRUE;
		}
		if(line.find("Mu")!=std::string::npos){
		  val_assign(&mu,line);
		  muset = TRUE;
		}
		if(line.find("csw")!=std::string::npos) {
		  val_assign(&csw,line);
		  cswset = TRUE;
		}
		if(line.find("Csw")!=std::string::npos) {
		  val_assign(&csw,line);
		  cswset = TRUE;
		}
		if(line.find("beta")!=std::string::npos) val_assign(&beta,line);
		if(line.find("tau")!=std::string::npos) val_assign(&tau,line);
		if(line.find("Beta")!=std::string::npos) val_assign(&beta,line);
		if(line.find("cgmax")!=std::string::npos) val_assign(&cgmax,line);
		if(line.find("CGmax")!=std::string::npos) val_assign(&cgmax,line);
		if(line.find("Cgmax")!=std::string::npos) val_assign(&cgmax,line);
		if(line.find("writefrequency")!=std::string::npos) val_assign(&writefrequency,line);
		if(line.find("savefrequency")!=std::string::npos) val_assign(&savefrequency,line);
		if(line.find("saveconfigs")!=std::string::npos) savecond_assign(&saveconfigs,line);
		if(line.find("prec")!=std::string::npos) val_assign(&prec,line);
		if(line.find("Prec")!=std::string::npos) val_assign(&prec,line);
		if(line.find("readsource")!=std::string::npos) cond_assign(&startcondition,line);
		if(line.find("startcondition")!=std::string::npos) cond_assign(&startcondition,line);
		if(line.find("sourcefile")!=std::string::npos) {
			val_assign(&sourcefile,line);
			sourcefilenumber_assign(&sourcefilenumber);
		}
		if(line.find("thermalizationsteps")!=std::string::npos) val_assign(&thermalizationsteps,line);
		if(line.find("heatbathsteps")!=std::string::npos) val_assign(&heatbathsteps,line);
		if(line.find("thermsteps")!=std::string::npos) val_assign(&thermalizationsteps,line);
		if(line.find("thermalization")!=std::string::npos) val_assign(&thermalizationsteps,line);
		if(line.find("overrelaxsteps")!=std::string::npos) val_assign(&overrelaxsteps,line);
		if(line.find("overrelax")!=std::string::npos) val_assign(&overrelaxsteps,line);
		if(line.find("oversteps")!=std::string::npos) val_assign(&overrelaxsteps,line);
		if(line.find("hmcsteps")!=std::string::npos) val_assign(&hmcsteps,line);
		if(line.find("integrationtseps1")!=std::string::npos) val_assign(&integrationsteps1,line);
		if(line.find("integrationsteps2")!=std::string::npos) val_assign(&integrationsteps2,line);

		if(line.find("fermaction")!=std::string::npos) fermact_assign(&fermact,line);
		if(line.find("fermionaction")!=std::string::npos) fermact_assign(&fermact,line);
		if(line.find("fermact")!=std::string::npos) fermact_assign(&fermact,line);

	}

	if(muset==TRUE && fermact != TWISTEDMASS) {
	  cout<<"Setting a value for mu is not allowed for fermion action other than twisted mass. Aborting..."<<endl;
	  exit(HMC_STDERR);
	}
	if(cswset==TRUE && fermact != CLOVER) {
	  cout<<"Setting a value for csw is not allowed for fermion action other than clover. Aborting..."<<endl;
	  exit(HMC_STDERR);
	}
	if(cswset==TRUE && muset==TRUE) {
	  cout<<"Setting values for both csw and mu is currently not allowed. Aborting..."<<endl;
	  exit(HMC_STDERR);
	}

	return HMC_SUCCESS;
}

void inputparameters::val_assign(hmc_float * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value=line.substr(pos+1);
	(*out) = atof(value.c_str());
	return;
}

void inputparameters::val_assign(int * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value=line.substr(pos+1);
	(*out) = atoi(value.c_str());
	return;
}

void inputparameters::sourcefilenumber_assign(std::string * out)
{
	//it is supposed that the file is called conf.xxxxx
	size_t length;
	char buffer[20];
	string str ("Test string...");
	length=sourcefile.copy(buffer,5,5);
	buffer[length]='\0';
	//atoi should neglect any letters left in the string
	int tmp = atoi(buffer);
	char buffer2[20];
	//there is no check if the number is bigger than 99999!!
	sprintf(buffer2, "%.5i", tmp+1);
	(*out) = buffer2;
}

void inputparameters::cond_assign(int * out, std::string line)
{
	if(std::strstr(line.c_str(),"cold")!=NULL) {
		(*out)=COLD_START;
		return;
	}
	if(std::strstr(line.c_str(),"hot")!=NULL) {
		(*out)=HOT_START;
		return;
	}
	if(std::strstr(line.c_str(),"source")!=NULL) {
		(*out)=START_FROM_SOURCE;
		return;
	}
	if(std::strstr(line.c_str(),"continue")!=NULL) {
		(*out)=START_FROM_SOURCE;
		return;
	}
	printf("invalid startcondition\n");
	exit(HMC_STDERR);
	return;
}

void inputparameters::fermact_assign(int * out, std::string line)
{
	if(std::strstr(line.c_str(),"TWISTEDMASS")!=NULL) {
		(*out)=TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(),"twistedmass")!=NULL) {
		(*out)=TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(),"Twistedmass")!=NULL) {
		(*out)=TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(),"TwistedMass")!=NULL) {
		(*out)=TWISTEDMASS;
		return;
	}
	if(std::strstr(line.c_str(),"clover")!=NULL) {
		(*out)=CLOVER;
		return;
	}
	if(std::strstr(line.c_str(),"CLOVER")!=NULL) {
		(*out)=CLOVER;
		return;
	}
	if(std::strstr(line.c_str(),"Clover")!=NULL) {
		(*out)=CLOVER;
		return;
	}
	if(std::strstr(line.c_str(),"WILSON")!=NULL) {
		(*out)=WILSON;
		return;
	}
	if(std::strstr(line.c_str(),"Wilson")!=NULL) {
		(*out)=WILSON;
		return;
	}
	if(std::strstr(line.c_str(),"wilson")!=NULL) {
		(*out)=WILSON;
		return;
	}
	if(std::strstr(line.c_str(),"unimproved")!=NULL) {
		(*out)=WILSON;
		return;
	}
	printf("invalid fermion action\n");
	exit(HMC_STDERR);
	return;
}


void inputparameters::savecond_assign(int * out, std::string line)
{
	if(std::strstr(line.c_str(),"yes")!=NULL) {
		(*out)=TRUE;
		return;
	}
	if(std::strstr(line.c_str(),"true")!=NULL) {
		(*out)=TRUE;
		return;
	}
	if(std::strstr(line.c_str(),"TRUE")!=NULL) {
		(*out)=TRUE;
		return;
	}
	if(std::strstr(line.c_str(),"True")!=NULL) {
		(*out)=TRUE;
		return;
	}
	if(std::strstr(line.c_str(),"no")!=NULL) {
		(*out)=FALSE;
		return;
	}
	if(std::strstr(line.c_str(),"false")!=NULL) {
		(*out)=FALSE;
		return;
	}
	if(std::strstr(line.c_str(),"FALSE")!=NULL) {
		(*out)=FALSE;
		return;
	}
	if(std::strstr(line.c_str(),"False")!=NULL) {
		(*out)=FALSE;
		return;
	}
	printf("invalid save condition\n");
	exit(HMC_STDERR);
	return;
}

void inputparameters::val_assign(std::string * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value=line.substr(pos+1);
	(*out) = value.c_str();
	return;
}

hmc_float inputparameters::get_kappa()
{
	return kappa;
}

hmc_float inputparameters::get_theta_fermion()
{
	return theta_fermion;
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
