#include "host_input.h"

hmc_error inputparameters::set_defaults(){
  kappa = 0.125;
  beta = 4.0;
  mu = 0.006;
  cgmax = 1000;
  prec = 64;
  startcondition = COLD_START;
  thermalizationsteps = 0;
  heatbathsteps = 1000;
  writefrequency = 1;
  savefrequency = 100;
  //sourcefile = "\0";
  sourcefilenumber = "00000";
  return HMC_SUCCESS;
}

hmc_error inputparameters::readfile(char* ifn){
  std::ifstream infile;
  infile.open(ifn);
  if(!infile.is_open()) {
    printf("Could not open input file: %s\n",ifn);
    exit(HMC_FILEERROR);
  }
  while (infile.good()) {
    std::string line;
    infile>>line;
    if(line.find("#")!=std::string::npos) continue; //allow comments
    if(line.find("kappa")!=std::string::npos) val_assign(&kappa,line);
    if(line.find("Kappa")!=std::string::npos) val_assign(&kappa,line);
    if(line.find("mu")!=std::string::npos) val_assign(&mu,line);
    if(line.find("Mu")!=std::string::npos) val_assign(&mu,line);
    if(line.find("beta")!=std::string::npos) val_assign(&beta,line);
    if(line.find("Beta")!=std::string::npos) val_assign(&beta,line);
    if(line.find("cgmax")!=std::string::npos) val_assign(&cgmax,line);
    if(line.find("CGmax")!=std::string::npos) val_assign(&cgmax,line);
    if(line.find("Cgmax")!=std::string::npos) val_assign(&cgmax,line);
    if(line.find("writefrequency")!=std::string::npos) val_assign(&writefrequency,line);
    if(line.find("savefrequency")!=std::string::npos) val_assign(&savefrequency,line);
    if(line.find("prec")!=std::string::npos) val_assign(&prec,line);
    if(line.find("Prec")!=std::string::npos) val_assign(&prec,line);
    if(line.find("readsource")!=std::string::npos) cond_assign(&startcondition,line);
    if(line.find("startcondition")!=std::string::npos) cond_assign(&startcondition,line);
    if(line.find("sourcefile")!=std::string::npos){
      val_assign(&sourcefile,line);
      sourcefilenumber_assign(&sourcefilenumber);
    }
    if(line.find("thermalizationsteps")!=std::string::npos) val_assign(&thermalizationsteps,line);
    if(line.find("heatbathsteps")!=std::string::npos) val_assign(&heatbathsteps,line);
  }
  return HMC_SUCCESS;
}

void inputparameters::val_assign(hmc_float * out, std::string line) {
  size_t pos = line.find("=");
  std::string value=line.substr(pos+1);
  (*out) = atof(value.c_str());
  return;
}

void inputparameters::val_assign(int * out, std::string line) {
  size_t pos = line.find("=");
  std::string value=line.substr(pos+1);
  (*out) = atoi(value.c_str());
  return;
}

void inputparameters::sourcefilenumber_assign(std::string * out){
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

void inputparameters::cond_assign(int * out, std::string line) {
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

void inputparameters::val_assign(std::string * out, std::string line) {
  size_t pos = line.find("=");
  std::string value=line.substr(pos+1);
  (*out) = value.c_str();
  return;
}

hmc_float inputparameters::get_kappa(){
  return kappa;
}

hmc_float inputparameters::get_beta(){
  return beta;
}

hmc_float inputparameters::get_mu(){
  return mu;
}

int inputparameters::get_cgmax(){
  return cgmax;
}

int inputparameters::get_prec(){
  return prec;
}

int inputparameters::get_thermalizationsteps(){
  return thermalizationsteps;
}

int inputparameters::get_heatbathsteps(){
  return heatbathsteps;
}

int inputparameters::get_writefrequency(){
  return writefrequency;
}

int inputparameters::get_savefrequency(){
  return savefrequency;
}

int inputparameters::get_startcondition(){
  return startcondition;
}

void inputparameters::display_sourcefile(){
  cout << sourcefile;
}

void inputparameters::display_sourcefilenumber(){
  cout << sourcefilenumber;
}
