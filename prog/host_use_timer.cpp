#include "host_use_timer.h"

void usetimer::reset() {
  timer.reset();
  return;
}

void usetimer::getTimeAndReset(){
  time_measurement = timer.getTimeAndReset();
  //!!
  num_meas = 0;
  return;
}

void usetimer::add(){
  time_measurement +=  timer.getTimeAndReset();
  num_meas ++;
  return;
}

void usetimer::zero(){
  time_measurement = 0;
  //!!
  num_meas = 0;
  return;
}

uint64_t usetimer::getTime(){
  return time_measurement;
}

int usetimer::getNumMeas(){
  return num_meas;
}

uint64_t divide(uint64_t a, int b){
  return uint64_t ( ( (float) a )/ ((float) b) );
}

float percent(uint64_t a, int b){
  return (( (float) a) / ( (float )b )) *100;
}

void time_output(
	usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy
#ifdef _FERMIONS_
, usetimer * inittimer, usetimer* singletimer, usetimer *Mtimer, usetimer *copytimer, usetimer *scalarprodtimer, usetimer *latimer, usetimer * solvertimer
#endif	
	) {

  uint64_t totaltime = (*total).getTime();
  uint64_t inittime = (*init).getTime();
  uint64_t polytime = (*poly).getTime();
  uint64_t plaqtime = (*plaq).getTime();
  uint64_t updatetime = (*update).getTime();
  uint64_t overrelaxtime = (*overrelax).getTime();
  uint64_t copytime = (*copy).getTime();

#ifdef _FERMIONS_
  uint64_t init_ferm = (*inittimer).getTime();
  uint64_t single_ferm = (*singletimer).getTime();
  uint64_t Mtime = (*Mtimer).getTime();
  uint64_t copy_ferm = (*copytimer).getTime();
  uint64_t scalprod = (*scalarprodtimer).getTime();
  uint64_t latime = (*latimer).getTime();
	uint64_t solvertime = (*solvertimer).getTime();
#endif
	
  int polysteps;
  int plaqsteps;
  int updatesteps;
  int overrelaxsteps;
  int copysteps;
#ifdef _FERMIONS_
	int single_ferm_steps;
  int M_steps;
  int copy_ferm_steps;
  int scalprod_steps;
  int la_steps;
	int solver_steps;
#endif

	uint64_t poly_avgtime_site;
  uint64_t plaq_avgtime_site;
  uint64_t update_avgtime_site;
  uint64_t overrelax_avgtime_site;
  uint64_t poly_avgtime;
  uint64_t plaq_avgtime;
  uint64_t update_avgtime;
  uint64_t overrelax_avgtime;
  uint64_t copy_avgtime;
#ifdef _FERMIONS_ 
  uint64_t M_avgtime_site;
  uint64_t scalprod_avgtime_site;
  uint64_t la_avgtime_site;
  uint64_t M_avgtime;
  uint64_t scalprod_avgtime;
  uint64_t la_avgtime;
  uint64_t copy_ferm_avgtime;	
	uint64_t single_ferm_avgtime;
	uint64_t solver_avgtime;
#endif	
	
  polysteps = (*poly).getNumMeas();
  plaqsteps = (*plaq).getNumMeas();
  updatesteps = (*update).getNumMeas();
  overrelaxsteps = (*overrelax).getNumMeas();
  copysteps = (*copy).getNumMeas();
#ifdef _FERMIONS_ 
	single_ferm_steps = (*singletimer).getNumMeas();
  M_steps = (*Mtimer).getNumMeas();
  copy_ferm_steps = (*copytimer).getNumMeas();
  scalprod_steps = (*scalarprodtimer).getNumMeas();
  la_steps = (*latimer).getNumMeas();
	solver_steps = (*solvertimer).getNumMeas();
#endif		
	
  if(polysteps!=0){
    poly_avgtime_site = divide(polytime, VOL4D*polysteps);
    poly_avgtime = divide(polytime, polysteps);
  }
  else{
    poly_avgtime_site = 0;
    poly_avgtime = 0;
  }
  if(plaqsteps!=0){
    plaq_avgtime_site = divide(plaqtime, VOL4D*plaqsteps);
    plaq_avgtime = divide(plaqtime, plaqsteps);
  }
  else{
    plaq_avgtime_site = 0;
    plaq_avgtime = 0;
  }
  if(updatesteps!=0){
    update_avgtime_site = divide(updatetime, VOL4D*updatesteps);
    update_avgtime = divide(updatetime, updatesteps);
  }
  else{
    update_avgtime_site = 0;
    update_avgtime = 0;
  }
  if(overrelaxsteps!=0){
    overrelax_avgtime_site = divide(overrelaxtime, VOL4D*overrelaxsteps);
    overrelax_avgtime = divide(overrelaxtime, overrelaxsteps);
  }
  else{
    overrelax_avgtime_site = 0;
    overrelax_avgtime = 0;
  }
  if(copysteps != 0) copy_avgtime = divide(copytime, copysteps);
  else copy_avgtime = 0;
#ifdef _FERMIONS_
	if(single_ferm_steps!=0){
    single_ferm_avgtime = divide(single_ferm, single_ferm_steps);
  }
  else{
    single_ferm_avgtime = 0;
  }
  if(copy_ferm_steps!=0){
    copy_ferm_avgtime = divide(copy_ferm, copy_ferm_steps);
  }
  else{
    copy_ferm_avgtime = 0;
  }
  if(M_steps!=0){
    M_avgtime = divide(Mtime, M_steps);
		M_avgtime_site = divide(Mtime, M_steps*VOL4D);
  }
  else{
    M_avgtime = 0;
		M_avgtime_site = 0;
  }
  if(scalprod_steps!=0){
    scalprod_avgtime = divide(scalprod, scalprod_steps);
		scalprod_avgtime_site = divide(scalprod, scalprod_steps*VOL4D);
  }
  else{
    scalprod_avgtime = 0;
		scalprod_avgtime_site = 0;
  }
	if(la_steps!=0){
    la_avgtime = divide(latime, la_steps);
		la_avgtime_site = divide(latime, la_steps*VOL4D);
  }
  else{
    la_avgtime = 0;
		la_avgtime_site = 0;
  }
	if(solver_steps!=0){
    solver_avgtime = divide(solvertime, solver_steps);
  }
  else{
    solver_avgtime = 0;
  }
#endif
  
  printf("\n");
  printf("**************************************************************\n");
  printf("total runtime:\t\t%lu\n\n", totaltime );
  printf("Times:\t\t tot\t\t avg\t\tsite\tperc\n");
  
  printf("Init.:\t%12lu\t%12lu\t%12lu\t%.3f\n", inittime, inittime, divide(inittime, VOL4D),percent (inittime, totaltime) );
  printf("Copy.:\t%12lu\t%12lu\t%12lu\t%.3f\n", copytime, copy_avgtime, copytime,percent (copytime, totaltime) );
  printf("Plaq.:\t%12lu\t%12lu\t%12lu\t%.3f\n", plaqtime, plaq_avgtime, plaq_avgtime_site, percent (plaqtime, totaltime));
  printf("Poly.:\t%12lu\t%12lu\t%12lu\t%.3f\n", polytime, poly_avgtime, poly_avgtime_site, percent (polytime, totaltime));
  printf("Updt.:\t%12lu\t%12lu\t%12lu\t%.3f\n", updatetime, update_avgtime, update_avgtime_site, percent (updatetime, totaltime));
  printf("Over.:\t%12lu\t%12lu\t%12lu\t%.3f\n", overrelaxtime, overrelax_avgtime, overrelax_avgtime_site, percent (overrelaxtime, totaltime));
#ifdef _FERMIONS_
	printf("Fermion Times:\t tot\t\t avg\t\tsite\tperc\n");
	printf("Init.:\t%12lu\t%12lu\t%12lu\t%.3f\n", init_ferm, init_ferm, divide(init_ferm, VOL4D),percent (init_ferm, totaltime) );
	printf("Solve:\t%12lu\t%12lu\t%12lu\t%.3f\n", solvertime, solver_avgtime, solver_avgtime,percent (solvertime, totaltime) );
	printf("Copy.:\t%12lu\t%12lu\t%12lu\t%.3f\n", copy_ferm, copy_ferm_avgtime, copy_ferm,percent (copy_ferm, totaltime) );
  printf("Sngle:\t%12lu\t%12lu\t%12lu\t%.3f\n", single_ferm, single_ferm_avgtime, single_ferm_avgtime, percent (single_ferm, totaltime));
	printf("ScPr.:\t%12lu\t%12lu\t%12lu\t%.3f\n", scalprod, scalprod_avgtime, scalprod_avgtime_site, percent (scalprod, totaltime));
	printf("BLAS.:\t%12lu\t%12lu\t%12lu\t%.3f\n", latime, la_avgtime, la_avgtime_site, percent (latime, totaltime));
	printf("Mferm:\t%12lu\t%12lu\t%12lu\t%.3f\n", Mtime, M_avgtime, M_avgtime_site, percent (Mtime, totaltime));
#endif
  printf("**************************************************************\n");
  
  
  //save some data to file
  ofstream out;
  stringstream str_filename;
  str_filename<<"time_H_";
#ifdef _USEGPU_
  str_filename<<"G_";
#else
  str_filename<<"C_";
#endif
#ifdef _USEDOUBLEPREC_
  str_filename<<"D_";
#else
  str_filename<<"S_";
#endif
#ifdef _RECONSTRUCT_TWELVE_
  str_filename<<"R_";
#else
  str_filename<<"N_";
#endif
  out.open(str_filename.str().c_str(), fstream::app); 
  if (out.is_open())
  {
    //output:
    //NTIME   NSPACE   VOL4D   totaltime   inittimer   polytime   plaqtime   updatetime   overrelaxtime  (all times average per time-measurement)
    out << NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t" 
    << inittime << "\t" << copy_avgtime << "\t" << poly_avgtime << "\t" << plaq_avgtime << "\t" << update_avgtime << "\t" << overrelax_avgtime << endl;
    out.close();
  }
  else cout << "Unable to open file for output" << endl;

#ifdef _FERMIONS_
	str_filename<<"ferm";
	out.open(str_filename.str().c_str(), fstream::app); 
  if (out.is_open())
  {    
		//output:
    //NTIME   NSPACE   VOL4D   totaltime   inittime   solvertime copytime singletime   scalarproducttime   latime  Mtime (all times average time-measurement)
    out << NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t" 
    << init_ferm << "\t" << solvertime << "\t" << copy_ferm_avgtime << "\t" << single_ferm_avgtime << "\t" << scalprod_avgtime << "\t" << la_avgtime << "\t" << M_avgtime << endl;
    out.close();

  }
  else cout << "Unable to open file for output" << endl;
#endif
  return;
}
