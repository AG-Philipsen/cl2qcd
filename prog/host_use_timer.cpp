#include "host_use_timer.h"

#include "logger.hpp"

void usetimer::reset()
{
	timer.reset();
	return;
}

void usetimer::getTimeAndReset()
{
	time_measurement = timer.getTimeAndReset();
	//!!
	num_meas = 0;
	return;
}

void usetimer::add()
{
	time_measurement +=  timer.getTimeAndReset();
	num_meas ++;
	return;
}

void usetimer::zero()
{
	time_measurement = 0;
	//!!
	num_meas = 0;
	return;
}

uint64_t usetimer::getTime()
{
	return time_measurement;
}

int usetimer::getNumMeas()
{
	return num_meas;
}

uint64_t divide(uint64_t a, int b)
{
	return uint64_t ( ( (float) a ) / ((float) b) );
}

float percent(uint64_t a, int b)
{
	return (( (float) a) / ( (float )b )) * 100;
}


/*
void time_output(
  usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy
#ifdef _FERMIONS_
  , usetimer * inittimer, usetimer* singletimer, usetimer *Mtimer, usetimer *copytimer, usetimer *scalarprodtimer, usetimer *latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer
#endif
#ifdef _PERFORM_BENCHMARKS_
  , int steps
#endif
#ifdef _USEHMC_
  , usetimer * hmctimer, usetimer * leapfrogtimer, usetimer * hmcinittimer, usetimer * metropolistimer
#endif
)
{

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
	uint64_t dslashtime = (*dslashtimer).getTime();
	uint64_t Mdiagtime = (*Mdiagtimer).getTime();
#endif
#ifdef _USEHMC_
	uint64_t hmctime = (*hmctimer).getTime();
	uint64_t leapfrogtime = (*leapfrogtimer).getTime();
	uint64_t hmcinittime = (*hmcinittimer).getTime();
	uint64_t metropolistime = (*metropolistimer).getTime();
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
	int dslash_steps;
	int Mdiag_steps;
#endif
#ifdef _USEHMC_
	int hmc_steps;
	int hmcinit_steps;
	int leapfrog_steps;
	int metropolis_steps;
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
	uint64_t dslash_avgtime_site;
	uint64_t Mdiag_avgtime_site;
	uint64_t M_avgtime;
	uint64_t scalprod_avgtime;
	uint64_t la_avgtime;
	uint64_t copy_ferm_avgtime;
	uint64_t single_ferm_avgtime;
	uint64_t solver_avgtime;
	uint64_t dslash_avgtime;
	uint64_t Mdiag_avgtime;
#endif
#ifdef _USEHMC_
	uint64_t hmc_avgtime;
	uint64_t hmcinit_avgtime;
	uint64_t leapfrog_avgtime;
	uint64_t metropolis_avgtime;
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
	Mdiag_steps = (*Mdiagtimer).getNumMeas();
	dslash_steps = (*dslashtimer).getNumMeas();
#endif
#ifdef _USEHMC_
	hmc_steps = (*hmctimer).getNumMeas();
	hmcinit_steps = (*hmcinittimer).getNumMeas();
	leapfrog_steps = (*leapfrogtimer).getNumMeas();
	metropolis_steps = (*metropolistimer).getNumMeas();
#endif

	if(polysteps != 0) {
		poly_avgtime_site = divide(polytime, VOL4D * polysteps);
		poly_avgtime = divide(polytime, polysteps);
	} else {
		poly_avgtime_site = 0;
		poly_avgtime = 0;
	}
	if(plaqsteps != 0) {
		plaq_avgtime_site = divide(plaqtime, VOL4D * plaqsteps);
		plaq_avgtime = divide(plaqtime, plaqsteps);
	} else {
		plaq_avgtime_site = 0;
		plaq_avgtime = 0;
	}
	if(updatesteps != 0) {
		update_avgtime_site = divide(updatetime, VOL4D * updatesteps);
		update_avgtime = divide(updatetime, updatesteps);
	} else {
		update_avgtime_site = 0;
		update_avgtime = 0;
	}
	if(overrelaxsteps != 0) {
		overrelax_avgtime_site = divide(overrelaxtime, VOL4D * overrelaxsteps);
		overrelax_avgtime = divide(overrelaxtime, overrelaxsteps);
	} else {
		overrelax_avgtime_site = 0;
		overrelax_avgtime = 0;
	}
	if(copysteps != 0) copy_avgtime = divide(copytime, copysteps);
	else copy_avgtime = 0;
#ifdef _FERMIONS_
	if(single_ferm_steps != 0) {
		single_ferm_avgtime = divide(single_ferm, single_ferm_steps);
	} else {
		single_ferm_avgtime = 0;
	}
	if(copy_ferm_steps != 0) {
		copy_ferm_avgtime = divide(copy_ferm, copy_ferm_steps);
	} else {
		copy_ferm_avgtime = 0;
	}
	if(M_steps != 0) {
		M_avgtime = divide(Mtime, M_steps);
		M_avgtime_site = divide(Mtime, M_steps * VOL4D);
	} else {
		M_avgtime = 0;
		M_avgtime_site = 0;
	}
	if(Mdiag_steps != 0) {
		Mdiag_avgtime = divide(Mdiagtime, M_steps);
		Mdiag_avgtime_site = divide(Mdiagtime, M_steps * VOL4D);
	} else {
		Mdiag_avgtime = 0;
		Mdiag_avgtime_site = 0;
	}
	if(dslash_steps != 0) {
		dslash_avgtime = divide(dslashtime, dslash_steps);
		dslash_avgtime_site = divide(dslashtime, dslash_steps * VOL4D);
	} else {
		dslash_avgtime = 0;
		dslash_avgtime_site = 0;
	}
	if(scalprod_steps != 0) {
		scalprod_avgtime = divide(scalprod, scalprod_steps);
		scalprod_avgtime_site = divide(scalprod, scalprod_steps * VOL4D);
	} else {
		scalprod_avgtime = 0;
		scalprod_avgtime_site = 0;
	}
	if(la_steps != 0) {
		la_avgtime = divide(latime, la_steps);
		la_avgtime_site = divide(latime, la_steps * VOL4D);
	} else {
		la_avgtime = 0;
		la_avgtime_site = 0;
	}
	if(solver_steps != 0) {
		solver_avgtime = divide(solvertime, solver_steps);
	} else {
		solver_avgtime = 0;
	}
#endif
#ifdef _USEHMC_
	if(hmc_steps != 0) {
		hmc_avgtime = divide(hmctime, hmc_steps);
	} else {
		hmc_avgtime = 0;
	}
	if(hmcinit_steps != 0) {
		hmcinit_avgtime = divide(hmcinittime, hmcinit_steps);
	} else {
		hmcinit_avgtime = 0;
	}
	if(metropolis_steps != 0) {
		metropolis_avgtime = divide(metropolistime, metropolis_steps);
	} else {
		metropolis_avgtime = 0;
	}
	if(leapfrog_steps != 0) {
		leapfrog_avgtime = divide(leapfrogtime, leapfrog_steps);
	} else {
		leapfrog_avgtime = 0;
	}
#endif

	printf("\n");
	printf("**************************************************************\n");
	cout<<"total runtime:\t\t"<<totaltime<<"\n"<<endl;
	cout<<"Times:\t\t tot\t\t avg\t\tsite\tperc\n"<<endl;
	cout_times("Init.:",inittime,inittime,divide(inittime,VOL4D),percent(inittime,totaltime));
	cout_times("Copy.:",copytime,copy_avgtime,copytime,percent(copytime,totaltime));
	cout_times("Plaq.:", plaqtime, plaq_avgtime, plaq_avgtime_site, percent (plaqtime, totaltime));
	cout_times("Poly.:", polytime, poly_avgtime, poly_avgtime_site, percent (polytime, totaltime));
	cout_times("Updt.:", updatetime, update_avgtime, update_avgtime_site, percent (updatetime, totaltime));
	cout_times("Over.:", overrelaxtime, overrelax_avgtime, overrelax_avgtime_site, percent (overrelaxtime, totaltime));

#ifdef _FERMIONS_
	cout<<"Fermion Times:\t tot\t\t avg\t\tsite\tperc"<<endl;;
	cout_times("Init.:", init_ferm, init_ferm, divide(init_ferm, VOL4D), percent (init_ferm, totaltime) );
	cout_times("Solve:", solvertime, solver_avgtime, solver_avgtime, percent (solvertime, totaltime) );
	cout_times("Copy.:", copy_ferm, copy_ferm_avgtime, copy_ferm, percent (copy_ferm, totaltime) );
	cout_times("Sngle:", single_ferm, single_ferm_avgtime, single_ferm_avgtime, percent (single_ferm, totaltime));
	cout_times("ScPr.:", scalprod, scalprod_avgtime, scalprod_avgtime_site, percent (scalprod, totaltime));
	cout_times("BLAS.:", latime, la_avgtime, la_avgtime_site, percent (latime, totaltime));
	cout_times("Mferm:", Mtime, M_avgtime, M_avgtime_site, percent (Mtime, totaltime));
	cout_times("Mdiag:", Mdiagtime, Mdiag_avgtime, Mdiag_avgtime_site, percent (Mdiagtime, totaltime));
	cout_times("Dslas:", dslashtime, dslash_avgtime, dslash_avgtime_site, percent (dslashtime, totaltime));
#endif
#ifdef _USEHMC_
	cout<<"HybridMC Times:\t tot\t\t avg\t\tsite\tperc"<<endl;
	cout_times("HMC:", hmctime, hmc_avgtime, hmc_avgtime, percent (hmctime, totaltime));
	cout_times("Init:", hmcinittime, hmcinit_avgtime, hmcinit_avgtime, percent (hmcinittime, totaltime));
	cout_times("Leap:", leapfrogtime, leapfrog_avgtime, leapfrog_avgtime, percent (leapfrogtime, totaltime));
	cout_times("Metr:", metropolistime, metropolis_avgtime, metropolis_avgtime, percent (metropolistime, totaltime));
#endif
	printf("**************************************************************\n");


	//save some data to file
	ofstream out;
	stringstream str_filename;
	//CP: this H is for heatbath benchmarking, it has to be replaced meaningfully for other occasions
	str_filename << "time_";
#ifdef _PERFORM_BENCHMARKS_
	str_filename << "B_";
#endif
#ifdef _USEGPU_
	str_filename << "G_";
#else
	str_filename << "C_";
#endif
#ifdef _USEDOUBLEPREC_
	str_filename << "D_";
#else
	str_filename << "S_";
#endif
#ifdef _RECONSTRUCT_TWELVE_
	str_filename << "R_";
#else
	str_filename << "N_";
#endif
#ifdef _PERFORM_BENCHMARKS_
	str_filename << benchmark_id;
	out.open(str_filename.str().c_str(), fstream::trunc);
#else
	out.open(str_filename.str().c_str(), fstream::app);
#endif

	//gauge output
	if (out.is_open()) {
		//output:
		//(benchmark_id) NTIME   NSPACE   VOL4D   totaltime   inittimer  copytime  polytime   plaqtime   updatetime   overrelaxtime  (all times average per time-measurement)
		out <<
#ifdef _PERFORM_BENCHMARKS_
		    benchmark_id << "\t" <<  steps << "\t" <<
#endif
		    NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t"
		    << inittime << "\t" << copy_avgtime << "\t" << poly_avgtime << "\t" << plaq_avgtime << "\t" << update_avgtime << "\t" << overrelax_avgtime << endl;
		out.close();
	} else logger.error() << "Unable to open file for output";

	//fermion output
#ifdef _FERMIONS_
	str_filename << "ferm";
	out.open(str_filename.str().c_str(), fstream::app);
	if (out.is_open()) {
		//output:
		//(benchmark_id) NTIME   NSPACE   VOL4D   totaltime   inittime  solvertime copytime singletime   scalarproducttime   latime  Mtime Mdiagtime dslashtime(all times average time-measurement)
		out <<
#ifdef _PERFORM_BENCHMARKS_
		    benchmark_id << "\t"  <<  steps << "\t" <<
#endif
		    NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t"
		    << init_ferm << "\t" << solvertime << "\t" << copy_ferm_avgtime << "\t" << single_ferm_avgtime << "\t" << scalprod_avgtime << "\t" << la_avgtime << "\t" << M_avgtime << "\t" << Mdiag_avgtime << "\t" << dslash_avgtime << endl;
		out.close();

	} else logger.error() << "Unable to open file for output";
#endif

	//hybrid monte carlo output
#ifdef _USEHMC_
	str_filename << "ferm";
	out.open(str_filename.str().c_str(), fstream::app);
	if (out.is_open()) {
		//output:
		//(benchmark_id) NTIME   NSPACE   VOL4D   totaltime  hmc_totaltime  hmcinittime  leapfrogtime   metropolistime (all times average time-measurement)
		out <<
#ifdef _PERFORM_BENCHMARKS_
		    benchmark_id << "\t"  <<  steps << "\t" <<
#endif
		    NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t"
		    << hmc_avgtime << "\t" << hmcinit_avgtime << "\t" << leapfrog_avgtime << "\t" << metropolis_avgtime << "\t" << endl;
		out.close();
	} else logger.error() << "Unable to open file for output";
#endif

	return;
}
void cout_times(const char* name, uint64_t timeA, uint64_t timeB, uint64_t timeC, float perc){
  size_t coutprec = cout.precision();
  cout<<name<<"\t"<<setw(12)<<timeA<<"\t"<<setw(12)<<timeB<<"\t"<<setw(12)<<timeC<<"\t";
  cout.precision(3);
  cout<<perc<<endl;
  cout.precision(coutprec);
  return;
}
*/

void time_output_heatbath(
  usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy
)
{

	uint64_t totaltime = (*total).getTime();
	uint64_t inittime = (*init).getTime();
	uint64_t polytime = (*poly).getTime();
	uint64_t plaqtime = (*plaq).getTime();
	uint64_t updatetime = (*update).getTime();
	uint64_t overrelaxtime = (*overrelax).getTime();
	uint64_t copytime = (*copy).getTime();

	int polysteps;
	int plaqsteps;
	int updatesteps;
	int overrelaxsteps;
	int copysteps;

	uint64_t poly_avgtime_site;
	uint64_t plaq_avgtime_site;
	uint64_t update_avgtime_site;
	uint64_t overrelax_avgtime_site;
	uint64_t poly_avgtime;
	uint64_t plaq_avgtime;
	uint64_t update_avgtime;
	uint64_t overrelax_avgtime;
	uint64_t copy_avgtime;

	polysteps = (*poly).getNumMeas();
	plaqsteps = (*plaq).getNumMeas();
	updatesteps = (*update).getNumMeas();
	overrelaxsteps = (*overrelax).getNumMeas();
	copysteps = (*copy).getNumMeas();

	if(polysteps != 0) {
		poly_avgtime_site = divide(polytime, VOL4D * polysteps);
		poly_avgtime = divide(polytime, polysteps);
	} else {
		poly_avgtime_site = 0;
		poly_avgtime = 0;
	}
	if(plaqsteps != 0) {
		plaq_avgtime_site = divide(plaqtime, VOL4D * plaqsteps);
		plaq_avgtime = divide(plaqtime, plaqsteps);
	} else {
		plaq_avgtime_site = 0;
		plaq_avgtime = 0;
	}
	if(updatesteps != 0) {
		update_avgtime_site = divide(updatetime, VOL4D * updatesteps);
		update_avgtime = divide(updatetime, updatesteps);
	} else {
		update_avgtime_site = 0;
		update_avgtime = 0;
	}
	if(overrelaxsteps != 0) {
		overrelax_avgtime_site = divide(overrelaxtime, VOL4D * overrelaxsteps);
		overrelax_avgtime = divide(overrelaxtime, overrelaxsteps);
	} else {
		overrelax_avgtime_site = 0;
		overrelax_avgtime = 0;
	}
	if(copysteps != 0) copy_avgtime = divide(copytime, copysteps);
	else copy_avgtime = 0;

	logger.trace() << "*******************************************************************";
	logger.trace() << "total runtime:" << setfill(' ') << setw(12) << totaltime;
	logger.trace() << "Times:\t" << setfill(' ') << setw(12) << "tot" << '\t' << setw(12) << "avg" << '\t' << setw(12) << "site" << '\t' << setw(5) << "perc";

	logger.trace() << "Init.:\t" << setfill(' ') << setw(12) << inittime << '\t' << setw(12) << inittime << '\t' << setw(12) << divide(inittime, VOL4D) << '\t' << setw(5) << fixed << setw(5) << setprecision(1) << percent(inittime, totaltime) ;
	logger.trace() << "Copy.:\t" << setfill(' ') << setw(12) << copytime << '\t' << setw(12) << copy_avgtime << '\t' << setw(12) << copytime << '\t' << fixed << setw(5) << setprecision(1) << percent(copytime, totaltime);
	logger.trace() << "Plaq.:\t" << setfill(' ') << setw(12) << plaqtime << '\t' << setw(12) << plaq_avgtime << '\t' << setw(12) << plaq_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(plaqtime, totaltime);
	logger.trace() << "Poly.:\t" << setfill(' ') << setw(12) << polytime << '\t' << setw(12) << poly_avgtime << '\t' << setw(12) << poly_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(polytime, totaltime);
	logger.trace() << "Updt.:\t" << setfill(' ') << setw(12) << updatetime << '\t' << setw(12) << update_avgtime << '\t' << setw(12) << update_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(updatetime, totaltime);
	logger.trace() << "Over.:\t" << setfill(' ') << setw(12) << overrelaxtime << '\t' << setw(12) << overrelax_avgtime << '\t' << setw(12) << overrelax_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(overrelaxtime, totaltime);
	logger.trace() << "*******************************************************************";


	//save some data to file
	ofstream out;
	stringstream str_filename;
	//CP: this H is for heatbath benchmarking, it has to be replaced meaningfully for other occasions
	str_filename << "time_";
#ifdef _USEGPU_
	str_filename << "G_";
#else
	str_filename << "C_";
#endif
#ifdef _USEDOUBLEPREC_
	str_filename << "D_";
#else
	str_filename << "S_";
#endif
#ifdef _RECONSTRUCT_TWELVE_
	str_filename << "R_";
#else
	str_filename << "N_";
#endif
	out.open(str_filename.str().c_str(), fstream::app);

	//gauge output
	if (out.is_open()) {
		//output:
		//(benchmark_id) NTIME   NSPACE   VOL4D   totaltime   inittimer  copytime  polytime   plaqtime   updatetime   overrelaxtime  (all times average per time-measurement)
		out <<
		    NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t"
		    << inittime << "\t" << copy_avgtime << "\t" << poly_avgtime << "\t" << plaq_avgtime << "\t" << update_avgtime << "\t" << overrelax_avgtime << endl;
		out.close();
	} else logger.error() << "Unable to open file for output";

	return;
}



void time_output_inverter(
  usetimer * total, usetimer * init, usetimer * poly, usetimer * plaq, usetimer * update, usetimer * overrelax, usetimer * copy
  , usetimer * inittimer, usetimer* singletimer, usetimer *Mtimer, usetimer *copytimer, usetimer *scalarprodtimer, usetimer *latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer)
{

	uint64_t totaltime = (*total).getTime();

	uint64_t init_ferm = (*inittimer).getTime();
	uint64_t single_ferm = (*singletimer).getTime();
	uint64_t Mtime = (*Mtimer).getTime();
	uint64_t copy_ferm = (*copytimer).getTime();
	uint64_t scalprod = (*scalarprodtimer).getTime();
	uint64_t latime = (*latimer).getTime();
	uint64_t solvertime = (*solvertimer).getTime();
	uint64_t dslashtime = (*dslashtimer).getTime();
	uint64_t Mdiagtime = (*Mdiagtimer).getTime();

	int single_ferm_steps;
	int M_steps;
	int copy_ferm_steps;
	int scalprod_steps;
	int la_steps;
	int solver_steps;
	int dslash_steps;
	int Mdiag_steps;

	uint64_t M_avgtime_site;
	uint64_t scalprod_avgtime_site;
	uint64_t la_avgtime_site;
	uint64_t dslash_avgtime_site;
	uint64_t Mdiag_avgtime_site;
	uint64_t M_avgtime;
	uint64_t scalprod_avgtime;
	uint64_t la_avgtime;
	uint64_t copy_ferm_avgtime;
	uint64_t single_ferm_avgtime;
	uint64_t solver_avgtime;
	uint64_t dslash_avgtime;
	uint64_t Mdiag_avgtime;

	single_ferm_steps = (*singletimer).getNumMeas();
	M_steps = (*Mtimer).getNumMeas();
	copy_ferm_steps = (*copytimer).getNumMeas();
	scalprod_steps = (*scalarprodtimer).getNumMeas();
	la_steps = (*latimer).getNumMeas();
	solver_steps = (*solvertimer).getNumMeas();
	Mdiag_steps = (*Mdiagtimer).getNumMeas();
	dslash_steps = (*dslashtimer).getNumMeas();

	if(single_ferm_steps != 0) {
		single_ferm_avgtime = divide(single_ferm, single_ferm_steps);
	} else {
		single_ferm_avgtime = 0;
	}
	if(copy_ferm_steps != 0) {
		copy_ferm_avgtime = divide(copy_ferm, copy_ferm_steps);
	} else {
		copy_ferm_avgtime = 0;
	}
	if(M_steps != 0) {
		M_avgtime = divide(Mtime, M_steps);
		M_avgtime_site = divide(Mtime, M_steps * VOL4D);
	} else {
		M_avgtime = 0;
		M_avgtime_site = 0;
	}
	if(Mdiag_steps != 0) {
		Mdiag_avgtime = divide(Mdiagtime, M_steps);
		Mdiag_avgtime_site = divide(Mdiagtime, M_steps * VOL4D);
	} else {
		Mdiag_avgtime = 0;
		Mdiag_avgtime_site = 0;
	}
	if(dslash_steps != 0) {
		dslash_avgtime = divide(dslashtime, dslash_steps);
		dslash_avgtime_site = divide(dslashtime, dslash_steps * VOL4D);
	} else {
		dslash_avgtime = 0;
		dslash_avgtime_site = 0;
	}
	if(scalprod_steps != 0) {
		scalprod_avgtime = divide(scalprod, scalprod_steps);
		scalprod_avgtime_site = divide(scalprod, scalprod_steps * VOL4D);
	} else {
		scalprod_avgtime = 0;
		scalprod_avgtime_site = 0;
	}
	if(la_steps != 0) {
		la_avgtime = divide(latime, la_steps);
		la_avgtime_site = divide(latime, la_steps * VOL4D);
	} else {
		la_avgtime = 0;
		la_avgtime_site = 0;
	}
	if(solver_steps != 0) {
		solver_avgtime = divide(solvertime, solver_steps);
	} else {
		solver_avgtime = 0;
	}

	time_output_heatbath(total,init,poly,plaq,update,overrelax,copy);

	logger.trace() << "*******************************************************************";
	logger.trace() << "Fermion times:\t" << setfill(' ') << setw(12) << "tot" << '\t' << setw(12) << "avg" << '\t' << setw(12) << "site" << '\t' << setw(5) << "perc";

	logger.trace() << "Init.:\t" << setfill(' ') << setw(12) << init_ferm << '\t' << setw(12) << init_ferm << '\t' << setw(12) << divide(init_ferm, VOL4D) << '\t' << setw(5) << fixed << setw(5) << setprecision(1) << percent(init_ferm, totaltime) ;
	logger.trace() << "Solve.:\t" << setfill(' ') << setw(12) << solvertime << '\t' << setw(12) << solver_avgtime << '\t' << setw(12) << solver_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(solvertime, totaltime);
	logger.trace() << "Copy.:\t" << setfill(' ') << setw(12) << copy_ferm << '\t' << setw(12) << copy_ferm_avgtime << '\t' << setw(12) << copy_ferm << '\t' << fixed << setw(5) << setprecision(1) << percent(copy_ferm, totaltime);
	logger.trace() << "Sngle.:\t" << setfill(' ') << setw(12) << single_ferm << '\t' << setw(12) << single_ferm_avgtime << '\t' << setw(12) << single_ferm_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(single_ferm, totaltime);
	logger.trace() << "ScPr.:\t" << setfill(' ') << setw(12) << scalprod << '\t' << setw(12) << scalprod_avgtime << '\t' << setw(12) << scalprod_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(scalprod, totaltime);
	logger.trace() << "BLAS:\t" << setfill(' ') << setw(12) << latime << '\t' << setw(12) << la_avgtime << '\t' << setw(12) << la_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(latime, totaltime);
	logger.trace() << "Mferm:\t" << setfill(' ') << setw(12) << Mtime << '\t' << setw(12) << M_avgtime << '\t' << setw(12) << M_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(Mtime, totaltime);
	logger.trace() << "Mdiag:\t" << setfill(' ') << setw(12) << Mdiagtime << '\t' << setw(12) << Mdiag_avgtime << '\t' << setw(12) << Mdiag_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(Mdiagtime, totaltime);
	logger.trace() << "Dslas:\t" << setfill(' ') << setw(12) << dslashtime << '\t' << setw(12) << dslash_avgtime << '\t' << setw(12) << dslash_avgtime_site << '\t' << fixed << setw(5) << setprecision(1) << percent(dslashtime, totaltime);
	logger.trace() << "*******************************************************************";

	//save some data to file
	ofstream out;
	stringstream str_filename;
	//CP: this H is for heatbath benchmarking, it has to be replaced meaningfully for other occasions
	str_filename << "time_";
#ifdef _USEGPU_
	str_filename << "G_";
#else
	str_filename << "C_";
#endif
#ifdef _USEDOUBLEPREC_
	str_filename << "D_";
#else
	str_filename << "S_";
#endif
#ifdef _RECONSTRUCT_TWELVE_
	str_filename << "R_";
#else
	str_filename << "N_";
#endif

	str_filename << "inverter";
	out.open(str_filename.str().c_str(), fstream::app);
	if (out.is_open()) {
		//output:
		//NTIME   NSPACE   VOL4D   totaltime   inittime  solvertime copytime singletime   scalarproducttime   latime  Mtime Mdiagtime dslashtime(all times average time-measurement)
		out <<
		    NTIME << "\t" << NSPACE << "\t" << VOL4D << "\t" << totaltime << "\t"
		    << init_ferm << "\t" << solvertime << "\t" << copy_ferm_avgtime << "\t" << single_ferm_avgtime << "\t" << scalprod_avgtime << "\t" << la_avgtime << "\t" << M_avgtime << "\t" << Mdiag_avgtime << "\t" << dslash_avgtime << endl;
		out.close();

	} else logger.error() << "Unable to open file for output";
	return;
}
