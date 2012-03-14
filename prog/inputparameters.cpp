#include "inputparameters.h"

#include "host_geometry.h"

#include "gitcommitid.h"

void inputparameters::set_defaults()
{
	//global parameters
	// default precision is 32 for float, 64 for double
	// this is only to ease tests. real input files should
	// always explicitly specify the precision to use
	prec = sizeof(hmc_float) * 8;

	use_rec12 = false;
	use_gpu = false;
	use_aniso = false;
	use_chem_pot_re = false;
	use_chem_pot_im = false;
	use_smearing = false;
	use_mp = false;
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

	mat_size = 9;
	float_size = 8;
	//complex_add is 2 fl. operations. There is no variable for that!
	//( a_re * b_re - b_im * a_im , a_re * b_im + a_im * b_re )
	flop_complex_mult = 6;
	//  1 entry: NC * complex mults and NC-1 complex adds
	//  NC*NC entries total
	flop_su3_su3 = (flop_complex_mult * NC + (NC - 1) * 2) * NC * NC;
	//  1 entry: NC * complex mults and NC-1 complex adds
	//  NC entries total
	flop_su3_su3vec = (flop_complex_mult * NC + (NC - 1) * 2) * NC;
	//  NC * complex adds
	flop_su3trace = (NC - 1) * 2;
	//  in prinicple, a gamma_matrix is a complex NDIM * NDIM matrix acting on NDIM vectors,
	//    which then again have NC entries for each of the NDIM entries.
	//    the matrix-mult alone are (flop_complex_mult * NDIM + (NDIM-1) * 2) * NDIM
	//    but for each such operation one has NC complex multiplications in addition.
	//    @TODO: gamma real??
	//    @TODO: simplifications possible!!
	flop_gamma_spinor = (flop_complex_mult * NDIM + (NDIM - 1) * 2) * NDIM * NC * flop_complex_mult;
	//  NDIM * flop_su3_su3vec
	flop_su3_spinor = NDIM * flop_su3_su3vec;
	//  NDIM * NC * complex_mult + ( NDIM * NC -1 ) complex adds
	flop_spinor_spinor = NDIM * NC * flop_complex_mult + (NDIM * NC - 1) * 2;
	//  NDIM * NC * 0.5 complex_mult + ( NDIM * NC -1 ) real adds
	flop_spinor_sqnorm = NDIM * NC * flop_complex_mult * 0.5 + (NC * NDIM - 1);

	//gaugefield parameters
	beta = 4.0;
	theta_gaugefield = 0.;
	rho = 0.;
	rho_iter = 0;
	gaugeact = WILSON;
	use_rectangles = false;
	c0 = 1.;
	c1 = 0.;
	c0_default_wilson = 1.;
	c1_default_tlsym = -0.083333333;
	c1_default_iwasaki = -0.331;
	c1_default_dbw2 = -1.4069;

	//heatbath parameters
	thermalizationsteps = 0;
	heatbathsteps = 1000;
	overrelaxsteps = 1;
	xi = 1;

	//fermionic parameters
	fermact = WILSON;
	kappa = 0.125;
	mu = 0.006;
	csw = 0.;
	kappa_mp = 0.125;
	mu_mp = 0.006;
	csw_mp = 0.;
	cgmax = 1000;
	cgmax_mp = cgmax;
	theta_fermion_spatial = 0.;
	theta_fermion_temporal = 0.;
	chem_pot_re = 0.;
	chem_pot_im = 0.;
	use_eo = true;
	//at the moment, only 2 solvers are implemented..
	use_cg = false;
	use_cg_mp = use_cg;
	use_bicgstab_save = false;
	use_bicgstab_save_mp = use_bicgstab_save;
	use_pointsource = true;
	num_sources = 12;
	pointsource_x = 0;
	pointsource_y = 0;
	pointsource_z = 0;
	pointsource_t = 0;
#ifdef _USEDOUBLEPREC_
	solver_prec = 1e-23;
	force_prec = 1e-12;
#else
	solver_prec = 1e-16;
	force_prec = 1e-8;
#endif
	iter_refresh = 100;
	iter_refresh_mp = iter_refresh;

	//HMC specific parameters
	tau = 0.5;
	reversibility_check = false;
	integrationsteps0 = 10;
	integrationsteps1 = integrationsteps0;
	hmcsteps = 10;
	num_timescales = 1;
	integrator0 = LEAPFROG;
	integrator1 = LEAPFROG;
	integrator2 = LEAPFROG;
	//this is the optimal value...
	lambda0 = 0.1931833275037836;
	lambda1 = 0.1931833275037836;
	lambda2 = 0.1931833275037836;

	//direction for the correlator
	corr_dir = 3;

	use_same_rnd_numbers = false;
	profile_solver = false;
	
	return;
}

void inputparameters::set_gauge_norm_factors() {
	
	//set normalization factors for gauge observables
	plaq_norm = (this->get_vol4d() * NDIM * (NDIM - 1) * NC) / 2.;
	tplaq_norm = (this->get_vol4d() * NC * (NDIM - 1));
	splaq_norm = (this->get_vol4d() * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	rect_norm = NDIM * (NDIM-1) * NC *  this->get_vol4d();
	poly_norm = (NC * this->get_volspace());
	
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
		bool csw_mpset = false;
		bool mu_mpset = false;
		bool c1set = false;

		while (infile.good()) {
			char linebuf[256];
			infile.getline(linebuf, 256);
			std::string line(linebuf);
			if(line.length() == 255) {
				/// @todo Handle such lines properly
				logger.fatal() << "The file contains a line longer than 255 characters - bailing out";
			}
			if(line.find("#") != std::string::npos) continue; //allow comments

			if(	line.find("kappa") != std::string::npos ||
					line.find("Kappa") != std::string::npos ) {
					if(	line.find("mp") != std::string::npos || 
							line.find("MP") != std::string::npos || 
							line.find("Mp") != std::string::npos )
						val_assign(&kappa_mp, line);
					else
						val_assign(&kappa, line);
			}
			if(	line.find("mu") != std::string::npos || 
					line.find("Mu") != std::string::npos) {
					if(	line.find("mp") != std::string::npos || 
							line.find("MP") != std::string::npos || 
							line.find("Mp") != std::string::npos ) {
						val_assign(&mu_mp, line);
						//calculate mubar
						this->calc_mubar_mp();
						mu_mpset = true;
					} else{
						val_assign(&mu, line);
						//calculate mubar
						this->calc_mubar();
						muset = true;
					}
			}
			if(	line.find("csw") != std::string::npos ||
					line.find("Csw") != std::string::npos) {
					if(	line.find("mp") != std::string::npos || 
							line.find("MP") != std::string::npos || 
							line.find("Mp") != std::string::npos ) {
						val_assign(&csw_mp, line);
						csw_mpset = true;
					} else {
						val_assign(&csw, line);
						cswset = true;
					}
			}
			
			if(	line.find("beta") != std::string::npos ||
					line.find("Beta") != std::string::npos) val_assign(&beta, line);

			if(	line.find("theta_fermion_spatial") != std::string::npos || 
					line.find("theta_spatial") != std::string::npos ||
					line.find("theta_s") != std::string::npos ||
					line.find("thetas") != std::string::npos ||
					line.find("ThetaS") != std::string::npos	) val_assign(&theta_fermion_spatial, line);
			
			if(	line.find("theta_fermion_temporal") != std::string::npos ||
					line.find("theta_temporal") != std::string::npos ||
					line.find("theta_t") != std::string::npos ||
					line.find("ThetaT") != std::string::npos ) val_assign(&theta_fermion_temporal, line);

			if(	line.find("cgmax") != std::string::npos || 
					line.find("CGmax") != std::string::npos ||
					line.find("Cgmax") != std::string::npos	) {
				if(	line.find("mp") != std::string::npos || 
						line.find("MP") != std::string::npos || 
						line.find("Mp") != std::string::npos ){
					val_assign(&cgmax_mp, line);
				} else {
					val_assign(&cgmax, line);
				}
			}

			if(	line.find("Solver") != std::string::npos ||
					line.find("solver") != std::string::npos ||
					line.find("SOLVER") != std::string::npos	) {
				if(	line.find("mp") != std::string::npos || 
						line.find("MP") != std::string::npos || 
						line.find("Mp") != std::string::npos ){
					cout << "here" << endl;
					solver_assign(&use_cg_mp, line, true);
				} else {
					solver_assign(&use_cg, line, false);
				}
			}
		
			if(line.find("sol_pr") != std::string::npos) val_assign(&solver_prec, line);
			if(line.find("force_pr") != std::string::npos) val_assign(&force_prec, line);

			if(line.find("writefrequency") != std::string::npos) val_assign(&writefrequency, line);
			if(line.find("savefrequency") != std::string::npos) val_assign(&savefrequency, line);
			if(line.find("saveconfigs") != std::string::npos) bool_assign(&saveconfigs, line);
			
			if(	line.find("prec") != std::string::npos ||
					line.find("Prec") != std::string::npos ) val_assign(&prec, line);

			if(line.find("readsource") != std::string::npos) startcond_assign(&startcondition, line);
			if(line.find("startcondition") != std::string::npos) startcond_assign(&startcondition, line);
			if(line.find("sourcefile") != std::string::npos) {
				val_assign(&sourcefile, line);
				sourcefilenumber_assign(&sourcefilenumber);
			}

			if(line.find("heatbathsteps") != std::string::npos) val_assign(&heatbathsteps, line);

			if(	line.find("thermalizationsteps") != std::string::npos ||
					line.find("thermsteps") != std::string::npos ||
					line.find("thermalization") != std::string::npos ) val_assign(&thermalizationsteps, line);

			if(	line.find("overrelaxsteps") != std::string::npos ||
					line.find("overrelax") != std::string::npos ||
					line.find("oversteps") != std::string::npos ) val_assign(&overrelaxsteps, line);

			if(line.find("iter_refresh") != std::string::npos) val_assign(&iter_refresh, line);

			if(line.find("num_dev") != std::string::npos) val_assign(&num_dev, line);

			if(	line.find("fermaction") 		!= std::string::npos	||
					line.find("fermionaction")	!= std::string::npos	||
					line.find("fermact") 				!= std::string::npos 	) {
				if(	line.find("mp") != std::string::npos || 
							line.find("MP") != std::string::npos || 
							line.find("Mp") != std::string::npos ) {
					fermact_assign(&fermact_mp, line);
				} else{
					fermact_assign(&fermact, line);
				}
			}
			
			if(	line.find("fermaction_mp") 		!= std::string::npos	||
					line.find("fermionaction_mp")	!= std::string::npos	||
					line.find("fermact_mp") 				!= std::string::npos 	) fermact_assign(&fermact_mp, line);

			if(	line.find("gaugeaction") 	!= std::string::npos ||
					line.find("gaugeact") 		!= std::string::npos ) {
				gaugeact_assign(&gaugeact, line, c1set);
			}
			if(line.find("c1") != std::string::npos) {
				val_assign(&c1, line);
				c1set = true;
				this->calc_c0_tlsym(this->get_c1());
			}
			
			if(line.find("hmcsteps") != std::string::npos) val_assign(&hmcsteps, line);
			if(	line.find("tau") != std::string::npos) val_assign(&tau, line);
			if(line.find("integrationsteps0") != std::string::npos) val_assign(&integrationsteps0, line);
			if(line.find("integrationsteps1") != std::string::npos) val_assign(&integrationsteps1, line);
			if(line.find("integrationsteps2") != std::string::npos) val_assign(&integrationsteps2, line);
			if(line.find("num_timescales") != std::string::npos) val_assign(&num_timescales, line);
			if(line.find("integrator0") != std::string::npos) integrator_assign(&integrator0, line);
			if(line.find("integrator1") != std::string::npos) integrator_assign(&integrator1, line);
			if(line.find("integrator2") != std::string::npos) integrator_assign(&integrator2, line);
			if(line.find("lambda0") != std::string::npos) val_assign(&lambda0, line);
			if(line.find("lambda1") != std::string::npos) val_assign(&lambda1, line);
			if(line.find("lambda2") != std::string::npos) val_assign(&lambda2, line);

			if(	line.find("evenodd") != std::string::npos	||
					line.find("even-odd") != std::string::npos ||
					line.find("even_odd") != std::string::npos ||
					line.find("even-odd-preconditioning") != std::string::npos ||
					line.find("use_eo") != std::string::npos ||
					line.find("use_evenodd") != std::string::npos ) bool_assign(&use_eo, line);
			
			if(	line.find("use_mp") != std::string::npos	||
					line.find("use_mass-preconditioning") != std::string::npos ||
					line.find("use_MP") != std::string::npos ||
					line.find("use_massprec") != std::string::npos ||
					line.find("use_mass-prec") != std::string::npos  ) bool_assign(&use_mp, line);

			if(	line.find("use_rec12") != std::string::npos ||
					line.find("REC12") != std::string::npos) bool_assign(&use_rec12, line);
			
			if(	line.find("use_gpu") != std::string::npos ||
					line.find("GPU") != std::string::npos  ) bool_assign(&use_gpu, line);

			if(	line.find("NS") != std::string::npos ||	
					line.find("NSPACE") != std::string::npos) val_assign(&nspace, line);
			if(	line.find("NT") != std::string::npos ||
					line.find("NTIME") != std::string::npos) val_assign(&ntime, line);

			if(	line.find("XI") != std::string::npos ||
					line.find("xi") != std::string::npos ||
					line.find("Xi") != std::string::npos) val_assign(&xi, line);
			if(	line.find("anisotropy") != std::string::npos ||
					line.find("Anisotropy") != std::string::npos ||
					line.find("use_aniso") != std::string::npos) bool_assign(&use_aniso, line);

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

			if(line.find("use_same_rnd_numbers") != std::string::npos) bool_assign(&use_same_rnd_numbers, line);
			if(line.find("taketimeofsol") != std::string::npos) {
				cout << line << endl;
				bool_assign(&profile_solver, line);
			}
		}

		//check for wrong settings in fermionaction
		if( ( (muset  == true) && (fermact != TWISTEDMASS) ) ||
		    ( (cswset == true) && (fermact != CLOVER     ) ) ||
		    ( (cswset == true) && (muset   == true       ) )  )
			throw Invalid_Fermact(fermact, muset, cswset);
		if( ( (mu_mpset  == true) && (fermact_mp != TWISTEDMASS) ) ||
		    ( (csw_mpset == true) && (fermact_mp != CLOVER     ) ) ||
		    ( (csw_mpset == true) && (mu_mpset   == true       ) )  )
			throw Invalid_Fermact(fermact_mp, mu_mpset, csw_mpset);
		//check for wrong settings in gaugeaction
		if( ( (use_rectangles  == true) && (gaugeact != TLSYM) ) ||
		    ( (use_rectangles  == false) && (gaugeact == TLSYM) ) )
			throw Invalid_Gaugeact();
		
		//check the read-in values against the compile time values
		this->set_settings_global();
		this->check_settings_global();
		this->set_gauge_norm_factors();

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
	} catch (Invalid_Gaugeact& ) {
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

	if(	value.find("TWISTEDMASS") != std::string::npos	||
			value.find("twistedmass") != std::string::npos	||
			value.find("Twistedmass") != std::string::npos	||
			value.find("TwistedMass") != std::string::npos	) {
		(*out) = TWISTEDMASS;
		return;
	}
	if(	value.find("clover") != std::string::npos	||
			value.find("CLOVER") != std::string::npos ||
			value.find("Clover") != std::string::npos ) {
		(*out) = CLOVER;
		return;
	}
	if(	value.find("WILSON") != std::string::npos ||
			value.find("Wilson") != std::string::npos ||
			value.find("wilson") != std::string::npos ||
			value.find("unimproved") != std::string::npos ) {
		(*out) = WILSON;
		return;
	}
	throw line;
	return;
}

void inputparameters::gaugeact_assign(int * out, std::string line, bool mu1set)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);

	if(	value.find("TLSYM") != std::string::npos	||
			value.find("tlsym") != std::string::npos	||
			value.find("Tlsym") != std::string::npos	||
			value.find("TLsym") != std::string::npos	) {
		(*out) = TLSYM;
		bool_assign(&use_rectangles, "1");
		if(!mu1set){
			c1 = c1_default_tlsym;
			calc_c0_tlsym(c1_default_tlsym);
		}
		return;
	}
	if(	value.find("IWASAKI") != std::string::npos	||
			value.find("iwasaki") != std::string::npos	||
			value.find("Iwasaki") != std::string::npos 	) {
		(*out) = IWASAKI;
		bool_assign(&use_rectangles, "1");
		if(!mu1set){
			c1 = c1_default_iwasaki;
			calc_c0_tlsym(c1_default_iwasaki);
		}
		return;
	}
	if(	value.find("DBW2") != std::string::npos	||
			value.find("dbw2") != std::string::npos	||
			value.find("Dbw2") != std::string::npos	) {
		(*out) = DBW2;
		bool_assign(&use_rectangles, "1");
		if(!mu1set){
			c1 = c1_default_dbw2;
			calc_c0_tlsym(c1_default_dbw2);
		}
		return;
	}
	if(	value.find("WILSON") != std::string::npos	||
			value.find("Wilson") != std::string::npos	||
			value.find("wilson") != std::string::npos	||
			value.find("unimproved") != std::string::npos	) {
		(*out) = WILSON;
		return;
	}
	throw line;
	return;
}

void inputparameters::solver_assign(bool * out, std::string line, bool mp)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);
	//LZ: note that the ordering of false/true is crucial here
	//    as any "bicgstab" hit will also be a "cg" hit...
	//    and any "bicgstab_save" hit will also be a "bicgcgstab" hit...
	cout << value << endl;
	if(	value.find("BiCGStab_save") != std::string::npos ||
			value.find("bicgstab_save") != std::string::npos ||	
			value.find("BICGSTAB_save") != std::string::npos ||	
			value.find("BICGSTAB_SAVE") != std::string::npos ) {
		(*out) = false;
		if(mp)
			this->use_bicgstab_save_mp = true;
		else
			this->use_bicgstab_save = true;
		return;
	}
	if(	value.find("BiCGStab") != std::string::npos ||
			value.find("bicgstab") != std::string::npos ||
			value.find("BICGSTAB") != std::string::npos ) {
		(*out) = false;
		return;
	}
	if(	value.find("Cg") != std::string::npos ||
			value.find("cg") != std::string::npos ||
			value.find("CG") != std::string::npos ) {
		(*out) = true;
	cout << "!" << endl;
		return;
	}
	cout << "throw" << endl;
	throw line;
	return;
}

void inputparameters::integrator_assign(int * out, std::string line)
{
	size_t pos = line.find("=");
	std::string value = line.substr(pos + 1);

	if(	value.find("LEAPFROG") != std::string::npos ||
			value.find("leapfrog") != std::string::npos) {
		(*out) = LEAPFROG;
		return;
	}
	if(	value.find("2MN") != std::string::npos ||
			value.find("2mn") != std::string::npos) {
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

	if(	value.find("1") != std::string::npos ||
			value.find("on") != std::string::npos ||
			value.find("ON") != std::string::npos ||
			value.find("yes") != std::string::npos ||
			value.find("TRUE") != std::string::npos ||
			value.find("true") != std::string::npos ||
			value.find("True") != std::string::npos ) {
		(*out) = true;
		return;
	}
	if(	value.find("0") != std::string::npos ||
			value.find("off") != std::string::npos ||
			value.find("OFF") != std::string::npos ||
			value.find("no") != std::string::npos ||
			value.find("false") != std::string::npos ||
			value.find("FALSE") != std::string::npos ||
			value.find("False") != std::string::npos ) {
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

hmc_float inputparameters::get_kappa_mp() const
{
	return kappa_mp;
}

void inputparameters::calc_c0_tlsym(hmc_float c1)
{
	c0 =  1. - 8.*c1;
	return;
}

hmc_float inputparameters::get_c0() const
{
	return c0;
}

hmc_float inputparameters::get_c1() const
{
	return c1;
}

void inputparameters::set_mubar_negative()
{
	mubar *= -1.;
}

void inputparameters::calc_mubar()
{
	mubar = 2.*kappa * mu;
}

void inputparameters::calc_mubar_mp()
{
	mubar_mp = 2.*kappa_mp * mu_mp;
}

hmc_float inputparameters::get_mubar() const
{
	return mubar;
}

hmc_float inputparameters::get_mubar_mp() const
{
	return mubar_mp;
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

hmc_float inputparameters::get_mu_mp() const
{
	return mu_mp;
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

hmc_float inputparameters::get_csw_mp() const
{
	return csw_mp;
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

hmc_float inputparameters::get_lambda(int number) const
{
	switch (number) {
		case 0:
			return lambda0;
		case 1:
			return lambda1;
		case 2:
			return lambda2;
		default:
			logger.fatal() << "Error requesting lambda for timescale " << number;
			logger.fatal() << "Aborting.";
			exit(EXIT_INPUTPARAMETERS);
	}
}

int inputparameters::get_cgmax() const
{
	return cgmax;
}

int inputparameters::get_cgmax_mp() const
{
	return cgmax_mp;
}

int inputparameters::get_ns() const
{
	return nspace;
}

int inputparameters::get_nt() const
{
	return ntime;
}

int inputparameters::get_xi() const
{
	int tmp;
	if (use_aniso == 0)
		tmp = 1;
	else
		tmp = xi;
	return tmp;
}

hmc_float inputparameters::get_xi_0() const
{
	hmc_float aniso = hmc_float (get_xi());
	hmc_float eta = (1.002503 * aniso * aniso * aniso + .39100 * aniso * aniso + 1.47130 * aniso - 0.19231) /
	                (aniso * aniso * aniso + 0.26287 * aniso * aniso + 1.59008 * aniso - 0.18224);
	return aniso / (1. + (1. - 1. / aniso) * eta / 6. * (1 - 0.55055 * 2 * NC / beta) / (1 - 0.77810 * 2 * NC / beta) * 2 * NC / beta );
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

int inputparameters::get_integrationsteps(int number) const
{
	switch (number) {
		case 0:
			return integrationsteps0;
		case 1:
			return integrationsteps1;
		case 2:
			return integrationsteps2;
		default:
			logger.fatal() << "Error requesting integrationsteps for timescale " << number;
			logger.fatal() << "Aborting.";
			exit(EXIT_INPUTPARAMETERS);
	}
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

int inputparameters::get_fermact_mp() const
{
	return fermact_mp;
}

int inputparameters::get_gaugeact() const
{
	return gaugeact;
}

int inputparameters::get_integrator(int which) const
{
	switch (which) {
		case 0:
			return integrator0;
		case 1:
			return integrator1;
		case 2:
			return integrator2;
		default:
		logger.fatal() << "cant make sense out of desired integrator number...\nAborting...";
		exit(EXIT_INPUTPARAMETERS);
	}	
}

int inputparameters::get_num_timescales() const
{
	return num_timescales;
}

bool inputparameters::get_saveconfigs() const
{
	return saveconfigs;
}

bool inputparameters::get_profile_solver() const
{
	return profile_solver;
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

bool inputparameters::get_use_mp() const
{
	return use_mp;
}

bool inputparameters::get_use_rec12() const
{
	return use_rec12;
}

bool inputparameters::get_use_rectangles() const
{
	return use_rectangles;
}

bool inputparameters::get_use_gpu() const
{
	return use_gpu;
}

bool inputparameters::get_use_aniso() const
{
	return use_aniso;
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

bool inputparameters::get_use_cg_mp() const
{
	return use_cg_mp;
}

bool inputparameters::get_use_bicgstab_save() const
{
	return use_bicgstab_save;
}

bool inputparameters::get_use_bicgstab_save_mp() const
{
	return use_bicgstab_save_mp;
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

bool inputparameters::get_use_same_rnd_numbers() const
{
	return use_same_rnd_numbers;
}

int inputparameters::get_mat_size() const
{
	return mat_size;
}
int inputparameters::get_float_size() const
{
	return float_size;
}
int inputparameters::get_flop_su3_su3() const
{
	return flop_su3_su3;
}
int inputparameters::get_flop_su3_su3vec() const
{
	return flop_su3_su3vec;
}
int inputparameters::get_flop_su3trace() const
{
	return flop_su3trace;
}
int inputparameters::get_flop_complex_mult() const
{
	return flop_complex_mult;
}
int inputparameters::get_flop_spinor_spinor() const
{
	return flop_spinor_spinor;
}
int inputparameters::get_flop_su3_spinor() const
{
	return flop_su3_spinor;
}
int inputparameters::get_flop_gamma_spinor() const
{
	return flop_gamma_spinor;
}
int inputparameters::get_flop_spinor_sqnorm() const
{
	return flop_spinor_sqnorm;
}

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

hmc_float inputparameters::get_force_prec() const
{
	return force_prec;
}


int inputparameters::get_iter_refresh() const
{
	return iter_refresh;
}

int inputparameters::get_iter_refresh_mp() const
{
	return iter_refresh_mp;
}

bool inputparameters::get_reversibility_check() const
{
	return reversibility_check;
}

int inputparameters::get_plaq_norm() const
{
	return plaq_norm;
}

int inputparameters::get_tplaq_norm() const
{
	return tplaq_norm;
}

int inputparameters::get_splaq_norm() const
{
	return splaq_norm;
}

int inputparameters::get_poly_norm() const
{
	return poly_norm;
}

int inputparameters::get_rect_norm() const
{
	return rect_norm;
}

void inputparameters::print_info_global() const
{
	logger.info() << "## Build based on commit: " << GIT_COMMIT_ID;
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
	if(this->get_use_aniso() == true) {
		logger.info() << "## USE ANISOTROPY: ON";
	} else {
		logger.info() << "## USE ANISOTROPY: OFF";
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
	*os  << "## Build based on commit: " << GIT_COMMIT_ID << endl;
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
	if(this->get_use_aniso() == true) {
		*os << "## USE ANISOTROPY: ON";
	} else {
		*os << "## USE ANISOTROPY: OFF";
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
	logger.info() << "## xi             = " << this->get_xi();
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
	*os  << "## xi             = " << this->get_xi();
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
	*os  << "## xi             = " << this->get_xi();
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
	logger.info() << "## xi             = " << this->get_xi();
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
	logger.info() << "## precision for inversions = " << this->get_solver_prec();
	if(this->get_use_eo() == true)
		logger.info() << "## Use even-odd preconditioning" ;
	if(this->get_use_eo() == false)
		logger.info() << "## Do NOT use even-odd preconditioning";
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
	if(this->get_use_cg() == false) {
		if(this->get_use_bicgstab_save() == false)
			logger.info() << "## Use BiCGStab for inversions";
		else
			logger.info() << "## Use BiCGStab-SAVE for inversions";
	}
	logger.info() << "## cgmax  = " << this->get_cgmax();
	logger.info() << "## iter_refresh  = " << this->get_iter_refresh();

	if(this->get_profile_solver() == true)
		logger.warn() << "## Profiling of solver activated. This may influence the overall performance time!";

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
	*os << "## precision for inversions = " << this->get_solver_prec() << endl;
	if(this->get_use_eo() == true)
		*os  << "## Use even-odd preconditioning" << endl;
	if(this->get_use_eo() == false)
		*os  << "## Do NOT use even-odd preconditioning" << endl;
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
	if(this->get_use_cg() == false) {
		if(this->get_use_bicgstab_save() == false)
			*os << "## Use BiCGStab for inversions" << endl;
		else
			*os << "## Use BiCGStab-SAVE for inversions" << endl;
	}
	*os << "## cgmax  = " << this->get_cgmax() << endl;
	*os << "## iter_refresh  = " << this->get_iter_refresh() << endl;

	if(this->get_profile_solver() == true)
		*os << "## Profiling of solver activated. This may influence the overall performance time!" << endl;

	//print extra warning if BC are set to default since this is a serious source of errors...
	if ( this->get_theta_fermion_spatial() == 0. && this->get_theta_fermion_temporal() == 0.) {
		*os << "\nNOTE: BCs have been set to periodic values by default!!\nTo change this use e.g. ThetaT/ThetaS in the input-file.\n";
	}
}

void inputparameters::print_info_gauge(ostream* os) const
{
	*os<< "## **********************************************************"<< endl;
	*os<< "## Gauge parameters:"<< endl;
	*os<< "##" << endl;
	if(this->get_gaugeact() == WILSON) {
		*os<<  "## gauge action: unimproved Wilson"<< endl;
	}
	if(this->get_gaugeact() == TLSYM) {
		*os<<  "## gauge action: tree level Symanzik"<< endl;
		*os<< "## c0  = " << this->get_c0()<< endl;
		*os<< "## c1  = " << this->get_c1()<< endl;
	}
}

void inputparameters::print_info_gauge() const
{
	logger.info() << "## **********************************************************";
	logger.info() << "## Gauge parameters:";
	logger.info() << "##" ;
	if(this->get_gaugeact() == WILSON) {
		logger.info() <<  "## gauge action: unimproved Wilson";
	}
	if(this->get_gaugeact() == TLSYM) {
		logger.info() <<  "## gauge action: tree level Symanzik";
		logger.info() << "## c0  = " << this->get_c0();
		logger.info() << "## c1  = " << this->get_c1();
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

void inputparameters::print_info_integrator(int number) const {
	string integrator_name;
	bool print_lambda = false;
	if(this->get_integrator(number) == LEAPFROG)
		integrator_name = "LEAPFROG";
	else if (this->get_integrator(number) == TWOMN) {
		integrator_name = "2MN";
		print_lambda = true;
	}
	else {
		logger.fatal() << "Fail in getting integrator information!";
		logger.fatal() << "Aborting...";
		exit(EXIT_INPUTPARAMETERS);
	}
	logger.info() << "## integrator" << number << " = " << integrator_name;
	logger.info() << "## integrationsteps" << number << " = " << this->get_integrationsteps(number);
	if(print_lambda) logger.info() << "## lambda" << number << " = " << get_lambda(number);
}

void inputparameters::print_info_integrator(ostream* os, int number) const {
	string integrator_name;
	bool print_lambda = false;
	if(this->get_integrator(number) == LEAPFROG)
		integrator_name = "LEAPFROG";
	else if (this->get_integrator(number) == TWOMN) {
		integrator_name = "2MN";
		print_lambda = true;
	}
	else {
		logger.fatal() << "Fail in getting integrator information!";
		logger.fatal() << "Aborting...";
		exit(EXIT_INPUTPARAMETERS);
	}
	*os<< "## integrator" << number << " = " << integrator_name << endl;
	*os<< "## integrationsteps" << number << " = " << this->get_integrationsteps(number) << endl;
	if(print_lambda) *os<< "## lambda" << number << " = " << get_lambda(number) << endl;
}

void inputparameters::print_info_hmc(char* progname) const
{

	logger.info() << "## Starting hmc program, executable name: " << progname ;
	this->print_info_global();
	this->print_info_fermion();
	this->print_info_gauge();
	logger.info() << "## **********************************************************";
	logger.info() << "## HMC parameters: " ;
	logger.info() << "##  ";
	logger.info() << "## tau  = " << this->get_tau();
	logger.info() << "## HMC steps  = " << this->get_hmcsteps();
	logger.info() << "## precision used in HMC-inversions = " << this->get_force_prec();
	logger.info() << "##  ";
	logger.info() << "## # Timescales  = " << this->get_num_timescales();
	//integrator infos
	for(int i = 0; i< this->get_num_timescales(); i++){
		print_info_integrator(i);
	}
	if(this->get_use_mp() == true){
		logger.info() << "##  ";
		logger.info() <<  "## use mass preconditioning:";
		if(this->get_fermact_mp() == WILSON) {
			logger.info() <<  "## mp action: unimproved Wilson";
			logger.info() << "## kappa_mp  = " << this->get_kappa_mp();
		}
		if(this->get_fermact_mp() == TWISTEDMASS) {
			logger.info() <<  "## mp action: twisted mass Wilson";
			logger.info() << "## kappa_mp  = " << this->get_kappa_mp();
			logger.info() << "## mu_mp     = " << this->get_mu_mp();
		}
		if(this->get_fermact_mp() == CLOVER) {
			logger.info() <<  "## mp action: clover Wilson";
			logger.info() << "## kappa_mp  = " << this->get_kappa_mp();
			logger.info() << "## csw_mp   = " << this->get_csw_mp();
		}
		logger.info() << "##" ;
		if(this->get_use_cg_mp() == true)
		logger.info() << "## Use CG-solver for mp inversions" ;
		if(this->get_use_cg_mp() == false) {
			if(this->get_use_bicgstab_save_mp() == false)
				logger.info() << "## Use BiCGStab for mp inversions";
			else
				logger.info() << "## Use BiCGStab-SAVE for mp inversions";
		}
		logger.info() << "## cgmax_mp  = " << this->get_cgmax_mp();
		logger.info() << "## iter_refresh_mp  = " << this->get_iter_refresh_mp();
		logger.info() << "##" ;
	}
	logger.info() << "## **********************************************************";
	return;
}

void inputparameters::print_info_hmc(char* progname, ostream* os) const
{
	*os << "## Starting hmc program, executable name: " << progname << endl;
	this->print_info_global(os);
	this->print_info_fermion(os);
	this->print_info_gauge(os);
	*os << "## **********************************************************" << endl;
	*os << "## HMC parameters: "  << '\n';
	*os << "##  " << '\n';
	*os << "## tau  = " << this->get_tau() << '\n';
	*os << "## HMC steps  = " << this->get_hmcsteps() << '\n';
	*os << "## precision used HMC-inversions = " << this->get_force_prec() << '\n';
	*os << "##  " << '\n';
	*os << "## # Timescales  = " << this->get_num_timescales() << '\n';
	//integrator infos
	for(int i = 0; i< this->get_num_timescales(); i++){
		print_info_integrator(os, i);
	}
	if(this->get_use_mp() == true){
		*os << "##  " << '\n';
		*os<<  "## use mass preconditioning:"  << '\n';
		if(this->get_fermact_mp() == WILSON) {
			*os<<  "## mp action: unimproved Wilson"  << '\n';
			*os<< "## kappa_mp  = " << this->get_kappa_mp()  << '\n';
		}
		if(this->get_fermact_mp() == TWISTEDMASS) {
			*os<<  "## mp action: twisted mass Wilson"  << '\n';
			*os<< "## kappa_mp  = " << this->get_kappa_mp()  << '\n';
			*os<< "## mu_mp     = " << this->get_mu_mp()  << '\n';
		}
		if(this->get_fermact_mp() == CLOVER) {
			*os<<  "## mp action: clover Wilson";
			*os<< "## kappa_mp  = " << this->get_kappa_mp()  << '\n';
			*os<< "## csw_mp   = " << this->get_csw_mp()  << '\n';
		}
		*os<< "##"  << endl;
		if(this->get_use_cg_mp() == true)
		*os<< "## Use CG-solver for mp inversions"  << '\n';
		if(this->get_use_cg_mp() == false) {
			if(this->get_use_bicgstab_save_mp() == false)
				*os<< "## Use BiCGStab for mp inversions"  << '\n';
			else
				*os<< "## Use BiCGStab-SAVE for mp inversions"  << '\n';
		}
		*os<< "## cgmax_mp  = " << this->get_cgmax_mp()  << '\n';
		*os<< "## iter_refresh_mp  = " << this->get_iter_refresh_mp()  << '\n';
		*os<< "##"  << '\n';
	}
	*os << "## **********************************************************" << '\n';
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
