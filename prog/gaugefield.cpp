#include "gaugefield.h"

hmc_error gaugefield::init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer)
{
	gf = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));

	parameters = input_parameters;

	init_gaugefield(timer);


	num_ocl_devices = numdevs;

	if(numdevs != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(num_ocl_devices > 0)
		devices = new opencl[num_ocl_devices];

	for(int n = 0; n < num_ocl_devices; n++) {
		cout << "init device #" << n << endl;
		devices[n].init(devicetypes[n], local_work_size, global_work_size, timer, parameters);
	}

	return HMC_SUCCESS;
}

hmc_error gaugefield::init_gaugefield(usetimer* timer)
{
	sourcefileparameters parameters_source;
	if(parameters->get_startcondition() == START_FROM_SOURCE) {
		int err;
		timer->reset();
		//tmp gauge field
		hmc_float * gaugefield_tmp;
		gaugefield_tmp = (hmc_float*) malloc(sizeof(hmc_float) * NDIM * NC * NC * NTIME * VOLSPACE);
		err = parameters_source.readsourcefile(&(parameters->sourcefile)[0], parameters->get_prec(), &gaugefield_tmp);
		err = copy_gaugefield_from_ildg_format(gf, gaugefield_tmp, parameters_source.num_entries_source);
		free(gaugefield_tmp);
		timer->add();
		if (err == 0) {
			print_info_source(&parameters_source);
		} else {
			printf("error in setting vals from source!!! check global settings!!!\n\n");
			return HMC_XMLERROR;
		}
	}
	if(parameters->get_startcondition() == COLD_START) {
		timer->reset();
		set_gaugefield_cold(gf);
		timer->add();
	}
	if(parameters->get_startcondition() == HOT_START) {
		timer->reset();
		set_gaugefield_cold(gf);
		timer->add();
	}
	return HMC_SUCCESS;
}

void gaugefield::print_info_source(sourcefileparameters* params)
{
	printf("**********************************************************\n");
	printf("Sourcefile parameters: (list not complete)\n");
	printf("field:  %s\n", params->field_source);
	printf("LX:  \t%d\n", params->lx_source);
	printf("LY:  \t%d\n", params->ly_source);
	printf("LZ:  \t%d\n", params->lz_source);
	printf("LT:  \t%d\n", params->lt_source);
	printf("entries: %d\n", params->num_entries_source);
	printf("beta:  \t%f\n", params->beta_source);
	printf("mu:  \t%f\n", params->mu_source);
	printf("kappa:  %f\n", params->kappa_source);
	printf("mubar:  %f\n", params->mubar_source);
	printf("plaq: \t%f\n", params->plaquettevalue_source);
	printf("**********************************************************\n");
	printf("\n");
	return;
}


hmc_error gaugefield::save(int number)
{

	ildg_gaugefield * gaugefield_buf;
	gaugefield_buf = (ildg_gaugefield*) malloc(sizeof(ildg_gaugefield));
	int gaugefield_buf_size = sizeof(ildg_gaugefield) / sizeof(hmc_float);

	//these are not yet used...
	hmc_float c2_rec = 0, epsilonbar = 0, mubar = 0;
	//LZ: generalize the following to larger numbers, if necessary...
	stringstream strnumber;
	strnumber.fill('0');
	strnumber.width(5);
	strnumber << right << number;
	stringstream outfilename;
	outfilename << "conf." << strnumber.str();
	string outputfile = outfilename.str();

	copy_gaugefield_to_ildg_format(gaugefield_buf, gf);

	hmc_float plaq = plaquette();

	write_gaugefield ( gaugefield_buf, gaugefield_buf_size , NSPACE, NSPACE, NSPACE, NSPACE, parameters->get_prec(), number, plaq, parameters->get_beta(), parameters->get_kappa(), parameters->get_mu(), c2_rec, epsilonbar, mubar, version.c_str(), outputfile.c_str());

	free(gaugefield_buf);

	return HMC_SUCCESS;
}

hmc_error gaugefield::copy_gaugefield_to_devices(usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	hmc_error err = devices[0].copy_gaugefield_to_device(gf, timer);
	return err;
}

hmc_error gaugefield::sync_gaugefield(usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	hmc_error err = devices[0].get_gaugefield_from_device(gf, timer);
	return err;
}

hmc_error gaugefield::copy_rndarray_to_devices(hmc_rndarray host_rndarray,  usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err = devices[0].copy_rndarray_to_device(host_rndarray, timer);
	return err;
}

hmc_error gaugefield::copy_rndarray_from_devices(hmc_rndarray rndarray, usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err = devices[0].copy_rndarray_from_device(rndarray, timer);
	return err;
}

void gaugefield::print_gaugeobservables(usetimer * timer, usetimer * timer2)
{
	timer->reset();
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	timer->add();
	timer2->reset();
	hmc_complex pol = polyakov();
	printf("%f\t%f\t%f\t%f\t%f\t%f\n", plaq, tplaq, splaq, pol.re, pol.im, sqrt(pol.re * pol.re + pol.im * pol.im));
	timer2->add();
	return;
}

void gaugefield::print_gaugeobservables(usetimer * timer, usetimer * timer2, int iter)
{
	timer->reset();
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	timer->add();
	timer2->reset();
	hmc_complex pol = polyakov();
	printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\n", iter, plaq, tplaq, splaq, pol.re, pol.im, sqrt(pol.re * pol.re + pol.im * pol.im));
	timer2->add();
	return;
}

void gaugefield::print_gaugeobservables(usetimer * timer, usetimer * timer2, int iter, std::string filename)
{
	timer->reset();
	hmc_float tplaq = 0;
	hmc_float splaq = 0;
	hmc_float plaq = plaquette(&tplaq, &splaq);
	timer->add();
	timer2->reset();
	hmc_complex pol = polyakov();
	timer2->add();
	//printf("%d\t%f\t%f\t%f\t%f\t%f\n",iter,plaq,tplaq,splaq,pol.re,pol.im);
	std::fstream gaugeout;
	gaugeout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!gaugeout.is_open()) exit(HMC_FILEERROR);
	gaugeout.width(8);
	gaugeout << iter;
	gaugeout << "\t";
	gaugeout.precision(15);
	gaugeout << plaq << "\t" << tplaq << "\t" << splaq << "\t" << pol.re << "\t" << pol.im << "\t" << sqrt(pol.re * pol.re + pol.im * pol.im) << std::endl;
	gaugeout.close();
	return;
}

void gaugefield::print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter, std::string filename)
{
	//printf("%d\t%f\t%f\t%f\t%f\t%f\n",iter,plaq,tplaq,splaq,pol.re,pol.im);
	std::fstream gaugeout;
	gaugeout.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!gaugeout.is_open()) exit(HMC_FILEERROR);
	gaugeout.width(8);
	gaugeout << iter;
	gaugeout << "\t";
	gaugeout.precision(15);
	gaugeout << plaq << "\t" << tplaq << "\t" << splaq << "\t" << pol.re << "\t" << pol.im << "\t" << sqrt(pol.re * pol.re + pol.im * pol.im) << std::endl;
	gaugeout.close();
	return;
}

hmc_float gaugefield::plaquette()
{
	hmc_float tdummy;
	hmc_float sdummy;
	//LZ: calculation of tdummy and sdummy is unnecessary for this function, chance to speed up a little bit...
	return plaquette(&tdummy, &sdummy);
}

hmc_float gaugefield::plaquette(hmc_float* tplaq, hmc_float* splaq)
{
	hmc_float plaq = 0;
	*tplaq = 0;
	*splaq = 0;
	for(int t = 0; t < NTIME; t++) {
		for(int n = 0; n < VOLSPACE; n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				for(int nu = 0; nu < mu; nu++) {
					hmc_su3matrix prod;
					local_plaquette(gf, &prod, n, t, mu, nu );
					hmc_float tmpfloat = trace_su3matrix(&prod).re;
					plaq += tmpfloat;
					if(mu == 0 || nu == 0) {
						*tplaq += tmpfloat;
					} else {
						*splaq += tmpfloat;
					}
				}
			}
		}
	}
	*tplaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1));
	*splaq /= static_cast<hmc_float>(VOL4D * NC * (NDIM - 1) * (NDIM - 2)) / 2. ;
	return plaq * 2.0 / static_cast<hmc_float>(VOL4D * NDIM * (NDIM - 1) * NC);
}


hmc_complex gaugefield::polyakov()
{
	hmc_complex res;
	res.re = 0;
	res.im = 0;
	for(int n = 0; n < VOLSPACE; n++) {
		hmc_su3matrix prod;
		local_polyakov(gf, &prod, n);
		hmc_complex tmpcomplex = trace_su3matrix(&prod);
		complexaccumulate(&res, &tmpcomplex);
	}
	res.re /= static_cast<hmc_float>(NC * VOLSPACE);
	res.im /= static_cast<hmc_float>(NC * VOLSPACE);
	return res;
}


hmc_complex gaugefield::spatial_polyakov(int dir)
{
	//assuming dir=1,2, or 3
	hmc_complex res;
	res.re = 0;
	res.im = 0;
	for(int x1 = 0; x1 < NSPACE; x1++) {
		for(int x2 = 0; x2 < NSPACE; x2++) {
			for(int t = 0; t < NTIME; t++) {
				hmc_su3matrix prod;
				unit_su3matrix(&prod);
				for(int xpol = 0; xpol < NSPACE; xpol++) {
					hmc_su3matrix tmp;
					int coord[NDIM];
					coord[0] = t;
					coord[dir] = xpol;
					int next = (dir % (NDIM - 1)) + 1;
					coord[next] = x1;
					int nnext = (next % (NDIM - 1)) + 1;
					coord[nnext] = x2;
					int pos = get_nspace(coord);
					get_su3matrix(&tmp, gf, pos, t, dir);
					accumulate_su3matrix_prod(&prod, &tmp);
				}
				hmc_complex tmpcomplex = trace_su3matrix(&prod);
				complexaccumulate(&res, &tmpcomplex);
			}
		}
	}
	res.re /= static_cast<hmc_float>(NC * NSPACE * NSPACE * NTIME);
	res.im /= static_cast<hmc_float>(NC * NSPACE * NSPACE * NTIME);
	return res;
}

hmc_error gaugefield::print_gaugeobservables_from_devices(const size_t local_work_size, const size_t global_work_size, hmc_float* plaq, hmc_float* tplaq, hmc_float* splaq, hmc_complex* pol, usetimer* plaqtime, usetimer* polytime, int i, string gaugeoutname)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	devices[0].gaugeobservables(local_work_size, global_work_size, plaq, tplaq, splaq, pol, plaqtime, polytime);
	print_gaugeobservables(*plaq, *tplaq, *splaq, *pol, i, gaugeoutname);

	return HMC_SUCCESS;
}

hmc_error gaugefield::print_gaugeobservables_from_devices(const size_t local_work_size, const size_t global_work_size, usetimer* plaqtime, usetimer* polytime, int i, string gaugeoutname)
{
	hmc_float plaq, tplaq, splaq;
	hmc_complex pol;

	hmc_error err = print_gaugeobservables_from_devices(local_work_size, global_work_size, &plaq, &tplaq, &splaq, &pol, plaqtime, polytime, i, gaugeoutname);

	return err;
}

hmc_error gaugefield::heatbath(const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err = devices[0].run_heatbath(parameters->get_beta(), local_work_size, global_work_size, timer);
	return err;
}

hmc_error gaugefield::overrelax(const size_t local_work_size, const size_t global_work_size, usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	hmc_error err = devices[0].run_overrelax(parameters->get_beta(), local_work_size, global_work_size, timer);
	return err;
}

hmc_error gaugefield::heatbath(const size_t local_work_size, const size_t global_work_size, int nheat, int nover, usetimer* timer_heat, usetimer* timer_over)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(local_work_size, global_work_size, timer_heat);
	for(int i = 0; i < nover; i++) err |= overrelax(local_work_size, global_work_size, timer_over);
	return err;
}

hmc_error gaugefield::heatbath(const size_t local_work_size, const size_t global_work_size, int nheat, usetimer* timer_heat)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(local_work_size, global_work_size, timer_heat);
	return err;
}

hmc_error gaugefield::finalize()
{
	free(gf);
	if(num_ocl_devices > 0)
		delete [] devices;
	return HMC_SUCCESS;
}


void gaugefield::kappa_karsch (hmc_float & kappa)
{
  //Initialize kappa
  kappa = 0.0;
  //Initializing beta
  const hmc_float beta = parameters->get_beta();
  //Compute diagonal spatial components of the energy-momentum-tensor
  hmc_float tdiag_11 [VOL4D];
  hmc_float tdiag_22 [VOL4D];
  hmc_float tdiag_33 [VOL4D];
  //a = 2T_11 - T_22 _ T_33
  hmc_float a [VOL4D];
  //b = 2T_22 - T_11 _ T_33
  hmc_float b [VOL4D];
  //c = 2T_33 - T_22 _ T_11
  hmc_float c [VOL4D];
    
  for (int t=0; t<NTIME; t++){
	for (int n=0; n<VOLSPACE; n++){
	      //Compute required plaquettes
	      hmc_su3matrix temp;
	      
	      local_plaquette(gf, & temp, n, t, 1, 0);
	      //faster, to take real first and then trace, but no method
	      hmc_float plaq_10 = trace_su3matrix(&temp).re;
	      local_plaquette(gf, & temp, n, t, 2, 0);
	      hmc_float plaq_20 = trace_su3matrix(&temp).re;
	      local_plaquette(gf, & temp, n, t, 3, 0);
	      hmc_float plaq_30 = trace_su3matrix(&temp).re;
	      local_plaquette(gf, & temp, n, t, 1, 2);
	      hmc_float plaq_12 = trace_su3matrix(&temp).re;
	      local_plaquette(gf, & temp, n, t, 1, 3);
	      hmc_float plaq_13 = trace_su3matrix(&temp).re;
	      local_plaquette(gf, & temp, n, t, 3, 2);
	      hmc_float plaq_32 = trace_su3matrix(&temp).re;

	      int point = n + VOLSPACE * t;
	      
	      tdiag_11 [point] = plaq_10 + plaq_12 + plaq_13 - plaq_20 - plaq_30 - plaq_32;
	      tdiag_22 [point] = plaq_20 + plaq_12 + plaq_32 - plaq_10 - plaq_30 - plaq_13;
	      tdiag_33 [point] = plaq_30 + plaq_13 + plaq_32 - plaq_10 - plaq_20 - plaq_12;

	      a[point] = 2.0*tdiag_11[point] - tdiag_22[point] - tdiag_33[point];
	      b[point] = 2.0*tdiag_22[point] - tdiag_11[point] - tdiag_33[point];
	      c[point] = 2.0*tdiag_33[point] - tdiag_22[point] - tdiag_11[point];
  }}

  hmc_float deltak = 2*PI/NSPACE;
  hmc_float result = 0.0;

  
  for (int x_3=0; x_3<NSPACE; x_3++){
    for (int y_3=0; y_3<x_3; y_3++){
      hmc_float factor = 1.0 - cos(deltak*hmc_float(x_3-y_3));
      for (int x_t=0; x_t<NTIME; x_t++){
	for (int y_t=0; y_t<NTIME; y_t++){
	  for (int x_1=0; x_1<NSPACE; x_1++){
	    for (int y_1=0; y_1<NSPACE; y_1++){
	      for (int x_2=0; x_2<NSPACE; x_2++){
	        for (int y_2=0; y_2<NSPACE; y_2++){
		  
		  int coord_x[NDIM];
		  coord_x[0]= x_t;
		  coord_x[1]= x_1;
		  coord_x[2]= x_2;
		  coord_x[3]= x_3;
		  //new method get_tnspace which gives n_x+VOLSPACE*x_t
		  int n_x = get_nspace (coord_x);
		  int point_x = n_x + VOLSPACE*x_t;
		  int coord_y[NDIM];
		  coord_y[0]= y_t;
		  coord_y[1]= y_1;
		  coord_y[2]= y_2;
		  coord_y[3]= y_3;
		  int n_y = get_nspace (coord_y);
		  int point_y = n_y + VOLSPACE*y_t;

		  result += factor * ( tdiag_11 [point_x] * a[point_y]
		                     + tdiag_22 [point_x] * b[point_y]
		                     + tdiag_33 [point_x] * c[point_y]);
  }}}}}}}}

  //Correlator, 2 by Def, 2 by T_12+T_21  3 by T_12+T_13+T_23 -->/12,
  //Volume for y /VOL4D, /2/pi^2*L_z^2 for derivation, *beta^2 / Nc^2 for T_munu, *2 for y_3<x_3
  
  hmc_float norm = hmc_float (NSPACE* NSPACE) / hmc_float (VOL4D) / hmc_float(NC * NC) /12.0 /PI /PI  * beta * beta;

  kappa = norm * result;
}

void gaugefield::kappa_clover (hmc_float & kappa){
  //Initializing result
  kappa = .0;
  //Initializing beta
  const hmc_float beta = parameters->get_beta();
  //Energy-momentum-tensor in clover-discretization
  hmc_float t_12 [VOL4D];
  hmc_float t_13 [VOL4D];
  hmc_float t_23 [VOL4D];
  
  for (int t=0; t<NTIME; t++){
    for (int n=0; n<VOLSPACE; n++){
      //Compute required plaquettes
      hmc_3x3matrix Q_22;
      local_Q_plaquette(&Q_22, gf, n, t, 2, 2);
      hmc_3x3matrix Q_10;
      local_Q_plaquette(&Q_10, gf, n, t, 1, 0);
      hmc_3x3matrix Q_20;
      local_Q_plaquette(&Q_20, gf, n, t, 2, 0);
      hmc_3x3matrix Q_02;
      adjoint_3x3matrix (&Q_02, &Q_20);
// 	      local_Q_plaquette(gf, &Q_02, n, t, 0, 2);
      hmc_3x3matrix Q_21;
      local_Q_plaquette(&Q_21,gf, n, t, 2, 1);
      hmc_3x3matrix Q_12;
      adjoint_3x3matrix (&Q_12, &Q_21);
// 	      local_Q_plaquette(gf, &Q_12, n, t, 1, 2);
      hmc_3x3matrix Q_03;
      local_Q_plaquette(&Q_03, gf, n, t, 0, 3);
      hmc_3x3matrix Q_30;
      adjoint_3x3matrix (&Q_30, &Q_03);
      hmc_3x3matrix Q_13;
      local_Q_plaquette( &Q_13, gf, n, t, 1, 3);
      hmc_3x3matrix Q_31;
      adjoint_3x3matrix (&Q_31, &Q_13);
      hmc_3x3matrix Q_23;
      local_Q_plaquette( &Q_23, gf, n, t, 2, 3);
      hmc_3x3matrix Q_32;
      adjoint_3x3matrix (&Q_32, &Q_23);
// 	      local_Q_plaquette(gf, &Q_32, n, t, 3, 2);
      hmc_3x3matrix Q_11;
      local_Q_plaquette( &Q_11, gf, n, t, 1, 1);

      int point = n + VOLSPACE * t;
      
      hmc_3x3matrix tmp;
      hmc_complex tmp_cmp;
	      
      //T_12
      //alpha=0
      subtract_3x3matrix (&tmp, &Q_20, &Q_02);
      multiply_3x3matrix (&tmp, &Q_10, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_12 [point] = tmp_cmp.re;
      //alpha=1
      subtract_3x3matrix (&tmp, &Q_21, &Q_12);
      multiply_3x3matrix (&tmp, &Q_11, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_12 [point] += tmp_cmp.re;
      //alpha=2, vanishes
// 	      subtract_3x3matrix (&tmp, &Q_22, &Q_22);
// 	      multiply_3x3matrix (&tmp, &Q_12, &tmp);
// 	      trace_3x3matrix (tmp_cmp, &tmp);
// 	      t_12 [point] += tmp_cmp.re;
      //alpha=3
      subtract_3x3matrix (&tmp, &Q_23, &Q_32);
      multiply_3x3matrix (&tmp, &Q_13, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_12 [point] += tmp_cmp.re;

      //T_13
      //alpha=0
      subtract_3x3matrix (&tmp, &Q_30, &Q_03);
      multiply_3x3matrix (&tmp, &Q_10, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_13 [point] = tmp_cmp.re;
      //alpha=1
      subtract_3x3matrix (&tmp, &Q_31, &Q_13);
      multiply_3x3matrix (&tmp, &Q_11, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_13 [point] += tmp_cmp.re;
      //alpha=2, vanishes
      subtract_3x3matrix (&tmp, &Q_32, &Q_23);
      multiply_3x3matrix (&tmp, &Q_12, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_13 [point] += tmp_cmp.re;
      //alpha=3, vanishes
// 	      subtract_3x3matrix (&tmp, &Q_33, &Q_33);
// 	      multiply_3x3matrix (&tmp, &Q_13, &tmp);
// 	      trace_3x3matrix (&tmp_cmp, &tmp);
// 	      t_13 [point] += tmp_cmp.re;

      //T_23
      //alpha=0
      subtract_3x3matrix (&tmp, &Q_30, &Q_03);
      multiply_3x3matrix (&tmp, &Q_20, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_23 [point] = tmp_cmp.re;
      //alpha=1
      subtract_3x3matrix (&tmp, &Q_31, &Q_13);
      multiply_3x3matrix (&tmp, &Q_21, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_23 [point] += tmp_cmp.re;
      //alpha=2, vanishes
      subtract_3x3matrix (&tmp, &Q_32, &Q_23);
      multiply_3x3matrix (&tmp, &Q_22, &tmp);
      trace_3x3matrix (&tmp_cmp, &tmp);
      t_23 [point] += tmp_cmp.re;
      //alpha=3, vanishes
// 	      subtract_3x3matrix (&tmp, &Q_33, &Q_33);
// 	      multiply_3x3matrix (&tmp, &Q_23, &tmp);
// 	      trace_3x3matrix (&tmp_cmp, &tmp);
// 	      t_23 [point] += tmp_cmp.re;
      }
    }
    
    //Momentum
    const hmc_float deltak = 2.0*PI / NSPACE;
    hmc_float result = 0.0;

  
    for (int x_3=0; x_3<NSPACE; x_3++){
      for (int y_3=0; y_3<x_3; y_3++){
	hmc_float factor = 1.0 - cos(deltak*hmc_float(x_3-y_3));
	for (int x_t=0; x_t<NTIME; x_t++){
	  for (int y_t=0; y_t<NTIME; y_t++){
	    for (int x_1=0; x_1<NSPACE; x_1++){
	      for (int y_1=0; y_1<NSPACE; y_1++){
		for (int x_2=0; x_2<NSPACE; x_2++){
		  for (int y_2=0; y_2<NSPACE; y_2++){
		  
		    int coord_x[NDIM];
		    coord_x[0]= x_t;
		    coord_x[1]= x_1;
		    coord_x[2]= x_2;
		    coord_x[3]= x_3;
		    //new method get_tnspace which gives n_x+VOLSPACE*x_t
		    int n_x = get_nspace (coord_x);
		    int point_x = n_x + VOLSPACE*x_t;
		    int coord_y[NDIM];
		    coord_y[0]= y_t;
		    coord_y[1]= y_1;
		    coord_y[2]= y_2;
		    coord_y[3]= y_3;
		    int n_y = get_nspace (coord_y);
		    int point_y = n_y + VOLSPACE*y_t;

		    //(T_12(x) T_12(y) + T_21(x) T_21(y) + T_13(x) T_13(y)) * factor
		    result += factor* ( t_12[point_x]*t_12[point_y]
				      + t_13[point_x]*t_13[point_y]
				      + t_23[point_x]*t_23[point_y]);
  }}}}}}}}
  //Normalization
  // / (-96)^2=/ 9216, * beta^2, /3 für T_12T_12+T_13T_13+T_23T_23, /VOL4D für Summe y, /2 /pi^2 * L_z^2 für Ableitung * 2 für forschleife
  hmc_float norm = hmc_float (NSPACE* NSPACE) / hmc_float (VOL4D) /PI /PI  * beta * beta /9216. /3.;
    
  kappa = norm * result * beta;
}

