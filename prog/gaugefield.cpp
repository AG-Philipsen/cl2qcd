#include "gaugefield.h"

hmc_error Gaugefield::init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer)
{
  //allocate memory for private gaugefield
	hmc_gaugefield* gftmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	set_gf(gftmp);


	set_parameters(input_parameters);

	init_gaugefield(timer);

	set_num_ocl_devices(numdevs);

	init_devices(devicetypes,timer);

	return HMC_SUCCESS;
}

hmc_error Gaugefield::init_devices(cl_device_type* devicetypes, usetimer* timer){
	if(get_num_ocl_devices() != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(get_num_ocl_devices() > 0) {
	  Opencl* dev_tmp = new Opencl[get_num_ocl_devices()];
	  set_devices(dev_tmp);
	}


	for(int n = 0; n < num_ocl_devices; n++) {
		cout << "init device #" << n << endl;
		get_devices()[n].init(devicetypes[n], local_work_size, global_work_size, timer, get_parameters());
		
	}
	return HMC_SUCCESS;
}


hmc_error Gaugefield::init_gaugefield(usetimer* timer)
{
	sourcefileparameters parameters_source;
	if((get_parameters())->get_startcondition() == START_FROM_SOURCE) {
		int err;
		timer->reset();
		//tmp gauge field
		hmc_float * gaugefield_tmp;
		gaugefield_tmp = (hmc_float*) malloc(sizeof(hmc_float) * NDIM * NC * NC * NTIME * VOLSPACE);
		err = parameters_source.readsourcefile(&(get_parameters()->sourcefile)[0], get_parameters()->get_prec(), &gaugefield_tmp);
		err = copy_gaugefield_from_ildg_format(get_gf(), gaugefield_tmp, parameters_source.num_entries_source);
		free(gaugefield_tmp);
		timer->add();
		if (err == 0) {
			print_info_source(&parameters_source);
		} else {
			printf("error in setting vals from source!!! check global settings!!!\n\n");
			return HMC_XMLERROR;
		}
	}
	if(get_parameters()->get_startcondition() == COLD_START) {
		timer->reset();
		set_gaugefield_cold(get_gf());
		timer->add();
	}
	if(get_parameters()->get_startcondition() == HOT_START) {
		timer->reset();
		set_gaugefield_cold(get_gf());
		timer->add();
	}
	return HMC_SUCCESS;
}



hmc_error Gaugefield::copy_gaugefield_to_devices(usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
  hmc_error err = get_devices()[0].copy_gaugefield_to_device(get_gf(), timer);
	return err;
}

hmc_error Gaugefield::sync_gaugefield(usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
  hmc_error err = get_devices()[0].get_gaugefield_from_device(get_gf(), timer);
	return err;
}

hmc_error Gaugefield::copy_rndarray_to_devices(hmc_rndarray host_rndarray,  usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
  hmc_error err = get_devices()[0].copy_rndarray_to_device(host_rndarray, timer);
	return err;
}

hmc_error Gaugefield::copy_rndarray_from_devices(hmc_rndarray rndarray, usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
  hmc_error err = get_devices()[0].copy_rndarray_from_device(rndarray, timer);
	return err;
}


void Gaugefield::print_info_source(sourcefileparameters* params)
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


hmc_error Gaugefield::save(int number)
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

	copy_gaugefield_to_ildg_format(gaugefield_buf, get_gf());

	hmc_float plaq = plaquette();

	write_gaugefield ( gaugefield_buf, gaugefield_buf_size , NSPACE, NSPACE, NSPACE, NTIME, get_parameters()->get_prec(), number, plaq, get_parameters()->get_beta(), get_parameters()->get_kappa(), get_parameters()->get_mu(), c2_rec, epsilonbar, mubar, version.c_str(), outputfile.c_str());

	free(gaugefield_buf);

	return HMC_SUCCESS;
}
void Gaugefield::print_gaugeobservables(usetimer * timer, usetimer * timer2)
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

void Gaugefield::print_gaugeobservables(usetimer * timer, usetimer * timer2, int iter)
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

void Gaugefield::print_gaugeobservables(usetimer * timer, usetimer * timer2, int iter, std::string filename)
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

void Gaugefield::print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter, std::string filename)
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

hmc_error Gaugefield::print_gaugeobservables_from_devices(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

  get_devices()[0].gaugeobservables(plaq, tplaq, splaq, pol, plaqtime, polytime);
	print_gaugeobservables(*plaq, *tplaq, *splaq, *pol, i, gaugeoutname);

	return HMC_SUCCESS;
}

hmc_error Gaugefield::print_gaugeobservables_from_devices(usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname)
{
	hmc_float plaq, tplaq, splaq;
	hmc_complex pol;

	hmc_error err = print_gaugeobservables_from_devices(&plaq, &tplaq, &splaq, &pol, plaqtime, polytime, i, gaugeoutname);

	return err;
}


hmc_float Gaugefield::plaquette()
{
	hmc_float tdummy;
	hmc_float sdummy;
	//LZ: calculation of tdummy and sdummy is unnecessary for this function, chance to speed up a little bit...
	return plaquette(&tdummy, &sdummy);
}

hmc_float Gaugefield::plaquette(hmc_float* tplaq, hmc_float* splaq)
{
	hmc_float plaq = 0;
	*tplaq = 0;
	*splaq = 0;
	for(int t = 0; t < NTIME; t++) {
		for(int n = 0; n < VOLSPACE; n++) {
			for(int mu = 0; mu < NDIM; mu++) {
				for(int nu = 0; nu < mu; nu++) {
					hmc_su3matrix prod;
					local_plaquette(get_gf(), &prod, n, t, mu, nu );
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


hmc_complex Gaugefield::polyakov()
{
	hmc_complex res;
	res.re = 0;
	res.im = 0;
	for(int n = 0; n < VOLSPACE; n++) {
		hmc_su3matrix prod;
		local_polyakov(get_gf(), &prod, n);
		hmc_complex tmpcomplex = trace_su3matrix(&prod);
		complexaccumulate(&res, &tmpcomplex);
	}
	res.re /= static_cast<hmc_float>(NC * VOLSPACE);
	res.im /= static_cast<hmc_float>(NC * VOLSPACE);
	return res;
}


hmc_complex Gaugefield::spatial_polyakov(int dir)
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
					get_su3matrix(&tmp, get_gf(), pos, t, dir);
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


hmc_error Gaugefield::heatbath(usetimer * const timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
  hmc_error err = get_devices()[0].run_heatbath(get_parameters()->get_beta(), timer);
	return err;
}

hmc_error Gaugefield::overrelax(usetimer * const timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

  hmc_error err = get_devices()[0].run_overrelax(get_parameters()->get_beta(), timer);
	return err;
}

hmc_error Gaugefield::heatbath(const int nheat, const int nover, usetimer * const timer_heat, usetimer * const timer_over)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(timer_heat);
	for(int i = 0; i < nover; i++) err |= overrelax(timer_over);
	return err;
}

hmc_error Gaugefield::heatbath(const int nheat, usetimer * const timer_heat)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(timer_heat);
	return err;
}

hmc_error Gaugefield::finalize()
{
  free(get_gf());
  if(num_ocl_devices > 0)
    delete [] get_devices();
  return HMC_SUCCESS;
}



hmc_error Gaugefield::set_gf (hmc_gaugefield * gf_val){
  gf = gf_val;
  return HMC_SUCCESS;
}

hmc_gaugefield * Gaugefield::get_gf ()
{
	return  gf;
}

hmc_error Gaugefield::set_devices (Opencl * devices_val){
  devices = devices_val;
  return HMC_SUCCESS;
}

Opencl * Gaugefield::get_devices ()
{
	return  devices;
}


hmc_error Gaugefield::set_num_ocl_devices (int num){
  num_ocl_devices = num;
  return HMC_SUCCESS;
}

int Gaugefield::get_num_ocl_devices ()
{
	return num_ocl_devices;
}


hmc_error Gaugefield::set_parameters (inputparameters * parameters_val){
  parameters = parameters_val; 
  return HMC_SUCCESS;
}

inputparameters * Gaugefield::get_parameters ()
{
	return  parameters;
}

	

	


	


