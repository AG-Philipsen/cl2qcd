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
