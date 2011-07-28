#include "gaugefield_k_hybrid.h"

hmc_error Gaugefield_k_hybrid::init_devices(cl_device_type* devicetypes, usetimer* timer){

	if(get_num_ocl_devices() != 1) {
	  logger.fatal() << "Gaugefiel_k needs exactly one device per type.";
	  exit(HMC_OCLERROR);
	}

	//this allocates two Opencl pointers in **devices
	alloc_devicetypes();


	Opencl_heatbath* dev_tmp_h = new Opencl_heatbath[get_num_ocl_devices()];
	set_devices(dev_tmp_h,DEV_HEATBATH);

	Opencl_k_hybrid* dev_tmp_k = new Opencl_k_hybrid[get_num_ocl_devices()];
	set_devices(dev_tmp_k,DEV_KAPPA);

	logger.trace()<<"init device for heatbath";
	get_devices_heatbath()[0].init(devicetypes[DEV_HEATBATH], timer, get_parameters());

	logger.trace()<<"init device for kappa calculation";
	get_devices_k()[0].init(devicetypes[DEV_KAPPA], timer, get_parameters());

	return HMC_SUCCESS;
}

hmc_error Gaugefield_k_hybrid::free_devices(){
  if(get_num_ocl_devices() > 0)
    delete [] get_devices_k();
  return HMC_SUCCESS;
}

hmc_float Gaugefield_k_hybrid::get_kappa_karsch (){
	return kappa_karsch_val;
}
	
hmc_float Gaugefield_k_hybrid::get_kappa_clover (){
	return kappa_clover_val;
}

hmc_error Gaugefield_k_hybrid::set_kappa_karsch (hmc_float in){
	kappa_karsch_val = in;
	return HMC_SUCCESS;
}

hmc_error Gaugefield_k_hybrid::set_kappa_clover (hmc_float in){
	kappa_clover_val = in;
	return HMC_SUCCESS;
  
}

hmc_error Gaugefield_k_hybrid::kappa_karsch_gpu (usetimer* timer){
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err;
 	err = get_devices_k()[0].run_kappa_karsch_gpu(get_parameters()->get_beta(), timer, &kappa_karsch_val);
	return err;
}

hmc_error Gaugefield_k_hybrid::kappa_clover_gpu (usetimer* timer){
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err;
 	err = get_devices_k()[0].run_kappa_clover_gpu(get_parameters()->get_beta(), timer, &kappa_clover_val);
	return err;
}

Opencl_k_hybrid * Gaugefield_k_hybrid::get_devices_k (){
  return  (Opencl_k_hybrid*)get_devices(DEV_KAPPA);
}

Opencl_heatbath * Gaugefield_k_hybrid::get_devices_heatbath (){
  return  (Opencl_heatbath*)get_devices(DEV_HEATBATH);
}


hmc_error Gaugefield_k_hybrid::copy_gaugefield_to_devices(usetimer* timer)
{
	hmc_error err = get_devices_k()[0].copy_gaugefield_to_device(get_sgf(), timer);
	err |= get_devices_heatbath()[0].copy_gaugefield_to_device(get_sgf(), timer);
	return err;
}


hmc_error Gaugefield_k_hybrid::copy_rndarray_to_devices(hmc_rndarray host_rndarray,  usetimer* timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err = get_devices_heatbath()[0].copy_rndarray_to_device(host_rndarray, timer);
	return err;
}

hmc_error Gaugefield_k_hybrid::copy_rndarray_from_devices(hmc_rndarray rndarray, usetimer* timer)
{
	hmc_error err = get_devices_heatbath()[0].copy_rndarray_from_device(rndarray, timer);
	return err;
}


hmc_error Gaugefield_k_hybrid::heatbath(usetimer * const timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
        hmc_error err = get_devices_heatbath()[0].run_heatbath(get_parameters()->get_beta(), timer);
	return err;
}

hmc_error Gaugefield_k_hybrid::overrelax(usetimer * const timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	hmc_error err = get_devices_heatbath()[0].run_overrelax(get_parameters()->get_beta(), timer);
	return err;
}

hmc_error Gaugefield_k_hybrid::heatbath(const int nheat, const int nover, usetimer * const timer_heat, usetimer * const timer_over)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(timer_heat);
	for(int i = 0; i < nover; i++) err |= overrelax(timer_over);
	return err;
}

hmc_error Gaugefield_k_hybrid::heatbath(const int nheat, usetimer * const timer_heat)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(timer_heat);
	return err;
}

hmc_error Gaugefield_k_hybrid::print_gaugeobservables_from_devices(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	get_devices_heatbath()[0].gaugeobservables(plaq, tplaq, splaq, pol, plaqtime, polytime);
	print_gaugeobservables(*plaq, *tplaq, *splaq, *pol, i, gaugeoutname);

	return HMC_SUCCESS;
}

hmc_error Gaugefield_k_hybrid::print_gaugeobservables_from_devices(usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname)
{
	hmc_float plaq, tplaq, splaq;
	hmc_complex pol;

	hmc_error err = print_gaugeobservables_from_devices(&plaq, &tplaq, &splaq, &pol, plaqtime, polytime, i, gaugeoutname);

	return err;
}

hmc_error Gaugefield_k_hybrid::print_TK(const string kappa_karsch_out,const string kappa_clover_out){
  std::fstream karschout;
  karschout.open(kappa_karsch_out.c_str(), std::ios::out | std::ios::app);
  if(!karschout.is_open()) exit(HMC_FILEERROR);
  //  karschout.width(8);
  //  karschout << iter;
  //  karschout << "\t";
  karschout.precision(15);
  karschout << get_kappa_karsch() << std::endl;
  karschout.close();

  std::fstream cloverout;
  cloverout.open(kappa_clover_out.c_str(), std::ios::out | std::ios::app);
  if(!cloverout.is_open()) exit(HMC_FILEERROR);
  //  cloverout.width(8);
  //  cloverout << iter;
  //  cloverout << "\t";
  cloverout.precision(15);
  cloverout << get_kappa_clover() << std::endl;
  cloverout.close();

  return HMC_SUCCESS;
}



hmc_error Gaugefield_k_hybrid::sync_gaugefield(usetimer* timer)
{
	hmc_error err = get_devices_heatbath()[0].get_gaugefield_from_device(get_sgf(), timer);
	err |= get_devices_k()[0].copy_gaugefield_to_device(get_sgf(), timer);

	return err;
}
