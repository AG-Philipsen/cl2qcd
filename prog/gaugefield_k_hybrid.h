/** @file
 * Provides a class for gauge fields, including calculation of transport coefficients
 *
 */
#ifndef _GAUGEFIELDKHYBRIDH_
#define _GAUGEFIELDKHYBRIDH_

#include "gaugefield.h"
#include "opencl_k_hybrid.h"
#include "opencl_heatbath.h"

#define DEV_KAPPA 1
#define DEV_HEATBATH 0

class Gaugefield_k_hybrid : public Gaugefield  {

  public:
    
  /**
   * Initializes the devices, to be called by init()
   * @return Error code as defined in hmcerrs.h
   * @param devicetypes array of cl_device_type handles
   * @param[in,out] timer timer for initialization
   */
  virtual hmc_error init_devices(cl_device_type* devicetypes, int nstates);

	/**
	 * Free device, called by finalize
	 */
	virtual hmc_error free_devices();
	
	/**
	 * Returns the transport coefficient kappa computed by Karsch&Wyld's method
	 * @return Result for the transport coefficient kappa
	 */	
	hmc_float get_kappa_karsch ();
	
	/**
	 * Returns the transport coefficient kappa computed by Clover method
	 * @return Result for the transport coefficient kappa
	 */	
	hmc_float get_kappa_clover ();
	
	/**
	 * Set the transport coefficient kappa computed by Karsch&Wyld's method
	 * @param[in] in Result for the transport coefficient kappa
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error set_kappa_karsch (hmc_float in);
	
	/**
	 * Set the transport coefficient kappa computed by Clover method
	 * @param[in] in Result for the transport coefficient kappa
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error set_kappa_clover (hmc_float in);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld on GPU
	 * @param[in,out] timer time for measurement
	 * @return Error code as defined in hmcerrs.h
	 */	
	 hmc_error kappa_karsch_gpu (usetimer* timer);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld  on GPU
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error kappa_clover_gpu (usetimer* timer);
  
	/**
	 * Returns private member opencl_k * devices
	 * @return devices of type opencl_k
	 * @todo LZ: CHECK IF THIS MAKES SENSE AT ALL!!!!
	 */
 	Opencl_k_hybrid * get_devices_k ();

	/**
	 * Returns private member opencl_heatbath * devices
	 * @return devices of type opencl_k
	 * @todo LZ: CHECK IF THIS MAKES SENSE AT ALL!!!!
	 */
 	Opencl_heatbath * get_devices_heatbath ();

	/**
	 * Copy gaugefield to devices.
	 */
	virtual hmc_error copy_gaugefield_to_devices();

	/**
	 * Perform a number of heatbath and (afterwards) overrelaxation steps.
	 * @param[in] nheat number of heatbath steps
	 * @param[in] nover number of overrelaxation steps
	 * @param[in,out] timer_heat time for heatbath steps
	 * @param[in,out] timer_over time for overrelaxation steps
	 */
	hmc_error heatbath(const int nheat, const int nover, usetimer * const timer_heat, usetimer * const timer_over);
	/**
	 * Perform a number of heatbath steps.
	 * @param[in] nheat number of heatbath steps
	 * @param[in,out] timer_heat time for heatbath steps
	 */
	hmc_error heatbath(const int nheat, usetimer * const timer_heat);
	/**
	 * Perform one heatbath step.
	 * @param[in,out] timer time for heatbath step
	 */
	hmc_error heatbath(usetimer * const timer);
	/**
	 * Perform one overrelaxation step.
	 * @param[in,out] timer time for overrelaxation step
	 */
	hmc_error overrelax(usetimer * const timer);

	/**
	 * Print gauge observables calculated on device, add iteration number, return values to program.
	 * @param[in,out] plaq pointer to plaquette value
	 * @param[in,out] tplaq pointer to timelike plaquette value
	 * @param[in,out] splaq pointer to spatial plaquette value
	 * @param[in,out] pol pointer to Polyakov loop value
	 * @param[in,out] plaqtime time to calculate plaquette
	 * @param[in,out] polytime time to calculate Polyakov loop
	 * @param[in] i integer number that accompanies output
	 * @param[in] gaugeoutname name of output file
	 * @return Error code as defined in hmcerrs.h
	 */
	virtual hmc_error print_gaugeobservables_from_devices(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname);

	/**
	 * Print transport coefficients to files.
	 * @param[in] kappa_karsch_out filename for Karsch TK
	 * @param[in] kappa_clover_out filename for clover TK
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error print_TK(const string kappa_karsch_out,const string kappa_clover_out);


	/**
	 * Print gauge observables calculated on device, add iteration number.
	 * @param[in,out] plaqtime time to calculate plaquette
	 * @param[in,out] polytime time to calculate Polyakov loop
	 * @param[in] i integer number that accompanies output
	 * @param[in] gaugeoutname name of output file
	 * @return Error code as defined in hmcerrs.h
	 */
	virtual hmc_error print_gaugeobservables_from_devices(usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname);

	/**
	 * Communicate gaugefield from heatbath device to everyone else.
	 */
	virtual hmc_error sync_gaugefield();

  private:
  hmc_float kappa_karsch_val;
  hmc_float kappa_clover_val;

};

#endif //_GAUGEFIELDKHYBRIDH_
