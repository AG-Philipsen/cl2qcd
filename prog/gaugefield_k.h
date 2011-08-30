/** @file
 * Provides a class for gauge fields, including calculation of transport coefficients
 *
 */
#ifndef _GAUGEFIELDKH_
#define _GAUGEFIELDKH_

#include "gaugefield.h"
#include "gaugefield_heatbath.h"
#include "opencl_heatbath.h"
#include "opencl_k.h"

class Gaugefield_k : public Gaugefield_heatbath  {

  public:
    
  /**
   * Initializes the devices, to be called by init()
   * @param devicetypes array of cl_device_type handles
   * @param[in,out] timer timer for initialization
   */
  virtual void init_devices(cl_device_type* devicetypes);

	/**
	 * Free device, called by finalize
	 */
	virtual void free_devices();

  /**
   * Free gaugefield and device allocations.
   */
  virtual void finalize();

	
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
	 */	
	void set_kappa_karsch (hmc_float in);
	
	/**
	 * Set the transport coefficient kappa computed by Clover method
	 * @param[in] in Result for the transport coefficient kappa
	 */	
	void set_kappa_clover (hmc_float in);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld on GPU
	 * @param[in,out] timer time for measurement
	 */	
	 void kappa_karsch_gpu (usetimer* timer);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld  on GPU
	 */	
	void kappa_clover_gpu (usetimer* timer);
  
	/**
	 * Returns private member opencl_k * devices
	 * @return devices of type opencl_k
	 * @todo LZ: CHECK IF THIS MAKES SENSE AT ALL!!!!
	 */
 	Opencl_k * get_devices_k ();

  private:
  hmc_float kappa_karsch_val;
  hmc_float kappa_clover_val;
  

};

#endif //_GAUGEFIELDKH_
