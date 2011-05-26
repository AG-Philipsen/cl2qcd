/** @file
 * Provides class Gaugefield_inversion, inherited from Gaugefield, adds solver capabilities
 *
 */
#ifndef _GAUGEFIELDINVERSIONH_
#define _GAUGEFIELDINVERSIONH_

#include "gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"
#include "opencl_fermions.h"

class Gaugefield_inversion : public Gaugefield {
  public:

  /**
   * Initializes the devices, to be called by init()
   * @return Error code as defined in hmcerrs.h
   * @param devicetypes array of cl_device_type handles
   * @param[in,out] timer timer for initialization
   */
  virtual hmc_error init_devices(cl_device_type* devicetypes, usetimer* timer);

  /**
   * Free gaugefield and device allocations.
   */
  virtual hmc_error finalize();

  /**
   * Perform inversion on host.
   */
  hmc_error perform_inversion_on_host();

  /**
   * Free device, called by finalize
   */
  virtual hmc_error free_devices();
  
  /**
   * Returns private member opencl_k * devices
   * @return devices of type opencl_k
   * @todo LZ: CHECK IF THIS MAKES SENSE AT ALL!!!!
   */
  Opencl_fermions * get_devices_fermions ();

 private:

};



#endif //_GAUGEFIELDINVERSIONH_
