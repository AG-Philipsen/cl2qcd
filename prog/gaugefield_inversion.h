/** @file
 * Provides class Gaugefield_inversion, inherited from Gaugefield, adds solver capabilities
 *
 * @todo global- and local work sizes should be moved inside Opencl_fermions (LZ)
 */
#ifndef _GAUGEFIELDINVERSIONH_
#define _GAUGEFIELDINVERSIONH_

#include "gaugefield.h"
//CP: these dont seem to be needed anymore
//#include "host_operations_spinor.h"
//#include "host_operations_spinorfield.h"
//#include "host_operations_fermionmatrix.h"
#include "opencl_fermions.h"

class Gaugefield_inversion : public Gaugefield {
  public:

  /**
   * Initializes the devices, to be called by init()
   * @param devicetypes array of cl_device_type handles
   */
  virtual void init_devices(cl_device_type* devicetypes);

  /**
   * Free gaugefield and device allocations.
   */
  virtual void finalize();

  /**
   * Perform inversion with point source
   * @return Pointer to OpenCL buffer with solution
   */
  cl_mem* perform_inversion_pointsource(usetimer* solvertimer);
  /**
   * Calculate set of correlators and print to file
   * Use 
   * @param[in] filename where to store the correlators
   */
  void calculate_correlators_from_pointsource();

  /**
   * Perform inversion on device and print pseudo-scalar correlator to std.
   * Use point sources.
   */
  void perform_inversion_pointsource_ps_corr_devices(usetimer* solvertimer);

  /**
   * Free device, called by finalize
   */
  virtual void free_devices();
  
  /**
   * Returns private member opencl_k * devices
   * @return devices of type opencl_k
   * @todo LZ: CHECK IF THIS MAKES SENSE AT ALL!!!!
   */
  Opencl_fermions * get_devices_fermions ();

 private:

};



#endif //_GAUGEFIELDINVERSIONH_
