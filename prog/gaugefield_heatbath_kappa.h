/** @file
 * Provides a class for gauge fields
 *
 */
#ifndef _GAUGEFIELDHEATBATHKAPPAH_
#define _GAUGEFIELDHEATBATHKAPPAH_

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "globaldefs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "inputparameters.h"
#include "host_readgauge.h"
#include "host_writegaugefield.h"
#include "host_use_timer.h"
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "exceptions.h"
#include "gaugefield_hybrid.h"
#include "opencl_module.h"
#include "opencl_module_heatbath.h"

#include "logger.hpp"


/**
 * Version number.
 */
extern string const version;

/**
 * Class for the gaugefield. Includes initialization, device management for multiple devices.
 *
 * @class Gaugefield
 */
class Gaugefield_heatbath_kappa : public Gaugefield_hybrid {
public:
  //init/finalize functions

  virtual void init_devices();


  virtual void delete_variables();
  virtual void finalize_opencl();

private:

};

#endif /* _GAUGEFIELDHEATBATHKAPPAH_ */
