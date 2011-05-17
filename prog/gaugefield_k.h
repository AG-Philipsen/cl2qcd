#include "gaugefield.h"

class Gaugefield_k : public Gaugefield{

  public:
  	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error kappa_karsch ();
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor built by a Clover discretization
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error kappa_clover ();
	
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
	 * @param[in] local_work_size OpenCL local_work_size
	 * @param[in] global_work_size OpenCL global_work_size
	 * @param[in,out] timer time for measurement
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error kappa_karsch_gpu (const size_t local_work_size, const size_t global_work_size, usetimer* timer_karsch);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld  on GPU
	 * @param[in] local_work_size OpenCL local_work_size
	 * @param[in] global_work_size OpenCL global_work_size
	 * @param[in,out] timer time for measurement
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error kappa_clover_gpu (const size_t local_work_size, const size_t global_work_size, usetimer* timer_karsch);
  
  
  private:
  hmc_float kappa_karsch_val;
  hmc_float kappa_clover_val;
  

};