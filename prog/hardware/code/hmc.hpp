/** @file
 * Basic OpenCL functionality
 */
#ifndef _HARDWARE_CODE_HMC_
#define _HARDWARE_CODE_HMC_

#include "../../host_use_timer.h"
#include "opencl_module.hpp"
#include "../../types.h"
#include "../../types_fermions.h"
#include "../../types_hmc.h"

#include "../../meta/counter.hpp"
#include "../hardware/buffers/plain.hpp"
#include "../hardware/buffers/prng_buffer.hpp"
#include "../hardware/buffers/su3.hpp"
#include "../hardware/buffers/spinor.hpp"
#include "../hardware/buffers/gaugemomentum.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Hmc : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Hmc();

	////////////////////////////////////////////////////
	//get members
	const hardware::buffers::Gaugemomentum * get_clmem_p();
	const hardware::buffers::Gaugemomentum * get_clmem_new_p();
	const hardware::buffers::SU3 * get_new_u();
	const hardware::buffers::Plain<spinor> * get_clmem_phi();
	const hardware::buffers::Plain<spinor> * get_clmem_phi_inv();
	const hardware::buffers::Spinor * get_clmem_phi_eo();
	const hardware::buffers::Spinor * get_clmem_phi_inv_eo();
	const hardware::buffers::Plain<spinor> * get_clmem_phi_mp();
	const hardware::buffers::Spinor * get_clmem_phi_mp_eo();
	hardware::buffers::Plain<hmc_float> * get_clmem_s_fermion_init();
	hardware::buffers::Plain<hmc_float> * get_clmem_s_fermion_mp_init();

	////////////////////////////////////////////////////
	//Methods needed for the HMC-algorithm
	void md_update_spinorfield(const hardware::buffers::SU3 * gaugefield, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void md_update_spinorfield_mp(usetimer * solvertimer, const hardware::buffers::SU3 * gaugefield);
	void generate_spinorfield_gaussian(const hardware::buffers::PRNGBuffer * prng);
	hmc_observables metropolis(hmc_float rnd, hmc_float beta, const hardware::buffers::SU3 * gaugefield);
	void calc_spinorfield_init_energy(hardware::buffers::Plain<hmc_float> * dest);
	void calc_gauge_force();
	void calc_fermion_force(usetimer * solvertimer, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void calc_fermion_force_detratio(usetimer * solvertimer, const hardware::buffers::SU3 * gaugefield);
	void md_update_gaugemomentum_device(hmc_float eps);
	void md_update_gaugefield_device(hmc_float eps);
	void set_zero_clmem_force_device();
	void gauge_force_device();
	void gauge_force_tlsym_device();
	void fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, hmc_float kappa = ARG_DEF);
	void fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, int evenodd, hmc_float kappa = ARG_DEF);
	hmc_float calc_s_fermion();
	hmc_float calc_s_fermion_mp(const hardware::buffers::SU3 * gaugefield);

	void virtual print_profiling(const std::string& filename, int number) const override;

protected:

	/**
	 * These functions can be removed once this module is turned into a different object
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;
	virtual size_t get_read_write_size(const std::string& in) const override;
	virtual uint64_t get_flop_size(const std::string& in) const override;

private:
	/**
	 * Constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Hmc(const meta::Inputparameters& params, hardware::Device * device);

	//variables
	//initial energy of the (gaussian) spinorfield
	hardware::buffers::Plain<hmc_float> clmem_s_fermion_init;
	hardware::buffers::Plain<hmc_float> clmem_s_fermion_mp_init;
	//squarenorm temps
	hardware::buffers::Plain<hmc_float> clmem_p2;
	hardware::buffers::Plain<hmc_float> clmem_new_p2;
	hardware::buffers::Plain<hmc_float> clmem_s_fermion;
	//new and old gaugemomentum, new gaugefield
	const hardware::buffers::Gaugemomentum clmem_p;
	const hardware::buffers::Gaugemomentum clmem_new_p;
	const hardware::buffers::SU3 new_u;
	//force field
	const hardware::buffers::Gaugemomentum clmem_force;
	//inverted spinorfield
	const hardware::buffers::Plain<spinor> clmem_phi_inv;
	const hardware::buffers::Spinor clmem_phi_inv_eo;
	//D(gaussian spinorfield)
	const hardware::buffers::Plain<spinor> clmem_phi;
	const hardware::buffers::Plain<spinor> clmem_phi_mp;
	const hardware::buffers::Spinor clmem_phi_eo;
	const hardware::buffers::Spinor clmem_phi_mp_eo;
};

}

}

#endif // _HARDWARE_CODE_HMC_
