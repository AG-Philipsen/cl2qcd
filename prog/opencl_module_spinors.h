/** @file
 * Heatbath for OpenCL
 */
#ifndef _OPENCLMODULSPINORSH_
#define _OPENCLMODULSPINORSH_

#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "host_geometry.h"
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_random.h"
#include "inputparameters.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"
#include "opencl_module_ran.h"

#include "exceptions.h"

#include "types_fermions.h"

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Spinors : public Opencl_Module_Ran {
public:
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(stringstream* collect_options);

	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual void fill_buffers();

	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();

	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual void clear_kernels();

	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual void clear_buffers();

	/**
	 * Add specific work_size determination for this child class
	 */
	virtual void get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups);


	/////////////////////////////////////////
	// device operations

	//    linear Algebra operations
	void convert_from_eoprec_device(cl_mem in1, cl_mem in2, cl_mem out);
	void convert_to_eoprec_device(cl_mem out1, cl_mem out2, cl_mem in);

	void set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out);
	void set_complex_to_scalar_product_eoprec_device(cl_mem a, cl_mem b, cl_mem out);
	void set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out);
	void set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out);
	void set_float_to_global_squarenorm_device(cl_mem a, cl_mem out);
	void set_float_to_global_squarenorm_eoprec_device(cl_mem a, cl_mem out);
	void set_zero_spinorfield_device(cl_mem x);
	void set_zero_spinorfield_eoprec_device(cl_mem x);
	void saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out);
	void sax_device(cl_mem x, cl_mem alpha, cl_mem out);
	void saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out);
	void saxpy_eoprec_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out);
	void sax_eoprec_device(cl_mem x, cl_mem alpha, cl_mem out);
	void saxsbypz_eoprec_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out);
	void create_point_source_device(cl_mem inout, int i, int spacepos, int timepos);
	void create_point_source_eoprec_device(cl_mem inout_even, cl_mem inout_odd, cl_mem gf, int i, int spacepos, int timepos);
	void set_spinorfield_cold_device(cl_mem inout);
	void set_eoprec_spinorfield_cold_device(cl_mem inout);

#ifdef _PROFILING_

	//BLAS
	usetimer timer_set_spinorfield_cold;
	usetimer timer_set_eoprec_spinorfield_cold;
	usetimer timer_convert_from_eoprec;
	usetimer timer_convert_to_eoprec;
	usetimer timer_saxpy;
	usetimer timer_sax;
	usetimer timer_saxsbypz;
	usetimer timer_set_zero_spinorfield;
	usetimer timer_create_point_source;
	usetimer timer_saxpy_eoprec;
	usetimer timer_sax_eoprec;
	usetimer timer_saxsbypz_eoprec;
	usetimer timer_set_zero_spinorfield_eoprec;
	usetimer timer_create_point_source_eoprec;

	//Scalar Product
	usetimer timer_scalar_product;
	usetimer timer_scalar_product_reduction;
	usetimer timer_global_squarenorm;
	usetimer timer_global_squarenorm_reduction;
	usetimer timer_scalar_product_eoprec;
	usetimer timer_global_squarenorm_eoprec;

	//Single
	usetimer timer_ratio;
	usetimer timer_product;

	//Observables
	usetimer timer_ps_correlator;

	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(const char * in);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(std::string filename, int number);

#endif

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual int get_read_write_size(const char * in);

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual int get_flop_size(const char * in);


protected:

	cl_mem clmem_scalar_product_buf_glob;
	cl_mem clmem_global_squarenorm_buf_glob;


	//BLAS
	cl_kernel set_spinorfield_cold;
	cl_kernel saxpy;
	cl_kernel sax;
	cl_kernel saxsbypz;
	cl_kernel set_zero_spinorfield;

	cl_kernel convert_from_eoprec;
	cl_kernel convert_to_eoprec;
	cl_kernel set_eoprec_spinorfield_cold;
	cl_kernel set_zero_spinorfield_eoprec;
	cl_kernel saxpy_eoprec;
	cl_kernel sax_eoprec;
	cl_kernel saxsbypz_eoprec;

	//Scalar Product
	cl_kernel scalar_product;
	cl_kernel scalar_product_reduction;
	cl_kernel global_squarenorm;
	cl_kernel global_squarenorm_reduction;
	cl_kernel scalar_product_eoprec;
	cl_kernel global_squarenorm_eoprec;

	//Single
	cl_kernel ratio;
	cl_kernel product;

	ClSourcePackage basic_fermion_code;

private:

};

#endif //OPENCLMODULSPINORSH
