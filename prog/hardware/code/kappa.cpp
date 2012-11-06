#include "kappa.hpp"

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

static std::string collect_build_options(hardware::Device *, const meta::Inputparameters& params)
{
	std::ostringstream options;
	options <<  "-D BETA=" << params.get_beta();
	options <<  " -D XI_0=" << meta::get_xi_0(params);

	return options.str();
}


void hardware::code::Kappa::fill_kernels()
{
	ClSourcePackage sources = get_device()->get_gaugefield_code()->get_sources() << ClSourcePackage(collect_build_options(get_device(), get_parameters()));

	cout << "Create TK clover kernels..." << endl;
	kappa_clover_gpu = createKernel("kappa_clover_gpu") << sources << "opencl_tk_kappa.cl";
}


void hardware::code::Kappa::run_kappa_clover(const hardware::buffers::SU3 * gaugefield, const hmc_float beta)
{
	//variables
	cl_int clerr = CL_SUCCESS;

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;

	get_work_sizes(kappa_clover_gpu, &local_work_size, &global_work_size, &num_groups);

	hardware::buffers::Plain<hmc_float> clmem_kappa_clover_buf_glob(num_groups, get_device());

	clerr = clSetKernelArg(kappa_clover_gpu, 0, sizeof(cl_mem), gaugefield->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 2, sizeof(cl_mem), clmem_kappa_clover_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(kappa_clover_gpu, global_work_size, local_work_size);

	// wait for results to have been read back
	//don't do that anymore ;-)
	//  clFinish(queue);
	//  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr,"clFinish",__FILE__,__LINE__);
}

hmc_float hardware::code::Kappa::get_kappa_clover()
{
	hmc_float kappa_clover;
	clmem_kappa_clover.dump(&kappa_clover);
	return kappa_clover;
}

void hardware::code::Kappa::clear_kernels()
{
	cl_int clerr = clReleaseKernel(kappa_clover_gpu);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}

hardware::code::Kappa::Kappa(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device), clmem_kappa_clover(1, device)
{
	fill_kernels();
}

hardware::code::Kappa::~Kappa()
{
	clear_kernels();
}
