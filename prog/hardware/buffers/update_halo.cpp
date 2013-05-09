#include "halo_update.hpp"
#include "../system.hpp"

hardware::buffers::DeviceAccessibleMemory::DeviceAccessibleMemory(const size_t bytes, const hardware::Device * device)
 : queue(device->command_queue)
{
	cl_int clerr;
	this->buf = clCreateBuffer(device->context, CL_MEM_ALLOC_HOST_PTR, bytes, nullptr, &clerr);
	if(clerr) {
		throw hardware::OpenclException(clerr, "clCreateBuffer", __FILE__, __LINE__);
	}
	this->mem = reinterpret_cast<char*>(clEnqueueMapBuffer(queue, this->buf, CL_TRUE, CL_MAP_WRITE_INVALIDATE_REGION, 0, bytes, 0, 0, 0, &clerr));
	if(clerr) {
		throw hardware::OpenclException(clerr, "clEnqueueMapBuffer", __FILE__, __LINE__);
	}
}

hardware::buffers::DeviceAccessibleMemory::~DeviceAccessibleMemory()
{
	cl_int clerr;
	clerr = clEnqueueUnmapMemObject(queue, this->buf, this->mem, 0, 0, 0);
	if(clerr) {
		throw hardware::OpenclException(clerr, "clEnqueueUnmapBuffer", __FILE__, __LINE__);
	}
	clerr = clReleaseMemObject(this->buf);
	if(clerr) {
		throw hardware::OpenclException(clerr, "clReleaseMemObject", __FILE__, __LINE__);
	}
}
