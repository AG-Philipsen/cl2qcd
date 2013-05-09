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

hardware::buffers::HostBufferCache::HostBufferCache()
	: cache()
{ }

hardware::buffers::HostBufferCache::~HostBufferCache()
{
	for(auto entry: cache) {
		for(auto buffer: entry.second) {
			delete buffer;
		}
	}
}

hardware::buffers::HostBufferCache& hardware::buffers::HostBufferCache::getInstance()
{
	static HostBufferCache bufferCache;
	return bufferCache;
}

const std::vector<hardware::buffers::DeviceAccessibleMemory*>& hardware::buffers::HostBufferCache::getBuffers(size_t num, size_t bytes, const hardware::Device* primary)
{
	auto& buffers = this->cache[std::make_pair(num, bytes)];
	if(buffers.size() == 0) {
		buffers.resize(num);
		for(size_t i = 0; i < num; ++i) {
			buffers[i] = new DeviceAccessibleMemory(bytes, primary);
		}
	}
	return buffers;
}
