#include "halo_update.hpp"
#include "../system.hpp"

hardware::buffers::DeviceAccessibleMemory::DeviceAccessibleMemory(const size_t bytes, hardware::Device * device)
 : Buffer(bytes, device, false)
{
	// nothing to do
}

hardware::buffers::ProxyBufferCache::ProxyBufferCache()
	: cache()
{ }

hardware::buffers::ProxyBufferCache::~ProxyBufferCache()
{
// leave it to the runtime to clean up. we cannot be sure we are executed before cl-finish, and in that case we segfault
// TODO make sure this is run before clfinish
//	for(auto entry: cache) {
//		for(auto buffer: entry.second) {
//			delete buffer;
//		}
//	}
}

hardware::buffers::ProxyBufferCache& hardware::buffers::ProxyBufferCache::getInstance()
{
	static ProxyBufferCache bufferCache;
	return bufferCache;
}

const std::vector<hardware::buffers::DeviceAccessibleMemory*>& hardware::buffers::ProxyBufferCache::getBuffers(size_t rows, size_t bytes, const std::vector<hardware::Device*>& devices)
{
	auto primary = devices[0];
	auto& buffers = this->cache[std::make_pair(primary->context,std::make_pair(rows, bytes))];
	if(buffers.size() == 0) {
		const size_t num_devs = devices.size();
		buffers.resize(rows * num_devs);
		for(size_t row = 0; row < rows; ++row) {
			for(size_t i = 0; i < num_devs; ++i) {
				buffers[row * num_devs + i] = new DeviceAccessibleMemory(bytes, devices[i]);
			}
		}
	}
	return buffers;
}
