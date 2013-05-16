#include "halo_update.hpp"
#include "../system.hpp"

hardware::buffers::DeviceAccessibleMemory::DeviceAccessibleMemory(const size_t bytes, hardware::Device * device)
 : Buffer(bytes, device, true)
{
	// nothing to do
}

hardware::buffers::HostBufferCache::HostBufferCache()
	: cache()
{ }

hardware::buffers::HostBufferCache::~HostBufferCache()
{
// leave it to the runtime to clean up. we cannot be sure we are executed before cl-finish, and in that case we segfault
// TODO make sure this is run before clfinish
//	for(auto entry: cache) {
//		for(auto buffer: entry.second) {
//			delete buffer;
//		}
//	}
}

hardware::buffers::HostBufferCache& hardware::buffers::HostBufferCache::getInstance()
{
	static HostBufferCache bufferCache;
	return bufferCache;
}

const std::vector<hardware::buffers::DeviceAccessibleMemory*>& hardware::buffers::HostBufferCache::getBuffers(size_t num, size_t bytes, hardware::Device* primary)
{
	auto& buffers = this->cache[std::make_pair(primary->context,std::make_pair(num, bytes))];
	if(buffers.size() == 0) {
		buffers.resize(num);
		for(size_t i = 0; i < num; ++i) {
			buffers[i] = new DeviceAccessibleMemory(bytes, primary);
		}
	}
	return buffers;
}
