/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

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
