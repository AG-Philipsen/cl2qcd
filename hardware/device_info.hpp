/** @file
 * Declaration of the hardware::DeviceInfo class
 *
 * Copyright (c) 2013 Matthias Bach
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _HARDWARE_DEVICE_INFO_HPP_
#define _HARDWARE_DEVICE_INFO_HPP_

#ifdef __APPLE__
#    include <OpenCL/cl.h>
#else
#    include <CL/cl.h>
#endif
#include <string>

namespace hardware {

    /**
     * Wrapper to represent some basic information off an OpenCL device
     */
    class DeviceInfo {
      public:
        /**
         * Initialize the info object.
         *
         * \param device_id id of the device to initialize
         */
        DeviceInfo(const cl_device_id device_id);

        /**
         * Initialize the info object.
         *
         * \param device_id id of the device to initialize
         */
        DeviceInfo(const DeviceInfo& other);

        /**
         * Get the number of compute units of the device
         */
        size_t get_num_compute_units() const noexcept;

        /**
         * Get the type of the OpenCL device
         */
        cl_device_type get_device_type() const noexcept;

        /**
         * Checks whether the device supports double precision.
         */
        bool is_double_supported() const noexcept;

        /**
         * Get the name of this device
         */
        std::string get_name() const noexcept;

        /**
         * Get the id of this device
         */
        cl_device_id get_id() const noexcept;

        /**
         * Get the prefered local work size of this device
         */
        size_t get_preferred_local_thread_num() const noexcept;

        /**
         * Get the default number of threads on this device
         */
        size_t get_preferred_global_thread_num() const noexcept;

        /**
         * Whether this device prefers blocked loops
         */
        bool get_prefers_blocked_loops() const noexcept;

        /**
         * Whether this device prefers SOA storage
         *
         * @todo This should be datatype specific
         */
        bool get_prefers_soa() const noexcept;

        /**
         * Check if the given extension is supported by the device.
         */
        bool check_extension(std::string) const;

      private:
        /**
         * The OpenCL device id of the device
         */
        cl_device_id device_id;

        /**
         * The prefered local work size of this device
         */
        size_t preferred_local_thread_num;

        /**
         * The preferred global work size of this device
         */
        size_t preferred_global_thread_num;

        /**
         * The number of compute units of the device
         */
        size_t num_compute_units;

        /**
         * The type of the device
         */
        cl_device_type device_type;

        /**
         * Whether this device prefers blocked loops
         */
        bool prefers_blocked_loops;

        /**
         * Whether this device prefers SOA storage
         *
         * @todo This should be datatype specific
         */
        bool prefers_soa;

        /**
         * The name of this device
         */
        std::string name;
    };

} /* namespace hardware */

#endif /* _HARDWARE_DEVICE_INFO_HPP_ */
