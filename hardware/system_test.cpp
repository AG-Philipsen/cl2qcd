/** @file
 * Testcases for the hardware::System class
 *
 * Copyright (c) 2012,2013 Matthias Bach
 * Copyright (c) 2015,2016 Christopher Pinke
 * Copyright (c) 2015 Francesca Cuteri
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::System
#include "system.hpp"

#include "device.hpp"
#include "hardwareTestUtilities.hpp"
#include "interfaceMockups.hpp"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(initialization)

    BOOST_AUTO_TEST_CASE(fromHardwareParameters)
    {
        const hardware::HardwareParametersMockup hardwareParameters(4, 4);
        const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
        BOOST_CHECK_NO_THROW(hardware::System system(hardwareParameters, kernelParameters));
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(systemSanity)

    void allDevicesMustSupportDoublePrecisionForSanityOfSystem(const hardware::System* system)
    {
        for (hardware::Device* device : system->get_devices()) {
            BOOST_REQUIRE_EQUAL(device->is_double_supported(), true);
        }
    }

    BOOST_AUTO_TEST_CASE(staticCastToPlatform)
    {
        const hardware::HardwareParametersMockup hardwareParameters(4, 4);
        const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
        hardware::System system(hardwareParameters, kernelParameters);
        BOOST_REQUIRE(static_cast<const cl_context>(system.getContext()));
    }

    BOOST_AUTO_TEST_CASE(enoughDevicesExist)
    {
        const hardware::HardwareParametersMockup hardwareParameters(4, 4);
        const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
        hardware::System system(hardwareParameters, kernelParameters);
        atLeastOneDeviceMustExistForSanityOfSystem(&system);
    }

    BOOST_AUTO_TEST_CASE(allDevicesHaveDoubleSupport)
    {
        const hardware::HardwareParametersMockup hardwareParameters(4, 4);
        const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
        hardware::System system(hardwareParameters, kernelParameters);
        allDevicesMustSupportDoublePrecisionForSanityOfSystem(&system);
    }

BOOST_AUTO_TEST_SUITE_END()

std::string getDeviceTypeAsName(const cl_device_type device_type)
{
    // https://www.khronos.org/registry/OpenCL/sdk/1.0/docs/man/xhtml/enums.html
    switch (device_type) {
        case CL_DEVICE_TYPE_CPU:
            return "CPU";
            break;
        case CL_DEVICE_TYPE_GPU:
            return "GPU";
            break;
        case CL_DEVICE_TYPE_ACCELERATOR:
            return "ACCELERATOR";
            break;
        case CL_DEVICE_TYPE_DEFAULT:
        case CL_DEVICE_TYPE_ALL:
            throw Print_Error_Message("Unexpected cl_device_type to be translated!");
            break;
        default:
            throw Print_Error_Message("Suspicious cl_device_type to be translated, please investigate!");
            break;
    }
}

void enableSpecificDeviceTypeOnly(const cl_device_type device_type, const hardware::HardwareParametersInterface& hI,
                                  const bool failTestIfMissingDeviceType = false)
{
    const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
    try {
        hardware::System system(hI, kernelParameters);
        for (hardware::Device* device : system.get_devices()) {
            BOOST_REQUIRE_EQUAL(device->get_device_type(), device_type);
        }
    } catch (hardware::OpenclException& exception) {
        if (checkIfNoOpenCLDevicesWereFound(exception)) {
            broadcastMessage_error("System does not seem to contain \"" + getDeviceTypeAsName(device_type) +
                                   "\" devices!");
            if (failTestIfMissingDeviceType)
                failTest();
        } else {
            broadcastMessage_fatal("Got unknown error code. Aborting...");
            failTest();
        }
    }
}

BOOST_AUTO_TEST_SUITE(devices)

    void checkThatOnlySpecifiedNumberOfDevicesIsInitialized(const int totalNumberOfDevicesInSystem)
    {
        for (int desiredNumberOfDevices = 1; desiredNumberOfDevices <= totalNumberOfDevicesInSystem;
             desiredNumberOfDevices++) {
            const hardware::HardwareParametersMockupForDeviceSelection hardwareParameters(4, 4, desiredNumberOfDevices,
                                                                                          std::vector<int>{0});
            const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
            hardware::System system(hardwareParameters, kernelParameters);
            BOOST_REQUIRE_EQUAL(system.get_devices().size(), desiredNumberOfDevices);
        }
    }

    BOOST_AUTO_TEST_CASE(setNumberOfDevicesByCommandLine)
    {
        const hardware::HardwareParametersMockup hardwareParameters(4, 4);
        const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
        hardware::System system(hardwareParameters, kernelParameters);
        checkThatOnlySpecifiedNumberOfDevicesIsInitialized(system.get_devices().size());
    }

    void checkThatOnlySpecifiedDeviceIsInitialized(const int totalNumberOfDevicesInSystem)
    {
        for (int desiredDevice = 0; desiredDevice < totalNumberOfDevicesInSystem; desiredDevice++) {
            const hardware::HardwareParametersMockupForDeviceSelection hardwareParameters(4, 4, 1,
                                                                                          std::vector<int>{
                                                                                              desiredDevice});
            const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);

            try {
                hardware::System system(hardwareParameters, kernelParameters);
                BOOST_REQUIRE_EQUAL(system.get_devices().size(), 1);
            } catch (std::invalid_argument) {
                // device might not support double-precision
            }
        }
    }

    BOOST_AUTO_TEST_CASE(setDeviceByCommandLine)
    {
        const hardware::HardwareParametersMockup hardwareParameters(4, 4);
        const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
        hardware::System system(hardwareParameters, kernelParameters);
        checkThatOnlySpecifiedDeviceIsInitialized(system.get_devices().size());
    }

    BOOST_AUTO_TEST_CASE(enableCpusOnly)
    {
        const hardware::HardwareParametersMockupWithCpusOnly hardwareParameters(4, 4);
        enableSpecificDeviceTypeOnly(CL_DEVICE_TYPE_CPU, hardwareParameters);
    }

    BOOST_AUTO_TEST_CASE(enableGpusOnly)
    {
        const hardware::HardwareParametersMockupWithGpusOnly hardwareParameters(4, 4);
        enableSpecificDeviceTypeOnly(CL_DEVICE_TYPE_GPU, hardwareParameters);
    }

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(checkDeviceAvailability)

    /*
     * NOTE: The following two tests are to create explicit tests in ctest in order to make the user easily discover
     *       which type of device is available on his/her architecture.
     */
    BOOST_AUTO_TEST_CASE(areCpusAvailable)
    {
        const hardware::HardwareParametersMockupWithCpusOnly hardwareParameters(4, 4);
        enableSpecificDeviceTypeOnly(CL_DEVICE_TYPE_CPU, hardwareParameters, true);
    }

    BOOST_AUTO_TEST_CASE(areGpusAvailable)
    {
        const hardware::HardwareParametersMockupWithGpusOnly hardwareParameters(4, 4);
        enableSpecificDeviceTypeOnly(CL_DEVICE_TYPE_GPU, hardwareParameters, true);
    }

BOOST_AUTO_TEST_SUITE_END()

void checkOnProperEnvironmentSettings()
{
    BOOST_REQUIRE_EQUAL(std::string("3"), std::string(getenv("GPU_DUMP_DEVICE_KERNEL")));
    BOOST_REQUIRE_NE(std::string(getenv("AMD_OCL_BUILD_OPTIONS_APPEND")).find("-save-temps"), -1);
}

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(dump_source_if_debugging, 1)
BOOST_AUTO_TEST_CASE(dump_source_if_debugging)
{
    /**
     * @todo: if debug mode is not activated at compile time, this switch does not seem to have any effect. Investigate!
     */
    switchLogLevel("debug");
    const hardware::HardwareParametersMockup hardwareParameters(4, 4);
    const hardware::code::OpenClKernelParametersMockup kernelParameters(4, 4);
    hardware::System system(hardwareParameters, kernelParameters);

    if (logger.beDebug()) {
        checkOnProperEnvironmentSettings();
    } else {
        broadcastMessage_fatal("Something went wrong, logger not in debug mode...");
        failTest();
    }
}
