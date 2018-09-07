/** @file
 * Utilities for tests for the hardware classes
 *
 * Copyright (c) 2015,2016 Christopher Pinke
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

#pragma once

#include "system.hpp"

void broadcastMessage_warn(const std::string message);
void broadcastMessage_error(const std::string message);
void broadcastMessage_fatal(const std::string message);
void failTest();
void atLeastOneDeviceMustExistForSanityOfSystem(const hardware::System* system);
bool checkIfNoOpenCLDevicesWereFound(const hardware::OpenclException exception);
bool checkBoostRuntimeArgumentsForGpuUsage();
bool checkBoostRuntimeArgumentsForRec12Usage();
void handleExceptionInTest(hardware::OpenclException& exception);
