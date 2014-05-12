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

#ifndef KERNELTESTER_H_
#define KERNELTESTER_H_

#include "../physics/prng.hpp"
#include "../hardware/system.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../host_functionality/logger.hpp"
#include "testUtilities.hpp"

class KernelTester
{
public:
  KernelTester(std::string kernelNameIn, std::string inputfileIn); 
  virtual void callSpecificKernel() {};

protected:
  const hardware::buffers::SU3* getGaugefieldBuffer()
  {
    return gaugefield->get_buffers()[0];
  }

  double testPrecision;
  meta::Inputparameters * parameters;
  hardware::System * system;
  hardware::Device * device;

  physics::PRNG * prng;  
  physics::lattices::Gaugefield * gaugefield;
};

class KernelTesterDouble: public KernelTester
{
public:
  KernelTesterDouble(std::string kernelNameIn, std::string inputfileIn);
  ~KernelTesterDouble();
protected:
  double referenceValue;
  double kernelResult;
};

class KernelTesterComplex: public KernelTester
{
public:
  KernelTesterComplex(std::string kernelNameIn, std::string inputfileIn);
  ~KernelTesterComplex();
protected:
  hmc_complex referenceValue;
  hmc_complex kernelResult;
};

#endif /* KERNELTESTER_H_ */
