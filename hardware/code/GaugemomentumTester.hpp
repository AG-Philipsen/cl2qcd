/*
 * Copyright (c) 2014,2015 Christopher Pinke
 * Copyright (c) 2016 Francesca Cuteri
 * Copyright (c) 2018,2021 Alessandro Sciarra
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

#ifndef GAUGEMOMENTUM_TESTER_HPP_
#define GAUGEMOMENTUM_TESTER_HPP_

#include "../../host_functionality/host_random.hpp"
#include "SpinorTester.hpp"
#include "gaugemomentum.hpp"
#include "kernelTester.hpp"
#include "prng.hpp"

enum class GaugeMomentumFilltype { One, Zero, Ascending };

int calculateGaugemomentumSize(LatticeExtents latticeExtentsIn) noexcept;
int calculateAlgebraSize(LatticeExtents latticeExtentsIn) noexcept;

double count_gm(ae* ae_in, int size);
double calc_var_gm(ae* ae_in, int size, double sum);

struct GaugemomentumTestParameters : public TestParameters {
    GaugemomentumTestParameters(const LatticeExtents latticeExtendsIn)
        : GaugemomentumTestParameters(latticeExtendsIn, GaugeMomentumFilltype::One, 1.0, 1000){};
    GaugemomentumTestParameters(const LatticeExtents latticeExtendsIn, const GaugeMomentumFilltype fillTypesIn,
                                const double alpha = 1.0)
        : GaugemomentumTestParameters(latticeExtendsIn, fillTypesIn, alpha, 1000){};
    GaugemomentumTestParameters(const LatticeExtents latticeExtendsIn, const unsigned int iter)
        : GaugemomentumTestParameters(latticeExtendsIn, GaugeMomentumFilltype::One, 1.0, iter){};

    const GaugeMomentumFilltype fillType;
    const double coefficient;       // Used in saxpy tests
    const unsigned int iterations;  // Used in gaussian tests

  private:
    GaugemomentumTestParameters(const LatticeExtents latticeExtendsIn, const GaugeMomentumFilltype fillTypesIn,
                                const double c, const unsigned int iter)
        : TestParameters(latticeExtendsIn), fillType(fillTypesIn), coefficient(c), iterations(iter){};
};

class GaugemomentumTester : public KernelTester {
  public:
    GaugemomentumTester(const std::string kernelName, const ParameterCollection pC, const ReferenceValues rV,
                        const GaugemomentumTestParameters tP);
    virtual ~GaugemomentumTester();

  protected:
    const hardware::code::Gaugemomentum* code;
    hardware::buffers::Plain<double>* doubleBuffer;
    void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Gaugemomentum* in, int index = 0);
};

struct GaugemomentumCreator {
    GaugemomentumCreator(const LatticeExtents lE) : numberOfElements(calculateAlgebraSize(lE)){};
    ae* createGaugemomentumBasedOnFilltype(const GaugeMomentumFilltype filltype = GaugeMomentumFilltype::One);
    void fill_with_one(ae* in);
    void fill_with_zero(ae* in);
    void fill_with_ascending(ae* in);

    size_t numberOfElements;
};

#endif
