/** @file
 *
 * Copyright (c) 2018 Francesca Cuteri
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

#include "parametersIntegrator.hpp"

#include "../executables/exceptions.hpp"
#include "../host_functionality/logger.hpp"

#include <boost/algorithm/string.hpp>
#include <stdexcept>

using namespace meta;

static common::integrator translateIntegratorToEnum(std::string);

double ParametersIntegrator::get_tau() const noexcept
{
    return tau;
}
bool ParametersIntegrator::get_reversibility_check() const noexcept
{
    return reversibility_check;
}
int ParametersIntegrator::get_integrationsteps(size_t timescale) const noexcept
{
    switch (timescale) {
        case 0:
            return integrationsteps0;
        case 1:
            return integrationsteps1;
        case 2:
            return integrationsteps2;
        default:
            throw std::out_of_range("No such timescale");
    }
}

int ParametersIntegrator::get_num_timescales() const noexcept
{
    return num_timescales;
}
common::integrator ParametersIntegrator::get_integrator(size_t timescale) const noexcept
{
    switch (timescale) {
        case 0:
            return integrator0;
        case 1:
            return integrator1;
        case 2:
            return integrator2;
        default:
            throw std::out_of_range("No such timescale");
    }
}
double ParametersIntegrator::get_lambda(size_t timescale) const noexcept
{
    switch (timescale) {
        case 0:
            return lambda0;
        case 1:
            return lambda1;
        case 2:
            return lambda2;
        default:
            throw std::out_of_range("No such timescale");
    }
}

ParametersIntegrator::ParametersIntegrator()
    : tau(0.5)
    , reversibility_check(false)
    , integrationsteps0(10)
    , integrationsteps1(10)
    , integrationsteps2(10)
    , num_timescales(1)
    , lambda0(0.1931833275037836)
    , lambda1(0.1931833275037836)
    , lambda2(0.1931833275037836)
    , options("Integrator options")
    , integrator0String("leapfrog")
    , integrator1String("leapfrog")
    , integrator2String("leapfrog")
    , integrator0(common::integrator::leapfrog)
    , integrator1(common::integrator::leapfrog)
    , integrator2(common::integrator::leapfrog)
{
    // clang-format off
    options.add_options()
    ("tau", po::value<double>(&tau)->default_value(tau),"The total time of integration.")
    ("reversibilityCheck", po::value<bool>(&reversibility_check)->default_value(reversibility_check),"Whether to check reversibility violation in the integrator by integrating back in time.")
    ("nTimeScales", po::value<int>(&num_timescales)->default_value(num_timescales),"The number of time scales (timescale0 for the gauge-part, timescale1 for the fermion, timescale2 for mass preconditioning). Consider that different timescales must use the same integrator.")
    ("integrator0", po::value<std::string>(&integrator0String)->default_value(integrator0String),"The integration scheme for timescale 0 (one among leapfrog and twomn).")
    ("integrator1", po::value<std::string>(&integrator1String)->default_value(integrator1String),"The integration scheme for timescale 1 (one among leapfrog and twomn).")
    ("integrator2", po::value<std::string>(&integrator2String)->default_value(integrator2String),"The integration scheme for timescale 2 (one among leapfrog and twomn).")
    ("integrationSteps0", po::value<int>(&integrationsteps0)->default_value(integrationsteps0),"The number of integration steps for timescale 0.")
    ("integrationSteps1", po::value<int>(&integrationsteps1)->default_value(integrationsteps1),"The number of integration steps for timescale 1.")
    ("integrationSteps2", po::value<int>(&integrationsteps2)->default_value(integrationsteps2),"The number of integration steps for timescale 2.")
    // this is the optimal value...
    ("lambda0", po::value<double>(&lambda0)->default_value(lambda0),"The lambda parameter for timescale 0 (for the twomn integrator).")
    ("lambda1", po::value<double>(&lambda1)->default_value(lambda1),"The lambda parameter for timescale 1 (for the twomn integrator).")
    ("lambda2", po::value<double>(&lambda2)->default_value(lambda2),"The lambda parameter for timescale 2 (for the twomn integrator).");
    // clang-format on
}

static common::integrator translateIntegratorToEnum(std::string s)
{
    boost::algorithm::to_lower(s);
    std::map<std::string, common::integrator> m;
    m["leapfrog"] = common::leapfrog;
    m["twomn"]    = common::twomn;
    m["2mn"]      = common::twomn;

    common::integrator a = m[s];
    if (a) {
        return a;
    } else {
        throw Invalid_Parameters("Unkown integrator!", "leapfrog, twomn or 2mn", s);
    }
}

void meta::ParametersIntegrator::makeNeededTranslations()
{
    integrator0 = translateIntegratorToEnum(integrator0String);
    integrator1 = translateIntegratorToEnum(integrator1String);
    integrator2 = translateIntegratorToEnum(integrator2String);
}
