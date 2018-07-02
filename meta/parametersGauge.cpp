/** @file
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2015,2018 Francesca Cuteri
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

#include "parametersGauge.hpp"

#include "../executables/exceptions.hpp"

#include <boost/algorithm/string.hpp>

// gaugefield parameters
double meta::ParametersGauge::get_beta() const noexcept
{
    return beta;
}
double meta::ParametersGauge::get_rho() const noexcept
{
    return rho;
}
int meta::ParametersGauge::get_rho_iter() const noexcept
{
    return rho_iter;
}
common::action meta::ParametersGauge::get_gaugeact() const noexcept
{
    return translateGaugeActionToEnum();
}

bool meta::ParametersGauge::get_use_smearing() const noexcept
{
    return use_smearing;
}

meta::ParametersGauge::ParametersGauge() : options("Gaugefield options")
{
    // clang-format off
    options.add_options()
    ("beta", po::value<double>(&beta)->default_value(4.0),"The beta-coupling in the gauge action.")
    ("useSmearing", po::value<bool>(&use_smearing)->default_value(false),"Whether to apply stout smearing to the gaugefield.")
    ("smearingFactor", po::value<double>(&rho)->default_value(0.),"The weight factor associated with the staples in stout smearing.")
    ("nSmearingSteps", po::value<int>(&rho_iter)->default_value(0),"The number of stout-smearing steps.")
    ("gaugeAction", po::value<std::string>(&gaugeact)->default_value("wilson"),"Which type of gauge action to use (e.g. wilson, tlsym, iwasaki, dbw2).");
    // clang-format on
}

common::action meta::ParametersGauge::translateGaugeActionToEnum() const
{
    std::string s = gaugeact;
    boost::algorithm::to_lower(s);
    std::map<std::string, common::action> m;
    m["wilson"]  = common::action::wilson;
    m["tlsym"]   = common::action::tlsym;
    m["iwasaki"] = common::action::iwasaki;
    m["dbw2"]    = common::action::dbw2;

    common::action a = m[s];
    if (a) {  // map returns 0 if element is not found
        return a;
    } else {
        throw Invalid_Parameters("Unkown gauge action!", "wilson, tlsym, iwasaki, dbw2", s);
    }
}
