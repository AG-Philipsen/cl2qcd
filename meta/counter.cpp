/*
 * Copyright (c) 2012,2013 Matthias Bach
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

#include "counter.hpp"

meta::Counter::Counter() noexcept : value(0)
{
    // already initialized
}

meta::Counter& meta::Counter::operator+=(const unsigned& inc) noexcept
{
    value += inc;
    return *this;
}

meta::Counter& meta::Counter::operator++() noexcept
{
    value++;
    return *this;
}

meta::Counter::operator unsigned() const noexcept
{
    return value;
}

void meta::Counter::reset() noexcept
{
    value = 0;
}
