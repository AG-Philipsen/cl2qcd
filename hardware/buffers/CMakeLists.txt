# Build specifications for the contents of the buffers package
#
# Copyright (c) 2012,2013 Matthias Bach
# Copyright (c) 2013,2018 Alessandro Sciarra
# Copyright (c) 2014 Christopher Pinke
#
# This file is part of CL2QCD.
#
# CL2QCD is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# CL2QCD is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.

add_library(buffers
    buffer.cpp
    prng_buffer.cpp
    su3.cpp
    spinor.cpp
    gaugemomentum.cpp
    3x3.cpp
    halo_update.cpp
    su3vec.cpp
)

target_link_libraries(buffers
    hardware
    code
)

#
# Definition of tests
#
add_unit_test(NAME hardware/buffers/buffer        LIBRARIES buffers crypto)
add_unit_test(NAME hardware/buffers/prng_buffer   LIBRARIES buffers crypto)
add_unit_test(NAME hardware/buffers/su3           LIBRARIES buffers crypto)
add_unit_test(NAME hardware/buffers/plain         LIBRARIES buffers crypto)
add_unit_test(NAME hardware/buffers/spinor        LIBRARIES buffers crypto)
add_unit_test(NAME hardware/buffers/gaugemomentum LIBRARIES buffers crypto)
add_unit_test(NAME hardware/buffers/3x3           LIBRARIES buffers crypto)
add_unit_test(NAME hardware/buffers/su3vec        LIBRARIES buffers crypto)
