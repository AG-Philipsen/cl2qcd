# Build specifications for the contents of the physics package
#
# Copyright (c) 2012,2013 Matthias Bach
# Copyright (c) 2012,2014,2015 Christopher Pinke
# Copyright (c) 2016 Francesca Cuteri
# Copyright (c) 2018 Alessandro Sciarra
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

add_library(prng
    prng.cpp
    utilities.cpp
)
target_link_libraries(prng
    ranlux
    hardware
    host_functionality
)

add_library(physics
    sources.cpp
)
add_modules(physics
    lattices
    fermionmatrix
    algorithms
    observables
)
target_link_libraries(physics
    hardware
    prng
    host_functionality
    ildg_io
)

#
# Definition of tests
#
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/prngstate_brokenTag prngstate_brokenTag COPYONLY )

add_unit_test(NAME physics/prng      LIBRARIES prng meta hardware)
add_unit_test(NAME physics/sources   LIBRARIES physics)
add_unit_test(NAME physics/utilities LIBRARIES physics)
