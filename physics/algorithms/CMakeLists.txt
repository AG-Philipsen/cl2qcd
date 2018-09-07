# Build specifications for the contents of the algorithms package
#
# Copyright (c) 2012,2013 Matthias Bach
# Copyright (c) 2013,2014,2018 Alessandro Sciarra
# Copyright (c) 2014 Christopher Pinke
# Copyright (c) 2017 Francesca Cuteri
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

add_library(algorithms
    su3heatbath.cpp
    inversion.cpp
    molecular_dynamics.cpp
    metropolis.cpp
    fermion_force.cpp
    fermion_force_staggered.cpp
    forces.cpp
    integrator.cpp
    hmc.cpp
    alg_remez.cpp
    find_minmax_eigenvalue.cpp
    rational_approximation.cpp
    solver_shifted.cpp
    rhmc.cpp
)
add_modules(algorithms
    solvers
)

target_link_libraries(algorithms
    lattices
    fermionmatrix
    physics
    observables
    optimal
)

#
# Definition of tests
#
add_unit_test(NAME physics/algorithms/fermion_force           LIBRARIES algorithms)
add_unit_test(NAME physics/algorithms/fermion_force_staggered LIBRARIES algorithms)
add_unit_test(NAME physics/algorithms/forces                  LIBRARIES algorithms)
add_unit_test(NAME physics/algorithms/molecular_dynamics      LIBRARIES algorithms)
add_unit_test(NAME physics/algorithms/rational_approximation  LIBRARIES algorithms COMMAND_LINE_OPTIONS ${CMAKE_CURRENT_SOURCE_DIR}/rational_approximation_test_input)
add_unit_test(NAME physics/algorithms/solver_shifted          LIBRARIES algorithms)
add_unit_test(NAME physics/algorithms/find_minmax_eigenvalue  LIBRARIES algorithms)
add_unit_test(NAME physics/algorithms/metropolis              LIBRARIES algorithms)
