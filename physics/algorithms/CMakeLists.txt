# Build specifications for the contents of the algorithms package
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.

#
# Collect all code in this part of the source tree
#

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
	alg_remez.C
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

add_unit_test(physics::algorithms::<fermion_force> fermion_force_test)
add_unit_test(physics::algorithms::<fermion_force_staggered> fermion_force_staggered_test)
add_unit_test(physics::algorithms::<forces> forces_test)
add_unit_test(physics::algorithms::<molecular_dynamics> molecular_dynamics_test)
add_unit_test(physics::algorithms::<rational_approximation> rational_approximation_test ${CMAKE_CURRENT_SOURCE_DIR}/rational_approximation_test_input)
add_unit_test(physics::algorithms::solver::<solver_shifted> solver_shifted_test)
add_unit_test(physics::algorithms::<find_minmax_eigenvalue> find_minmax_eigenvalue_test)
add_unit_test(physics::algorithms::<metropolis> metropolis_test)
