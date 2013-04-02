# Build specifications for the contents of the algorithms package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(algorithms
	solver.cpp
	heatbath.cpp
	kappa_clover.cpp
	flavour_doublet.cpp
	inversion.cpp
	molecular_dynamics.cpp
	metropolis.cpp
	fermion_force.cpp
	forces.cpp
	integrator.cpp
	hmc.cpp
)
target_link_libraries(algorithms
		lattices
		physics
)

#
# Definition of tets
#

add_unit_test(physics::algorithms::flavour_doublet flavour_doublet_test)
add_unit_test(physics::algorithms::<fermion_force> fermion_force_test)
add_unit_test(physics::algorithms::<forces> forces_test)
