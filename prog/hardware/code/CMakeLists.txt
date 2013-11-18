# Build specifications for the contents of the code package
#
# Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(code
		opencl_module.cpp
		gaugefield.cpp
		prng.cpp
		heatbath.cpp
		kappa.cpp
		correlator.cpp
		complex.cpp
		spinors.cpp
		spinors_staggered.cpp
		fermions.cpp
		fermions_staggered.cpp
		gaugemomentum.cpp
		molecular_dynamics.cpp
		buffer.cpp
)
target_link_libraries(code
		meta
		buffers
		klepsydra
		optimal # todo remove me!
		einhard
)

#
# Definition of tets
#

