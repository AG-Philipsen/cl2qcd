# Build specifications for the contents of the code package
#
# Copyright (c) 2012,2013 Matthias Bach
# Copyright (c) 2013-2016 Christopher Pinke
# Copyright (c) 2014,2018 Alessandro Sciarra
# Copyright (c) 2016 Francesca Cuteri
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

add_library(code
    opencl_module.cpp
    gaugefield.cpp
    prng.cpp
    heatbath.cpp
    kappa.cpp
    correlator.cpp
    correlator_staggered.cpp
    real.cpp
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
    buffers
    klepsydra
    einhard
    hardwareTestUtilities
)

add_library(kernelTestUtilities
    testUtilities.cpp
    kernelTester.cpp
    GaugefieldTester.cpp
    SpinorTester.cpp
    GaugemomentumTester.cpp
    SpinorStaggeredTester.cpp
)

target_link_libraries(kernelTestUtilities
    ${Boost_LIBRARIES}
    host_functionality
    hardware
    geometry
)

add_subdirectory(miscellaneousTests)

# default random number seed
configure_file( ${CMAKE_SOURCE_DIR}/rand_seeds rand_seeds COPYONLY )

#
# Definition of tests
#
add_unit_test(            NAME hardware/code/KernelTester            LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/gaugefield              LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/spinors                 LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/spinors_merged_kernels  LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/fermions                LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/fermions_merged_kernels LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/spinors_staggered       LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/correlator              LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/real                    LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/gaugemomentum           LIBRARIES kernelTestUtilities) #<--these contain only dummy tests
add_unit_test(CREATE_ONLY NAME hardware/code/molecular_dynamics      LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/fermions_staggered      LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/correlator_staggered    LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/complex                 LIBRARIES kernelTestUtilities)

add_unit_test(ADD_ONLY    NAME hardware/code/gaugefield_CPU                    COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/gaugefield_CPU_REC12              COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/gaugefield_GPU                    COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/gaugefield_GPU_REC12              COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/spinors_CPU                       COMMAND_LINE_OPTIONS -- --useGPU=false                         )
add_unit_test(ADD_ONLY    NAME hardware/code/spinors_GPU                       COMMAND_LINE_OPTIONS -- --useGPU=true                          )
add_unit_test(ADD_ONLY    NAME hardware/code/spinors_merged_kernels_CPU        COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/spinors_merged_kernels_GPU        COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_CPU                      COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_CPU_REC12                COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_GPU                      COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_GPU_REC12                COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_merged_kernels_CPU       COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_merged_kernels_CPU_REC12 COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_merged_kernels_GPU       COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_merged_kernels_GPU_REC12 COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/spinors_staggered_CPU             COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/spinors_staggered_GPU             COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/correlator_CPU                    COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/correlator_GPU                    COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/real_CPU                          COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/real_GPU                          COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/gaugemomentum_CPU                 COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/gaugemomentum_CPU_REC12           COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/gaugemomentum_GPU                 COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/gaugemomentum_GPU_REC12           COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/molecular_dynamics_CPU            COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/molecular_dynamics_CPU_REC12      COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/molecular_dynamics_GPU            COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/molecular_dynamics_GPU_REC12      COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_staggered_CPU            COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_staggered_CPU_REC12      COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_staggered_GPU            COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/fermions_staggered_GPU_REC12      COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true )
add_unit_test(ADD_ONLY    NAME hardware/code/correlator_staggered_CPU          COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/correlator_staggered_GPU          COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/complex_CPU                       COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(ADD_ONLY    NAME hardware/code/complex_GPU                       COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
