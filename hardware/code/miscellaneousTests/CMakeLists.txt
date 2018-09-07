# Miscellaneous tests
#
# Copyright (c) 2014-2016 Christopher Pinke
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

add_unit_test(CREATE_ONLY NAME hardware/code/miscellaneousTests/Calc_staple        LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/miscellaneousTests/gaugefield_convert LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/miscellaneousTests/local_Q            LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/miscellaneousTests/su2su3extend       LIBRARIES kernelTestUtilities)
add_unit_test(CREATE_ONLY NAME hardware/code/miscellaneousTests/geometry           LIBRARIES kernelTestUtilities)

add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/Calc_staple_CPU        COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/Calc_staple_CPU_REC12  COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/Calc_staple_GPU        COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/Calc_staple_GPU_REC12  COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/gaugefield_convert_CPU COMMAND_LINE_OPTIONS -- --useGPU=false )
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/gaugefield_convert_GPU COMMAND_LINE_OPTIONS -- --useGPU=true )
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/local_Q_CPU            COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=false)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/local_Q_CPU_REC12      COMMAND_LINE_OPTIONS -- --useGPU=false --useReconstruct12=true)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/local_Q_GPU            COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=false)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/local_Q_GPU_REC12      COMMAND_LINE_OPTIONS -- --useGPU=true  --useReconstruct12=true)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/su2su3extend_CPU       COMMAND_LINE_OPTIONS -- --useGPU=false)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/su2su3extend_GPU       COMMAND_LINE_OPTIONS -- --useGPU=true)
add_unit_test(   ADD_ONLY NAME hardware/code/miscellaneousTests/geometry) # <-- this contains only dummy test and need to be worked on
