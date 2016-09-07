# Miscellaneous tests
#
# Copyright 2014 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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

add_executable(su2su3extend su2su3extend.cpp)
target_link_libraries(su2su3extend kernelTestUtilities )
add_test(Su2su3extend_CPU su2su3extend --use_gpu=false)
add_test(Su2su3extend_GPU su2su3extend --use_gpu=true)

add_executable(gaugefield_convert gaugefield_convert.cpp)
target_link_libraries(gaugefield_convert kernelTestUtilities)
add_test(Gaugefield_Conversion_CPU gaugefield_convert --use_gpu=false )
add_test(Gaugefield_Conversion_GPU gaugefield_convert --use_gpu=true )

#@todo: these tests contain only dummy test and need to be worked on
add_executable(hardware_code_geometry geometry.cpp)
target_link_libraries(hardware_code_geometry 	kernelTestUtilities )
add_test(geometry hardware_code_geometry)

add_executable(calc_staple Calc_staple_test.cpp)
target_link_libraries(calc_staple kernelTestUtilities )
add_test(Calc_staple_1_CPU       calc_staple --use_gpu=false --use_rec12=false)
add_test(Calc_staple_1_CPU_REC12 calc_staple --use_gpu=false --use_rec12=true)
add_test(Calc_staple_1_GPU       calc_staple --use_gpu=true  --use_rec12=false)
add_test(Calc_staple_1_GPU_REC12 calc_staple --use_gpu=true  --use_rec12=true)

add_executable(localQ local_Q_test.cpp)
target_link_libraries(localQ kernelTestUtilities)
add_test(LocalQ_1_CPU       localQ --use_gpu=false --use_rec12=false)
add_test(LocalQ_1_CPU_REC12 localQ --use_gpu=false --use_rec12=true)
add_test(LocalQ_1_GPU       localQ --use_gpu=true  --use_rec12=false)
add_test(LocalQ_1_GPU_REC12 localQ --use_gpu=true  --use_rec12=true)
