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
target_link_libraries(su2su3extend
	optimal
	${Boost_LIBRARIES}
)
add_test(Su2su3extend su2su3extend)

add_executable(geometry geometry.cpp)
target_link_libraries(geometry
	optimal
	${Boost_LIBRARIES}
)
add_test(geometry geometry)

add_executable(gaugefield_convert gaugefield_convert.cpp)
target_link_libraries(gaugefield_convert
	optimal
	${Boost_LIBRARIES}
)
add_test(Gaugefield_Conversion gaugefield_convert)

configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/conf.00200 conf.00200 COPYONLY )

add_executable(calc_staple Calc_staple_test.cpp)
target_link_libraries(calc_staple
	optimal
	${Boost_LIBRARIES}
)

add_test(Calc_staple_1_CPU       calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=false --use_rec12=false)
add_test(Calc_staple_1_CPU_REC12 calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=false --use_rec12=true)
add_test(Calc_staple_1_GPU       calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=true  --use_rec12=false)
add_test(Calc_staple_1_GPU_REC12 calc_staple ${CMAKE_CURRENT_SOURCE_DIR}/staple_input_1 --use_gpu=true  --use_rec12=true)


add_executable(localQ local_Q_test.cpp)
target_link_libraries(localQ
	optimal
	${Boost_LIBRARIES}
)
add_test(LocalQ_1_CPU       localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=false --use_rec12=false)
add_test(LocalQ_1_CPU_REC12 localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=false --use_rec12=true)
add_test(LocalQ_1_GPU       localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=true  --use_rec12=false)
add_test(LocalQ_1_GPU_REC12 localQ ${CMAKE_CURRENT_SOURCE_DIR}/local_Q_input_1 --use_gpu=true  --use_rec12=true)
