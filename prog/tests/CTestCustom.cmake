# Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
# Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

# most tests require double precision
# some will explicitly require single precision
set( USE_DOUBLE_PRECISION @USE_DOUBLE_PRECISION@ )
if( USE_DOUBLE_PRECISION )
	set( CTEST_CUSTOM_TESTS_IGNORE
	     ${CTEST_CUSTOM_TESTS_IGNORE}
	     Heatbath5_CPU_SP
	     Heatbath6_CPU_SP_REC12
	     Heatbath7_GPU_SP
	     Heatbath8_GPU_SP_REC12
	   )
else( USE_DOUBLE_PRECISION )
	set( CTEST_CUSTOM_TESTS_IGNORE
	     ${CTEST_CUSTOM_TESTS_IGNORE}
	     Inverter1_CPU_TM_EO
	     Inverter2_CPU_TM_EO_REC12
	     Inverter3_GPU_TM_EO
	     Inverter4_GPU_TM_EO_REC12
	     Inverter5_CPU_TM
	     Inverter6_CPU_TM_REC12
	     Inverter7_GPU_TM
	     Inverter8_GPU_TM_REC12

	     Heatbath1_CPU_DP
	     Heatbath2_CPU_DP_REC12
	     Heatbath3_GPU_DP
	     Heatbath4_GPU_DP_REC12
	   )
endif( USE_DOUBLE_PRECISION )
