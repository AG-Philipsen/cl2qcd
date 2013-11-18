## This file should be placed in the root directory of your project.
## Then modify the CMakeLists.txt file in the root directory of your
## project to incorporate the testing dashboard.
## # The following are required to uses Dart and the Cdash dashboard
##   ENABLE_TESTING()
##   INCLUDE(CTest)
##
## Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
## Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
##
## This file is part of CL2QCD.
##
## CL2QCD is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## CL2QCD is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.

set(CTEST_PROJECT_NAME "clhmc")
set(CTEST_NIGHTLY_START_TIME "01:00:00 UTC")

set(CTEST_DROP_METHOD "http")
set(CTEST_DROP_SITE "code.compeng.uni-frankfurt.de")
set(CTEST_DROP_LOCATION "/dashboard/submit.php?project=clhmc")
set(CTEST_DROP_SITE_CDASH TRUE)
