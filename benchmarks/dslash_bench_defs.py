# Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
# Matthias Bach, Christian Schaefer, Stefano Lottini, Alessandro Sciarra
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

#name of the executable
executable = 'dslash_benchmark'

#debugging options:
#print more information
debug = 0
#print executable-output to stdout
stdout = 1
#save result file under a different name for backup
backup = 1

default_space_dims = [16, 24, 32, 48]
default_time_dims = [4, 8, 12, 16, 20, 24, 28, 32]

input_glob = """#global settings
prec=64
use_gpu=true
num_dev=2

startcondition=cold

#fermion settings
fermact=TWISTEDMASS
kappa=0.05
mu=0.2
corr_dir=3
ThetaT=1.
use_evenodd=yes

#this controls hows many bencharking steps are performed
#   NOTE: times 2, for each step EVEN and ODD is performed!!
hmcsteps=5000

#variable settings depending on test
NS={0}

#variable settings depending on test
NT={1}

"""
