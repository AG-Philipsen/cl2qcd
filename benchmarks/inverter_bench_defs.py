#!/usr/bin/env python
# coding=utf8
#
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

#name of the executable
executable = 'inverter'

#debugging options:
#print more information
debug = 0
#print executable-output to stdout
stdout = 1
#save result file under a different name for backup
backup = 0

input_glob = """#global settings
prec=64
use_gpu=true
num_dev=2
enable_profiling=true

startcondition=cold

#fermion settings
fermact=TWISTEDMASS
kappa=0.2
mu=0.02
corr_dir=3
ThetaT=1.

cgmax=100
startcondition=cold
savefrequency=10
fermact=TWISTEDMASS
use_evenodd=yes

hmcsteps=1

solver=CG
# solver=BICGSTAB
"""

#arrays for the different tests, this is not nice, but a quick workaround
#the programm will perform tests with all members of this list if no input-file is given
input_var1 = [#"""
#variable settings depending on test
#NS=16
#"""#,
"""
#variable settings depending on test
NS=24
"""#,
#"""
#variable settings depending on test
#NS=32
#""",
#"""
##variable settings depending on test
#NS=48
#"""
]

input_var2 = [#"""
#variable settings depending on test
#NT=4
#""",
#"""
#variable settings depending on test
#NT=8
#""",
"""
#variable settings depending on test
NT=12
"""#,
#"""
#variable settings depending on test
#NT=16
#""",
#"""
#variable settings depending on test
#NT=20
#""",
#"""
#variable settings depending on test
#NT=24
#""",
#"""
#variable settings depending on test
#NT=28
#""",
#"""
#variable settings depending on test
#NT=32
#"""
]
