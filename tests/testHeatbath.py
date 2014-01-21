#!/usr/bin/env python
# coding=utf8
#
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

from subprocess import *
import re
import sys

# use a rather large tolerance as we don't want to have to do too many iterations
tolerance = .3

def compareFloat(reference, actual, precision):
	return actual > reference * (1.0 - precision / 100) and actual < reference * (1.0 + precision / 100)


def main():
	# Expect first arg to be reference file, remaining args are passed to application

	reference = open(sys.argv[1])
	refval = float(reference.readline())

	subject = Popen(['../heatbath'] + sys.argv[2:])

	subject.wait()
	if subject.returncode == 0:
		#print "Program completed successfully"
		pass
	else:
		print "Program terminated with exit code %i" % ( subject.returncode )
		return subject.returncode

	# Read the output file and check whether the plaquette in the last line matches the reference
	for line in open('gaugeObs.dat'):
		pass
	actual = float(line.split()[1])
	if not compareFloat(actual, refval, tolerance):
		print "Invalid Result! Expected: %f  Got: %f (at %.1f %% tolerance)" % (refval, actual, tolerance)
		return 127

	# no check failed
	return 0

if __name__ == '__main__':
	sys.exit(main())
