#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Author: Matthias Bach <bach@compeng.uni-frankfurt.de>
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


import matplotlib.pyplot as plt
from pylab import recfromtxt
import optparse
import operator

#
# MAIN
#
if __name__ == "__main__":

	parser = optparse.OptionParser( usage='Usage: %prog [options] datafile reference', description='Compare two sets of heatbath output data' )
	parser.add_option("-f","--file", dest="filename", help="write comparison results to files")
	parser.add_option("-e","--format", dest="format", help="format for the result file", default="pdf")

	( args, files ) = parser.parse_args()

	if len( files ) != 2:
		parser.error( 'Please specify two data files' )
		exit -1

	data = recfromtxt(files[0], usecols=(0,1), dtype=[('iter',int),('kappa',float)])
	reference = recfromtxt(files[1], usecols=(0,1), dtype=[('iter',int),('kappa',float)])



	plt.figure(1)
	plt.title('Comparison of Gauge Variables')
	plt.ylabel('Kappa')
	plt.xlabel('Iteration')

	plt.plot( data.iter, data.kappa, label='Data' )
	plt.plot( reference.iter, reference.kappa, label='Reference' )

	plt.legend()

	if args.filename:
		plt.savefig( args.filename + '.comparison.' + args.format )



	difference = data.kappa - reference.kappa

	plt.figure(2)
	plt.title('Difference of Gauge Variables')
	plt.ylabel('Kappa')
	plt.xlabel('Iteration')

	plt.plot( data.iter, difference, label='Difference' )

	plt.legend()

	if args.filename:
		plt.savefig( args.filename + '.difference.' + args.format )
	else:
		plt.show()

