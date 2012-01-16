#!/usr/bin/env python
# coding=utf8
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

import matplotlib.pyplot as plt
import numpy as np
import sys
import re

def main(datafile):
	runs = []

	for line in open(datafile, 'r'):
		match = re.match(r'## NSPACE:\s+(\d+)', line)
		if match:
			nspace = int(match.group(1))
		match = re.match(r'## NTIME:\s+(\d+)', line)
		if match:
			ntime = int(match.group(1))
		floatpattern = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
		match = re.match(r'\s*dslash_eoprec\s+\d+\s+\d+\s+\d+\s+\d+\s+(' + floatpattern + ')\s(' + floatpattern + ')', line)
		if match:
			print match
			bandwidth = float(match.group(1))
			gflops = float(match.group(5))
			runs.append((nspace, ntime, bandwidth, gflops))

	# dump data to cli
	print 'NTIME   NSPACE    GB/S   GFLOPS'
	for run in runs:
		print '{0[0]:>5}   {0[1]:>5}   {0[2]:>5}   {0[3]:>5}'.format(run)

	runs = np.array(runs)

	xpos = map(lambda (x,y): int(x)**3 * int(y), runs[:,0:2])
	xtics = map(lambda (x,y): '{0}^3x{1}'.format(int(x), int(y)), runs[:,0:2])

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	line1 = ax1.plot(xpos, runs[:,2], 'r.', label='Bandwidth')
	line2 = ax2.plot(xpos, runs[:,3], 'b.', label='Gflops')
	ax1.set_title('Dslash Performance')
	ax1.set_xticks(xpos)
	ax1.set_xticklabels(xtics, rotation=90)
	ax1.set_xlabel('Lattice Size')
	ax1.set_ylabel('Bandwidth GB/s')
	ax2.set_ylabel('Gflops')
	ax2.set_ylim(ax1.get_ylim())
	fig.legend((line1, line2),
	           ('Bandwidth', 'Gflops'))

	plt.show()

if __name__ == '__main__':
	if len(sys.argv) > 1:
		datafile = sys.argv[1]
	else:
		datafile = 'dslash_benchmark_profiling_data'
	main(datafile)
