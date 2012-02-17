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

	raw_times = []
	raw_gflops = []
	nextRun = False

	for line in open(datafile, 'r'):
		match = re.match(r'^\[\d\d:\d\d:\d\d\] INFO: ## NSPACE:\s+(\d+)', line)
		if match:
			nspace = int(match.group(1))
			nextRun = True
		match = re.match(r'^\[\d\d:\d\d:\d\d\] INFO: ## NTIME:\s+(\d+)', line)
		if match:
			ntime = int(match.group(1))
			nextRun = True

		if nextRun and len(raw_times) > 1:
			# always ignore first run TODO
			times = np.array(raw_times[1:]);
			gflops = np.array(raw_gflops[1:]);

			runs.append((nspace, ntime, np.mean(times), np.std(times), np.mean(gflops), np.std(gflops)));

			raw_times = []
			raw_gflops = []
			nextRun = False

		floatpattern = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
		match = re.match(r'^\[\d\d:\d\d:\d\d\] INFO: \w+ completed in (\d+) ms @ (' + floatpattern + ') Gflops.', line)
		if match:
			raw_times.append(int(match.group(1)))
			raw_gflops.append(float(match.group(2)))

	# dump data to cli
	print 'NTIME   NSPACE     ms (std)  GFLOPS (std)'
	for run in runs:
		print '{0[0]:>5}   {0[1]:>5}   {0[2]:>5.4} {0[3]:>5.2}   {0[4]:>5.4} {0[5]:>5.2}'.format(run)

	runs = np.array(runs)

	xpos = map(lambda (x,y): int(x)**3 * int(y), runs[:,0:2])
	xtics = map(lambda (x,y): '{0}^3x{1}'.format(int(x), int(y)), runs[:,0:2])

	fig = plt.figure()
	ax1 = fig.add_subplot(111)
	line1 = ax1.plot(xpos, runs[:,4], 'b.', label='Gflops')
	ax1.set_title('Inverter Performance')
	ax1.set_xticks(xpos)
	ax1.set_xticklabels(xtics, rotation=90)
	ax1.set_xlabel('Lattice Size')
	ax1.set_ylabel('Gflops')
	ax1.set_ylim(bottom=1)
	fig.legend((line1),
	           ('Gflops'))

	plt.show()

if __name__ == '__main__':
	datafile = sys.argv[1]
	main(datafile)

