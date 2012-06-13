#!/usr/bin/env python
# coding=utf8
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

import matplotlib.pyplot as plt
import numpy as np
import sys
import re
import argparse
from collections import namedtuple

linestyles = ['r.', 'b.', 'r*', 'b*', 'g.', 'k.', 'r,', 'b,', 'g,', 'k,', 'g*', 'k*']

FileData = namedtuple('FileData', ['label', 'runs', 'xpos'])

def main(datafiles, filelabels, output=None):

	filedatas = []

	for i in xrange(len(datafiles)):
		runs = []

		for line in open(datafiles[i], 'r'):
			match = re.match(r'## NSPACE:\s+(\d+)', line)
			if match:
				nspace = int(match.group(1))
			match = re.match(r'## NTIME:\s+(\d+)', line)
			if match:
				ntime = int(match.group(1))
			floatpattern = r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
			match = re.match(r'\s*dslash_(eo|eoprec)\s+\d+\s+\d+\s+\d+\s+\d+\s+(' + floatpattern + ')\s(' + floatpattern + ')', line)
			if match:
				print match
				bandwidth = float(match.group(2))
				gflops = float(match.group(6))
				runs.append((nspace, ntime, bandwidth, gflops))

		# dump data to cli
		print 'NSPACE   NTIME    GB/S   GFLOPS'
		for run in runs:
			print '{0[0]:>5}   {0[1]:>5}   {0[2]:>5}   {0[3]:>5}'.format(run)

		runs = np.array(sorted(runs, key=lambda p: p[0]**3 * p[1]))
		xpos = map(lambda (x,y): int(x)**3 * int(y), runs[:,0:2])
		filedatas.append(FileData(filelabels[i], runs, xpos))

	xtic_pos = []
	for data in filedatas:
		xtic_pos.extend(map(lambda p, (x, y): (p, x, y), data.xpos, data.runs[:,0:2]))
	xtic_pos.sort(key=lambda e: e[0])
	print xtic_pos
	# make sure tics don't overlapp
	min_delta = 10000
	i = len(xtic_pos) - 1
	while i > 0:
		if xtic_pos[i][0] < xtic_pos[i-1][0] + min_delta:
			del xtic_pos[i-1]
		i -= 1
	print xtic_pos
	xtic_label = map(lambda (p, x, y): '{0}^3x{1}'.format(int(x), int(y)), xtic_pos)
	xtic_pos = map(lambda (p, x, y): p, xtic_pos)

	fig = plt.figure(figsize=(10,4))
	fig.subplots_adjust(bottom=0.28)
	ax1 = fig.add_subplot(111)
	ax2 = ax1.twinx()
	lines = []
	labels = []
	linestyle = 0
	for data in filedatas:
		lines.append(ax1.plot(data.xpos, data.runs[:,2], linestyles[linestyle]))
		labels.append(data.label + ' Bandwidth')
		linestyle += 1
		lines.append(ax2.plot(data.xpos, data.runs[:,3], linestyles[linestyle]))
		labels.append(data.label + ' Gflops')
		linestyle += 1
	ax1.set_title('Dslash Performance')
	ax1.set_xticks(xtic_pos)
	ax1.set_xticklabels(xtic_label, rotation=90)
	ax1.set_xlabel('Lattice Size')
	ax1.set_ylabel('Bandwidth GB/s')
	ax2.set_ylabel('Gflops')
	ax1.set_ylim(bottom=1)
	ax2.set_ylim(ax1.get_ylim())
	fig.legend(lines, labels)

	if output:
		plt.savefig(output)
	else:
		plt.show()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot output from dslash_bench.py runs')
	parser.add_argument('--labels', metavar='LABEL', nargs='*', help='Labels to mark the line from each input file.')
	parser.add_argument('files', metavar='FILE', nargs='*')
	parser.add_argument('-o', '--output', metavar='FILE', default=None, help='File to dump the plot to')
	args = parser.parse_args()

	if args.labels and len(args.files) != len(args.labels):
		print 'Please specify exactly one label per file or none at all'

	if args.files:
		datafiles = args.files
	else:
		datafiles = ['dslash_benchmark_profiling_data']

	if args.labels:
		labels = args.labels
	else:
		labels = args.files

	main(datafiles, labels, args.output)
