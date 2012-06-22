#!/usr/bin/env python
# coding=utf8
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

import sys
import argparse
from bench_plot import main

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot output from heatpath_bench.py runs')
	parser.add_argument('--labels', metavar='LABEL', nargs='*', help='Labels to mark the line from each input file.')
	parser.add_argument('files', metavar='FILE', nargs='*')
	parser.add_argument('-o', '--output', metavar='FILE', default=None, help='File to dump the plot to')
	parser.add_argument('--metric', default='both', help='Output gflops, gbytes or both')
	parser.add_argument('--notitle', default=False, action='store_true', help='Suppress plot title')
	parser.add_argument('--maxSize', type=int, help='Maximum lattice size to plot');
	args = parser.parse_args()

	if args.labels and len(args.files) != len(args.labels):
		print 'Please specify exactly one label per file or none at all'

	if args.files:
		datafiles = args.files
	else:
		datafiles = ['heatbath_benchmark_profiling_data']

	if args.labels:
		labels = args.labels
	else:
		labels = args.files

	if not args.metric in ('gflops', 'gbytes', 'both'):
		print 'Metric must be gflops, gbytes or both.'
		sys.exit(1)

	main(datafiles, labels, r'heatbath_(odd|even)', args.output, args.metric, False if args.notitle else 'Heatbath Performance', args.maxSize)
