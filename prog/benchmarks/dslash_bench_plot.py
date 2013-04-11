#!/usr/bin/env python
# coding=utf8
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

import sys
import argparse
from bench_plot import main

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot output from dslash_bench.py runs')
	parser.add_argument('--labels', metavar='LABEL', nargs='*', help='Labels to mark the line from each input file.')
	parser.add_argument('files', metavar='FILE', nargs='*')
	parser.add_argument('-o', '--output', metavar='FILE', default=None, help='File to dump the plot to')
	parser.add_argument('--metric', default='both', help='Output gflops, gbytes or both')
	parser.add_argument('--notitle', default=False, action='store_true', help='Suppress plot title')
	parser.add_argument('--title', default='Overrelax Performance', help='Title to use for the plot')
	parser.add_argument('--maxSize', type=int, help='Maximum lattice size to plot')
	parser.add_argument('--legend-pos', help='Position for the plot legend')
	parser.add_argument('--offset-lines', action='store_true', default=False, help='Offset lines in x direction to avoid overlapping points')
	parser.add_argument('--grayscale', default=False, action='store_true', help='Create a grayscale plot')
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

	if not args.metric in ('gflops', 'gbytes', 'both'):
		print 'Metric must be gflops, gbytes or both.'
		sys.exit(1)

	main(datafiles, labels, r'dslash_(eo|eoprec)', args.output, args.metric, False if args.notitle else args.title, args.maxSize, args.legend_pos, args.offset_lines, args.grayscale)
