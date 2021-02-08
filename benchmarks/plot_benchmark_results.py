#!/usr/bin/env python3
#
# Copyright (c) 2012-2013 Matthias Bach
# Copyright (c) 2021 Alessandro Sciarra
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.

import argparse
import os
import re
from collections import namedtuple
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
import itertools

use_tex = False
plot_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
plot_markers = ['o', 's', 'p', '*', 'h', 'H', 'D', 'd']
plot_markersizes = [6, 6, 8, 8, 7, 7, 6, 6]
plot_lines = ['solid', 'dotted', 'dashed', 'dashdot']
FileData = namedtuple('FileData', ['label', 'runs', 'xlabels'])

def set_tex_usage(use):
	global use_tex
	use_tex = use
	if use_tex:
		from matplotlib import rc
		rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
		## for Palatino and other serif fonts use:
		#rc('font',**{'family':'serif','serif':['Palatino']})
		rc('text', usetex=True)
		#matplotlib.rcParams['text.usetex'] = use

def get_xtic_label(lattice_sizes):
	if use_tex:
		label = rf'${int(lattice_sizes[0])}^3\times{int(lattice_sizes[1])}$'
	else:
		label = f'{int(lattice_sizes[0])}^3x{int(lattice_sizes[1])}'
	return label

def get_performance_plot_xtics_dictionary(datafiles):
	all_lattice_sizes = []
	for data in datafiles:
		for run in data.runs:
			all_lattice_sizes.append([run[0], run[1]])
	# Remove duplicates and keep sorted by volume
	all_lattice_sizes = sorted(all_lattice_sizes, key=lambda p: p[0]**3 * p[1])
	all_lattice_sizes = list(all_lattice_sizes for all_lattice_sizes,_ in itertools.groupby(all_lattice_sizes))
	# Build up dictionary
	xtics_dictionary = {}
	counter = 1
	for lattice in all_lattice_sizes:
		new_xtic = get_xtic_label(lattice)
		if new_xtic not in xtics_dictionary:
			xtics_dictionary[new_xtic] = counter
			counter += 1
	return xtics_dictionary

def get_plot_style(style):
	color=plot_colors[style%len(plot_colors)]
	marker=plot_markers[style%len(plot_markers)]
	markersize=plot_markersizes[style%len(plot_markersizes)]
	return color, plot_lines[0], marker, markersize

def make_performance_plot(axis, datafiles, plot_title):
	style = 0
	x_positions = get_performance_plot_xtics_dictionary(datafiles)

	for data in datafiles:
		x_values = [x_positions[label] for label in data.xlabels]
		color, line, marker, size = get_plot_style(style)
		axis.plot(x_values, data.runs[:,3], label=data.label, color=color, linestyle=line, marker=marker, markerfacecolor=color, markersize=size)
		style += 1

	right_axis = axis.twinx()
	axis.set_title(plot_title)
	axis.set_xticks(list(x_positions.values()))
	axis.set_xticklabels(list(x_positions.keys()), rotation=60, ha='right', rotation_mode="anchor")
	if use_tex:
		axis.set_xlabel(r'\textsc{Lattice Size}')
		axis.set_ylabel(r'\textsc{GB/s}')
		right_axis.set_ylabel(r'\textsc{GFLOPS}')
	else:
		axis.set_xlabel('Lattice Size')
		axis.set_ylabel('GB/s')
		right_axis.set_ylabel('GFLOPS')

	axis.yaxis.labelpad = 15
	right_axis.yaxis.labelpad = 15
	axis.set_ylim(top=axis.get_ylim()[1]*1.05)
	flops_memory_conversion_factor = datafiles[0].runs[0][2]/datafiles[0].runs[0][3]
	right_axis.set_ylim(top=axis.get_ylim()[1]/flops_memory_conversion_factor)

	legend_ncol = len(datafiles) if len(datafiles)<=7 else 7
	figure = plt.gcf()
	figure.legend(ncol=legend_ncol, loc='upper center')

def make_bandwidth_histogram_plot(axis, all_data, bandwidths):
	from statistics import mean
	average_bandwidths = [mean(data.runs[:,3]) for data in all_data]
	x_values = range(len(all_data))
	y_values = []
	for average, maximum in zip(average_bandwidths, bandwidths):
		y_values.append(average/maximum*100)
	axis.bar(x_values, y_values, color=plot_colors)

	axis.yaxis.tick_right()
	axis.set_ylim(bottom=0, top=100)
	axis.set_yticks([val for val in range(0,101,20)])
	if use_tex:
		axis.set_title(r'{\ttfamily <BW>/BW$_\mathtt{max}$}')
		axis.set_yticklabels([rf'{val}\%' for val in range(0,101,20)])
		xtic_labels = [rf'\textsc{{{xtic}}}' for xtic in xtic_labels]
	else:
		axis.set_title('<BW>/BW_max')
		axis.set_yticklabels([f'{val}%' for val in range(0,101,20)])
		xtic_labels = [data.label for data in all_data]
	axis.set_xticks(x_values)
	axis.set_xticklabels(xtic_labels, rotation=20, ha='right', rotation_mode="anchor")

def check_datafiles_existence(datafiles):
	filenames = []
	for file in datafiles:
		if not os.path.exists(file):
			print(f'File "{file}" was not found, aborting.')
			raise SystemExit
		else:
			filenames.append(file)
	return filenames

def check_provided_BWs(files, bw):
	if args.maxBWs != None:
		if len(args.maxBWs) != len(filenames):
			print(f'Number of BW-values provided mismatch number of files: {len(args.maxBWs)} != {len(filenames)}')
			raise SystemExit

def check_provided_labels(files, labels):
	if args.labels != None:
		if len(args.labels) != len(filenames):
			print(f'Number of labels provided mismatch number of files: {len(args.labels)} != {len(filenames)}')
			raise SystemExit

def read_data_from_all_files(filenames, labels):
	read_data = []
	label_counter = 0
	for filename in filenames:
		runs = []
		with open(filename) as file:
			next(file) #Skip header line
			for line in file:
				row = line.split()
				if row:
					runs.append([float(x) for x in row])

		runs = np.array(sorted(runs, key=lambda p: p[0]**3 * p[1]))
		xlabels = [get_xtic_label([run[0], run[1]]) for run in runs]
		read_data.append(FileData(labels[label_counter], runs, xlabels))
		label_counter += 1
	return read_data

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Plot benchmark results')
	parser.add_argument('files', metavar='FILE', nargs='+')
	parser.add_argument('--title', default='Benchmark performance', help='Title to use for the plot')
	parser.add_argument('--maxBWs', nargs='+', type=float, help='The maximum bandwidth of the device used in each benchmark datafile. Specifying this option, a histogram will be produced next to the plot. It requires as many values as the number of files.')
	parser.add_argument('--labels', nargs='+', type=str, help='The labels to be used in the plot e.g. in the legend. It requires as many values as the number of files.')
	parser.add_argument('--TeX', default=False, action='store_true', help='Use TeX style in plots. WARNING: Slow but nice.')
	args = parser.parse_args()

	set_tex_usage(args.TeX)
	filenames = check_datafiles_existence(args.files)
	check_provided_BWs(filenames, args.maxBWs)
	check_provided_labels(filenames, args.labels)
	data_read_from_files = read_data_from_all_files(filenames, args.labels if args.labels != None else filenames)

	if args.maxBWs == None:
		figure, axis = plt.subplots(figsize=(18, 8))
		figure.subplots_adjust(left=0.1, bottom=0.23, right=0.9, top=0.8, wspace=None, hspace=None)
		make_performance_plot(axis, data_read_from_files, args.title)
	else:
		figure, (axis_left, axis_right) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [4, 1]}, figsize=(18, 8))
		figure.subplots_adjust(left=0.06, bottom=0.15, right=0.94, top=0.88, wspace=0.22, hspace=None)
		make_performance_plot(axis_left, data_read_from_files, args.title)
		make_bandwidth_histogram_plot(axis_right, data_read_from_files, args.maxBWs)

	plt.show()
