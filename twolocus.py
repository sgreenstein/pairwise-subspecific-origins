"""
File: twolocus.py
Authors: Seth Greenstein, Andrew P Morgan
Purpose:
		Treat the genome as a series of intervals, some of which have a known subspecific origin.
		Over subsets of all samples, find the frequency of all pairwise combinations of origins for intervals.
		E.g. find the number of samples that are domesticus at interval A and musculus at interval B
"""

import os
import sys
import numpy as np
import csv
import string
import glob

from time import clock
from collections import OrderedDict

## --- CONSTANTS --- ##
# integer representations of subspecies
NO_SPECIES = 0
SPECIES_TO_INT = OrderedDict({'None': NO_SPECIES, 'dom': 0o001, 'mus': 0o010, 'cas': 0o100})
# identify 'homosubspecific' allele combinations
INTRA_SUBSPECIFIC = [ v+v for k,v in SPECIES_TO_INT.iteritems() ] + [ -1*(v+v) for k,v in SPECIES_TO_INT.iteritems() ]

# integer representations of chromosomes
CHROMO_TO_INT = OrderedDict()
for x in range(1,21):
	CHROMO_TO_INT[ str(x) ] = x
#CHROMO_TO_INT['X'] = 20
CHROMO_TO_INT['Y'] = 21
CHROMO_TO_INT['MT'] = 22

# chromosome sizes (mm9 coordinates)
CHROMO_SIZES = [	197195432, 181748087, 159599783, 155630120, 152537259,
					149517037, 152524553, 131738871, 124076172, 129993255,
					121843856, 121257530, 120284312, 125194864, 103494974,
					98319150,	95272651,  90772031,  61342430, 166650296,
					91744698,		16299 ]
CHROMO_OFFSETS = np.cumsum([0] + CHROMO_SIZES)

class TwoLocus:

	def __init__(self, in_path = None, chrom_sizes = None):
		""" Load a databse of pairwise labels for a collection of samples.
		:param in_path: default path to database of pre-computed intervals
		"""
		self.path = in_path if in_path else os.getcwd()
		self.available = self.list_available_strains()
		self.sizes = chrom_sizes if chrom_sizes else CHROMO_SIZES
		self.offsets = np.cumsum([0] + chrom_sizes) if chrom_sizes else CHROMO_OFFSETS

	def genome_index(self, chromosome, position):
		""" Converts chromosome and position to a single position in a coordinate system that covers
		the whole genome.  Chromosomes are just concanteated in karyotype order.
		:param chromosome: integer representation of chromosome
		:param position: position on chromosome
		:return: integer denoting chromosome and position
		"""
		return self.offsets[chromosome-1] + position


	def list_available_strains(self, in_path = None):
		""" Given path to a databse of pre-computed subspecies origin intervals,
		list the strains that are available.
		:param in_path: path to directory containing pre-computed intervals
		:return: list of available strains
		 """
		in_path = self.path if not in_path else in_path
		sources = []
		for s in glob.glob(os.path.join(in_path, "*_sources.npy")):
			sources.append( os.path.basename(s).replace("_sources.npy","") )
		intervals = []
		for s in glob.glob(os.path.join(in_path, "*_intervals.npy")):
			intervals.append( os.path.basename(s).replace("_intervals.npy","") )
		avail = set(intervals) & set(sources)
		self.available = avail
		return avail


	def preprocess(self, file_list):
		""" Creates the preprocessed matrix for each sample. Run once.
		:param file_list: list of file names of csv files with subspecific origins
		:return: list of the strain names that were processed
		"""
		strain_names = []
		# remove chars from strain names that can't be in filenames (e.g. /)
		valid_chars = "-_.()%s%s" % (string.ascii_letters, string.digits)
		for strain_name, chromosomes in self.parse_csvs(file_list).iteritems():
			strain_name = ''.join(c for c in strain_name if c in valid_chars)
			self.make_incidence_matrix(strain_name, chromosomes)
			strain_names.append(strain_name)
		self.available = self.available | set(strain_names)
		return strain_names


	def parse_csvs(self, file_list):
		""" Parses downloaded haplotypes from the Mouse Phylogeny Viewer.
		:param file_list: list of filenames
		:return: dictionary of {strain name: {chromosome number: list of (subspecies id, interval end)}}
		"""
		# csv headers
		STRAIN = 'strain'
		CHROMOSOME = 'chrom'
		START = 'start'
		END = 'end'
		SUBSPECIES = 'subspecies'
		strains = {}
		for filename in file_list:
			with open(filename) as csvfile:
					reader = csv.DictReader(csvfile)
					for row in reader:
						chromosomes = strains.setdefault(row[STRAIN], OrderedDict())
						intervals = chromosomes.setdefault(CHROMO_TO_INT[row[CHROMOSOME]], [])
						lastEnd = 0 if not intervals else intervals[-1][1]
						# add null interval if there is a gap between intervals with assigned subspecies
						if not int(row[START]) - 1 <= lastEnd <= int(row[START]) + 1:
								intervals.append((NO_SPECIES, int(row[START])))
						intervals.append((SPECIES_TO_INT[row[SUBSPECIES]], int(row[END])))
						lastEnd = row[END]
		# add null interval to end of each chromosome
		for chromosomes in strains.itervalues():
			for chromosome, intervals in chromosomes.iteritems():
				if intervals[-1] < self.sizes[chromosome-1]:
					intervals.append((NO_SPECIES, self.sizes[chromosome-1]))
		return strains


	def make_incidence_matrix(self, strain_name, chromosomes):
		""" Constructs and saves to disk the pairwise combinations of subspecific origins of intervals.
		:param strain_name: string, name of the strain/sample
		:param chromosomes: dictionary with the intervals for each chromosome
		:return: void
		"""
		in_path = self.path if self.path else os.getcwd()
		num_intervals = sum([len(ints) for ints in chromosomes.itervalues()])
		intervals = np.empty(num_intervals, dtype=np.uint32)
		source_array = np.empty(num_intervals, dtype=np.int8)
		sources = np.empty([num_intervals, num_intervals], dtype=np.int8)
		interval_num = 0
		for chromosome, interval_list in sorted(chromosomes.iteritems(), key=lambda x: x[0]):
			for species, end in interval_list:
				intervals[interval_num] = self.genome_index(chromosome, end)
				source_array[interval_num] = species
				interval_num += 1
		for interval_num in xrange(num_intervals):
			sources[interval_num, :] = source_array
		for interval_num in xrange(num_intervals):
			sources[:, interval_num] = (source_array + sources[:, interval_num])
		## prox::dist is upper triangle; dist::prox is lower triangle
		sources_new = np.empty([num_intervals, num_intervals], dtype=np.int8)
		sources_new = np.triu(sources) + -1*np.tril(sources, k = -1)
		np.save(os.path.join(in_path, strain_name + "_intervals.npy"), intervals)
		np.save(os.path.join(in_path, strain_name + "_sources.npy"), sources_new)


	@classmethod
	def make_elementary_intervals(self, interval_lists):
		""" Given a list of lists of interval endpoints, find minimal set of intervals
		which can cover them all without breaking any in two.
		:param interval_lists: list of lists; elements of inner list are endpoints of genomic intervals
		:return: endpoints of the 'elementary intervals' induced by the input intervals
		"""
		elem_intervals = []
		while interval_lists:
			elem_intervals.append( min(intervals[0] for intervals in interval_lists) )
			i = 0
			while i < len(interval_lists):
			#for i, intervals in enumerate(interval_lists):
				intervals = interval_lists[i]
				if intervals[0] == elem_intervals[-1]:
					intervals.pop(0)
				if not intervals:
					interval_lists.pop(i)
				else:
					i += 1
		return elem_intervals


	def calculate_pairwise_frequencies(self, strain_names = None, verbose = False):
		""" For every locus pair and every label pair, count the number of strains which have those
		labels at those pairs of loci.
		:param strain_names: list of strain names to analyze (must be a subset of the output from preprocess())
		:param verbose: spit out progress messages to stderr
		"""
		if verbose:
			sys.stderr.write("Loading intervals...\n")

		in_path = self.path if self.path else os.getcwd()
		interval_lists = []
		interval_dict = {}
		for strain_name in strain_names:
			interval_lists.append( list(np.load(os.path.join(in_path, strain_name + "_intervals.npy"))) )
			interval_dict[strain_name] = np.copy(interval_lists[-1])

		# compute elementary intervals
		if verbose:
			sys.stderr.write("Computing elementary intervals...")
		elem_intervals = self.make_elementary_intervals(interval_lists)
		if verbose:
			sys.stderr.write(" {} total intervals\n".format(len(elem_intervals)))

		# dictionary of matrices, one for each combination of two sources. i.e. 1 matrix has the counts for the dom-mus combo
		# elementary intervals are along the axes
		# Each element contains the number of samples that have that matrix's combo at that pair of intervals
		source_counts = {}
		for species1, source1 in SPECIES_TO_INT.iteritems():
			for species2, source2 in SPECIES_TO_INT.iteritems():
				source_counts[ source1+source2 ] = np.zeros([len(elem_intervals), len(elem_intervals)],
															dtype=np.int16)
				source_counts[ -1*(source1+source2) ] = np.zeros([len(elem_intervals), len(elem_intervals)],
																dtype=np.int16)
		## fill in matrices
		if verbose:
			sys.stderr.write("Counting incidence of subspecies pairs...\n")
		for strain_name in strain_names:
			sys.stderr.write("\t-- {}\n".format(strain_name))
			start = clock()
			intervals = interval_dict[strain_name]
			sources = np.load(os.path.join(in_path, strain_name + "_sources.npy"))
			elem_row = 0
			# map this strain's intervals onto the elementary intervals
			breaks = np.searchsorted(elem_intervals, intervals)
			for row, row_end in enumerate(intervals):
				row_ind = range(breaks[row], breaks[row]+1)
				for col, col_end in enumerate(intervals):
					source = sources[row,col]
					col_ind = range(breaks[col], breaks[col]+1)
					if source in source_counts:
						for i in col_ind:
							source_counts[source][ row_ind,i ] += 1
			end = clock()
			if verbose:
				sys.stderr.write("\t\t... in {} seconds\n".format(end-start))

		return ( source_counts, elem_intervals )


	def calculate_genomic_area(self, counts, intervals):
		"""
		Compute the total genomic 'area' occupied by each combination of subspecies.
		:param counts: dictionary of incidence matrices, one per subspecies combo
		:param intervals: the 'elementary intervals' over which the counts were computed
		"""
		# compute area of each cell in the interval grid
		intervals = np.array([0] + intervals, dtype = np.float32)/1.0e6
		areas = np.zeros([ len(intervals)-1, len(intervals)-1 ], dtype = np.float32)
		for row in range(1, len(intervals)):
			for col in range(row, len(intervals)):
				areas[row-1,col-1] = (intervals[row]-intervals[row-1]) * (intervals[col]-intervals[col-1])
				if col > row:
					areas[col-1,row-1] = areas[row-1,col-1]
					
		areas_masked = OrderedDict()
		# get denominator for genomic area
		genome_len = sum(float(x)/1.0e6 for x in self.sizes)-(3.0*2)
		denom = genome_len**2
		for combo, vals in counts.iteritems():
			factor = 2.0 if combo in INTRA_SUBSPECIFIC else 1.0
			areas_masked.update({ combo: np.sum((vals > 0)*areas*factor)/denom })
		return areas_masked


def main():
	""" Run some tests with a dummy file, overriding chromosome lengths locally for sake of testing. """

	x = TwoLocus(chrom_sizes = [20e6, 20e6])
	x.preprocess(["test.csv"])
	rez = x.calculate_pairwise_frequencies(["A","B","C"])

	combos = OrderedDict()
	for s1, i1 in SPECIES_TO_INT.iteritems():
		for s2, i2 in SPECIES_TO_INT.iteritems():
			combos[ s1 + "_" + s2 ] = -1*(i1+i2)
			combos[ s2 + "_" + s1 ] = i1+i2

	areas = x.calculate_genomic_area(rez[0], rez[1])
	total = 0.0
	for combo, code in combos.iteritems():
		if code in areas:
			print "\t{} ({}): {: 1.5f}".format(combo, code, areas[code])
			total += areas[code]
	print total

	for combo, code in combos.iteritems():
		print "\n", rez[1]
		print "\t{} ({}):\n{}".format(combo, code, rez[0][code])


if __name__ == "__main__":
	main()
