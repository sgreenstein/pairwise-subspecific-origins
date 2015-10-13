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
import logging

from time import clock
from collections import OrderedDict, Counter

# --- CONSTANTS --- #
# integer representations of subspecies
NO_SPECIES = 0
SPECIES_TO_INT = OrderedDict({'None': NO_SPECIES, 'dom': 0o001, 'mus': 0o010, 'cas': 0o100})
# identify 'homosubspecific' allele combinations
INTRA_SUBSPECIFIC = [v + v for v in SPECIES_TO_INT.itervalues()] + [-1 * (v + v) for v in SPECIES_TO_INT.itervalues()]

# integer representations of chromosomes
CHROMO_TO_INT = OrderedDict()
for x in xrange(1, 21):
    CHROMO_TO_INT[str(x)] = x
# CHROMO_TO_INT['X'] = 20
CHROMO_TO_INT['Y'] = 21
CHROMO_TO_INT['MT'] = 22

INT_TO_CHROMO = {v: k for k, v in CHROMO_TO_INT.iteritems()}

# chromosome sizes (mm9 coordinates)
CHROMO_SIZES = [197195432, 181748087, 159599783, 155630120, 152537259,
                149517037, 152524553, 131738871, 124076172, 129993255,
                121843856, 121257530, 120284312, 125194864, 103494974,
                98319150, 95272651, 90772031, 61342430, 166650296,
                91744698, 16299]


class TwoLocus:
    def __init__(self, in_path=None, chrom_sizes=None):
        """ Load a databse of pairwise labels for a collection of samples.
        :param in_path: default path to database of pre-computed intervals
        """
        self.path = in_path if in_path else os.getcwd()
        self.available = self.list_available_strains()
        self.sizes = chrom_sizes if chrom_sizes else CHROMO_SIZES
        self.offsets = np.cumsum([0] + chrom_sizes)

    def genome_index(self, chromosome, position):
        """ Converts chromosome and position to a single position in a coordinate system that covers
        the whole genome.  Chromosomes are just concatenated in karyotype order.
        :param chromosome: integer representation of chromosome
        :param position: position on chromosome
        :return: integer denoting chromosome and position
        """
        return self.offsets[chromosome - 1] + position

    def chrom_and_pos(self, index):
        """ Converts a single genome position to a chromosome and position
        :param index: integer denoting chromosome and position
        :return: string representation of chromosome, position on chromosome
        """
        chromo_num = 1
        while index > self.sizes[chromo_num]:
            index -= self.sizes[chromo_num]
            chromo_num += 1
        return INT_TO_CHROMO[chromo_num], index

    def list_available_strains(self, in_path=None):
        """ Given path to a databse of pre-computed subspecies origin intervals,
        list the strains that are available.
        :param in_path: path to directory containing pre-computed intervals
        :return: list of available strains
        """
        in_path = self.path if not in_path else in_path
        sources = []
        for s in glob.glob(os.path.join(in_path, "*_sources.npy")):
            sources.append(os.path.basename(s).replace("_sources.npy", ""))
        intervals = []
        for s in glob.glob(os.path.join(in_path, "*_intervals.npy")):
            intervals.append(os.path.basename(s).replace("_intervals.npy", ""))
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
        # add null interval to end of each chromosome
        for chromosomes in strains.itervalues():
            for chromosome, intervals in chromosomes.iteritems():
                if intervals[-1] < self.sizes[chromosome - 1]:
                    intervals.append((NO_SPECIES, self.sizes[chromosome - 1]))
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
        # prox::dist is upper triangle; dist::prox is lower triangle
        sources_new = np.triu(sources) + -1 * np.tril(sources, k=-1).astype(np.int8)
        np.save(os.path.join(in_path, strain_name + "_intervals.npy"), intervals)
        np.save(os.path.join(in_path, strain_name + "_sources.npy"), sources_new)

    @classmethod
    def make_elementary_intervals(cls, interval_lists):
        """ Given a list of lists of interval endpoints, find minimal set of intervals
        which can cover them all without breaking any in two.
        :param interval_lists: list of lists; elements of inner list are endpoints of genomic intervals
        :return: endpoints of the 'elementary intervals' induced by the input intervals
        """
        elem_intervals = []
        while interval_lists:
            elem_intervals.append(min(intervals[0] for intervals in interval_lists))
            i = 0
            while i < len(interval_lists):
                # for i, intervals in enumerate(interval_lists):
                intervals = interval_lists[i]
                if intervals[0] == elem_intervals[-1]:
                    intervals.pop(0)
                if not intervals:
                    interval_lists.pop(i)
                else:
                    i += 1
        return elem_intervals

    def calculate_pairwise_frequencies(self, strain_names=None, verbose=False):
        """ For every locus pair and every label pair, count the number of strains which have those
        labels at those pairs of loci.
        :param strain_names: list of strain names to analyze (must be a subset of the output from preprocess())
        :param verbose: spit out progress messages to stderr
        """
        if verbose:
            logging.basicConfig(level=logging.INFO)
        logging.info('Loading intervals')

        in_path = self.path if self.path else os.getcwd()
        interval_lists = []
        interval_dict = {}
        for strain_name in strain_names:
            interval_lists.append(list(np.load(os.path.join(in_path, strain_name + "_intervals.npy"))))
            interval_dict[strain_name] = np.copy(interval_lists[-1])

        # compute elementary intervals
        logging.info("Computing elementary intervals...")
        elem_intervals = self.make_elementary_intervals(interval_lists)
        logging.info("%d total intervals", len(elem_intervals))

        # dictionary of matrices, one for each combo of two sources. i.e. 1 matrix has the counts for the dom-mus combo
        # elementary intervals are along the axes
        # Each element contains the number of samples that have that matrix's combo at that pair of intervals
        source_counts = {}
        for species1, source1 in SPECIES_TO_INT.iteritems():
            for species2, source2 in SPECIES_TO_INT.iteritems():
                source_counts[source1 + source2] = np.zeros([len(elem_intervals), len(elem_intervals)],
                                                            dtype=np.int16)
                source_counts[-1 * (source1 + source2)] = np.zeros([len(elem_intervals), len(elem_intervals)],
                                                                   dtype=np.int16)
        # fill in matrices
        logging.info("Counting incidence of subspecies pairs...")
        for strain_name in strain_names:
            logging.info("\t-- %s", strain_name)
            start = clock()
            intervals = interval_dict[strain_name]
            sources = np.load(os.path.join(in_path, strain_name + "_sources.npy"))
            # map this strain's intervals onto the elementary intervals
            breaks = np.searchsorted(elem_intervals, intervals)
            for row, row_end in enumerate(intervals):
                row_ind = range(breaks[row], breaks[row] + 1)
                for col, col_end in enumerate(intervals):
                    source = sources[row, col]
                    col_ind = range(breaks[col], breaks[col] + 1)
                    if source in source_counts:
                        for i in col_ind:
                            source_counts[source][row_ind, i] += 1
            end = clock()
            logging.info("\t\t... in %.2f seconds", end - start)

        return source_counts, elem_intervals

    def calculate_genomic_area(self, counts, intervals):
        """
        Compute the total genomic 'area' occupied by each combination of subspecies.
        :param counts: dictionary of incidence matrices, one per subspecies combo
        :param intervals: the 'elementary intervals' over which the counts were computed
        """
        # compute area of each cell in the interval grid
        intervals = np.array([0] + intervals, dtype=np.float32) / 1.0e6
        areas = np.zeros([len(intervals) - 1, len(intervals) - 1], dtype=np.float32)
        for row in range(1, len(intervals)):
            for col in range(row, len(intervals)):
                areas[row - 1, col - 1] = (intervals[row] - intervals[row - 1]) * (intervals[col] - intervals[col - 1])
                if col > row:
                    areas[col - 1, row - 1] = areas[row - 1, col - 1]

        areas_masked = OrderedDict()
        denom = np.sum(np.array(self.sizes) / 1.0e6) ** 2
        for combo, vals in counts.iteritems():
            factor = 1
            areas_masked.update({combo: np.sum((vals > 0) * areas * factor) / denom})
        return areas_masked

    def sources_at_point_pair(self, chrom1, pos1, chrom2, pos2, strain_names):
        """ Prints the range of the 2D interval and the counts of subspecific combos at 2 loci in the genome
        :param chrom1: chromosome of one locus
        :param pos1: position of one locus
        :param chrom2: chromosome of another locus
        :param pos2: position of another locus
        :param strain_names: list of strain names to analyze
        """
        coords = [self.genome_index(CHROMO_TO_INT[str(chrom1)], pos1), self.genome_index(CHROMO_TO_INT[str(chrom2)], pos2)]
        mins = [0] * 2
        maxes = [np.sum(self.sizes)] * 2
        coords.sort()
        source_counts = Counter()
        for strain_name in strain_names:
            intervals = np.load(strain_name + '_intervals.npy')
            sources = np.load(strain_name + '_sources.npy')
            print intervals
            print sources
            # find interval containing each location
            i = 0
            interval_indices = [None, None]
            print 'coords', coords
            for loc_num in xrange(2):
                while intervals[i] < coords[loc_num]:
                    i += 1
                if i > 0:
                    mins[loc_num] = max(mins[loc_num], intervals[i-1])
                maxes[loc_num] = min(maxes[loc_num], intervals[i])
                interval_indices[loc_num] = i
            source_counts[sources[interval_indices[0], interval_indices[1]]] += 1
        print 'Proximal range: Chr %s:%d - Chr %s:%d'\
              % tuple(cp for cp in self.chrom_and_pos(mins[0]) + self.chrom_and_pos(maxes[0]))
        print 'Distal range: Chr %s:%d - Chr %s:%d'\
              % tuple(cp for cp in self.chrom_and_pos(mins[1]) + self.chrom_and_pos(maxes[1]))
        print source_counts
        for species1, species1_int in SPECIES_TO_INT.iteritems():
            for species2, species2_int in SPECIES_TO_INT.iteritems():
                if species1_int + species2_int in source_counts:
                    print species1, '+', species2, source_counts[species1_int + species2_int]


def main():
    """ Run some tests with a dummy file, overriding chromosome lengths locally for sake of testing. """
    x = TwoLocus(chrom_sizes=[200000000] * len(CHROMO_TO_INT))
    x.sources_at_point_pair(19, 3500000, 19, 10000000, ['CASTEiJ', 'AJ'])
    exit()

    x = TwoLocus(chrom_sizes=[20e6, 20e6])
    x.preprocess(["test.csv"])
    rez = x.calculate_pairwise_frequencies(["A"])

    combos = OrderedDict()
    for s1, i1 in SPECIES_TO_INT.iteritems():
        for s2, i2 in SPECIES_TO_INT.iteritems():
            combos[-1 * (i1 + i2)] = s1 + " :: " + s2
            combos[i1 + i2] = s2 + " :: " + s1

    areas = x.calculate_genomic_area(rez[0], rez[1])
    total = 0.0
    for code, combo in combos.iteritems():
        if code in areas:
            print "\t{:15s}({:4d}):{:1.5f}".format(combo, code, areas[code])
            total += areas[code]
    print "\t{:21s}:{:1.5f}".format("Total", total)

    sys.exit(1)
    for code, combo in combos.iteritems():
        print "\n", rez[1]
        print "\t{} ({}):\n{}".format(combo, code, rez[0][code])


if __name__ == "__main__":
    main()
