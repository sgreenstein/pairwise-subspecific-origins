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
import subspecies
# import pyximport
# pyximport.install()
# import subspeciesCython as subspecies
import pickle
from scipy import stats

from time import clock
from collections import OrderedDict, Counter

# integer representations of chromosomes
CHROMO_TO_INT = {}
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
        """ Load a database of pairwise labels for a collection of samples.
        :param in_path: default path to database of pre-computed intervals
        """
        self.path = in_path or os.getcwd()
        self._sample_dict_path = os.path.join(self.path, 'sample_dict.p')
        if not os.path.exists(self._sample_dict_path):
            with open(self._sample_dict_path, 'w+') as fp:
                pickle.dump({}, fp)
        with open(self._sample_dict_path) as fp:
            self.sample_dict = pickle.load(fp)
        self.sizes = chrom_sizes or CHROMO_SIZES
        self.offsets = np.cumsum([0] + self.sizes, dtype=int)

    def genome_index_to_dict(self, index):
        """ Converts a genome position to a dictionary of chromosome and position
        :param index: position in genome according to our coordinates
        :return: {'Chromosome': chrom_num, 'Position': pos_on_chrom}
        """
        chrom_pos = self.chrom_and_pos(index)
        return {'Chromosome': chrom_pos[0], 'Position': chrom_pos[1]}

    def save_sample_dict(self):
        """ Saves an updated sample dictionary to disk
        :param new_dict: {sample name: (interval list, origin list)}
        """
        with open(self._sample_dict_path, 'w+') as fp:
            pickle.dump(self.sample_dict, fp)

    def genome_index(self, chromosome, position):
        """ Converts chromosome and position to a single position in a coordinate system that covers
        the whole genome.  Chromosomes are just concatenated in karyotype order.
        :param chromosome: integer or string representation of chromosome
        :param position: position on chromosome
        :return: integer denoting chromosome and position
        :raises: ValueError if position is invalid
        """
        if type(chromosome) is str:
            try:
                chromosome = CHROMO_TO_INT[chromosome]
            except KeyError:
                raise ValueError('Invalid chromosome')
        if position > self.offsets[chromosome]:
            raise ValueError('Position exceeds chromosome length')
        return self.offsets[chromosome - 1] + position

    def chrom_and_pos(self, index):
        """ Converts a single genome position to a chromosome and position
        :param index: integer denoting chromosome and position
        :return: string representation of chromosome, position on chromosome
        """
        chromo_num = 1
        while index > self.sizes[chromo_num - 1]:
            index -= self.sizes[chromo_num - 1]
            chromo_num += 1
        return INT_TO_CHROMO[chromo_num], index

    def list_available_strains(self):
        """ Given path to a databse of pre-computed subspecies origin intervals,
        list the strains that are available.
        :return: list of available strains
        """
        return [strain for strain in self.sample_dict]

    def is_available(self, strain):
        return strain in self.sample_dict

    def preprocess(self, file_list):
        """ Parses and saves the subspecific origins
        :param file_list: list of csv files with subspecific origin information
        """
        for strain_name, chromosomes in self.parse_csvs(file_list).iteritems():
            self.sample_dict[strain_name] = self.intervals_and_sources(chromosomes)
        self.save_sample_dict()

    def parse_csvs(self, file_list):
        """ Parses downloaded haplotypes from the Mouse Phylogeny Viewer
        :param file_list: list of filenames
        :return: dictionary of {strain name: {chromosome number: list of (subspecies id, interval end)}}
        """
        # csv headers
        STRAIN = 'strain'
        CHROMOSOME = 'chrom'
        START = 'start'
        END = 'end'
        SUBSPECIES = 'subspecies'
        COLOR = 'color'
        COLOR_TO_NAME = ['mus', 'cas', 'dom']
        strains = {}
        for filename in file_list:
            with open(os.path.join(self.path, filename)) as csvfile:
                reader = csv.DictReader(csvfile)
                if SUBSPECIES in reader.fieldnames:
                    def subspec_int(r):
                        return subspecies.to_int(r[SUBSPECIES])
                else:
                    def subspec_int(r):
                        return subspecies.to_int(COLOR_TO_NAME[np.argmax([int(v) for v in r[COLOR].split(' ')])])
                for row in reader:
                    chromosomes = strains.setdefault(row[STRAIN], OrderedDict())
                    intervals = chromosomes.setdefault(CHROMO_TO_INT[row[CHROMOSOME]], [])
                    lastEnd = 0 if not intervals else intervals[-1][1]
                    # add null interval if there is a gap between intervals with assigned subspecies
                    if not int(row[START]) - 1 <= lastEnd <= int(row[START]) + 1:
                        intervals.append((subspecies.NONE, int(row[START])))
                    intervals.append((subspec_int(row), int(row[END])))
        # add null interval to end of each chromosome
        for chromosomes in strains.itervalues():
            for chromosome, intervals in chromosomes.iteritems():
                if intervals[-1] < self.sizes[chromosome - 1]:
                    intervals.append((subspecies.NONE, self.sizes[chromosome - 1]))
        return strains

    def intervals_and_sources(self, chromosomes):
        """ Converts dictionary to lists of intervals and sources
        :param chromosomes: dictionary of {strain name: {chromosome number: list of (subspecies id, interval end)}}
        :return: list of ends of intervals, list of int representations of sources
        """
        num_intervals = sum([len(ints) for ints in chromosomes.itervalues()])
        intervals = np.empty(num_intervals, dtype=np.uint32)
        sources = np.empty(num_intervals, dtype=np.uint8)
        interval_num = 0
        for chromosome, interval_list in sorted(chromosomes.iteritems(), key=lambda x: x[0]):
            for species, end in interval_list:
                intervals[interval_num] = self.genome_index(chromosome, end)
                sources[interval_num] = species
                interval_num += 1
        return intervals, sources

    @staticmethod
    def make_elementary_intervals(interval_lists):
        """ Given a list of lists of interval endpoints, find minimal set of intervals
        which can cover them all without breaking any in two.
        :param interval_lists: list of lists; elements of inner list are endpoints of genomic intervals
        :return: endpoints of the 'elementary intervals' induced by the input intervals
        """
        # copy so we don't destroy the original lists
        interval_lists = [[i for i in il] for il in interval_lists]
        elem_intervals = []
        while interval_lists:
            elem_intervals.append(min(intervals[0] for intervals in interval_lists))
            i = 0
            while i < len(interval_lists):
                intervals = interval_lists[i]
                if intervals[0] == elem_intervals[-1]:
                    intervals.pop(0)
                if not intervals:
                    interval_lists.pop(i)
                else:
                    i += 1
        return elem_intervals

    @profile
    def pairwise_frequencies(self, strain_names, include_unknown=False, verbose=False):
        """ For every locus pair and every label pair, count the number of strains which have those
        labels at those pairs of loci.
        :param strain_names: list of strain names to analyze (must be a subset of the output from preprocess())
        :param verbose: log progress
        :returns {strain combo: matrix of counts}, list of elementary interval ends
        """
        if verbose:
            logging.basicConfig(level=logging.INFO)
        logging.info('Loading intervals...')

        interval_lists = []
        for strain_name in strain_names:
            interval_lists.append(list(self.sample_dict[strain_name][0]))

        # compute elementary intervals
        logging.info("Computing elementary intervals...")
        elem_intervals = self.make_elementary_intervals(interval_lists)
        logging.info("%d total intervals", len(elem_intervals))

        # dictionary of matrices, one for each combo of two sources. i.e. 1 matrix has the counts for the dom-mus combo
        # elementary intervals are along the axes
        # Each element contains the number of samples that have that matrix's combo at that pair of intervals
        source_counts = {}
        for combo in subspecies.iter_combos(include_unknown=include_unknown):
            source_counts[combo] = np.zeros([len(elem_intervals), len(elem_intervals)], dtype=np.int16)
        # fill in matrices
        logging.info("Counting incidence of subspecies pairs...")
        for strain_name in strain_names:
            logging.info("\t-- %s", strain_name)
            start = clock()
            intervals, sources = self.sample_dict[strain_name]
            # map this strain's intervals onto the elementary intervals
            breaks = np.searchsorted(elem_intervals, intervals)
            for row, row_end in enumerate(intervals):
                row_ind = range(breaks[row], breaks[row] + 1)
                for col, col_end in enumerate(intervals):
                    source = subspecies.combine(sources[row], sources[col])
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
        for row in xrange(1, len(intervals)):
            for col in xrange(row, len(intervals)):
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
        coords = [self.genome_index(chrom1, pos1), self.genome_index(chrom2, pos2)]
        mins = [0] * 2
        maxes = [np.sum(self.sizes)] * 2
        coords.sort()
        source_counts = Counter()
        output = {}
        for strain_name in strain_names:
            intervals = self.sample_dict[strain_name][0]
            sources = self.sample_dict[strain_name][1]
            # find interval containing each location
            i = 0
            interval_indices = [None, None]
            # output['coords'] = coords
            for loc_num in xrange(2):
                while intervals[i] < coords[loc_num]:
                    i += 1
                if i > 0:
                    mins[loc_num] = max(mins[loc_num], intervals[i - 1])
                maxes[loc_num] = min(maxes[loc_num], intervals[i])
                interval_indices[loc_num] = i
            source_counts[subspecies.combine(sources[interval_indices[0]], sources[interval_indices[1]])] += 1
        output['Proximal'] = [self.genome_index_to_dict(mins[0]), self.genome_index_to_dict(maxes[0])]
        output['Distal'] = [self.genome_index_to_dict(mins[1]), self.genome_index_to_dict(maxes[1])]
        for proximal in subspecies.iter_subspecies():
            for distal in subspecies.iter_subspecies():
                if subspecies.combine(proximal, distal) in source_counts:
                    output.setdefault(subspecies.to_string(proximal), {})[subspecies.to_string(distal)] = \
                        source_counts[subspecies.combine(proximal, distal)]
        return output

    def interlocus_dependence(self, strain_names):
        """ Performs a chi square test to find interval pairs whose origins are interdependent
        :param strain_names: list of strain names to analyze
        :return: elementary intervals, matrix of chi square values, matrix of p values (both upper triangular)
        """
        combo_count_dict, intervals = self.pairwise_frequencies(strain_names)
        # convert source_counts to matrix combo_counts
        combo_counts = np.empty([len(intervals), len(intervals), subspecies.NUM_SUBSPECIES ** 2], dtype=np.uint16)
        species_counts = np.zeros([len(intervals), subspecies.NUM_SUBSPECIES])
        for i, prox_species in enumerate(subspecies.iter_subspecies()):
            for j, dist_species in enumerate(subspecies.iter_subspecies()):
                counts = combo_count_dict[subspecies.combine(prox_species, dist_species)]
                species_counts[:, i] += np.diag(counts)
                combo_counts[:, :, i * subspecies.NUM_SUBSPECIES + j] = counts
        # compute expected combo frequencies from source frequencies
        combo_expectations = np.zeros([len(intervals), len(intervals), subspecies.NUM_SUBSPECIES ** 2])
        for i in xrange(subspecies.NUM_SUBSPECIES):
            for j in xrange(subspecies.NUM_SUBSPECIES):
                combo_expectations[:, :, i * subspecies.NUM_SUBSPECIES + j] += \
                    np.outer(species_counts[:, i], species_counts[:, j])
        # normalize expectations using the actual total frequency
        sums = np.sum(combo_counts, axis=2)
        old_settings = np.seterr(invalid='ignore')  # ignore division by 0 errors for intervals with no assigned origin
        for i in xrange(subspecies.NUM_SUBSPECIES ** 2):
            combo_expectations[:, :, i] = np.true_divide(combo_expectations[:, :, i], sums)
        np.seterr(**old_settings)
        combo_expectations = np.nan_to_num(combo_expectations)
        # do chi-square test
        chi_sq = np.zeros_like(combo_expectations[:, :, 0])
        p_values = np.ones_like(chi_sq)
        for i in xrange(len(intervals)):
            # only upper triangle is meaningful
            for j in xrange(i + 1, len(intervals)):
                nonzero_expectations = np.where(combo_expectations[i, j])
                chi_sq[i, j], p_values[i, j] = stats.chisquare(
                    combo_counts[i, j][nonzero_expectations], combo_expectations[i, j][nonzero_expectations])
        return intervals, chi_sq, p_values

    @staticmethod
    def _find_interval_bounds(intervals1, index1, intervals2, index2):
        """
        :param intervals1: list of interval ends
        :param index1: index into intervals1
        :param intervals2: list of interval ends
        :param index2: index into intervals2
        :return: start of proximal interval, end of distal interval
        """
        if intervals1[index1] < intervals2[index2]:
            lo = 0 if index1 == 0 else intervals1[index1 - 1]
            hi = intervals2[index2]
        else:
            lo = 0 if index2 == 0 else intervals2[index2 - 1]
            hi = intervals1[index1]
        return lo, hi

    @profile
    def unique_combos(self, strain_names1, strain_names2):
        """ finds combinations at interval pairs that are unique to one set of samples
        :param strain_names1: list of strain names in first set
        :param strain_names2: list of strain names in second set
        :return: json object containing the interval pairs unique to each set of samples
        """
        counts1, intervals1 = self.pairwise_frequencies(strain_names1)
        counts2, intervals2 = self.pairwise_frequencies(strain_names2)
        row1 = 0
        row2 = 0
        output = {}
        while row1 < len(intervals1) and row2 < len(intervals2):
            # only do upper triangle
            col1 = row1 + 1
            col2 = row2 + 1
            while col1 < len(intervals1) and col2 < len(intervals2):
                for combo in subspecies.iter_combos():
                    if counts1[combo][row1, col1] and not counts2[combo][row2, col2]:
                        unique_dict = output.setdefault('A', {})
                    elif not counts1[combo][row1, col1] and counts2[combo][row2, col2]:
                        unique_dict = output.setdefault('B', {})
                    else:
                        unique_dict = None
                    if unique_dict is not None:
                        prox_interval = self._find_interval_bounds(intervals1, row1, intervals2, row2)
                        dist_interval = self._find_interval_bounds(intervals1, col1, intervals2, col2)
                        unique_dict.setdefault(subspecies.to_string(combo), []).append({
                            'Proximal': [
                                self.genome_index_to_dict(prox_interval[0]),
                                self.genome_index_to_dict(prox_interval[1])
                            ],
                            'Distal': [
                                self.genome_index_to_dict(dist_interval[0]),
                                self.genome_index_to_dict(dist_interval[1])
                            ]
                        })
                if intervals1[col1] < intervals2[col2]:
                    col1 += 1
                elif intervals1[col1] > intervals2[col2]:
                    col2 += 1
                else:
                    col1 += 1
                    col2 += 1
            if intervals1[row1] < intervals2[row2]:
                row1 += 1
            elif intervals1[row1] > intervals2[row2]:
                row2 += 1
            else:
                row1 += 1
                row2 += 1
        return output


def main():
    """ Run some tests with a dummy file, overriding chromosome lengths locally for sake of testing.
    """
    tl = TwoLocus(in_path='/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    classical = [s for s in
                 ["129P1/ReJ",# "129P3/J", "129S1SvlmJ", "129S6", "129T2/SvEmsJ", "129X1/SvJ", "A/J", "A/WySnJ",
                  # "AEJ/GnLeJ", "AEJ/GnRk", "AKR/J", "ALR/LtJ", "ALS/LtJ", "BALB/cByJ", "BALB/cJ", "BDP/J", "BPH/2J",
                  # "BPL/1J", "BPN/3J", "BTBR T<+>tf/J", "BUB/BnJ", "BXSB/MpJ", "C3H/HeJ", "C3HeB/FeJ", "C57BL/10J",
                  # "C57BL/10ScNJ", "C57BL/10SAAAJ", "C57BL/6CR", "C57BL/6J", "C57BL/6NCI", "C57BL/6Tc", "C57BLKS/J",
                  # "C57BR/cdJ", "C57L/J", "C58/J", "CBA/CaJ", "CBA/J", "CE/J", "CHMU/LeJ", "DBA/1J", "DBA/1LacJ",
                  # "DBA/2DeJ", "DBA/2HaSmnJ", "DBA/2J", "DDK/Pas", "DDY/JclSidSeyFrkJ", "DLS/LeJ", "EL/SuzSeyFrkJ",
                  # "FVB/NJ", "HPG/BmJ", "I/LnJ", "IBWSP2", "IBWSR2", "ICOLD2", "IHOT1", "IHOT2", "ILS", "ISS", "JE/LeJ",
                  # "KK/HlJ", "LG/J", "LP/J", "LT/SvEiJ", "MRL/MpJ", "NOD/ShiLtJ", "NON/ShiLtJ", "NONcNZO10/LtJ",
                  # "NONcNZO5/LtJ", "NOR/LtJ", "NU/J", "NZB/BlNJ", "NZL/LtJ", "NZM2410/J", "NZO/HlLtJ", "NZW/LacJ", "P/J",
                  # "PL/J", "PN/nBSwUmabJ", "RF/J", "RHJ/LeJ", "RIIIS/J", "RSV/LeJ", "SB/LeJ", "SEA/GnJ", "SEC/1GnLeJ",
                  # "SEC/1ReJ", "SH1/LeJ", "SI/Col Tyrp1 Dnahc11/J", "SJL/Bm", "SJL/J", "SM/J", "SSL/LeJ", "ST/bJ",
                  "STX/Le", ]#"SWR/J", "TALLYHO/JngJ", "TKDU/DnJ", "TSJ/LeJ", "YBR/EiJ", "ZRDCT Rax<+>ChUmdJ"]
                 if tl.is_available(s)]
    wild_derived = [s for s in
                    ['22MO',# 'BIK/g', 'BULS', 'BUSNA', 'BZO', 'CALB/RkJ', 'CASA/RkJ', 'CAST/EiJ', 'CIM', 'CKN', 'CKS',
                     # 'CZECHI/EiJ', 'CZECHII/EiJ', 'DCA', 'DCP', 'DDO', 'DEB', 'DGA', 'DIK', 'DJO', 'DKN', 'DMZ', 'DOT',
                     # 'IS/CamRkJ', 'JF1/Ms', 'LEWES/EiJ', 'MBK', 'MBS', 'MCZ', 'MDG', 'MDGI', 'MDH', 'MGA', 'MH',
                     # 'MOLD/RkJ', 'MOLF/EiJ', 'MOLG/DnJ', 'MOR/RkJ', 'MPB', 'MSM/Ms', 'PERA/EiJ', 'PERC/EiJ', 'POHN/Deh',
                     # 'PWD/PhJ', 'PWK/PhJ', 'RBA/DnJ', 'RBB/DnJ', 'RBF/DnJ', 'SF/CamEiJ', 'SKIVE/EiJ', 'SOD1/EiJ',
                     # 'STLT', 'STRA', 'STRB', 'STUF', 'STUP', 'STUS', 'TIRANO/EiJ', 'WLA', 'WMP', 'WSB/EiJ',
                     'ZALENDE/EiJ'] if tl.is_available(s)]
    tl.unique_combos(classical, wild_derived)
    # tl.preprocess(['subspecific_origins.csv'])
    exit()
    x = TwoLocus(chrom_sizes=[20e6, 20e6])
    x.preprocess(["test2.csv"])
    x.unique_combos(['A', 'B', 'D'], ['C', 'E'])
    print type(x.sources_at_point_pair('1', 1, '1', 10000000, ['A'])['Distal']['Start']['Position'])
    # x.interlocus_dependence([chr(c) for c in xrange(ord('A'), ord('J')+1)])
    exit()

    x = TwoLocus(chrom_sizes=[20 * 10 ** 6, 20 * 10 ** 6])
    x.preprocess(["test.csv"])
    rez = x.pairwise_frequencies(["A"], include_unknown=True)

    areas = x.calculate_genomic_area(rez[0], rez[1])
    total = 0.0

    for combo in subspecies.iter_combos(include_unknown=True):
        print "\t{:15s}({:4d}):{:1.5f}".format(subspecies.to_string(combo), combo, areas[combo])
        total += areas[combo]
    print "\t{:21s}:{:1.5f}".format("Total", total)

    sys.exit(1)
    for code, combo in combos.iteritems():
        print "\n", rez[1]
        print "\t{} ({}):\n{}".format(combo, code, rez[0][code])


if __name__ == "__main__":
    main()
