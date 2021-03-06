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
import glob
import logging
# import subspecies
import pyximport

pyximport.install()
import subspeciesCython as subspecies
import pickle
from scipy import stats
from time import clock
from collections import OrderedDict, Counter

INT_TO_CHROMO = [str(integer) for integer in range(20)] + ['X', 'Y', 'MT']
# integer representations of chromosomes
CHROMO_TO_INT = {string: integer for integer, string in enumerate(INT_TO_CHROMO)}

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

    def chrom_and_pos(self, index, index2=None):
        """ Converts genome position to chromosome and position
        :param index: integer denoting chromosome and position
        :param index2: second integer denoting chromosome and position (optional)
        :return: string representation of chromosome, position on chromosome
        """
        if index2 is None:
            return self._chrom_and_pos(index)
        else:
            start_chrom, start_pos = self._chrom_and_pos(index)
            end_chrom, end_pos = self._chrom_and_pos(index2)
            if start_chrom != end_chrom:
                start_pos = 0  # was end of previous chromosome, so change to beginning of current
            return end_chrom, int(start_pos), int(end_pos)  # convert from numpy int type

    def _chrom_and_pos(self, index):
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
            # TODO: remove stuff about "bad" (residual heterozygosity)
            bad = False
            for intervals in chromosomes.itervalues():
                for ss, _ in intervals:
                    if ss == -999:
                        bad = True
            if not bad:
                print 'good', strain_name
                self.sample_dict[strain_name] = self.intervals_and_sources(chromosomes)
            else:
                print 'bad', strain_name
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
                        intervals.append((subspecies.UNKNOWN, int(row[START])))
                    intervals.append((subspec_int(row), int(row[END])))
        # add null interval to end of each chromosome
        for chromosomes in strains.itervalues():
            for chromosome, intervals in chromosomes.iteritems():
                if intervals[-1] < self.sizes[chromosome - 1]:
                    intervals.append((subspecies.UNKNOWN, self.sizes[chromosome - 1]))
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

    # @profile
    def build_pairwise_matrix(self, strain_names, elem_intervals):
        # 3d matrix. First index is combo, remaining 2d matrices are counts for pairwise intervals
        source_counts = np.zeros([(subspecies.NUM_SUBSPECIES + 1) ** 2, len(elem_intervals), len(elem_intervals)],
                                 dtype=np.int16)
        for strain_name in strain_names:
            intervals, sources = self.sample_dict[strain_name]
            # map this strain's intervals onto the elementary intervals
            breaks = np.insert(np.searchsorted(elem_intervals, intervals), 0, -1)
            for row in xrange(len(intervals)):
                for col in xrange(row, len(intervals)):  # only upper triangle
                    source = subspecies.combine(sources[row], sources[col])
                    source_ordinate = subspecies.to_ordinal(source)
                    source_counts[source_ordinate, breaks[row] + 1:breaks[row + 1] + 1,
                    breaks[col] + 1:breaks[col + 1] + 1] += 1
        return source_counts

    def pairwise_frequencies(self, strain_names):
        """ For every locus pair and every label pair, count the number of strains which have those
        labels at those pairs of loci.
        :param strain_names: list of strain names to analyze (must be a subset of the output from preprocess())
        """
        output = [[[], [], [], []] for _ in xrange(subspecies.NUM_SUBSPECIES**2)]
        for strain_name in strain_names:
            intervals, sources = self.sample_dict[strain_name]
            for i in xrange(len(intervals)):
                # only upper triangle is meaningful
                if subspecies.is_known(sources[i]):
                    for j in xrange(i, len(intervals)):
                        if subspecies.is_known(sources[j]):
                            combo_output = output[subspecies.to_ordinal(subspecies.combine(sources[i], sources[j]))]
                            combo_output[0].append(intervals[i-1])
                            combo_output[1].append(intervals[i])
                            combo_output[2].append(intervals[j-1])
                            combo_output[3].append(intervals[j])
        return output, [subspecies.to_color(i, True) for i in xrange(subspecies.NUM_SUBSPECIES**2)]

    def absent_regions(self, strain_names):
        """ finds regions in which no samples have a certain combo
        :param strain_names: list of strain names to analyze (must be a subset of the output from preprocess())
        """
        elem_intervals = self.make_elementary_intervals(
            [self.sample_dict[sn][0] for sn in strain_names])
        background = self.build_pairwise_matrix(strain_names, elem_intervals)
        output = [[[], [], [], []] for _ in xrange(subspecies.NUM_SUBSPECIES**2)]
        for combo in xrange(subspecies.NUM_SUBSPECIES**2):
            for i in xrange(len(elem_intervals)):
                for j in xrange(i, len(elem_intervals)):
                    if not background[combo, i, j]:
                        output[combo][0].append(elem_intervals[i-1])
                        output[combo][1].append(elem_intervals[i])
                        output[combo][2].append(elem_intervals[j-1])
                        output[combo][3].append(elem_intervals[j])
        return output

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
        for combo, vals in enumerate(counts):
            factor = 1
            areas_masked.update({str(subspecies.to_string(combo, True)): np.sum((vals > 0) * areas * factor) / denom})
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
        output = {}
        samples = [[[] for _ in subspecies.iter_subspecies(True)] for _ in subspecies.iter_subspecies(True)]
        key = [subspecies.to_string(s) for s in subspecies.iter_subspecies(True)]
        for strain_name in strain_names:
            intervals = self.sample_dict[strain_name][0]
            sources = self.sample_dict[strain_name][1]
            # find interval containing each location
            i = 0
            interval_indices = [None, None]
            for loc_num in xrange(2):
                while intervals[i] < coords[loc_num]:
                    i += 1
                if i > 0:
                    mins[loc_num] = max(mins[loc_num], intervals[i - 1])
                maxes[loc_num] = min(maxes[loc_num], intervals[i])
                interval_indices[loc_num] = i
            samples[subspecies.to_ordinal(sources[interval_indices[0]])][
                subspecies.to_ordinal(sources[interval_indices[1]])].append(strain_name)
        output['Key'] = key
        output['Samples'] = samples
        output['Intervals'] = [
            self.chrom_and_pos(mins[0], maxes[0]),
            self.chrom_and_pos(mins[1], maxes[1])
        ]
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
        output = []
        for i in xrange(len(intervals)):
            # only upper triangle is meaningful
            for j in xrange(i + 1, len(intervals)):
                nonzero_expectations = np.where(combo_expectations[i, j])
                chi_sq, p_value = stats.chisquare(
                    combo_counts[i, j][nonzero_expectations], combo_expectations[i, j][nonzero_expectations])
                output.append([
                        chi_sq,
                        p_value,
                        # proximal interval
                        intervals[i-1],
                        intervals[i],
                        # distal interval
                        intervals[j],
                        intervals[j-1]
                        ])
        return output

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

    def not_in_background(self, background_strains, foreground_strains):
        """ finds combinations at interval pairs that are present in 1+ fg strains but is absent from the background
        :param background_strains: list of strain names
        :param foreground_strains: list of strain names
        :return: json object containing interval pairs
        """
        output = [[[], [], [], [], []] for _ in xrange(subspecies.NUM_SUBSPECIES**2)]
        for strain in foreground_strains:
            elem_intervals = self.make_elementary_intervals(
                [self.sample_dict[sn][0] for sn in background_strains + [strain]])
            background_absent = np.logical_not(self.build_pairwise_matrix(background_strains, elem_intervals))
            foreground = self.build_pairwise_matrix([strain], elem_intervals)
            uniquities = np.logical_and(foreground, background_absent)
            for combo in xrange(subspecies.NUM_SUBSPECIES**2):
                combo_uniquities = np.where(uniquities[combo])
                for i, j in zip(combo_uniquities[0], combo_uniquities[1]):
                    output[combo][0].append(elem_intervals[i-1])
                    output[combo][1].append(elem_intervals[i])
                    output[combo][2].append(elem_intervals[j-1])
                    output[combo][3].append(elem_intervals[j])
                    output[combo][4].append(strain)
        return output, [subspecies.to_color(combo, ordinal=True) for combo in xrange(subspecies.NUM_SUBSPECIES**2)]

    # @profile
    def unique_combos(self, background_strains, foreground_strains):
        """ finds combinations at interval pairs that is absent from the background but shared by all foreground samples
        :param background_strains: list of strain names
        :param foreground_strains: list of strain names
        :return: json object containing interval pairs
        """
        elem_intervals = self.make_elementary_intervals(
            [self.sample_dict[sn][0] for sn in background_strains + foreground_strains])
        background = self.build_pairwise_matrix(background_strains, elem_intervals)
        foreground = self.build_pairwise_matrix(foreground_strains, elem_intervals)
        output = []
        uniquities = np.logical_and(foreground == len(foreground_strains), np.logical_not(background))
        for combo in xrange(subspecies.NUM_SUBSPECIES**2):
            combo_uniquities = np.where(uniquities[combo])
            combo_color = subspecies.to_color(combo, ordinal=True)
            for i, j in zip(combo_uniquities[0], combo_uniquities[1]):
                output.append([
                    # proximal interval start, end
                    elem_intervals[i - 1],
                    elem_intervals[i],
                    # distal interval start, end
                    elem_intervals[j - 1],
                    elem_intervals[j],
                    combo_color
                ])
        return output

    def contingency_table(self, dead_strains, live_strains, output_file):
        elem_intervals = self.make_elementary_intervals(
            [self.sample_dict[sn][0] for sn in dead_strains + live_strains]
        )
        num_dead = len(dead_strains)
        num_live = len(live_strains)
        dead_observed = self.build_pairwise_matrix(dead_strains, elem_intervals)
        live_observed = self.build_pairwise_matrix(live_strains, elem_intervals)
        with open(output_file, 'w+') as fp:
            writer = csv.writer(fp)
            writer.writerow(['Proximal chromosome', 'Proximal start', 'Proximal end',
                             'Distal chromosome', 'Distal start', 'Distal end',
                             'Proximal origin', 'Distal origin', 'chi squared', 'p-value'])
            elem_intervals.insert(0, 0)
            for combo in xrange(subspecies.NUM_SUBSPECIES**2):
                for i in xrange(len(elem_intervals)-1):
                    for j in xrange(i+1, len(elem_intervals)-1):
                        if dead_observed[combo, i, j] and live_observed[combo, i, j]:
                            contingency = np.array([[dead_observed[combo, i, j], live_observed[combo, i, j]],
                                                    [num_dead-dead_observed[combo, i, j],
                                                     num_live-live_observed[combo, i, j]]])
                            chi_squared, p, _, _ = stats.chi2_contingency(contingency)
                            proximal_pos = self.chrom_and_pos(elem_intervals[i], elem_intervals[i+1])
                            distal_pos = self.chrom_and_pos(elem_intervals[j], elem_intervals[j+1])
                            writer.writerow(proximal_pos + distal_pos +
                                            (subspecies.proximal(combo), subspecies.distal(combo), chi_squared, p))


def main():
    """ Run some tests with a dummy file, overriding chromosome lengths locally for sake of testing.
    """
    tl = TwoLocus(in_path='/csbiodata/public/www.csbio.unc.edu/htdocs/sgreens/pairwise_origins/')
    # tl = TwoLocus()
    # tl.preprocess(glob.glob('OR_ss_origins/*.hap'))
    print len(tl.list_available_strains())
    exit()
    # print len(tl.list_available_strains())
    # tl.preprocess(['cc_origins.csv'])
    # tl.preprocess(['ccv_origins.csv'])
    classical = [s for s in
                 ["129P1/ReJ",  # "129P3/J", "129S1SvlmJ", "129S6", "129T2/SvEmsJ", "129X1/SvJ", "A/J", "A/WySnJ",
                  "AEJ/GnLeJ", "AEJ/GnRk", "AKR/J", "ALR/LtJ", "ALS/LtJ", "BALB/cByJ", "BALB/cJ", "BDP/J", "BPH/2J",
                  # "BPL/1J", "BPN/3J", "BTBR T<+>tf/J", "BUB/BnJ", "BXSB/MpJ", "C3H/HeJ", "C3HeB/FeJ", "C57BL/10J",
                  # "C57BL/10ScNJ", "C57BL/10SAAAJ", "C57BL/6CR", "C57BL/6J", "C57BL/6NCI", "C57BL/6Tc", "C57BLKS/J",
                  # "C57BR/cdJ", "C57L/J", "C58/J", "CBA/CaJ", "CBA/J", "CE/J", "CHMU/LeJ", "DBA/1J", "DBA/1LacJ",
                  # "DBA/2DeJ", "DBA/2HaSmnJ", "DBA/2J", "DDK/Pas", "DDY/JclSidSeyFrkJ", "DLS/LeJ", "EL/SuzSeyFrkJ",
                  # "FVB/NJ", "HPG/BmJ", "I/LnJ", "IBWSP2", "IBWSR2", "ICOLD2", "IHOT1", "IHOT2", "ILS", "ISS", "JE/LeJ",
                  # "KK/HlJ", "LG/J", "LP/J", "LT/SvEiJ", "MRL/MpJ", "NOD/ShiLtJ", "NON/ShiLtJ", "NONcNZO10/LtJ",
                  # "NONcNZO5/LtJ", "NOR/LtJ", "NU/J", "NZB/BlNJ", "NZL/LtJ", "NZM2410/J", "NZO/HlLtJ", "NZW/LacJ", "P/J",
                  # "PL/J", "PN/nBSwUmabJ", "RF/J", "RHJ/LeJ", "RIIIS/J", "RSV/LeJ", "SB/LeJ", "SEA/GnJ", "SEC/1GnLeJ",
                  # "SEC/1ReJ", "SH1/LeJ", "SI/Col Tyrp1 Dnahc11/J", "SJL/Bm", "SJL/J", "SM/J", "SSL/LeJ", "ST/bJ",
                  "STX/Le", ]  # "SWR/J", "TALLYHO/JngJ", "TKDU/DnJ", "TSJ/LeJ", "YBR/EiJ", "ZRDCT Rax<+>ChUmdJ"]
                 if tl.is_available(s)]
    wild_derived = [s for s in
                    ['22MO',
                     # 'BIK/g', 'BULS', 'BUSNA', 'BZO', 'CALB/RkJ', 'CASA/RkJ', 'CAST/EiJ', 'CIM', 'CKN', 'CKS',
                     'CZECHI/EiJ', 'CZECHII/EiJ', 'DCA', 'DCP', 'DDO', 'DEB', 'DGA', 'DIK', 'DJO', 'DKN', 'DMZ', 'DOT',
                     # 'IS/CamRkJ', 'JF1/Ms', 'LEWES/EiJ', 'MBK', 'MBS', 'MCZ', 'MDG', 'MDGI', 'MDH', 'MGA', 'MH',
                     # 'MOLD/RkJ', 'MOLF/EiJ', 'MOLG/DnJ', 'MOR/RkJ', 'MPB', 'MSM/Ms', 'PERA/EiJ', 'PERC/EiJ', 'POHN/Deh',
                     # 'PWD/PhJ', 'PWK/PhJ', 'RBA/DnJ', 'RBB/DnJ', 'RBF/DnJ', 'SF/CamEiJ', 'SKIVE/EiJ', 'SOD1/EiJ',
                     # 'STLT', 'STRA', 'STRB', 'STUF', 'STUP', 'STUS', 'TIRANO/EiJ', 'WLA', 'WMP', 'WSB/EiJ',
                     'ZALENDE/EiJ'] if tl.is_available(s)]
    tl.contingency_table(classical, wild_derived, '/csbiohome01/sgreens/Projects/intervals/contingency.csv')
    exit()
    x = TwoLocus(chrom_sizes=[20e6, 20e6])
    x.preprocess(["test2.csv"])
    x.unique_combos(['A', 'B', 'D'], ['C', 'E'])
    x.sources_at_point_pair('1', 1, '1', 10000000, ['A'])
    # x.interlocus_dependence([chr(c) for c in xrange(ord('A'), ord('J')+1)])
    # exit()

    x = TwoLocus(chrom_sizes=[20 * 10 ** 6, 20 * 10 ** 6])
    x.preprocess(["test.csv"])
    rez = x.pairwise_frequencies(["A"])

    areas = x.calculate_genomic_area(rez[0], rez[1])
    total = 0.0

    for combo in subspecies.iter_combos():
        print "\t{:15s}({:4d}):{:1.5f}".format(subspecies.to_string(combo), combo,
                                               areas[str(subspecies.to_string(combo))])
        total += areas[str(subspecies.to_string(combo))]
    print "\t{:21s}:{:1.5f}".format("Total", total)

    sys.exit(1)
    # for code, combo in combos.iteritems():
    #     print "\n", rez[1]
    #     print "\t{} ({}):\n{}".format(combo, code, rez[0][code])


if __name__ == "__main__":
    main()
