"""
File: intervals.py
Author: Seth Greenstein
Purpose: Treat the genome as a series of intervals, some of which have a known subspecific origin.
        Over subsets of all samples, find the frequency of all pairwise combinations of origins for intervals.
        E.g. find the number of samples that are domesticus at interval A and musculus at interval B
"""

import numpy as np
import csv
import string
from collections import Counter

# integer representations of subspecies
NO_SPECIES = 0
SPECIES_TO_INT = {'None': NO_SPECIES, 'Dom': 0o001, 'Mus': 0o010, 'Cast': 0o100}
MAX_POSITION = 200000000  # 200m, bigger than the largest mouse chromosome
CHROMO_TO_INT = {str(x): x for x in xrange(1, 20)}
CHROMO_TO_INT['X'] = 21
CHROMO_TO_INT['Y'] = 22
CHROMO_TO_INT['MT'] = 23
INT_TO_CHROMO = [str(x) for x in xrange(0, 20)] + ['X', 'Y', 'MT']
CHROMO_SIZES = [200000000] * len(CHROMO_TO_INT)  # upper bound on size of each chromosome


def genome_index(chromosome, position):
    """ Converts chromosome and position to a single number
    :param chromosome: integer representation of chromosome
    :param position: position on chromosome
    :return: integer denoting chromosome and position
    """
    return ((chromosome - 1) * MAX_POSITION) + position


def index_to_chrom_pos(index):
    chromo_num = 1
    while index > CHROMO_SIZES[chromo_num]:
        index -= CHROMO_SIZES[chromo_num]
        chromo_num += 1
    return INT_TO_CHROMO[chromo_num], index


def make_source_matrix(strain_name, chromosomes):
    """ Constructs and saves the pairwise combinations of subspecific origins of intervals
    :param strain_name: string, name of the strain/sample
    :param chromosomes: dictionary with the intervals for each chromosome
    """
    num_intervals = sum([len(ints) for ints in chromosomes.itervalues()])
    intervals = np.empty(num_intervals, dtype=np.uint32)
    source_array = np.empty(num_intervals, dtype=np.uint8)
    sources = np.empty([num_intervals, num_intervals], dtype=np.uint8)
    interval_num = 0
    for chromosome, interval_list in sorted(chromosomes.iteritems(), key=lambda x: x[0]):
        for species, end in interval_list:
            intervals[interval_num] = genome_index(chromosome, end)
            source_array[interval_num] = species
            interval_num += 1
    for interval_num in xrange(num_intervals):
        sources[interval_num, :] = source_array
    for interval_num in xrange(num_intervals):
        sources[:, interval_num] = (source_array + sources[:, interval_num])
    np.save(strain_name + '_intervals.npy', intervals)
    np.save(strain_name + '_sources.npy', sources)


def rgb_to_subspecies(file_list, suffix='_ss_names'):
    """ replaces color values with subspecies names
    :param file_list: list of files to convert
    :param suffix: the converted files are saved as oldfilename_suffix.extension
    :return: list of filenames of converted files
    """
    new_filenames = []
    for file_name in file_list:
        dot_pos = file_name.find('.')
        with open(file_name) as in_file, open(file_name[:dot_pos] + suffix + file_name[dot_pos:], 'w+') as out_file:
            for row in in_file:
                row = row.replace('255 0 0', 'Mus')
                row = row.replace('0 255 0', 'Cast')
                row = row.replace('0 0 255', 'Dom')
                out_file.write(row)
    return new_filenames


def parse_csvs(file_list):
    """ parses downloaded haplotypes from the Mouse Phylogeny Viewer
    :param file_list: list of filenames
    :return: dictionary of {strain name: {chromosome number: list of (subspecies id, interval end)}}
    """
    # csv headers
    STRAIN = 'strain'
    CHROMOSOME = 'chromosome'
    START = 'start_position'
    END = 'end_position'
    SUBSPECIES = 'subspecies'
    strains = {}
    for filename in file_list:
        with open(filename) as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                chromosomes = strains.setdefault(row[STRAIN], {})
                intervals = chromosomes.setdefault(CHROMO_TO_INT[row[CHROMOSOME]], [])
                lastEnd = 0 if not intervals else intervals[-1]
                # add null interval if there is a gap between intervals with assigned subspecies
                if not int(row[START]) - 1 <= lastEnd <= int(row[START]) + 1:
                    intervals.append((NO_SPECIES, int(row[START])))
                intervals.append((SPECIES_TO_INT[row[SUBSPECIES]], int(row[END])))
    # add null interval to end of each chromosome
    for chromosomes in strains.itervalues():
        for chromosome, intervals in chromosomes.iteritems():
            intervals.append((NO_SPECIES, CHROMO_SIZES[chromosome]))
    return strains


def preprocess(file_list):
    """ creates the preprocessed matrix for each sample. Run once
    :param file_list: list of file names of csv files with subspecific origins
    :return: list of the strain names that were processed
    """
    strain_names = []
    # remove chars from strain names that can't be in filenames (e.g. /)
    valid_chars = "-_.() %s%s" % (string.ascii_letters, string.digits)
    for strain_name, chromosomes in parse_csvs(file_list).iteritems():
        strain_name = ''.join(c for c in strain_name if c in valid_chars)
        make_source_matrix(strain_name, chromosomes)
        strain_names.append(strain_name)
    return strain_names


def source_frequencies(strain_names):
    """ creates matrices of the pairwise subspecific origin frequencies
    Andrew: this is the only function you should need to edit
    :param strain_names: list of strain names to analyze (must be a subset of the output from preprocess())
    """
    interval_lists = []
    for strain_name in strain_names:
        interval_lists.append(list(np.load(strain_name + '_intervals.npy')))
    # compute elementary intervals
    elem_intervals = []
    while interval_lists:
        elem_intervals.append(min(intervals[0] for intervals in interval_lists))
        for i, intervals in enumerate(interval_lists):
            if intervals[0] == elem_intervals[-1]:
                intervals.pop(0)
            if not intervals:
                interval_lists.pop(i)
    # dictionary of matrices, one for each combination of two sources. i.e. 1 matrix has the counts for the dom-mus combo
    # elementary intervals are along the axes
    # Each element contains the number of samples that have that matrix's combo at that pair of intervals
    source_counts = {}
    for source1 in SPECIES_TO_INT.itervalues():
        if source1 != NO_SPECIES:  # don't count combinations for which one source is null
            for source2 in SPECIES_TO_INT.itervalues():
                if source2 != NO_SPECIES:
                    source_counts[source1 + source2] = np.zeros([len(elem_intervals), len(elem_intervals)],
                                                                dtype=np.uint16)
    # fill in matrices
    for strain_name in strain_names:
        intervals = np.load(strain_name + '_intervals.npy')
        sources = np.load(strain_name + '_sources.npy')
        elem_row = 0
        for row, row_end in enumerate(intervals):
            while elem_row < len(elem_intervals) and elem_intervals[elem_row] <= row_end:
                elem_col = 0
                for source, col_end in zip(sources[row], intervals):
                    while elem_col < len(elem_intervals) and elem_intervals[elem_col] <= col_end:
                        if source in source_counts:
                            source_counts[source][elem_row, elem_col] += 1
                        elem_col += 1
                elem_row += 1
    print sum([source for source in source_counts.itervalues()])
    print elem_intervals
    print source_counts


def sources_at_point_pair(chrom1, pos1, chrom2, pos2, strain_names):
    """ Prints the range of the 2D interval and the counts of subspecific combos at 2 loci in the genome
    :param chrom1: chromosome of one locus
    :param pos1: position of one locus
    :param chrom2: chromosome of another locus
    :param pos2: position of another locus
    :param strain_names: list of strain names to analyze
    """
    coords = [genome_index(CHROMO_TO_INT[str(chrom1)], pos1), genome_index(CHROMO_TO_INT[str(chrom2)], pos2)]
    mins = [0] * 2
    maxes = [np.sum(CHROMO_SIZES)] * 2
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
    print 'Proximal range: Chr %s:%d - Chr %s:%d' % tuple(x for x in index_to_chrom_pos(mins[0]) + index_to_chrom_pos(maxes[0]))
    print 'Distal range: Chr %s:%d - Chr %s:%d' % tuple(x for x in index_to_chrom_pos(mins[1]) + index_to_chrom_pos(maxes[1]))
    print source_counts
    for species1, species1_int in SPECIES_TO_INT.iteritems():
        for species2, species2_int in SPECIES_TO_INT.iteritems():
            if species1_int + species2_int in source_counts:
                print species1, '+', species2, source_counts[species1_int + species2_int]


def main():
    # you should only need to run these once
    # converted_file_names = rgb_to_subspecies(['put', 'file', 'names', 'here'])
    # file_names = converted_file_names + ['put files', 'you do not need to convert', 'here']
    # strain_names = preprocess(file_names)
    # print strain_names
    # the following you'll need to run multiple times
    # for any subset you want of strain names, run something like:
    # source_frequencies(strain_names)
    sources_at_point_pair(19, 3500000, 19, 10000000, ['CASTEiJ', 'AJ'])


if __name__ == '__main__':
    main()