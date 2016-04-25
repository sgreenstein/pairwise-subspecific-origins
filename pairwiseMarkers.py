import scipy.stats
import numpy as np
import csv


def main():
    genotype_dump = 'CC_3.csv'
    outfile = 'outfile.csv'
    with open(genotype_dump) as p_fp, open(genotype_dump) as d_fp, open(outfile, 'w+') as output_fp:
        p_reader = csv.reader(p_fp)
        d_reader = csv.reader(d_fp)
        writer = csv.writer(output_fp)
        writer.writerow(['Proximal Marker', 'Proximal ChrB37', 'Proximal PosB37',
                         'Distal Marker', 'Distal ChrB37', 'Distal PosB37',
                         'Chi-squared', 'p', 'corrected p-value'])
        d_reader.next()  # skip header
        table = np.empty([2, 2])
        for d_row in d_reader:
            filter_distal = False
            p_fp.seek(0)
            p_reader.next()  # skip header
            d_present_alleles = []
            for p_row in p_reader:
                filter_proximal = False
                if p_row[1] >= d_row[1]:
                    break
                table.fill(0)
                p_present_alleles = []
                for p_allele, d_allele in zip(p_row[4::2], d_row[4::2]):
                    if p_allele not in p_present_alleles:
                        p_present_alleles.append(p_allele)
                        if p_allele == 'H' or len(p_present_alleles) > 2:
                            filter_proximal = True
                            break
                    if d_allele not in d_present_alleles:
                        d_present_alleles.append(d_allele)
                        if d_allele == 'H' or len(d_present_alleles) > 2:
                            filter_distal = True
                            break
                    table[p_present_alleles.index(p_allele), d_present_alleles.index(d_allele)] += 1
                if filter_distal:
                    break
                elif not filter_proximal:
                    try:
                        writer.writerow(tuple(p_row[:3]) + tuple(d_row[:3]) + scipy.stats.chi2_contingency(table)[:2])
                    except ValueError:
                        pass  # had zero as expected


if __name__ == '__main__':
    main()
