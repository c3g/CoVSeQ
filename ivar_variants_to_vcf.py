#!/usr/bin/env python

"""
Convert iVar variants tsv file to vcf format. Based on nf-core-viralrecon https://github.com/drpatelh/nf-core-viralrecon/blob/0ca056e2e4b2d6a4f2c5ce08915b19e972bf5ab1/bin/ivar_variants_to_vcf.py#L57 script
"""

import os
import sys
import re
import errno
import argparse

def parse_args(args=None):
    """
    Arguments parsing
    """
    description = "Convert iVar variants tsv file to vcf format."
    epilog = """Example usage: python ivar_variants_to_vcf.py <FILE_IN> <FILE_OUT>"""

    parser = argparse.ArgumentParser(description=description, epilog=epilog)
    parser.add_argument('FILE_IN', help="Input tsv file.")
    parser.add_argument('FILE_OUT', help="Full path to output vcf file")

    return parser.parse_args(args)

def make_dir(path):
    """
    Creating directory of output file if required
    """
    if path:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

def ivar_variants_to_vcf(input_file, output_file):
    """
    Converting ivar variants tsv to vcf
    """
    filename = os.path.splitext(input_file)[0]
    header = ('##fileformat=VCFv4.2\n'
              '##source=iVar\n'
              '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">\n'
              '##FILTER=<ID=PASS,Description="Result of p-value <= 0.05">\n'
              '##FILTER=<ID=FAIL,Description="Result of p-value > 0.05">\n'
              '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
              '##FORMAT=<ID=ref_DP,Number=1,Type=Integer,Description="Depth of reference base">\n'
              '##FORMAT=<ID=ref_RV,Number=1,Type=Integer,Description="Depth of reference base on reverse reads">\n'
              '##FORMAT=<ID=ref_QUAL,Number=1,Type=Integer,Description="Mean quality of reference base">\n'
              '##FORMAT=<ID=alt_DP,Number=1,Type=Integer,Description="Depth of alternate base">\n'
              '##FORMAT=<ID=alt_RV,Number=1,Type=Integer,Description="Deapth of alternate base on reverse reads">\n'
              '##FORMAT=<ID=alt_QUAL,Number=1,Type=String,Description="Mean quality of alternate base">\n'
              '##FORMAT=<ID=alt_FREQ,Number=1,Type=String,Description="Frequency of alternate base">\n')
    header += '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + filename + '\n'

    variants_count_dict = {'SNP':0, 'INS':0, 'DEL':0}

    make_dir(os.path.dirname(output_file))

    with open(input_file, 'r') as in_file, open(output_file, 'w') as out_file:
        out_file.write(header)
        next(in_file)
        for line in in_file:
            line = re.split("\t", line)
            chrom = line[0]
            pos = line[1]
            id = '.'
            ref = line[2]
            alt = line[3]
            if alt[0] == '+':
                alt = ref + alt[1:]
                variants_count_dict['INS'] += 1
            elif alt[0] == '-':
                ref += alt[1:]
                alt = line[2]
                variants_count_dict['DEL'] += 1
            else:
                variants_count_dict['SNP'] += 1
            qual = '.'
            pass_test = line[13]
            if pass_test == 'TRUE':
                filter = 'PASS'
            else:
                filter = 'FAIL'
            info = 'DP=' + line[11]
            format = 'GT:ref_DP:ref_RV:ref_QUAL:alt_DP:alt_RV:alt_QUAL:alt_FREQ'
            sample = '1:' + line[4] + ':' + line[5] + ':' + line[6] + ':' + line[7] + ':' + line[8] + ':' + line[9] + ':' + line[10]
            line = chrom + '\t' + pos + '\t' + id + '\t' + ref + '\t' + alt + '\t' + qual + '\t' + filter + '\t' + info + '\t' + format + '\t' + sample + '\n'
            out_file.write(line)

    ## Print variant counts
    variants_count_list = [(variant, str(count)) for variant, count in sorted(variants_count_dict.items())]
    print('\t'.join(['sample'] + [x[0] for x in variants_count_list]))
    print('\t'.join([filename] + [x[1] for x in variants_count_list]))


def main(args=None):
    """
    main function
    """
    args = parse_args(args)
    ivar_variants_to_vcf(args.FILE_IN, args.FILE_OUT)


if __name__ == '__main__':
    sys.exit(main())
