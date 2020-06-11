#!/usr/bin/env python

import os
import sys
import re
import argparse

import vcf
from Bio import SeqIO

def parse_args():
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(description="""1. Checks if the consensus has the same length as the regference and print a warning if not.
2. If no Warning from first step, checks if all the variants in the vcf are found in the fasta file and return the variants not found.""")

    parser.add_argument('-c',
                        '--consensus',
                        help="consensus fasta file",
                        required=True)
    parser.add_argument('-v',
                        '--vcf',
                        help="variants vcf file",
                        required=True)

    parser.add_argument('-r',
                        '--reference',
                        help="reference fasta file",
                        required=True)

    parser.add_argument('-w',
                        '--window',
                        help="Number of bases to add on each side of the variant position to check against consensus (Default: 6)",
                        default=6,
                        type=int,
                        required=False)

    return parser.parse_args()

def check_consensus_size(fasta_file, reference_index):
    """
    Returns a warning if the len of the consensus != len of reference
    """
    ret = True
    with open(reference_index, "r") as reference:
        for line in reference:
            ref_len = int(line.split("\t")[1])
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) != ref_len:
            sys.stderr.write("""WARNING: The length of the consensus (%i) is diferent from the length of the reference (%i).
The variants can't be checked.\n""" % (len(record.seq), ref_len))
            ret = False
    return ret

def check_variants(vcf_file, consensus_file, reference_file, window):
    """
    Returns the positions in the vcf where the fasta is different
    """
    out_list = []
    insertion = 0
    deletion = 0
    indel_shift = 0
    for variant in vcf.Reader(open(vcf_file, 'r')):
        for sample in variant.samples:
            # print(sample['alt_FREQ'])
            # print(len(variant.REF), len(variant.ALT))
            # exit()
            variant.ALT = "".join(str(v) for v in variant.ALT)
            for consensus in SeqIO.parse(consensus_file, "fasta"):
                for reference in SeqIO.parse(reference_file, "fasta"):
                    # print(len(reference.seq), len(consensus.seq))
                    # for indice, nt_r in enumerate(reference.seq):
                    #     if nt_r != consensus.seq[indice]:
                    #         print(nt_r, consensus.seq[indice])
                    # print(len(variant.REF) - len(variant.ALT))

                    # deletion within variant
                    if len(variant.REF) - len(variant.ALT) > 0:
                        if consensus.seq[variant.POS - (window + 1 + deletion - insertion) : variant.POS] == reference.seq[variant.POS - (window + 1) : variant.POS] and consensus.seq[variant.POS : variant.POS + window - (deletion + insertion)] == reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))]:
                            print(
                                variant.POS,
                                variant.ALT,
                                sample['alt_FREQ'],
                                sample['alt_DP'],
                                consensus.seq[variant.POS - (window + deletion - insertion) : variant.POS] + "." + variant.ALT + "." + consensus.seq[variant.POS : variant.POS + window - (deletion + insertion)],
                                reference.seq[variant.POS - window : variant.POS] + "." + variant.REF + "." + reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))]
                                )
                            deletion += (len(variant.REF) - len(variant.ALT))

                    # insertion within variant
                    elif len(variant.REF) - len(variant.ALT) < 0:
                        print(variant.POS, variant.REF, variant.ALT, consensus.seq[variant.POS - (window + 1 + deletion - insertion) : variant.POS + window - (deletion + insertion) + abs(len(variant.REF) - len(variant.ALT))], reference.seq[variant.POS - (window + 1) : variant.POS + window], sample['alt_FREQ'], sample['alt_DP'])
                        insertion += abs(len(variant.REF) - len(variant.ALT))

                    else:
                        print(variant.POS, variant.REF, variant.ALT, consensus.seq[variant.POS - (window + 1 + deletion - insertion) : variant.POS + window - (deletion - insertion)], reference.seq[variant.POS - (window + 1) : variant.POS + window], sample['alt_FREQ'], sample['alt_DP'])

                    # if record.seq[variant.POS - 1] != "".join(str(v) for v in variant.ALT):
                    #     # print(str(variant.POS), str(variant.REF), "".join(str(v) for v in variant.ALT), str(record.seq[variant.POS - 1]))
                    #     out_list.append([str(variant.POS), str(variant.REF), str(variant.ALT), str(sample['alt_FREQ']), str(record.seq[variant.POS - 1])])
    return out_list

def main():
    """
    main
    """
    args = parse_args()

    check_variants(args.vcf, args.consensus, args.reference, args.window)
    # print("\n".join("\t".join(variant) for variant in check_variants(args.vcf, args.consensus)))
    exit()
    if check_consensus_size(args.consensus, args.index_reference):
        not_found_variants = check_variants(args.vcf, args.consensus)
        if not_found_variants:
            print("POS\tALT\tALT_FREQ\tALT_DEPTH\tCONSENSUS\tREF")
            print("\n".join("\t".join(variant) for variant in not_found_variants))


if __name__ == "__main__":
    main()
