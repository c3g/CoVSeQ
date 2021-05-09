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
    description = """1. Checks if the consensus has the same length as the regference and print a warning if not.
2. If no Warning from first step, checks if all the variants in the vcf are found in the fasta file and return the variants not found."""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('-c',
                        '--consensus',
                        help="Consensus fasta file",
                        required=True)

    parser.add_argument('-v',
                        '--vcf',
                        help="Variants vcf file",
                        required=True)

    parser.add_argument('-r',
                        '--reference',
                        help="Reference fasta file",
                        required=True)

    parser.add_argument('-w',
                        '--window',
                        help="Number of bases to add on each side of the variant position to check against consensus (Default: 6)",
                        default=6,
                        type=int,
                        required=False)

    parser.add_argument('-o',
                        '--output',
                        help="tsv file with variants found within vcf file (Default: stdout)",
                        type=argparse.FileType('w'),
                        default='-')

    return parser.parse_args()

def check_consensus_size(fasta_file, reference):
    """
    Returns a warning if the len of the consensus != len of reference
    """
    reference_index = reference + ".fai"
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

def store_current_out(variant, consensus, reference, sample, window, deletion, insertion, variant_type, variant_match, context_match, out_list):
    ref_interval = str(reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + variant.REF + "." + reference.seq[variant.POS : variant.POS + window])
    cons_interval = str(consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] + "." + consensus.seq[variant.POS - 1 - deletion + insertion: variant.POS - deletion + insertion] + "." + consensus.seq[variant.POS - deletion + insertion : variant.POS + window - deletion + insertion])
    if variant_type == "deletion":
        ref_interval = str(reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + reference.seq[variant.POS - 1 : variant.POS + (len(variant.REF) - len(variant.ALT))] + "." + reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))])
    elif variant_type == "insertion":
        cons_interval = str(consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] + "." + consensus.seq[variant.POS - 1 - deletion + insertion : variant.POS - deletion + insertion + (len(variant.ALT) - len(variant.REF))] + "." + consensus.seq[variant.POS - deletion + insertion + (len(variant.ALT) - len(variant.REF)) : variant.POS + window - deletion + insertion + (len(variant.ALT) - len(variant.REF))])
    out_list.append([
        str(variant.POS),
        variant.REF,
        variant.ALT,
        str(sample['alt_FREQ']),
        str(sample['alt_DP']),
        ref_interval,
        cons_interval,
        str(variant_match),
        str(context_match)
        ])

def check_variants(vcf_file, consensus_file, reference_file, window):
    """
    Returns the positions in the vcf where the fasta is different
    """
    out_list = []
    insertion = 0
    deletion = 0
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

                    variant_match = True

                    # deletion within variant
                    if len(variant.REF) - len(variant.ALT) > 0:
                        context_match = consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] == reference.seq[variant.POS - window - 1 : variant.POS - 1] and consensus.seq[variant.POS - deletion : variant.POS + window - deletion + insertion] == reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))]
                        variant_match = variant.ALT == consensus.seq[variant.POS - 1 - deletion + insertion : variant.POS - deletion + insertion] and consensus.seq[variant.POS - deletion : variant.POS + window - deletion + insertion] == reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))]
                        # if consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] == reference.seq[variant.POS - window - 1 : variant.POS - 1] and consensus.seq[variant.POS - deletion : variant.POS + window - deletion + insertion] == reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))]:
                        #     context_match = True
                        # else:
                        #     context_match = False
                            # out_list.append([
                            #     str(variant.POS),
                            #     variant.REF,
                            #     variant.ALT,
                            #     str(sample['alt_FREQ']),
                            #     str(sample['alt_DP']),
                            #     str(reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + variant.REF + "." + reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))]),
                            #     str(consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] + "." + variant.ALT + "." + consensus.seq[variant.POS - deletion + insertion : variant.POS + window - deletion + insertion])
                            #     ])
                            # print(
                            #     variant.POS,
                            #     variant.REF,
                            #     variant.ALT,
                            #     sample['alt_FREQ'],
                            #     sample['alt_DP'],
                            #     reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + variant.REF + "." + reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))],
                            #     consensus.seq[variant.POS - (window + deletion - insertion) - 1 : variant.POS - 1] + "." + variant.ALT + "." + consensus.seq[variant.POS : variant.POS + window - (deletion + insertion)]
                            #     )
                        store_current_out(variant, consensus, reference, sample, window, deletion, insertion, "deletion", variant_match, context_match, out_list)
                        deletion += (len(variant.REF) - len(variant.ALT))
                            # exit()

                    # insertion within variant
                    elif len(variant.REF) - len(variant.ALT) < 0:
                        context_match = consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] == reference.seq[variant.POS - window - 1 : variant.POS - 1] and consensus.seq[variant.POS - deletion + insertion + (len(variant.ALT) - len(variant.REF)) : variant.POS + window - deletion + insertion + (len(variant.ALT) - len(variant.REF))] == reference.seq[variant.POS : variant.POS + window]
                        variant_match = variant.ALT == consensus.seq[variant.POS - 1 - deletion + insertion : variant.POS - deletion + insertion + (len(variant.ALT) - len(variant.REF))]
                        # if consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] == reference.seq[variant.POS - window - 1 : variant.POS - 1] and consensus.seq[variant.POS - deletion + insertion + (len(variant.ALT) - len(variant.REF)) : variant.POS + window - deletion + insertion + (len(variant.ALT) - len(variant.REF))] == reference.seq[variant.POS : variant.POS + window]:
                        #     context_match = True
                        # else:
                        #     context_match = False
                            # out_list.append([
                            #     str(variant.POS),
                            #     variant.REF,
                            #     variant.ALT,
                            #     str(sample['alt_FREQ']),
                            #     str(sample['alt_DP']),
                            #     str(reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + variant.REF + "." + reference.seq[variant.POS : variant.POS + window]),
                            #     str(consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] + "." + variant.ALT + "." + consensus.seq[variant.POS - deletion + insertion + (len(variant.ALT) - len(variant.REF)) : variant.POS + window - deletion + insertion + (len(variant.ALT) - len(variant.REF))])
                            #     ])
                            # print(
                            #     variant.POS,
                            #     variant.REF,
                            #     variant.ALT,
                            #     sample['alt_FREQ'],
                            #     sample['alt_DP'],
                            #     reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + variant.REF + "." + reference.seq[variant.POS : variant.POS + window],
                            #     consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] + "." + variant.ALT + "." + consensus.seq[variant.POS - deletion + insertion + (len(variant.ALT) - len(variant.REF)) : variant.POS + window - deletion + insertion + (len(variant.ALT) - len(variant.REF))]
                            #     )
                        store_current_out(variant, consensus, reference, sample, window, deletion, insertion, "insertion", variant_match, context_match, out_list)
                        insertion += (len(variant.ALT) - len(variant.REF))
                        # print(variant.POS, variant.REF, variant.ALT, consensus.seq[variant.POS - (window + 1 + deletion - insertion) : variant.POS + window - (deletion + insertion) + abs(len(variant.REF) - len(variant.ALT))], reference.seq[variant.POS - (window + 1) : variant.POS + window], sample['alt_FREQ'], sample['alt_DP'])

                    else:
                        context_match = consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] == reference.seq[variant.POS - window - 1 : variant.POS - 1] and consensus.seq[variant.POS - deletion + insertion : variant.POS + window - deletion + insertion] == reference.seq[variant.POS : variant.POS + window]
                        variant_match = variant.ALT == consensus.seq[variant.POS - 1 - deletion + insertion : variant.POS - deletion + insertion]
                        # if consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] == reference.seq[variant.POS - window - 1 : variant.POS - 1] and consensus.seq[variant.POS - deletion + insertion : variant.POS + window - deletion + insertion] == reference.seq[variant.POS : variant.POS + window]:
                        #     context_match = True
                        # else:
                        #     context_match = False
                        store_current_out(variant, consensus, reference, sample, window, deletion, insertion, None, variant_match, context_match, out_list)
                            # out_list.append([
                            #     str(variant.POS),
                            #     variant.REF,
                            #     variant.ALT,
                            #     str(sample['alt_FREQ']),
                            #     str(sample['alt_DP']),
                            #     str(reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + variant.REF + "." + reference.seq[variant.POS : variant.POS + window]),
                            #     str(consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] + "." + variant.ALT + "." + consensus.seq[variant.POS - deletion + insertion : variant.POS + window - deletion + insertion])
                            #     ])
                        # print(
                        #     variant.POS,
                        #     variant.REF,
                        #     variant.ALT,
                        #     sample['alt_FREQ'],
                        #     sample['alt_DP'],
                        #     reference.seq[variant.POS - window - 1 : variant.POS - 1] + "." + variant.REF + "." + reference.seq[variant.POS + (len(variant.REF) - len(variant.ALT)) : variant.POS + window + (len(variant.REF) - len(variant.ALT))],
                        #     consensus.seq[variant.POS - window - deletion + insertion - 1 : variant.POS - 1 - deletion + insertion] + "." + variant.ALT + "." + consensus.seq[variant.POS - deletion + insertion : variant.POS + window - deletion + insertion]
                        #     )
                        # exit()
                        # print(variant.POS, variant.REF, variant.ALT, consensus.seq[variant.POS - (window + 1 + deletion - insertion) : variant.POS + window - (deletion - insertion)], reference.seq[variant.POS - (window + 1) : variant.POS + window], sample['alt_FREQ'], sample['alt_DP'])

                    # if record.seq[variant.POS - 1] != "".join(str(v) for v in variant.ALT):
                    #     # print(str(variant.POS), str(variant.REF), "".join(str(v) for v in variant.ALT), str(record.seq[variant.POS - 1]))
                    #     out_list.append([str(variant.POS), str(variant.REF), str(variant.ALT), str(sample['alt_FREQ']), str(record.seq[variant.POS - 1])])
    return out_list

def main():
    """
    main
    """
    args = parse_args()

    check_consensus_size(args.consensus, args.reference)

    variants_list = check_variants(args.vcf, args.consensus, args.reference, args.window)
    print("POS\tREF\tALT\tALT_FREQ\tALT_DEPTH\tREF+-{window_size}\tCONSENSUS+-{window_size}\tVARIANT_MATCH_CONSENSUS\tCONTEXT_MATCH+-{window_size}".format(window_size=args.window))
    print("\n".join("\t".join(variant) for variant in variants_list))

    # TODO: add a flag for context match and variant match

    # print("\n".join("\t".join(variant) for variant in check_variants(args.vcf, args.consensus)))
    # exit()
    # if check_consensus_size(args.consensus, args.index_reference):
    #     not_found_variants = check_variants(args.vcf, args.consensus)
    #     if not_found_variants:
    #         print("POS\tALT\tALT_FREQ\tALT_DEPTH\tCONSENSUS\tREF")
    #         print("\n".join("\t".join(variant) for variant in not_found_variants))


if __name__ == "__main__":
    main()
