#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 5/25/24
'''

import gzip
import os
import sys

AUTOSOMESXY = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                  "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", 'chrY']

def do_parse(in_vcf, sample):
    headers = []
    sv_records = []

    with gzip.open(in_vcf, 'rt') as f:
        for line in f:
            if '##' in line or '#CHROM' in line:
                if '#CHROM' not in line:
                    headers.append(line.strip())
                continue

            entries = line.strip().split("\t")
            chrom, start, sv_id = entries[0], int(entries[1]), entries[2]

            if chrom not in AUTOSOMESXY:
                continue

            ## Reference is a Gap or no alternative sequence
            if 'N' in entries[3]:
                continue

            info_dict = parse_vcf_info_column(entries[7])

            svtype = info_dict['SVTYPE']

            if svtype == 'SNV':
                continue

            svlen = abs(int(info_dict['SVLEN']))

            if svlen < 50 or svlen > 100000:
                continue

            end = start + svlen

            info_strs = ''
            for key, val in info_dict.items():
                info_strs += f'{key}={val};'
            info_strs += f'END={end}'

            str1 = '\t'.join(entries[0:7])
            str2 = '\t'.join(entries[8:])
            new_vcf_str = f'{str1}\t{info_strs}\t{str2}'

            sv_records.append(new_vcf_str)

    for header_line in headers:
        if 'ID=SVLEN' in header_line:
            continue
        print(header_line, file=sys.stdout)

    print('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Variant length">', file=sys.stdout)
    print('##INFO=<ID=END,Number=1,Type=Integer,Description="Variant end">', file=sys.stdout)
    print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}', file=sys.stdout)

    for sv in sv_records:
        print(sv, file=sys.stdout)

def parse_vcf_info_column(info_str):
    info_tokens = info_str.split(";")
    info_dict = {}

    for token in info_tokens:
        if "=" not in token:
            continue
        info_dict[token.split('=')[0]] = token.split('=')[1]

    return info_dict

if __name__ == '__main__':

    do_parse(sys.argv[1], sys.argv[2])
