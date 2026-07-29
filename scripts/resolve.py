#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 5/6/24
'''

import sys
import pysam
import truvari

REF_CLEAN = True  # Set to false if you're working with the right reference
MAX_SV = 100_000_000  # Filter things smaller than this

RC = str.maketrans("ATCG", "TAGC")


def parse_vcf_info_column(info_str):
    info_tokens = info_str.split(";")
    info_dict = {}

    for token in info_tokens:
        if "=" not in token:
            continue
        info_dict[token.split('=')[0]] = token.split('=')[1]

    return info_dict
def do_rc(s):
    """
    Reverse complement a sequence
    """
    return s.translate(RC)[::-1]


def resolve(entry, ref):
    """
    """

    if entry.start > ref.get_reference_length(entry.chrom):
        return None
    if entry.alts[0] in ['<CNV>', '<INS>']:
        return None

    seq_prime = ref.fetch(entry.chrom, entry.start, entry.stop)
    # FC> Replacing non-standard DNA characters with an N
    seq = ""
    for i in range(len(seq_prime)):
        c = seq_prime[i].upper()
        if c != 'A' and c != 'C' and c != 'G' and c != 'T':
            seq = seq + 'N'
        else:
            seq = seq + seq_prime[i]

    if entry.alts[0] == '<DEL>':
        entry.ref = seq
        entry.alts = [seq[0]]
    elif entry.alts[0] == '<INV>':
        entry.ref = seq
        entry.alts = [do_rc(seq)]
    elif entry.alts[0] == '<DUP>':
        entry.info['SVTYPE'] = 'INS'
        entry.ref = seq[0]
        entry.alts = [seq]
        entry.stop = entry.start + 1
    entry.qual = 1

    return entry

if __name__ == '__main__':
    default_quals = {"pbsv": 6,
                     "sniffles": 5,
                     "pav": 7, 'svimasm': 3, 'cutesv': 4, 'delly': 2, 'dip': 8}

    VALIDCHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                   "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",
                   'chrY']

    caller = sys.argv[2]
    d_qual = 1 if caller not in default_quals else default_quals[caller]


    if caller == 'svimasm':
        headers = []
        sv_records = []
        for line in sys.stdin:
            if line.startswith('#'):
                headers.append(line.strip())
                # print(line.strip(), file=sys.stdout)
                continue

            entries = line.strip().split('\t')
            info_dict = parse_vcf_info_column(entries[7])

            if entries[0] not in VALIDCHROMS:
                continue
            if info_dict['SVTYPE'] != 'BND':
                end = int(info_dict['END'])
                svlen = abs(int(info_dict['SVLEN'])) if 'SVLEN' in info_dict else end - int(entries[1]) + 1
                info_dict['SVIMASM_ID'] = entries[2]
                if 1000000 > svlen >= 50:

                    new_info_str = ''
                    for key, val in info_dict.items():
                        new_info_str += f'{key}={val};'
                    rec = f'{entries[0]}\t{entries[1]}\t{entries[2]}\t{entries[3]}\t{entries[4]}\t{default_quals[caller]}\t{entries[6]}\t{new_info_str[:-1]}\t{entries[8]}\t{entries[9]}'
                    sv_records.append(rec)
                    # print(out_line, file=sys.stdout)


        for header in headers[:-2]:
            print(header, file=sys.stdout)
        print('##INFO=<ID=SVIMASM_ID,Number=1,Type=String,Description="SV ID">', file=sys.stdout)
        print(f'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample', file=sys.stdout)

        for rec in sv_records:
            print(rec, file=sys.stdout)


    elif caller == 'pav':

        vcf = pysam.VariantFile(sys.stdin)
        ref = pysam.FastaFile(sys.argv[3])
        n_header = vcf.header.copy()

        if REF_CLEAN:
            for ctg in vcf.header.contigs.keys():
                if ctg not in ref.references:
                    n_header.contigs[ctg].remove_header()

        out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
        seen = set()
        for entry in vcf:

            key = truvari.entry_to_hash(entry)
            if key in seen:
                continue
            seen.add(key)

            if REF_CLEAN and entry.chrom not in VALIDCHROMS:
                continue

            if truvari.entry_size(entry) >= MAX_SV or truvari.entry_size(entry) < 50:
                continue

            entry.qual = d_qual
            if entry.alts[0].startswith("<"):
                entry = resolve(entry, ref)

            if entry is None:
                continue

            if entry.info['SVTYPE'] == 'DEL' and (set(entry.alts[0]) == {'N'} or entry.alts[0] == '.'):
                entry.alts = entry.ref[0]
                if len(entry.ref) < truvari.entry_size(entry):
                    entry = resolve(entry, ref)

            if entry.info['SVTYPE'] != 'INV':
                entry.info['SVLEN'] = abs(len(entry.ref) - len(entry.alts[0]))
            else:
                entry.info['SVLEN'] = len(entry.ref)
            # No more blank genotypes
            n_gt = tuple([_ if _ is not None else 0 for _ in entry.samples[0]['GT']])
            # Preserve phasing informatino
            is_phased = entry.samples[0].phased
            entry.samples[0]['GT'] = n_gt
            entry.samples[0].phased = is_phased

            entry.translate(n_header)
            try:
                out.write(entry)
            except Exception:
                sys.stderr.write(f"{entry}\n{type(entry)}\n")

    else:
        vcf = pysam.VariantFile(sys.stdin)
        ref = pysam.FastaFile(sys.argv[3])
        n_header = vcf.header.copy()
        filt = pysam.VariantFile('excluded_SVs.vcf', 'w', header=n_header)
        if REF_CLEAN:
            for ctg in vcf.header.contigs.keys():
                if ctg not in ref.references:
                    n_header.contigs[ctg].remove_header()

        out = pysam.VariantFile("/dev/stdout", 'w', header=n_header)
        seen = set()
        for entry in vcf:

            if entry.info['SVTYPE'] == 'BND':
                continue

            if entry.filter.keys()[0] != 'PASS':
                continue

            if entry.alts[0] == '.':
                continue

            if entry.start < 10000:
                continue

            # if 'IMPRECISE' in entry.info and entry.info['IMPRECISE']:
            #     filt.write(entry)
            #     continue

            key = truvari.entry_to_hash(entry)
            if key in seen:
                continue
            seen.add(key)

            if REF_CLEAN and entry.chrom not in VALIDCHROMS:
                continue

            if truvari.entry_size(entry) >= MAX_SV:
                continue

            entry.qual = d_qual
            if entry.alts[0].startswith("<"):
                entry = resolve(entry, ref)

            if entry is None:
                continue

            if entry.info['SVTYPE'] == 'DEL' and set(entry.alts[0]) == {'N'}:
                entry.alts = entry.ref[0]

            if entry.info['SVTYPE'] != 'INV':
                entry.info['SVLEN'] = abs(len(entry.ref) - len(entry.alts[0]))
            else:
                entry.info['SVLEN'] = len(entry.ref)
            # No more blank genotypes
            n_gt = tuple([_ if _ is not None else 0 for _ in entry.samples[0]['GT']])
            # Preserve phasing informatino
            is_phased = entry.samples[0].phased
            entry.samples[0]['GT'] = n_gt
            entry.samples[0].phased = is_phased

            entry.translate(n_header)
            try:
                out.write(entry)
            except Exception:
                sys.stderr.write(f"{entry}\n{type(entry)}\n")