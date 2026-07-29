#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 5/6/24
'''
import gzip
import sys
import pysam

if __name__ == '__main__':
    # in_file = sys.argv[1]
    ref = pysam.FastaFile(sys.argv[2])
    # in_file = sys.stdin
    # svout = sys.stdout

    headers = []

    sv_list = []
    for line in sys.stdin:
        if '#' in line:
            headers.append(line.strip())
            continue
        entries = line.strip().split('\t')
        delta = len(entries[4]) - len(entries[3])
        if 1000000 > abs(delta) >= 50:
            svtype = 'INS' if delta > 0 else 'DEL'
            entries[6] = "PASS"
            entries[5] = '3'
            entries[2] = f'dipcall_{entries[0]}_{entries[1]}_{svtype}_{abs(delta)}'
            end = int(entries[1]) + 1 if svtype == 'INS' else int(entries[1]) + abs(delta)
            ## DEL
            # if entries[4] == '*':
                # seq_prime = ref.fetch(entries[0], int(entries[1]), end)
                # seq = ""
                # for i in range(len(seq_prime)):
                #     c = seq_prime[i].upper()
                #     if c != 'A' and c != 'C' and c != 'G' and c != 'T':
                #         seq = seq + 'N'
                #     else:
                #         seq = seq + seq_prime[i]
                # entries[4] = seq
                # print(entries[4])
            ## INS
            if entries[3] == '*' or entries[4] == '*':
                continue
                # seq_prime = ref.fetch(entries[0], int(entries[1]), end)
                # seq = ""
                # for i in range(len(seq_prime)):
                #     c = seq_prime[i].upper()
                #     if c != 'A' and c != 'C' and c != 'G' and c != 'T':
                #         seq = seq + 'N'
                #     else:
                #         seq = seq + seq_prime[i]
                # entries[3] = seq

            entries[7] = f'SVTYPE={svtype};SVLEN={delta};END={end}'

            sv_list.append('\t'.join(entries))
    header_last = headers[-1]
    headers.pop()
    headers.append("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">")
    headers.append("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">")
    headers.append("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variation\">")
    headers.append(header_last)
    for val in headers:
        print(val, file=sys.stdout)
    for sv in sv_list:
        print(sv, file=sys.stdout)