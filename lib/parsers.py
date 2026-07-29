#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 1/11/25
'''
import pandas as pd
import numpy as np

def parse_vcf_info_column(info_str):
    info_tokens = info_str.split(";")
    info_dict = {}

    for token in info_tokens:
        if "=" not in token:
            continue
        info_dict[token.split('=')[0]] = token.split('=')[1]

    return info_dict
def decode_truvari_supp(n, bits):
    return ''.join(str(1 & int(n) >> i) for i in range(bits))

def create_size_df(ins_length_info, ins_values, del_length_info, del_values, ref, sample):
    df_data_ins = pd.DataFrame()
    df_data_del = pd.DataFrame()

    num = 0
    for x, y in zip(ins_values,[ins_length_info[i] for i in ins_values]):
        num += 1
        df_data_ins.loc[num, "x"] = x
        df_data_ins.loc[num, "y"] = np.log10(y) if y != 0 else 0
        df_data_ins.loc[num, "type"] = "ins"
        df_data_ins.loc[num, "ref"] = ref
        df_data_ins.loc[num, "sample"] = sample
    for x, y in zip(del_values,[del_length_info[i] for i in del_values]):
        num += 1
        df_data_del.loc[num, "x"] = -x
        df_data_del.loc[num, "y"] = np.log10(y) if y != 0 else 0
        df_data_del.loc[num, "type"] = "del"
        df_data_del.loc[num, "ref"] = ref
        df_data_del.loc[num, "sample"] = sample
    df_data = pd.concat([df_data_ins, df_data_del])
    return df_data


def read_svsize_bed(bed_df, ref, sample):
    ins_values1 = list(range(100,1000,10))
    del_values1 = [i * -1 for i in list(range(100,1000,10))]

    ins_values2 = list(range(1000,10000,100))
    del_values2 = [i * -1 for i in list(range(1000,10000,100))]

    ins_length_info1 = {i: 0 for i in ins_values1}
    del_length_info1 = {i: 0 for i in del_values1}

    ins_length_info2 = {i: 0 for i in ins_values2}
    del_length_info2 = {i: 0 for i in del_values2}

    for idx, row in bed_df.iterrows():
        sv_len, svtype = row['SVLEN'], row['SVTYPE']
        if abs(sv_len) > 100000 or abs(sv_len) < 50: continue

        if svtype == 'INS':
            if 1000 >= abs(sv_len) >= 100:
                bin_id = sv_len // 10 * 10
                if bin_id in ins_length_info1:
                    ins_length_info1[bin_id] += 1

            if 10000 >= abs(sv_len) >= 1000:
                bin_id = sv_len // 100 * 100
                if bin_id in ins_length_info2:
                    ins_length_info2[bin_id] += 1

        if svtype == 'DEL':
            sv_len = -sv_len if sv_len > 0 else sv_len
            if 1000 >= abs(sv_len) >= 100:
                bin_id = (sv_len // 10 + 1) * 10
                if bin_id in del_length_info1:
                    del_length_info1[bin_id] += 1

            if 10000 >= abs(sv_len) >= 1000:
                bin_id = (sv_len // 100 + 1) * 100
                if bin_id in del_length_info2:
                    del_length_info2[bin_id] += 1

    size100_df = create_size_df(ins_length_info1,ins_values1,del_length_info1,del_values1,ref,sample)
    size1k_df = create_size_df(ins_length_info2,ins_values2,del_length_info2,del_values2,ref,sample)
    return size100_df, size1k_df
