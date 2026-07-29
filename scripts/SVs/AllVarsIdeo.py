#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jdlin@uw.edu

@time: 7/6/26
'''
import json
import gzip
import math
from matplotlib.colors import Normalize

from intervaltree import IntervalTree
import joblib as jl
import pysam
import numpy as np
import plotly.io as pio
import pandas as pd
import seaborn as sns
from brokenaxes import brokenaxes
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from plotly.subplots import make_subplots
import pickle
from adjustText import adjust_text
import matplotlib as mpl
from plotly import graph_objects as go
import plotly.express as px
import re
from collections import Counter

from helpers.Constants import *
from helpers.VcfFunc import read_svsize_bed
from helpers.Functions import nonoverlap_intervals

new_rc_params = {'text.usetex': False,
"svg.fonttype": 'none',
"font.family": "sans-serif",
"font.sans-serif": "Arial"
}
mpl.rcParams.update(new_rc_params)

def load_ranges(range_tbl, color=''):
    xranges = []
    kwargs = {'color': [], 'edgecolor': [], 'facecolor': []}
    for idx, row in range_tbl.iterrows():
        this_start = row['start']
        this_end = row['end'] - row['start']
        xranges.append((this_start, this_end))
        if color == '':
            kwargs['facecolor'].append(row['color'])
            kwargs['edgecolor'].append(row['color'])
        else:
            kwargs['facecolor'].append(color)
            kwargs['edgecolor'].append(color)
    return xranges, kwargs

def load_cytoband(range_tbl):
    xranges = []
    kwargs = {'color': [], 'edgecolor': [], 'facecolor': []}
    for idx, row in range_tbl.iterrows():
        this_start = row['start']
        this_end = row['end'] - row['start']
        xranges.append((this_start, this_end))
        # r, g, b = CYTO_COLOR[row['band']]
        # color = rgb_to_hex(int(r), int(g), int(b))
        kwargs['facecolor'].append(CYTO_COLOR[row['band']])
        kwargs['edgecolor'].append('k')
    # legends = [Patch(label=key, color=val) for key, val in legend_dict.items()]

    return xranges, kwargs



def all_var_ideogram():
    chrom_size = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13_size.txt', sep='\t', usecols=[0, 2], names=['chrom', 'size'])
    chrom_cens = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13v2.0_cytobands_CEN.txt', sep='\t', usecols=[0, 1, 2],
                             names=['chrom', 'start', 'end'])
    chrom_cyto = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13v2.0_cytobands_allchrs.txt', sep='\t', usecols=[0, 1, 2, 4],
                             names=['chrom', 'start', 'end', 'band'])

    # chrom_sds = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13v2.0_SD.bed', sep='\t', usecols=[0, 1, 2],
    #                         names=['chrom', 'start', 'end'])

    ## Load SV hotspots
    hotspot = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/pct09_dist500/tables/hotspots.tsv', sep='\t', usecols=[0, 1, 2], skiprows=1, names=['chrom', 'start', 'end'])

    ## Load INV
    dp_inv = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/dp_inv/invCallset.bed', names=['chrom', 'start', 'end'], usecols=[0,1,2], sep='\t')

    ## Load SD variation
    sd_annot = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/DongAhnSD_Var/sd_var_stats2.tsv', header=0, sep='\t')
    # sd_annot = sd_annot.loc[sd_annot['max_min_diff']<5]
    other_sds = sd_annot.loc[sd_annot['max_min_diff'] < 1]
    sd_annot = sd_annot.loc[sd_annot['max_min_diff'] >= 1]

    ## Load contig breaks
    # brk_tbl = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/break_synteny/all_sample_brks_100kb.bed', sep='\t', names=['chrom', 'start', 'end', 'cnt'])
    # brk_tbl = brk_tbl.loc[brk_tbl['cnt'] > 0]

    # sd_chroms = list(sd_annot['chrom'].unique())

    # cmap = mpl.colormaps.get_cmap("RdBu")
    cmap = mpl.colormaps.get_cmap("YlOrRd")
    normalized_indices = np.linspace(0, 1, 50)
    rgba_colors = cmap(normalized_indices)
    # hex_colors = [mpl.colors.to_hex(color) for color in rgba_colors[::-1]]
    hex_colors = [mpl.colors.to_hex(color) for color in rgba_colors]

    bins = np.linspace(1, sd_annot['max_min_diff'].max(), 51)
    sd_annot['color'] = pd.cut(sd_annot['max_min_diff'], bins, labels=hex_colors)


    cmap_brk = mpl.colormaps.get_cmap("coolwarm")

    rgba_colors = cmap_brk(normalized_indices)
    hex_colors = [mpl.colors.to_hex(color) for color in rgba_colors]

    # bins = np.linspace(0, brk_tbl['cnt'].max(), 51)
    # brk_tbl['color'] = pd.cut(brk_tbl['cnt'], bins, labels=hex_colors)

    var_start = 0
    ylabels = []
    yticks = []
    fig, ax = plt.subplots(figsize=(24, 12))
    # fig, ax = plt.subplots(figsize=(20, 3))
    # fig, ax = plt.subplots()

    for idx, row in chrom_size.iterrows():
        chrom, size = row['chrom'], row['size']
        if chrom in ['chrX', 'chrY']:
            continue

        # if chrom != 'chr10':
        #     continue

        yticks.append(var_start + 0.01)

        # brk_ranges, brk_kwargs = load_ranges(brk_tbl.loc[brk_tbl['chrom'] == chrom])
        inv_ranges, inv_kwargs = load_ranges(dp_inv.loc[dp_inv['chrom'] == chrom], '#66bd63')
        hotspot_ranges, hotspot_kwargs = load_ranges(hotspot.loc[hotspot['chrom'] == chrom], '#e78ac3')

        # cen_ranges, cen_kwargs = load_ranges(chrom_cens.loc[chrom_cens['chrom'] == chrom], '#f0f0f0')
        cyto_ranges, cyto_kwargs = load_cytoband(chrom_cyto.loc[chrom_cyto['chrom']==chrom])
        # acro_ranges, acro_kwargs = load_ranges(acro_pos.loc[acro_pos['chrom'] == chrom], '#f0f0f0')
        sd_ranges, sd_kwargs = load_ranges(sd_annot.loc[sd_annot['chrom'] == chrom])
        other_sd_ranges, other_sd_kwargs = load_ranges(other_sds.loc[other_sds['chrom'] == chrom], '#8da0cb')

        # ax.broken_barh([(0, size)], (var_start, 0.02), color='#f0f0f0')


        # ax.broken_barh(acro_ranges, (var_start, 0.02), **acro_kwargs)
        # inv_ranges, inv_kwargs= load_sv_ranges(inv_tbl.loc[inv_tbl['CHR']==chrom], 10000, 10000000, SVCOLOR['INV'])
        # ax.broken_barh(inv_ranges, (var_start, 0.02), **inv_kwargs)

        # ax.broken_barh(hotspot_ranges, (var_start, 0.02), **hotspot_kwargs)

        # ax.broken_barh(svnum_ranges, (var_start, 0.02), **svnum_kwargs)
        # ax.broken_barh(other_sd_ranges, (var_start, 0.02), **other_sd_kwargs)
        ax.broken_barh(cyto_ranges, (var_start, 0.015), **cyto_kwargs)

        ax.broken_barh(sd_ranges, (var_start + 0.016, 0.025), **sd_kwargs)
        ax.broken_barh(other_sd_ranges, (var_start + 0.016, 0.025), **other_sd_kwargs)

        ax.broken_barh(inv_ranges, (var_start + 0.042, 0.025), **inv_kwargs)


        # ax.broken_barh(brk_ranges, (var_start + 0.081, 0.02), **brk_kwargs)

        var_start += 0.08

        ylabels.append(row['chrom'])

    # ax1 = ax.inset_axes([0.7, 0.2, 0.01, 0.45])
    # new_cmap = mpl.colors.ListedColormap(hex_colors)
    # norm = mpl.colors.BoundaryNorm(normalized_indices, new_cmap.N)
    # cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=new_cmap), cax=ax1, orientation='vertical')
    # cbar.set_label('Difference between max and min (Mbp)', size=16)
    # cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    # cbar.set_ticklabels([1, 2, 3, 4, 5, 6], fontsize=14)

    ranges = np.linspace(0, sd_annot['max_min_diff'].max(), 6)

    # cbar.set_ticklabels([int(ele) for ele in ranges], fontsize=14)


    ax.set_yticks(yticks, labels=ylabels, fontsize=22)
    # ax.set_yticks([])
    # ax.legend(handles=[Patch(label='INV', color='#a6d854'), Patch(label='SD range < 1Mb', color='#8da0cb'), Patch(label='Contig breaks', color='#3690c0'),
    #                    Patch(label='INS/DEL hotspot', color='#e78ac3')], loc='lower right', fontsize=16)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    # ax.legend(handles=legends, loc='lower right', fontsize=14)
    # ax.set_xlim([0, chrom_size['size'].max() + 5000000])

    plt.margins(x=0)
    plt.margins(y=0.01)
    ax.xaxis.set_label_position('top')
    ax.set_xticks([])
    ax.spines[['bottom', 'top', 'left', 'right']].set_visible(False)
    fig.tight_layout()
    # fig.savefig(f'{VOL28}/SVREF/hq_asm/CHM13/DongAhnSD_Var/var_ideogram.png', dpi=300)
    fig.savefig(f'{VOL28}/SVREF/hq_asm/Figures/var_ideogram_mainPanel2.svg')
    # fig.savefig(f'{VOL28}/SVREF/hq_asm/Figures/chr10_var_ideogram.png', dpi=300)
    plt.show()


def selected_chrom_ideogram():
    chrom_size = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13_size.txt', sep='\t', usecols=[0, 2], names=['chrom', 'size'])
    chrom_cens = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13v2.0_cytobands_CEN.txt', sep='\t', usecols=[0, 1, 2],
                             names=['chrom', 'start', 'end'])
    chrom_cyto = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13v2.0_cytobands_allchrs.txt', sep='\t', usecols=[0, 1, 2, 4],
                             names=['chrom', 'start', 'end', 'band'])

    # chrom_sds = pd.read_csv(f'{VOL28}/ref_diff/chm13/chm13v2.0_SD.bed', sep='\t', usecols=[0, 1, 2],
    #                         names=['chrom', 'start', 'end'])

    ## Load SV hotspots
    hotspot = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/pct09_dist500/tables/hotspots.tsv', sep='\t', usecols=[0, 1, 2], skiprows=1, names=['chrom', 'start', 'end'])

    ## Load INV
    dp_inv = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/dp_inv/invCallset.bed', names=['chrom', 'start', 'end'], usecols=[0,1,2], sep='\t')

    ## Load SD variation
    sd_annot = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/DongAhnSD_Var/sd_var_stats2.tsv', header=0, sep='\t')
    # sd_annot = sd_annot.loc[sd_annot['max_min_diff']<5]
    other_sds = sd_annot.loc[sd_annot['max_min_diff'] < 1]
    sd_annot = sd_annot.loc[sd_annot['max_min_diff'] >= 1]

    ## Load contig breaks
    # brk_tbl = pd.read_csv(f'{VOL28}/SVREF/hq_asm/CHM13/break_synteny/all_sample_brks_100kb.bed', sep='\t', names=['chrom', 'start', 'end', 'cnt'])
    # brk_tbl = brk_tbl.loc[brk_tbl['cnt'] > 0]

    # sd_chroms = list(sd_annot['chrom'].unique())

    # cmap = mpl.colormaps.get_cmap("RdBu")
    cmap = mpl.colormaps.get_cmap("YlOrRd")
    normalized_indices = np.linspace(0, 1, 50)
    rgba_colors = cmap(normalized_indices)
    # hex_colors = [mpl.colors.to_hex(color) for color in rgba_colors[::-1]]
    hex_colors = [mpl.colors.to_hex(color) for color in rgba_colors]

    bins = np.linspace(1, sd_annot['max_min_diff'].max(), 51)
    sd_annot['color'] = pd.cut(sd_annot['max_min_diff'], bins, labels=hex_colors)


    cmap_brk = mpl.colormaps.get_cmap("coolwarm")

    rgba_colors = cmap_brk(normalized_indices)
    hex_colors = [mpl.colors.to_hex(color) for color in rgba_colors]

    # bins = np.linspace(0, brk_tbl['cnt'].max(), 51)
    # brk_tbl['color'] = pd.cut(brk_tbl['cnt'], bins, labels=hex_colors)

    var_start = 0
    ylabels = []
    yticks = []
    # fig, ax = plt.subplots(figsize=(24, 12))
    fig, ax = plt.subplots(figsize=(20, 5))
    # fig, ax = plt.subplots()

    for idx, row in chrom_size.iterrows():
        chrom, size = row['chrom'], row['size']
        if chrom not in ['chr1', 'chr10', 'chr16']:
            continue

        # if chrom != 'chr10':
        #     continue

        yticks.append(var_start + 0.01)

        # brk_ranges, brk_kwargs = load_ranges(brk_tbl.loc[brk_tbl['chrom'] == chrom])
        inv_ranges, inv_kwargs = load_ranges(dp_inv.loc[dp_inv['chrom'] == chrom], '#66bd63')
        hotspot_ranges, hotspot_kwargs = load_ranges(hotspot.loc[hotspot['chrom'] == chrom], '#e78ac3')

        # cen_ranges, cen_kwargs = load_ranges(chrom_cens.loc[chrom_cens['chrom'] == chrom], '#f0f0f0')
        cyto_ranges, cyto_kwargs = load_cytoband(chrom_cyto.loc[chrom_cyto['chrom']==chrom])
        # acro_ranges, acro_kwargs = load_ranges(acro_pos.loc[acro_pos['chrom'] == chrom], '#f0f0f0')
        sd_ranges, sd_kwargs = load_ranges(sd_annot.loc[sd_annot['chrom'] == chrom])
        other_sd_ranges, other_sd_kwargs = load_ranges(other_sds.loc[other_sds['chrom'] == chrom], '#8da0cb')

        # ax.broken_barh([(0, size)], (var_start, 0.02), color='#f0f0f0')


        # ax.broken_barh(acro_ranges, (var_start, 0.02), **acro_kwargs)
        # inv_ranges, inv_kwargs= load_sv_ranges(inv_tbl.loc[inv_tbl['CHR']==chrom], 10000, 10000000, SVCOLOR['INV'])
        # ax.broken_barh(inv_ranges, (var_start, 0.02), **inv_kwargs)

        # ax.broken_barh(hotspot_ranges, (var_start, 0.02), **hotspot_kwargs)

        # ax.broken_barh(svnum_ranges, (var_start, 0.02), **svnum_kwargs)
        # ax.broken_barh(other_sd_ranges, (var_start, 0.02), **other_sd_kwargs)
        ax.broken_barh(cyto_ranges, (var_start, 0.015), **cyto_kwargs)

        ax.broken_barh(sd_ranges, (var_start + 0.016, 0.025), **sd_kwargs)
        ax.broken_barh(other_sd_ranges, (var_start + 0.016, 0.025), **other_sd_kwargs)

        ax.broken_barh(inv_ranges, (var_start + 0.042, 0.025), **inv_kwargs)


        # ax.broken_barh(brk_ranges, (var_start + 0.081, 0.02), **brk_kwargs)

        var_start += 0.08

        ylabels.append(row['chrom'])

    # ax1 = ax.inset_axes([0.7, 0.2, 0.01, 0.45])
    # new_cmap = mpl.colors.ListedColormap(hex_colors)
    # norm = mpl.colors.BoundaryNorm(normalized_indices, new_cmap.N)
    # cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=new_cmap), cax=ax1, orientation='vertical')
    # cbar.set_label('Difference between max and min (Mbp)', size=16)
    # cbar.set_ticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    # cbar.set_ticklabels([1, 2, 3, 4, 5, 6], fontsize=14)

    ranges = np.linspace(0, sd_annot['max_min_diff'].max(), 6)

    # cbar.set_ticklabels([int(ele) for ele in ranges], fontsize=14)


    ax.set_yticks(yticks, labels=ylabels, fontsize=22)
    # ax.set_yticks([])
    # ax.legend(handles=[Patch(label='INV', color='#a6d854'), Patch(label='SD range < 1Mb', color='#8da0cb'), Patch(label='Contig breaks', color='#3690c0'),
    #                    Patch(label='INS/DEL hotspot', color='#e78ac3')], loc='lower right', fontsize=16)
    ax.invert_yaxis()
    ax.xaxis.tick_top()
    # ax.legend(handles=legends, loc='lower right', fontsize=14)
    # ax.set_xlim([0, chrom_size['size'].max() + 5000000])

    plt.margins(x=0)
    plt.margins(y=0.01)
    ax.xaxis.set_label_position('top')
    ax.set_xticks([])
    ax.spines[['bottom', 'top', 'left', 'right']].set_visible(False)
    fig.tight_layout()
    # fig.savefig(f'{VOL28}/SVREF/hq_asm/CHM13/DongAhnSD_Var/var_ideogram.png', dpi=300)
    # fig.savefig(f'{VOL28}/SVREF/hq_asm/Figures/var_ideogram_mainPanel2.svg')
    fig.savefig(f'{VOL28}/SVREF/hq_asm/Figures/selected_ideogram.png', dpi=300)
    fig.savefig(f'{VOL28}/SVREF/hq_asm/Figures/selected_ideogram.svg')
    plt.show()



def main():

    # sd_ideogram()
    # all_var_ideogram()

    selected_chrom_ideogram()

if __name__ == '__main__':
    main()
