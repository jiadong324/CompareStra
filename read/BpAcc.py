#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/10/5
'''
import pysam
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

# from hgsvc.Match import *
from helpers.Functions import *
from helpers.Constant import *
from helpers.Annot import *

# def find_nearest_matched_pavs(workdir, platforms, samples):
#
#     long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']
#
#     for sample in samples:
#         sample_dir = f'{workdir}/{sample}'
#         matched_out = f'{sample_dir}/intersect/nrst_match'
#         if not os.path.exists(matched_out):
#             os.mkdir(matched_out)
#
#         for platform in platforms:
#             for caller in long_read_callers:
#                 print(f'Matching {sample} {caller} with PAV calls')
#                 source_bed = f'{sample_dir}/pav_v112/pav_{sample}_CCS_svs.bed'
#                 ngmlr_bed = f'{sample_dir}/{caller}/bed/{sample}.{platform}.ngmlr.{caller}.s5.filtered.bed'
#                 minimap2_bed = f'{sample_dir}/{caller}/bed/{sample}.{platform}.minimap2.{caller}.s5.filtered.bed'
#                 print(f'Producing {platform} com-ngmlr')
#                 matched_file1 = f'{matched_out}/{sample}.{caller}.{platform}.com-ngmlr.nrst.match.tsv'
#                 nearest_by_svlen_overlap(source_bed, ngmlr_bed, 500, 0.5, matched_file1)
#                 print(f'Producing {platform} com-minimap2')
#                 matched_file2 = f'{matched_out}/{sample}.{caller}.{platform}.com-minimap2.nrst.match.tsv'
#                 nearest_by_svlen_overlap(source_bed, minimap2_bed, 500, 0.5, matched_file2)

                # if platform == 'hifi' and caller != 'nanovar':
                #     print(f'Producing {platform} com-pbmm2')
                #     pbmm2_bed = f'{sample_dir}/{caller}/bed/{sample}.{platform}.pbmm2.{caller}.s5.filtered.bed'
                #     matched_file3 = f'{matched_out}/{sample}.{caller}.{platform}.com-pbmm2.nrst.match.tsv'
                #     nearest_by_svlen_overlap(source_bed, pbmm2_bed, 500, 0.5, matched_file3)


def samples_pav_bpshift_stackplot(workdir, samples):

    # caller_colors = [TOOLCOLORS[caller] for caller in long_read_callers]

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    shift_label_colors = ['#DA4302', '#009975', '#9B3675', '#F7B71D']
    shift_legends = [Patch(facecolor=shift_label_colors[::-1][idx], label=shift_labels[::-1][idx]) for idx in range(len(shift_labels))]

    aligners = ['minimap2', 'ngmlr']
    platforms = ['hifi', 'ont']

    fig = plt.figure(figsize=(10, 8))
    # fig.suptitle(f'Breakpoint shift compared to PAV', fontsize=14)

    gs = fig.add_gridspec(2, 2, hspace=0.13, wspace=0.05)
    axes = gs.subplots(sharex='col', sharey='row')

    for aligner_idx in range(len(aligners)):
        aligner = aligners[aligner_idx]
        for platform_idx in range(len(platforms)):
            platform = platforms[platform_idx]

            callers_bpshift = []

            for caller in CALLERS:
                shift_by_cutoff = {'0,10': 0, '10,50': 0, '50,100': 0, '>100': 0}
                this_caller_total_matches = 0
                for sample_idx in range(len(samples)):

                    sample = samples[sample_idx]
                    matched_dir = f'{workdir}/{sample}/intersect/nrst_match'
                    matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'

                    df_matches = pd.read_csv(matched_file, sep='\t', header=[0])
                    this_caller_total_matches += len(df_matches)

                    for idx, row in df_matches.iterrows():
                        bpshift = int(row['BP_DIFF'])
                        if bpshift <= 10:
                            shift_by_cutoff['0,10'] += 1
                        elif bpshift <= 50:
                            shift_by_cutoff['10,50']+= 1
                        elif bpshift <= 100:
                            shift_by_cutoff['50,100']+= 1
                        else:
                            shift_by_cutoff['>100']+= 1

                for shift_label, count in shift_by_cutoff.items():
                    callers_bpshift.append((shift_label, count * 100 / this_caller_total_matches, caller))

            df_caller_bpshift = pd.DataFrame(callers_bpshift, columns=['shift_label', 'pcrt', 'caller'])
            # df_caller_bpshift_groupby_sample = df_caller_bpshift.groupby('sample')
            ax = axes[aligner_idx][platform_idx]

            if aligner_idx == 0:
                ax.set_title(PLATDICT[platform])

            ax.set_ylabel(aligner)
            ax.legend(handles=shift_legends, loc='upper right', title='BpShift')
            each_caller_shifts = [df_caller_bpshift[df_caller_bpshift['shift_label']==shift_label]['pcrt'].values for shift_label in shift_labels]
            ax.stackplot(CALLERS, each_caller_shifts, labels = shift_labels, colors=shift_label_colors)

            ax.set_ylim(0, 100)
            ax.set_yticks(np.linspace(0, 100, 6))
            ax.set_yticklabels([0, '20%', '40%', '60%', '80%', '100%'])

    gs.tight_layout(fig)
    # plt.show()
    # fig.savefig(f'{workdir}/thesis_figures/Figure3/samples_pav_bpshift_stackplot.pdf')
    fig.savefig(f'{workdir}/Figures/Figure5/samples_pav_bpshift_stackplot.pdf')

def samples_pav_bpshift_stackplot2(workdir, samples, aligners, platforms):


    long_read_callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']


    fig = plt.figure(figsize=(11, 4))

    bar_width = 0.35
    bar_gap = 0.04

    r1 = np.arange(len(long_read_callers))
    xticks = [x + (bar_width + bar_gap) / 2 for x in r1]
    r2 = [x + bar_width + bar_gap for x in r1]
    rs = [r1, r2]

    shift_labels = ['0,10', '10,50', '50,100', '>100']


    all_legend = [Patch(facecolor='white', edgecolor='black', linestyle='--', label='minimap2'),
                  Patch(facecolor='white', edgecolor='black', linestyle='-', label='ngmlr')]

    for idx in range(len(shift_labels)):
        all_legend.append(Patch(facecolor=BPSHIFTCOLORS[idx], label=shift_labels[idx]))

    plat_fig_idx = 1
    for platform in platforms:
        aligners_shift = {}
        ax = fig.add_subplot(1, 2, plat_fig_idx)
        for aligner in aligners:
            shift_by_cutoff = {'0,10': np.zeros(len(long_read_callers)), '10,50': np.zeros(len(long_read_callers)),
                               '50,100': np.zeros(len(long_read_callers)), '>100': np.zeros(len(long_read_callers))}

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]
                for sample_fig_idx in range(len(samples)):
                    sample = samples[sample_fig_idx]
                    matched_dir = f'{workdir}/{sample}/intersect/nrst_match'

                    matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'

                    if not os.path.exists(matched_file):
                        continue
                    df_matches = pd.read_csv(matched_file, sep='\t', header=[0])

                    for idx, row in df_matches.iterrows():
                        bpshift = int(row['BP_DIFF'])
                        if bpshift <= 10:
                            shift_by_cutoff['0,10'][caller_idx] += 1
                        elif bpshift <= 50:
                            shift_by_cutoff['10,50'][caller_idx] += 1
                        elif bpshift <= 100:
                            shift_by_cutoff['50,100'][caller_idx] += 1
                        else:
                            shift_by_cutoff['>100'][caller_idx] += 1

            aligners_shift[aligner] = pd.DataFrame(shift_by_cutoff)


        for x_axis_idx in range(len(aligners)):
            aligner = aligners[x_axis_idx]
            df_bpshifts = aligners_shift[aligner]
            totals = [i+j+k+m for i,j,k,m in zip(df_bpshifts['0,10'], df_bpshifts['10,50'],df_bpshifts['50,100'],df_bpshifts['>100'])]
            for shift_label_idx in range(len(shift_labels)):
                label = shift_labels[shift_label_idx]
                this_label_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[label], totals)]
                if shift_label_idx == 0:
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, hatch='//', edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])

                else:
                    bottom_sum = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[0]], totals)]
                    for bottom_idx in range(1, shift_label_idx):
                        bottom_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[bottom_idx]], totals)]
                        for bottom_val_idx in range(len(bottom_sum)):
                            bottom_sum[bottom_val_idx] += bottom_pcrt[bottom_val_idx]
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, hatch='//', edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])

        # x_axis_ticks = [x + (bar_width + bar_gap) / 2 for x in r1]

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 6))
        if plat_fig_idx == 1:
            ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 6)])
        else:
            ax.set_yticklabels('')

        ax.set_xticks(xticks)
        ax.set_xticklabels(long_read_callers)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(color='grey', linestyle=(0, (5, 5)), lw=1, axis='y')

        plat_fig_idx += 1

    plt.tight_layout()
    # plt.show()
    fig.savefig(f'{workdir}/Figures/Figure5/persample_pav_bpshift_stackplot.pdf')

def samples_pav_inaccbp(workdir, samples, simple_reps_path, rmsk_path):
    simple_reps = pysam.Tabixfile(simple_reps_path, 'r')
    rmsk = pysam.Tabixfile(rmsk_path, 'r')

    # long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']
    long_read_callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']

    aligners = ['minimap2', 'ngmlr']

    platforms = ['hifi', 'ont']

    inaccurate_calls = []
    accurate_calls = []
    rptypes = ['VNTR', 'STR', 'LTR', 'LINE', 'SINE', 'None']

    for aligner_idx in range(len(aligners)):
        aligner = aligners[aligner_idx]
        for platform_idx in range(len(platforms)):
            platform = platforms[platform_idx]

            # callers_bpshift = []

            for caller in long_read_callers:
                for sample_idx in range(len(samples)):

                    shift_by_cutoff = {'0,10': 0, '10,50': 0, '50,100': 0, '>100': 0}

                    sample = samples[sample_idx]
                    matched_dir = f'{workdir}/{sample}/intersect/nrst_match'
                    matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'

                    # inaccbp_events = open(f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.inaccbp.tsv', 'w')

                    df_matches = pd.read_csv(matched_file, sep='\t', header=[0])

                    for idx, row in df_matches.iterrows():
                        chrom, start, end = row['#CHROM'], int(row['POS']), int(row['END'])
                        rptype, pcrt = rep_annotation(chrom, start, end, simple_reps, rmsk)

                        bpshift = int(row['BP_DIFF'])
                        if bpshift <= 10:
                            shift_by_cutoff['0,10'] += 1
                            accurate_calls.append((row['SOURCE_ID'], rptype, pcrt, caller, platform, aligner, sample))

                        elif bpshift <= 50:
                            shift_by_cutoff['10,50']+= 1
                            accurate_calls.append((row['SOURCE_ID'], rptype, pcrt, caller, platform, aligner, sample))

                        elif bpshift <= 100:
                            shift_by_cutoff['50,100']+= 1
                            inaccurate_calls.append((row['SOURCE_ID'], rptype, pcrt, caller, platform, aligner, sample))

                        else:
                            shift_by_cutoff['>100']+= 1
                            inaccurate_calls.append((row['SOURCE_ID'], rptype, pcrt, caller, platform, aligner, sample))


                    # for shift_label, count in shift_by_cutoff.items():
                    #     callers_bpshift.append((shift_label, count * 100 / len(df_matches), sample, caller))


    df_inaccurate = pd.DataFrame(inaccurate_calls, columns=['pav_id', 'rptype', 'pcrt', 'caller', 'platform', 'aligner', 'sample'])
    df_accurate = pd.DataFrame(accurate_calls, columns=['pav_id', 'rptype', 'pcrt', 'caller', 'platform', 'aligner', 'sample'])

    all_calls_rptypes = []
    inaccurate_caller_rptypes = []
    accurate_caller_rptypes = []

    inaccurate_rptype_sum = {}
    accurate_rptype_sum = {}

    for caller in long_read_callers:
        inaccurate_rptype_dict = {}
        inaccurate_total = 0

        for idx, row in df_inaccurate[df_inaccurate['caller']==caller].iterrows():
            inaccurate_total += 1
            rptype = row['rptype']
            if rptype not in rptypes:
                rptype = 'Others'
            if rptype in inaccurate_rptype_dict:
                inaccurate_rptype_dict[rptype] += 1
            else:
                inaccurate_rptype_dict[rptype] = 1

        for rptype, count in inaccurate_rptype_dict.items():
            inaccurate_caller_rptypes.append((rptype, count / inaccurate_total, caller))
            if rptype not in accurate_rptype_sum:
                inaccurate_rptype_sum[rptype] = 0
            else:
                inaccurate_rptype_sum[rptype] += count / inaccurate_total
            # all_calls_rptypes.append((rptype, count / inaccurate_total, caller, 'inaccurate'))

        accurate_rptype_dict = {}
        accurate_total = 0

        for idx, row in df_accurate[df_accurate['caller']==caller].iterrows():
            accurate_total += 1
            rptype = row['rptype']
            if rptype not in rptypes:
                rptype = 'Others'
            if rptype in accurate_rptype_dict:
                accurate_rptype_dict[rptype] += 1
            else:
                accurate_rptype_dict[rptype] = 1

        for rptype, count in accurate_rptype_dict.items():
            accurate_caller_rptypes.append((rptype, count / accurate_total, caller))
            if rptype not in accurate_rptype_sum:
                accurate_rptype_sum[rptype] = 0
            else:
                accurate_rptype_sum[rptype] += count / accurate_total
            # all_calls_rptypes.append((rptype, count / accurate_total, caller, 'accurate'))


    # df_inaccurate_caller_rptypes = pd.DataFrame(inaccurate_caller_rptypes, columns=['rptype', 'pcrt', 'caller'])
    # df_accurate_caller_rptypes = pd.DataFrame(accurate_caller_rptypes, columns=['rptype', 'pcrt', 'caller'])
    # df_all_calls_rptypes = pd.DataFrame(all_calls_rptypes, columns=['rptype', 'pcrt', 'caller', 'class'])

    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1,1,1)

    # ax1 = fig.add_subplot(1, 2, 1)
    # sns.barplot(data=df_accurate_caller_rptypes, x='rptype', y='pcrt', hue='caller', ax=ax1, edgecolor='w', ci=None, palette=TOOLCOLORS)
    # ax1.set_title('Repeats of accurate calls')
    #
    # ax2 = fig.add_subplot(1, 2, 2)
    # sns.barplot(data=df_inaccurate_caller_rptypes, x='rptype', y='pcrt', hue='caller', ax=ax2, edgecolor='w', ci=None, palette=TOOLCOLORS)
    # ax2.set_title('Repeats of inaccurate calls')
    #
    # for ax in [ax1, ax2]:
    #     ax.set_ylim(0, 1)
    #     ax.set_yticks(np.linspace(0, 1, 6))
    #     ax.set_yticklabels([0, '20%', '40%', '60%', '80%', '100%'])

    bar_width = 0.35
    bar_gap = 0.04

    xlabels = []
    accurate_rptype_avg = []
    inaccurate_rptype_avg = []

    for rptype, pcrt in accurate_rptype_sum.items():
        if rptype not in xlabels:
            xlabels.append(rptype)
        accurate_rptype_avg.append(pcrt / len(long_read_callers))
        inaccurate_rptype_avg.append(inaccurate_rptype_sum[rptype] / len(long_read_callers))

    r1 = np.arange(len(xlabels))
    xticks = [x + (bar_width + bar_gap) / 2 for x in r1]
    r2 = [x + bar_width + bar_gap for x in r1]

    print(xlabels)
    print(accurate_rptype_avg)
    print(inaccurate_rptype_avg)

    ax.bar(r1, accurate_rptype_avg,  width=bar_width)
    ax.bar(r2, inaccurate_rptype_avg,  width=bar_width)
    ax.plot(r1, accurate_rptype_avg)
    ax.plot(r2, inaccurate_rptype_avg)

    ax.set_ylim(0, 0.6)
    ax.set_yticks(np.linspace(0, 0.6, 4))

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    plt.tight_layout()
    plt.show()

    # fig.savefig(f'{workdir}/thesis_figures/Figure5/samples_bpshift_rptypes.pdf')
    fig.savefig(f'{workdir}/Figures/Figure5/samples_bpshift_rptypes.pdf')

def pav_bpshift_of_sample(workdir, sample):
    long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']

    aligners = ['minimap2', 'ngmlr']

    platforms = ['hifi', 'ont']

    fig = plt.figure(figsize=(11, 4))
    fig.suptitle(f'{sample} breakpoint compared to PAV', fontsize=14)

    gs = fig.add_gridspec(1, 2, hspace=0.13, wspace=0.05)
    axes = gs.subplots(sharex='col', sharey='row')

    shift_labels = ['0,10', '10,50', '50,100', '>100']
    shift_label_colors = ['#DA4302', '#009975', '#9B3675', '#F7B71D']

    all_legend = [Patch(facecolor='white', edgecolor='black', linestyle='--', label='minimap2'), Patch(facecolor='white', edgecolor='black', label='ngmlr')]

    for idx in range(len(shift_labels)):
        all_legend.append(Patch(facecolor=shift_label_colors[idx], label=shift_labels[idx]))

    for plat_fig_idx in range(len(platforms)):
        platform = platforms[plat_fig_idx]
        matched_dir = f'{workdir}/{sample}/intersect/nrst_match'
        aligners_shift = {}

        for aligner in aligners:

            shift_by_cutoff = {'0,10': np.zeros(len(long_read_callers)), '10,50': np.zeros(len(long_read_callers)),
                               '50,100': np.zeros(len(long_read_callers)), '>100': np.zeros(len(long_read_callers))}

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]

                matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'

                if not os.path.exists(matched_file):
                    continue
                df_matches = pd.read_csv(matched_file, sep='\t', header=[0])

                for idx, row in df_matches.iterrows():
                    bpshift = int(row['BP_DIFF'])
                    if bpshift <= 10:
                        shift_by_cutoff['0,10'][caller_idx] += 1
                    elif bpshift <= 50:
                        shift_by_cutoff['10,50'][caller_idx] += 1
                    elif bpshift <= 100:
                        shift_by_cutoff['50,100'][caller_idx] += 1
                    else:
                        shift_by_cutoff['>100'][caller_idx] += 1

            aligners_shift[aligner] = pd.DataFrame(shift_by_cutoff)

        ax = axes[plat_fig_idx]
        bar_width = 0.4
        bar_gap = 0.05
        r1 = np.arange(len(long_read_callers))
        r2 = [x + bar_width + bar_gap  for x in r1]
        rs = [r1, r2]
        for x_axis_idx in range(len(aligners)):
            aligner = aligners[x_axis_idx]
            df_bpshifts = aligners_shift[aligner]
            totals = [i+j+k+m for i,j,k,m in zip(df_bpshifts['0,10'], df_bpshifts['10,50'],df_bpshifts['50,100'],df_bpshifts['>100'])]
            for shift_label_idx in range(len(shift_labels)):
                label = shift_labels[shift_label_idx]
                this_label_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[label], totals)]
                if shift_label_idx == 0:
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, edgecolor='black', linestyle="--", width=bar_width, color=shift_label_colors[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, edgecolor='black', width=bar_width, color=shift_label_colors[shift_label_idx])
                else:
                    bottom_sum = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[0]], totals)]
                    for bottom_idx in range(1, shift_label_idx):
                        bottom_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[bottom_idx]], totals)]
                        for bottom_val_idx in range(len(bottom_sum)):
                            bottom_sum[bottom_val_idx] += bottom_pcrt[bottom_val_idx]
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, edgecolor='black', linestyle="--", width=bar_width, color=shift_label_colors[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, edgecolor='black', width=bar_width, color=shift_label_colors[shift_label_idx])

        x_axis_ticks = [x + (bar_width + bar_gap) / 2 for x in r1]
        ax.set_title(PLATDICT[platform])
        if plat_fig_idx == 0:
            ax.set_ylabel('Percent')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 6))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 6)])
        ax.set_xticks(x_axis_ticks)
        ax.set_xticklabels(long_read_callers)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(color='grey', linestyle=(0, (5, 5)), lw=1, axis='y')

        if plat_fig_idx == 1:
            ax.legend(handles=all_legend, loc='upper left', bbox_to_anchor=(1,1), ncol=1)

    gs.tight_layout(fig)
    # plt.show()
    fig.savefig(f'{workdir}/{sample}/figures/{sample}.callers.com.bpshift.barplot.pdf')

def find_match_ngs(workdir, platforms, samples):
    long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']
    for sample in samples:
        for platform in platforms:
            for caller in long_read_callers:

                manta_bed = f'{workdir}/{sample}/manta/manta.passonly.exbnd.bed'

                bed_dir = f'{workdir}/{sample}/{caller}/bed'

                minimap2_bed = f'{bed_dir}/{sample}.{platform}.minimap2.{caller}.s5.filtered.bed'
                matched_file1 = f'{bed_dir}/{sample}.{platform}.{caller}.minimap2-manta.nrst.match.tsv'
                print(f'Producing {caller}-{platform} minimap2-manta matched calls')
                short_long_nearest_match(manta_bed, minimap2_bed, 500, 0.5, matched_file1)

                ngmlr_bed = f'{bed_dir}/{sample}.{platform}.ngmlr.{caller}.s5.filtered.bed'
                matched_file2 = f'{bed_dir}/{sample}.{platform}.{caller}.ngmlr-manta.nrst.match.tsv'
                print(f'Producing {caller}-{platform} ngmlr-manta matched calls')

                short_long_nearest_match(manta_bed, ngmlr_bed, 500, 0.5, matched_file2)

def sample_ngs_bpshift_stackplot(workdir, samples, platforms):

    long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']

    aligners = ['minimap2', 'ngmlr']
    shift_labels = ['0,10', '10,50', '50,100', '>100']
    shift_label_colors = ['#DA4302', '#009975', '#9B3675', '#F7B71D']
    shift_legends = [Patch(facecolor=shift_label_colors[::-1][idx], label=shift_labels[::-1][idx]) for idx in
                     range(len(shift_labels))]


    fig = plt.figure(figsize=(10, 8))

    gs = fig.add_gridspec(2, 2, hspace=0.13, wspace=0.05)
    axes = gs.subplots(sharex='col', sharey='row')


    for aligner_idx in range(len(aligners)):
        aligner = aligners[aligner_idx]
        for platform_idx in range(len(platforms)):
            platform = platforms[platform_idx]

            callers_bpshift = []
            for caller_idx in range(len(long_read_callers)):
                shift_by_cutoff = {'0,10': 0, '10,50': 0, '50,100': 0, '>100': 0}
                this_caller_total_matches = 0
                caller = long_read_callers[caller_idx]

                for sample_idx in range(len(samples)):
                    sample = samples[sample_idx]
                    matched_dir = f'{workdir}/{sample}/{caller}/bed'
                    matched_file = f'{matched_dir}/{sample}.{platform}.{caller}.{aligner}-manta.nrst.match.tsv'

                    df_matches = pd.read_csv(matched_file, sep='\t', header=[0])
                    this_caller_total_matches += len(df_matches)

                    for idx, row in df_matches.iterrows():
                        bpshift = int(row['BP_DIFF'])
                        if bpshift <= 10:
                            shift_by_cutoff['0,10'] += 1
                        elif bpshift <= 50:
                            shift_by_cutoff['10,50'] += 1
                        elif bpshift <= 100:
                            shift_by_cutoff['50,100'] += 1
                        else:
                            shift_by_cutoff['>100'] += 1

                for shift_label, count in shift_by_cutoff.items():
                    callers_bpshift.append((shift_label, count * 100 / this_caller_total_matches, caller))

            df_caller_bpshift = pd.DataFrame(callers_bpshift, columns=['shift_label', 'pcrt', 'caller'])
            # df_caller_bpshift_groupby_sample = df_caller_bpshift.groupby('sample')
            ax = axes[aligner_idx][platform_idx]

            if aligner_idx == 0:
                ax.set_title(PLATDICT[platform])

            ax.set_ylabel(aligner)
            ax.legend(handles=shift_legends, loc='upper right', title='BpShift')
            each_caller_shifts = [df_caller_bpshift[df_caller_bpshift['shift_label'] == shift_label]['pcrt'].values for
                                  shift_label in shift_labels]
            ax.stackplot(long_read_callers, each_caller_shifts, labels=shift_labels, colors=shift_label_colors)

            ax.set_ylim(0, 100)
            ax.set_yticks(np.linspace(0, 100, 6))
            ax.set_yticklabels([0, '20%', '40%', '60%', '80%', '100%'])

    gs.tight_layout(fig)
    # plt.show()
    fig.savefig(f'{workdir}/figures/Figure3/samples_ngs_bpshift_stackplot.pdf')


def sample_ngs_bpshift_stackplot2(workdir, samples, platforms):
    long_read_callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']

    aligners = ['minimap2', 'ngmlr']
    shift_labels = ['0,10', '10,50', '50,100', '>100']



    fig = plt.figure(figsize=(11, 4))

    bar_width = 0.35
    bar_gap = 0.04

    r1 = np.arange(len(long_read_callers))
    xticks = [x + (bar_width + bar_gap) / 2 for x in r1]
    r2 = [x + bar_width + bar_gap for x in r1]
    rs = [r1, r2]

    plat_fig_idx = 1

    for platform in platforms:
        aligners_shift = {}
        ax = fig.add_subplot(1, 2, plat_fig_idx)
        for aligner in aligners:
            shift_by_cutoff = {'0,10': np.zeros(len(long_read_callers)), '10,50': np.zeros(len(long_read_callers)),
                               '50,100': np.zeros(len(long_read_callers)), '>100': np.zeros(len(long_read_callers))}

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]
                for sample_idx in range(len(samples)):
                    sample = samples[sample_idx]
                    matched_dir = f'{workdir}/{sample}/{caller}/bed'
                    matched_file = f'{matched_dir}/{sample}.{platform}.{caller}.{aligner}-manta.nrst.match.tsv'

                    df_matches = pd.read_csv(matched_file, sep='\t', header=[0])

                    for idx, row in df_matches.iterrows():
                        bpshift = int(row['BP_DIFF'])
                        if bpshift <= 10:
                            shift_by_cutoff['0,10'] += 1
                        elif bpshift <= 50:
                            shift_by_cutoff['10,50'] += 1
                        elif bpshift <= 100:
                            shift_by_cutoff['50,100'] += 1
                        else:
                            shift_by_cutoff['>100'] += 1

            aligners_shift[aligner] = pd.DataFrame(shift_by_cutoff)

        for x_axis_idx in range(len(aligners)):
            aligner = aligners[x_axis_idx]
            df_bpshifts = aligners_shift[aligner]
            totals = [i+j+k+m for i,j,k,m in zip(df_bpshifts['0,10'], df_bpshifts['10,50'],df_bpshifts['50,100'],df_bpshifts['>100'])]
            for shift_label_idx in range(len(shift_labels)):
                label = shift_labels[shift_label_idx]
                this_label_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[label], totals)]
                if shift_label_idx == 0:
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, hatch='//', edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])

                else:
                    bottom_sum = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[0]], totals)]
                    for bottom_idx in range(1, shift_label_idx):
                        bottom_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[bottom_idx]], totals)]
                        for bottom_val_idx in range(len(bottom_sum)):
                            bottom_sum[bottom_val_idx] += bottom_pcrt[bottom_val_idx]
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, hatch='//', edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, edgecolor='black', width=bar_width, color=BPSHIFTCOLORS[shift_label_idx])

        # x_axis_ticks = [x + (bar_width + bar_gap) / 2 for x in r1]

        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 6))
        if plat_fig_idx == 1:
            ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 6)])
        else:
            ax.set_yticklabels('')

        ax.set_xticks(xticks)
        ax.set_xticklabels(long_read_callers)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(color='grey', linestyle=(0, (5, 5)), lw=1, axis='y')

        plat_fig_idx += 1

    plt.tight_layout()
    # plt.show()
    fig.savefig(f'{workdir}/Figures/Figure5/persample_ngs_bpshift_stackplot.pdf')


def merge_samples(vcf_dir, merge_dir, samples):
    long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']

    platforms = ['hifi', 'ont']
    aligners = ['minimap2', 'ngmlr']

    for caller in long_read_callers:
        this_caller_merge_dir = f'{merge_dir}/{caller}'
        if not os.path.exists(this_caller_merge_dir):
            os.mkdir(this_caller_merge_dir)

        for platform in platforms:
            for aligner in aligners:
                sample_file = f'{this_caller_merge_dir}/{caller}-{platform}-{aligner}_sammples_file'
                sample_file_writer = open(sample_file, 'w')
                for sample in samples:
                    vcf_file = f'{vcf_dir}/{sample}/{caller}/{sample}.{platform}.{aligner}.{caller}.s5.vcf'
                    print(vcf_file, file=sample_file_writer)
                sample_file_writer.close()

                merged_out_vcf = f'{this_caller_merge_dir}/{platform}.{aligner}.samples_merged.vcf'
                print(f'Producing {caller}-{aligner}-{platform} samples overlaps ...')

                cmd = f'~/Biotools/SURVIVOR/Debug/SURVIVOR merge {sample_file} 500 3 0 0 0 50 {merged_out_vcf}'
                os.system(cmd)

                # os.remove(sample_file)


def reccurent_sv_bpdiff(merge_dir):

    # long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']
    long_read_callers = ['pbsv', 'svim', 'cutesv', 'sniffles', 'nanovar']

    platforms = ['hifi', 'ont']
    aligners = ['minimap2', 'ngmlr']

    bpdiff_acc = {}
    bpdiff_inacc = {}

    for platform in platforms:
        for aligner in aligners:
            for idx in range(len(long_read_callers)):
                caller = long_read_callers[idx]
                this_caller_merge_dir = f'{merge_dir}/{caller}'
                merged_out_vcf = f'{this_caller_merge_dir}/{platform}.{aligner}.samples_merged.vcf'

                for line in open(merged_out_vcf, 'r'):
                    if '#' in line:
                        continue
                    entries = line.strip().split('\t')
                    chrom, start, id = entries[0], int(entries[1]), entries[2]

                    info_tokens = entries[7].split(';')
                    info_dict = {}
                    for token in info_tokens:
                        info_dict[token.split('=')[0]] = token.split('=')[1]

                    if info_dict['SVTYPE'] == 'TRA':
                        continue
                    if 'chr' not in chrom:
                        chrom = f'chr{chrom}'

                    if chrom not in VALID_CHROMS:
                        continue

                    starts = []
                    ends = []

                    sample1_co_str_tokens = entries[9].split(':')[-1].split(',')
                    for sample1_co_str in sample1_co_str_tokens:
                        sample1_start, sample1_end = int(sample1_co_str.split('-')[0].split('_')[1]), int(sample1_co_str.split('-')[1].split('_')[1])
                        starts.append(sample1_start)
                        ends.append(sample1_end)

                    sample2_co_str_tokens = entries[10].split(':')[-1].split(',')
                    for sample2_co_str in sample2_co_str_tokens:
                        sample2_start, sample2_end = int(sample2_co_str.split('-')[0].split('_')[1]), int(sample2_co_str.split('-')[1].split('_')[1])

                        starts.append(sample2_start)
                        ends.append(sample2_end)

                    sample3_co_str_tokens = entries[11].split(':')[-1].split(',')
                    for sample3_co_str in sample3_co_str_tokens:
                        sample3_start, sample3_end = int(sample3_co_str.split('-')[0].split('_')[1]), int(sample3_co_str.split('-')[1].split('_')[1])
                        starts.append(sample3_start)
                        ends.append(sample3_end)

                    start_std = abs(np.std(starts))
                    end_std = abs(np.std(ends))

                    p_a_key = f'{platform}-{aligner}'

                    if start_std <= 50 and end_std <= 50:
                        if p_a_key in bpdiff_acc:
                            bpdiff_acc[p_a_key][idx] += 1
                        else:
                            bpdiff_acc[p_a_key] = np.zeros(len(long_read_callers))
                            bpdiff_acc[p_a_key][idx] += 1
                    else:
                        if p_a_key in bpdiff_inacc:
                            bpdiff_inacc[p_a_key][idx] += 1
                        else:
                            bpdiff_inacc[p_a_key] = np.zeros(len(long_read_callers))
                            bpdiff_inacc[p_a_key][idx] += 1



    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)
    x_ticks = np.arange(len(long_read_callers))

    platform_ls = {'hifi': '-', 'ont': '--'}
    platform_mk = {'hifi': 'o', 'ont': 'X'}

    platform_legend = [Line2D([], [], label=plat, color='black', ls=platform_ls[plat], marker=platform_mk[plat]) for plat in platforms]

    for p_a_key, acc_count in bpdiff_acc.items():
        platform = p_a_key.split('-')[0]
        aligner = p_a_key.split('-')[1]

        inacc_count = bpdiff_inacc[p_a_key]
        total = [i + j for i, j in zip(acc_count, inacc_count)]
        acc_pcrt = [i/j for i, j in zip(acc_count, total)]
        # inacc_pcrt = [i/j for i, j in zip(inacc_count, total)]
        ax.plot(x_ticks, acc_pcrt, color=ALIGNERCOLOR[aligner], ls=platform_ls[platform], marker=platform_mk[platform], lw=2, markersize=10)


    # ax.bar(x_ticks, inacc_pcrt, bottom=acc_pcrt)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(long_read_callers)
    ax.set_ylim(0.4, 1)

    # ax.grid(linewidth=1.5, linestyle='--', axis='y')

    ax.set_yticks(np.linspace(0.4, 1, 4))
    ax.set_yticklabels([round(val, 2) for val in np.linspace(0.4 , 1, 4)])

    ax.legend(handles=platform_legend)

    plt.tight_layout()
    # plt.show()
    # fig.savefig(f'{MAC}/thesis_figures/Figure6/samples_recurrent_sv_bpshift.pdf')
    fig.savefig(f'{MAC}/Figures/Figure5/samples_recurrent_sv_bpshift.pdf')


'''

def pieplot_pav_bpshift_count_platform(workdir, samples, platform):
    long_read_callers = ['pbsv', 'svision', 'svim', 'cutesv', 'sniffles', 'nanovar']

    aligners_dict = {'hifi': ['minimap2', 'ngmlr'], 'ont': ['minimap2', 'ngmlr']}
    aligners = aligners_dict[platform]

    for sample in samples:
        print(f'Plot {sample}-{platform} ...')
        matched_dir = f'{workdir}/{sample}/intersect/nrst_match'
        aligner_bpshift_dict = {}
        for aligner in aligners:
            # matched_calls = {}
            bpshift_count_dict = {}
            for caller in long_read_callers:
                matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'
                shift_by_cutoff = {'0,10': 0, '10,50': 0, '50,100': 0, '>100': 0}
                if not os.path.exists(matched_file):
                    continue
                df_matches = pd.read_csv(matched_file, sep='\t', header=[0])
                # matched_calls[caller] = len(df_matches)

                for idx, row in df_matches.iterrows():
                    bpshift = int(row['BP_DIFF'])
                    if bpshift <= 10:
                        shift_by_cutoff['0,10'] += 1
                    elif bpshift <= 50:
                        shift_by_cutoff['10,50'] += 1
                    elif bpshift <= 100:
                        shift_by_cutoff['50,100'] += 1
                    else:
                        shift_by_cutoff['>100'] += 1

                bpshift_count_dict[caller] = shift_by_cutoff

            aligner_bpshift_dict[aligner] = bpshift_count_dict


        fig = plt.figure(figsize=(20, 6))
        fig.suptitle(f'Comparison of {sample}-{platform} breakpoint with PAV', fontsize=13)
        fig_dix = 1
        for aligner, bpshift_dict in aligner_bpshift_dict.items():
            for caller in long_read_callers:
                count_dict = bpshift_dict[caller]

                values = []
                labels = []
                for label, count in count_dict.items():
                    values.append(count)
                    labels.append(label)

                ax = fig.add_subplot(len(aligners), len(long_read_callers), fig_dix)
                ax.pie(values, labels=labels, autopct=make_autopct(values), startangle=90)
                ax.set_title(caller)
                fig_dix += 1

        plt.tight_layout()
        # plt.show()
        fig.savefig(f'{workdir}/{sample}/figures/{sample}.{platform}.nrst.match.minbpshift.pieplot.pdf')

def plot_pav_bpshift_density(matched_dir, sample, platform, aligner):

    all_match_bpshift = []
    for matched_file in os.listdir(matched_dir):
        if '.tsv' in matched_file:
            caller = matched_file.split('.')[1]
            df_matches = pd.read_csv(os.path.join(matched_dir, matched_file), sep='\t', header=[0])
            for idx, row in df_matches.iterrows():
                bpshift = int(row['BP_DIFF'])
                all_match_bpshift.append((bpshift, float(row['TARGET_SIZE_SIM']), caller))

    df_bpshifts = pd.DataFrame(all_match_bpshift, columns=['BP_DIFF', 'TARGET_SIZE_SIM', 'caller'])

    df_bpshifts1 = df_bpshifts[df_bpshifts['BP_DIFF'] <= 100]
    df_bpshifts2 = df_bpshifts[df_bpshifts['BP_DIFF'] <= 10]

    fig = plt.figure(figsize=(11, 4))
    fig.suptitle(f'{sample}-{platform}-{aligner}', fontsize=13)

    ax1 = fig.add_subplot(1, 2, 1)
    sns.kdeplot(data=df_bpshifts1, x='BP_DIFF', hue='caller', ax=ax1)
    ax2 = fig.add_subplot(1, 2, 2)
    sns.kdeplot(data=df_bpshifts2, x='BP_DIFF', hue='caller', ax=ax2)

    plt.tight_layout()
    # plt.show()

    fig.savefig(os.path.join(matched_dir, f'{sample}.{platform}.{aligner}.minbpshift.density.pdf'))

def pieplot_pav_bpshift_count_aligner(workdir, samples, aligner):
    long_read_callers = ['pbsv', 'svision', 'svim', 'cutesv', 'sniffles', 'nanovar']

    for sample in samples:
        platform_bpshift_dict = {}
        matched_dir = f'{workdir}/{sample}/intersect/nrst_match'

        for platform in ['hifi', 'ont']:
            bpshift_count_dict = {}
            for caller in long_read_callers:
                matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'

                shift_by_cutoff = {'0,10': 0, '10,20': 0, '20,50': 0, '50,100': 0, '>100': 0}
                if not os.path.exists(matched_file):
                    continue
                df_matches = pd.read_csv(matched_file, sep='\t', header=[0])
                # matched_calls[caller] = len(df_matches)

                for idx, row in df_matches.iterrows():
                    bpshift = int(row['BP_DIFF'])
                    if bpshift <= 10:
                        shift_by_cutoff['0,10'] += 1
                    elif bpshift <= 20:
                        shift_by_cutoff['10,20'] += 1
                    elif bpshift <= 50:
                        shift_by_cutoff['20,50'] += 1
                    elif bpshift <= 100:
                        shift_by_cutoff['50,100'] += 1
                    else:
                        shift_by_cutoff['>100'] += 1

                bpshift_count_dict[caller] = shift_by_cutoff

            platform_bpshift_dict[platform] = bpshift_count_dict

        fig = plt.figure(figsize=(20, 6))
        fig_dix = 1
        for platform, bpshift_dict in platform_bpshift_dict.items():
            for caller in long_read_callers:
                count_dict = bpshift_dict[caller]
                values = []
                labels = []
                for label, count in count_dict.items():
                    values.append(count)
                    labels.append(label)

                ax = fig.add_subplot(2, len(long_read_callers), fig_dix)
                ax.pie(values, labels=labels, autopct=make_autopct(values), startangle=90)
                ax.set_title(caller)
                fig_dix += 1

        plt.tight_layout()
        # plt.show()
        fig.savefig(f'{workdir}/{sample}/figures/{sample}.{aligner}.hifi-ont.nrst.match.minbpshift.pieplot.pdf')

def barplot_recall_precision(sample_dir, sample, aligner):
    benchmark_bed = f'{sample_dir}/pav_v112/pav_{sample}_CCS_svs.bed'
    df_benchmark = pd.read_csv(benchmark_bed, sep='\t')

    num_truth_calls = len(df_benchmark)

    long_read_callers = ['svision', 'pbsv', 'svim', 'cutesv', 'sniffles']

    rp_dict = {}
    for platform in ['ont', 'hifi']:
        this_platform_rp = []
        for caller in long_read_callers:
            caller_bed = f'{sample_dir}/{caller}/bed/{sample}.{platform}.{aligner}.{caller}.s5.filtered.bed'
            df_caller = pd.read_csv(caller_bed, sep='\t')
            num_caller_calls = len(df_caller)
            matched_file = f'{sample_dir}/intersect/nrst_match/{platform}_{aligner}/{sample}.{caller}.nrst.match.tsv'
            df_matches = pd.read_csv(matched_file, sep='\t', header=[0])
            num_matches = len(df_matches)

            recall = num_matches / num_truth_calls
            precision = num_matches / num_caller_calls
            this_platform_rp.append((recall, precision, fmeasure(precision, recall)))

        rp_dict[platform] = this_platform_rp

    recalls = {}
    precisions = {}
    fscores = {}

    for platform, rp_list in rp_dict.items():
        recalls[platform] = [val[0] for val in rp_list]
        precisions[platform] = [val[1] for val in rp_list]
        fscores[platform] = [val[2] for val in rp_list]
        # ax.plot(x_axis, fscores, label=platform, marker='o')

    width = 0.4
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(1,1,1)
    x_axis = np.arange(len(long_read_callers))

    ax.bar(x_axis, recalls['hifi'], width=width, read='center', edgecolor='white', label='HiFi')
    ax.bar(x_axis + width, recalls['ont'], width=width, read='center', edgecolor='white', label='ONT')
    ax.set_xticks(x_axis + width/2)
    ax.set_xticklabels(long_read_callers)

    ax.set_ylim([0, 1])
    ax.set_yticks(np.arange(0, 1.2, 0.2))
    ax.set_yticklabels([0, '20%', '40%', '60%', '80%','100%'])

    plt.legend()
    plt.tight_layout()
    plt.show()

def samples_hifi_recall_bpshift(workdir, samples):
    platform = 'hifi'

    long_read_callers = ['pbsv', 'svision', 'svim', 'cutesv']

    aligners = ['minimap2', 'pbmm2', 'ngmlr']

    fig = plt.figure(figsize=(18, 4))
    fig.suptitle(f'{platform} breakpoint shift', fontsize=14)

    gs = fig.add_gridspec(1, 3, hspace=0.13, wspace=0.05)
    axes = gs.subplots(sharex='col', sharey='row')

    shift_labels = ['0,10', '10,20', '20,50', '50,100', '>100']
    shift_label_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#a65628']
    aligner_ls = ['--', '-.', '-']

    all_legend = [Patch(facecolor='white', edgecolor='black', linestyle='--', label='minimap2'),
                  Patch(facecolor='white', edgecolor='black', linestyle='-.', label='pbmm2'),
                  Patch(facecolor='white', edgecolor='black', linestyle='-', label='ngmlr')]

    for idx in range(len(shift_labels)):
        all_legend.append(Patch(facecolor=shift_label_colors[idx], label=shift_labels[idx]))

    for sample_fig_idx in range(len(samples)):
        sample = samples[sample_fig_idx]
        print(f'Plot {sample} ...')
        matched_dir = f'{workdir}/{sample}/intersect/nrst_match'
        df_pavs = pd.read_csv(f'{workdir}/{sample}/pav_v112/pav_{sample}_CCS_svs.bed', sep='\t')
        pav_calls = len(df_pavs)

        aligners_shift = {}

        for aligner in aligners:

            shift_by_cutoff = {'0,10': np.zeros(len(long_read_callers)), '10,20': np.zeros(len(long_read_callers)),
                               '20,50': np.zeros(len(long_read_callers)), '50,100': np.zeros(len(long_read_callers)),
                               '>100': np.zeros(len(long_read_callers))}

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]

                matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'

                if not os.path.exists(matched_file):
                    continue
                df_matches = pd.read_csv(matched_file, sep='\t', header=[0])

                for idx, row in df_matches.iterrows():
                    bpshift = int(row['BP_DIFF'])
                    if bpshift <= 10:
                        shift_by_cutoff['0,10'][caller_idx] += 1
                    elif bpshift <= 20:
                        shift_by_cutoff['10,20'][caller_idx] += 1
                    elif bpshift <= 50:
                        shift_by_cutoff['20,50'][caller_idx] += 1
                    elif bpshift <= 100:
                        shift_by_cutoff['50,100'][caller_idx] += 1
                    else:
                        shift_by_cutoff['>100'][caller_idx] += 1

            aligners_shift[aligner] = pd.DataFrame(shift_by_cutoff)

        ax = axes[sample_fig_idx]
        bar_width = 0.25
        bar_gap = 0.02
        r1 = np.arange(len(long_read_callers))
        r2 = [x + bar_width + bar_gap for x in r1]
        r3 = [x + 2 * (bar_width + bar_gap) for x in r1]
        rs = [r1, r2, r3]

        for x_axis_idx in range(len(aligners)):
            aligner = aligners[x_axis_idx]
            df_bpshifts = aligners_shift[aligner]
            # totals = [i+j+k+m+n for i,j,k,m,n in zip(df_bpshifts['0,10'], df_bpshifts['10,20'], df_bpshifts['20,50'],df_bpshifts['50,100'],df_bpshifts['>100'])]
            totals = [pav_calls, pav_calls, pav_calls, pav_calls, pav_calls, pav_calls]
            for shift_label_idx in range(len(shift_labels)):
                label = shift_labels[shift_label_idx]
                this_label_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[label], totals)]
                if shift_label_idx == 0:
                    ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, edgecolor='black', linestyle=aligner_ls[x_axis_idx], width=bar_width, color=shift_label_colors[shift_label_idx])
                else:
                    bottom_sum = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[0]], totals)]
                    for bottom_idx in range(1, shift_label_idx):
                        bottom_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[bottom_idx]], totals)]
                        for bottom_val_idx in range(len(bottom_sum)):
                            bottom_sum[bottom_val_idx] += bottom_pcrt[bottom_val_idx]
                    ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, edgecolor='black', linestyle=aligner_ls[x_axis_idx], width=bar_width, color=shift_label_colors[shift_label_idx])

        # x_axis_ticks = [x + (bar_width + bar_gap) / 2 for x in r1]
        ax.set_title(sample)
        if sample_fig_idx == 0:
            ax.set_ylabel('Percent')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 6))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 6)])
        ax.set_xticks(r2)
        ax.set_xticklabels(long_read_callers)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(color='grey', linestyle=(0, (5, 5)), lw=1, axis='y')

        if sample_fig_idx == 2:
            ax.legend(handles=all_legend, loc='upper left', bbox_to_anchor=(1,1), ncol=1)

    gs.tight_layout(fig)
    # plt.show()
    fig.savefig(f'{workdir}/figures/{platform}.samples.minimap2-pbmm2-ngmlr.recall.bpshift.barplot.pdf')


    
def barplot_bpshift_count(matched_dir, sample, aligner):
    long_read_callers = ['pbsv', 'svision', 'svim', 'cutesv', 'sniffles']
    bpshift_dict = {}

    for caller in long_read_callers:
        caller_bpshifts = {}
        for platform in ['ont', 'hifi']:
            matched_file = f'{matched_dir}/{platform}_{aligner}/{sample}.{caller}.nrst.match.tsv'

            shift_by_cutoff = {'0,10': 0, '10,20': 0, '20,50': 0, '50,100': 0, '>100': 0}

            df_matches = pd.read_csv(matched_file, sep='\t', header=[0])
            for idx, row in df_matches.iterrows():
                bpshift = int(row['BP_DIFF'])
                if bpshift <= 10:
                    shift_by_cutoff['0,10'] += 1
                elif 10 < bpshift <= 20:
                    shift_by_cutoff['10,20'] += 1
                elif 20 < bpshift <= 50:
                    shift_by_cutoff['20,50'] += 1
                elif 50 < bpshift <= 100:
                    shift_by_cutoff['50,100'] += 1
                else:
                    shift_by_cutoff['>100'] += 1

            caller_bpshifts[platform] = shift_by_cutoff

        bpshift_dict[caller] = caller_bpshifts

    fig = plt.figure(figsize=(19, 4))
    fig.suptitle(f'{sample}-{aligner}', fontsize=13)

    i = 1
    for caller, shifts in bpshift_dict.items():
        y_val_pcrt = {}
        y_val = {}
        x_val = []
        for aligner, count_dict in shifts.items():
            x_val = list(count_dict.keys())
            y_val[aligner] = list(count_dict.values())
            y_sum = sum(list(count_dict.values()))
            y_val_pcrt[aligner] = [(ele / y_sum) * 100 for ele in list(count_dict.values())]

        width = 0.4

        x_axis = np.arange(0, len(x_val))

        ax = fig.add_subplot(1, len(bpshift_dict), i)

        ax.set_ylim([0, 100])
        ax.set_yticks(np.linspace(0, 100, 6))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 100, 6)])

        ax.set_ylabel("Percentage")

        ax.bar(x_axis, y_val_pcrt['minimap2'], width=width, read='center', edgecolor='white', label='minimap2')
        ax.bar(x_axis + width, y_val_pcrt['ngmlr'], width=width, read='center', edgecolor='white',label='ngmlr')

        ax.set_title(caller, fontsize=14)
        ax.set_xticks(x_axis + width/2)
        ax.set_xticklabels(x_val)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # for j in range(len(x_axis)):
        #     ax.text(x=x_axis[j], y=y_val_pcrt['minimap2'][j], ha='center', va='bottom', s=y_val['minimap2'][j])
        #     ax.text(x=x_axis[j] + width, y=y_val_pcrt['ngmlr'][j], ha='center', va='bottom', s=y_val['ngmlr'][j])

        i += 1

    plt.legend()
    plt.tight_layout()
    plt.show()
    # fig.savefig(os.path.join(matched_dir, f'{sample}.{platform}.minbpshift.count.pdf'))
    
def samples_recall_bpshift_of_platform(workdir, samples, platform):
    long_read_callers = ['pbsv', 'svision', 'svim', 'cutesv', 'sniffles', 'nanovar']

    aligners = ['minimap2', 'ngmlr']

    fig = plt.figure(figsize=(18, 4))
    fig.suptitle(f'{platform} breakpoint shift', fontsize=14)

    gs = fig.add_gridspec(1, 3, hspace=0.13, wspace=0.05)
    axes = gs.subplots(sharex='col', sharey='row')

    shift_labels = ['0,10', '10,20', '20,50', '50,100', '>100']
    shift_label_colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#a65628']

    all_legend = [Patch(facecolor='white', edgecolor='black', linestyle='--', label='minimap2'), Patch(facecolor='white', edgecolor='black', label='ngmlr')]

    for idx in range(len(shift_labels)):
        all_legend.append(Patch(facecolor=shift_label_colors[idx], label=shift_labels[idx]))

    for sample_fig_idx in range(len(samples)):
        sample = samples[sample_fig_idx]
        print(f'Plot {sample} ...')
        matched_dir = f'{workdir}/{sample}/intersect/nrst_match'
        df_pavs = pd.read_csv(f'{workdir}/{sample}/pav_v112/pav_{sample}_CCS_svs.bed', sep='\t')
        pav_calls = len(df_pavs)

        aligners_shift = {}

        for aligner in aligners:

            shift_by_cutoff = {'0,10': np.zeros(len(long_read_callers)), '10,20': np.zeros(len(long_read_callers)),
                               '20,50': np.zeros(len(long_read_callers)), '50,100': np.zeros(len(long_read_callers)),
                               '>100': np.zeros(len(long_read_callers))}

            for caller_idx in range(len(long_read_callers)):
                caller = long_read_callers[caller_idx]

                matched_file = f'{matched_dir}/{sample}.{caller}.{platform}.com-{aligner}.nrst.match.tsv'

                if not os.path.exists(matched_file):
                    continue
                df_matches = pd.read_csv(matched_file, sep='\t', header=[0])

                for idx, row in df_matches.iterrows():
                    bpshift = int(row['BP_DIFF'])
                    if bpshift <= 10:
                        shift_by_cutoff['0,10'][caller_idx] += 1
                    elif bpshift <= 20:
                        shift_by_cutoff['10,20'][caller_idx] += 1
                    elif bpshift <= 50:
                        shift_by_cutoff['20,50'][caller_idx] += 1
                    elif bpshift <= 100:
                        shift_by_cutoff['50,100'][caller_idx] += 1
                    else:
                        shift_by_cutoff['>100'][caller_idx] += 1

            aligners_shift[aligner] = pd.DataFrame(shift_by_cutoff)

        ax = axes[sample_fig_idx]
        bar_width = 0.4
        bar_gap = 0.05
        r1 = np.arange(len(long_read_callers))
        r2 = [x + bar_width + bar_gap  for x in r1]
        rs = [r1, r2]
        for x_axis_idx in range(len(aligners)):
            aligner = aligners[x_axis_idx]
            df_bpshifts = aligners_shift[aligner]
            # totals = [i+j+k+m+n for i,j,k,m,n in zip(df_bpshifts['0,10'], df_bpshifts['10,20'], df_bpshifts['20,50'],df_bpshifts['50,100'],df_bpshifts['>100'])]
            totals = [pav_calls, pav_calls, pav_calls, pav_calls, pav_calls, pav_calls]
            for shift_label_idx in range(len(shift_labels)):
                label = shift_labels[shift_label_idx]
                this_label_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[label], totals)]
                if shift_label_idx == 0:
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, edgecolor='black', linestyle="--", width=bar_width, color=shift_label_colors[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, label=label, edgecolor='black', width=bar_width, color=shift_label_colors[shift_label_idx])
                else:
                    bottom_sum = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[0]], totals)]
                    for bottom_idx in range(1, shift_label_idx):
                        bottom_pcrt = [i / j * 100 for i, j in zip(df_bpshifts[shift_labels[bottom_idx]], totals)]
                        for bottom_val_idx in range(len(bottom_sum)):
                            bottom_sum[bottom_val_idx] += bottom_pcrt[bottom_val_idx]
                    if aligner == 'minimap2':
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, edgecolor='black', linestyle="--", width=bar_width, color=shift_label_colors[shift_label_idx])
                    else:
                        ax.bar(rs[x_axis_idx], this_label_pcrt, bottom=bottom_sum, label=label, edgecolor='black', width=bar_width, color=shift_label_colors[shift_label_idx])

        x_axis_ticks = [x + (bar_width + bar_gap) / 2 for x in r1]
        ax.set_title(sample)
        if sample_fig_idx == 0:
            ax.set_ylabel('Percent')
        ax.set_ylim(0, 100)
        ax.set_yticks(np.linspace(0, 100, 6))
        ax.set_yticklabels([f'{int(val)}%' for val in np.linspace(0, 100, 6)])
        ax.set_xticks(x_axis_ticks)
        ax.set_xticklabels(long_read_callers)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(color='grey', linestyle=(0, (5, 5)), lw=1, axis='y')

        if sample_fig_idx == 2:
            ax.legend(handles=all_legend, loc='upper left', bbox_to_anchor=(1,1), ncol=1)

    gs.tight_layout(fig)
    # plt.show()
    fig.savefig(f'{workdir}/figures/{platform}.samples.recall.bpshift.barplot.pdf')
    
'''