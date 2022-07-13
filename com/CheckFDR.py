#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2021/10/25

'''
import pandas as pd
import pysam
from statistics import mean

from helpers.Annot import *
from helpers.Functions import *
from helpers.Reader import *



def merge_assm_uniques(workdir, datasets):

    simple_reps = pysam.Tabixfile(iMACSIMPLE, 'r')
    rmsk = pysam.Tabixfile(iMACRMSK, 'r')
    sds = pysam.Tabixfile(iMACSD, 'r')

    plat_assembler = {'HiFi': 'hifiasm', 'ONT': 'flye'}

    for dataset in datasets:
        plat = 'HiFi'
        if 'ont' in dataset:
            plat = 'ONT'

        output_dir = f'{workdir}/{plat}/minimap2_{dataset}/filtered/comstra'

        assm_tmp_file = open(f'{output_dir}/asm_tmp.txt', 'w')
        read_tmp_file = open(f'{output_dir}/read_tmp.txt', 'w')

        for asm_caller in ASMCALLERS:
            for read_caller in CALLERS:
                merged_vcf = f'{output_dir}/{read_caller}-minimap2.{asm_caller}-minimap2-{plat_assembler[plat]}.{dataset}.jasmine.merged.vcf'
                asm_unique_vcf = f'{output_dir}/{read_caller}-minimap2.{asm_caller}-minimap2-{plat_assembler[plat]}.{dataset}.{asm_caller}.unique.vcf'
                read_unique_vcf = f'{output_dir}/{read_caller}-minimap2.{asm_caller}-minimap2-{plat_assembler[plat]}.{dataset}.{read_caller}.unique.vcf'

                print(f'{dataset}-{read_caller}-{asm_caller}')

                separate_stra_merged_vcfs(merged_vcf, asm_unique_vcf, read_unique_vcf)

                print(asm_unique_vcf, file=assm_tmp_file)
                print(read_unique_vcf, file=read_tmp_file)

        asm_unique_merged = f'{output_dir}/assembly_based.uniques.jasmine.merged.vcf'
        read_unique_merged = f'{output_dir}/read_based.uniques.jasmine.merged.vcf'

        assm_tmp_file.close()
        read_tmp_file.close()

        # cmd1 = f'{JASMINE} file_list={output_dir}/asm_tmp.txt out_file={asm_unique_merged} max_dist=1000 spec_len=50 spec_reads=1'
        # os.system(cmd1)
        # cmd2 = f'{JASMINE} file_list={output_dir}/read_tmp.txt out_file={read_unique_merged} max_dist=1000 spec_len=50 spec_reads=1'
        # os.system(cmd2)

        # asm_uniques = open(f'{output_dir}/assembly_based.uniques.annot.tsv', 'w')
        #
        # for rec in vcf.Reader(open(asm_unique_merged, 'r')):
        #     chrom, start = rec.CHROM, int(rec.POS)
        #     end, svlen = int(rec.INFO['END']), int(rec.INFO['SVLEN'][0])
        #     svtype = rec.INFO['SVTYPE']
        #     supp_vec = rec.INFO['SUPP_VEC']
        #     supp = rec.INFO['SUPP']
        #
        #     if start == end:
        #         end += svlen
        #
        #     region_label, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simple_reps, rmsk, sds)
        #     print(f'{chrom}\t{start}\t{rec.ID}\t{svlen}\t{svtype}\t{supp_vec}\t{supp}\t{region_label}\t{rptype}\t{round(pcrt, 2)}', file=asm_uniques)
        #
        # read_uniques = open(f'{output_dir}/read_based.uniques.annot.tsv', 'w')
        # for rec in vcf.Reader(open(read_unique_merged, 'r')):
        #     chrom, start = rec.CHROM, int(rec.POS)
        #     end, svlen = int(rec.INFO['END']), int(rec.INFO['SVLEN'][0])
        #     svtype = rec.INFO['SVTYPE']
        #     supp_vec = rec.INFO['SUPP_VEC']
        #     supp = rec.INFO['SUPP']
        #
        #     if start == end:
        #         end += svlen
        #
        #     region_label, rptype, pcrt = annotate_sv_region(chrom, start, end, 0, simple_reps, rmsk, sds)
        #     print(f'{chrom}\t{start}\t{rec.ID}\t{svlen}\t{svtype}\t{supp_vec}\t{supp}\t{region_label}\t{rptype}\t{round(pcrt, 2)}', file=read_uniques)

    # os.remove(f'{output_dir}/*.unique.vcf')
    # os.remove(f'{output_dir}/*.txt')

def annotate_uniques(workdir, dataset, bam_path, flank):

    columns = ['chrom', 'start', 'id', 'svlen', 'svtype', 'suppvec', 'supp', 'region', 'reptype', 'pcrt']

    bam_file = pysam.AlignmentFile(bam_path, 'rb')

    plat = 'HiFi'
    if 'ont' in dataset:
        plat = 'ONT'

    for stra in ['read', 'assembly']:
        unique_file = f'{workdir}/{plat}/minimap2_{dataset}/filtered/comstra/{stra}_based.uniques.annot.tsv'
        classified_out = open(f'{workdir}/{plat}/minimap2_{dataset}/filtered/comstra/{stra}_based.uniques.annot.sigs.tsv', 'w')

        df_unique = pd.read_csv(unique_file, sep='\t', names=columns)

        print(f'Process {len(df_unique)} unique {stra} calls')
        counter = len(df_unique)
        for idx, row in df_unique.iterrows():
            chrom, start, svtype, svlen, svregion, rptype, pcrt = str(row['chrom']), int(row['start']), \
                                                                       row['svtype'], abs(int(row['svlen'])), row['region'], row['reptype'], float(row['pcrt'])

            suppvec, supp = row['suppvec'], row['supp']
            end = start + svlen

            mean_mapq, num_signature_reads, nearest_sigs = check_sv_signature_in(chrom, start, end, bam_file, 50, 100000, flank)

            print(f'{chrom}\t{start}\t{end}\t{svtype}\t{svlen}\t{svregion}\t{rptype}\t{pcrt}\t{suppvec}\t{supp}\t{mean_mapq}\t{num_signature_reads}\t{nearest_sigs}', file=classified_out)
            counter -= 1

            if counter % 1000 == 0:
                print(f'{counter} remaining ...')



def check_sv_signature_in(chrom, start, end, aln_file, min_sv_size, max_sv_size, flank):

    aligns = aln_file.fetch(chrom, start - flank, end + flank)
    mapping_quality = []
    qnames = []
    raw_signatures = []
    while True:
        try:
            current_alignment = next(aligns)
            qname = current_alignment.query_name

            qnames.append(qname)
            mapping_quality.append(current_alignment.mapping_quality)

            if current_alignment.is_unmapped or current_alignment.is_secondary or current_alignment.mapping_quality < 20:
                continue
            if current_alignment.is_supplementary:
                sigs = analyze_alignment_indel(current_alignment, aln_file, current_alignment.query_name,)
                raw_signatures.extend(sigs)
            else:

                supplementary_alignments = retrieve_other_alignments(current_alignment, aln_file)
                good_suppl_alns = [aln for aln in supplementary_alignments if
                                   not aln.is_unmapped and aln.mapping_quality >= 20]
                sigs = analyze_alignment_indel(current_alignment, aln_file, current_alignment.query_name)
                raw_signatures.extend(sigs)

                sigs = analyze_read_segments(current_alignment, good_suppl_alns, aln_file, min_sv_size, max_sv_size)
                raw_signatures.extend(sigs)

        except StopIteration:
            break

    if len(mapping_quality) == 0:
        ## Cannot find mapped reads in this region
        mean_mapq = -1
    else:
        mean_mapq = mean(mapping_quality)


    sv_sig = []
    nearest_sv_sig = []

    for sig in raw_signatures:
        sig_start, sig_end = sig[1], sig[2]
        sig_size = sig_end - sig_start + 1

        sig_bpshift = min_breakpoint_shift(start, end, sig_start, sig_end)
        size_sim = size_similarity(sig_size, end - start + 1)

        if sig_bpshift < 500:
            sv_sig.append(sig)

        if sig_bpshift < 500 and size_sim >= 0.5:
            nearest_sv_sig.append(sig)

    return mean_mapq, len(sv_sig), len(nearest_sv_sig)


def analyze_alignment_indel(alignment, bam, query_name):
    sv_signatures = []
    ref_chr = bam.getrname(alignment.reference_id)
    ref_start = alignment.reference_start
    indels = analyze_cigar_indel(alignment.cigartuples, 50)
    for pos_ref, pos_read, length, typ in indels:
        if typ == "DEL":
            sv_signatures.append((ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, "cigar", query_name, 'DEL'))
        elif typ == "INS":
            sv_signatures.append((ref_chr, ref_start + pos_ref, ref_start + pos_ref + length, "cigar", query_name, 'INS'))
    return sv_signatures

def analyze_cigar_indel(tuples, min_length):
    """Parses CIGAR tuples (op, len) and returns Indels with a length > minLength"""
    pos_ref = 0
    pos_read = 0
    indels = []
    for operation, length in tuples:
        if operation == 0:                     # alignment match
            pos_ref += length
            pos_read += length
        elif operation == 1:                   # insertion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "INS"))
            pos_read += length
        elif operation == 2:                   # deletion
            if length >= min_length:
                indels.append((pos_ref, pos_read, length, "DEL"))
            pos_ref += length
        elif operation == 4:                   # soft clip
            pos_read += length
        elif operation == 7 or operation == 8:        # match or mismatch
            pos_ref += length
            pos_read += length

    return indels

def analyze_read_segments(primary, supplementaries, bam, min_sv_size, max_sv_size, segment_overlap=50, segment_gap=50):
    read_name = primary.query_name

    alignments = [primary] + supplementaries
    alignment_list = []
    for alignment in alignments:
        #correct query coordinates for reversely mapped reads
        if alignment.is_reverse:
            q_start = alignment.infer_read_length() - alignment.query_alignment_end
            q_end = alignment.infer_read_length() - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        new_alignment_dict = {  'q_start': q_start,
                                'q_end': q_end,
                                'ref_id': alignment.reference_id,
                                'ref_start': alignment.reference_start,
                                'ref_end': alignment.reference_end,
                                'is_reverse': alignment.is_reverse  }

        alignment_list.append(new_alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))
    #inferred_read_length = alignments[0].infer_read_length()

    sv_signatures = []
    tandem_duplications = []

    for index in range(len(sorted_alignment_list) - 1):
        alignment_current = sorted_alignment_list[index]
        alignment_next = sorted_alignment_list[index + 1]


        distance_on_read = alignment_next['q_start'] - alignment_current['q_end']
        #Same chromosome
        if alignment_current['ref_id'] == alignment_next['ref_id']:
            ref_chr = bam.getrname(alignment_current['ref_id'])
            #Same orientation
            if alignment_current['is_reverse'] == alignment_next['is_reverse']:
                #Compute distance on reference depending on orientation
                if alignment_current['is_reverse']:
                    distance_on_reference = alignment_current['ref_start'] - alignment_next['ref_end']
                else:
                    distance_on_reference = alignment_next['ref_start'] - alignment_current['ref_end']
                #No overlap on read
                if distance_on_read >= -segment_overlap:
                    #No overlap on reference
                    if distance_on_reference >= -segment_overlap:
                        deviation = distance_on_read - distance_on_reference
                        #INS candidate
                        if deviation >= min_sv_size:
                            #No gap on reference
                            if not alignment_current['is_reverse']:
                                sv_signatures.append((ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] + deviation, "suppl", read_name, 'INS'))
                            else:
                                sv_signatures.append((ref_chr, alignment_current['ref_start'], alignment_current['ref_start'] + deviation, "suppl", read_name, 'INS'))
                        #DEL candidate
                        elif -max_sv_size <= deviation <= -min_sv_size:
                            #No gap on read
                            if not alignment_current['is_reverse']:
                                sv_signatures.append((ref_chr, alignment_current['ref_end'], alignment_current['ref_end'] - deviation, "suppl", read_name, 'DEL'))
                            else:
                                sv_signatures.append((ref_chr, alignment_next['ref_end'], alignment_next['ref_end'] - deviation, "suppl", read_name, 'DEL'))

                    #overlap on reference
                    else:
                        #Tandem Duplication
                        if distance_on_reference <= -min_sv_size:
                            if not alignment_current['is_reverse']:
                                #Tandem Duplication
                                if alignment_next['ref_end'] > alignment_current['ref_start']:
                                    tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_current['ref_end'], True, True))
                                #Large tandem duplication
                                elif distance_on_reference >= -max_sv_size:
                                    tandem_duplications.append((ref_chr, alignment_next['ref_start'], alignment_current['ref_end'], False, True))
                            else:
                                #Tandem Duplication
                                if alignment_next['ref_start'] < alignment_current['ref_end']:
                                    tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_next['ref_end'], True, False))
                                #Large tandem duplication
                                elif distance_on_reference >= -max_sv_size:
                                    tandem_duplications.append((ref_chr, alignment_current['ref_start'], alignment_next['ref_end'], False, False))
            #Different orientations
            else:
                #Normal to reverse
                if not alignment_current['is_reverse'] and alignment_next['is_reverse']:
                    if alignment_next['ref_start'] - alignment_current['ref_end'] >= -segment_overlap: # Case 1
                        #INV candidate
                        if min_sv_size <= alignment_next['ref_end'] - alignment_current['ref_end'] <= max_sv_size:
                            sv_signatures.append((ref_chr, alignment_current['ref_end'], alignment_next['ref_end'], "suppl", read_name, "INV"))

                    elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -segment_overlap: # Case 3
                        #INV candidate
                        if min_sv_size <= alignment_current['ref_end'] - alignment_next['ref_end'] <= max_sv_size:
                            sv_signatures.append((ref_chr, alignment_next['ref_end'], alignment_current['ref_end'], "suppl", read_name, "INV"))

                #Reverse to normal
                if alignment_current['is_reverse'] and not alignment_next['is_reverse']:
                    if alignment_next['ref_start'] - alignment_current['ref_end'] >= -segment_overlap: # Case 2
                        #INV candidate
                        if min_sv_size <= alignment_next['ref_start'] - alignment_current['ref_start'] <= max_sv_size:
                            sv_signatures.append((ref_chr, alignment_current['ref_start'], alignment_next['ref_start'], "suppl", read_name, "INV"))
                    elif alignment_current['ref_start'] - alignment_next['ref_end'] >= -segment_overlap: # Case 4
                        #INV candidate
                        if min_sv_size <= alignment_current['ref_start'] - alignment_next['ref_start'] <= max_sv_size:
                            sv_signatures.append((ref_chr, alignment_next['ref_start'], alignment_current['ref_start'], "suppl", read_name, "INV"))

        # Handle tandem duplications
    current_chromosome = None
    current_starts = []
    current_ends = []
    current_copy_number = 0
    current_fully_covered = []
    for tandem_duplication in tandem_duplications:
        if current_chromosome == None:
            current_chromosome = tandem_duplication[0]
            current_starts.append(tandem_duplication[1])
            current_ends.append(tandem_duplication[2])
            current_copy_number = 1
            current_fully_covered.append(tandem_duplication[3])
            current_direction = tandem_duplication[4]
        else:
            if is_similar(current_chromosome, mean(current_starts), mean(current_ends), tandem_duplication[0],
                          tandem_duplication[1], tandem_duplication[2]) and current_direction == tandem_duplication[4]:
                current_starts.append(tandem_duplication[1])
                current_ends.append(tandem_duplication[2])
                current_copy_number += 1
                current_fully_covered.append(tandem_duplication[3])
            else:
                fully_covered = True if sum(current_fully_covered) else False
                sv_signatures.append((current_chromosome, int(mean(current_starts)), int(mean(current_ends)),
                                               current_copy_number, fully_covered, "suppl", read_name))
                current_chromosome = tandem_duplication[0]
                current_starts = [tandem_duplication[1]]
                current_ends = [tandem_duplication[2]]
                current_copy_number = 1
                current_fully_covered = [tandem_duplication[3]]
    if current_chromosome != None:
        sv_signatures.append((current_chromosome, int(mean(current_starts)), int(mean(current_ends)), "suppl", read_name, 'DUP'))


    return sv_signatures


def is_similar(chr1, start1, end1, chr2, start2, end2):
    if chr1 == chr2 and abs(start1 - start2) < 20 and abs(end1 - end2) < 20:
        return True
    else:
        return False

def retrieve_other_alignments(main_alignment, bam):
    """Reconstruct other alignments of the same read for a given alignment from the SA tag"""
    #reconstructing other alignments from SA tag does not work if sequence of main_alignment is hard-clipped
    if main_alignment.get_cigar_stats()[0][5] > 0:
        return []
    try:
        sa_tag = main_alignment.get_tag("SA").split(";")
    except KeyError:
        return []
    other_alignments = []
    # For each other alignment encoded in the SA tag
    for element in sa_tag:
        # Read information from the tag
        fields = element.split(",")
        if len(fields) != 6:
            continue
        rname = fields[0]
        pos = int(fields[1])
        strand = fields[2]
        # CIGAR string encoded in SA tag is shortened
        cigar = fields[3]
        mapq = int(fields[4])
        nm = int(fields[5])

        # Generate an aligned segment from the information
        a = pysam.AlignedSegment()
        a.query_name = main_alignment.query_name
        a.query_sequence= main_alignment.query_sequence
        if strand == "+":
            a.flag = 2048
        else:
            a.flag = 2064
        a.reference_id = bam.get_tid(rname)
        a.reference_start = pos - 1
        try:
            a.mapping_quality = mapq
        except OverflowError:
            a.mapping_quality = 0
        a.cigarstring = cigar
        a.next_reference_id = -1
        a.next_reference_start = -1
        a.template_length = 0
        a.query_qualities = main_alignment.query_qualities
        a.set_tags([("NM", nm, "i")])

        other_alignments.append(a)
    return other_alignments

# if __name__ == '__main__':
#     low_cov_datasets = ['hifi_10kb', 'hifi_11kb']
#     compare_merged_insdel_to_pav(iMACDIR, low_cov_datasets, 'giab_ctg', iMACBAM, iMACSIMPLE, iMACRMSK)