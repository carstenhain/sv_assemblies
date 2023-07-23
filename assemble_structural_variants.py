import pysam
import read_methods
import os
import numpy as np
import mappy as mp
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from Bio import SeqIO
import warnings
import vcf

warnings.filterwarnings("ignore")

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def write_generic_vcf_header (fo, generic_vcf_path):
    for line in open(generic_vcf_path, "r"):
        fo.write(line)
    fo.flush()

def get_sv_reads_from_segment (segments, samfile, external_blacklist, breakpoint_distance_tolerance):

    # array with reads of interest
    sv_reads = []

    # list of reads (=read names) already in sv_reads to not add duplicates
    blacklist = []

    # query all segments
    for idx, segment in enumerate(segments):
        chrom = segment[0]
        start = segment[1]
        end = segment[2]

        small_internal_segment = False
        if end - start < 500:
            small_internal_segment = True

        # all reads spanning this segment
        for read in samfile.fetch(chrom, start, end):

            # only reads with supplementary alignment that are not yet in blacklist (=sv_reads)
            if read.has_tag("SA") and not (read.query_name in blacklist) and not (read.query_name in external_blacklist):

                write_read = False
                # add read to sv_reads if at least one end is close to the segment border
                if np.abs(read.reference_start - start) < breakpoint_distance_tolerance:
                    write_read = True
                if np.abs(read.reference_end - end) < breakpoint_distance_tolerance:
                    write_read = True

                if small_internal_segment:
                    write_read = False
                    if np.abs(read.reference_start - start) < breakpoint_distance_tolerance and np.abs(read.reference_end - end) < breakpoint_distance_tolerance:
                        write_read = True

                if write_read:
                    sv_reads.append(read)
                    blacklist.append(read.query_name)

    return sv_reads

def get_assembly_mappings_mappy (alignments, query_length, min_length, min_mapq):

    mappings = []

    for aln in alignments:

        direction = "+"
        if aln.strand < 0:
            direction = "-"

        if aln.mlen > min_length and aln.mapq > min_mapq:

            continue_with_mapping = True

            if aln.mapq == 1 and aln.mlen < 5000:
                continue_with_mapping = False

            if continue_with_mapping:
                if direction == "-":

                    mappings.append(
                        [
                            aln.q_st,
                            aln.q_en,
                            direction,
                            aln.ctg,
                            aln.r_st,
                            aln.r_en,
                            aln.mapq,
                            aln.cigar
                        ])
                else:

                    mappings.append(
                        [
                            aln.q_st,
                            aln.q_en,
                            direction,
                            aln.ctg,
                            aln.r_st,
                            aln.r_en,
                            aln.mapq,
                            aln.cigar
                        ])

    sorted_mappings = sorted(mappings)
    """
    0   query_start
    1   query_end
    2   direction
    3   ref_name
    4   ref_start
    5   ref_end
    6   mapq
    7   cigar
    """
    return sorted_mappings

def get_adjacent_mappings (sorted_mappings, chrom, start, end):
    for i in range(0, len(sorted_mappings)):
        if sorted_mappings[i][3] == chrom and np.abs(start - sorted_mappings[i][4]) < 10 and np.abs(end - sorted_mappings[i][5]) < 10:
            adjacent_mappings = []
            if i > 0:
                adjacent_mappings.append(i - 1)
            if i < len(sorted_mappings) - 1:
                adjacent_mappings.append(i + 1)
            return adjacent_mappings

# removes terminal streched with coverage lower than 2 from the assembly
# uses bedtools, samtools and mappy
# return True if terminal bases are pruned
def prune_final_assembly (assembly_name):

    assembly_sequence = list(SeqIO.parse("assemblies/" + assembly_name + ".fasta", "fasta"))[0].seq

    #print("\tPruning assembly assemblies/" + assembly_name + ".fasta with length " + str(len(assembly_sequence)))

    # map reads to assembly and write to bed
    tmp_reads_bed = open("tmp.prune.reads.bed", "w")
    assembly_aligner = mp.Aligner("assemblies/" + assembly_name + ".fasta", preset="map-ont")

    if not assembly_aligner: raise Exception("ERROR: failed to load/build index")
    for read in SeqIO.parse("reads/" + assembly_name + ".fastq", "fastq"):
        for aln in list(assembly_aligner.map(str(read.seq))):
            if aln.mapq > 10 and (aln.r_en - aln.r_st) > 100:
                tmp_reads_bed.write("ctg\t" + str(aln.r_st) + "\t" + str(aln.r_en) + "\n")

    tmp_reads_bed.close()

    # build single base bed file
    tmp_per_base_bed = open("tmp.prune.perbase.bed", "w")

    # build per base bed
    for i in range(0, len(assembly_sequence)):
        tmp_per_base_bed.write("ctg\t" + str(i) + "\t" + str(1 + i) + "\n")
    tmp_per_base_bed.close()

    # count coverage using bedtools
    os.system("bedtools intersect -c -a tmp.prune.perbase.bed -b tmp.prune.reads.bed > tmp.prune.perbase.count.bed")

    # load coverage data
    coverage_data = []
    for line in open("tmp.prune.perbase.count.bed"):
        coverage_data.append([int(line.split("\t")[1]), int(line.split("\t")[-1])])

    # prune start
    prune_start_position = 0
    for i in range(0, len(coverage_data)):
        if coverage_data[i][1] < 2:
            prune_start_position += 1
        else:
            break

    # prune end
    prune_end_position = len(assembly_sequence)
    for i in range(0, len(coverage_data)):
        if coverage_data[len(coverage_data) - i - 1][1] < 2:
            prune_end_position = len(coverage_data) - i
        else:
            break

    if np.any([prune_start_position != 0, prune_end_position != len(assembly_sequence)]):

        print(f"{bcolors.WARNING}" + "\t\t\tPruning assembly" + f"{bcolors.ENDC}")
        if prune_start_position != 0:
            print(f"{bcolors.WARNING}" + "\t\t\t\tPruning " + str(prune_start_position) + " bp from start" + f"{bcolors.ENDC}")
        if prune_end_position != len(assembly_sequence):
            print(f"{bcolors.WARNING}" + "\t\t\t\tPruning " + str(len(assembly_sequence) - prune_end_position) + " bp from end" + f"{bcolors.ENDC}")

        # extract and save unpruned sequence using bedtools getfasta

        tmp_pruning_bed = open("pruning.bed", "w")
        tmp_pruning_bed.write("\t".join(
            [
                assembly_name,
                str(prune_start_position),
                str(prune_end_position)
            ]) + "\n")
        tmp_pruning_bed.close()

        os.system("samtools faidx assemblies/" + assembly_name + ".fasta")
        os.system("bedtools getfasta -fi assemblies/" + assembly_name + ".fasta -bed pruning.bed > assemblies/" + assembly_name + ".pruned.fasta")
        os.system("mv assemblies/" + assembly_name + ".pruned.fasta assemblies/" + assembly_name + ".fasta")

        return True

    return False

def plot_query_mappings (sorted_mappings, axis, query_length, segments):

    counter = 0

    axis.set_xlim(0, query_length)

    for sm in sorted_mappings:

        counter += 0.01

        direction = sm[2]
        head_length = np.amin([1000, np.abs(sm[1] - sm[0])])

        col = "blue"
        for seg in segments:
            if seg[0] == sm[3] and np.abs(seg[1] - sm[4]) < 10 and np.abs(seg[2] - sm[5]) < 10:
                col = "red"

        if direction == "+":

            axis.arrow(x=sm[0],
                       y=counter,
                       dx=sm[1] - sm[0],
                       dy=0,
                       head_length=head_length,
                       length_includes_head=True,
                       color=col)
            axis.text(sm[0] + 0.5 * np.abs(sm[1] - sm[0]),
                      counter + 0.002,
                      sm[3] + ":" + str(sm[4]) + "-" + str(sm[5]))
        else:

            axis.arrow(x=sm[1],
                       y=counter,
                       dx=sm[0] - sm[1],
                       dy=0,
                       head_length=head_length,
                       length_includes_head=True,
                       color=col)
            axis.text(sm[0] + 0.5 * np.abs(sm[1] - sm[0]),
                      counter + 0.002,
                      sm[3] + ":" + str(sm[4]) + "-" + str(sm[5]))

        axis.set_ylim(0, counter + 0.02)

def plot_qc_plot (assembly_name, sorted_mappings, ax_reads, ax_cov):

    try:

        assembly_sequence = list(SeqIO.parse("assemblies/" + assembly_name + ".fasta", "fasta"))[0].seq

        # build bed files
        # tmp.bed reads
        tmp_bed = open("tmp.bed", "w")
        # tmp.perbase.bed coverage per base
        tmp_per_base_bed = open("tmp.perbase.bed", "w")
        # tmp.perregion.bed coverage per region
        tmp_per_region_bed = open("tmp.perregion.bed", "w")

        # build per base bed
        for i in range(0, len(assembly_sequence)):
            tmp_per_base_bed.write("ctg\t" + str(i) + "\t" + str(1 + i) + "\n")
        tmp_per_base_bed.close()

        """
        # build per region bed
        for i in range(0, len(sorted_mappings)):
            tmp_per_region_bed.write("ctg\t" + str(sorted_mappings[i][0]) + "\t" + str(sorted_mappings[i][1]) + "\n")
        tmp_per_region_bed.close()
        """

        # axes
        ax_reads.set_xlim(0, len(assembly_sequence))
        ax_reads.set_ylabel("Reads")

        ax_cov.set_xlim(0, len(assembly_sequence))
        ax_cov.set_ylabel("Coverage")

        # map reads to assembly
        assembly_aligner = mp.Aligner("assemblies/" + assembly_name + ".fasta", preset="map-ont")
        if not assembly_aligner: raise Exception("ERROR: failed to load/build index")

        pos = 0
        for read in SeqIO.parse("reads/" + assembly_name + ".fastq", "fastq"):
            alignments_list = list(assembly_aligner.map(str(read.seq)))
            for aln in alignments_list:
                if aln.mapq > 10 and (aln.r_en - aln.r_st) > 100:
                    ax_reads.add_patch(Rectangle(
                        (aln.r_st, pos),
                        aln.r_en - aln.r_st,
                        1,
                        facecolor=(0.9, 0.9, 0.9)))
                    tmp_bed.write("ctg\t" + str(aln.r_st) + "\t" + str(aln.r_en) + "\n")
                    pos += 1
                    ax_reads.set_ylim(0, pos)
        tmp_bed.close()

        # get coverage
        os.system("bedtools intersect -c -a tmp.perbase.bed -b tmp.bed > tmp.perbase.count.bed")
        #os.system("bedtools intersect -c -a tmp.perregion.bed -b tmp.bed > tmp.perregion.count.bed")

        cov_x_data = []
        cov_y_data = []
        for line in open("tmp.perbase.count.bed"):
            pos = int(line.split("\t")[1])
            cov = int(line.split("\t")[-1])
            cov_x_data.append(pos)
            cov_y_data.append(cov)

        ax_cov.plot([0, len(assembly_sequence)], [2, 2], c="red")
        ax_cov.plot(cov_x_data, cov_y_data)
        ax_cov.set_ylim(0, np.amax(cov_y_data) + 2)

        scaling = np.floor(np.amax(cov_y_data) / 2)
        for idx, sm in enumerate(sorted_mappings):
            y_pos = 0
            if idx % 2 == 0:
                y_pos = 1

            color = (0.9, 0.9, 0.9)
            if sm[6] == 1:
                color = (0.9, 0.6, 0.6)

            ax_cov.add_patch(Rectangle(
                (sm[0], y_pos * scaling),
                sm[1] - sm[0],
                scaling,
                facecolor=color))
    except FileNotFoundError:
        print("\t\t" + f"{bcolors.FAIL}" + "Assembly not found for QC Plot" + f"{bcolors.ENDC}")
    except IndexError:
        print("\t\t" + f"{bcolors.FAIL}" + "Cannot load assembly sequence" + f"{bcolors.ENDC}")

def assembly_to_vcf_record (sorted_mappings, id_name, expanded_fo, collapsed_fo):

    if len(sorted_mappings) == 1:

        has_long_insertion = False
        insertion_pos = 0

        start_pos = sorted_mappings[0][4]
        for cig in sorted_mappings[0][7]:
            if cig[1] in [0, 2]:
                start_pos += cig[0]
            if cig[1] == 1 and cig[0] > 500:
                insertion_pos = start_pos
                has_long_insertion = True


        if has_long_insertion:
            vcf_record = "\t".join(
                [
                    sorted_mappings[0][3],
                    str(insertion_pos),
                    id_name,
                    "N",
                    "<INS>",
                    ".",
                    "PASS",
                    "CHR2=" + sorted_mappings[0][3] + ";END=" + str(insertion_pos) + ";SVTYPE=INS",
                    "GT:DR:DV",
                    "0/1:2:1"
                ])
        else:
            vcf_record = "\t".join(
                [
                    sorted_mappings[0][3],
                    str(sorted_mappings[0][4]),
                    id_name,
                    "N",
                    "<UNKNOWN>",
                    ".",
                    "PASS",
                    "CHR2=" + sorted_mappings[0][3] + ";END=" + str(sorted_mappings[0][5]) + ";SVTYPE=UNKNOWN",
                    "GT:DR:DV",
                    "0/1:2:1"
                ])

        expanded_fo.write(vcf_record + "\n")
        collapsed_fo.write(vcf_record + "\n")

    # expanded form
    for i in range(0, len(sorted_mappings) - 1):
        expanded_fo.write(two_segments_to_vcf_record(sorted_mappings[i], sorted_mappings[i + 1], id_name + "_" + str(i)))
        expanded_fo.write("\n")

    # collapsed form
    vcf_record = two_segments_to_vcf_record(sorted_mappings[0], sorted_mappings[-1], id_name)
    # vcf_record # hier irgendwie noch die insertions in INFO schreiben
    #print(sorted_mappings)
    ins_segments = []
    for i in range(1, len(sorted_mappings) - 1):
        #print(sorted_mappings[i])
        ins_segments.append(sorted_mappings[i][3] + ":" + str(sorted_mappings[i][4]) + "-" + str(sorted_mappings[i][5]))
    segment_info = "SEGMENTS=" + ";".join(ins_segments)
    vcf_record_splits = vcf_record.rstrip().split("\t")
    vcf_record_splits[7] += ";" + segment_info
    collapsed_fo.write("\t".join(vcf_record_splits))
    collapsed_fo.write("\n")

def two_segments_to_vcf_record (segment_a, segment_b, name):

    chr1 = segment_a[3]
    pos1 = -1
    out_pos1 = -1
    chr2 = segment_b[3]
    pos2 = -1
    out_pos2 = -1

    if segment_a[2] == "+":
        pos1 = segment_a[5]
        out_pos1 = segment_a[4]
    else:
        pos1 = segment_a[4]
        out_pos1 = segment_a[5]

    if segment_b[2] == "+":
        pos2 = segment_b[4]
        out_pos2 = segment_b[5]
    else:
        pos2 = segment_b[5]
        out_pos2 = segment_b[4]

    if chr1 == chr2:
        # same orientation - either DUP or DEL
        if segment_a[2] == segment_b[2]:

            segment_a_in_sv = False
            if out_pos1 > np.amin([pos1, pos2]) and out_pos1 < np.amax([pos1, pos2]):
                segment_a_in_sv = True

            segment_b_in_sv = False
            if out_pos2 > np.amin([pos1, pos2]) and out_pos2 < np.amax([pos1, pos2]):
                segment_b_in_sv = True

            if segment_a_in_sv and segment_b_in_sv:
                svtype = "DUP"
            elif not segment_a_in_sv and not segment_b_in_sv:
                svtype = "DEL"
            else:
                svtype = "DELDUP"

        else:
            svtype = "INV"
    else:
        svtype = "BND"



    reference_allele = ""
    if chr1 == chr2:
        if svtype == "INV":
            reference_allele = "<INV>"
        if svtype == "DEL":
            reference_allele = "<DEL>"
        if svtype == "DUP":
            reference_allele = "<DUP>"
        if svtype == "DELDUP":
            reference_allele = "<DELDUP>"
        if pos1 > pos2:
            tmp = pos2
            pos2 = pos1
            pos1 = tmp
    else:
        reference_allele = "]" + chr2 + ":" + str(pos2) + "]N"

    record = "\t".join(
        [
            chr1,
            str(pos1),
            name,
            "N",
            reference_allele,
            ".",
            "PASS",
            "CHR2=" + chr2 + ";END=" + str(pos2) + ";SVTYPE=" + svtype,
            "GT:DR:DV",
            "0/1:2:1"
        ])

    return record

    """
    chr1	8595285	33212	N	]chr14:87478047]N	.	PASS	CHR2=chr14;END=87478047;RE=7;IMPRECISE;SVLEN=1;SVMETHOD=Snifflesv1.0.12;SVTYPE=BND;STRANDS=--;AF=0.112903;Kurtosis_quant_start=0.297521;Kurtosis_quant_stop=0.499986;REF_strand=27,35;RNAMES=0ffb80e6-1f38-40d6-bae2-59dfd240bb43,135bb061-88e3-4220-bd98-0efdeda8c917,2e4147ce-3207-4977-a741-bbb401c0a030,5a770d4a-dd37-4715-a6eb-15845c2f7bc8,5efdb766-63b0-4e06-99d7-711ae70c26fb_2,8a87336d-8c86-4d9f-b471-6fccbc6ab424,bfc1b46d-668f-4128-aa57-537d1693d921;STD_quant_start=2.17124;STD_quant_stop=371.493;STRANDS2=2,5,5,2;SUPTYPE=SR;Strandbias_pval=0.690481	GT:DR:DV	0/0:55:7
    chr1	12524813	326	N	<DEL>	.	PASS	CHR2=chr1;END=12554640;RE=9;PRECISE;SVLEN=-29827;SVMETHOD=Snifflesv1.0.12;SVTYPE=DEL;STRANDS=+-;AF=0.113924;Kurtosis_quant_start=-0.75;Kurtosis_quant_stop=1.1095;REF_strand=43,36;RNAMES=04a86973-cb0b-4535-95b7-9b90e900da6c,346fe074-56ce-427e-ac5c-75337d7d00c7,4cb9d22b-06f9-48b9-8146-a562a42f6df0,6acec029-5346-40cf-b179-f1c0d18060dc,79b23b34-51b4-48bc-a6cb-03f48e6ae938,8eb7dd1b-232a-4d0a-a698-23a70d1befcf,a2cb8928-6ae3-4c8e-9e47-327f6af3143a,e5f1ec9e-1a44-4b9d-b8c6-289778853a93,e872f4df-37ef-46ed-8c08-14ef95b56d47;STD_quant_start=0.666667;STD_quant_stop=2.21108;STRANDS2=5,4,5,4;SUPTYPE=SR;Strandbias_pval=1.0	GT:DR:DV	0/0:70:9
    """

def write_deletion_record (record, expanded_fo, collapsed_fo):
    record_string = "\t".join(
        [
            record.CHROM,
            str(record.POS),
            record.ID,
            "N",
            "<DEL>",
            ".",
            "PASS",
            "CHR2=" + record.CHROM + ";END=" + str(record.INFO["END"]) + ";SVTYPE=DEL",
            "GT:DR:DV",
            "0/1:2:1"
        ])

    expanded_fo.write(record_string + "\n")
    collapsed_fo.write(record_string + "\n")

def assembly_to_bed (sorted_mappings, id_name, bed_fo):
    for idx, sm in enumerate(sorted_mappings):
        bed_fo.write("\t".join([sm[3], str(sm[4]), str(sm[5]), id_name + "_" + str(idx)]) + "\n")
    bed_fo.flush()

def assemble_sv (chrom, start, end, id_name, num_seed, manual_blacklist, mappy_aligner, samfile, expanded_vcf_out, collapsed_vcf_out, bed_out, blacklist):

    segments = []
    segments.append([chrom, start, end])
    segments_added = True

    j = 0

    sorted_mappings = []

    print("\tRunning assemble_sv for SV " + id_name + " seed " + str(num_seed + 1))

    breakpoint_distance_tolerance = 10

    while segments_added:

        # maximum of 100 cycles
        j += 1
        if j == 100:
            break

        """
        for sm in segments:
            if sm[0] == "chr8":
                breakpoint_distance_tolerance = 1
        """

        segments_added = False
        name = id_name + ".assembly." + str(j)

        print("\t\tRound " + str(j))
        # get reads from segments
        sv_reads = get_sv_reads_from_segment(segments, samfile, blacklist, breakpoint_distance_tolerance)
        print("\t\t\tGathered " + str(len(sv_reads)) + " reads from " + str(len(segments)) + " segments")
        for sm in segments:
            print("\t\t\t\t" + str(sm))

        if len(sv_reads) < 2 or len(sv_reads) > 250:
            print("\t\t" + f"{bcolors.FAIL}" + "Aborting assemble_sv for " + name + " due to insufficient read count" + f"{bcolors.ENDC}")
            break

        # write reads to fastq file
        tmp_fastq_fo = open("reads/" + name + ".fastq", "w")
        for read in sv_reads:
            if not (read.query_name) in manual_blacklist:
                read_methods.write_fastq_original(read, tmp_fastq_fo)
        tmp_fastq_fo.close()

        # assemble reads
        read_methods.assemble_reads_lamassemble_s("reads/" + name + ".fastq", "lamassemble",
                                                  "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/train/promethion.mat",
                                                  name, "assemblies/" + name + ".fasta", "2")


        assembly_name = list(SeqIO.parse("assemblies/" + name + ".fasta", "fasta"))[0].id
        assembly_sequence = list(SeqIO.parse("assemblies/" + name + ".fasta", "fasta"))[0].seq

        print("\t\t\tProduced assembly " + assembly_name + " with " + str(len(assembly_sequence)) + " bases")

        # prune assembly, reload assembly if necessary
        pruning = prune_final_assembly(name)
        if pruning:
            try:
                assembly_name = list(SeqIO.parse("assemblies/" + name + ".fasta", "fasta"))[0].id
                assembly_sequence = list(SeqIO.parse("assemblies/" + name + ".fasta", "fasta"))[0].seq

                print("\t\t\tProduced pruned assembly " + assembly_name + " with " + str(len(assembly_sequence)) + " bases")
            except IndexError:
                return

        # map assembly to reference
        alignments_list = list(mappy_aligner.map(str(assembly_sequence)))
        num_alns_with_mq_larger_1 = 0
        for aln in alignments_list:
            if aln.mapq > 1:
                num_alns_with_mq_larger_1 += 1
        print("\t\t\tFound " + str(len(alignments_list)) + " alignments for the assembly, " + str(num_alns_with_mq_larger_1) + " alignments have MAPQ > 1")

        # get assembly mappings
        sorted_mappings = get_assembly_mappings_mappy(alignments_list, len(assembly_sequence), 100, 1)

        if len(alignments_list) == 1:
            print("\t\t" + f"{bcolors.WARNING}" + "Assembly maps in one part -> aborting" + f"{bcolors.ENDC}")
            break

        if len(alignments_list) == 0:
            print("\t\t" + f"{bcolors.FAIL}" + "Aborting assemble_sv" + f"{bcolors.ENDC}")
            break

        # plot
        fig, ax = plt.subplots()
        plot_query_mappings(sorted_mappings, ax, len(assembly_sequence), segments)
        plt.savefig("plots/" + name + ".png")
        plt.close()

        # print identified segments
        for sm in sorted_mappings:
            in_segments = False
            half_in_segments = False
            segment_id = -1

            for idx, seg in enumerate(segments):
                if seg[0] == sm[3] and np.abs(seg[1] - sm[4]) < 10 and np.abs(seg[2] - sm[5]) < 10:
                    in_segments = True
                    segment_id = idx
                if seg[0] == sm[3]:
                    if np.abs(seg[1] - sm[4]) < 10:
                        half_in_segments = True
                        segment_id = idx
                    if np.abs(seg[2] - sm[5]) < 10:
                        half_in_segments = True
                        segment_id = idx

            if in_segments:
                print("\t\t\t\t" + f"{bcolors.OKGREEN}" + str(sm[:-1]) + f"{bcolors.ENDC}")
            elif half_in_segments:
                print("\t\t\t\t" + f"{bcolors.WARNING}" + str(sm[:-1]))  # direkt hier ersetzen?
                print("\t\t\t\t\tSegment fits only on one side to the original segment")
                print("\t\t\t\t\tChanging segment " + segments[segment_id][0] + ":" + str(segments[segment_id][1]) + "-" + str(segments[segment_id][2]) + " to " + sm[3] + ":" + str(sm[4]) + "-" + str(sm[5]) + f"{bcolors.ENDC}")
                segments[segment_id] = [sm[3], sm[4], sm[5]]
            else:
                print("\t\t\t\t" + str(sm[:-1]))


        # add adjacent segments to segments

        print("\t\t\tAdding adjacent segments")

        new_segments = []
        for seg in segments:

            new_segments.append(seg)

            # get adjacent mappings to seg and add them to "new_segments"
            adj_mappings = get_adjacent_mappings(sorted_mappings, seg[0], seg[1], seg[2])

            num_added = 0
            if adj_mappings == None:
                print("\t\t\t\tFound no adjacent mappings for mapping " + seg[0] + ":" + str(seg[1]) + "-" + str(seg[2]))
            else:
                for x in adj_mappings:
                    mapping = sorted_mappings[x]

                    # add mapping only if it is not yet present in "segments"
                    add_mapping = True
                    for all_old_seg in segments:
                        if all_old_seg[0] == mapping[3] and np.abs(all_old_seg[1] - mapping[4]) < 10 and np.abs(
                                all_old_seg[2] - mapping[5]) < 10:
                            add_mapping = False

                    if add_mapping:
                        new_segments.append([mapping[3], mapping[4], mapping[5]])
                        # print("\t\t\tAdded segment " + mapping[3] + ":" + str(mapping[4]) + "-" + str(mapping[5]))
                        num_added += 1
                        segments_added = True

                print("\t\t\t\tFound " + str(len(adj_mappings)) + " adjacent mappings for mapping " + seg[0] + ":" + str(
                    seg[1]) + "-" + str(seg[2]) + ", added " + str(num_added) + " segments")

        print("\t\t\tSV assembly spans " + str(len(new_segments)) + " segments")

        segments = new_segments

    # manage output if assembly maps
    if len(sorted_mappings) > 0:

        # if sum of mapping lengths is much smaller than total assembly length
        distance_start = sorted_mappings[0][0]
        distance_end = len(assembly_sequence) - sorted_mappings[-1][1]
        long_unmapped_assembly_strech = False
        if np.any([distance_start > 5000, distance_end > 5000, distance_start / len(assembly_sequence) > 0.4,
                   distance_end / len(assembly_sequence) > 0.4]):
            long_unmapped_assembly_strech = True
        if long_unmapped_assembly_strech:
            sorted_mappings = get_assembly_mappings_mappy(alignments_list, len(assembly_sequence), 100, 0)

        # QC plot
        ax_qc_1 = plt.subplot(2, 1, 1)
        ax_qc_2 = plt.subplot(2, 1, 2)
        plot_qc_plot(id_name + ".assembly." + str(j), sorted_mappings, ax_qc_1, ax_qc_2)
        plt.savefig("final/qcplots/" + id_name + ".png")
        plt.close()

        # Copy results to new folder
        os.system("cp plots/" + id_name + ".assembly." + str(j) + ".png final/plots/" + id_name + ".assembly." + str(j) + ".png")
        os.system("cp assemblies/" + id_name + ".assembly." + str(j) + ".fasta final/assemblies/" + id_name + ".assembly." + str(j) + ".fasta")
        os.system("cp reads/" + id_name + ".assembly." + str(j) + ".fastq final/reads/" + id_name + ".assembly." + str(j) + ".fastq")

        # Assembly to VCF record
        assembly_to_vcf_record(sorted_mappings, id_name, expanded_vcf_out, collapsed_vcf_out)
        expanded_vcf_out.flush()
        collapsed_vcf_out.flush()

        # Assembly to bed record
        assembly_to_bed(sorted_mappings, id_name, bed_out)
        bed_out.flush()

        #os.system("minimap2 -ax map-ont -t 26 /home/ubuntu/seq/2022_02_WGS_LSK/reference/chm13v2.0.mmi final/reads/" + id_name + "*.fastq -Y | samtools view -b - | samtools sort > final/" + id_name + ".bam")
        #os.system("samtools index final/" + id_name + ".bam")
        #os.system("rm tmp.fasta")

    else:
        print(f"{bcolors.FAIL}" + "ERROR for record " + str(id_name) + f"{bcolors.ENDC}")

def map_assemblies (reference):
    os.system("cat final/assemblies/*.fasta > tmp.fasta")
    os.system("minimap2 -ax map-ont -t 26 " + reference + " tmp.fasta -Y | samtools view -b - | samtools sort > final/assemblies_to_reference.bam")
    os.system("samtools index final/assemblies_to_reference.bam")
    os.system("rm tmp.fasta")

def assemble_record (vcf_record, vcf_sv_assemblies_folder, mappy_aligner, samfile, expanded_vcf_out, collapsed_vcf_out, bed_out, blacklist):

    # seeding

    print("Running assembly method for record " + vcf_record.ID)



    try:
        sv_assembly_seq = str(list(SeqIO.parse(vcf_sv_assemblies_folder + vcf_record.ID + ".sv.fasta", "fasta"))[0].seq)
        aln = list(mappy_aligner.map(sv_assembly_seq))
    except FileNotFoundError:
        print("\tSV assembly not found")
        aln = []




    seed_mappings = []
    for a in aln:
        if a.mapq > 0:

            long_deletion = False

            if not (vcf_record.INFO["SVTYPE"] == "BND"):

                # check for long insertion

                pos = a.r_st

                for cig_op in a.cigar:
                    if cig_op[1] == 0:
                        pos += cig_op[0]
                    if cig_op[1] == 2:
                        if cig_op[0] > 5000:
                            if np.abs(cig_op[0] - np.abs(record.INFO["SVLEN"])) < 100:
                                long_deletion = True
                        pos += cig_op[0]


            seed_mappings.append(
                [
                    a.mapq,
                    a.ctg,
                    a.r_st,
                    a.r_en,
                    long_deletion
                ])

    print("\tFound " + str(len(seed_mappings)) + " seeds")
    for i in range(0, np.amin([3, len(seed_mappings)])):
        print("\t\t" + seed_mappings[i][1] + ":" + str(seed_mappings[i][2]) + "-" + str(seed_mappings[i][3]) + ", MQ: " + str(seed_mappings[i][0]) + ", LENGTH: " + str(seed_mappings[i][3] - seed_mappings[i][2]) + ", DEL: " + str(seed_mappings[i][-1]))

    # only one mapping and that mapping carries a long deletion
    
    has_long_insertion = False
    if len(seed_mappings) <= 2:
        for seed in seed_mappings:
            if seed[-1]:
                has_long_insertion = True
    

    if len(seed_mappings) < 2:# and has_long_insertion:
        print("\tDetected assembly with only one mapping --> writing original record")
        write_deletion_record(record, expanded_vcf_out, collapsed_vcf_out)

    else:

        for num, seed in enumerate(seed_mappings):

            assemble_sv(seed[1], seed[2], seed[3], record.ID + "_" + str(num), num, [], mappy_aligner, samfile, expanded_vcf_out, collapsed_vcf_out, bed_out, blacklist)

            # only the 2 best seeds
            if num == 1:
                break
