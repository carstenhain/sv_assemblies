import vcf
import mappy as mp
import pysam
import read_methods
from Bio import SeqIO
import os
import numpy as np
import assemble_structural_variants


# removes terminal streched with coverage lower than 2 from the assembly
# uses bedtools, samtools and mappy
# return True if terminal bases are pruned
def prune_assembly (asm_file_path, fastq_file_path, prune_cutoff):

    asm_name = list(SeqIO.parse(asm_file_path, "fasta"))[0].id
    asm_sequence = list(SeqIO.parse(asm_file_path, "fasta"))[0].seq
    assembly_aligner = mp.Aligner(asm_file_path, preset="map-ont")
    if not assembly_aligner: raise Exception("ERROR: failed to load/build index")

    # map reads to assembly and write reads to bed
    tmp_reads_bed = open("tmp.prune.reads.bed", "w")
    for read in SeqIO.parse(fastq_file_path, "fastq"):
        for aln in list(assembly_aligner.map(str(read.seq))):
            if aln.mapq > 10 and (aln.r_en - aln.r_st) > 100:
                tmp_reads_bed.write("ctg\t" + str(aln.r_st) + "\t" + str(aln.r_en) + "\n")
    tmp_reads_bed.close()

    # build single base bed file
    tmp_per_base_bed = open("tmp.prune.perbase.bed", "w")
    for i in range(0, len(asm_sequence)):
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
        if coverage_data[i][1] < prune_cutoff:
            prune_start_position += 1
        else:
            break

    # prune end
    prune_end_position = len(asm_sequence)
    for i in range(0, len(coverage_data)):
        if coverage_data[len(coverage_data) - i - 1][1] < prune_cutoff:
            prune_end_position = len(coverage_data) - i
        else:
            break

    if np.any([prune_start_position != 0, prune_end_position != len(asm_sequence)]):

        # extract and save unpruned sequence using bedtools getfasta

        tmp_pruning_bed = open("pruning.bed", "w")
        tmp_pruning_bed.write("\t".join(
            [
                asm_name,
                str(prune_start_position),
                str(prune_end_position)
            ]) + "\n")
        tmp_pruning_bed.close()

        os.system("samtools faidx " + asm_file_path)
        os.system("bedtools getfasta -fi " + asm_file_path + " -bed pruning.bed > " + asm_file_path.replace(".fasta", ".pruned.fasta"))

        return True

    return False

def assemble_record (record, mappy_aligner, samfile, workdir, expanded_vcf_fo, collapsed_vcf_fo, bed_fo, read_blacklist):

    # gather SV supporting reads

    reads_at_start = {}
    for read in samfile.fetch(record.CHROM, record.POS - 150, record.POS + 150):
        if read.has_tag("SA"):
            if not (read.query_name in reads_at_start):
                reads_at_start[read.query_name] = read

    sv_reads = []
    for read in samfile.fetch(record.INFO["CHR2"], record.INFO["END"] - 150, record.INFO["END"] + 150):
        if read.query_name in reads_at_start:
            sv_reads.append(read)

    fastq_fo = open(workdir + record.ID + "_0.fastq", "w")
    for read in sv_reads:
        read_methods.write_fastq_original(read, fastq_fo)
    fastq_fo.close()

    read_methods.assemble_reads_lamassemble_s(workdir + record.ID + "_0.fastq", "lamassemble", "/home/ubuntu/seq/Miltenyi_ONT/lamassemble/lamassemble/train/promethion.mat", record.ID + "_0", workdir + record.ID + "_0.asm.fasta", "3")

    try:
        asm_name = list(SeqIO.parse(workdir + record.ID + "_0.asm.fasta", "fasta"))[0].id
    except IndexError:
        return

    prune_assembly(workdir + record.ID + "_0.asm.fasta", workdir + record.ID + "_0.fastq", 3)

    pruned_assembly_seq = list(SeqIO.parse(workdir + record.ID + "_0.asm.fasta", "fasta"))[0].seq

    aln = list(mappy_aligner.map(str(pruned_assembly_seq)))

    seed_mappings = []
    for a in aln:
        if a.mapq > 20 and (a.r_en - a.r_st) > 1000:
            seed_mappings.append(
                [
                    a.mapq,
                    a.ctg,
                    a.r_st,
                    a.r_en,
                    False
                ])

    expanded_vcf_out = open(workdir + "exp.vcf", "w")
    expanded_vcf_out = open(workdir + "exp.vcf", "w")
    if len(seed_mappings) > 0:
        assemble_structural_variants.assemble_sv\
            (
                seed_mappings[0][1],
                seed_mappings[0][2],
                seed_mappings[0][3],
                record.ID, 1, [],
                mappy_aligner,
                samfile,
                expanded_vcf_fo,
                collapsed_vcf_fo,
                bed_fo,
                read_blacklist
            )


if __name__ == '__main__':

    os.system("rm -r final")
    os.system("rm -r reads")
    os.system("rm -r assemblies")
    os.system("rm -r plots")

    os.system("mkdir final")
    os.system("mkdir reads")
    os.system("mkdir assemblies")
    os.system("mkdir plots")

    os.system("mkdir final/assemblies")
    os.system("mkdir final/reads")
    os.system("mkdir final/plots")
    os.system("mkdir final/qcplots")

    ids_to_assemble = []
    read_blacklist = []

    sample_id = ""
    vcf_file_path = ".vcf"
    mapping_file_path = ".bam"
    reference_file_path = "reference/chm13v2.0.mmi"

    mappy_aligner = mp.Aligner(reference_file_path, preset="map-ont")
    if not mappy_aligner: raise Exception("ERROR: failed to load/build index")

    samfile = pysam.AlignmentFile(mapping_file_path, "rb")
    vcf_reader = vcf.Reader(open(vcf_file_path, "r"))

    workdir = ""

    expanded_vcf_out = open(workdir + "expanded.vcf", "w")
    collapsed_vcf_out = open(workdir + "collapsed.vcf", "w")
    bed_out = open(workdir + "segments.bed", "w")

    assemble_structural_variants.write_generic_vcf_header(expanded_vcf_out, "/home/ubuntu/seq/2022_02_WGS_LSK/Python/generic.vcf")
    assemble_structural_variants. write_generic_vcf_header(collapsed_vcf_out, "/home/ubuntu/seq/2022_02_WGS_LSK/Python/generic.vcf")

    for record in vcf_reader:

        if record.ID in ids_to_assemble:#record.INFO["SVTYPE"] != "INS" and record.FILTER == []:#
            assemble_record(record, mappy_aligner, samfile, workdir, expanded_vcf_out, collapsed_vcf_out, bed_out, read_blacklist)

    expanded_vcf_out.close()
    collapsed_vcf_out.close()
    bed_out.close()

    # mapping der assemblies auf das Refenzgenome
    assemble_structural_variants.map_assemblies("/home/ubuntu/seq/2022_02_WGS_LSK/reference/chm13v2.0.mmi")
