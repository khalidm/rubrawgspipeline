#!/bin/env python

"""
GATK-based variant-calling pipeline, exome version.

Authors: Bernie Pope, Clare Sloggett, Gayle Philip.
Thanks to Dmitri Mouradov and Maria Doyle for input on the initial
analysis design.
Thanks to Matt Wakefield for contributions to Rubra
(https://github.com/bjpop/rubra) during pipeline development.

Description:

This program implements a workflow pipeline for next generation
sequencing variant detection using the Broad Institute's GATK for
variant calling and using ENSEMBL for basic annotation.

It uses Rubra (https://github.com/bjpop/rubra) based on the
Ruffus library.

It supports parallel evaluation of independent pipeline stages,
and can run stages on a cluster environment.

The pipeline is configured by an options file in a python file,
including the actual commands which are run at each stage.
"""


import sys
import re
import os.path
import os
import string
from collections import defaultdict
from glob import *
# from __future__ import print_function
import shutil
from ruffus import *
from rubra.utils import pipeline_options
from rubra.utils import (runStageCheck, mkLogFile, mkDir, mkForceLink, mkLink)
from input_fastq import parse_and_link

def make_metadata_string(metadata):
    return r'-r"@RG\tID:%s\tPU:%s\tSM:%s\tPL:%s"' % (metadata['ID'], metadata['PU'], metadata['SM'], metadata['PL'])

# Shorthand access to options
ref_files = pipeline_options.ref_files
working_files = pipeline_options.working_files
logDir = pipeline_options.pipeline['logDir']

# Data setup process and input organisation and metadata functions

#Metadata holding structures
fastq_metadata = defaultdict(dict)
#sample_metadata = defaultdict(dict)

original_fastq_files = []
for fastq_dir in working_files['fastq_dirs']:
    #original_fastq_files += glob(os.path.join(fastq_dir, '*.fastq.gz'))
    #original_fastq_files += glob(os.path.join(fastq_dir, '*_sequence.txt.gz'))
    original_fastq_files += glob(os.path.join(fastq_dir, '*.fastq.gz'))

if len(original_fastq_files)==0:
    print "No input files found. Do the filenames follow the naming convention?"
    print "Directories searched:"
    print "\n".join(working_files['fastq_dirs'])
    sys.exit(1)

# Parse metadata out of input file names and construct symlinks
# Metadata is put into a dict (for the rest of ruffus) and some of it also into symlinks (for filename uniqueness)
# currently parsing by assuming AGRF naming structure and paired-end reads
mkDir(working_files['fastq_symlink_dir'])
all_fastq_files = []
for file in original_fastq_files:
    name = os.path.basename(file)
    #print name
    #pu,remaining = name.split("_")
    #barcode = pu.split(".")[0]
    #lane = pu.split(".")[1]
    #id = pu[:5] + "." + lane
    #sm = os.path.basename(os.path.dirname(file))
    #sm = string.replace(sm,"_","-")

    #print "sm = " + sm + "\npu = " + pu + "\nlane = " + lane + "\nid = " + id + "\nbarcode = " + barcode

    ##file = sm + "_" + name
    #file2 = sm + "_" + barcode + ".L" + lane + "_" + remaining
    #print file + " -> " + file2
    symlink = parse_and_link(file, working_files['fastq_symlink_dir'], fastq_metadata)
    all_fastq_files.append(symlink)

# Make a list of files we will actually use
if pipeline_options.pipeline['restrict_samples']:
    allowed_samples = set(pipeline_options.pipeline['allowed_samples'])
    fastq_files = [file for file in sorted(all_fastq_files)
            if (fastq_metadata[os.path.basename(file)]['sample'] in allowed_samples)]
else:
    fastq_files = sorted(all_fastq_files)

# print "Symlinked files that will be used:"
# for file in fastq_files:
#     print file
# print
# print "Output dir is %s" % working_files['output_dir']
# print "Log dir is %s" % logDir
# print

# Create output subdirectories

output_dir = working_files['output_dir']
mkDir(output_dir)

# directory for final summary tables
ref_dir = os.path.join(output_dir, "reference")
mkDir(ref_dir)

fastqc_dir = os.path.join(output_dir, "FastQC")
mkDir(fastqc_dir)

sambam_dir = os.path.join(output_dir, "alignments")
mkDir(sambam_dir)

variant_dir = os.path.join(output_dir, "variant_calls")
mkDir(variant_dir)

varscan_temp_dir = os.path.join(variant_dir, "varscan_temp")
mkDir(varscan_temp_dir)
# varscan_temp_file = os.path.join(varscan_temp_dir, "ALL_SAMPLENAME")
# mkLogFile(varscan_temp_file)

coverage_dir = os.path.join(output_dir, "coverage")
mkDir(coverage_dir)

ensembl_dir = os.path.join(output_dir, "ensembl")
mkDir(ensembl_dir)

snpeff_dir = os.path.join(output_dir, "snpeff")
mkDir(snpeff_dir)

# directory for final summary tables
results_dir = os.path.join(output_dir, "results")
mkDir(results_dir)

# Pipeline declarations

# Making references

#Reference file setup
RefName=ref_files['fasta_reference'].split("/")[-1]
print RefName
fasta_reference=os.path.join(ref_dir, RefName)
mkLink(ref_files['fasta_reference'], fasta_reference)

RefDict=ref_files['fasta_dict'].split("/")[-1]
print RefDict
# print ref_files['fasta_dict']
fasta_dict=os.path.join(ref_dir, RefDict)
# print fasta_dict
mkLink(ref_files['fasta_dict'], fasta_dict)

# sys.exit(0)

# RefDict_temp=ref_files['fasta_reference'].split("/")[-1]
# RefDict=RefDict_temp.split(".")[0]+".dict"
# ref_dir=ref_dir+"/"
# print RefDict
# print ref_dir
# fasta_dict=os.path.join(ref_dir, RefDict)
# mkLink(ref_files['fasta_reference'], fasta_dict)
# sys.exit(0)

@files(fasta_reference, fasta_reference+".indexReferenceBWA.Success")
def indexRefbwa(inputs, outputs):
    """
    Index reference file for use with BWA.
    """
    ref = inputs
    flagFile = outputs
    runStageCheck('indexReferenceBWA', flagFile, ref)

@files(fasta_reference, [fasta_reference+".fai", fasta_reference+".indexReferenceSAM.Success"])
def indexRefSamtools(inputs, outputs):
    """
    Index reference file for use with samtools.
    """
    ref = inputs
    out, flagFile = outputs
    runStageCheck('indexReferenceSAM', flagFile, ref)

# Gzip fastq files
# @files(fastq_files, regex('(.+\/)?(.+?)\.fastq\.gz')
# @files(fastq_files, regex('(.+\/)?(.+?)\.fastq\.gz'))
#def gzip(inputs, outputs):
#    """
#    Gzip all the split fastq files.
#    """
#    print "INPUT =" + inputs
#    _GZIPFastq_success, file = inputs
#    flagFile, zipped_file = outputs
#    runStageCheck('gzip', flagFile, file)

# Alignment and correction steps

@transform(fastq_files, regex('(.+\/)?(.+?)\.fastq\.gz'),
        [r'%s/\2_fastqc' % fastqc_dir, r'%s/\2.fastqc.Success' % fastqc_dir])
def fastqc(inputs, outputs):
    """
    Run FastQC on each fastq file.
    """
    # print "HERE IN FastQC"
    sequence = inputs
    # print "INPUTS = " + inputs
    fastqc_dest, flagFile = outputs
    # print "Success = " + flagFile
    runStageCheck('fastqc', flagFile, fastqc_dir, sequence)

@follows(indexRefbwa, indexRefSamtools)
@transform(fastq_files, regex(r".*?(([^/]+)(_1|_2))\.fastq.gz"),
        [r"%s/\1.sai" % sambam_dir, r"%s/\1.alignBwa.Success" % sambam_dir])
def alignBWA(inputs, outputs):
    """
    Align sequence reads to the reference genome. This is bwa's first stage, bwa aln.
    Use -I for _sequence.txt files.
    """
    seq = inputs
    output, flag_file = outputs
    encodingflag = ''
    if fastq_metadata[os.path.basename(seq)]['encoding'] == 'I':
        encodingflag = '-I'
    print "bwa aln on %s" % os.path.basename(seq)
    runStageCheck('alignBWA', flag_file, encodingflag, fasta_reference, seq, output)

# Convert alignments to SAM format.
# This assumes paired-end; if we have single end we should wrap in a conditional and in the other case
#   define with @transform not @collate, and call SamSE not SamPE
# @collate(alignBWA, regex(r"(.*?)([^/]+)(_1|_2)\.sai"),
@collate(alignBWA, regex(r"(.*?)([^/]+)(_1|_2)\.sai"),
# @collate(alignBWA, regex(r"(.*?)([/]+)(.*)(_.*)(_1|_2)\.sai"),
        add_inputs(r"%s/\2\3.fastq.gz" % working_files['fastq_symlink_dir']),
        # add_inputs(r"%s/\2\3\4\5.fastq.gz" % working_files['fastq_symlink_dir']),
        [r"\1\2.sam", r"\1\2.alignToSam.Success"])
def alignToSam(inputs, outputs):
    """
    Turn two paired-end bwa "sai" alignments into a sam file.
    """
    output,flag_file = outputs
    [sai1, seq1], [sai2, seq2] = [[sai, seq] for [[sai, _flag_file], seq] in inputs]
    fastq_name = os.path.splitext(os.path.basename(sai1))[0] + ".fastq.gz"
    sample = fastq_metadata[fastq_name]['sample']
    puID = fastq_metadata[fastq_name]['run_id']
    runID = fastq_metadata[fastq_name]['run_id'][:5]
    # runID = fastq_metadata[fastq_name]['run_id']
    lane = fastq_metadata[fastq_name]['run_id'][-1:]
    readgroup_metadata = { 'PL': 'ILLUMINA',
                           'SM': sample,
                           'ID': "%s.%s" % (runID, lane),
                           'PU': puID }
                           #'PU': "%s_%s_Lane%d" % (sample, runID, lane) }
    metadata_str = make_metadata_string(readgroup_metadata)
    # print "bwa sampe on %s,%s" % (os.path.basename(sai1), os.path.basename(sai2))
    runStageCheck('alignToSamPE', flag_file, fasta_reference, metadata_str, sai1, sai2, seq1, seq2, output)

@transform(alignToSam, suffix(".sam"),
            [".bam", ".samToBam.Success"])
def samToBam(inputs, outputs):
    """
    Convert sam to bam and sort, using Picard.
    """
    output, flag_file = outputs
    sam, _success = inputs
    print "converting to sorted bam: %s" % os.path.basename(sam)
    runStageCheck('samToSortedBam', flag_file, sam, output)

# r"(.*?)([^/]+)(_1|_2)\.sai"
@collate(samToBam, regex(r"(.*?)([^/]+)_([^/]+)\.bam"),
            [r"\1\3/\3.merged.bam", r'\1\3/\3.mergeBams.Success'])
            #[r"\1\3/\3.merged.bam", r'\1\3/\3.mergeBams.Success'])
def mergeBams(inputs, outputs):
    """
    Merge the sorted bams together for each sample.
    Picard should cope correctly if there is only one input.
    """
    bams = [bam for [bam, _success] in inputs]
    output_bam, flag_file = outputs
    baminputs = ' '.join(["INPUT=%s" % bam for bam in bams])
    print "merging %s into %s" % (",".join([os.path.basename(bam) for bam in bams]), os.path.basename(output_bam))
    runStageCheck('mergeBams', flag_file, baminputs, output_bam)

@follows('indexMergedBams')
@transform(mergeBams, suffix('.bam'),
            ['.dedup.bam', '.bam.dedup.Success'])
def dedup(inputs, outputs):
    """
    Remove apparent duplicates from merged bams using Picard MarkDuplicates.
    """
    input_bam, _success = inputs
    output_bam, flag_file = outputs
    logFile = mkLogFile(logDir, input_bam, '.dedup.log')
    print "de-duping %s" % os.path.basename(input_bam)
    runStageCheck('dedup', flag_file, input_bam, logFile, output_bam)

@follows('indexDedupedBams')
@transform(dedup, suffix('.bam'),
            ['.realigner.intervals', '.bam.realignIntervals.Success'])
def realignIntervals(inputs, outputs):
    """
    Run GATK RealignTargetCreator to find suspect intervals for realignment.
    """
    bam, _success = inputs
    output_intervals, flag_file = outputs
    logFile = mkLogFile(logDir, bam, '.realignIntervals.log')
    print "calculating realignment intervals for %s" % os.path.basename(bam)
    runStageCheck('realignIntervals', flag_file, fasta_reference, bam, ref_files['indels_realign_goldstandard'], ref_files['indels_realign_1000G'], ref_files['exon_bed'], logFile, output_intervals)
    # runStageCheck('realignIntervals', flag_file, fasta_reference, bam, ref_files['indels_realign_goldstandard'], ref_files['indels_realign_1000G'], ref_files['exon_bed'], logFile, output_intervals)

def remove_GATK_bai(bamfile):
    """
    A bug in some versions of GATK cause it to create an x.bai file, and this gets in the way of using the properly named x.bam.bai file. If the given file exists, delete it.
    """
    bad_bai = os.path.splitext(bamfile)[0] + ".bai"
    try:
        os.remove(bad_bai)
    except OSError, e:
        # Ignore error only if it is OSError #2, ie File Not Found
        if e.errno != 2:
            raise e

@transform(realignIntervals, regex(r"(.*?)([^/]+)\.realigner\.intervals"),
            add_inputs([r'\1\2.bam']),
            [r'\1\2.realigned.bam', r'\1\2.bam.realign.Success'])
def realign(inputs, outputs):
    """
    Run GATK IndelRealigner for local realignment, using intervals found by realignIntervals.
    """
    [intervals, _success], [input_bam] = inputs
    output_bam, flag_file = outputs
    logFile = mkLogFile(logDir, input_bam, '.realign.log')
    print "realigning %s" % os.path.basename(input_bam)
    runStageCheck('realign', flag_file, fasta_reference, input_bam, ref_files['indels_realign_goldstandard'], ref_files['indels_realign_1000G'], intervals, logFile, output_bam)
    # runStageCheck('realign', flag_file, fasta_reference, input_bam, intervals, logFile, output_bam)
    remove_GATK_bai(output_bam)

@transform(realign, suffix('.bam'), ['.bam.bai', '.indexRealignedBams.Success'])
def indexRealignedBam(inputs, outputs):
    """
    Index the locally realigned bams using samtools.
    """
    bam, _realign_success = inputs
    output, flagFile = outputs
    print "samtools index on %s" % bam
    runStageCheck('indexBam', flagFile, bam)

@transform([realign, indexRealignedBam], suffix('.bam'), ['.leftAligned.bam', '.leftAlignIndels.Success'])
def leftAlign(inputs,outputs):
    """
    GATK LeftAlignIndels is a tool that takes a bam file and left-aligns any indels inside it
    'command': "java -Xmx22g -jar " + GATK_HOME + "GenomeAnalysisTK.jar -allowPotentiallyMisencodedQuals -T LeftAlignIndels -I %input -R %ref -o %output"
    """
    bam, _realign_success = inputs
    output_bam, flagFile = outputs
    runStageCheck('leftalignindels', flagFile, bam, fasta_reference, output_bam)
    remove_GATK_bai(output_bam)

@transform(leftAlign, suffix('.bam'), ['.bam.bai', '.indexleftAlignBam.Success'])
def indexleftAlignBam(inputs, outputs):
    """
    Index the locally realigned bams using samtools.
    """
    bam, _leftAlign_success = inputs
    output, flagFile = outputs
    print "samtools index on %s" % bam
    runStageCheck('indexBam', flagFile, bam)


@transform([leftAlign,indexleftAlignBam], suffix('.bam'), ['.recal_data.grp', '.baseQualRecal.Success'])
def baseQualRecal(inputs, outputs):
    """
    GATK BaseRecalibrator, first step of base quality score recalibration.
    'command': "java -Xmx22g -jar /vlsci/VR0245/shared/charlotte-working/programs/GenomeAnalysisTKLite-2.3-9-gdcdccbb/GenomeAnalysisTKLite.jar -T BaseRecalibrator -I %bam -R %ref --knownSites %dbsnp -nt 8 -log %log -o %out"
    """
    bam, _leftAlign_success = inputs
    output_grp, flagFile = outputs
    logFile = mkLogFile(pipeline_options.pipeline['logDir'], bam, '.baseQualRecal.log')
    print "Base Quality recal using GATK for: %s" % bam
    runStageCheck('baseQualRecal', flagFile, bam, fasta_reference, ref_files['dbsnp'], ref_files['indels_realign_goldstandard'], logFile, output_grp)

@transform([baseQualRecal, realign], regex(r'(.*?)([^/]+)\.recal_data\.grp'),
            add_inputs([r'\1\2.bam']),
            [r'\1\2.recal.bam', r'\1\2.baseQualRecalPrint.Success'])
def baseQualRecalPrint(inputs, outputs):
    """
    GATK TableRecalibration: write reads after base quality scores using the output of baseQualRecal.
    'command': "java -Xmx7g -jar /vlsci/VR0245/shared/charlotte-working/programs/GenomeAnalysisTKLite-2.3-9-gdcdccbb/GenomeAnalysisTKLite.jar -T PrintReads -I %bam -R %ref -BQSR %csvfile -log %log -o %out"
    """
    [[input_grp, _baseQualRecal_success], [input_bam]] = inputs
    output_bam, flagFile = outputs
    logFile = mkLogFile(pipeline_options.pipeline['logDir'], input_bam, '.baseQualRecalTabulate.log')
    print "recalibrate base quality scores using GATK on %s" % input_bam
    runStageCheck('baseQualRecalPrintReads', flagFile, input_bam, fasta_reference, input_grp, logFile, output_bam)
    remove_GATK_bai(output_bam)

@transform(baseQualRecalPrint, suffix('.bam'), ['.bam.bai', '.indexRealignedBams.Success'])
def indexBaseQualRecalBam(inputs, outputs):
    """
    Index the locally realigned bams using samtools.
    """
    bam, _baseRecalBam_success = inputs
    output, flagFile = outputs
    print "samtools index on %s" % bam
    runStageCheck('indexBam', flagFile, bam)

# Variant calling steps

#@follows(indexBaseQualRecalBam)
#@transform(baseQualRecalPrint,
#            regex(r'(.*?)([^/]+)\.recal\.bam'),
#            [r'%s/\2.recal.UGT.SNP.vcf' % variant_dir,
#             r'%s/\2.recal.UGT.SNP.vcf.idx' % variant_dir,
#             r'%s/\2.recal.UGT.callSNPs.Success' % variant_dir])
#def callSNPs(inputs, outputs):
#    """
#    Use GATK UnifiedGenotyper to call SNPs from recalibrated bams.
#    """
#    bam, _success = inputs
#    output_vcf, _idx, flag_file = outputs
#    logFile = mkLogFile(logDir, bam, '.callSNPs.log')
#    print "calling SNPs from %s" % bam
#    # runStageCheck('callSNPs', flag_file, fasta_reference, bam, ref_files['exon_bed_extended'], ref_files['dbsnp'], logFile, output_vcf)
#    runStageCheck('callSNPs', flag_file, fasta_reference, bam, ref_files['exon_bed'], ref_files['dbsnp'], logFile, output_vcf)

@follows(indexBaseQualRecalBam)
@transform(baseQualRecalPrint,
            regex(r'(.*?)([^/]+)\.recal\.bam'),
            [r'%s/\2.recal.hap.raw.vcf' % variant_dir,
             r'%s/\2.recal.hap.raw.vcf.idx' % variant_dir,
             r'%s/\2.callHAPs.Success' % variant_dir])
def callHAP(inputs, outputs):
    """
    Use GATK HaplotypeCaller to call SNPs/Indels from recalibrated bams.
    """
    bam, _success = inputs
    output_vcf, _idx, flag_file = outputs
    logFile = mkLogFile(logDir, bam, '.callHAPs.log')
    #print "calling haplotypes from %s" % bam
    runStageCheck('callHAP', flag_file, fasta_reference, bam, ref_files['dbsnp'], logFile, output_vcf)
    #runStageCheck('callHAP', flag_file, fasta_reference, bam, ref_files['exon_bed'], ref_files['dbsnp'], logFile, output_vcf)

'''
@follows(indexBaseQualRecalBam)
@merge(baseQualRecalPrint,
             [r'%s/All.hapmerged.raw.vcf' % variant_dir,
              r'%s/All.hapmerged.raw.vcf.idx' % variant_dir,
              r'%s/All.callHAPMerged.Success' % variant_dir])
             # regex(r'(.*?)([^/]+)\.recal\.bam'),
             # [r'%s/\2.HAP.merged.raw.vcf' % variant_dir,
             #  r'%s/\2.HAP.merged.raw.vcf.idx' % variant_dir,
             #  r'%s/\2.callHAPMerged.Success' % variant_dir])
def callHAPMerged(inputs, outputs):
     """
     Use GATK HaplotypeCaller to call SNPs/Indels from recalibrated bams.
     """
     bams = ""
     success = []
     for i in inputs:
         bams += " -I " + i[0] + " "
         success.append(i[1])
     # bams = [" -I " + i for i in temp_bam]
     output_vcf, _idx, flag_file = outputs
     logFile = mkLogFile(logDir, bams, '.callHAPMerged.log')
     print "calling haplotypes from %s" % bams
     runStageCheck('callHAPMerged', flag_file, fasta_reference, bams, ref_files['dbsnp'], logFile, output_vcf)
     #runStageCheck('callHAPMerged', flag_file, fasta_reference, bams, ref_files['exon_bed_extended'], ref_files['dbsnp'], logFile, output_vcf)
'''
'''
# CALL VARSCAN SAMPLE
@follows(indexBaseQualRecalBam)
@transform(baseQualRecalPrint,
             regex(r'(.*?)([^/]+)\.recal\.bam'),
             [r'%s/\2.recal.varscan.raw.vcf' % variant_dir,
             # r'%s/\2.recal.varscan.raw.vcf.idx' % variant_dir,
             r'%s/\2.recal.varscan.callVARSCAN.Success' % variant_dir])
def callVARSCAN(inputs, outputs):
     """
     Use GATK VARSCAN to call SNPs/Indels from recalibrated bams.
     """
     bam, _success = inputs
     # output_vcf, _idx, flag_file = outputs
     output_vcf, flag_file = outputs
     pathtobam, bamfilename = os.path.split(bam)
     sample_file = varscan_temp_dir + "/" + bamfilename.split(".")[0]+"_SAMPLENAME"
     with open(sample_file, "w") as text_file:
         text_file.write("{}".format(bamfilename.split(".")[0]))
     logFile = mkLogFile(logDir, bam, '.callVARSCAN.log')
     # LOG print "calling varscan from %s" % bam
     runStageCheck('callVARSCAN', flag_file, fasta_reference, working_files['varscanjar'], sample_file, output_vcf, bam)
'''
# @follows(indexBaseQualRecalBam)
# @merge(baseQualRecalPrint,
#             [r'%s/All.varscanmerged.raw.vcf' % variant_dir,
#              # r'%s/All.varscanmerged.raw.vcf.idx' % variant_dir,
#              r'%s/All.callVARSCANMerged.Success' % variant_dir])
# def callVARSCANMerged(inputs, outputs):
#     """
#     Use GATK HaplotypeCaller to call SNPs/Indels from recalibrated bams.
#     """
#     samples = []
#     bams = ""
#     success = []
#     sample_file_name = varscan_temp_dir + "/" + "ALL_SAMPLENAME"
#     sample_file = open(sample_file_name, "w")
#     for i in inputs:
#         bams += i[0] + " "
#         success.append(i[1])
#         pathtobam, bamfilename = os.path.split(i[0])
#         # with open(sample_file_name, "w") as text_file:
#             # text_file.write("{}".format(bamfilename.split(".")[0]))
#         sample_file.write("{}".format(bamfilename.split(".")[0]))
#         sample_file.write("{}".format("\n"))
#         # write to file <- bamfilename.split(".")[0]
#     sample_file.close()
#     # print "BAMS = " + bams
#     # bams = [" -I " + i for i in temp_bam]
#     # output_vcf, _idx, flag_file = outputs
#     output_vcf, flag_file = outputs
#     logFile = mkLogFile(logDir, bams, '.VARSCANMerged.log')
#     print "calling varscan merged from %s" % bams
#     runStageCheck('callVARSCANMerged', flag_file, fasta_reference, working_files['varscanjar'], sample_file_name, output_vcf, bams)
#     # runStageCheck('callVARSCANMerged', flag_file, fasta_reference, bams, working_files['varscanjar'], output_vcf)

# os.system('mv *_SAMPLENAME tempdir/.')

# @follows(indexBaseQualRecalBam)
# @transform(baseQualRecalPrint,
#             regex(r'(.*?)([^/]+)\.recal\.bam'),
#             [r'%s/\2.recal.UGT.INDEL.vcf' % variant_dir,
#              r'%s/\2.recal.UGT.INDEL.vcf.idx' % variant_dir,
#              r'%s/\2.recal.UGT.callIndels.Success' % variant_dir])
# def callIndels(inputs, outputs):
#     """
#     Use GATK UnifiedGenotyper to call indels from recalibrated bams.
#     """
#     bam, _success = inputs
#     output_vcf, _idx, flag_file = outputs
#     logFile = mkLogFile(logDir, bam, '.callIndels.log')
#     print "calling Indels from %s" % bam
#     runStageCheck('callIndels', flag_file, fasta_reference, bam, ref_files['exon_bed'], ref_files['dbsnp'], logFile, output_vcf)

# # add VQSR here
# regex(r'(.*?)([^/]+)\.recal\.bam')
@transform(callHAP, suffix('recal.hap.raw.vcf'),
            ['recal.hap.raw.recal', 'recal.hap.raw.tranches', 'recal.hap.raw.recal.R' ,'VarRecal.Success'])
def callVariantRecalibrator(inputs, outputs):
    """
    Use GATK VariantFiltration to filter raw SNP calls.
    """
    input_vcf, _idx, _success = inputs
    output_recal, output_tranches, output_R, flag_file = outputs
    logFile = mkLogFile(logDir, input_vcf, '.VarRecal.log')
    print "VariantRecalibrator -> %s" % input_vcf
    runStageCheck('callVariantRecalibrator', flag_file, fasta_reference, input_vcf, ref_files['hapmap'],ref_files['omnimap'], ref_files['1kghc'], ref_files['dbsnp'], output_recal, output_tranches, output_R, logFile)

#recal.hap.raw.vcf
@transform([callVariantRecalibrator, callHAP], regex(r'(.*?)([^/]+)\.recal\.hap\.raw\.recal'),
            add_inputs(r'\1\2.recal.hap.raw.vcf'),
            # add_inputs([r'\1\2.recal.hap.raw.vcf']), add_inputs([r'\1\2.recal.hap.raw.tranches']),
            [r'\1\2.recal.hap.raw.recal.vcf', r'\1\2.ApplyVarRecal.Success'])
def callApplyRecalibration(inputs, outputs):
    """
    Use GATK VariantFiltration to filter raw SNP calls.
    """
    #[[input_recal, _VarRecal_success], [input_vcf], [input_recal] , [input_trances]] = inputs
    [input_recal, input_tranches, input_R, _VarRecal_success], input_vcf = inputs
    output_vcf, flag_file = outputs
    logFile = mkLogFile(logDir, input_recal, '.ApplyVarRecal.log')
    print "ApplyRecalibration SNP -> %s" % input_vcf
    runStageCheck('callApplyRecalibration', flag_file, fasta_reference, input_vcf, input_recal, input_tranches, output_vcf, logFile)

# # add VQSR here

#@transform(callSNPs, suffix('.SNP.vcf'),
#            ['.SNP.filtered.vcf', '.SNP.filtered.vcf.idx', '.filterSNPs.Success'])
#def filterSNPs(inputs, outputs):
#    """
#    Use GATK VariantFiltration to filter raw SNP calls.
#    """
#    input_vcf, _idx, _success = inputs
#    output_vcf, _idxout, flag_file = outputs
#    logFile = mkLogFile(logDir, input_vcf, '.filterSNPs.log')
#    print "filtering SNPs from %s" % input_vcf
#    runStageCheck('filterSNPs', flag_file, fasta_reference, input_vcf, logFile, output_vcf)

# @transform(callVARSCAN, suffix('.varscan.raw.vcf'),
#             ['.varscan.filtered.vcf', '.varscan.filtered.vcf.idx', '.filterVarscan.Success'])
# def filterVarscan(inputs, outputs):
#     """
#     Use GATK VariantFiltration to filter raw SNP calls.
#     """
#     input_vcf, _idx, _success = inputs
#     print input_vcf
#     output_vcf, _idxout, flag_file = outputs
#     logFile = mkLogFile(logDir, input_vcf, '.filterVarscan.log')
#     print "filtering varcsan SNPs from %s" % input_vcf
#     runStageCheck('filterVarscan', flag_file, fasta_reference, input_vcf, logFile, output_vcf)

'''
# Use Varscan filter
@transform(callVARSCAN, suffix('.varscan.raw.vcf'),
            ['.varscan.filtered.vcf', '.filterVarscan.Success'])
def filterVarscan(inputs, outputs):
    """
    Use Varscan filter to filter raw SNP calls.
    """
    # input_vcf, _idx, _success = inputs
    input_vcf, _success = inputs
    # print input_vcf
    # output_vcf, _idxout, flag_file = outputs
    output_vcf, flag_file = outputs
    logFile = mkLogFile(logDir, input_vcf, '.filterVarscan.log')
    # print "filtering varcsan SNPs (using varscan filter) from %s" % input_vcf
    runStageCheck('filterVarscan', flag_file, fasta_reference, working_files['varscanjar'], input_vcf, output_vcf)
    #runStageCheck('callVARSCAN', flag_file, fasta_reference, working_files['varscanjar'], sample_file, output_vcf, bam)
# sys.exit(0)
'''

@transform(callHAP, suffix('.raw.vcf'),
            ['.filtered.vcf', '.filtered.vcf.idx', '.filterHapVcfs.Success'])
def filterHapVcfs(inputs, outputs):
    """
    Use GATK VariantFiltration to filter raw sample HAP calls.
    """
    input_vcf, _idx, _success = inputs
    output_vcf, _idxout, flag_file = outputs
    logFile = mkLogFile(logDir, input_vcf, '.filterSNPs.log')
    # print "filtering haplotyper vcf from %s" % input_vcf
    runStageCheck('filterHapVcfs', flag_file, fasta_reference, input_vcf, logFile, output_vcf)

@transform(callApplyRecalibration, suffix('.raw.recal.vcf'),
            ['.recal.filtered.vcf', '.recal.filtered.vcf.idx', '.recal.filterHapVcfs.Success'])
def filterHapVcfsAfterVQSR(inputs, outputs):
    """
    Use GATK VariantFiltration to filter raw sample HAP calls.
    """
    #input_vcf, _idx, _success = inputs
    input_vcf, _success = inputs
    output_vcf, _idxout, flag_file = outputs
    logFile = mkLogFile(logDir, input_vcf, '.recal.filterSNPs.log')
    # print "filtering haplotyper vcf from %s" % input_vcf
    runStageCheck('filterHapVcfs', flag_file, fasta_reference, input_vcf, logFile, output_vcf)

'''
@transform(callHAPMerged, suffix('.raw.vcf'),
             ['.filtered.vcf', '.filtered.vcf.idx', '.filterHapVcfs.Success'])
def filterMergedHapVcfs(inputs, outputs):
     """
     Use GATK VariantFiltration to filter raw merged HAP call.
     """
     input_vcf, _idx, _success = inputs
     output_vcf, _idxout, flag_file = outputs
     logFile = mkLogFile(logDir, input_vcf, '.filterSNPs.log')
     print "filtering haplotyper vcf from %s" % input_vcf
     runStageCheck('filterMergedHapVcfs', flag_file, fasta_reference, input_vcf, logFile, output_vcf)
'''

# @transform(callIndels, suffix('.INDEL.vcf'),
#             ['.INDEL.filtered.vcf', '.INDEL.filtered.vcf.idx', '.filterIndels.Success'])
# def filterIndels(inputs, outputs):
#     """
#     Use GATK VariantFiltration to filter raw INDEL calls.
#     """
#     input_vcf, _idx, _success = inputs
#     output_vcf, _idxout, flag_file = outputs
#     logFile = mkLogFile(logDir, input_vcf, '.filterIndels.log')
#     print "filtering indels from %s" % input_vcf
#     runStageCheck('filterIndels', flag_file, fasta_reference, input_vcf, logFile, output_vcf)

# @transform([filterSNPs, filterIndels], regex(r'.*?([^/]+)\.vcf'),
#     [r'%s/\1.ensembl.vcf' % ensembl_dir,r'%s/\1.getEnsemblAnnotations.Success' % ensembl_dir])
# def getEnsemblAnnotations(inputs, outputs):
#     """
#     Annotate vcf using ENSEMBL variant effect predictor.
#     """
#     vcf, _idx, _success = inputs
#     output, flag_file = outputs
#     logFile = mkLogFile(logDir, vcf, '.EnsemblAnnotation.log')
#     print "Annotating %s with ENSEMBL variant effect predictor" % os.path.basename(vcf)
#     runStageCheck('annotateEnsembl', flag_file, vcf, output, logFile)

# @transform([filterSNPs, filterIndels, filterVarscan, filterHapVcfs], regex(r'.*?([^/]+)\.vcf'),
# @transform([filterSNPs, filterIndels, filterHapVcfs], regex(r'.*?([^/]+)\.vcf'),
# @transform([filterIndels, filterHapVcfs], regex(r'.*?([^/]+)\.vcf'),
#@transform([filterHapVcfs, filterVarscan], regex(r'.*?([^/]+)\.vcf'),
@transform([filterHapVcfs, filterHapVcfsAfterVQSR], regex(r'.*?([^/]+)\.vcf'),
    [r'%s/\1.snpeff.vcf' % snpeff_dir,r'%s/\1.getSnpeffAnnotations.Success' % snpeff_dir])
def getSnpeffAnnotations(inputs, outputs):
    """
    Annotate vcf using snpeff variant effect predictor.
    """
    print "Inputs: %s" % inputs
    vcf, _idx, _success = inputs
    #vcf, _success = inputs
    output, flag_file = outputs
    logFile = mkLogFile(logDir, vcf, '.snpEffAnnotation.log')
    config = working_files['snpeff'] + "snpEff.config"
    # print "Annotating %s with snpeff variant effect predictor" % os.path.basename(vcf)
    runStageCheck('annotateSNPEff', flag_file, working_files['snpeff'], config, vcf, output)
    # runStageCheck('annotateSNPEff', flag_file, vcf, output, logFile)

# FILTER HQ VARIANTS
@transform(getSnpeffAnnotations, suffix('.merged.dedup.realigned.leftAligned.recal.hap.recal.filtered.snpeff.vcf'),
            ['.hq.vcf', '.filterHQVar.Success'])
            #['.hq.vcf', '.hq.vcf.idx', '.filterHQVar.Success'])
def filterSnpSift(inputs, outputs):
    """
    Filter Recalibrated variants.
    """
    #print "Inputs: %s" % inputs
    #print "Inputs: %s" % outputs
    vcf, _success = inputs
    #output, _idxout, flag_file = outputs
    output, flag_file = outputs
    logFile = mkLogFile(logDir, vcf, '.snpSiftFilter.log')
    # config = working_files['snpeff'] + "snpEff.config"
    # print "Annotating %s with snpeff variant effect predictor" % os.path.basename(vcf)
    #runStageCheck('filterSnpSift', flag_file, working_files['snpeff'], config, vcf, output)
    runStageCheck('filterSnpSift', flag_file, working_files['snpeff'], vcf, output)

#@transform(filterHapVcfs, suffix('.vcf'),
@transform(filterSnpSift, suffix('.vcf'),
            ['.vcf.gz', '.vcf.gz.tbi', '.vcfindexSNPs.Success'])
def vcfIndexSNPs(inputs, outputs):
    """
    Use bgzip and tabix to prepare raw SNPs vcf for vcftools handling.
    """
    #vcf, _idx, _success = inputs
    vcf, _success = inputs
    zipfile, tabix_index, flag_file = outputs
    print "bgzip and tabix (for vcftools) on %s" % vcf
    runStageCheck('indexVCF', flag_file, vcf)

# MERGE ALL SAMPLE VCF FILES TO A SINGLE MERGED.VCF
#@collate(filterSnpSift, regex(r".*.hq.vcf"), r'%s/\1.snpeff.vcf' % snpeff_dir
@collate(vcfIndexSNPs, regex(r".*.hq.vcf.gz"),
            [r'%s/merged.vcf' % snpeff_dir, r'%s/mergeVCFs.Success' % snpeff_dir])
def mergeVCFs(inputs, outputs):
    """
    Merge the sorted bams together for each sample.
    Picard should cope correctly if there is only one input.
    """
    vcfs = [vcf for [vcf, _idx, _success] in inputs]
    output_vcf, flag_file = outputs
    vcfinputs = ' '.join([" %s" % vcf for vcf in vcfs])
    print "merging %s into %s" % (",".join([os.path.basename(vcf) for vcf in vcfs]), os.path.basename(output_vcf))
    runStageCheck('mergeVCFs', flag_file, vcfinputs, output_vcf)

#@transform([filterVarscan], regex(r'.*?([^/]+)\.vcf'),
#    [r'%s/\1.snpeff.vcf' % snpeff_dir,r'%s/\1.geitSnpeffAnnotations.Success' % snpeff_dir])
#def getSnpeffAnnotationsVarscan(inputs, outputs):
#    """
#    Annotate vcf using snpeff variant effect predictor.
#    """
#    print "Inputs: %s" % inputs
#    #vcf, _idx, _success = inputs
#    vcf, _success = inputs
#    output, flag_file = outputs
#    logFile = mkLogFile(logDir, vcf, '.snpEffAnnotation.log')
#    config = working_files['snpeff'] + "snpEff.config"
#    # print "Annotating %s with snpeff variant effect predictor" % os.path.basename(vcf)
#    runStageCheck('annotateSNPEff', flag_file, working_files['snpeff'], config, vcf, output)
#    # runStageCheck('annotateSNPEff', flag_file, vcf, output, logFile)

# # Indexing steps

@transform(mergeBams, suffix('.bam'),
            ['.bam.bai', '.bam.indexMergedBams.Success'])
def indexMergedBams(inputs, outputs):
    """
    Index the merged bams using samtools.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "samtools index on %s" % os.path.basename(bam)
    runStageCheck('indexBam', flag_file, bam)

@transform(dedup, suffix('.bam'),
            ['.bam.bai', '.bam.indexDedupedBams.Success'])
def indexDedupedBams(inputs, outputs):
    """
    Index the de-duplicated bams using samtools. Note that this actually goes from the fixMate-ed bams.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "samtools index on %s" % os.path.basename(bam)
    runStageCheck('indexBam', flag_file, bam)

@transform(realign, suffix('.bam'),
            ['.bam.bai', '.bam.indexRealignedBams.Success'])
def indexRealignedBams(inputs, outputs):
    """
    Index the locally realigned bams using samtools.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "samtools index on %s" % os.path.basename(bam)
    runStageCheck('indexBam', flag_file, bam)

@transform(mergeBams, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountMergedBams.Success'])
def igvcountMergedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the merged bam files, to improve viewing of the bam coverage in igv.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

@transform(realign, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountRealignedBams.Success'])
def igvcountRealignedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the merged bam files, to improve viewing of the bam coverage in igv.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

@transform(dedup, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountDedupedBams.Success'])
def igvcountDedupedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the deduped bam files, to improve viewing of the bam coverage in igv. Note that this actually goes from the fixMate-ed bams.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

@transform(baseQualRecalPrint, suffix('.bam'),
            ['.bam.tdf', '.bam.igvcountRecalibratedBams.Success'])
def igvcountRecalibratedBams(inputs, outputs):
    """
    Use igvtools count to create a .tdf file for the recalibrated bam files, to improve viewing of the bam coverage in igv.
    """
    bam, _success = inputs
    outfile, flag_file = outputs
    print "igvtools count on %s" % os.path.basename(bam)
    runStageCheck('igvcount', flag_file, bam, outfile)

'''
#@transform(filterHapVcfs, suffix('.vcf'),
@transform(filterSnpSift, suffix('.vcf'),
            ['.vcf.gz', '.vcf.gz.tbi', '.vcfindexSNPs.Success'])
def vcfIndexSNPs(inputs, outputs):
    """
    Use bgzip and tabix to prepare raw SNPs vcf for vcftools handling.
    """
    vcf, _idx, _success = inputs
    zipfile, tabix_index, flag_file = outputs
    print "bgzip and tabix (for vcftools) on %s" % vcf
    runStageCheck('indexVCF', flag_file, vcf)
'''

#@transform(filterIndels, suffix('.vcf'),
#            ['.vcf.gz', '.vcf.gz.tbi', '.vcfindexIndels.Success'])
#def vcfIndexIndels(inputs, outputs):
#    """
#    Use bgzip and tabix to prepare raw indels vcf for vcftools handling.
#    """
#    vcf, _idx, _success = inputs
#    zipfile, tabix_index, flag_file = outputs
#    print "bgzip and tabix (for vcftools) on %s" % vcf
#    runStageCheck('indexVCF', flag_file, vcf)

# Coverage steps
#@follows(indexMergedBams)
#@transform(mergeBams,
#            regex(r'(.*?)([^/]+)\.bam'),
#            [r'%s/\2.early.DepthOfCoverage.sample_cumulative_coverage_counts' % coverage_dir,
#            r'%s/\2.early.DepthOfCoverage.sample_cumulative_coverage_proportions' % coverage_dir,
##            r'%s/\2.early.DepthOfCoverage.sample_gene_summary' % coverage_dir,
#            r'%s/\2.early.DepthOfCoverage.sample_interval_statistics' % coverage_dir,
#            r'%s/\2.early.DepthOfCoverage.sample_interval_summary' % coverage_dir,
#            r'%s/\2.early.DepthOfCoverage.sample_statistics' % coverage_dir,
#            r'%s/\2.early.DepthOfCoverage.sample_summary' % coverage_dir,
#            r'%s/\2.earlyDepthOfCoverage.Success' % coverage_dir])
#def earlyDepthOfCoverage(inputs, outputs):
#    """
#    Use GATK DepthOfCoverage to get a first pass at coverage statistics, after merging bams.
#    """
#    bam, _success = inputs
#    flag_file = outputs[-1]
#    output_example = outputs[0]
#    output_base = os.path.splitext(output_example)[0]
#    print "calculating coverage statistics using GATK DepthOfCoverage on %s" % bam
#    runStageCheck('depthOfCoverage', flag_file, fasta_reference, bam, ref_files['exon_bed'], output_base)

#@follows(indexDedupedBams)
#@transform(dedup,
#            regex(r'(.*?)([^/]+)\.dedup\.bam'),
#            [r'%s/\2.deduped.DepthOfCoverage.sample_cumulative_coverage_counts' % coverage_dir,
#            r'%s/\2.deduped.DepthOfCoverage.sample_cumulative_coverage_proportions' % coverage_dir,
##            r'%s/\2.early.DepthOfCoverage.sample_gene_summary' % coverage_dir,
#            r'%s/\2.deduped.DepthOfCoverage.sample_interval_statistics' % coverage_dir,
#            r'%s/\2.deduped.DepthOfCoverage.sample_interval_summary' % coverage_dir,
#            r'%s/\2.deduped.DepthOfCoverage.sample_statistics' % coverage_dir,
#            r'%s/\2.deduped.DepthOfCoverage.sample_summary' % coverage_dir,
#            r'%s/\2.dedupedDepthOfCoverage.Success' % coverage_dir])
#def dedupedDepthOfCoverage(inputs, outputs):
#    """
#    Use GATK DepthOfCoverage to get a coverage statistics as soon as duplicates are removed.
#    """
#    bam, _success = inputs
#    flag_file = outputs[-1]
#    output_example = outputs[0]
#    output_base = os.path.splitext(output_example)[0]
#    print "calculating coverage statistics using GATK DepthOfCoverage on %s" % bam
#    runStageCheck('depthOfCoverage', flag_file, fasta_reference, bam, ref_files['exon_bed'], output_base)

@follows(indexBaseQualRecalBam)
@transform(baseQualRecalPrint,
            regex(r'(.*?)([^/]+)\.recal\.bam'),
            [r'%s/\2.DepthOfCoverage.sample_cumulative_coverage_counts' % coverage_dir,
            r'%s/\2.DepthOfCoverage.sample_cumulative_coverage_proportions' % coverage_dir,
#            r'%s/\2.DepthOfCoverage.sample_gene_summary' % coverage_dir,
            r'%s/\2.DepthOfCoverage.sample_interval_statistics' % coverage_dir,
            r'%s/\2.DepthOfCoverage.sample_interval_summary' % coverage_dir,
            r'%s/\2.DepthOfCoverage.sample_statistics' % coverage_dir,
            r'%s/\2.DepthOfCoverage.sample_summary' % coverage_dir,
            r'%s/\2.depthOfCoverage.Success' % coverage_dir])
def finalDepthOfCoverage(inputs, outputs):
    """
    Use GATK DepthOfCoverage to get coverage statistics.
    """
    bam, _success = inputs
    flag_file = outputs[-1]
    output_example = outputs[0]
    output_base = os.path.splitext(output_example)[0]
    print "calculating coverage statistics using GATK DepthOfCoverage on %s" % bam
    runStageCheck('depthOfCoverage', flag_file, fasta_reference, bam, output_base)
    #runStageCheck('depthOfCoverage', flag_file, fasta_reference, bam, ref_files['exon_bed'], output_base)

'''
@follows(indexBaseQualRecalBam)
@transform(baseQualRecalPrint,
            regex(r'(.*?)([^/]+)\.recal\.bam'),
            [r'%s/\2.exon_coverage.txt' % coverage_dir,
            r'%s/\2.exonCoverage.Success' % coverage_dir])
def exonCoverage(inputs, outputs):
    """
    Use bedtools' coverageBed to get coverage for each exon.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "calculating exon coverage using coverageBed on %s" % bam
    runStageCheck('exonCoverage', flag_file, bam, ref_files['exon_bed'], output)
'''

# Read-counting steps

#@transform(dedup, suffix('.bam'),
#            ['.exon_intersections.bam', '.intersectBamExons.Success'])
#def intersectBamExons(inputs, outputs):
#    """
#    Use bedtools' intersectBed on deduped BAM file to get the reads which
#    map to exon capture array regions.
#    """
#    bam, _success = inputs
#    output, flag_file = outputs
#    print "Calculating intersection with exon regions (%s) of %s" % (ref_files['exon_bed'], bam)
#    runStageCheck('intersectBam', flag_file, bam, ref_files['exon_bed'], output)

#@transform(dedup, suffix('.bam'),
#            ['.exon_extended_intersections.bam', '.bam.intersectBamExonsExtended.Success'])
#def intersectBamExonsExtended(inputs, outputs):
#    """
#    Use bedtools' intersectBed on deduped BAM file to get the reads which
#    map to exon capture array regions.
#    """
#    bam, _success = inputs
#    output, flag_file = outputs
#    print "Calculating intersection with exon + buffer (%s) regions of %s" % (ref_files['exon_bed'], bam)
#    runStageCheck('intersectBam', flag_file, bam, ref_files['exon_bed'], output)

@transform(samToBam, suffix('.bam'),
            ['.bam.flagstat', '.bam.countRunBam.Success'])
def countRunBam(inputs, outputs):
    """
    Run samtools flagstat on the initial per-lane, per-run bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)

@transform(mergeBams, suffix('.bam'),
            ['.bam.flagstat', '.bam.countRunBam.Success'])
def countMergedBam(inputs, outputs):
    """
    Run samtools flagstat on the merged bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)

@transform(realign, suffix('.bam'),
            ['.bam.flagstat', '.bam.countRealignedBam.Success'])
def countRealignedBam(inputs, outputs):
    """
    Run samtools flagstat on the realigned bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)

@transform(dedup, suffix('.bam'),
            ['.bam.flagstat', '.bam.countDedupedBam.Success'])
def countDedupedBam(inputs, outputs):
    """
    Run samtools flagstat on the deduped bam file.
    """
    bam, _success = inputs
    output, flag_file = outputs
    print "Running samtools flagstat on %s" % bam
    runStageCheck('flagstat', flag_file, bam, output)

#@transform(intersectBamExons, suffix('.bam'),
#            ['.bam.flagstat', '.bam.countExonsBam.Success'])
#def countExonsBam(inputs, outputs):
#    """
#    Run samtools flagstat on the exon-regions-only bam file.
#    """
#    bam, _success = inputs
#    output, flag_file = outputs
#    print "Running samtools flagstat on %s" % bam
#    runStageCheck('flagstat', flag_file, bam, output)

#@transform(intersectBamExonsExtended, suffix('.bam'),
#            ['.bam.flagstat', '.bam.countExonsExtendedBam.Success'])
#def countExonsExtendedBam(inputs, outputs):
#    """
#    Run samtools flagstat on the extended-exon-regions-only bam file.
#    """
#    bam, _success = inputs
#    output, flag_file = outputs
#    print "Running samtools flagstat on %s" % bam
#    runStageCheck('flagstat', flag_file, bam, output)

# Data collation and plotting steps

#@merge([countExonsBam, countExonsExtendedBam, countDedupedBam, countMergedBam],
#@merge([countDedupedBam, countMergedBam],
#        ["%s/readcounts.txt" % results_dir, "%s/readcount_fractions.txt" % results_dir, "%s/collateReadcounts.Success" % results_dir])
#def collateReadCounts(inputs, outputs):
#    """
#    Run samtools flagstat on the extended-exon-regions-only bam file.
#    """
#    # Note expected input and output directories are effectively hard-coded
#    in_dir =  sambam_dir
#    out_dir = results_dir
#    flag_file = outputs[-1]
#    print "Collating read counts"
#    runStageCheck('collateReadcounts', flag_file, in_dir, out_dir)

