#!/bin/env python

import sys
import re
import os.path
import os
import string
from glob import *

working_files = {
    'fastq_dirs': ['/vlsci/VR0182/shared/Data/Lobb/fastq/'],
    'fastq_symlink_dir': '/vlsci/VR0182/shared/km/analysis/output/fastq_symlinks',
    'output_dir': '/vlsci/VR0182/shared/km/analysis/output/',
    'varscanjar': '/vlsci/VR0002/kmahmood/Programs/varscan/VarScan.v2.3.6.jar',
    'snpeff': '/vlsci/VR0002/kmahmood/Programs/snpEff3.4/snpEff/'
}

original_fastq_files = []
for fastq_dir in working_files['fastq_dirs']:
    original_fastq_files += glob(os.path.join(fastq_dir, '*/*_[1|2].fastq'))

#print original_fastq_files

for file in original_fastq_files:
    name = os.path.basename(file)
    pu,remaining = name.split("_")
    barcode = pu.split(".")[0]
    lane = pu.split(".")[1]
    id = pu[:5] + "." + lane
    sm = os.path.basename(os.path.dirname(file))
    sm = string.replace(sm,"_","-")
    sm_p = string.replace(sm,"-","_")

    #print "sm = " + sm + "\npu = " + pu + "\nlane = " + lane + "\nid = " + id + "\nbarcode = " + barcode

    #file = sm + "_" + name
    #file2 = sm + "_" + barcode + ".L" + lane + "_" + remaining

    #
    #file2 = "/vlsci/VR0182/shared/Data/Lobb/fastq_test/0015_99_002/" + barcode + "L" + lane + "_" + sm + "_" + remaining + ".gz"
    # OUTPUT TO SCRATCH
    file2 = "/scratch/VR0182/lobb/fastqz/" + sm_p + "/" + barcode + "L" + lane + "_" + sm + "_" + remaining + ".gz"

    #print "ln -s " + file + " " + file2
    #print "mkdir /scratch/VR0182/lobb/fastqz/" + sm_p
    # TEMP print "cat " + file + " | gzip -9 > " + file2
    print "in=\"" + file + "\"" + " out=\"" + file2 + "\"" + " sbatch --job-name=" + sm + " FastqtoGz.sh"


