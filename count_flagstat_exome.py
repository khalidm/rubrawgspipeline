
"""
Read all output flagstat files and store values in a text table.
Creates two tables in the output directory: readcounts.txt and readcounts_fractions.txt, which contains useful human-readable percentage statistics.
Each row will be a sample.
Values to store (columns) are:
    Total reads from SAMPLE.bam.flagstat
    Mapped reads from SAMPLE.bam.flagstat
    Mapped reads from SAMPLE.dedup.bam.flagstat
    Mapped reads from SAMPLE.dedup.exon_intersections.bam.flagstat
    Mapped reads from SAMPLE.dedup.exon_extended_intersections.bam.flagstat
There is an assumption that the bams were mapped with something like bwa, which includes unmapped reads.

Usage: python count_flagstat_wgs.py flagstat_directory output_directory
"""

import sys
import os
import re
import optparse
from collections import defaultdict

class FlagstatParseException (Exception):
    pass

def read_flagstat(filename):
    """
    Given a filename, parse the flagstat values and return as a hash.
    Relying on flagstat contents being in usual order.
    """
    numbers = re.compile(r'^(\d+)\s+\+\s+(\d+)\s+')
    values = {}
    f = open(filename)
    for field in ['total',
                    'duplicates',
                    'mapped',
                    'paired',
                    'read1',
                    'read2',
                    'properly_paired',
                    'both_mapped',
                    'singletons',
                    'mate_distant',
                    'mate_distant_goodqual']:
        line = f.readline().strip()
        match = numbers.match(line)
        if not match:
            raise FlagstatParseException
        # print field + "\t" + match.group(1) + "\t" + match.group(2)
        values[field] = int(match.group(1))
        values[field+'_QCfailed'] = int(match.group(2))
    f.close()
    return values

# -----
def printplus(obj):
    """
    Pretty-prints the object passed in.

    """
    print "\t\t-- "
    # Dict
    if isinstance(obj, dict):
        for k, v in sorted(obj.items()):
            print u'{0}: {1}'.format(k, v)

    # List or tuple            
    elif isinstance(obj, list) or isinstance(obj, tuple):
        for x in obj:
            print x

    # Other
    else:
        print obj
# -----

# Get arguments and input filenames
parser = optparse.OptionParser(usage=__doc__)
#parser.add_option()
(options, args) = parser.parse_args()
if len(args) != 2:
    parser.error("Wrong number of arguments - see usage info")
in_dir = args[0]
out_dir = args[1]
if not (os.path.exists(out_dir) and os.path.isdir(out_dir)):
    sys.exit("There does not seem to be a directory %s , exiting" % out_dir)

filenames = os.listdir(in_dir)

#print ', '.join(filenames)

alignedname = re.compile('^([^_\.]+).merged.bam.flagstat')
dedupname = re.compile('^([^_\.]+).merged.dedup.bam.flagstat')
exonintname = re.compile('^([^_\.]+).merged.dedup.exon_intersections.bam.flagstat')
exonplusname = re.compile('^([^_\.]+).merged.dedup.exon_extended_intersections.bam.flagstat')

samples = defaultdict(dict)
for filename in filenames:    
    if exonplusname.match(filename):        
        match = exonplusname.match(filename)
        name = match.group(1)
        # print "filename " + filename + " - " + name
        values = read_flagstat( os.path.join(in_dir, filename) )
        samples[name]['exonplus'] = values['mapped']
        # printplus(values)
    elif exonintname.match(filename):        
        match = exonintname.match(filename)
        name = match.group(1)
        # print "filename " + filename + " - " + name
        values = read_flagstat( os.path.join(in_dir, filename) )
        samples[name]['exon'] = values['mapped']    
        # printplus(values)
    elif dedupname.match(filename):        
        match = dedupname.match(filename)
        name = match.group(1)
        # print "filename " + filename + " - " + name
        values = read_flagstat( os.path.join(in_dir, filename) )
        samples[name]['deduped'] = values['mapped']            
        # printplus(values)
    elif alignedname.match(filename):        
        match = alignedname.match(filename)
        name = match.group(1)
        # print "filename " + filename + " - " + name
        values = read_flagstat( os.path.join(in_dir, filename) )
        samples[name]['mapped'] = values['mapped']
        samples[name]['total'] = values['total']
        # printplus(values)
        
# print ', '.join(samples.keys())

tablefile = os.path.join(out_dir, "readcounts.txt")
tablefile_plus = os.path.join(out_dir, "readcounts_fractions.txt")

OUT_TABLE = open(tablefile, 'w')
OUT_TABLEPLUS = open(tablefile_plus, 'w')

OUT_TABLE.write("Sample\tTotal\tMapped\tDeduped\tOn-target\tNear-target\n") 
OUT_TABLEPLUS.write("Sample\tTotal\tMapped\t%\tDeduped\t%\tOn-target\t%\tNear-target\t%\n")

# debug
# for sample in sorted(samples.keys()):
#     values = samples[sample]
#     printplus(values)
#     fraction_mapped = float(values['mapped'])/values['total']
#     fraction_deduped = float(values['deduped'])/values['mapped']
#     fraction_target = float(values['exon'])/values['deduped']  
#     print "values length = " + str(len(values)) + " - " + sample + " %mapped " + str(fraction_mapped) + " %deduped " + str(fraction_deduped) + " %fraction_target " + str(fraction_target)
# exit(0)
# debug

for sample in sorted(samples.keys()):
    values = samples[sample]
    fraction_mapped = float(values['mapped'])/values['total']    
    # if int(values['deduped']) > 0:
    fraction_deduped = float(values['deduped'])/values['mapped']
    fraction_target = float(values['exon'])/values['deduped']
    fraction_targetplus = float(values['exonplus'])/values['deduped']
    # else:
    #     fraction_deduped = 0
    #     fraction_target = 0
    #     fraction_targetplus = 0

    OUT_TABLE.write( "%s\t%d\t%d\t%d\t%d\t%d\n" % (sample, values['total'], values['mapped'], values['deduped'], values['exon'], values['exonplus']) )
    # OUT_TABLE.write( "%s\t%d\t%d\n" % (sample, values['total'], values['mapped']) )
    
    OUT_TABLEPLUS.write( "%s\t%d\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\t%d\t%.3f\n" % (sample, values['total'], values['mapped'], fraction_mapped, values['deduped'], fraction_deduped, values['exon'], fraction_target, values['exonplus'], fraction_targetplus) )
    # OUT_TABLEPLUS.write( "%s\t%d\t%d\t%.3f\n" % (sample, values['total'], values['mapped'], fraction_mapped) )

OUT_TABLE.close()
OUT_TABLEPLUS.close()
