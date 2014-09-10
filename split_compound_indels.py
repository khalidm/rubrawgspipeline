#!/usr/bin/env python2.7
# encoding: utf-8
"""
split_compound_indels.py

Created by Matthew Wakefield on 2013-10-17.
Copyright (c) 2013  Matthew Wakefield and The University of Melbourne. All rights reserved.
"""
from __future__ import print_function, unicode_literals, division
import sys, os
import re
from io import StringIO
import unittest
import argparse

__author__ = "Matthew Wakefield"
__copyright__ = "Copyright 2013,  Matthew Wakefield and The University of Melbourne"
__credits__ = ["Matthew Wakefield",]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Matthew Wakefield"
__email__ = "wakefield@wehi.edu.au"
__status__ = "Development"


def indel_lines(mpileupfile):
    inputlines = 0
    outputlines = 0
    for line in mpileupfile:
        inputlines +=1
        if int(line.split()[3]) != len(line.split()[4])-(line.count('$')+(line.count('^')*2)):
            outputlines += 1
            yield line.strip('\n')
    print('Processed {0} indel containing positions from a total of {1} positions'.format(outputlines,inputlines), file=sys.stderr)

def split_out_indels(pileupstring, errorstring):
    snpstring = ''
    snperrors = ''
    indelstring = ''
    indelerrors = ''
    pileupstring_position = 0
    errorstring_position = 0
    while pileupstring_position < len(pileupstring):
        if pileupstring[pileupstring_position] in '.,*':
            snpstring += pileupstring[pileupstring_position]
            snperrors += errorstring[errorstring_position]
            indelstring += pileupstring[pileupstring_position]
            indelerrors += errorstring[errorstring_position]
            pileupstring_position += 1
            errorstring_position += 1
        elif pileupstring[pileupstring_position] == '$':
            snpstring += pileupstring[pileupstring_position]
            indelstring += pileupstring[pileupstring_position]
            pileupstring_position += 1
        elif pileupstring[pileupstring_position] == '^':
            snpstring += pileupstring[pileupstring_position:pileupstring_position+2]
            indelstring += pileupstring[pileupstring_position:pileupstring_position+2]
            pileupstring_position += 2
        elif pileupstring[pileupstring_position] == '-':
            feature_size = int(re.match(r"\d+",pileupstring[pileupstring_position+1:]).group(0))
            description_length = 1+len(str(feature_size))+feature_size # eg -3CGG = 1+1+3 = 5
            indelstring += pileupstring[pileupstring_position:pileupstring_position+description_length]
            pileupstring_position += description_length
        elif pileupstring[pileupstring_position] == '+':
            feature_size = int(re.match(r"\d+",pileupstring[pileupstring_position+1:]).group(0))
            description_length = 1+len(str(feature_size))+feature_size # eg -3CGG = 1+1+3 = 5
            indelstring += pileupstring[pileupstring_position:pileupstring_position+description_length]
            indelerrors += errorstring[errorstring_position:errorstring_position+feature_size]
            pileupstring_position += description_length
        else:
            snpstring += pileupstring[pileupstring_position]
            snperrors += errorstring[errorstring_position]
            indelstring += '.'
            indelerrors += errorstring[errorstring_position]
            pileupstring_position += 1
            errorstring_position += 1
    return (snpstring, snperrors), (indelstring, indelerrors)

def command_line_interface(*args,**kw):
    parser = argparse.ArgumentParser(description='A script for to split positions with compound \
                                                indels and SNVs into separate entries.\
                                                Takes a single sample mpileup file and produces \
                                                a mpileup file with two entries for each indel \
                                                containing position - one with only SNVs and one \
                                                with only indels.  This can be used with single \
                                                position indel callers to produce either SNV or \
                                                indel calls.\
                                                This script may have unintended consequences on \
                                                calling and should ONLY be used as a secondary \
                                                screen to rescue missed compound mutations.')
    parser.add_argument('mpileup',
                        type=argparse.FileType('U'),
                        default=sys.stdin,
                        help='a single sample mpileup or pileup file')
    parser.add_argument('--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='a mpileup file containing a SNP and indel entry for each indel position')
    return parser.parse_args(*args,**kw)

def main():
    args = command_line_interface()
    for line in indel_lines(args.mpileup):
        chromosome, position, reference, readcount, bases, qualities = line.split('\t')
        (snpstring, snperrors), (indelstring, indelerrors) = split_out_indels(bases, qualities)
        print(chromosome, position, reference, readcount, snpstring, snperrors, sep='\t', file=args.output)
        print(chromosome, position, reference, readcount, indelstring, indelerrors, sep='\t', file=args.output)

class test_split_compound_indels(unittest.TestCase):
    def setUp(self):
        self.mpileupfile = StringIO(''.join(["1\t27023256\tC\t35\t...........T-3CGG.................,,.,,,,\tFGG@GFEGGGGGGGGGGGGGGGGEGGGGFCGDFFF\n",
                                    "1\t27023257\tC\t35\t....$......*.................,,.,,,,\tGGFEGGGGGGGGGGGGGGGGGGGEGGGGFAGAA>C\n",
                                    "1\t27023258\tG\t34\t.........*.................,,.,,,,\tGG?FGGEG?GEGGGGGGGGGGGEGEGGFFCFDDD\n",
                                    "1\t27023259\tG\t34\t..^I.......*.................,,.,,,,\t?G>GGGGGEGFGFGGGGGGGCG?EDEGFAGB=@9\n",
                                    ])
                               )

        pass
    
    def test_indel_lines(self):
        self.assertEqual(list(indel_lines(self.mpileupfile)),['1\t27023256\tC\t35\t...........T-3CGG.................,,.,,,,\tFGG@GFEGGGGGGGGGGGGGGGGEGGGGFCGDFFF\n'])
        pass
    
    def test_split_out_indels(self):
        pileupstring = "...........T-3CGG.................,,.,,,,"
        errorstring =  "ABCDEFGHIJKL     MNOPQRSTUVWXYZ!@#$%^&*()".replace(' ','') #spaces only for clarity
        snpstring =    "...........T.................,,.,,,,"
        snperrors =    "ABCDEFGHIJKLMNOPQRSTUVWXYZ!@#$%^&*()".replace(' ','') #spaces only for clarity
        indelstring =  "............-3CGG.................,,.,,,,"
        indelerrors =  "ABCDEFGHIJKL     MNOPQRSTUVWXYZ!@#$%^&*()".replace(' ','') #spaces only for clarity
        self.assertEqual(split_out_indels(pileupstring, errorstring), ((snpstring, snperrors), (indelstring, indelerrors)))


if __name__ == '__main__':
    #unittest.main()
    main()

