# Created on 02.03.2017
# Edited on 25.04.2017

#To run: python ~/Desktop/Pysam_AlignmentStatistics_TestLoop.py ~/Documents/Singapore_StemCells/Script_Output/BWA_OutputFiles/Sumba/YamagishiReads_Human/aligned_human_${file}.sam

import sys
import pysam
from collections import Counter
import csv

if len(sys.argv) < 2:
    sys.exit("usage: alnstat.py in.sam")

fname = sys.argv[1]

#using Pysam to view samfile
samfile = pysam.AlignmentFile(fname)

stats = Counter()

for read in samfile:
    stats['total'] += 1
    stats['qcfail'] += int(read.is_qcfail)
    if read.is_unmapped:
        stats['unmapped'] += 1
        continue # other flags don't apply
    # record if mapping quality <= 30
    stats["mapping quality <= 30"] += int(read.mapping_quality <= 30)
    stats['mapped'] += 1

# specify the output order, since dicts don't have order

output_order = ("total",
                "mapped",
                "unmapped",
                "qcfail",
                "mapping quality <= 30")

# format output and print to standard out
for key in output_order:
    format_args = (stats[key])
    sys.stdout.write("{}\n".format(format_args))

