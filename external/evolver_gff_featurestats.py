#!/usr/bin/python
# Re-assign gene_index
# Designed for interchr evolver where each chr has gene_index's
# starting from 0 but ice needs unique.

import sys
import gff
import re

def Die(s):
	print >> sys.stderr, sys.argv[0], "***ERROR***", s
	sys.exit(1)

def DoRec():
	global RecordCounts, BaseCounts
	
	Length = gff.End - gff.Start + 1

	if gff.Feature not in RecordCounts.keys():
		RecordCounts[gff.Feature] = {}
		BaseCounts[gff.Feature] = {}
	
	if gff.Label not in RecordCounts[gff.Feature].keys():
#		print >> sys.stderr, "SET", gff.Label, gff.Feature, "keys=", RecordCounts[gff.Feature].keys()
		RecordCounts[gff.Feature][gff.Label] = 1
		BaseCounts[gff.Feature][gff.Label] = Length
	else:
		RecordCounts[gff.Feature][gff.Label] += 1
		BaseCounts[gff.Feature][gff.Label] += Length
	
FileName = sys.argv[1]

RecordCounts = {}
BaseCounts = {}

gff.GetRecs(FileName, DoRec)

print "             Seq     Feature        Recs       Bases"
print "================  ==========  ==========  =========="
for Feature in RecordCounts.keys():
	for Label in RecordCounts[Feature].keys():
		RecordCount = RecordCounts[Feature][Label]
		BaseCount = BaseCounts[Feature][Label]
		print "%16.16s  %10.10s  %10u  %10u" % (Label, Feature, RecordCount, BaseCount)
	print ""
