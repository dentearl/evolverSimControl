#!/usr/bin/python

import sys
import gff

#########################################################################
# Set SORT=0 if input is known to be sorted, this can be significantly
# faster for large files.
# WARNING: this script does NOT check that input is sorted; if you set
# SORT=0 and give it incorrectly sorted input the output will be GARBAGE.
#########################################################################
SORT=1

FileName = sys.argv[1]

def Die(s):
	print >> sys.stderr, sys.argv[0], "***ERROR***", s
	sys.exit(1)
	
def IsExonFeature(Feature):
	return Feature == "CDS" or Feature == "UTR"

def DoRec():
	global LastGeneIndex
	global ExonStart
	global LastLabel
	global LastEnd
	global Attr

	if not IsExonFeature(gff.Feature):
		return
	
	GeneIndex = gff.GetRequiredIntAttr("gene_index")
	Start = gff.Start
	End = gff.End
	Feature = gff.Feature
	Label = gff.Label

	if GeneIndex != LastGeneIndex or gff.Label != LastLabel or gff.Start != LastEnd + 1:
		if ExonStart != -1:
			gff.Label = LastLabel
			gff.Start = ExonStart
			gff.End = LastEnd
			gff.Source = "cdsutr2exons"
			gff.Feature = "exon"
			gff.Frame = "."
			gff.Strand = "."
			gff.Attrs = "gene_index %u; ces %s;" % (LastGeneIndex, Attr)
			gff.WriteRec(sys.stdout)
			Attr = ""
		
		ExonStart = Start

	LastGeneIndex = GeneIndex
	LastStart = Start
	LastLabel = Label
	LastEnd = End
	s = "%s:%u-%u" % (Feature, Start, End)
	if Attr == "":
		Attr = s
	else:
		Attr += "," + s

ExonStart = -1
LastGeneIndex = -1
LastStart = -1
LastEnd = -1
LastFeature = ""
LastLabel = ""
Attr = ""

if SORT:
	gff.GetSortedRecs(FileName, DoRec)
else:
	gff.GetSortedRecs(FileName, DoRec)

if ExonStart != -1:
	gff.Start = ExonStart
	gff.End = LastEnd
	gff.Source = "cdsutr2exons"
	gff.Feature = "exon"
	gff.Frame = "."
	gff.Strand = "."
	gff.Attrs = "gene_index %u; ces %s;" % (LastGeneIndex, Attr)
	gff.WriteRec(sys.stdout)
