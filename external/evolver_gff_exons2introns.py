#!/usr/bin/env python
# Copyright (C) 2008-2011 by
# George Asimenos, Robert C. Edgar, Serafim Batzoglou and Arend Sidow.
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# 
##############################
import sys
import evolverSimControl.lib.evolver_gff as gff

#########################################################################
# Set SORT=0 if input is known to be sorted, this can be significantly
# faster for large files.
# WARNING: this script does NOT check that input is sorted; if you set
# SORT=0 and give it incorrectly sorted input the output will be GARBAGE.
#########################################################################
SORT=1

def Die(s):
	print >> sys.stderr, sys.argv[0], "***ERROR***", s
	sys.exit(1)
	
FileName = sys.argv[1]

def DoRec():
	global LastGeneIndex
	global LastExonStart
	global LastLabel
	global LastExonEnd
	global IntronCounts
	global ExonCounts

	if gff.Feature != "exon":
		return

	GeneIndex = gff.GetRequiredIntAttr("gene_index")
	Key = gff.Label + "%%%" + str(GeneIndex)
	if Key not in ExonCounts.keys():
		ExonCounts[Key] = 1
	else:
		ExonCounts[Key] += 1
	# print "ExonCounts[%s] = %u" % (Key, ExonCounts[Key])
	
	Start = gff.Start
	End = gff.End
	Label = gff.Label
	
	if GeneIndex == LastGeneIndex and gff.Label == LastLabel and LastExonStart != -1:
		if LastExonStart != -1:
			gff.Label = LastLabel
			gff.Start = LastExonEnd + 1
			gff.End = Start - 1
			gff.Source = "exons2introns"
			gff.Feature = "intron"
			gff.Frame = "."
			gff.Strand = "."
			gff.Attrs = "gene_index %u; exons %u-%u,%u-%u;" % (GeneIndex, LastExonStart, LastExonEnd, Start, End)
			gff.WriteRec(sys.stdout)
			if Key not in IntronCounts.keys():
				IntronCounts[Key] = 1
			else:
				IntronCounts[Key] += 1
			# print >> sys.stderr, "IntronCounts[%s] = %u" % (Key, IntronCounts[Key])
		
	LastGeneIndex = GeneIndex
	LastExonStart = Start
	LastExonEnd = End
	LastLabel = Label

IntronStart = -1
LastGeneIndex = -1
LastExonStart = -1
LastExonEnd = -1
LastLabel = ""

ExonCounts = {}
IntronCounts = {}

if SORT:
	gff.GetSortedRecs(FileName, DoRec)
else:
	gff.GetRecs(FileName, DoRec)

print >> sys.stderr, "           Label        Gene       Exons     Introns"
print >> sys.stderr, "================  ==========  ==========  =========="
for Key in ExonCounts.keys():
	Fields = Key.split("%%%")

	Label = Fields[0]
	GeneIndex = int(Fields[1])

	ExonCount = ExonCounts[Key]
	IntronCount = 0
	if ExonCount == 0:
		Die("%s gene index %u has no exons" % (Label, GeneIndex, ExonCount))
	if ExonCount == 1 and Key in IntronCounts.keys():
		Die("%s gene index %u has %u exons and %u introns" % (Label, GeneIndex, ExonCounts[Key], IntronCounts[Key]))
	if ExonCount > 1:
		if Key not in IntronCounts.keys():
			Die("%s gene index %u has %u exons but no introns" % (Label, GeneIndex, ExonCounts[Key]))
		IntronCount = IntronCounts[Key]
		if IntronCount != ExonCount - 1:
			Die("Gene index %u has %u exons and %u introns" % (GeneIndex, ExonCounts[Key], IntronCounts[Key]))
	print >> sys.stderr, "%16.16s  %10u  %10u  %10u" % (Label, GeneIndex, ExonCount, IntronCount)
