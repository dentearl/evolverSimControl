#!/usr/bin/python
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

FileName = sys.argv[1]

def Die(s):
	print >> sys.stderr, "**ERROR**", s, sys.argv
	sys.exit(1)

def DoRec():
	global GeneLos, GeneHis, GeneStrands, Label, Lo, Hi
	
	if gff.Feature != "CDS" and gff.Feature != "UTR" and gff.Feature != "exon":
		return
		
	if Label == "":
		Label = gff.Label
	elif Label != gff.Label:
		Die("More than one label in EVOLVER_GFF")
	
	if Hi == -1:
		Hi = gff.End
	else:
		if gff.End > Hi:
			Hi = gff.End
	
	GeneIndex = gff.GetRequiredIntAttr("gene_index")
	if GeneIndex in GeneLos.keys():
		if gff.Start < GeneLos[GeneIndex]:
			GeneLos[GeneIndex] = gff.Start
		if gff.End > GeneHis[GeneIndex]:
			GeneHis[GeneIndex] = gff.End
		if GeneStrands[GeneIndex] != gff.Strand:
			Die("Gene on both strands")
	else:
		GeneLos[GeneIndex] = gff.Start
		GeneHis[GeneIndex] = gff.End
		GeneStrands[GeneIndex] = gff.Strand

Hi = -1
GeneLos = {}
GeneHis = {}
GeneStrands = {}

Label = ""
gff.GetRecs(FileName, DoRec)

gff.Source = "gene_lengths"
gff.Feature = "gene"
gff.Score = 0
gff.Frame = "."

TotGeneLength = 0
for GeneIndex in GeneLos.keys():
	gff.Start = GeneLos[GeneIndex]
	gff.End = GeneHis[GeneIndex]
	gff.Strand = GeneStrands[GeneIndex]
	gff.Attrs = "gene_index %u;" % GeneIndex
	gff.WriteRec(sys.stdout)
	TotGeneLength += gff.End - gff.Start + 1

print >> sys.stderr, "Max annot end     %10u" % Hi
print >> sys.stderr, "Total gene length %10u" % TotGeneLength
print >> sys.stderr, "Pct               %10.1f%%" % (float(TotGeneLength)*100/Hi)
