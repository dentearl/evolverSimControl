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
	
FileName1 = sys.argv[1]
FileName2 = sys.argv[2]
Name1 = FileName1
Name2 = FileName2
GenomeLength1 = -1
GenomeLength2 = -1
if len(sys.argv) > 3:
	Name1 = sys.argv[3]
if len(sys.argv) > 4:
	Name2 = sys.argv[4]
if len(sys.argv) > 6:
	GenomeLength1 = int(sys.argv[5])
	GenomeLength2 = int(sys.argv[6])
	

ConstrainedFeatures = [ "CDS", "UTR", "NXE", "NGE" ]

def Die(s):
	print >> sys.stderr, sys.argv[0], "***ERROR***", s
	sys.exit(1)

def DoRec():
	global Counts, Bases

	Length = gff.End - gff.Start + 1
	Feature = gff.Feature

	if Feature not in Bases.keys():
		Bases[Feature] = Length
		Counts[Feature] = 1
	else:
		Bases[Feature] += Length
		Counts[Feature] += 1
	
	if Feature in ConstrainedFeatures:
		Bases["Constrained"] += Length
		Counts["Constrained"] += 1

def Get(L, k):
	if k in L.keys():
		return L[k]
	return 0

def PctChg(x, y):
	if x == 0:
		if y == 0:
			return "100"
		else:
			return "--"
	else:
		return str(100*(y-x)/x)

def GetCounts(FileName):
	global Bases, Counts
	Bases = {}
	Counts = {}
	Bases["Constrained"] = 0
	Counts["Constrained"] = 0
	gff.GetRecs(FileName, DoRec)
	return Counts, Bases

Counts1, Bases1 = GetCounts(FileName1)
Counts2, Bases2 = GetCounts(FileName2)

Features = [ "CDS", "UTR", "NXE", "NGE", "island", "tandem", "Constrained" ]
Keys = Counts1.keys()
Keys.extend(Counts2.keys())
for Feature in Keys:
	if Feature not in Features:
		Features.append(Feature)

if GenomeLength1 != -1:
	Features.append("Neutral")
	Features.append("Total")
	Counts1["Neutral"] = 0
	Counts2["Neutral"] = 0
	Counts1["Total"] = 0
	Counts2["Total"] = 0
	Bases1["Neutral"] = GenomeLength1 -  Bases1["Constrained"]
	Bases2["Neutral"] = GenomeLength2 -  Bases2["Constrained"]
	Bases1["Total"] = GenomeLength1
	Bases2["Total"] = GenomeLength2

print "         Feature  1=%8.8s  2=%8.8s       Nr2-1   2-1 Pct      Bases1      Bases2    Bases2-1   2-1 Pct" % (Name1, Name2)
print "================  ==========  ==========  ==========  ========  ==========  ==========  ==========  ========"
for Feature in Features:
	n1 = Get(Counts1, Feature)
	n2 = Get(Counts2, Feature)
	dn = n2 - n1
	b1 = Get(Bases1, Feature)
	b2 = Get(Bases2, Feature)
	db = b2 - b1
	pn = PctChg(n1, n2)
	pb = PctChg(b1, b2)
	s = ""
	s += "%16.16s" % Feature
	s += "  %10u" % n1
	s += "  %10u" % n2
	s += "  %+10d" % (n2 - n1)
	s += "  %7.7s%%" % pn
	s += "  %10u" % b1
	s += "  %10u" % b2
	s += "  %+10d" % (b2-b1)
	s += "  %7.7s%%" % pb
	print s
