#!/usr/bin/python
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
# Re-assign gene_index
# Designed for interchr evolver where each chr has gene_index's
# starting from 0 but ice needs unique.

import sys
import evolver_gff # -dae, modification for simCtrl
import re

def Die(s):
	print >> sys.stderr, sys.argv[0], "***ERROR***", s
	sys.exit(1)

def DoRec():
	global GeneIndexes

	AttrDict = evolver_gff.GetAttrDict()
	CurrGeneIndex = evolver_gff.GetIntAttr("gene_index", -1)
	if CurrGeneIndex == -1:
		evolver_gff.WriteRec(sys.stdout)
		return
	
	Key = evolver_gff.Label + "." + str(CurrGeneIndex)
	if Key in GeneIndexes.keys():
		NewGeneIndex = GeneIndexes[Key]
	else:
		NewGeneIndex = len(GeneIndexes)
		GeneIndexes[Key] = NewGeneIndex
	
	AttrDict["gene_index"] = NewGeneIndex
	evolver_gff.SetAttrsFromDict(AttrDict)

	evolver_gff.WriteRec(sys.stdout)	

FileName = sys.argv[1]
GeneIndexes = {}

evolver_gff.GetRecs(FileName, DoRec)

print >> sys.stderr, len(GeneIndexes), "genes found"
