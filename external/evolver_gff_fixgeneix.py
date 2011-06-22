#!/usr/bin/python
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
