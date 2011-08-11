#!/usr/bin/python

import sys
import evolver_gff

FileName = sys.argv[1]

def DoRec():
	evolver_gff.WriteRec(sys.stdout)

evolver_gff.GetSortedRecs(FileName, DoRec)
