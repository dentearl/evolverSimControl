#!/usr/bin/python

import sys
import evolverSimControl.lib.evolver_gff as gff

FileName = sys.argv[1]

def DoRec():
	gff.WriteRec(sys.stdout)

gff.GetSortedRecs(FileName, DoRec)
