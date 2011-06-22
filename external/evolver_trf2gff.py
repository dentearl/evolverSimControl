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

MAX_MOTIF_ATTR = 32
FEATURE = "tandem"

FileName = sys.argv[1]

def Quit(s):
	print >> sys.stderr, "*** ERROR ***", sys.argv[0],s
	sys.exit(0)

File = open(FileName)

# Skip 6-line header
for i in range(0, 6):
	Line = File.readline()

while 1:
	Line = File.readline()
	if len(Line) == 0:
		break
	Line = Line.strip()
	if len(Line) == 0:
		continue
	elif Line.startswith("Sequence: "):
		Label = Line[10:]
		continue
	elif Line.startswith("Parameters: "):
		continue

# 1 59 7 8.0 7 83 9 64 37 38 1 22 1.64 CCCTAAA CCCTAAACCCTAAACCCTAAACCCTAAACCTCTGAATCCTTAATCCCTAAATCCCTAAA
# 0  1 2   3 4  5 6  7  8  9 19 11  12 13      14

# 0, 1	Indices of the repeat relative to the start of the sequence. 
# 2		Period size of the repeat. 
# 3		Number of copies aligned with the consensus pattern. 
# 4		Size of consensus pattern (may differ slightly from the period size). 
# 5		Percent of matches between adjacent copies overall. 
# 6		Percent of indels between adjacent copies overall. 
# 7		Alignment score. 
# 8-11	Percent composition for each of the four nucleotides. 
# 12	Entropy measure based on percent composition. 
# 13	Consensus motif
# 14	Sequence

	Fields = Line.split()
	if len(Fields) != 15:
		Quit("Expected 15 fields, got: " + Line)
	
	Start = int(Fields[0])
	End = int(Fields[1])
	Copies = float(Fields[3])
	Score = int(Fields[7])
	Motif = Fields[13]
	Length = len(Motif)
	Attrs = "replen %u;" % Length
	Attrs += " copies %.1f;" % Copies
	if len(Motif) > MAX_MOTIF_ATTR:
		Attrs += " cons \"%s...\";" % Motif[0:MAX_MOTIF_ATTR]
	else:
		Attrs += " cons \"%s\";" % Motif

	s = "%s\ttrf\t%s\t%d\t%d\t%d\t+\t.\t%s" % (Label, FEATURE, Start, End, Score, Attrs)
	print s
