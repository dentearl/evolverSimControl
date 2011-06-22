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

def Die(s):
	print >> sys.stderr, "**ERROR**", s, sys.argv
	sys.exit(1) 

Values = {}
for i in range(1, len(sys.argv)):
	FileName = sys.argv[i]
	File = open(FileName)
	while 1:
		Line = File.readline()
		if len(Line) == 0:
			break
		Line = Line.strip()
		if len(Line) == 0:
			continue
		if Line[0] == "#":
			continue
		Fields = Line.split(";")
		N = len(Fields)
		if N == 1:
			Die("File '%s', no semi-colons in line: %s"  % (FileName, Line))
		Value = int(Fields[N-1])
		Key = ""
		for i in range(0, N-1):
			Key += Fields[i] + ";"
		if Key in Values.keys():
			Values[Key] += Value
		else:
			Values[Key] = Value

for Key in Values.keys():
	print "%s%u" % (Key, Values[Key])
