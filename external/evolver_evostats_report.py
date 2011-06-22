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

FileName = sys.argv[1]

Events = [
  "Substitute", \
  "Delete", \
  "Invert", \
  "Move", \
  "Copy", \
  "Tandem", \
  "TandemExpand", \
  "TandemContract", \
  "Insert", \
  "MEInsert", \
  "InterChrCopy", \
  "InterChrMove", \
  "ChrSplit", \
  "ChrFuse", \
  "RecipTransloc", \
  "NonRecipTransloc", \
  "DeleteUTR", \
  "CreateNontermUTR", \
  "CreateTermUTR", \
  "CreateCDS", \
  "DeleteCDS", \
  "MoveStartCodonIntoUTR", \
  "MoveStartCodonIntoCDS", \
  "MoveStopCodonIntoUTR", \
  "MoveStopCodonIntoCDS", \
  "MoveUTRTerm", \
  "MoveCDSDonorIntoIntron", \
  "MoveUTRDonorIntoIntron", \
  "MoveDonorIntoCDS", \
  "MoveDonorIntoUTR", \
  "MoveCDSAcceptorIntoIntron", \
  "MoveUTRAcceptorIntoIntron", \
  "MoveAcceptorIntoCDS", \
  "MoveAcceptorIntoUTR", \
  "ChangeGeneSpeed", \
  "ChangeNGESpeed", \
  "DeleteNXE", \
  "CreateNXE", \
  "MoveNXE", \
  "DeleteNGE", \
  "CreateNGE", \
  "MoveNGE", \
  "DeleteIsland", \
  "CreateIsland", \
  "MoveIsland" \
  ]
 
ConstraintChangeEvents = [  \
  "DeleteUTR", \
  "CreateNontermUTR", \
  "CreateTermUTR", \
  "CreateCDS", \
  "DeleteCDS", \
  "MoveStartCodonIntoUTR", \
  "MoveStartCodonIntoCDS", \
  "MoveStopCodonIntoUTR", \
  "MoveStopCodonIntoCDS", \
  "MoveUTRTerm", \
  "MoveCDSDonorIntoIntron", \
  "MoveUTRDonorIntoIntron", \
  "MoveDonorIntoCDS", \
  "MoveDonorIntoUTR", \
  "MoveCDSAcceptorIntoIntron", \
  "MoveUTRAcceptorIntoIntron", \
  "MoveAcceptorIntoCDS", \
  "MoveAcceptorIntoUTR", \
  "ChangeGeneSpeed", \
  "ChangeNGESpeed", \
  "DeleteNXE", \
  "CreateNXE", \
  "MoveNXE", \
  "DeleteNGE", \
  "CreateNGE", \
  "MoveNGE", \
  "DeleteIsland", \
  "CreateIsland", \
  "MoveIsland" \
  ]

# fprintf(g_fStats, "EV_FAILS;%s;%u\n", EvType, FailCount);
# fprintf(g_fStats, "EV_REJECTS;%s;%u\n", EvType, RejectCount);
# fprintf(g_fStats, "EV_ACCEPTS;%s;%u\n", EvType, AcceptCount);
# fprintf(g_fStats, "EV_DETAIL;%s;%s;%s;%u\n", EvType, ResultType, Subtype, Count);
# fprintf(g_fStats, "EV_RANGE_FAIL;%s;%s;%u\n", RateType, BinStr, FailCount);
# fprintf(g_fStats, "EV_RANGE_REJECT;%s;%s;%u\n", RateType, BinStr, RejectCount);
# fprintf(g_fStats, "EV_RANGE_ACCEPT;%s;%s;%u\n", RateType, BinStr, AcceptCount);
# fprintf(g_fStats, "EV_ACC_BASES;%s;%u\n", RateType, Acc);
# fprintf(g_fStats, "EV_REJ_BASES;%s;%u\n", RateType, Rej);

def Die(s):
	print >> sys.stderr, "**ERROR**", s, sys.argv
	sys.exit(1)

Fails = {}
Rejects = {}
Accepts = {}
Details = {}
RangeFails = {}
RangeRejects = {}
RangeAccepts = {}
AccBases = {}
RejBases = {}

def Cmp(Ev1, Ev2):
	global Accepts
	N1 = 0
	N2 = 0
	if Ev1 in Accepts.keys():
		N1 = Accepts[Ev1]
	if Ev2 in Accepts.keys():
		N2 = Accepts[Ev2]

	if N1 > N2:
		return 1
	if N1 < N2:
		return -1
	return 0

def Cmp2(b1, b2):
	b1 = b1.replace("-", "")
	b1 = b1.replace(",", "")
	b1 = b1.replace(" ", "")
	b2 = b2.replace("-", "")
	b2 = b2.replace(",", "")
	b2 = b2.replace(" ", "")
	i1 = int(b1)
	i2 = int(b2)
	if i1 > i2:
		return 1
	if i1 < i2:
		return -1
	return 0

def Add(Dict, Key, N):
	if Key in Dict.keys():
		Dict[Key] += int(N)
	else:
		Dict[Key] = int(N)

File = open(sys.argv[1])
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
	Rec = Fields[0]
	Ev = Fields[1]
	if Ev not in Events:
		Events.append(Ev)

	if Rec == "EV_FAILS":
		Add(Fails, Ev, Fields[2])

	elif Rec == "EV_REJECTS":
		Add(Rejects, Ev, Fields[2])

	elif Rec == "EV_ACCEPTS":
		Add(Accepts, Ev, Fields[2])

	elif Rec == "EV_DETAIL":
		Key = Ev + ";" + Fields[2] + ";" + Fields[3]
		Add(Details, Key, Fields[4])

	elif Rec == "EV_RANGE_FAIL":
		Key = Ev + ";" + Fields[2]
		Add(RangeFails, Key, Fields[3])

	elif Rec == "EV_RANGE_REJECT":
		Key = Ev + ";" + Fields[2]
		Add(RangeRejects, Key, Fields[3])

	elif Rec == "EV_RANGE_ACCEPT":
		Key = Ev + ";" + Fields[2]
		Add(RangeAccepts, Key, Fields[3])

	elif Rec == "EV_ACC_BASES":
		Add(AccBases, Ev, Fields[2])

	elif Rec == "EV_REJ_BASES":
		Add(RejBases, Ev, Fields[2])

	else:
		Die("Unknown record type: " + Rec)

Total = 0
EvTotals = {}
for Ev in Accepts.keys():
	R = 0
	if Ev in Rejects.keys():
		R = Rejects[Ev]
	N = Accepts[Ev] + R
	EvTotals[Ev] = N
	Total += N
if Total == 0:
	Total = 1

SortedEvents = Accepts.keys()
SortedEvents.sort(Cmp)

print r"               Event       Fails     Rejects     Accepts        Acc%      Event%"
print r"--------------------  ----------  ----------  ----------  ----------  ----------"
for Ev in Events:
	s = "%20.20s" % Ev

	T = 0
	if Ev in Fails.keys():
		N = Fails[Ev]
		T += N
		s += "  %10u" % N
	else:
		s += "            "

	RejectCount = 0
	if Ev in Rejects.keys():
		RejectCount = Rejects[Ev]
		T += RejectCount
		s += "  %10u" % RejectCount
	else:
		s += "            "

	AcceptCount = 0
	if Ev in Accepts.keys():
		AcceptCount = Accepts[Ev]
		T += AcceptCount
		s += "  %10u" % AcceptCount
	else:
		s += "            "
	
	AcceptPct = 0
	N = AcceptCount + RejectCount
	if N > 0:
		AcceptPct = float(AcceptCount)*100.0/float(N)
	s += "  %9.1f%%" % AcceptPct

	EventPct = float(N)*100.0/float(Total)
	s += "  %9.3g%%" % EventPct

	if T > 0:
		print s

print ""
print r"               Event       Range       Fails     Rejects     Accepts       Total     Acc%"
print "--------------------  ----------  ----------  ----------  ----------  ----------  -------"

for Ev in Events:
	Any = 0	

	Bins = []
	BinFails = {}
	BinRejects = {}
	BinAccepts = {}

	for Key in RangeFails.keys():
		Fields = Key.split(";")
		if Fields[0] != Ev:
			continue
		Bin = Fields[1]
		if Bin not in Bins:
			Bins.append(Bin)
		BinFails[Bin] = RangeFails[Key]

	for Key in RangeRejects.keys():
		Fields = Key.split(";")
		if Fields[0] != Ev:
			continue
		Bin = Fields[1]
		if Bin not in Bins:
			Bins.append(Bin)
		BinRejects[Bin] = RangeRejects[Key]

	for Key in RangeAccepts.keys():
		Fields = Key.split(";")
		if Fields[0] != Ev:
			continue
		Bin = Fields[1]
		if Bin not in Bins:
			Bins.append(Bin)
		BinAccepts[Bin] = RangeAccepts[Key]

	Bins.sort(Cmp2)
	for Bin in Bins:
		Fails = 0
		Rejects = 0
		Accepts = 0
		
		if Bin in BinFails.keys():
			Fails = BinFails[Bin]
		if Bin in BinRejects.keys():
			Rejects = BinRejects[Bin]
		if Bin in BinAccepts.keys():
			Accepts = BinAccepts[Bin]

		Total = Fails + Rejects + Accepts
		if Total == 0:
			continue
		
		s = "%20.20s  %10.10s" % (Ev, Bin)
		
		if Fails > 0:
			s += "  %10u" % Fails
		else:
			s += "            "

		if Rejects > 0:
			s += "  %10u" % Rejects
		else:
			s += "            "

		if Accepts > 0:
			s += "  %10u" % Accepts
		else:
			s += "            "
		
		s += "  %10u" % Total
		
		if Total > 0:
			AccPct = (100.0*Accepts)/Total
			s += "  %6.1f%%" % AccPct
		
		print s
		Any = 1
	if Any:
		print ""

def GetSubtype(s):
	Fields = s.split(";")
	return Fields[2]
	
def CmpSTS(x, y):
	sx = GetSubtype(x)
	sy = GetSubtype(y)
	if sx > sy:
		return 1
	elif sx < sy:
		return -1
	else:
		return 0

print ""
print "               Event        Type           Subtype       Total      Pct"
print "--------------------  ----------  ----------------  ----------  -------"

for Ev in Events:
	if Ev == "sub_ins":
		print ""
	Any = 0
	Keys = Details.keys()
	Types = [ "Accept", "Reject", "Fail", "" ]
	for Type in Types:
		Total = 0
		SubtypeStrings = []
		for Key in Keys:
			Fields = Key.split(";")
			if Fields[0] != Ev:
				continue
			SubtypeStrings.append(Key)
			Count = Details[Key]
			Total += Count
		if Total == 0:
			Total = 1

		SubtypeStrings.sort(CmpSTS)
		for Key in SubtypeStrings:
			Fields = Key.split(";")
			Type2 = Fields[1]
			if Type2 != Type:
				# print "Type=%s Type2=%s" % (Type, Type2)
				continue

		for Key in SubtypeStrings:
			Fields = Key.split(";")
			Type2 = Fields[1]
			if Type2 != Type:
				# print "Type=%s Type2=%s" % (Type, Type2)
				continue
			Count = Details[Key]
			if Count == 0:
				continue
			Subtype = Fields[2]
			print "%20.20s  %10.10s  %16.16s  %10u  %6.1f%%" % (Ev, Type2, Subtype, Count, 100.0*float(Count)/Total)
			Any = 1
	if Any and not Ev in ConstraintChangeEvents:
		print ""

print ""
print "               Event   Acc.Bases   Rej.Bases"
print "--------------------  ----------  ----------"

Any = 0
for Ev in Events:
	Acc = 0
	Rej = 0
	if Ev in AccBases.keys():
		Acc = AccBases[Ev]
	if Ev in RejBases.keys():
		Rej = RejBases[Ev]
	if Acc == 0 and Rej == 0:
		continue
	Any = 1
	print "%20.20s  %10u  %10u" % (Ev, Acc, Rej)
if Any:
	print ""
