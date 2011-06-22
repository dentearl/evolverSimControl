import re
import sys

MaxError = -1
TargetLabel = ""
TargetStart = -1
TargetEnd = -1
Start = -1
End = -1
Score = -1
Strand = ""
Frame = "."
Attrs = ""

def Quit(s):
	print >> sys.stderr, "*** ERROR ***", s, sys.argv
	sys.exit(1)

def ProgressFile(File, FileSize):
#	if not sys.stderr.isatty():
	return
	Pos = File.tell()
	Pct = (100.0*Pos)/FileSize
	Str = "%5.1f%%\r" % Pct
	sys.stderr.write(Str)

def Progress(i, N):
#	if not sys.stderr.isatty():
	return
	Pct = (100.0*i)/N
	Str = "%5.1f%%\r" % Pct
	sys.stderr.write(Str)

def AppendAttr(a):
	global Attrs

	if len(Attrs) > 0:
		Attrs += " ; "
	Attrs += a

def AppAttr(Name, Value):
	global Attrs
	if len(Attrs) > 0:
		Attrs += " "
	Attrs += "%s %s;" % (Name, Value)

def RepeatFullLength():
	global RepeatStart
	global RepeatEnd
	global RepeatLeft

	return RepeatEnd + RepeatLeft

def RepeatHitLength():
	global RepeatStart
	global RepeatEnd
	global RepeatLeft

	return RepeatEnd - RepeatStart + 1

def RepeatFract():
	return float(RepeatHitLength())/float(RepeatFullLength())

def HasTarget():
	return Attrs.find("Target ") >= 0

def HasComboCCType():
	return Attrs.find("ComboCCType[") >= 0

def HasRepeat():
	return Attrs.find("Repeat ") >= 0

def HasBundle():
	return Attrs.find("Bundle[") >= 0

def HitLength():
	return End - Start + 1

def HasSP():
	return Attrs.find("SP") >= 0

def HasDeleted():
	return Attrs.find("Deleted") >= 0
	
def HasDeletedMerged():
	return Attrs.find("DeletedMerged") >= 0

def HasCRCC():
	return Attrs.find("CRCC") >= 0
	
def GetAttrDict():
	AttrDict = {}
	a = Attrs
	while 1:
		Match = re.match("[ \t;]*([A-Za-z_][A-Za-z_0-9]*)[ \t]+([^;]*);", a)
		if Match == None:
			break
		AttrName = Match.groups()[0]
		AttrValue = Match.groups()[1]
		if AttrValue[0] == '"':
			AttrValue = AttrValue[1:-1]
		AttrDict[AttrName] = AttrValue
		Match = re.match("([ \t;]*[A-Za-z_][A-Za-z_0-9]*[ \t]+([^;]*);)", a)
		if Match == None:
			print "ERROR"
			break			
		x = Match.groups()[0]
		n = len(x)
		a = a[n:]
		#print ""
		#print "x=%s" % x
		#print "n=%d" % n
		#print "a=%s" % a
		#print ""

	return AttrDict

def GetAttrDictFromStr(a):
	AttrDict = {}
	while 1:
		Match = re.match("[ \t;]*([A-Za-z_][A-Za-z_0-9]*)[ \t]+([^;]*);", a)
		if Match == None:
			break
		AttrName = Match.groups()[0]
		AttrValue = Match.groups()[1]
		if AttrValue[0] == '"':
			AttrValue = AttrValue[1:-1]
		AttrDict[AttrName] = AttrValue
		Match = re.match("([ \t;]*[A-Za-z_][A-Za-z_0-9]*[ \t]+([^;]*);)", a)
		if Match == None:
			print "ERROR"
			break			
		x = Match.groups()[0]
		n = len(x)
		a = a[n:]
		#print ""
		#print "x=%s" % x
		#print "n=%d" % n
		#print "a=%s" % a
		#print ""

	return AttrDict

def GetStrAttr(Name, DefaultValue):
	AD = GetAttrDict()
	if Name in AD.keys():
		return AD[Name]
	else:
		return DefaultValue

def GetIntAttr(Name, DefaultValue):
	AD = GetAttrDict()
	if Name in AD.keys():
		return int(AD[Name])
	else:
		return DefaultValue

def GetRequiredStrAttr(Name):
	Value = GetStrAttr(Name, None)
	if Value == None:
		Quit("Required attr '%s' not found" % Name)
	return Value

def GetRequiredIntAttr(Name):
	Value = GetStrAttr(Name, None)
	if Value == None:
		Quit("Required attr '%s' not found" % Name)
	return int(Value)

def SetAttrsFromDict(AttrDict):
	global Attrs
	Attrs = ''
	for Name in AttrDict.keys():
		Value = AttrDict[Name]
		if len(Attrs) > 0:
			Attrs += ' '
		if type(Value) == int:
			Attrs += Name + ' ' + str(Value) + ';'
		else:
			Attrs += Name + ' ' + '"' + str(Value) + '";'

def PokeRepeatAttrs():
	global Attrs
	global RepeatName
	global RepeatFamily
	global RepeatStart
	global RepeatEnd
	global RepeatLeft

	if RepeatStart == -1:
		Attrs = "Repeat %s %s . . ." % (RepeatName, RepeatFamily)
	else:
		Attrs = "Repeat %s %s %d %d %d" % (RepeatName, RepeatFamily, RepeatStart, RepeatEnd, RepeatLeft)

def PokeTargetAttrs():
	global Attrs
	global TargetLabel
	global TargetStart
	global TargetEnd
	global MaxError

	if MaxError == -1:
		NewAttrs = "Target %s %d %d" % (TargetLabel, TargetStart, TargetEnd)
	else:
		NewAttrs = "Target %s %d %d ; maxe %.3g" % (TargetLabel, TargetStart, TargetEnd, MaxError)

	Attrs = re.sub("Target[^;]*;", "", Attrs)
	if len(Attrs) == 0:
		Attrs = NewAttrs
	else:
		Attrs = NewAttrs + " ; " + Attrs 	

def GetBundleIndex():
	Match = re.search("Bundle\[(\d+)\]", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (Bundle)"
		return 0
	str = Match.groups()[0]
	return int(str)	

def GetNumericAttr(Name):
	Match = re.search(Name + "\[([0-9.]+)\]", Attrs)
	if Match == None:
		print >> sys.stderr, "No match GetNumericAttr", Name
		return 0
	str = Match.groups()[0]
	return float(str)
	
def GetCCIndex():
	Match = re.search("CC\[(\d+)\]", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (CC)"
		return 0
	str = Match.groups()[0]
	return int(str)	
	
def GetGeneIndex():
	a = GetAttr("gene_index")
	if a == None:
		print >> sys.stderr, "gene_index not found in " + Attrs
		return -1
	return int(a)	

def HasAttr(a):
	return GetAttr(a) != None
	
def GetGeneId():
	Match = re.search('gene_id "([^"]+)"', Attrs)
	if Match == None:
		print >> sys.stderr, "No match (GeneId)"
		return 0
	return Match.groups()[0]
	
def GetSPIndex():
	Match = re.search("SP *(\d+)", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (SP)"
		return 0
	str = Match.groups()[0]
	return int(str)	

# CRCC[46]
def GetCRCC():
	Match = re.search("CRCC\[(\d+)]", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (CRCC): " + Attrs
		return 0

	strIndex = Match.groups()[0]
	return int(strIndex)
	
# ComboCCType[Clique]
def GetComboCCType():
	Match = re.search("ComboCCType\[([^\]]+)]", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (ComboCCType): " + Attrs
		return 0

	Type = Match.groups()[0]
	return Type

# Pyramid 46 
def ParsePyramidAttrs():
	global PyramidIndex

	Match = re.search("Pyramid (\d+)", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (Pyramid): " + Attrs
		return 0

	strIndex = Match.groups()[0]
	PyramidIndex = int(strIndex)

# TESFam 46 
def ParseTESAttrs():
	global TESFam
	global TESRev
	global TESCliffIndex

	Match = re.search("TESFam (\d+) (\d+) ([-+])", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (TESFam): " + Attrs
		return 0

	strFam = Match.groups()[0]
	strCliffIndex = Match.groups()[1]
	strRev = Match.groups()[2]
	
	TESFam = int(strFam)
	TESCliffIndex = int(strCliffIndex)
	TESRev = (strRev == "-")

def ParseBlastLengthAttrs():
	global BlastTargetLength
	global BlastQueryLength
	
	Match = re.search("QueryLength (\d+)", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (QueryLength): " + Attrs
		return 0
	BlastQueryLength = int(Match.groups()[0])

	Match = re.search("TargetLength (\d+)", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (TargetLength): " + Attrs
		return 0
	BlastTargetLength = int(Match.groups()[0])

# Local seq1 1235995 1236048 seq2 62 114
def ParseLocalAttrs():
	global LocalLabel
	global LocalStart
	global LocalEnd
	global LocalTargetLabel
	global LocalTargetStart
	global LocalTargetEnd

	Match = re.search("Local ([^ ]+) (\d+)\ (\d+) ([^ ]+) (\d+)\ (\d+)", Attrs)
	if Match == None:
		return 0
	
	Fields = Match.groups()
	LocalLabel = Fields[0]
	LocalStart = int(Fields[1])
	LocalEnd = int(Fields[2])
	LocalTargetLabel = Fields[3]
	LocalTargetStart = int(Fields[4])
	LocalTargetEnd = int(Fields[5])
	
	return 1

# Family 46.35 ; Pile 2
def ParseTRSAttrs():
	global TRSFam
	global TRSSuperFam
	global TRSPile

	Match = re.search("Family (\d+)\.(\d+)", Attrs)
	if Match == None:
		Match = re.search("Family (\d+)", Attrs)
		if Match == None:
			print "No match (Family): " + Attrs
			return 0
		else:
			strFam = Match.groups()[0]
			strSuperFam = "-1"
	else:
		strFam = Match.groups()[0]
		strSuperFam = Match.groups()[1]
	
	TRSFam = int(strFam)
	TRSSuperFam = int(strSuperFam)
	
	Match = re.search("Pile (\d+)", Attrs)
	if Match == None:
		print "No match (Pile): " + Attrs
		return 0

	strPile = Match.groups()[0]
	TRSPile = int(strPile)

def HasAnnot():
	Match = re.search('Annot "([^"]*)"', Attrs)
	return Match != None

# Annot "11111111111111111111  AT9NMU1(100%)"
def ParseAnnotAttrs():
	global AnnotString

	Match = re.search('Annot "([^"]*)"', Attrs)
	if Match == None:
		print "No match (Annot): " + Attrs
		return 0

	AnnotString = Match.groups()[0]

def GetAttr(Name):
	global Attrs

# '@' is used to make sure Name is not a prefix of some other name.
# Integer attr
	Match = re.search("[^A-Za-z]" + Name + r"[ \t]+(\d+)", "@" + Attrs)
	if Match != None:
		return Match.groups()[0]
	
# String attr
# Replace \" by \v to deal with embedded "s.
	AttrsCopy = Attrs.replace(r'\"', "\v")
	Match = re.search("[^A-Za-z]" + Name + r'[ \t]+"([^"]*)"', "@" + AttrsCopy)
	if Match != None:
		Value = Match.groups()[0]
		return Value.replace("\v", '"')
	
	return None

def GetRequiredAttr(Name):
	Value = GetAttr(Name)
	if Value == None:
		Quit("Required attribute " + Name + " not found in GFF record")
	return Value

def GetAttrs():
	global Attrs
	
	Values = {}

	AttrsCopy = Attrs	
	while 1:
		SpaceNameSpacePatt = r"[; \t]*([A-Za-z][A-Za-z0-9_]*)[ \t]*"
		Match = re.search(SpaceNamePatt, AttrsCopy)
		if Match == None:
			return Values

		Name = Match.groups()[0]		
		AttrsCopy = re.sub(SpaceNamePatt, "", AttrsCopy)

		IntPatt = r"(\d)+"
		Match = re.search(IntPatt, AttrsCopy)
		if Match == None:
			StrPatt = r'"([^"]*"'
			Match = re.search(StrPatt, AttrsCopy)
			if Match == None:
				return Values
			else:
				AttrsCopy = re.sub(StrPatt, "", AttrsCopy)
		else:
			AttrsCopy = re.sub(IntPatt, "", AttrsCopy)			
		Value = Match.groups()[0]
		Values[Name] = Value

# Target <label> <begin> <end>
def ParseTargetAttrs():
	global TargetLabel
	global TargetStart
	global TargetEnd
	global MaxError

	Match = re.search("Target ([^;]*)", Attrs)
	if Match == None:
		print >> sys.stderr, "No match (Target) : " + Attrs
		return 0

	a = Match.groups()[0]
	
	Match2 = re.match("([^ ]+) (\d+) (\d+)", a)
	if Match2 == None:
		Quit("Invalid Target attribute in GFF file: " + a)

	g = Match2.groups()
	TargetLabel = g[0]
	TargetStart = int(g[1])
	TargetEnd = int(g[2])
	
	Match3 = re.search("maxe ([^;]*)", Attrs)
	if Match3 == None:
		MaxError = -1
	else:	
		a = Match3.groups()[0]
		MaxError = float(a)
	return 1

def AppendAttr(a):
	global Attrs
	if len(Attrs) == 0:
		Attrs = a
		return
	if not Attrs.endswith(";"):
		Attrs += "; "
	else:
		Attrs += " "
	Attrs += a
	if not a.endswith(";"):
		Attrs += ";"

def AppendStrAttr(Name, Value):
	AppendAttr(Name + ' "' + Value + '"')

def AppendIntAttr(Name, Value):
	AppendAttr(Name + " " + str(Value))

def ParseRepeatAttrs():
	global Attrs
	global RepeatName
	global RepeatFamily
	global RepeatStart
	global RepeatEnd
	global RepeatLeft

	Match = re.search("Repeat ([^;]*)", Attrs)
	if Match == None:
		print "No match (Repeat): " + Attrs
		return 0

	a = Match.groups()[0]
	
	Match2 = re.match("([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+) ([^ ]+)", a)
	if Match2 == None:
		Quit("Invalid Repeat attribute in GFF file: " + a)
	
	g = Match2.groups()

	RepeatName = g[0]
	RepeatFamily = g[1]
	strRepeatStart = g[2]
	strRepeatEnd = g[3]
	strRepeatLeft = g[4]
	if strRepeatStart == ".":
		if strRepeatEnd != "." or strRepeatLeft != ".":
			Quit("Invalid Repeat attribute in GFF file: " + a)
		RepeatStart = -1
		RepeatEnd = -1
		RepeatLeft = -1
	else:
		RepeatStart = int(strRepeatStart)
		RepeatEnd = int(strRepeatEnd)
		RepeatLeft = int(strRepeatLeft)

	return 1

# GFF record:
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
#     0         1         2        3     4       5       6       7         8           9

# Attrributes: 
# Repeat <RepeatName> <RepeatFamily> <RepeatFrom> <RepeatTo> <RepeatLeft>
#             0             1            2             3            4

def ParseRec(Line):
	global Label
	global Source
	global Feature
	global Start
	global End
	global Score
	global Strand
	global Frame
	global Attrs

	Fields = Line.split("\t")
	if len(Fields) < 8:
		Quit("Expected 8 fields in GFF record, got: " + Line)
	
	Label = Fields[0]
	Source = Fields[1]
	Feature = Fields[2]
	Start = int(Fields[3])
	End = int(Fields[4])
	if Fields[5] == ".":
		Score = 0
	else:
		Score = float(Fields[5])
	Strand = Fields[6]
	Frame = Fields[7]

	Attrs = ""
	if len(Fields) > 8:
		Attrs = Fields[8]
		
def GetRec(File, OnRecord):
	global Line
	Line = File.readline()
	if len(Line) == 0:
		return 0
	Line = Line.strip()
	if len(Line) == 0:
		return 1
	ParseRec(Line)
	OnRecord()
	return 1

def GetRecs(FileName, OnRecord):
	File = open(FileName)
	File.seek(0, 2)
	FileSize = File.tell()
	File.seek(0)
	Counter = 0
	while GetRec(File, OnRecord):
		if Counter%1000 == 0:
			ProgressFile(File, FileSize)
		Counter += 1
		continue
#	if sys.stderr.isatty():
#		sys.stderr.write("100.0%\n")
	
def CompareGFFLines(Line1, Line2):
	ParseRec(Line1)
	Label1 = Label
	Start1 = Start
	End1 = End

	ParseRec(Line2)
	Label2 = Label
	Start2 = Start
	End2 = End
	
	if Label1 == Label2:
		if Start1 == Start2:
			if End1 == End2:
				return 0
			elif End1 < End2:
				return -1
			else:
				return 1
			return 0
		elif Start1 < Start2:
			return -1
		else:
			return 1
	elif Label1 < Label2:
		return -1
	else:
		return 1

def GetSortedLines(FileName):
	SortedLines = []
	File = open(FileName)
	File.seek(0, 2)
	FileSize = File.tell()
	File.seek(0)
	Counter = 0
	while 1:
		if Counter%1000 == 0:
			ProgressFile(File, FileSize)
		Counter += 1
		Line = File.readline()
		if len(Line) == 0:
			break
		SortedLines.append(Line.strip())
	if sys.stderr.isatty():
		sys.stderr.write("100.0%\n")

	print >> sys.stderr, "Sorting %d recs" % len(SortedLines)
	SortedLines.sort(CompareGFFLines)
	return SortedLines

def GetSortedRecs(FileName, OnRecord):
	SortedLines = GetSortedLines(FileName)
	Counter = 0
	N = len(SortedLines)
	for Line in SortedLines:
		if Counter%1000 == 0:
			Progress(Counter, N)
		Counter += 1
		ParseRec(Line)
		OnRecord()
	if sys.stderr.isatty():
		sys.stderr.write("100.0%\n")

# GFF record:
# <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
#     0         1         2        3     4       5       6       7         8           9
def RecToStr():
	return "%s\t%s\t%s\t%d\t%d\t%.5g\t%s\t%s\t%s" % (Label, Source, Feature, Start, End, Score, Strand, Frame, Attrs)
	#        0   1   2   3   4   5   6   7   8       0       1       2       3     4     5       6      7      8

def WriteRec(File):
	s = RecToStr()
	print >> File, s

def GetRecAsDict():
	global Label
	global Source
	global Feature
	global Start
	global End
	global Score
	global Strand
	global Frame
	global Attrs

	r = {}
	r["Label"] = Label
	r["Source"] = Source
	r["Feature"] = Feature
	r["Start"] = Start
	r["End"] = End
	r["Score"] = Score
	r["Strand"] = Strand
	r["Frame"] = Frame
	r["Attrs"] = Attrs
	return r

def isinteger(s):
	if len(s) == 0:
		return 0
	for c in s:
		if not c.isdigit():
			return 0
	return 1

def GetAttrStrFromDict(r):
	s = ""
	for Name in r.keys():
		if s != "":
			s += " "
		Value = r[Name]
		if isinteger(Value):
			s += "%s %s;" % (Name, Value)
		else:
			s += "%s \"%s\";" % (Name, Value)
	return s

def PutRecAsDict(r):
	global Label
	global Source
	global Feature
	global Start
	global End
	global Score
	global Strand
	global Frame
	global Attrs

	Label = r["Label"]
	Source  = r["Source"]
	Feature= r["Feature"]
	Start = r["Start"]
	End = r["End"]
	Score = r["Score"]
	Strand = r["Strand"]
	Frame = r["Frame"]
	Attrs = r["Attrs"]
