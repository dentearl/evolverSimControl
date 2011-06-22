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
