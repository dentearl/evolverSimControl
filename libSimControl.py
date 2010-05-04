def standardOptions(parser):
    parser.add_option('-p', '--parent',dest='parentDir',
                      help='Parent directory.')
    parser.add_option('-j', '--jobFile',dest='jobFile',
                      help='jobFile, passed in by jobTree.py.')
    parser.add_option('-s', '--step',dest='stepSize', action="store",
                      type ='float', default=0.001,
                      help='stepSize for each cycle. [default: %default]')
    parser.add_option('-g', '--params',dest='gParamsDir',
                      help='Parameter directory.')
    parser.add_option('-e', '--seed',dest='seed',default='random',
                      type='string', help='Random seed. [default: %default]')

def standardOptionsCheck(options, usage):
    import os, sys, random
    if (options.parentDir == None):
        sys.stderr.write('%s: Error, specify parent dir.\n' % sys.argv[0])
        usage()
    if (options.jobFile == None):
        sys.stderr.write('%s: Error, specify jobFile.\n' % sys.argv[0])
        usage()
    if (options.gParamsDir == None):
        sys.stderr.write('%s: Error, specify params.\n' % sys.argv[0])
        usage()
    if (not os.path.isdir(options.parentDir)):
        sys.stderr.write('%s: Error, Parent dir "%s" not a directory!\n' % (sys.argv[0], options.parentDir))
        usage()
    if (not os.path.isdir(options.gParamsDir)):
        sys.stderr.write('%s: Error, Params dir "%s" not a directory!\n' % (sys.argv[0], options.gParamsDir))
        usage()
    if(not os.path.exists(options.jobFile)):
        sys.stderr.write('%s: Error, jobFile does not exist.\n' % sys.argv[0])
        usage()
    options.parentDir  = os.path.abspath(options.parentDir)
    options.gParamsDir = os.path.abspath(options.gParamsDir)
    options.theParent  = os.path.basename(options.parentDir)

def verifyPrograms(programs):
    """verifyPrograms(programs) takes a list of executable names, and acts on the list object
    to look up the full path to the executables, or if they are not found, to print an
    error message and sys.exit(1)
    """
    import sys
    c=-1
    for i in programs:
        c=c+1
        i = which(i)
        if(i == None):
            sys.stderr.write('ERROR: Could not locate "%s" in PATH.\n' %(programs[c]))
            sys.exit(1)
        else:
            programs[c] = i

def which(program):
    """which() acts like the unix utitily which, but is portable between os.
    If the program does not exist in the PATH then 'None' is returned. 
    """
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None

def commandPacker(c):
    """commandPacker() takes in a string and replaces all of the characters
    that need to be escaped when being passed to jobTree. The ultimate goal
    is for this command to be passed to simCtrl_commandEval.py where it will
    be unpacked and then executed.
    """
    a = c.replace('"', '&myQuot;')
    a = a.replace("'", "&myApos;")
    a = '&myCMD;'+a+'&myCMD;'
    return a


def discritizeTree(nt, ss):
    """discritizeTree() takes a newickTree (binaryTree object) and a step size
    and translates the tree branch distances into discrete steps using a ceiling
    function. So a branch length of 0.9 and a step size of 0.4 yields a new branch
    length of 3 (math.ceil(0.9/0.4) = 3.0). Trees are changed in place.
    """
    import math
    if(nt==None):
        return
    nt.distance=math.ceil(float(nt.distance)/ss)
    discritizeTree(nt.right, ss)
    discritizeTree(nt.left, ss)

def typeTimestamp(file, type, value):
    """file is a cycleInfo.xml file, type is in {main, stats}, value is in {Start, End}
    """
    import xml.etree.ElementTree as ET
    import time
    infoTree=ET.parse(file)
    root=infoTree.getroot()
    if (type == 'main' and value == 'Start'):
        # create new xml
        tObj=ET.SubElement(root, 'timestamps')
        tObj.attrib['startEpochUTC'] = str(time.time())
        cycle=ET.SubElement(tObj,type)
    elif value == 'Start':
        # create 'stats' section
        tObj=root.find('timestamps')
        cycle=ET.SubElement(tObj,type)
    else:
        # continue with the present type
        tObj=root.find('timestamps')
        cycle=tObj.find(type)

    timeStart=ET.SubElement(cycle,type+value)
    timeHuman=ET.SubElement(timeStart, 'humanUTC')
    timeHuman.text=str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime()))
    timeEpoch=ET.SubElement(timeStart, 'epochUTC')
    timeEpoch.text=str(time.time())
    if value== 'End' and type == 'stats':
        # this is the end of the cycle
        tObj=root.find('timestamps')
        tObj.attrib['endEpochUTC'] = str(time.time())
    info=ET.ElementTree(root)
    info.write(file)

def subTypeTimestamp(file, type, timeName):
    """file is a cycleInfo.xml file, type is in {main, stats}, timeName is something
    like 'cyleStep_1_cycleMain_1_start' or 'cyleStep_1_cycleMain_1_end'
    """
    import xml.etree.ElementTree as ET
    import time
    infoTree=ET.parse(file)
    root=infoTree.getroot()
    tObj=root.find('timestamps')
    cycle=tObj.find(type)
    timeObj=ET.SubElement(cycle, timeName)
    timeHuman=ET.SubElement(timeObj, 'humanUTC')
    timeHuman.text=str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC) ", time.gmtime()))
    timeEpoch=ET.SubElement(timeObj, 'epochUTC')
    timeEpoch.text=str(time.time())
    info=ET.ElementTree(root)
    info.write(file)

def stem_print(close, dist, ndigits):
   """/*
   *  R : A Computer Language for Statistical Data Analysis
   *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
   *  Copyright (C) 1997-2000   Robert Gentleman, Ross Ihaka and the
   *                            R Development Core Team
   #  This function ported to Python by Dent Earl, UCSC BME Dept. 2010
   *
   *  This program is free software; you can redistribute it and/or modify
   *  it under the terms of the GNU General Public License as published by
   *  the Free Software Foundation; either version 2 of the License, or
   *  (at your option) any later version.
   *
   *  This program is distributed in the hope that it will be useful,
   *  but WITHOUT ANY WARRANTY; without even the implied warranty of
   *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   *  GNU General Public License for more details.
   *
   *  You should have received a copy of the GNU General Public License
   *  along with this program; if not, a copy is available at
   *  http://www.r-project.org/Licenses/
   */
   """
   import sys
   if (close/10 == 0) and (dist < 0):
      sys.stdout.write( '  %*s | ' %(ndigits, '-0'))
   else:
      sys.stdout.write('  %*d | ' %(ndigits, int(close/10)))

def stem(data, scale=1, width=80, atom=1e-8):
   """/*
   *  R : A Computer Language for Statistical Data Analysis
   *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
   *  Copyright (C) 1997-2000   Robert Gentleman, Ross Ihaka and the
   *                            R Development Core Team
   #  This function ported to Python by Dent Earl, UCSC BME Dept. 2010
   *
   *  This program is free software; you can redistribute it and/or modify
   *  it under the terms of the GNU General Public License as published by
   *  the Free Software Foundation; either version 2 of the License, or
   *  (at your option) any later version.
   *
   *  This program is distributed in the hope that it will be useful,
   *  but WITHOUT ANY WARRANTY; without even the implied warranty of
   *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   *  GNU General Public License for more details.
   *
   *  You should have received a copy of the GNU General Public License
   *  along with this program; if not, a copy is available at
   *  http://www.r-project.org/Licenses/
   */
   """
   import math, sys
   if len(data) <=1:
      return False
   data = sorted(data)
   if data[-1] > data[0]:
      r = atom+(data[-1]-data[0])/scale
      c = 10.0**(11.0 - int(math.log10(r)+10))
      mm = min(2, max(0, int(r*c/25)))
      k = 3*mm + 2 - 150/(len(data)+50)
      if (k-1)*(k-2)*(k-5)==0:
         c = c*10.0
      x1 = abs(data[0])
      x2 = abs(data[-1]);
      if(x2 > x1):
         x1 = x2;
      while(x1*c > sys.maxint):
         c /= 10
      if (k*(k-4)*(k-8)==0):
         mu = 5
      if ((k-1)*(k-5)*(k-6)==0):
         mu = 20
   else:
      r = atom + abs(data[0])/scale;
      c = 10.0**(11.0-int(math.log10(r)+10))
      k = 2 #/* not important what   */
   
   mu = 10
   if (k*(k-4)*(k-8))==0:
      mu = 5
   if ((k-1)*(k-5)*(k-6)==0):
      mu = 20
      
   # Find and print width of the stem.
   
   lo = math.floor(data[0]*c/mu)*mu
   hi = math.floor(data[-1] *c/mu)*mu
   ldigits = int(math.floor(math.log10(-lo))+1) if(lo < 0) else 0
   hdigits = int(math.floor(math.log10(hi))) if(hi > 0) else 0
   ndigits = int(hdigits) if(ldigits < hdigits) else ldigits
   
   # Starting cell
   
   if (lo < 0) and (math.floor(data[0]*c) == lo):
      lo = lo - mu
   hi = lo + mu
   if math.floor(data[0]*c+0.5) > hi:
      lo = hi
      hi = lo+mu
   # Print decimal info
   pdigits = 1 - math.floor(math.log10(c)+0.5)
   decStr = '\n  The decimal point is '
   if pdigits == 0:
      decStr = decStr + 'at the |\n'
   else:
      direction = 'right' if pdigits > 0 else 'left'
      decStr = decStr + '%d digit(s) to the %s of the |\n'%(pdigits, direction)
   print decStr
   i=0
   while True:
      if lo < 0:
         stem_print(int(hi), int(lo), ndigits)
      else:
         stem_print(int(lo), int(hi), ndigits)
      j=0
      while i < len(data):
         if data[i] < 0:
            xi = data[i]*c - .5
         else:
            xi = data[i]*c + .5
         if (hi == 0 and data[i] >=0) or (lo<0 and xi > hi) or (lo >= 0 and xi >=hi):
            break
         j+=1
         if (j<=width-12):
            sys.stdout.write( '%d' %(abs(xi)%10))
         i+=1
      if j > width:
         sys.stdout.write( '+%d' %(j-width))
      sys.stdout.write('\n')
      if i >=len(data):
         break
      hi += mu
      lo += mu
   sys.stdout.write('\n')
