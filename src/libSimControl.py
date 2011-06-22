##################################################
# Copyright (C) 2009-2011 by
# Dent Earl (dearl@soe.ucsc.edu, dentearl@gmail.com)
# Benedict Paten (benedict@soe.ucsc.edu, benedictpaten@gmail.com)
# ... and other members of the Reconstruction Team of David Haussler's
# lab (BME Dept. UCSC)
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
##################################################
requiredPrograms = ['cat', 'cp', 'egrep', 'ln', 'mkdir', 
                    'mv', 'sed', 'touch', 'trf',
                    'evolver_codon_report.pl', 
                    'evolver_cvt', 
                    'evolver_drawrev', 
                    'evolver_evo', 
                    'evolver_evostats_report.py', 
                    'evolver_gene_deactivate.sh', 
                    'evolver_gff_cdsutr2exons.py', 
                    'evolver_gff_exons2introns.py', 
                    'evolver_gff_featurestats2.py', 
                    'evolver_gff_featurestats2.sh', 
                    'evolver_gff_fixgeneix.py', 
                    'evolver_handle_mobiles.pl', 
                    'evolver_merge_evostats.py', 
                    'evolver_mobile_report.pl', 
                    'evolver_transalign', 
                    'evolver_trf2gff.py']

def verifyPrograms(programs):
    """verifyPrograms(programs) takes a list of executable names, and acts on the list object
    to look up the full path to the executables, or if they are not found it raises an exeption
    """
    from libSimControlClasses import BadInputError
    from libSimControl import which
    if not isinstance(programs, list):
       raise TypeError('verifyPrograms takes a list of program '
                       'names, not %s.\n' % programs.__class__)
    c=-1
    for p in programs:
       if not isinstance(p, str):
          raise TypeError('verifyPrograms list members should all be strings, '
                          '"%s" not a string, is a %s.\n' %(str(p), p.__class__))
       c=c+1
       p = which(p)
       if p is None:
           raise RuntimeError('Error verifyPrograms(): Could not locate "%s" '
                              'in PATH.\n' %(programs[c]))
       else:
           programs[c] = p

def which(program):
    """which() acts like the unix utility which, but is portable between os.
    If the program does not exist in the PATH then 'None' is returned. 
    """
    import os
    def is_exe(fpath):
        return os.path.exists(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath != '':
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def discritizeTree(nt, ss):
    """discritizeTree() takes a newickTree (binaryTree object) and a step size
    and translates the tree branch distances into discrete steps using a ceiling
    function. So a branch length of 0.9 and a step size of 0.4 yields a new branch
    length of 3 (math.ceil(0.9/0.4) = 3.0). Trees are changed in place.
    """
    import math
    if nt is None :
        return
    nt.distance=math.ceil(float(nt.distance)/ss)
    discritizeTree(nt.right, ss)
    discritizeTree(nt.left, ss)

def typeTimestamp(dirname, typeTS, value):
    """dirname is a cycle directory, typeTS is in {cycle, stats, transalign}, value is in {start, end}
    """
    from libSimControlClasses import BadInputError
    from libSimControl import lockfile, unlockfile, addTimestampsTag
    import os
    import sys
    import time
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat
    value = value.lower()
    if typeTS not in ('cycle', 'stats', 'transalign'):
        raise BadInputError('typeTS must be either "cycle", "stats", or "transalign" not %s\n' % typeTS)
    if value not in ('start', 'end'):
        raise BadInputError('value must be either "start" or "end", not %s\n' % value)
    value = value[0].upper() + value[1:]
    fileMap = {'cycle':'cycle', 'stats':'stats', 'transalign':'transalign'}
    filename = os.path.join(dirname, 'xml', fileMap[typeTS] + '.xml')
    if value == 'Start':
        # new xml, needs new timestamp tag
        addTimestampsTag(filename)
    lockname = lockfile(filename)
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    root = infoTree.getroot()
    timeTag = root.find('timestamps')
    timeStart = ET.SubElement(timeTag, typeTS + value)
    timeHuman = ET.SubElement(timeStart, 'humanUTC')
    timeHuman.text = str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime()))
    timeEpoch = ET.SubElement(timeStart, 'epochUTC')
    timeEpoch.text = str(time.time())
    info = ET.ElementTree(root)
    info.write(lockname)
    unlockfile(lockname)

def addTimestampsTag(filename):
    """
    """
    from libSimControl import lockfile, unlockfile
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat
    import time
    lockname = lockfile(filename)
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    root = infoTree.getroot()
    if len(root.findall('timestamps')) > 0:
        raise RuntimeError('There should be no timestamps tag\n')
    timeTag = ET.SubElement(root, 'timestamps')
    timeTag.attrib['startEpochUTC'] = str(time.time())
    info = ET.ElementTree(root)
    info.write(lockname)
    unlockfile(lockname)

def subTypeTimestamp(dirname, typeTS, timeName, chrName = None):
    """dirname is the cycle directory, type is in {cycle, stats, transalign}, timeName is something
    like 'CycleStep1_start' or 'cycleStep1_end'
    """
    from libSimControlClasses import BadInputError
    from libSimControl import lockfile, unlockfile
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat
    import time
    if typeTS not in ('cycle', 'stats', 'transalign', 'cycleChr'):
        raise BadInputError('typeTS must be either "cycle", "stats", '
                             '"transalign", or "cycleChr" not %s\n' % typeTS)
    if chrName is not None:
        filename = os.path.join(dirname, 'xml', 'cycle.%s.xml' % chrName)
    else:
        filename = os.path.join(dirname, 'xml', typeTS + '.xml')
    lockname = lockfile(filename)
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    root = infoTree.getroot()
    timeTag = root.find('timestamps')
    timeObj = ET.SubElement(timeTag, timeName)
    timeHuman = ET.SubElement(timeObj, 'humanUTC')
    timeHuman.text = str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC) ", time.gmtime()))
    timeEpoch = ET.SubElement(timeObj, 'epochUTC')
    timeEpoch.text = str(time.time())
    info = ET.ElementTree(root)
    info.write(lockname)
    unlockfile(lockname)

def lockfile(filename):
    """ This is a fragile attempt at avoiding too many collisions on the xml info files
    Many processes will try to lock simultaneously.
    """
    import os
    import time
    timeout = 60
    for numFails in xrange(0, timeout):
        try:
            # shutil.move:
            # shutil.move is atomic in linux when both source and dest are on same filesystem
            # shutil.move(filename, filename + '.lock')
            # 
            # os.rename:
            # The operation may fail on some Unix flavors if src and dst are on different 
            # filesystems. If successful, the renaming will be an atomic operation (this 
            # is a POSIX requirement).
            os.rename(filename, filename + '.lock')
            break
        except:
            time.sleep(1)
    if not os.path.exists(filename + '.lock'):
        raise RuntimeError('Unable to lock file %s after %d seconds!\n' % (filename, timeout))
    return filename + '.lock'

def unlockfile(filename):
    """ There will be only one process seeking to unlock a file at a time.
    """
    import os
    import os
    if filename.endswith('.lock'):
        try:
            os.rename(filename, filename[:-5])
        except:
            raise RuntimeError('Database access rate is likely too high for filesystem, %s\n' % filename)
    else:                                                                                            
        raise RuntimeError('unlockfile: Supplied filename, %s, does not end in ".lock"\n' % filename)

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
        sys.stdout.write('  %*s | ' %(ndigits, '-0'))
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
        if x2 > x1:
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
    ldigits = int(math.floor(math.log10(-lo))+1) if (lo < 0) else 0
    hdigits = int(math.floor(math.log10(hi))) if (hi > 0) else 0
    ndigits = int(hdigits) if (ldigits < hdigits) else ldigits
   
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
            j += 1
            if (j<=width-12):
                sys.stdout.write('%d' %(abs(xi)%10))
            i += 1
        if j > width:
            sys.stdout.write('+%d' %(j-width))
        sys.stdout.write('\n')
        if i >=len(data):
            break
        hi += mu
        lo += mu
    sys.stdout.write('\n')

def nameTree(nt, reportDistance = True):
    """nameTree(nt) takes a newick tree and returns a str that can be used
    to name the cycle-step that the tree represents. Distance included in
    the name by default.
    """
    from libSimControl import sanitizeTreeName
    from sonLib.bioio import printBinaryTree
    if nt is None:
        print 'nt was none'
        return ''
    if nt.iD is not None:
        if nt.distance == 0.0 or not reportDistance:
            name = nt.iD
        else:
            name = nt.iD + str(nt.distance)
    else:
        name = printBinaryTree(nt, True)
    name = sanitizeTreeName(name)
    return name

def sanitizeTreeName(name):
    """sanitizeTreeName(name) takes all the nasty characters out of a newickTree and
    returns a str that is more amenable to being a file (or directory) name.
    """
    name = name.replace(' ','')
    name = name.replace(',','')
    name = name.replace(':','-')
    name = name.replace('.','_')
    name = name.replace(';','')
    name = name.replace('\'','')
    name = name.replace('"','')
    name = name.replace('(','_L_')
    name = name.replace(')','_R_')
    name = name.rstrip('0')
    name = name.rstrip('-0_')
    return name

def tree2str(nt):
    """ tree2str takes a newick tree object and returns an 
    unsanitized string containing node distances
    """
    from sonLib.bioio import printBinaryTree
    ts = printBinaryTree(nt, True)
    return ts.rstrip(':.0;') # necessary due to weird newickTree code

def takeNewickStep(thisNewickStr, options):
    """ takeNewickStep takes a newick tree and an options object
    and returns the name of the next step.
    """ 
    from libSimControl import tree2str
    from sonLib.bioio import newickTreeParser
    nt = newickTreeParser(thisNewickStr, 0.0)
    if nt.distance > options.stepLength:
        nt.distance -= options.stepLength
        step = options.stepLength
    else:
        step = nt.distance
        nt.distance = 0.0
    return (tree2str(nt), step)

def myLog(s):
    import os
    if os.path.exists('sc_log.log'):
        f = open('sc_log.log', 'a')
    else:
        f = open('sc_log.log', 'w')
        f.write(s)
    f.close()

def runCommands(cmds, localTempDir, inPipes = [], outPipes = [], mode='s', debug = False):
    """ runCommands is a wrapper function for the parallel and serial
    versions of runCommands. mode may either be s or p.
    """
    from libSimControl import runCommandsP, runCommandsS
    import os
    from sonLib.bioio import logger
    if not os.path.exists(localTempDir):
        raise ValueError('localTempDir "%s" does not exist!\n' % localTempDir)
    if not isinstance(cmds, list):
        raise TypeError('runCommands takes a list for the "cmds" '
                        'argument, not a %s.\n' % cmds.__class__)
    if mode not in ('s', 'p'):
        raise ValueError('runCommands "mode" argument must be either '
                         's or p, not %s.\n' % mode)
    if outPipes != []:
        if len(cmds) != len(outPipes):
            raise ValueError('runCommands length of outPipes list %d '
                             'not equal to cmds list %d!.\n' % (len(outPipes), len(cmds)))
    else:
        outPipes = [None] * len(cmds)
    if inPipes != []:
        if len(cmds) != len(inPipes):
            raise ValueError('runCommands length of inPipes list %d '
                             'not equal to cmds list %d!.\n' % (len(inPipes), len(cmds)))
    else:
        inPipes = [None] * len(cmds)
    if mode == 's':
        logger.info('Issuing serial commads %s %s %s.\n' % (str(cmds), str(inPipes), str(outPipes)))
        runCommandsS(cmds, localTempDir, inPipes = inPipes, outPipes = outPipes, debug = debug)
    else:
        logger.info('Issuing parallel commads %s %s %s.\n' % (str(cmds), str(inPipes), str(outPipes)))
        runCommandsP(cmds, localTempDir, inPipes = inPipes, outPipes = outPipes, debug = debug)

def runCommandsP(cmds, localTempDir, inPipes=[], outPipes=[], debug = False):
    """ runCommandsP uses the subprocess module
    to issue parallel processes from the cmds list.
    """
    import os
    from sonLib.bioio import logger
    import subprocess
    procs = []
    i = -1
    for c in cmds:
        i += 1
        if inPipes[i] is None:
            sin = None
        else:
            sin = subprocess.PIPE
        if outPipes[i] is None:
            sout = None
        else:
            sout = subprocess.PIPE
        logger.debug('Executing parallel %s < %s > %s\n' % (' '.join(c), inPipes[i], outPipes[i]))
        procs.append(subprocess.Popen(c, cwd=localTempDir, stdin = sin, stdout = sout))
    i = -1
    for p in procs:
        i += 1
        if inPipes[i] is None:
            sin = None
        else:
            if not os.path.exists(inPipes[i]):
                raise IOError('Unable to locate inPipe file: %s for command %s' % (inPipes[i], cmds[i]))
            sin = open(inPipes[i], 'r').read()
        if outPipes[i] is None:
            pout, perr = p.communicate(sin)
            handleReturnCode(p.returncode, cmds[i])
        else:
            f = open(outPipes[i], 'w')
            f.write(p.communicate(sin)[0])
            f.close()
            handleReturnCode(p.returncode, cmds[i])

def runCommandsS(cmds, localTempDir, inPipes=[], outPipes=[], debug = False):
    """ runCommandsS uses the subprocess module
    to issue serial processes from the cmds list.
    """
    import os
    from sonLib.bioio import logger
    import subprocess
    i = -1
    for c in cmds:
        i += 1
        if inPipes[i] is None:
            sin = None
        else:
            sin = subprocess.PIPE
        if outPipes[i] is None:
            sout = None
        else:
            sout = subprocess.PIPE
        logger.debug('Executing serial %s < %s > %s\n' % (' '.join(c), inPipes[i], outPipes[i]))
        p = subprocess.Popen(c, cwd=localTempDir, stdin = sin, stdout = sout)
            
        if inPipes[i] is None:
            sin = None
        else:
            if not os.path.exists(inPipes[i]):
                raise IOError('Unable to locate inPipe file: %s for command %s' % (inPipes[i], c))
            sin = open(inPipes[i], 'r').read()
        if outPipes[i] is None:
            pout, perr = p.communicate(sin)
            handleReturnCode(p.returncode, cmds[i])
        else:
            f = open(outPipes[i], 'w')
            f.write(p.communicate(sin)[0])
            f.close()
            handleReturnCode(p.returncode, cmds[i])

def handleReturnCode(retcode, cmd):
    if not isinstance(retcode, int):
        raise TypeError('handleReturnCode takes an integer for '
                        'retcode, not a %s.\n' % retcode.__class__)
    if not isinstance(cmd, list):
        raise TypeError('handleReturnCode takes a list for '
                        'cmd, not a %s.\n' % cmd.__class__)
    if retcode:
        if retcode < 0:
            raise RuntimeError('Experienced an error while trying to execute: '
                               '%s SIGNAL:%d\n' %(' '.join(cmd), -retcode))
        else:
            raise RuntimeError('Experienced an error while trying to execute: '
                               '%s retcode:%d\n' %(' '.join(cmd), retcode))

def createNewCycleXmls(directory, parentDir, stepLength, newickStr, options):
    """
    """
    from libSimControl import typeTimestamp, subTypeTimestamp, newInfoXml, tree2str, takeNewickStep
    import os
    from sonLib.bioio import newickTreeParser
    import xml.etree.ElementTree as ET
    if not os.path.exists(directory):
        raise RuntimeError('cycleNewCycleInfoXml: directory: %s does not exist!\n' % directory)
    if not os.path.isdir(directory):
        raise RuntimeError('cycleNewCycleInfoXml: directory: %s is not a directory!\n' % directory)
    if not os.path.exists(os.path.join(directory, 'xml', 'summary.xml')):
        root = ET.Element('info')
        e = ET.SubElement(root, 'parentDir')  
        e.text = parentDir
        e = ET.SubElement(root, 'thisDir')
        e.text = directory
        e = ET.SubElement(root, 'stepLength')
        e.text = str(stepLength).rstrip('0')
        nt = newickTreeParser(newickStr, 0.0)
        children = {}
        if nt.distance == 0:
            if nt.internal:
                branches = {'left' : tree2str(nt.left),
                            'right': tree2str(nt.right)}
                for b in branches:
                    children[b] = nameTree(newickTreeParser(takeNewickStep(branches[b], options)[0], 0.0))
        else:
            children['stem'] = nameTree(newickTreeParser(takeNewickStep(tree2str(nt), options)[0], 0.0))
        e = ET.SubElement(root, 'numberChildren')
        e.text = str(len(children))
        for c in children:
            e = ET.SubElement(root, 'child')
            e.text = children[c]
            e.attrib['type'] = c # left, right, stem

        info = ET.ElementTree(root)
        info.write(os.path.join(directory, 'xml', 'summary.xml'))
        addTimestampsTag(os.path.join(directory, 'xml', 'summary.xml'))
    if not os.path.exists(os.path.join(directory, 'xml', 'cycle.xml')):
        newInfoXml(os.path.join(directory, 'xml', 'cycle.xml'))
        typeTimestamp(directory, 'cycle', 'start')

def newInfoXml(filename):
    """
    """
    import xml.etree.ElementTree as ET
    root=ET.Element('info')
    info=ET.ElementTree(root)
    info.write(filename)

def createRootXmls(command, options):
    """
    """
    from libSimControl import createSimulationInfoXml
    import os
    import xml.etree.ElementTree as ET

    os.mkdir(os.path.join(options.outDir, options.rootName, 'xml'))
    root = ET.Element('info')
    e = ET.SubElement(root, 'cycleIsRoot')
    e.text = 'True'
    info = ET.ElementTree(root)
    info.write(os.path.join(options.outDir, options.rootName, 'xml', 'summary.xml'))
    createSimulationInfoXml(command, options)

def createSimulationInfoXml(command, options):
    """
    """
    import os
    import time
    import xml.etree.ElementTree as ET
    if(os.path.exists(os.path.join(options.outDir, 'simulationInfo.xml'))):
        raise RuntimeError('simulationInfo.xml already exists at %s' % options.outDir)
    root = ET.Element('info')
    tObj = ET.SubElement(root, 'sourceRootDir')
    tObj.text = str(options.rootInputDir)
    tObj = ET.SubElement(root, 'sourceParamsDir')
    tObj.text = str(options.paramsInputDir)
    tObj = ET.SubElement(root, 'rootDir')
    tObj.text = str(options.rootDir)
    tObj = ET.SubElement(root, 'paramsDir')
    tObj.text = str(options.paramsDir)
    tObj = ET.SubElement(root, 'rootName')
    tObj.text = str(options.rootName)
    tObj = ET.SubElement(root, 'tree')
    tObj.text = str(options.inputNewick)
    tObj = ET.SubElement(root, 'stepLength')
    tObj.text = str(options.stepLength)
    timeTag = ET.SubElement(root, 'timestamps')
    timeTag.attrib['startEpochUTC'] = str(time.time())
    timeStart = ET.SubElement(timeTag,'start')
    timeLocal = ET.SubElement(timeStart, 'humanLocal')
    timeLocal.text = str(time.strftime("%a, %d %b %Y %H:%M:%S (%Z) ", time.localtime()))
    timeHuman = ET.SubElement(timeStart, 'humanUTC')
    timeHuman.text = str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC) ", time.gmtime()))
    timeEpoch = ET.SubElement(timeStart, 'epochUTC')
    timeEpoch.text = str(time.time())
    cmd = ET.SubElement(root, 'command')
    cmd.text = ' '.join(command)
    info = ET.ElementTree(root)
    info.write(os.path.join(options.outDir,'simulationInfo.xml'))

def verifyDirExists(directory):
    """ Convenience function to verify the existence of a directory
    """
    import os
    if not os.path.exists(directory):
        raise RuntimeError('Error, unable to locate directory %s\n' % directory)
    if not os.path.isdir(directory):
        raise RuntimeError('Error, directory %s is not a directory\n' % directory)

def verifyFileExists(filename):
    """ Convenience function to verify the existence of a file
    """
    import os
    if not os.path.exists(filename):
        raise RuntimeError('Error, unable to locate file %s\n' % filename)

def evolverInterStepCmd(thisDir, thisParentDir, theChild, thisStepLength, seed, paramsDir):
    """ produces the command argument list needed to run an evolver inter step.
    Called by CycleStep1.
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir, paramsDir, os.path.join(thisDir, 'inter'),
              os.path.join(thisDir, 'stats'), os.path.join(thisDir, 'logs')]:
        verifyDirExists(d)
    for f in [os.path.join(thisParentDir,'seq.rev'), os.path.join(thisParentDir, 'annots.gff'),
              os.path.join(paramsDir, 'model.txt')]:
        verifyFileExists(f)

    outname = os.path.join(thisDir, 'inter','inter.outseq.rev')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-interchr')
        cmd.append(os.path.join(thisParentDir, 'seq.rev'))
        cmd.append('-inannots')
        cmd.append(os.path.join(thisParentDir, 'annots.gff'))
        cmd.append('-aln')
        cmd.append(os.path.join(thisDir, 'inter','inter.aln.rev'))
        cmd.append('-outchrnames')
        cmd.append(os.path.join(thisDir, 'inter','inter.chrnames.txt'))
        cmd.append('-outannots')
        cmd.append(os.path.join(thisDir, 'inter','inter.outannots.gff'))
        cmd.append('-outseq')
        cmd.append(outname + '.tmp')
        cmd.append('-outgenome')
        cmd.append(theChild + '.inter')
        cmd.append('-branchlength')
        cmd.append(str(thisStepLength))
        cmd.append('-statsfile')
        cmd.append(os.path.join(thisDir, 'stats', 'inter.stats.txt'))
        cmd.append('-model')
        cmd.append(os.path.join(paramsDir, 'model.txt'))
        cmd.append('-seed')
        cmd.append(str(seed))
        cmd.append('-logevents')
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'logs', 'inter.log'))
        followCmd = [which('mv')]
        followCmd.append(outname + '.tmp')
        followCmd.append(outname)
    return cmd, followCmd

def evolverInterStepMobilesCmd(thisDir, thisParentDir, theParent, thisStepLength, paramsDir):
    """ produces the command argument list needed to run an evolver inter step mobiles command.
    Called by runEvolverInterCmds(). This is only run when options.noMEs is False.
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir, paramsDir, os.path.join(thisDir, 'mobiles'),
              os.path.join(thisParentDir, 'mobiles')]:
        verifyDirExists(d)
    for f in [os.path.join(thisParentDir,'mobiles', 'ME.fa'), 
              os.path.join(thisParentDir,'mobiles', 'ME.gff'),
              os.path.join(thisParentDir,'mobiles', 'LTR.fa'), 
              os.path.join(paramsDir, 'mes.cfg'),
              os.path.join(paramsDir, 'model.mes.txt')]:
        verifyFileExists(f)

    cmd = []
    cmd.append(which('evolver_handle_mobiles.pl'))
    cmd.append('--evo')
    cmd.append(which('evolver_evo'))
    cmd.append('--cvt')
    cmd.append(which('evolver_cvt'))
    cmd.append('--py')
    cmd.append(os.path.dirname(which('evolver_evo')))
    cmd.append('--parentDir')
    cmd.append(thisParentDir)
    cmd.append('--genome')
    cmd.append(theParent)
    cmd.append('--stepSize')
    cmd.append(str(thisStepLength))
    cmd.append('--mefa')
    cmd.append(os.path.join(thisParentDir, 'mobiles', 'ME.fa'))
    cmd.append('--megff')
    cmd.append(os.path.join(thisParentDir, 'mobiles', 'ME.gff'))
    cmd.append('--ltr')
    cmd.append(os.path.join(thisParentDir, 'mobiles', 'LTR.fa'))
    cmd.append('--mescfg')
    cmd.append(os.path.join(paramsDir, 'mes.cfg'))
    cmd.append('--model')
    cmd.append(os.path.join(paramsDir, 'model.mes.txt'))
    return cmd

def evolverInterStepMobilesMoveCmd(thisLocalTempDir, thisDir):
    """ produces the command argument list needed to run an evolver inter step mobiles command.
    Called by CycleStep1.
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisLocalTempDir]:
        verifyDirExists(d)
    for f in [os.path.join(thisLocalTempDir, 'mes.fa'), 
               os.path.join(thisLocalTempDir, 'ME_output.fa'),
               os.path.join(thisLocalTempDir, 'ME_output.gff'), 
               os.path.join(thisLocalTempDir, 'ME_output_ltrs.fa')]:
        verifyFileExists(f)

    cmds = []
    cmds.append([which('mv'), os.path.join(thisLocalTempDir, 'mes.fa'),
                  os.path.join(thisDir, 'mobiles', 'mes.fa')])
    cmds.append([which('mv'), os.path.join(thisLocalTempDir, 'ME_output.fa'),
                  os.path.join(thisDir, 'mobiles', 'ME.fa')])
    cmds.append([which('mv'), os.path.join(thisLocalTempDir, 'ME_output.gff'),
                  os.path.join(thisDir, 'mobiles', 'ME.gff')])
    cmds.append([which('mv'), os.path.join(thisLocalTempDir, 'ME_output_ltrs.fa'),
                  os.path.join(thisDir, 'mobiles', 'LTR.fa')])
    return cmds

def evolverIntraStepCmd(thisDir, theChild, thisStepLength, thisChr, 
                        seed, paramsDir, localTempDir, options):
    """ produces the command argument list needed to run an evolver intra step.
    Called by CycleStep2Chromosome.
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, paramsDir, os.path.join(thisDir, 'intra'),
              os.path.join(thisDir, 'stats'), os.path.join(thisDir, 'logs'), localTempDir]:
        verifyDirExists(d)
    for f in [os.path.join(paramsDir, 'model.txt'), os.path.join(thisDir, 'inter','inter.outseq.rev'),
               os.path.join(thisDir, 'inter', 'inter.outannots.gff')]:
        verifyFileExists(f)

    cmd = []
    cmds = []
    outname = os.path.join(localTempDir, thisChr+'.outseq.rev')
    if not os.path.exists(outname):
        cmd.append(which('evolver_evo'))
        cmd.append('-inseq')
        cmd.append(os.path.join(thisDir, 'inter','inter.outseq.rev'))
        cmd.append('-chrname')
        cmd.append(thisChr)
        cmd.append('-branchlength')
        cmd.append(str(thisStepLength))
        cmd.append('-seed')
        cmd.append(str(seed))
        if not options.noMEs:
            cmd.append('-mes')
            cmd.append(os.path.join(thisDir, 'mobiles', 'mes.fa'))
        cmd.append('-inannots')
        cmd.append(os.path.join(thisDir, 'inter', 'inter.outannots.gff'))
        cmd.append('-statsfile')
        cmd.append(os.path.join(thisDir, 'stats', thisChr+'.stats.txt'))
        cmd.append('-codonsubs')
        cmd.append(os.path.join(thisDir, 'intra', thisChr+'.codonsubs.txt'))
        cmd.append('-outannots')
        cmd.append(os.path.join(thisDir, 'intra', thisChr+'.outannots.gff'))
        cmd.append('-outgenome')
        cmd.append(theChild)
        cmd.append('-model')
        cmd.append(os.path.join(paramsDir, 'model.txt'))
        cmd.append('-aln')
        cmd.append(os.path.join(localTempDir, thisChr+'.aln.rev'))
        cmd.append('-outseq')
        cmd.append(outname + '.tmp')
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'logs', 'intra.'+thisChr+'.log'))
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        cmds.append(cmd)
    return cmds

def evolverIntraStepToFastaCmd(thisDir, thisStepLength, thisChr, paramsDir, localTempDir):
    """ produces the command argument list needed to convert the .rev files into .fa files
    Called by CycleStep2Chromosome.
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, paramsDir, os.path.join(thisDir, 'intra'), localTempDir,
              os.path.join(thisDir, 'stats'), os.path.join(thisDir, 'logs')]:
        verifyDirExists(d)
    for f in [os.path.join(localTempDir, thisChr+'.outseq.rev')]:
        verifyFileExists(f)
    cmd = []
    outname = os.path.join(localTempDir, thisChr+'.outseq.fa')
    if not os.path.exists(outname):
        cmd.append(which('evolver_cvt'))
        cmd.append('-fromrev')
        cmd.append(os.path.join(localTempDir, thisChr+'.outseq.rev'))
        cmd.append('-tofasta')
        cmd.append(outname + '.tmp')
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'intra', 'intra.'+thisChr+'.tofasta.log'))
        cmds = [cmd]
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        cmds.append(cmd)
    return cmds

def callEvolverIntraStepTRFCmd(thisDir, thisChr, localTempDir):
    """ calls tandem repeats finder (trf) on the per chromosome .fa files.
    Called by CycleStep2Chromosome.
    """
    from libSimControl import which, verifyDirExists, verifyFileExists, runCommands
    import os
    import subprocess
    for d in [thisDir, localTempDir]:
        verifyDirExists(d)
    for f in [os.path.join(localTempDir, thisChr+'.outseq.fa')]:
        verifyFileExists(f)
    MAX_PERIOD_SIZE = 2000
    
    # trfBig
    # cmd = []
    # outname = os.path.join(localTempDir, thisChr+'.trf.bed')
    # if not os.path.exists(outname):
    #     cmd.append(which('trfBig'))
    #     cmd.append('-bedAt=%s' % (outname + '.tmp'))
    #     cmd.append('-tempDir=%s' % localTempDir)
    #     cmd.append('-maxPeriod=%d' % MAX_PERIOD_SIZE)
    #     cmd.append(os.path.join(localTempDir, thisChr+'.outseq.fa'))
    #     cmd.append(os.path.join(localTempDir, thisChr+'.outseq.trf.fa'))
    #     cmds = [cmd]
    #     cmd = [which('mv')]
    #     cmd.append(outname + '.tmp')
    #     cmd.append(outname)
    #     cmds.append(cmd)
    # runCommands(cmds, localTempDir)
    
    # trf
    outname = os.path.join(thisDir, 'logs', 'trf.' + thisChr + '.log')
    if not os.path.exists(outname):
        cmd = [which('trf')]
        cmd.append(os.path.join(localTempDir, thisChr+'.outseq.fa'))
        cmd += ['2', '7', '7', '80', '10', '50', str(MAX_PERIOD_SIZE), '-d', '-h']
        returncode = subprocess.call(cmd, cwd=localTempDir)
        # note that TRF's returncode is the number of successfully processed
        # sequences special wrapper.
        if returncode != 1:
            if returncode < 0:
                raise RuntimeError('callEvolverIntraStepTRFCmd: Experienced an error while trying to execute: '
                                   '%s SIGNAL:%d\n' %(' '.join(cmd), -returncode))
            else:
                raise RuntimeError('callEvolverIntraStepTRFCmd: Experienced an error while trying to execute: '
                                   '%s retcode:%d\n' %(' '.join(cmd), returncode))
        f=open(outname, 'w')
        f.close()
    

def evolverIntraStepMoveTRFCmd(thisDir, thisChr, localTempDir):
    """ calls tandem repeats finder (trf) on the per chromosome .fa files.
    Called by CycleStep2Chromosome.
    """
    import glob
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, os.path.join(thisDir, 'intra'), localTempDir]:
        verifyDirExists(d)
    for f in [os.path.join(localTempDir, thisChr + '.outseq.rev'),
               os.path.join(localTempDir, thisChr + '.aln.rev')]:
        verifyFileExists(f)
    
    cmds = []
    cmd = [which('mv')]
    cmd.append(os.path.join(localTempDir, thisChr + '.outseq.rev'))
    cmd.append(os.path.join(thisDir, 'intra', thisChr + '.outseq.rev'))
    cmds.append(cmd)
    cmd = [which('mv')]
    cmd.append(os.path.join(localTempDir, thisChr + '.aln.rev'))
    cmd.append(os.path.join(thisDir, 'intra', thisChr + '.aln.rev'))
    cmds.append(cmd)
    # trfBig
    # cmd = [which('mv')]
    # cmd.append(os.path.join(localTempDir, thisChr+'.trf.bed'))
    # cmd.append(os.path.join(thisDir, 'intra', thisChr+'.trf.bed'))
    # cmds.append(cmd)

    # trf
    files = glob.glob(os.path.join(localTempDir, thisChr + '.*.dat'))
    for f in files:
        cmd = [which('mv')]
        cmd.append(f)
        cmd.append(os.path.join(thisDir, 'intra', os.path.basename(f)))
        cmds.append(cmd)
    return cmds

def runMergeTrfBedsToGff(thisDir):
    """ Takes the bed output of trfBig and converts it to
    a gff that evolver would like to use.
    """
    import glob
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    MAX_MOTIF_ATTR = 32
    for d in [thisDir, os.path.join(thisDir, 'intra')]:
        verifyDirExists(d)
    
    # merge
    outname = os.path.join(thisDir, 'intra', 'trfannots.bed')
    if not os.path.exists(outname):
        files = glob.glob(os.path.join(thisDir, 'intra', '*.trf.bed'))
        out = open(outname + '.tmp', 'w')
        for aFile in files:
            f = open(aFile, 'r')
            for line in f:
                out.write(line)
            f.close()
        out.close()
        os.rename(outname + '.tmp', outname)
    
    # convert bed to evolver gff
    outname = os.path.join(thisDir, 'intra', 'trfannots.gff')
    if not os.path.exists(outname):
        b = open(os.path.join(thisDir, 'intra', 'trfannots.bed'), 'r')
        g = open(outname + '.tmp', 'w')
        for line in b:
            line = line.strip()
            (chrom, start, end, trf, periodSize, nCopies, 
             consensusSize, percentMatch, percentIndel, alignScore, 
             percentA, percentC, percentG, percentT, entropy, motif) = line.split('\t')
            if len(motif) > MAX_MOTIF_ATTR:
                motif = motif[0:MAX_MOTIF_ATTR] + '...'
            g.write('%s\ttrf\ttandem\t%d\t%s\t%s\t'
                    '+\t.\treplen %s; copies %s; cons"%s";\n' % (chrom, int(start)+ 1, end, alignScore,
                                                                 periodSize, nCopies, motif))
        b.close()
        g.close()
        os.rename(outname + '.tmp', outname)

def runEvolverInterCmds(thisDir, thisParentDir, theChild, theParent, thisStepLength, 
                        seed, paramsDir, localTempDir, options):
    from libSimControl import (which, evolverInterStepMobilesCmd, evolverInterStepMobilesMoveCmd, 
                               runCommands, handleReturnCode)
    import os
    import subprocess
    
    outname = os.path.join(thisDir, 'inter', 'inter.outseq.rev')
    if not os.path.exists(outname):
        # we split this up to run it in parallel with the mobiles command
        cmd1, followCmd1 = evolverInterStepCmd(thisDir, thisParentDir, theChild, 
                                               thisStepLength, seed, paramsDir)
        p1 = subprocess.Popen(cmd1, cwd=localTempDir) 
    else:
        p1 = None
        followCmd1 = None
        
    if not options.noMEs:
        outname = os.path.join(thisDir, 'logs', 'mobiles.log')
        if not os.path.exists(outname):
            cmd2 = evolverInterStepMobilesCmd(thisDir, thisParentDir, theParent, thisStepLength, paramsDir)
            p2 = subprocess.Popen(cmd2, cwd=localTempDir, stdout=subprocess.PIPE)
            f = open(outname + '.tmp', 'w')
            f.write(p2.communicate()[0])
            f.close()
            handleReturnCode(p2.returncode, cmd2)
            os.rename(outname + '.tmp', outname)
        mvCmds = evolverInterStepMobilesMoveCmd(localTempDir, thisDir)
        runCommands(mvCmds, localTempDir, mode = 'p')
    if p1 is not None:
        p1.wait()
        handleReturnCode(p1.returncode, cmd1)
        runCommands([followCmd1], localTempDir)

def evolverIntraMergeCmds(thisDir, theChild):
    """ Produces three lists of commands to be run in parallel
    by CycleStep3
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, os.path.join(thisDir, 'intra')]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'inter', 'inter.chrnames.txt')]:
        verifyFileExists(f)
    catCmd = [which('cat')]
    evoCmd = [which('evolver_evo')]
    cvtCmd = [which('evolver_cvt')]
    followCmds = [] # make things atomic
    firstLine = True
    f = open(os.path.join(thisDir, 'inter', 'inter.chrnames.txt'), 'r')
    evoChrStr = ''
    cvtChrStr = ''
    for chrom in f:
        chrom = chrom.strip()
        verifyFileExists(os.path.join(thisDir, 'intra', chrom + 'trfannots.gff'))
        verifyFileExists(os.path.join(thisDir, 'intra', chrom + '.outannots.gff'))
        catCmd.append(os.path.join(thisDir, 'intra', chrom + 'trfannots.gff'))
        catCmd.append(os.path.join(thisDir, 'intra', chrom + '.outannots.gff'))
        verifyFileExists(os.path.join(thisDir, 'intra', chrom + '.aln.rev'))
        verifyFileExists(os.path.join(thisDir, 'intra', chrom + '.outseq.rev'))
        if not firstLine:
            # these commands require comma separated values
            evoChrStr += ','+os.path.join(thisDir, 'intra', chrom + '.aln.rev')
            cvtChrStr += ','+os.path.join(thisDir, 'intra', chrom + '.outseq.rev')
        else:
            evoChrStr += os.path.join(thisDir, 'intra', chrom + '.aln.rev')
            cvtChrStr += os.path.join(thisDir, 'intra', chrom + '.outseq.rev')
            firstLine = False
    f.close()
    
    # this filename comes from a pipe in the CycleStep3() code
    followCmds.append([which('mv'), 
                       os.path.join(thisDir, 'intra', 'evannots.gff.tmp'), 
                       os.path.join(thisDir, 'intra', 'evannots.gff')])
    
    outname = os.path.join(thisDir, 'intra', 'intra.aln.rev')
    evoCmd.append('-mergechrs')
    evoCmd.append(evoChrStr)
    evoCmd.append('-outgenome')
    evoCmd.append(theChild)
    evoCmd.append('-out')
    evoCmd.append(outname + '.tmp')
    followCmds.append([which('mv'), outname + '.tmp', outname])
    
    outname = os.path.join(thisDir, 'seq.rev')
    cvtCmd.append('-mergerevseqs')
    cvtCmd.append(cvtChrStr)
    cvtCmd.append('-out')
    cvtCmd.append(outname + '.tmp')
    followCmds.append([which('mv'), outname + '.tmp', outname])
    
    return (catCmd, evoCmd, cvtCmd, followCmds)

def evolverGeneDeactivationStep(thisDir, thisParentDir):
    """ Produces a list of commands to run the final CycleStep.
    by CycleStep4
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir, os.path.join(thisDir, 'intra')]:
        verifyDirExists(d)
    for f in [os.path.join(thisParentDir, 'annots.gff'),
                   os.path.join(thisDir, 'intra', 'evannots.gff')]:
        verifyFileExists(f)
    cmd = [which('evolver_gene_deactivate.sh')]
    cmd.append(os.path.join(thisParentDir, 'annots.gff'))
    cmd.append(os.path.join(thisDir, 'intra', 'evannots.gff'))
    cmd.append(os.path.join(thisDir, 'annots.gff'))
    cmd.append(which('evolver_evo'))
    return cmd

def statsStep1CmdsP(thisDir, thisParentDir):
    """ Produces a list of commands to run the first stats step.
    called by StatsStep1
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir, os.path.join(thisDir, 'stats')]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'seq.rev'),
               os.path.join(thisDir, 'annots.gff')]:
        verifyFileExists(f)
    
    outPipes = []
    cmds  = []
    followCmds = [] # used for atomicity

    outname = os.path.join(thisDir, 'stats', 'annotstats.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-annotstats')
        cmd.append(os.path.join(thisDir, 'annots.gff'))
        cmd.append('-seq')
        cmd.append(os.path.join(thisDir, 'seq.rev'))
        cmd.append('-log')
        cmd.append(outname + '.tmp')
        cmds.append(cmd)
        outPipes.append(None)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        followCmds.append(cmd)

    outname = os.path.join(thisDir, 'stats', 'probstats.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-probstats')
        cmd.append(os.path.join(thisDir, 'annots.gff'))
        cmd.append('-log')
        cmd.append(outname + '.tmp')
        cmds.append(cmd)
        outPipes.append(None)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        followCmds.append(cmd)

    return cmds, followCmds, outPipes

def statsStep1CmdsS(thisDir, thisParentDir):
    """ Produces a list of commands to run the first stats step.
    called by StatsStep1
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir, os.path.join(thisDir, 'stats')]:
       verifyDirExists(d)
    for f in [os.path.join(thisDir, 'annots.gff'), 
               os.path.join(thisParentDir, 'stats', 'expanded_annots.gff'),
               os.path.join(thisDir, 'seq.rev'), os.path.join(thisParentDir, 'seq.rev')]:
       verifyFileExists(f)

    outPipes = []
    cmds  = []

    outname = os.path.join(thisDir, 'stats', 'cds_annots.gff')
    if not os.path.exists(outname):
        cmd = [which('egrep')]
        cmd.append('CDS|UTR')
        cmd.append(os.path.join(thisDir, 'annots.gff'))
        cmds.append(cmd)
        outPipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        outPipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'exons.gff')
    if not os.path.exists(outname):
        cmd = [which('evolver_gff_cdsutr2exons.py')]
        cmd.append(os.path.join(thisDir, 'stats', 'cds_annots.gff'))
        cmds.append(cmd)
        outPipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        outPipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'introns.gff')
    if not os.path.exists(outname):
        cmd = [which('evolver_gff_exons2introns.py')]
        cmd.append(os.path.join(thisDir, 'stats', 'exons.gff'))
        cmds.append(cmd)
        outPipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        outPipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'expanded_annots.gff')
    if not os.path.exists(outname):
        cmd = [which('cat')]
        cmd.append(os.path.join(thisDir, 'annots.gff'))
        cmd.append(os.path.join(thisDir, 'stats', 'exons.gff'))
        cmd.append(os.path.join(thisDir, 'stats', 'introns.gff'))
        cmds.append(cmd)
        outPipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        outPipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'tmpstats.cycle.diffannots.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_gff_featurestats2.sh')]
        cmd.append(which('evolver_gff_featurestats2.py'))
        cmd.append(os.path.join(thisParentDir, 'stats', 'expanded_annots.gff'))
        cmd.append(os.path.join(thisDir, 'stats', 'expanded_annots.gff'))
        cmd.append(os.path.basename(thisParentDir))
        cmd.append(os.path.basename(thisDir))
        cmd.append(os.path.join(thisParentDir, 'seq.rev'))
        cmd.append(os.path.join(thisDir, 'seq.rev'))
        cmd.append(which('evolver_cvt'))
        cmds.append(cmd)
        outPipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        outPipes.append(None)
        cmds.append(cmd)
    
    return cmds, outPipes

def getParentDir(thisDir):
    """ getParentDir inspects the summary.xml file contained in the supplied directory
    and returns the text of the parentDir tag.
    CAUTION, if thisDir is None, getParentDir() returns None, but if getParentDir() is 
    a string and is not a valid directory, getParentDir() throws an exception.
    The None in None out behavior allows us to easily run Cycle (step n) with Stats (n-1) 
    and Transalign (n-1) in parallel.
    """
    from libSimControl import verifyDirExists, verifyFileExists, lockfile, unlockfile, dirIsRoot
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat

    if thisDir is None:
       return None
    for d in [thisDir]:
       verifyDirExists(d)
    if dirIsRoot(thisDir):
       return None
    lockname = lockfile(os.path.join(thisDir, 'xml', 'summary.xml'))
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    unlockfile(lockname)
    root = infoTree.getroot()
    t = root.find('parentDir')
    if t is not None:
       parentDir = t.text
       return parentDir
    return None

def getBranchDir(thisDir):
    """ Returns the first ancestor that is either the root or a branch point
    """
    from libSimControl import verifyDirExists, isBranchOrRoot, getParentDir
    verifyDirExists(thisDir)
    p = getParentDir(thisDir)
    while not isBranchOrRoot(p):
        p = getParentDir(p)
    return p

def dirIsRoot(thisDir):
    """ dirIsRoot checks to see if the supplied directory is the root directory
    for the simulation. It uses the summary.xml file contained in the supplied
    directory to make the decision. Returns True or False
    """
    from libSimControl import verifyDirExists, lockfile, unlockfile
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat

    for d in [thisDir]:
       verifyDirExists(d)
    lockname = lockfile(os.path.join(thisDir, 'xml', 'summary.xml'))
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    unlockfile(lockname)
    root = infoTree.getroot()
    t = root.find('cycleIsRoot')
    if t is not None:
       if t.text == 'True':
          return True
    return False

def statsStep2Cmds(thisDir, thisParentDir, options):
    """ Produces a list of commands to run the a stats step,
    called by StatsStep2
    """
    import glob
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir, os.path.join(thisDir, 'stats')]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'seq.rev'), os.path.join(thisParentDir, 'seq.rev')]:
        verifyFileExists(f)

    pipes = []
    cmds  = []

    outname = os.path.join(thisDir, 'stats', 'merged_cycle.stats.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_merge_evostats.py')]
        for f in glob.glob(os.path.join(thisDir, 'stats', '*.stats.txt')):
            cmd.append(f)
        cmds.append(cmd)
        pipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'events_cycle.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evostats_report.py')]
        cmd.append(os.path.join(thisDir, 'stats', 'merged_cycle.stats.txt'))
        cmds.append(cmd)
        pipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    if not options.noMEs:
        outname = os.path.join(thisDir, 'stats', 'stats.mobiles.txt')
        if not os.path.exists(outname):
            cmd = [which('evolver_mobile_report.pl')]
            for f in glob.glob(os.path.join(thisDir, 'logs', 'intra.*.log')):
                cmd.append(f)
            cmd.append('--mobilesLog')
            cmd.append(os.path.join(thisDir, 'logs', 'mobiles.log'))
            cmds.append(cmd)
            pipes.append(outname + '.tmp')
            cmd = [which('mv')]
            cmd.append(outname + '.tmp')
            cmd.append(outname)
            pipes.append(None)
            cmds.append(cmd)

    outname = os.path.join(thisDir, 'stats', 'tmpstats.cycle.diffcompost.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-nologcmdlineandtime')
        cmd.append('-compost1')
        cmd.append(os.path.join(thisParentDir, 'seq.rev'))
        cmd.append('-compost2')
        cmd.append(os.path.join(thisDir, 'seq.rev'))
        cmd.append('-log')
        cmd.append(outname + '.tmp')
        cmds.append(cmd)
        pipes.append(None)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    return cmds, pipes

def statsStep3Cmds(thisDir, thisParentDir, options):
    """ Produces a list of commands to run a stats step,
    called by StatsStep3
    """
    from libSimControl import which, verifyDirExists, verifyFileExists, getBranchDir
    import os
    for d in [thisDir, thisParentDir, os.path.join(thisDir, 'stats'), options.rootDir]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'seq.rev'), os.path.join(options.rootDir, 'seq.rev')]:
        verifyFileExists(f)
        
    pipes = []
    cmds  = []

    outname = os.path.join(thisDir, 'stats', 'tmpstats.branch.diffcompost.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-nologcmdlineandtime')
        cmd.append('-compost1')
        cmd.append(os.path.join(getBranchDir(thisDir), 'seq.rev'))
        cmd.append('-compost2')
        cmd.append(os.path.join(thisDir, 'seq.rev'))
        cmd.append('-log')
        cmd.append(outname + '.tmp')
        cmds.append(cmd)
        pipes.append(None)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)

    outname = os.path.join(thisDir, 'stats', 'tmpstats.branch.diffannots.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_gff_featurestats2.sh')]
        cmd.append(which('evolver_gff_featurestats2.py'))
        cmd.append(os.path.join(getBranchDir(thisDir), 'stats', 'expanded_annots.gff'))
        cmd.append(os.path.join(thisDir, 'stats', 'expanded_annots.gff'))
        cmd.append(os.path.basename(options.rootDir))
        cmd.append(os.path.basename(thisDir))
        cmd.append(os.path.join(options.rootDir, 'seq.rev'))
        cmd.append(os.path.join(thisDir, 'seq.rev'))
        cmd.append(which('evolver_cvt'))
        cmds.append(cmd)
        pipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)

    return cmds, pipes

def statsStep4Cmds(thisDir, thisParentDir, options):
    """ Produces a list of commands to run a stats step,
    called by StatsStep4
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir, os.path.join(thisDir, 'stats'), options.rootDir]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'seq.rev'), os.path.join(options.rootDir, 'seq.rev')]:
        verifyFileExists(f)
        
    pipes = []
    cmds  = []
    
    outname = os.path.join(thisDir, 'stats', 'merged_root.stats.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_merge_evostats.py')]
        cmd.append(os.path.join(thisDir, 'stats', 'merged_cycle.stats.txt'))
        cmd.append(os.path.join(thisParentDir, 'stats', 'merged_root.stats.txt'))
        cmds = [cmd]
        pipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'events_root.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evostats_report.py')]
        cmd.append(os.path.join(thisDir, 'stats', 'merged_root.stats.txt'))
        cmds.append(cmd)
        pipes.append(outname + '.tmp')
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    return cmds, pipes

def isBranchOrRoot(thisDir):
    """ Checks the summary.xml file contained in thisDir
    for evidence that the directory is either a branch or the root.
    """
    from libSimControl import verifyFileExists, verifyDirExists, lockfile, unlockfile, dirIsRoot
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat
    
    for d in [thisDir]:
        verifyDirExists(d)
    if dirIsRoot(thisDir):
        return True
    lockname = lockfile(os.path.join(thisDir, 'xml', 'summary.xml'))
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    unlockfile(lockname)
    root = infoTree.getroot()
    nc = root.find('numberChildren')
    if nc is not None:
        # branches will have two children
        if nc.text == '2':
            return True
    return False

def isLeaf(thisDir):
    """ Checks the summary.xml file contained in thisDir
    for evidence that the directory is a leaf.
    """
    from libSimControl import verifyFileExists, verifyDirExists, lockfile, unlockfile, dirIsRoot
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat
    
    for d in [thisDir]:
        verifyDirExists(d)
    if dirIsRoot(thisDir):
        return True
    lockname = lockfile(os.path.join(thisDir, 'xml', 'summary.xml'))
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    unlockfile(lockname)
    root = infoTree.getroot()
    nc = root.find('numberChildren')
    if nc is not None:
        # leafs will have no children
        if nc.text == '0':
            return True
    return False

def nodeIsLeaf(nt):
    """Returns True if the newick tree supplied has 0 distance
    and is not a node.
    Input can either be a newickTree object or string.
    """
    from sonLib.bioio import newickTreeParser
    if isinstance(nt, str):
        nt = newickTreeParser(nt, 0.0)
    return not nt.internal and nt.distance == 0.0
    
def treeStr2Dir(treeStr, simDir):
    """ Takes a newick tree string and an options object and
    returns the directory path for this pair.
    """
    from libSimControl import nameTree
    import os
    from sonLib.bioio import newickTreeParser
    if not isinstance(treeStr, str):
        raise TypeError('treeStr should be a string, is %s\n' % treeStr.__class__)
    if not isinstance(simDir, str):
        raise TypeError('simDir should be a string, is %s\n' % treeStr.__class__)
    return os.path.abspath(os.path.join(simDir, 
                                        nameTree(newickTreeParser(treeStr, 0.0))
                                        ))

def lastOneOutTurnOffTheLightsCycle(thisDir):
    """ lastOneOutTurnOffTheLightsCycle() checks the .xml files in thisDir
    and if (1) both Stats and Transalign are finished and (2) it has not already
    been done, it adds an attribute (endEpochUTC) to the timestamps tag in the 
    xml file.
    This is called by both StatsStep4 and TransalignStep2
    """
    import os
    if statsAndTransAreComplete(thisDir):
        # this is the end of the cycle
        addEndTimeAttribute(os.path.join(thisDir, 'xml', 'summary.xml'))

def cycleIsComplete(thisDir):
    from libSimControl import lockfile, unlockfile, dirIsRoot
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat

    if not os.path.exists(thisDir):
        return False
    lockname = lockfile(os.path.join(thisDir, 'xml', 'summary.xml'))
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    unlockfile(lockname)
    root = infoTree.getroot()
    stamps = root.find('timestamps')
    if stamps is not None:
        return 'endEpochUTC' in stamps.attrib

def statsAndTransAreComplete(thisDir):
    """ Is the cycle done with both Stats and Transalign?
    """
    from libSimControl import lockfile, unlockfile, dirIsRoot
    import os
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat

    if not os.path.exists(thisDir):
        return False
    endings = { 'transalign': False, 'stats': False }
    for f in endings:
        lockname = lockfile(os.path.join(thisDir, 'xml', f+'.xml'))
        try:
            infoTree = ET.parse(lockname)
        except expat.ExpatError: # broken xml file
            sys.stderr.write('Bad xml: %s\n' % lockname)
            raise
        unlockfile(lockname)
        root = infoTree.getroot()
        stamps = root.find('timestamps')
        if stamps is not None:
            s = stamps.find('%sEnd' % f)
            if s:
                endings[f] = True
    
    return endings['transalign'] and endings['stats']

def addEndTimeAttribute(filename):
    """ opens the filename xml file, finds the timestamps tag, and adds an
    attribute for endEpochUTC. Used at the end of the Stats, Transalign and Cycle
    jobs.
    """
    from libSimControl import verifyFileExists, lockfile, unlockfile, dirIsRoot
    import os
    import time
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat
    verifyFileExists(filename)
    lockname = lockfile(filename)
    try:
        infoTree = ET.parse(lockname)
    except expat.ExpatError: # broken xml file
        sys.stderr.write('Bad xml: %s\n' % lockname)
        raise
    root = infoTree.getroot()
    t = root.find('timestamps')
    if t is None:
        raise RuntimeError('%s does not contain "timestamps" tag' % filename)
    now = time.time()
    t.attrib['endEpochUTC'] = str(now)
    t.attrib['elapsedTime'] = str(now - float(t.attrib['startEpochUTC']))
    info = ET.ElementTree(root)
    info.write(lockname)
    unlockfile(lockname)

def lastOneOutTurnOffTheLightsSimulation(thisDir, options):
    """ lastOneOutTurnOffTheLightsSimulation() checks the simulationInfo.xml file in thisDir
    and if (1) all of the leaf cycles for this sim are complete and (2) it has not already
    been recorded, it changes the simulationInfo.xml to record the end of the simulation.
    This is called by 
    """
    from libSimControl import verifyFileExists, verifyDirExists, lockfile, unlockfile, dirIsRoot
    import os
    import time
    import xml.etree.ElementTree as ET
    import xml.parsers.expat as expat
    for d in [thisDir]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'simulationInfo.xml')]:
        verifyFileExists(f)
    if allLeafsComplete(options):
        lockname = lockfile(os.path.join(thisDir, 'simulationInfo.xml'))
        try:
            infoTree = ET.parse(lockname)
        except expat.ExpatError: # broken xml file
            sys.stderr.write('Bad xml: %s\n' % lockname)
            raise
        root = infoTree.getroot()
        timeTag = root.find('timestamps')
        if timeTag is None:
            raise RuntimeError('%s does not contain "timestamps" tag' % 
                               os.path.join(thisDir, 'simulationInfo.xml'))
        now = time.time()
        timeTag.attrib['endEpochUTC'] = str(now)
        timeTag.attrib['elapsedTime'] = str(now - float(timeTag.attrib['startEpochUTC']))
        endTag = timeTag.find('end')
        if endTag is not None:
            return
        timeEnd        = ET.SubElement(timeTag, 'end')
        timeLocal      = ET.SubElement(timeEnd, 'humanLocal')
        timeLocal.text = str(time.strftime("%a, %d %b %Y %H:%M:%S (%Z) ", time.localtime()))
        timeHuman      = ET.SubElement(timeEnd, 'humanUTC')
        timeHuman.text = str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC) ", time.gmtime()))
        timeEpoch      = ET.SubElement(timeEnd, 'epochUTC')
        timeEpoch.text = str(time.time())
        info=ET.ElementTree(root)
        info.write(lockname)
        unlockfile(lockname)

def allLeafsComplete(options):
    """ allLeafsComplete checks all of the leaf directories for a simulation and looks in their
    xml/summary.xml files to see if they are finished running.
    """
    import os
    from sonLib.bioio import newickTreeParser
    leafs = extractLeafNames(newickTreeParser(options.inputNewick, 0.0))
    for l in leafs:
        if not cycleIsComplete(os.path.join(options.simDir, l)):
            return False
    return True

def extractLeafNames(nt):                                   
    """ given a newickTree BinaryTree object, returns a list
    of leaf names as given by their .ID values              
    """                                                     
    names = []
    if nt.right is not None:
        names += extractLeafNames(nt.right)
    if nt.left is not None:
        names += extractLeafNames(nt.left)
    if not nt.internal:
        names.append(nt.iD)
    return names

def transalignStep1Cmds_1(thisDir, thisParentDir, options):
    """ Produces a list of commands to run a stats step,
    called by StatsStep4
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for d in [thisDir, thisParentDir]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'seq.rev'), os.path.join(options.rootDir, 'seq.rev')]:
        verifyFileExists(f)
    
    DRAW_REV_BLOCK_SIZE=10000
    DRAW_REV_NT_PER_PIX=100000
    
    pipes = []
    cmds  = []
    
    outname = os.path.join(thisDir, 'inter-intra.aln.rev')
    if not os.path.exists(outname):
        cmd = [which('evolver_transalign')]
        cmd.append('-in1')
        cmd.append(os.path.join(thisDir, 'inter', 'inter.aln.rev'))
        cmd.append('-in2')
        cmd.append(os.path.join(thisDir, 'intra', 'intra.aln.rev'))
        cmd.append('-out')
        cmd.append(outname + '.tmp')
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'logs', 'transalign1.log'))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'aln.rev')
    if not os.path.exists(outname):
        # we have this check because if this is a restarted job
        # then we needn't recalculate these
        if isBranchOrRoot(thisParentDir):
            # In these cases the alignment above the branch point should not be carried
            # into the the descendant genomes. Alignments should only go back to the most
            # recent branch point.
            cmd = [which('ln')]
            cmd.append('-s')
            cmd.append(os.path.join(thisDir, 'inter-intra.aln.rev'))
            cmd.append(outname)
            pipes.append(None)
            cmds.append(cmd)
        else:
            cmd = [which('evolver_transalign')]
            cmd.append('-in1')
            cmd.append(os.path.join(thisParentDir, 'aln.rev'))
            cmd.append('-in2')
            cmd.append(os.path.join(thisDir, 'inter-intra.aln.rev'))
            cmd.append('-out')
            cmd.append(outname + '.tmp')
            cmd.append('-log')
            cmd.append(os.path.join(thisDir, 'logs', 'transalign2.log'))
            pipes.append(None)
            cmds.append(cmd)
            cmd = [which('mv')]
            cmd.append(outname + '.tmp')
            cmd.append(outname)
            pipes.append(None)
            cmds.append(cmd)

    outname = os.path.join(thisDir, 'stats', 'cds_aln.cycle.rev')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-cdsalns')
        cmd.append(os.path.join(thisDir, 'inter', 'inter.aln.rev'))
        cmd.append('-alns')
        cmd.append(outname + '.tmp')
        cmd.append('-annots1')
        cmd.append(os.path.join(thisDir, 'inter', 'inter.outannots.gff'))
        cmd.append('-annots2')
        cmd.append(os.path.join(thisDir, 'annots.gff'))
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'logs', 'cds_alns.cycle.log'))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)

    outname = os.path.join(thisDir, 'stats', 'img.cycle.cmap.pdf')
    if not os.path.exists(outname):
        cmd = [which('evolver_drawrev')]
        cmd.append('-fromrev')
        cmd.append(os.path.join(thisDir, 'inter', 'inter.aln.rev'))
        cmd.append('-tocmap')
        cmd.append(outname + '.tmp')
        cmd.append('-blocksize')
        cmd.append(str(DRAW_REV_BLOCK_SIZE))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)

    outname = os.path.join(thisDir, 'stats', 'img.cycle.lmap.pdf')
    if not os.path.exists(outname):
        cmd = [which('evolver_drawrev')]
        cmd.append('-fromrev')
        cmd.append(os.path.join(thisDir, 'inter', 'inter.aln.rev'))
        cmd.append('-tolmap')
        cmd.append(outname + '.tmp')
        cmd.append('-npp')
        cmd.append(str(DRAW_REV_NT_PER_PIX))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)

    outname = os.path.join(thisDir, 'stats', 'tmpstats.branch.difflength.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-nologcmdlineandtime')
        cmd.append('-ancstats')
        cmd.append(os.path.join(thisDir, 'aln.rev'))
        cmd.append('-log')
        cmd.append(outname + '.tmp')
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'tmpstats.cycle.difflength.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-nologcmdlineandtime')
        cmd.append('-ancstats')
        cmd.append(os.path.join(thisDir, 'intra', 'intra.aln.rev'))
        cmd.append('-log')
        cmd.append(outname + '.tmp')
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'img.branch.cmap.pdf')
    if not os.path.exists(outname):
        cmd = [which('evolver_drawrev')]
        cmd.append('-fromrev')
        cmd.append(os.path.join(thisDir, 'aln.rev'))
        cmd.append('-tocmap')
        cmd.append(outname + '.tmp')
        cmd.append('-blocksize')
        cmd.append(str(DRAW_REV_BLOCK_SIZE))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'img.branch.lmap.png')
    if not os.path.exists(outname):
        cmd = [which('evolver_drawrev')]
        cmd.append('-fromrev')
        cmd.append(os.path.join(thisDir, 'aln.rev'))
        cmd.append('-tolmap')
        cmd.append(outname + '.tmp')
        cmd.append('-npp')
        cmd.append(str(DRAW_REV_NT_PER_PIX))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'cds_aln.root.rev')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-cdsalns')
        cmd.append(os.path.join(thisDir, 'aln.rev'))
        cmd.append('-alns')
        cmd.append(outname + '.tmp')
        cmd.append('-annots1')
        cmd.append(os.path.join(options.rootDir, 'stats', 'cds_annots.gff'))
        cmd.append('-annots2')
        cmd.append(os.path.join(thisDir, 'annots.gff'))
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'logs', 'cds_alns.root.log'))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    outname = os.path.join(thisDir, 'stats', 'codonSubs.cycle.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-getcodonsubs')
        cmd.append(os.path.join(thisDir, 'stats', 'cds_aln.cycle.rev'))
        cmd.append('-out')
        cmd.append(outname + '.tmp')
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'logs', 'getCodonSubs.cycle.log'))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
        
    outname = os.path.join(thisDir, 'stats', 'codonSubs.root.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_evo')]
        cmd.append('-getcodonsubs')
        cmd.append(os.path.join(thisDir, 'stats', 'cds_aln.root.rev'))
        cmd.append('-out')
        cmd.append(outname + '.tmp')
        cmd.append('-log')
        cmd.append(os.path.join(thisDir, 'logs', 'getCodonSubs.root.log'))
        pipes.append(None)
        cmds.append(cmd)
        cmd = [which('mv')]
        cmd.append(outname + '.tmp')
        cmd.append(outname)
        pipes.append(None)
        cmds.append(cmd)
    
    return cmds, pipes

def runTransalignStep1Cmds_2(thisDir, thisParentDir, localTempDir, options):
    """ Produces a list of commands to run a stats step,
    called by StatsStep4
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    import subprocess
    for d in [thisDir, thisParentDir]:
        verifyDirExists(d)
    for f in [os.path.join(thisDir, 'seq.rev'), os.path.join(options.rootDir, 'seq.rev')]:
        verifyFileExists(f)
    
    inPipes  = []
    outPipes = []
    cmds     = []
    outnames = []

    outname = os.path.join(thisDir, 'stats', 'protStats.cycle.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_codon_report.pl')]
        cmd.append(os.path.basename(thisDir))
        cmd.append(os.path.basename(thisParentDir))
        cmds.append(cmd)
        inPipes.append(os.path.join(thisDir, 'stats', 'codonSubs.cycle.txt'))
        outPipes.append(outname + '.tmp')
        outnames.append(outname)
    
    outname = os.path.join(thisDir, 'stats', 'protStats.root.txt')
    if not os.path.exists(outname):
        cmd = [which('evolver_codon_report.pl')]
        cmd.append(os.path.basename(options.rootDir))
        cmd.append(os.path.basename(thisParentDir))
        cmds.append(cmd)
        inPipes.append(os.path.join(thisDir, 'stats', 'codonSubs.root.txt'))
        outPipes.append(outname + '.tmp')
        outnames.append(outname)
    
    outname = os.path.join(thisDir, 'stats', 'diffs.cycle.txt')
    if not os.path.exists(outname):
        cmd = [which('cat')]
        verifyFileExists(os.path.join(thisDir, 'stats', 'tmpstats.cycle.difflength.txt'))
        verifyFileExists(os.path.join(thisDir, 'stats', 'tmpstats.cycle.diffcompost.txt'))
        verifyFileExists(os.path.join(thisDir, 'stats', 'tmpstats.cycle.diffannots.txt'))
        cmd.append(os.path.join(thisDir, 'stats', 'tmpstats.cycle.difflength.txt'))
        cmd.append(os.path.join(thisDir, 'stats', 'tmpstats.cycle.diffcompost.txt'))
        cmd.append(os.path.join(thisDir, 'stats', 'tmpstats.cycle.diffannots.txt'))
        cmds.append(cmd)
        inPipes.append(None)
        outPipes.append(outname + '.tmp')
        outnames.append(outname)
    
    outname = os.path.join(thisDir, 'stats', 'diffs.branch.txt')
    if not os.path.exists(outname):
        cmd = [which('cat')]
        verifyFileExists(os.path.join(thisDir, 'stats', 'tmpstats.branch.difflength.txt'))
        verifyFileExists(os.path.join(thisDir, 'stats', 'tmpstats.branch.diffcompost.txt'))
        verifyFileExists(os.path.join(thisDir, 'stats', 'tmpstats.branch.diffannots.txt'))
        cmd.append(os.path.join(thisDir, 'stats', 'tmpstats.branch.difflength.txt'))
        cmd.append(os.path.join(thisDir, 'stats', 'tmpstats.branch.diffcompost.txt'))
        cmd.append(os.path.join(thisDir, 'stats', 'tmpstats.branch.diffannots.txt'))
        cmds.append(cmd)
        inPipes.append(None)
        outPipes.append(outname + '.tmp')
        outnames.append(outname)
    
    i = -1
    for c in cmds:
        i += 1
        if inPipes[i] is None:
            p = subprocess.Popen(c, cwd=localTempDir, stdout=subprocess.PIPE)
            f = open(outPipes[i], 'w')
            f.write(p.communicate()[0])
            f.close()
            handleReturnCode(p.returncode, c)
        else:
            p = subprocess.Popen(c, cwd=localTempDir, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
            f = open(outPipes[i], 'w')
            f.write(p.communicate(open(inPipes[i]).read())[0])
            f.close()
            handleReturnCode(p.returncode, c)
    i = -1
    for out in outnames:
        i += 1
        c = [which('mv'), outnames[i] + '.tmp', outnames[i]]
        p = subprocess.Popen(c)
        p.wait()
        handleReturnCode(p.returncode, c)
    
def extractLeafsFromNewick(nt, leafDict):
    """Given a newick tree object, it returns a dict of
    leaf objects. Operates recursively.
    """
    if nt is None:
        return None
    nt.distance = 0.0
    if nt.right is None and nt.left is None:
        leafDict[nt.iD] = True
    else:
        extractLeafsFromNewick(nt.right, leafDict = leafDict)
        extractLeafsFromNewick(nt.left , leafDict = leafDict)

def buildNodesListFromNewick(nt, nodesList, leafs, parentNode=''):
    """buildNodesListFromNewick() takes in a newick tree and the
    final output list and it creates a new Node object
    with appropriate information given the status of
    the newick passed in. Recursively calls itself to
    build the output list.
    """
    from libSimControl import nameTree
    from libSimControlClasses import Node
    if nt is None:
        return None
    nt.distance = 0
    n = Node()
    n.parent = parentNode
    n.name = nameTree(nt)
    leftName  = buildNodesListFromNewick(nt.right, nodesList = nodesList, 
                                         leafs = leafs, parentNode = n.name)
    rightName = buildNodesListFromNewick(nt.left,  nodesList = nodesList, 
                                         leafs = leafs, parentNode = n.name)
    if leftName is not None and rightName is not None:
        n.children.append(leftName)
        n.children.append(rightName)
        nodesList.append(n)
    if n.name in leafs:
        n.isLeaf = True
    if parentNode is not None:
        return n.name

def buildNodesDict(nl, leafsDict):
    """ takes a list of nodes and a dict of leafs
    and builds a dict containing the union of the two groups
    """
    from libSimControl import nameTree
    from libSimControlClasses import Node
    nd = {}
    parent = {}
    for n in nl:
        nd[n.name] = n
        for c in n.children:
            parent[c] = n.name
    for l in leafsDict:
        n = Node()
        n.name = l
        n.parent = parent[l]
        n.isLeaf = True
        nd[n.name] = n
    return nd

def buildNodeParentDict(nl):
    """ Takes a node list and builds a dict
    of parents with child names as keys and parent
    names as values
    """
    npd = {}
    for n in nl:
        for c in n.children:
            npd[c] = n.name
    return npd

def buildMergeCommand(maf1, maf2, out, treelessRootCmd, name, options, drop=None):
    """
    """
    from libSimControl import which, verifyDirExists, verifyFileExists
    import os
    for m in [maf1, maf2]:
        verifyFileExists(m)
    cmd = [which('mafJoin')]
    if len(treelessRootCmd) > 0:
        cmd += treelessRootCmd
    cmd.append(name)
    cmd.append('-maxBlkWidth=%d' % options.maxBlkWidth)
    # ' -maxBlkWidth=500' # avg120 mammal mouse-rat mouse-rat-human-tmp merge 
    # ' -maxBlkWidth=1000' avg120 mammal cow-dog merge.  
    # ' -maxBlkWidth=10000' default
    cmd.append('-maxInputBlkWidth=%d' % options.maxInputBlkWidth)
    # ' -maxInputBlkWidth=500'  avg120 mammal mouse-rat mouse-rat-human-tmp merge 
    # ' -maxInputBlkWidth=1000'
    if drop:
        cmd.append('-multiParentDropped=%s' % drop)
    cmd.append(maf1)
    cmd.append(maf2)
    cmd.append(out + '.tmp')
    cmds = [cmd]
    cmds.append([which('mv'), out + '.tmp', out])
    return cmds

def burninRootName(options):
    """returns a str that contains the name of the burnin's root genome
    """
    import os
    import re
    if not os.path.exists(os.path.join(options.rootDir, 'burnin.tmp.maf')):
        return ''
    pat = re.compile('^s (.*?)\.chr')
    f = open(os.path.join(options.rootDir, 'burnin.tmp.maf'))
    names = {}
    for line in f:
        line = line.strip()
        r = re.match(pat, line)
        if r:
            if r.group(1) not in names:
                names[r.group(1)] = True
    if len(names) == 2:
        for n in names:
            if n != options.rootName:
                return n
    return ''
