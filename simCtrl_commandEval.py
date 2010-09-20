#!/usr/bin/env python
"""
simCtrl_commandEval.py, 17 February 2009
 dent earl, dearl (a) soe ucsc edu

commandEval.py takes a JOB_TREE file and a quoted
string of individual commands, each of which is
separated by &myCMD; .

COMMAND MACROS:
&myQuot; \"
&myApos; \'
LOCAL_DIR the location of the local temp directory
JOB_XML the jobFile

We use this program not only to ensure that commands are executed
in the local_temp_dir but also to check the return code of a string
of commands that are executed in serial.

For example:
simCtrl_commandEval.py \" &myCMD;command1&myCMD; &myCMD;command2&myCMD; &myCMD;command3&myCMD; ... \"
internal quotes have to be replaced with &myQuot; and apostrophes with &myApos; .

"""
########################################
import os
import resource
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from optparse import OptionParser
import simulation.lib.libSimControl as LSC

def usage():
    print "USAGE: %s JOB_FILE \" &myCMD;command1&myCMD; &myCMD;command2&myCMD; &myCMD;command3&myCMD; ... \"" %(sys.argv[0])
    print __doc__
    sys.exit(2)

def initOptions(parser):
    """initOptions() initializes the options that will go to the
    parser object
    """
    parser.add_option('-s', '--statXML',dest='statXML',
                      help='file to write run stats to.')

def recordStats(file, i, name , t0, t1, r):
    """recordStats() writes out details of the command run
    to an xml file.
    """
    if os.path.exists(file):
        xmlTree = ET.parse(file)
        infotree = xmlTree.getroot()
        if infotree.find('timestamps'):
            ts = infotree.find('timestamps')
            if ts.find('micro'):
                tm = ts.find('micro')
            else:
                tm = ET.SubElement(ts, 'micro')
        else:
            ts = ET.SubElement(infotree, 'timestamps')
            tm = ET.SubElement(ts, 'micro')
    else:
        infotree = ET.Element('info')
        ts = ET.SubElement(infotree, 'timestamps')
        tm = ET.SubElement(ts, 'micro')
    # Step one, either create the file with the correct tag hierarchy, or if
    # it already exists then just get us to the starting point.
    ##############################
    nObj = ET.SubElement(tm, 'command_'+str(i))
    nObj.text = str(name)
    oObj = ET.SubElement(nObj, 'start')
    timeHuman=ET.SubElement(oObj, 'humanUTC')
    timeHuman.text=str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime(t0)))
    timeEpoch=ET.SubElement(oObj, 'epochUTC')
    timeEpoch.text=str(t0)
    oObj = ET.SubElement(nObj, 'end')
    timeHuman=ET.SubElement(oObj, 'humanUTC')
    timeHuman.text=str(time.strftime("%a, %d %b %Y %H:%M:%S (UTC)", time.gmtime(t1)))
    timeEpoch=ET.SubElement(oObj, 'epochUTC')
    timeEpoch.text=str(t1)
    oObj = ET.SubElement(nObj, 'elapsed')
    oObj.text=str(t1-t0)
    oObj = ET.SubElement(nObj, 'elapsedPretty')
    oObj.text=prettyTime(t1-t0)
    oObj = ET.SubElement(nObj, 'resources')
    pObj = ET.SubElement(oObj, 'ru_utime')
    pObj.text = str(r[0])
    pObj = ET.SubElement(oObj, 'ru_stime')
    pObj.text = str(r[1])
#     pObj = ET.SubElement(oObj, 'ru_maxrss')
#     pObj.text = str(r[2])
#     pObj = ET.SubElement(oObj, 'ru_ixrss')
#     pObj.text = str(r[3])
    pObj = ET.SubElement(oObj, 'myWait')
    pObj.text = str((t1-t0)-(r[0]+r[1]))
    info=ET.ElementTree(infotree)
    info.write(file)

def prettyTime(t):
    """ prettyTime(t): input t comes in as seconds. Output is a str
    with units attached (secs, mins, hrs, etc.)
    """
    if(t < 210):
        return '%.2f secs' %(t)
    t = t/60.0
    if(t < 121):
        return '%.2f mins' %(t)
    t = t/60.0
    if(t < 25):
        return '%.2f hrs' %(t)
    t = t/24.0
    if(t < 28):
        return '%.2f days' %(t)
    t = t/7
    return '%.2f weeks' %(t)

def tupleDiff(a,b):
    c=[]
    for i in range(0,len(a)):
        c.append(a[i] - b[i])
    return c

def main(argv):
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    jobFile = args[0]
    command = args[1]
    if not os.path.exists(jobFile):
        sys.stderr.write('%s: Unable to find jobFile, %s.\n' %(sys.argv[0], jobFile))
        usage()
    xmlTree  = ET.parse(jobFile)
    jobElm   = xmlTree.getroot()
    localDir = jobElm.attrib['local_temp_dir']
    if not os.path.isdir(localDir):
        sys.stderr.write('%s: Experienced an error, localDir, (%s), is not a directory!\n' %(sys.argv[0], localDir))
        sys.exit(1)
    commands = LSC.commandUnPacker(command)
    i = -1
    for c in commands:
        d=c.strip()
        if not (c and d):
            continue
        i=i+1
        c=c.replace('LOCAL_DIR', localDir)
        c=c.replace('JOB_XML', jobFile)
        t0=time.time()
        r = resource.getrusage(resource.RUSAGE_CHILDREN)
        oldR = r[0:4]
        p = subprocess.Popen(c, cwd=localDir, shell=True)
        p.wait()
        t1=time.time()
        r = resource.getrusage(resource.RUSAGE_CHILDREN)
        newR = tupleDiff(r[0:4], oldR)
        oldR=newR
        if p.returncode:
            if p.returncode < 0:
                sys.stderr.write('%s: Experienced an error while trying to execute: %s SIGNAL:%d\n' %(sys.argv[0], c, -(p.returncode)))
            else:
                sys.stderr.write('%s: Experienced an error while trying to execute: %s retcode:%d\n' %(sys.argv[0], c, p.returncode))
            sys.exit(1)
        if options.statXML:
            recordStats(options.statXML, i, c, t0, t1, newR)
        

if __name__ == "__main__":
    main(sys.argv)
