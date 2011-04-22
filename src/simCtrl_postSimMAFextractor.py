#!/usr/bin/env python
"""
postSimMAFextractor.py
dent earl, dearl (a) soe ucsc edu
15 april 2010
a script to be run following the completion of a
simulation. The script will check the simulationInfo.xml
file and then figure out the order to both
extract the relvant node's MAFs but then combine
the MAFs into larger alignments in the correct order.
"""
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
import os
import re
import sys
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = ['evolver_cvt', 'evolver_transalign',
            'simCtrl_commandEval.py', 'mafJoin']
LSC.verifyPrograms(programs)
(CVT_BIN, TRANS_BIN, CMD_EVAL_BIN,
 MAF_MERGE_BIN ) = programs

def usage():
    print 'USAGE: '+sys.argv[0]+' --simDir <dir> --jobFile JOB_FILE'
    print __doc__
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-p', '--simDir',dest='simDir',
                      help='Simulation directory.')
    parser.add_option('-j', '--jobFile',dest='jobFile',
                      help='jobFile, passed in by jobTree.')
    parser.add_option('-v', '--verbose',action='store_true', dest='isVerbose',
                      default=False,
                      help='if --debug is on, prints out more verbose debug statements.')
    parser.add_option('-m', '--mergeStep',action='store_true', dest='isMergeStep',
                      default=False,
                      help='the .aln.rev and .maf files have been created, now merge them.')
    parser.add_option('-d', '--debug', action='store_true', dest='isDebug',
                      default=False, help='Turns on debug output, does not issue jobs')
    parser.add_option('--noParalogy', action='store_true', dest='noParalogy',
                      default=False, help='adds a flag -noparalogy to the transalign call, switches paralogous blocks off. ')

class Node:
    """Nodes have one parent and two children,
    unless they are children in which case their
    children list is empty, or they are the root
    node in which case their parent is None.
    """
    def __init__(self):
        self.parent   = None
        self.name     = ''
        self.children = []
        self.isLeaf   = False

def findNode( nList, nName ):
    ans = None
    for n in nList:
        if n.name == nName:
            ans = n
            continue
    return ans

def checkOptions(options):
    if (options.simDir == None):
        sys.stderr.write('%s: Error, specify simulation dir.\n' % sys.argv[0])
        usage()
    if (options.jobFile == None) and (options.isDebug == False):
        sys.stderr.write('%s: Error, specify jobFile.\n' % sys.argv[0])
        usage()
    if (options.isDebug == False) and (not os.path.exists(options.jobFile)):
        sys.stderr.write('%s: Error, jobFile does not exist.\n' % sys.argv[0])
        usage()
    options.simDir  = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to find simulationInfo.xml.\n' % sys.argv[0])
        usage()
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    rootNameObj = infoTree.find('rootDir')
    options.rootName = os.path.basename( rootNameObj.text )
    options.rootDir    = os.path.abspath( rootNameObj.text )

def extractLeaves(nt, leafDict):
    """Given a newick tree object, it returns a dict of
    leaf objects. Operates recursively.
    """
    if nt == None:
        return None
    nt.distance=0
    if (nt.right == None) and (nt.left == None):
        leafDict[nt.iD] = 1
    else:
        extractLeaves(nt.right, leafDict=leafDict)
        extractLeaves(nt.left , leafDict=leafDict)

def buildNodesList(nt, nodesList, leaves, parentNode=''):
    """buildNodesList() takes in a newick tree and the
    final output list and it creates a new Node object
    with appropriate information given the status of
    the newick passed in. Recursively calls itself to
    build the output list.
    """
    if nt == None:
        return None
    nt.distance = 0
    n=Node()
    n.parent=parentNode
    n.name=LST.nameTree(nt)
    leftName  = (buildNodesList(nt.right, nodesList=nodesList, leaves=leaves, parentNode=n.name))
    rightName = (buildNodesList(nt.left,  nodesList=nodesList, leaves=leaves, parentNode=n.name))
    if leftName and rightName:
        n.children.append( leftName )
        n.children.append( rightName )
        nodesList.append( n )
    if n.name in leaves:
        n.isLeaf = True
    if parentNode:
        return n.name
    
def printNodesList(nl, leaves):
    if nl == None:
        return
    for n in nl:
        print 'Name: %s' %(n.name)
        for c in n.children:
            if c in leaves:
                print '  child %s (leaf)' %(c)
            else:
                print '  child %s' %(c)

def printCyclesDict( nd ):
    if nd == None:
        return
    for n in nd:
        print 'Name: %s (leaf:%s)' %( n, str( nd[ n ].isLeaf ) )
        for c in nd[ n ].children:
            print '  child %s (leaf:%s)' %(c, str( nd[ c ].isLeaf ))

def buildMAFpairs( options, nodesList, leaves ):
    """buldMAFpairs() takes the nodesList and creates all of the
    proper MAF pairs using the evolver scripts. eg
    ((a, b)E,(c,d)F)root;
    would decompose into the:
    aE.maf
    bE.maf
    cF.maf
    dF.maf
    Eroot.maf
    Froot.maf
    After this step is complete, the mafs will be merged.
    """
    if not options.isDebug:
        xmlTree = ET.parse(options.jobFile)
        jobElm=xmlTree.getroot()
        childrenElm = xmlTree.find( 'children' )
    runningJobs = 0
    # if the rootDir has an aln file, then it is considered to have
    # been a burnin genome. So the aln in the rootDir is back to the
    # root of the entire simulation. We will generate that maf.
    if os.path.exists( os.path.join(options.rootDir, 'aln.rev')) and not os.path.exists( os.path.join(options.rootDir, 'burnin.tmp.maf') ):
        transCMD = CMD_EVAL_BIN+ ' JOB_FILE "'
        transCMD += LSC.commandPacker(CVT_BIN+\
                                      ' -fromrev '+os.path.join(options.rootDir, 'aln.rev')+\
                                      ' -tomaf '+os.path.join(options.rootDir, 'burnin.tmp.maf'))
        transCMD += '"'
        runningJobs = runningJobs + 1
        if not options.isDebug:
            newChild = ET.SubElement( childrenElm, 'child' )
            newChild.attrib['command'] = transCMD
        else:
            mCMD = transCMD.replace('/cluster/home/dearl/sonTrace/bin/simCtrl_commandEval.py JOB_FILE "', '')
            mCMD = mCMD.replace('&myApos;', '\'')
            mCMD = mCMD.replace('&myCMD;&myCMD;', '\n')
            mCMD = mCMD.replace('&myCMD;', '\n')
            mCMD = mCMD.replace('"','')
            sys.stderr.write('\n%s' %( mCMD ))
    for n in nodesList:
        for c in n.children:
            if c in leaves:
                ext = '.maf'
            else:
                ext = '.tmp.maf'
            transCMD = CMD_EVAL_BIN+ ' JOB_FILE "'
            if not os.path.exists(os.path.join(options.simDir, c, n.name + ext)):
                if options.isVerbose and options.isDebug:
                    print 'well, I can\'t seem to see %s' %(os.path.join(options.simDir, c, n.name+'.maf'))
                transCMD += LSC.commandPacker(CVT_BIN+\
                                             ' -fromrev '+os.path.join(options.simDir, c, 'aln.rev')+\
                                             ' -tomaf '+os.path.join(options.simDir, c, n.name+ext))
                runningJobs = runningJobs + 1
            if os.path.exists(os.path.join(options.simDir, c, 'aln.rev')) and \
                   not os.path.exists(os.path.join(options.simDir, c, n.name+ext)):
                if options.isVerbose:
                    print '\nwell, I can\'t seem to see %s\n' %( os.path.join(options.simDir, c, n.name+ext) )
                runningJobs = runningJobs + 1
            transCMD += '"'
            if options.isDebug:
                if transCMD != CMD_EVAL_BIN+ ' JOB_FILE ""':
                    tcmd = transCMD.replace('&myCMD;&myCMD;','\n\t')
                    tcmd = tcmd.replace('"', '')
                    tcmd = tcmd.replace('/cluster/home/dearl/sonTrace/bin/simCtrl_commandEval.py JOB_FILE ', '')
                    tcmd = tcmd.replace('&myCMD;','')
                    if options.isVerbose:
                        print 'components of transCMD is(are): \n\t%s' % tcmd
                    else:
                        print '%s\n' %tcmd
            else:
                if transCMD != CMD_EVAL_BIN+ ' JOB_FILE ""':
                    newChild = ET.SubElement(childrenElm, 'child')
                    newChild.attrib['command'] = transCMD
    if runningJobs:
        followUpCommand = sys.argv[0] +\
                          ' --simDir '+options.simDir+\
                          ' --jobFile JOB_FILE '
    else:
        followUpCommand = sys.argv[0] +\
                          ' --simDir '+options.simDir+\
                          ' --jobFile JOB_FILE '+\
                          ' --mergeStep '
    if options.noParalogy:
        followUpCommand += ' --noParalogy '
    if not options.isDebug:
        jobElm.attrib['command'] = followUpCommand
        xmlTree.write(options.jobFile)
    else:
        sys.stderr.write('\nI wish to followUpCommand %s\n' %(followUpCommand))

def compareCommand( maf1, maf2, out, missing=None ):
    """ if isRoot, then childName is actually the nodeParent name.
    """
    compareStr  = MAF_COMPARE_BIN
    compareStr += ' --mAFFile1=' + maf1
    compareStr += ' --mAFFile2=' + maf2
    compareStr += ' --outputFile=' + out
    compareStr += ' --sampleNumber 100000000'
    if missing:
        compareStr += ' --ultraVerbose'
        compareStr += ' > ' + missing
    return compareStr

def mergeCommand( maf1, maf2, out, treelessRootStr, name, drop=None):
    mergeStr  = MAF_MERGE_BIN
    mergeStr += treelessRootStr 
    mergeStr += " '" + name + "'"
    mergeStr += ' -maxBlkWidth=250' #' -maxBlkWidth=500' # avg120 mammal mouse-rat mouse-rat-human-tmp merge # ' -maxBlkWidth=1000' avg120 mammal cow-dog merge.  #' -maxBlkWidth=10000' default
    mergeStr += ' -maxBlkWidth=250' #' -maxInputBlkWidth=500' # avg120 mammal mouse-rat mouse-rat-human-tmp merge # ' -maxInputBlkWidth=1000'
    if drop:
        mergeStr += ' -multiParentDropped=' + drop
    mergeStr += ' ' + maf1
    mergeStr += ' ' + maf2
    mergeStr += ' ' + out
    return mergeStr

def performMAFmerge( options, nodesList, leaves, nodeParentDict, nodesDict ):
    """performMAFmerge() takes the nodesList and goes through the
    steps of progressively merging the MAFs, starting from the leaves
    and working its way up the tree.
    ((a, b)E,(c,d)F)root;
    which was earlier decomposed into the:
    aE.maf
    bE.maf
    cF.maf
    dF.maf
    Eroot.maf
    Froot.maf
    would now be combined in the following pattern:
    aE.maf + bE.maf -> abE.maf
    cF.maf + dF.maf -> cdF.maf
    abE.maf + Eroot.maf -> abEroot.maf
    cdF.maf + Froot.maf -> cdFroot.maf
    abEroot.maf + cdFroot.maf -> abcdEFroot.maf
    ... and then we drink root beers and have highfives.
    OR, in the specific directory based data structure
    that we actually use,
    a/E.maf
    b/E.maf
    c/F.maf
    d/F.maf
    E/root.tmp.maf
    F/root.tmp.maf
    a/E.maf b/E.maf -> E/E.maf
    E/E.maf E/root.tmp.maf -> E/root.maf
    F/F.maf F/root.tmp.maf -> F/root.maf
    E/root.maf F/root.maf -> root/all.maf
    OR, most generally,
    child0/parent.maf child1/parent.maf -> parent/parent.maf
    parent/parent.maf parent/grandParent.tmp.maf -> parent/grandParent.maf
    """
    if not options.isDebug:
        xmlTree = ET.parse(options.jobFile)
        jobElm=xmlTree.getroot()
        childrenElm = xmlTree.find('children')
    runningJobs=0
    for n in nodesList:
        if n.name != options.rootName:
            nodeParent = nodeParentDict[n.name]
        else:
            nodeParent = options.rootName
        if (not os.path.exists(os.path.join(options.simDir,n.name, nodeParent+'.maf')) and 
            os.path.exists(os.path.join(options.simDir, n.children[0], n.name+'.maf')) and 
            os.path.exists(os.path.join(options.simDir, n.children[1], n.name+'.maf'))):
            mergeCMD = CMD_EVAL_BIN + ' JOB_FILE "'

            treelessRootStr = ''
            if n.children[0] in leaves :
                treelessRootStr += ' -treelessRoot1='+ str( n.name )
            if n.children[1] in leaves :
                treelessRootStr += ' -treelessRoot2='+ str( n.name )
            ##############################
            # the 'lookdown' aspect of the merge is performed for every node, including the root.
            maf1 = os.path.join(options.simDir, n.children[0], n.name+'.maf')
            maf2 = os.path.join(options.simDir, n.children[1], n.name+'.maf')
            mergeOut  = os.path.join(options.simDir, n.name, n.name+'.maf')
            drop = os.path.join(options.simDir, n.name, n.name+'.dropped.tab')
            mergeStr = mergeCommand( maf1, maf2, mergeOut, treelessRootStr, str( n.name ), drop )
            mergeCMD += LSC.commandPacker( mergeStr )

            ##############################
            # The 'lookup' aspect of the merge is only performed when we are not at the root
            # This merge merges the results of the 'lookdown' merge, that is to say the maf that contains
            # all descendant sequences including the node, with the node-parent maf, to produce a maf
            # that the parent can use to merge its children.
            if n.name != options.rootName:
                treelessRootStr = ' -treelessRoot2='+ str( nodeParent )
                maf1 = os.path.join(options.simDir, n.name, n.name+'.maf')
                maf2 = os.path.join(options.simDir, n.name, nodeParent+'.tmp.maf')
                mergeOut  = os.path.join(options.simDir, n.name, nodeParent+'.maf')
                drop = os.path.join(options.simDir, n.name, nodeParent+'.dropped.tab')
                mergeStr = mergeCommand( maf1, maf2, mergeOut, treelessRootStr, str( n.name ), drop)
                mergeCMD += LSC.commandPacker( mergeStr )

            mergeCMD += '"'
            runningJobs = runningJobs + 1
            if not options.isDebug:
                newChild = ET.SubElement(childrenElm, 'child')
                newChild.attrib['command'] = mergeCMD
            else:
                mCMD = mergeCMD.replace('/cluster/home/dearl/sonTrace/bin/simCtrl_commandEval.py JOB_FILE "', '')
                mCMD = mCMD.replace('&myApos;', '\'')
                mCMD = mCMD.replace('&myCMD;&myCMD;', '\n')
                mCMD = mCMD.replace('&myCMD;', '\n')
                mCMD = mCMD.replace('"','')
                if options.isVerbose:
                    sys.stderr.write('n: %s\n  I wish to perform:%s' %( n.name, mCMD ))
                else:
                    sys.stderr.write('%s' % mCMD )
        else:
            if options.isDebug and options.isVerbose:
                print 'n: %s\n  Not running cycle for the following reasons:' %(n.name)
                print 'exists: '+os.path.join(options.simDir, n.name, n.name+'.maf') +': '+ str(os.path.exists(os.path.join(options.simDir, n.name, n.name+'.maf')))
                print 'exists: '+os.path.join(options.simDir, n.children[0], n.name+'.maf') +': '+ str(os.path.exists(os.path.join(options.simDir, n.children[0], n.name+'.maf')))
                print 'exists: '+os.path.join(options.simDir, n.children[1], n.name+'.maf') +': '+ str(os.path.exists(os.path.join(options.simDir, n.children[1], n.name+'.maf')))
    if runningJobs:
        followUpCommand = sys.argv[0] +\
                          ' --simDir '+options.simDir+\
                          ' --jobFile JOB_FILE '+\
                          ' --mergeStep '
        if options.noParalogy:
            followUpCommand += ' --noParalogy '
    else:
        if os.path.exists( os.path.join( options.rootDir, 'burnin.tmp.maf') ) and os.path.exists( os.path.join( options.rootDir, options.rootName+'.maf') ) and not os.path.exists( os.path.join( options.rootDir, 'burnin.maf')):
            mergeCMD = CMD_EVAL_BIN + ' JOB_FILE "'
            treelessRootStr = ' -treelessRoot2='+ burninRootName( options ) # the burnin's root name
            maf1 = os.path.join(options.rootDir, options.rootName+'.maf' )
            maf2 = os.path.join(options.rootDir, 'burnin.tmp.maf' )
            mergeOut  = os.path.join(options.rootDir, 'burnin.maf' )
            drop = os.path.join(options.rootDir, 'burnin.dropped.tab' )
            mergeStr = mergeCommand( maf1, maf2, mergeOut, treelessRootStr, options.rootName, drop )
            mergeCMD += LSC.commandPacker( mergeStr )
            mergeCMD += '"'
            runningJobs = runningJobs + 1
            followUpCommand = sys.argv[0] +\
                              ' --simDir '+options.simDir+\
                              ' --jobFile JOB_FILE '+\
                              ' --mergeStep '
            if not options.isDebug:
                newChild = ET.SubElement( childrenElm, 'child' )
                newChild.attrib['command'] = mergeCMD
                xmlTree.write( options.jobFile )
            
    if options.isDebug:
        if runningJobs:
            mCMD = mergeCMD.replace('/cluster/home/dearl/sonTrace/bin/simCtrl_commandEval.py JOB_FILE "', '')
            mCMD = mCMD.replace('&myApos;', '\'')
            mCMD = mCMD.replace('&myCMD;&myCMD;', '\n')
            mCMD = mCMD.replace('&myCMD;', '\n')
            mCMD = mCMD.replace('"','')
            sys.stderr.write('\nI wish to perform %s\n%s\n' %( mCMD, followUpCommand ))
        else:
            sys.stderr.write('I think I\'m done -- I have no running jobs\n')
    else:
        if runningJobs:
            jobElm.attrib['command'] = followUpCommand
            xmlTree.write(options.jobFile)

def burninRootName( options ):
    """returns a str that contains the name of the burnin's root genome
    """
    if not os.path.exists( os.path.join(options.rootDir, 'burnin.tmp.maf')):
        return ''
    pat = re.compile('^s (.*?)\.chr')
    f = open( os.path.join(options.rootDir, 'burnin.tmp.maf') )
    names = {}
    for line in f:
        line = line.strip()
        r = re.match( pat, line )
        if r:
            if r.group(1) not in names:
                names[ r.group(1) ] = True
    if  len( names ) == 2:
        for n in names:
            if n != options.rootName:
                return n
    return ''

def nodeParentDictBuilder(nl):
    npd = {}
    for n in nl:
        for c in n.children:
            npd[c] = n.name
    return npd

def buildCyclesDict( nl, leaves ):
    """ this will be useful to lookup a node's leaf status.
    """
    nd = {}
    parent = {}
    for n in nl:
        nd[ n.name ] = n
        for c in n.children:
            parent[ c ] = n.name
    for l in leaves:
        n = Node()
        n.name = l
        n.parent = parent[ l ]
        n.isLeaf = True
        nd[ n.name ] = n
    return nd

def main():
    parser=OptionParser()
    initOptions( parser )
    ( options, args ) = parser.parse_args()
    checkOptions( options )
    nt = newickTreeParser( options.inputNewick, 0.0 )
    if nt.iD == None:
        nt.iD = options.rootName
    #####
    # step one, discover all the nodes and children
    leaves={}
    extractLeaves( nt, leaves )
    nodesList=[]
    cyclesDict={} # keyed by name 
    buildNodesList( nt, nodesList, leaves )
    cyclesDict = buildCyclesDict( nodesList, leaves )
    if options.isDebug and options.isVerbose:
        printNodesList( nodesList, leaves )
    if not options.isMergeStep:
        # step two, create the commands to build all pairwise MAFs
        buildMAFpairs( options, nodesList, leaves )
    else:
        # step three, create the commands to progressively combine MAFs
        nodeParentDict = nodeParentDictBuilder( nodesList )
        performMAFmerge( options, nodesList, leaves, nodeParentDict, cyclesDict )


if __name__ == "__main__":
    main()
