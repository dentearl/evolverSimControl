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
##############################
import os, sys
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import simualtion.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

programs = ['evolver_cvt', 'evolver_transalign',
            'simCtrl_commandEval.py', 'head']
LSC.verifyPrograms(programs)
(CVT_BIN, TRANS_BIN, CMD_EVAL_BIN,
 MAF_MERGE_BIN) = programs

def usage():
    print 'USAGE: '+sys.argv[0]+' --simDir <dir> --jobFile JOB_FILE'
    print __doc__
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-p', '--simDir',dest='simDir',
                      help='Simulation directory.')
    parser.add_option('-j', '--jobFile',dest='jobFile',
                      help='jobFile, passed in by jobTree.py.')
    parser.add_option('-m', '--mergeStep',action='store_true', dest='isMergeStep',
                      default=False,
                      help='the .aln.rev and .maf files have been created, now merge them.')
    parser.add_option('-d', '--debug', action='store_true', dest='isDebug',
                      default=False, help='Turns on debug output, does not issue jobs')

class Node:
    """Nodes have one parent and two children,
    unless they are children in which case their
    children list is None, or they are the root
    node in which case their parent is None.
    """
    def __init__(self):
        self.parent   = None
        self.name     = ''
        self.children = []
        self.isLeaf   = False

def checkOptions(options):
    if (options.simDir == None):
        sys.stderr.write('%s: Error, specify simulation dir.\n' % sys.argv[0])
        usage()
    if (options.jobFile == None) and (options.isDebug == False):
        sys.stderr.write('%s: Error, specify jobFile.\n' % sys.argv[0])
        usage()
    if (not os.path.exists(options.jobFile)) and (options.isDebug == False):
        sys.stderr.write('%s: Error, jobFile does not exist.\n' % sys.argv[0])
        usage()
    options.simDir  = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to find simutaltionInfo.xml.\n' % sys.argv[0])
        usage()
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    treeObj = infoTree.find('rootDir')
    options.rootDir=treeObj.text

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
        n.children.append(leftName)
        n.children.append(rightName)
        nodesList.append(n)
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

def buildMAFpairs(options, nodesList, leaves):
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
        childrenElm = xmlTree.find('children')
    runningJobs = 0
    for n in nodesList:
        for c in n.children:
            if c in leaves:
                ext = '.maf'
            else:
                ext = '.tmp.maf'
            transCMD = CMD_EVAL_BIN+ ' JOB_FILE "'
            if not os.path.exists(os.path.join(options.simDir, c, n.name+'.aln.rev')):
                print 'well, I can\'t seem to see %s' %(os.path.join(options.simDir, c, n.name+'.aln.rev'))
                transCMD += LSC.commandPacker(TRANS_BIN+\
                                              ' -in1 '+ os.path.join(options.simDir,c, 'root.aln.rev') + \
                                              ' -in2 '+ os.path.join(options.simDir,n.name, 'root.aln.rev')  + \
                                              ' -out '+ os.path.join(options.simDir,c, n.name+'.aln.rev') + \
                                              ' -log '+ os.path.join(options.simDir,c, 'logs', 'transalign.'+n.name+'.log'))+\
                           LSC.commandPacker(CVT_BIN+\
                                             ' -fromrev '+os.path.join(options.simDir, c, n.name+'.aln.rev')+\
                                             ' -tomaf '+os.path.join(options.simDir, c, n.name+ext))
                runningJobs = runningJobs + 1
            if os.path.exists(os.path.join(options.simDir, c, n.name+'.aln.rev')) and \
                   not os.path.exists(os.path.join(options.simDir, c, n.name+ext)):
                print '\nwell, I can\'t seem to see %s or %s\n' %(os.path.join(options.simDir, c, n.name+'.aln.rev'), os.path.join(options.simDir, c, n.name+ext))
                transCMD += LSC.commandPacker(CVT_BIN+\
                                              ' -fromrev '+os.path.join(options.simDir, c, n.name+'.aln.rev')+\
                                              ' -tomaf '+os.path.join(options.simDir, c, n.name+ext))
                runningJobs = runningJobs + 1
            transCMD += '"'
            if options.isDebug:
                print 'transCMD is '+transCMD
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
    if not options.isDebug:
        jobElm.attrib['command'] = followUpCommand
        xmlTree.write(options.jobFile)
    else:
        sys.stderr.write('I wish to followUpCommand %s\n' %(followUpCommand))

def performMAFmerge(options, nodesList, leaves, nodeParentDict):
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
    """
    if not options.isDebug:
        xmlTree = ET.parse(options.jobFile)
        jobElm=xmlTree.getroot()
        childrenElm = xmlTree.find('children')
    runningJobs=0
    for n in nodesList:
        if n.name != 'root':
            nodeParent = nodeParentDict[n.name]
        else:
            nodeParent = 'root'
        if (not os.path.exists(os.path.join(options.simDir,n.name, nodeParent+'.maf')) and 
            os.path.exists(os.path.join(options.simDir, n.children[0], n.name+'.maf')) and 
            os.path.exists(os.path.join(options.simDir, n.children[1], n.name+'.maf'))):
            mergeCMD = CMD_EVAL_BIN + ' JOB_FILE "'+\
                       LSC.commandPacker(MAF_MERGE_BIN +\
                                         ' '+os.path.join(options.simDir, n.children[0], n.name+'.maf')+\
                                         ' '+os.path.join(options.simDir, n.children[1], n.name+'.maf')+\
                                         ' > ' + os.path.join(options.simDir, n.name, n.name+'.maf'))
            if n.name != 'root':
                mergeCMD +=LSC.commandPacker(MAF_MERGE_BIN +\
                                             ' '+os.path.join(options.simDir, n.name, nodeParent+'.tmp.maf')+\
                                             ' '+os.path.join(options.simDir, n.name, n.name+'.maf')+\
                                             ' > '+os.path.join(options.simDir, n.name, nodeParent+'.maf'))
            mergeCMD += '"'
            runningJobs = runningJobs + 1
            if not options.isDebug:
                newChild = ET.SubElement(childrenElm, 'child')
                newChild.attrib['command'] = mergeCMD
            else:
                sys.stderr.write('n: %s\n  I wish to perform %s\n' %(n.name, mergeCMD))
        else:
            if options.isDebug:
                print 'n: %s\n  Not running cycle for the following reasons:' %(n.name)
                print os.path.join(options.simDir, n.name, n.name+'.maf') +': '+ str(os.path.exists(os.path.join(options.simDir, n.name, n.name+'.maf')))
                print os.path.join(options.simDir, n.children[0], n.name+'.maf') +': '+ str(os.path.exists(os.path.join(options.simDir, n.children[0], n.name+'.maf')))
                print os.path.join(options.simDir, n.children[1], n.name+'.maf') +': '+ str(os.path.exists(os.path.join(options.simDir, n.children[1], n.name+'.maf')))
    if runningJobs:
        followUpCommand = sys.argv[0] +\
                          ' --simDir '+options.simDir+\
                          ' --jobFile JOB_FILE '+\
                          ' --mergeStep '
    if options.isDebug:
        if runningJobs:
            sys.stderr.write('I wish to perform %s\n' %(followUpCommand))
        else:
            sys.stderr.write('I seem to think I\'m done I have no running jobs\n')
    else:
        if runningJobs:
            jobElm.attrib['command'] = followUpCommand
            xmlTree.write(options.jobFile)

def nodeParentDictBuilder(nl):
    npd = {}
    for n in nl:
        for c in n.children:
            npd[c] = n.name
    return npd

def main():
    parser=OptionParser()
    initOptions(parser)
    (options, args) = parser.parse_args()
    checkOptions(options)
    nt = newickTreeParser(options.inputNewick, 0.0)
    if nt.iD == None:
        nt.iD= 'root'
    #####
    # step one, discover all the nodes and children
    leaves={}
    extractLeaves(nt, leaves)
    nodesList=[]
    buildNodesList(nt, nodesList, leaves)
    if options.isDebug:
        printNodesList(nodesList, leaves)
    if not options.isMergeStep:
        # step two, create the command to build all MAFs
        buildMAFpairs(options, nodesList, leaves)
    else:
        # step three, create the command to progressively combine MAFs
        nodeParentDict = nodeParentDictBuilder(nodesList)
        performMAFmerge(options, nodesList, leaves, nodeParentDict)

    ####
    # finished
    if options.isDebug:
        return

if __name__ == "__main__":
    main()
