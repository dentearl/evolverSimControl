#!/usr/bin/env python
# simTreeSingle.py
# dent earl, dearl (a) soe ucsc edu
# 22 oct 2009
# a recursive script to control an entire tree's simulation.
# Give the parent, and a newick tree as a command line
# option and simTree launches an evolver cycle via the bash script
# cycle.sh
#
# "((Steve:.003, Zack:.005):.002);"
# "(Steve:0.003, (Zack:.005, Chris:.02):.001):.002;"
# "(  (Chimp: 0.052,   Human: 0.042): 0.007,  Gorilla: 0.060,  (Gibbon: 0.124,   Orangutan: 0.0971): 0.038);"
#
##############################
from optparse import OptionParser
#import xml.etree.ElementTree as ET
import copy, glob, os, re, shutil, subprocess, sys

CYCLE_SH = '/cluster/home/dearl/sonTrace/src/eval/cycleControl/cycleStandAlone.sh'
if(not os.path.exists(CYCLE_SH)):
    sys.stderr.write('ERROR: Could not locate "cycleStandAlone.sh" script at expected location: %s. Please change source code.\n' %(CYCLE_SH))
    sys.exit(1)

#########################################################
# BEGINNING OF BENEDICT PATEN'S BINARY TREE CODE.
#////////////////////////////////////////
#////////////////////////////////////////
DEFAULT_DISTANCE = 0.001
def newickTreeParser(newickTree, defaultDistance=DEFAULT_DISTANCE, \
                     sortNonBinaryNodes=False, reportUnaryNodes=False):
    """
    lax newick tree parser
    """
    newickTree = newickTree.replace("(", " ( ")
    newickTree = newickTree.replace(")", " ) ")
    newickTree = newickTree.replace(":", " : ")
    newickTree = newickTree.replace(";", "")
    newickTree = newickTree.replace(",", " , ")
    
    newickTree = re.compile("[\s]*").split(newickTree)
    while "" in newickTree:
        newickTree.remove("")
    def fn(newickTree, i):
        if i[0] < len(newickTree):
            if newickTree[i[0]] == ':':
                d = float(newickTree[i[0]+1])
                i[0] += 2
                return d
        return defaultDistance
    def fn2(newickTree, i):
        if i[0] < len(newickTree):
            j = newickTree[i[0]]
            if j != ':' and j != ')' and j != ',':
                i[0] += 1
                return j
        return None
    def fn3(newickTree, i):
        if newickTree[i[0]] == '(':
            #subTree1 = None
            subTreeList = []
            i[0] += 1
            k = []
            while newickTree[i[0]] != ')':
                if newickTree[i[0]] == ',':
                    i[0] += 1
                subTreeList.append(fn3(newickTree, i))
            i[0] += 1
            def cmp(i, j):
                if i.distance < j.distance:
                    return -1
                if i.distance > j.distance:
                    return 1
                return 0
            if sortNonBinaryNodes:
                subTreeList.sort(cmp)
            subTree1 = subTreeList[0]
            if len(subTreeList) > 1:
                for subTree2 in subTreeList[1:]:
                    subTree1 = BinaryTree(0.0, True, subTree1, subTree2, None)
                subTree1.iD = fn2(newickTree, i)
                subTree1.distance += fn(newickTree, i)
            elif reportUnaryNodes:
                subTree1 = BinaryTree(0.0, True, subTree1, None, None)
                subTree1.iD = fn2(newickTree, i)
                subTree1.distance += fn(newickTree, i)
            else:
                fn2(newickTree, i)
                subTree1.distance += fn(newickTree, i)
            return subTree1
        leafID = fn2(newickTree, i)
        return BinaryTree(fn(newickTree, i), False, None, None, leafID)
    return fn3(newickTree, [0])
def printBinaryTree(binaryTree, includeDistances, dontStopAtID=True, distancePrintFn=(lambda f : "%f" % f)):
    def fn(binaryTree):
        #print " tree Node ", binaryTree.left, binaryTree.right, binaryTree.distance, binaryTree.internal, binaryTree.iD 
        if binaryTree.iD != None:
            iD = str(binaryTree.iD)
        else:
            iD = ''
        if binaryTree.internal and (dontStopAtID or binaryTree.iD == None):
            if binaryTree.right != None:
                s = '(' + fn(binaryTree.left) + ',' + fn(binaryTree.right) + ')' + iD
            else:
                s = '(' + fn(binaryTree.left) + ')' + iD 
        else:
            s = iD
        if includeDistances:
            return s + ':' + distancePrintFn(binaryTree.distance)
        return s
    return fn(binaryTree) + ';'
class BinaryTree:
    def __init__(self, distance, internal, left, right, iD):
        self.distance = distance
        self.internal = internal
        self.left = left
        self.right = right
        self.iD = iD
class TraversalID:
    """
    tree traversal numbers, used as nodeIDs for identifying
    orders in the tree
    """
    def __init__(self, midStart, mid, midEnd):
        self.midStart = midStart
        self.mid = mid
        self.midEnd = midEnd
#////////////////////////////////////////
#////////////////////////////////////////
# END OF BENEDICT PATEN'S BINARY TREE CODE.
##################################################
def usage():
    print 'USAGE: '+sys.argv[0]+' --parent <dir> --out [optional dir] --tree <newick tree in quotes> --params <parameter dir>'
    sys.exit(2)
def fixName(name):
    name=name.replace(' ','')
    name=name.replace(',','')
    name=name.replace(':','-')
    name=name.replace('.','_')
    name=name.replace(';','')
    name=name.replace('\'','')
    name=name.replace('"','')
    name=name.replace('(','_L_')
    name=name.replace(')','_R_')
    return name
def branchLog(message):
    curr = os.curdir
    logPath = os.path.join(curr, 'branch_log.log')
    if(not os.path.exists(logPath)):
        FILE = open(logPath, 'w')
        FILE.write( '%s' % (message))
        FILE.close()
    else:
        FILE = open(logPath, 'a')
        FILE.write( '%s' % (message))
        FILE.close()
def main():
#    sys.stderr.write('\n%s\n\n' %(sys.argv))
    stepSize =0.001
    parser=OptionParser()
    parser.add_option('-p', '--parent',dest='parentDir',
                      help='Parent directory.')
    parser.add_option('-o', '--out',dest='outDir',
                      help='Out directory.')
    parser.add_option('-s', '--step',dest='stepSize', action="store",
                      type ='float', help='stepSize for each cycle.')
    parser.add_option('-g', '--params',dest='gParamsDir',
                      help='Parameter directory.')
    parser.add_option('-t', '--tree',dest='inputNewick',
                      help='Newick tree.')
    parser.add_option('-y', '--saveParent', action='store_true', dest='saveParent',
                      default=False, help='Save all cycles.')
    parser.add_option('-c', '--isContinue', action='store_true', dest='isContinue',
                      default=False, help='Is this a recursive run?')
    parser.add_option('-l', '--logBranch', action='store_true', dest='logBranch',
                      default=False, help='record all tree travel information.')
    parser.add_option('-T', '--testTree', action='store_true', dest='testTree',
                      default=False, help='Instead of performing a simulation, does dry run with empty dirs.')
    parser.add_option('-v', '--verbose', action='store_true', dest='verbose',
                      default=False, help='Turn on verbose output.')
    (options, args) = parser.parse_args()

    if (options.inputNewick == None):
        sys.stderr.write('%s: Error, specify newick.\n' % sys.argv[0])
        usage()
    if (options.parentDir == None):
        sys.stderr.write('%s: Error, specify parent dir.\n' % sys.argv[0])
        usage()
    parentDir  = os.path.abspath(options.parentDir)
    if (options.gParamsDir == None):
        sys.stderr.write('%s: Error, specify params.\n' % sys.argv[0])
        usage()
    gParamsDir = os.path.abspath(options.gParamsDir)
    if(options.outDir != None):
        if(not os.path.exists(options.outDir)):
            os.mkdir(options.outDir)
        workingDir = options.outDir
    else:
        (workingDir,tail) = os.path.split(options.parentDir)
    outDir = options.outDir
    if (not os.path.isdir(options.gParamsDir)):
        sys.stderr.write('%s: Error, Params directory "%s" not a directory!\n' %  (sys.argv[0], options.gParamsDir))
        usage()
    if (not os.path.isdir(options.parentDir)):
        sys.stderr.write('%s: Error, Parent directory "%s" not a directory!\n' % (sys.argv[0], options.parentDir))
        usage()

    inputNewick = options.inputNewick
    if(options.stepSize == None):
        options.stepSize = stepSize
    if(not options.testTree):
        if(len(glob.glob(parentDir+'/*rev')) == 0):
            sys.stderr.write('Error: no *.rev found in parent dir.\n')
            usage()
        
    ##############################
    # End of pre-processing.
    ##############################
    newickTree = newickTreeParser(inputNewick, 0.0)
    if(newickTree.distance <= 0):
        newickTree.distance = 0
        if(newickTree.internal):
            if(options.verbose):
                sys.stderr.write('at branch point with tree: %s\n' %(printBinaryTree(newickTree,1)))
            nextTree = newickTreeParser(inputNewick, 0.0)
            # branch point!
            ####################
            # LEFT BRANCH
            lTree = newickTree.left

            nextLTree = nextTree.left
            nextLTree.distance = nextLTree.distance - options.stepSize
            if(nextLTree.distance < 0):
                nextLTree.distance = 0
            if(nextLTree.iD):
                if(nextLTree.distance == 0):
                    name=nextLTree.iD
                else:
                    name=nextLTree.iD+str(nextLTree.distance)
            else:
                name=printBinaryTree(nextLTree,1)
            name=fixName(name)
            leftChildPath = os.path.join(workingDir, name)
            if (lTree.distance < options.stepSize):
                cycleStepSize = lTree.distance
            else:
                cycleStepSize = options.stepSize
            # left branch command
            leftChildCMD = CYCLE_SH +\
                       ' '+parentDir+\
                       ' '+leftChildPath+\
                       ' '+gParamsDir+\
                       ' '+str(cycleStepSize)

            ##################################################
            # TEST TREE COMMAND
            if(options.testTree):
                leftChildCMD = 'mkdir '+ leftChildPath
                sys.stderr.write('left: now planning on running \n\t%s\n' %leftChildCMD)
                if(options.logBranch):
                    (head, tail) = os.path.split(leftChildPath)
                    branchLog( '%25s: %s\n' % ('left mkdir',tail))
            ##################################################
            # spawn a parallel process to run the cycle.
            leftProcess = subprocess.Popen(leftChildCMD, shell=True)

            ####################
            # RIGHT BRANCH
            rTree = newickTree.right
            nextRTree = nextTree.right
            nextRTree.distance = nextRTree.distance - options.stepSize
            if(nextRTree.distance < 0):
                nextRTree.distance = 0
            if(nextRTree.iD):
                if(nextRTree.distance == 0):
                    name=nextRTree.iD
                else:
                    name=nextRTree.iD+str(nextRTree.distance)
            else:
                name=printBinaryTree(nextRTree,1)
            name=fixName(name)
            rightChildPath = os.path.join(workingDir, name)
            if (rTree.distance < options.stepSize):
                cycleStepSize = rTree.distance
            else:
                cycleStepSize = options.stepSize
            # right branch command
            rightChildCMD = CYCLE_SH +\
                       ' '+parentDir+\
                       ' '+rightChildPath+\
                       ' '+gParamsDir+\
                       ' '+str(cycleStepSize)
            ##################################################
            # TEST TREE COMMAND
            if(options.testTree):
                rightChildCMD = 'mkdir '+ rightChildPath
                sys.stderr.write('right: now planning on running \n\t%s\n' %rightChildCMD)
                if(options.logBranch):
                    (head, tail) = os.path.split(rightChildPath)
                    branchLog( '%25s: %s\n' % ('right mkdir',tail))
            ##################################################
            # spawn a parallel process to run the cycle
            rightProcess = subprocess.Popen(rightChildCMD, shell=True)

            ##########
            # wait for the subprocesses to return with no errors.
            # they should both take about the same amount of time to complete
            leftProcess.wait()
            rightProcess.wait()            
            if(leftProcess.returncode or rightProcess.returncode):
                sys.stderr.write('%s: Error in one of the branch cycle sub-processes,\n  left exit status: %s\n  right exit status: %s\n' %(sys.argv[0], leftChildCMD, rightChildCMD))
                sys.exit(1)

            #####
            # Left Recursion
            treeString = printBinaryTree(nextLTree, 1)
            leftRecursiveCommand= sys.argv[0] +\
                                  ' --parent ' + leftChildPath +\
                                  ' --tree "'  + treeString +'"'+ \
                                  ' --params '+ gParamsDir + \
                                  ' --step '   + str(options.stepSize)
            if (outDir != None):
                leftRecursiveCommand = leftRecursiveCommand + ' --out ' + outDir
            leftRecursiveCommand = leftRecursiveCommand + ' --isContinue '
            if (options.testTree):
                leftRecursiveCommand = leftRecursiveCommand + ' --testTree '
            if(options.saveParent):
                leftRecursiveCommand = leftRecursiveCommand + ' --saveParent '
            if(options.verbose):
                leftRecursiveCommand = leftRecursiveCommand + ' --verbose '
            if(options.logBranch):
                leftRecursiveCommand = leftRecursiveCommand + ' --logBranch '
            if (options.isContinue):
                ####################
                # Write this command to parent dir, in case of downstream failure
                commandPath = os.path.join(parentDir, 'nextCommand_left.txt')
                FILE = open(commandPath, 'w')
                FILE.write( '%s\n' % (leftRecursiveCommand))
                FILE.close()
            leftProcess = subprocess.Popen(leftRecursiveCommand, shell=True)
            
            #####
            # Right Recursion
            treeString = printBinaryTree(nextRTree, 1)
            rightRecursiveCommand= sys.argv[0] +\
                                   ' --parent ' + rightChildPath +\
                                   ' --tree "'  + treeString +'"' +\
                                   ' --params '+ gParamsDir + \
                                   ' --step '   + str(options.stepSize)
            if (outDir != None):
                rightRecursiveCommand = rightRecursiveCommand + ' --out ' + outDir
            rightRecursiveCommand = rightRecursiveCommand + ' --isContinue '
            if (options.testTree):
                rightRecursiveCommand = rightRecursiveCommand + ' --testTree '
            if(options.saveParent):
                rightRecursiveCommand = rightRecursiveCommand + ' --saveParent '
            if(options.verbose):
                rightRecursiveCommand = rightRecursiveCommand + ' --verbose '
            if(options.logBranch):
                rightRecursiveCommand = rightRecursiveCommand + ' --logBranch '
            if (options.isContinue):
                ####################
                # Write this command to parent dir, in case of downstream failure
                commandPath = os.path.join(parentDir, 'nextCommand_right.txt')
                FILE = open(commandPath, 'w')
                FILE.write( '%s\n' % (rightRecursiveCommand))
                FILE.close()
            rightProcess = subprocess.Popen(rightRecursiveCommand, shell=True)

            leftProcess.wait()
            rightProcess.wait()
            if(leftProcess.returncode):
                sys.stderr.write('%s: Error left sub-processes "%s", returncode: %s.\n' %(sys.argv[0], leftRecursiveCommand, leftProcess.returncode))
                sys.exit(1)
            if(rightProcess.returncode):
                sys.stderr.write('%s: Error right sub-processes "%s", returncode: %s.\n' %(sys.argv[0], rightRecursiveCommand, rightProcess.returncode))
                sys.exit(1)
                                            
    else:
        ####################
        # not a branch, traveling down stem
        ####################
        nextTree = newickTreeParser(inputNewick, 0.0)
        if (nextTree.distance < options.stepSize):
            cycleStepSize = nextTree.distance
            nextTree.distance = 0
        else:
            nextTree.distance = nextTree.distance - options.stepSize
            cycleStepSize = options.stepSize
        if(nextTree.iD):
            if(nextTree.distance == 0):
                name=nextTree.iD
            else:
                name=nextTree.iD+str(nextTree.distance)
        else:
            name=printBinaryTree(nextTree,1)
        name=fixName(name)
        
        if(options.verbose):
                sys.stderr.write('At stem, node %s, with tree: %s\n' %(name, printBinaryTree(nextTree,1)))
        childPath = os.path.join(workingDir, name)

        treeString=printBinaryTree(nextTree, 1)
        
        childCMD = CYCLE_SH +\
                   ' '+parentDir+\
                   ' '+childPath+\
                   ' '+gParamsDir+\
                   ' '+str(cycleStepSize)

        ##################################################
        # TEST TREE COMMAND
        if(options.testTree):
            childCMD = 'mkdir '+ childPath
            sys.stderr.write('stem: now planning on running \n\t%s\n' %childCMD)
            if(options.logBranch):
                (head, tail) = os.path.split(childPath)
                branchLog( '%25s: %s\n' % ('stem mkdir',tail))
        ##################################################
        returnCode = subprocess.Popen(childCMD, shell=True).wait()
        if(returnCode):
            sys.stderr.write('%s: Error in cycle sub-processes, exit status: %s\n' % (sys.argv[0], returnCode))
            sys.exit(1)

        recursiveCommand= sys.argv[0] +\
                          ' --parent ' + childPath+\
                          ' --tree "'  + treeString+'"'+\
                          ' --params '+ gParamsDir+\
                          ' --step '   + str(options.stepSize)
        if (outDir != None):
            recursiveCommand = recursiveCommand + ' --out ' + outDir
            recursiveCommand = recursiveCommand + ' --isContinue '
        if (options.testTree):
            recursiveCommand = recursiveCommand + ' --testTree '
        if (options.saveParent):
            recursiveCommand = recursiveCommand + ' --saveParent '
        if (options.verbose):
            recursiveCommand = recursiveCommand + ' --verbose '
            sys.stderr.write('PEPARING TO RECURSE FROM STEM: %s\n' %(recursiveCommand))
        if (options.logBranch):
            recursiveCommand = recursiveCommand + ' --logBranch '
        if (options.isContinue):
            ####################
            # Write this command to parent dir, in case of downstream failure
            commandPath = os.path.join(parentDir, 'nextCommand.txt')
            FILE = open(commandPath, 'w')
            FILE.write( '%s\n' % (recursiveCommand))
            FILE.close()
        subprocess.Popen(recursiveCommand, shell=True)

    newickTree = newickTreeParser(inputNewick, 0.0) # used to check if this is a leaf.
    if((not options.saveParent) and (options.isContinue) and newickTree.distance > 0): # settings check 
        shutil.rmtree(parentDir) # burn it to the ground
        if(options.logBranch):
            (head, tail) = os.path.split(parentDir)
            branchLog( '%25s: %s. saveParent:%d isContinue:%d distance:%f\n' % ('performing rmtree() on ', tail, options.saveParent, options.isContinue, newickTree.distance))
    else:
        if(options.logBranch):
            (head, tail) = os.path.split(parentDir)
            branchLog( 'will not delete parent (%s). saveParent:%d isContinue:%d distance:%f\n' % (tail, options.saveParent, options.isContinue, newickTree.distance))

if __name__ == "__main__":
    main()
