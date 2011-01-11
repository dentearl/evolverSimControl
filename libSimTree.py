def nameTree(nt, reportDistance=1):
    """nameTree(nt) takes a newick tree and returns a str that can be used
    to name the cycle-step that the tree represents. Distance included in
    the name by default.
    """
    from sonLib.bioio import printBinaryTree
    if nt == None:
        return ''
    if(nt.iD):
        if(nt.distance == 0 or not reportDistance):
            name=nt.iD
        else:
            name=nt.iD+str(nt.distance)
    else:
        name=printBinaryTree(nt,1)
    name=fixName(name)
    return name

def newickRootName( nt ):
    if nt.iD:
        return nt.iD
    else:
        return 'root'

def fixName(name):
    """fixName(name) takes all the nasty characters out of a newickTree and
    returns a str that is more amenable to being a file (or directory) name.
    """
    name=name.replace(' ','')
    name=name.replace(',','')
    name=name.replace(':','-')
    name=name.replace('.','_')
    name=name.replace(';','')
    name=name.replace('\'','')
    name=name.replace('"','')
    name=name.replace('(','_L_')
    name=name.replace(')','_R_')
    name=name.rstrip('0')
    name=name.rstrip('-0_')
    return name

def branchLog(message):
    """branchLog(message) sends a message to the branching log
    """
    import os
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

def standardOptions(parser):
    """Takes in an OptionParser parser object and adds a series of standard options.
    """
    parser.add_option('-o', '--out',dest='outDir',
                      help='Out directory.')
    parser.add_option('-t', '--tree',dest='inputNewick',
                      help='Newick tree.')
    parser.add_option('-y', '--removeParent', action='store_true', dest='removeParent',
                      default=False, help='Remove cycles that are not nodes in order to save disk space (not recommended).')
    parser.add_option('-b', '--isBranchChild', action='store_true', dest='isBranchChild',
                      default=False, help='Establishes the parent as un-deletable.')
    parser.add_option('-c', '--isContinue', action='store_true', dest='isContinue',
                      default=False, help='Is this a recursive run?')
    parser.add_option('-l', '--logBranch', action='store_true', dest='logBranch',
                      default=False, help='record all tree travel information.')
    parser.add_option('-T', '--testTree', action='store_true', dest='testTree',
                      default=False, help='Instead of performing a simulation, does dry run with empty dirs.')

def standardOptionsCheck(options, usage):
    """Takes in an OptionParser options object and checks the standard options.
    """
    import os, sys
    if (options.inputNewick == None):
        sys.stderr.write('%s: Error, specify newick.\n' % sys.argv[0])
        usage()
    if(options.outDir == None):
        sys.stderr.write('%s: Error, specify out directory.\n' % sys.argv[0])
        usage()
    if(not os.path.exists(options.outDir)):
        os.mkdir(options.outDir)

def cycleBranchCommandBuilder(nt, nxt, branchString, stepSize, workingDir, parentDir, gParamsDir, seed, logBranch, isTestTree):
    import simulation.lib.libSimTree as LST
    import simulation.lib.libSimControl as LSC
    import os
    programs = ['simCtrl_cycleMain_1.py']
    LSC.verifyPrograms(programs)
    CYCLEBEGIN_PY = programs[0]
    nxt.distance = nxt.distance - stepSize
    if(nxt.distance < 0):
        nxt.distance = 0
    name = LST.nameTree(nxt)
    childPath = os.path.join(workingDir, name)
    if (nt.distance < stepSize):
        cycleStepSize = nt.distance
    else:
        cycleStepSize = stepSize
    # branch command
    if seed != 'random':
        seed = int(seed)
        if branchString == 'left':
            seed += 1
        elif branchString == 'right':
            seed -= 1
        else:
            seed += 10
        seed = abs(seed)
    childCMD = CYCLEBEGIN_PY +\
               ' --parent ' + parentDir +\
               ' --child '  + childPath +\
               ' --params ' + gParamsDir +\
               ' --step '   + str(cycleStepSize) +\
               ' --seed '   + str(seed) +\
               ' --jobFile JOB_FILE'
    ##################################################
    # TEST TREE COMMAND
    if(isTestTree):
        childCMD = 'mkdir '+ childPath + ' && echo JOB_FILE  > '+childPath+'/rev'
        if(logBranch):
            (head, tail) = os.path.split(childPath)
            LST.branchLog( '%25s: %s\n' % (branchString+' mkdir',tail))
    return childCMD

def commandRecorder(cmd, dir):
    import xml.etree.ElementTree as ET
    import os, sys
    inputXML=os.path.join(dir, 'cycleInfo.xml')
    if not os.path.exists(inputXML):
        sys.stderr.write('Error, unable to locate %s.\n' % (inputXML))
        sys.exit(1)
    infoTree=ET.parse(inputXML)
    root=infoTree.getroot()
    tObj=ET.SubElement(root, 'followUpCommand')
    cmd = cmd.replace('"', '\'')
    tObj.text=cmd
    info=ET.ElementTree(root)
    info.write(inputXML)

def treeBranchCommandBuilder(nt, branchStr, options, commonParent, gParamsDir):
    from sonLib.bioio import printBinaryTree
    import simulation.lib.libSimTree as LST
    import simulation.lib.libSimControl as LSC
    programs = ['simCtrl_simTree.py']
    LSC.verifyPrograms(programs)
    SIMTREE_PY = programs[0]
    treeString = printBinaryTree(nt, 1)
    treeString = treeString.rstrip(':.0;') # necessary due to weird newickTree code
    if options.seed != 'random':
        options.seed = int(options.seed)
        if branchStr == 'Left Branch':
            options.seed += 1
        elif branchStr == 'Right Branch':
            options.seed -= 1
        else:
            options.seed += 10
        options.seed = abs(options.seed)
    childCMD = SIMTREE_PY +\
               ' --parent '+commonParent+\
               ' --tree "' +treeString+'"'+\
               ' --params '+gParamsDir+\
               ' --step '  +str(options.stepSize) +\
               ' --seed '  +str(options.seed)+\
               ' --jobFile JOB_FILE'
    if (options.outDir != None):
        childCMD = childCMD + ' --out ' + options.outDir
    if options.isContinue:
        childCMD = childCMD + ' --isContinue '
    if branchStr != 'Stem':
        childCMD = childCMD + ' --isBranchChild '
    if options.logBranch:
        childCMD = childCMD + ' --logBranch '
    if options.removeParent:
        childCMD = childCMD + ' --removeParent '
    if options.noMEs:
        childCMD = childCMD + ' --noMEs '
    if options.testTree:
        childCMD = childCMD + ' --testTree '
    if options.logBranch:
        LST.branchLog( '%25s: %s\n' % (branchStr+' to simTree.py',childCMD))
    if not options.testTree:
        LST.commandRecorder(childCMD, commonParent)
    return childCMD
