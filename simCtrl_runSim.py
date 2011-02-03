#!/usr/bin/env python
"""
runSim.py
dent earl, dearl (a) soe ucsc edu
19 April 2009
In order to run a simulation you call runSim.py
runSim.py checks to make sure that every single
piece of software called (eventually) by the simulation
will exist and it does some prep work before the first
call to simTree.py, which begins the recursive process
of performing the genome simulation.
"""
##############################
from sonLib.bioio import newickTreeParser
from sonLib.bioio import printBinaryTree
from optparse import OptionParser
import xml.etree.ElementTree as ET
import os
import sys
import time
from datetime import datetime
import simulation.lib.libSimControl as LSC
import simulation.lib.libSimTree as LST

# this is the first script to run in a simulation, so it will check for the
# existance of everything the entire simulation will end up calling, not
# just the scripts or external files used by runSim.py.
programs = ['simCtrl_simTree.py','simCtrl_rootCycleInfoCreator.py', 'cp',
            'mkdir', 'simCtrl_commandEval.py',
            'evolver_evo', 'evolver_cvt', 'evolver_transalign',
            'evolver_drawrev', 'evolver_gff_cdsutr2exons.py',
            'evolver_gff_exons2introns.py', 'evolver_gff_featurestats2.sh',
            'evolver_gff_featurestats2.py', 'evolver_codon_report.pl',
            'evolver_merge_evostats.py','evolver_mobile_report.pl',
            'touch', 'ln','egrep', 'cat',
            'simCtrl_cycleMain_1.py', 'simCtrl_cycleMain_2.py',
            'simCtrl_cycleMain_3.py', 'simCtrl_cycleMain_4.py',
            'simCtrl_cycleStats_1.py', 'simCtrl_cycleStats_2.py',
            'simCtrl_cycleStats_3.py', 'simCtrl_cycleStats_4.py',
            'simCtrl_simTreeFollowUp.py', 'simCtrl_wrapperTRF.py',
            'simCtrl_completeTimestamp.py']
LSC.verifyPrograms( programs )
( SIMTREE_PY, ROOT_CYCLEXML_MAKER,
  CP_BIN, MKDIR_BIN, CMD_EVAL_BIN ) = programs[ 0:5 ]

def usage():
    print 'USAGE: '+sys.argv[0]+' --root [dir] --out --tree [newick tree in quotes] --params [parameter dir] --stepSize [0.001] --jobFile JOB_FILE '
    print __doc__
    sys.exit(2)

def newickContainsReservedWord(nt):
    """newickContainsReservedWord() checks the newick to make sure that
    there are no reserved names used as IDs. At present the only reserved
    name is 'parameters'.
    """
    reservedWords = {'parameters':1}
    if nt == None:
        return False
    if nt.iD in reservedWords:
        return True
    left = newickContainsReservedWord(nt.right)
    right= newickContainsReservedWord(nt.left)
    if left or right:
        return True

"""
Okay, here's what we do: The first call to simTree.py will ALWAYS perform
the task of copying over the parent genome into a directory called 'root'
and will create the cycleInfo.xml file. Then, if the newick tree starts at a branch
then we will recall simTree.py using the root as the parent. there must be a variable
passed to branchcommandbuilder through the branches that represents this.
If instead we are on a stem to begin with then the only change will be to make
the new root directory the parent directory for the first real cycle.
"""

def initOptions(parser):
    """initOptions() initializes the options that will go to the
    parser object
    """
    parser.add_option('--isFollowUp',dest='isFollowUp', action='store_true',
                      default=False, help='file to write run-stats to.')
    parser.add_option('--rootName',dest='rootName', 
                      help='name of the root genome, to differentiate it from the input newick.')

def main():
    parser=OptionParser()
    LSC.standardOptions(parser)
    LST.standardOptions(parser)
    initOptions(parser)
    (options, args) = parser.parse_args()
    LSC.standardOptionsCheck(options, usage)
    LST.standardOptionsCheck(options, usage)
    # check newickTree for reserved words
    nt = newickTreeParser( options.inputNewick, 0.0 )
    if newickContainsReservedWord(nt):
        sys.stderr.write('%s: Error: newick tree contains reserved words.\n' %(sys.argv[0]))
        usage()
    if not options.rootName:
        rootName = LST.newickRootName( nt )
    else:
        rootName = options.rootName
    if not os.path.exists( os.path.join( options.parentDir, 'seq.rev')):
        sys.stderr.write('ERROR: Unable to find seq.rev in --parentDir [%s].\n' % options.parentDir)
        sys.exit(1)
    xmlTree = ET.parse(options.jobFile)
    childrenElm = xmlTree.find('children')
    if not options.isFollowUp:
        # if the root/ dir does not exist inside the simulation directory
        # then we need to create it and then re-call runSim
        childCMD = CMD_EVAL_BIN+' JOB_FILE "'
        childCMD+= LSC.commandPacker(CP_BIN+\
                                     ' -r '+options.parentDir+\
                                     ' '+os.path.join( options.outDir, rootName ))
        childCMD+= LSC.commandPacker( ROOT_CYCLEXML_MAKER +\
                                     ' --dir '+os.path.join(options.outDir, rootName ) )
        childCMD+= LSC.commandPacker(MKDIR_BIN+\
                                     ' '+os.path.join(options.outDir, 'parameters'))
        childCMD+= LSC.commandPacker(CP_BIN+\
                                     ' '+os.path.join(options.gParamsDir,'model.txt')+\
                                     ' '+os.path.join(options.outDir, 'parameters'))
        childCMD+= LSC.commandPacker(CP_BIN+\
                                     ' '+os.path.join(options.gParamsDir,'model.mes.txt')+\
                                     ' '+os.path.join(options.outDir, 'parameters'))
        if not options.noMEs:
            childCMD+= LSC.commandPacker(CP_BIN+\
                                         ' '+os.path.join(options.gParamsDir,'mes.cfg')+\
                                         ' '+os.path.join(options.outDir, 'parameters'))
        childCMD+= '"'
        options.gParamsDir = os.path.join(options.outDir, 'parameters')
        newChild = ET.SubElement(childrenElm, 'child')
        newChild.attrib['command']=childCMD
        options.parentDir = os.path.join(options.outDir, rootName)
        followUpCommand = sys.argv[0] +\
                          ' --parent '+options.parentDir+\
                          ' --tree "'+options.inputNewick + '"'+\
                          ' --params '+options.gParamsDir +\
                          ' --step '+str( options.stepSize ) +\
                          ' --seed '+options.seed +\
                          ' --jobFile JOB_FILE'+\
                          ' --isFollowUp '
        if options.outDir != None:
            followUpCommand = followUpCommand + ' --out ' + options.outDir
        if options.removeParent:
            followUpCommand = followUpCommand + ' --removeParent'
        if options.isBranchChild:
            followUpCommand = followUpCommand + ' --isBranchChild'
        if options.noMEs:
            followUpCommand = followUpCommand + ' --noMEs'
        if options.testTree:
            followUpCommand = followUpCommand + ' --testTree'
        if options.logBranch:
            followUpCommand = followUpCommand + ' --logBranch'
        jobElm=xmlTree.getroot()
        jobElm.attrib['command'] = followUpCommand
        if (options.logBranch):
            dt=datetime.utcnow()
            nowStr = dt.strftime("%a, %d %b %Y %H:%M:%S +0000")
            separator= '####################\n'
            LST.branchLog('%s%s : Starting new run with tree %s\n%s' %( separator, nowStr, options.inputNewick, separator ))

        if(os.path.exists(os.path.join(options.outDir, 'simulationInfo.xml'))):
            os.remove(os.path.join(options.outDir, 'simulationInfo.xml'))
        root=ET.Element( 'info' )
        tObj=ET.SubElement( root, 'rootDir' )
        tObj.text=str( options.parentDir )
        tObj=ET.SubElement( root, 'gParamsDir' )
        tObj.text=str( options.gParamsDir )
        tObj=ET.SubElement( root, 'tree')
        tObj.text=str( options.inputNewick )
        tObj=ET.SubElement( root, 'stepSize' )
        tObj.text=str( options.stepSize )
        tObj=ET.SubElement( root, 'removeParent' )
        tObj.text=str( options.removeParent )
        tObj=ET.SubElement( root, 'timestamps')
        timeStart      = ET.SubElement(tObj,'start')
        timeLocal      = ET.SubElement( timeStart, 'humanLocal' )
        timeLocal.text = str( time.strftime("%a, %d %b %Y %H:%M:%S (%Z) ", time.localtime()) )
        timeHuman      = ET.SubElement( timeStart, 'humanUTC' )
        timeHuman.text = str( time.strftime("%a, %d %b %Y %H:%M:%S (UTC) ", time.gmtime()) )
        timeEpoch      = ET.SubElement( timeStart, 'epochUTC' )
        timeEpoch.text = str(time.time())
        cmd = ET.SubElement( root, 'command')
        cmd.text = ' '.join(sys.argv)
        info=ET.ElementTree( root )
        info.write( os.path.join( options.outDir,'simulationInfo.xml' ) )
    else:
        # since the root/ dir exists, we can make the first call to simTree.py
        # WRITE INFO.XML IF this is the first run of the simulation.
        if nt.distance == 0:
            # branch point, create two child processes
            newChild = ET.SubElement(childrenElm, 'child')
            newChild.attrib['command'] = LST.treeBranchCommandBuilder(nt.left, 'Left Branch', options,
                                                              options.parentDir, options.gParamsDir)
            newChild = ET.SubElement(childrenElm, 'child')
            newChild.attrib['command'] = LST.treeBranchCommandBuilder(nt.right, 'Right Branch', options,
                                                              options.parentDir, options.gParamsDir)
            
        else:
            # stem, create one child process
            newChild = ET.SubElement(childrenElm, 'child')
            newChild.attrib['command'] = LST.treeBranchCommandBuilder(nt, 'Stem', options, options.parentDir,
                                                              options.gParamsDir)

    xmlTree.write(options.jobFile)


if __name__ == "__main__":
    main()
