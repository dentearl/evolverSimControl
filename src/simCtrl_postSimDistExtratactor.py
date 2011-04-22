#!/usr/bin/env python
"""
postSimDistExtractor.py
dent earl, dearl (a) soe ucsc edu
29 June 2010
a script to be run following the completion of a
simulation. The script will check the simulationInfo.xml
file and then open the leaves' and roots'
stats/annotstats.txt file and extract the distributions of
annotations, printing them out in an R readible format for
plotting and analysis.
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

def usage():
    print 'USAGE: '+sys.argv[0]+' --simDir <dir>'
    print __doc__
    sys.exit(2)

def initOptions(parser):
    parser.add_option('-d', '--simDir',dest='simDir',
                      help='Simulation directory.')
    parser.add_option('-i', '--internalBranches',dest='internalBranches',
                      action='store_true', default='false',
                      help='Turns on output of the internal braches. Otherwise it is just root and leaves.')

def checkOptions(options):
    if (options.simDir == None):
        sys.stderr.write('%s: Error, specify simulation dir.\n' % sys.argv[0])
        usage()
    options.simDir  = os.path.abspath(options.simDir)
    if not os.path.exists(os.path.join(options.simDir, 'simulationInfo.xml')):
        sys.stderr.write('%s: Error, unable to find simulationInfo.xml.\n' % sys.argv[0])
        usage()
    infoTree = ET.parse(os.path.join(options.simDir, 'simulationInfo.xml'))
    treeObj = infoTree.find('tree')
    options.inputNewick=treeObj.text
    treeObj = infoTree.find('rootDir')
    options.rootDir = treeObj.text
    options.rootName = os.path.basename( options.rootDir )

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

def extractLeavesAndIntBranches( nt, options, leafDict ):
    """Given a newick tree object, it returns a dict of
    leaf and internal branch objects. Operates recursively.
    """
    if nt == None:
        return None
    nt.distance = 0
    if nt.right == None and nt.left == None:
        # leaf
        leafDict[ nt.iD ] = True
    else:
        if options.internalBranches:
            if nt.right != None and nt.left != None and nt.iD != options.rootName:
                # internal branch
                leafDict[ nt.iD ] = True
        extractLeavesAndIntBranches( nt.right, options, leafDict = leafDict )
        extractLeavesAndIntBranches( nt.left , options, leafDict = leafDict )

def parseStats(options, leaves):
    results={}
    for l in leaves:
        results[l] = extractDist( options, l )

    return results

def extractDist(options, leaf):
    infile = open(os.path.join(options.simDir, leaf, 'stats', 'annotstats.txt'), 'r')
    typePat = re.compile('^Length dist. (\S+):')
    summPat = re.compile('^Avg (\S+), N=(\d+), total (\d+)')
    distPat = re.compile('^\s+(\d+)\s+-\s+(\d+)\s+(\d+)')
    results={}
    while infile:
        line = infile.readline()
        if not line:
            break
        p = typePat.match(line)
        if p != None:
            type = p.group(1)
            line = infile.readline()
            p = summPat.match(line)
            results[type] = {'ave':float(p.group(1)),
                             'n':int(p.group(2)),
                             'total':int(p.group(3)),
                             'rangeStart':[],
                             'rangeEnd':[],
                             'values':[]}
            infile.readline()
            infile.readline()
            line = infile.readline()
            while line != '\n':
                p = distPat.match(line)
                if p != None:
                    results[type]['rangeStart'].append(int(p.group(1)))
                    results[type]['rangeEnd'].append(int(p.group(2)))
                    results[type]['values'].append(int(p.group(3)))
                line=infile.readline()
    return results

def printStats(options, results):
    for r in results:
        dists = results[r]
        for d in dists:
            print '%s.%s.ave = %f' %(r, d, dists[d]['ave'])
            print '%s.%s.n = %d' %(r, d, dists[d]['n'])
            print '%s.%s.total = %d' %(r, d, dists[d]['total'])
            sys.stdout.write('%s.%s.range = c(' %(r, d))
            for i in range(len(dists[d]['rangeStart'])):
                if i == len(dists[d]['rangeStart']) - 1:
                    s = ')\n'
                else:
                    s = ', '
                sys.stdout.write("'%s - %s'%s" %(dists[d]['rangeStart'][i],
                                                 dists[d]['rangeEnd'][i],s))
            sys.stdout.write('%s.%s.value = c(' %(r, d))
            for i in range(len(dists[d]['values'])):
                if i == len(dists[d]['values']) - 1:
                    s = ')\n'
                else:
                    s = ', '
                sys.stdout.write('%s%s' %(dists[d]['values'][i], s))

def printScript(options, results):
    print 'pdf("smallMultiples.pdf", height=15, width=30)'
    print 'par(mfrow=c(%d, 4))' %len(results)
    order = ['CDS', 'UTR', 'introns', 'genes']
    lengths = {'CDS':'',
               'UTR':'',
               'introns':'',
               'genes':''}
    colors = {'CDS':'rgb(248/255, 118/255, 109/255)',
               'UTR':'rgb(124/255, 174/255, 0)',
               'introns':'rgb(0, 191/255, 196/255)',
               'genes':'rgb(199/255, 124/255, 255/255)'}
    
    
    for o in order:
        print "barplot(%s.%s.value%s, names.arg='', border=NA, las=2, main='%s.%s', col=%s)" %(
            options.rootName, o, lengths[o], options.rootName, o, colors[o])
    i = 0
    for r in results:
        i += 1
        if r != options.rootName:
            for o in order:
                if i == len(results):
                    print "barplot(%s.%s.value%s, names.arg=%s.%s.range%s, border=NA, las=2, main='%s.%s', col=%s)" % (
                        r, o, lengths[o], r, o, lengths[o], r, o, colors[o])
                else:
                    print "barplot(%s.%s.value%s, names.arg='', border=NA, las=2, main='%s.%s', col=%s)" % (
                        r, o, lengths[o], r, o, colors[o])
    print 'dev.off()'

def standardizeResults(options, results):
    """ standardize the counts so that the plots will look okay
    """
    types = results[ options.rootName ].keys()
    typesDict = {}
    for t in types:
        if t not in typesDict:
            typesDict[t] = { 'minName'       : 'paul',
                             'minRangeStart' : 100000000000,
                             'minRangeEnd'   : 100000000000,
                             'maxName'       : 'steven',
                             'maxRangeStart' : 0,
                             'maxRangeEnd'   : 0 }
        for r in results:
            if results[r][t]['rangeStart'][0] < typesDict[t]['minRangeStart']:
                typesDict[t]['minName'] = r
                typesDict[t]['minRangeStart'] = results[r][t]['rangeStart'][0]
                typesDict[t]['minRangeEnd'] = results[r][t]['rangeEnd'][0]
            if results[r][t]['rangeEnd'][-1] > typesDict[t]['maxRangeEnd']:
                typesDict[t]['maxName'] = r
                typesDict[t]['maxRangeStart'] = results[r][t]['rangeStart'][-1]
                typesDict[t]['maxRangeEnd'] = results[r][t]['rangeEnd'][-1]
    # standardize all of the ranges to be the same, for samples
    # without a value at a particular range, append a 0 into their 'values'
    for t in typesDict:
        #####
        # Do the low side of the range first,
        for r in results:
            newVals = []
            newStarts = []
            newEnds = []
            i = 0
            while results[typesDict[t]['minName']][t]['rangeStart'][i] != results[r][t]['rangeStart'][0]:
                newVals.append(0)
                newStarts.append(results[ typesDict[t]['minName'] ][t]['rangeStart'][i])
                newEnds.append(results[ typesDict[t]['minName'] ][t]['rangeEnd'][i])
                i += 1
            if i:
                newVals.extend(results[r][t]['values'])
                newStarts.extend(results[r][t]['rangeStart'])
                newEnds.extend(results[r][t]['rangeEnd'])
                results[r][t]['values'] = newVals
                results[r][t]['rangeStart'] = newStarts
                results[r][t]['rangeEnd'] = newEnds
        #####
        # Do the highside of the range second,
        for r in results:
            newVals = []
            newStarts = []
            newEnds = []
            i = 0
            if len(results[ typesDict[t]['maxName'] ][t]['rangeEnd']) != len(results[r][t]['rangeEnd']):
                if typesDict[t]['maxRangeStart'] == results[r][t]['rangeStart'][-1] and \
                       typesDict[t]['maxRangeEnd'] != results[r][t]['rangeEnd'][-1]:
                    results[r][t]['rangeEnd'][-1] = typesDict[t]['maxRangeEnd']
                else:
                    for i in range(len(results[r][t]['values']), len(results[ typesDict[t]['maxName'] ][t]['values'])):
                        results[r][t]['values'].append(0)
                        results[r][t]['rangeStart'].append(results[ typesDict[t]['maxName'] ][t]['rangeStart'][i])
                        results[r][t]['rangeEnd'].append(results[ typesDict[t]['maxName'] ][t]['rangeEnd'][i])
            else:
                results[r][t]['rangeEnd'][-1] = results[ typesDict[t]['maxName'] ][t]['rangeEnd'][ len(results[r][t]['values']) - 1 ]
                

def printScriptGG(options, results):
    """script to print out R commands to use ggplot2 library to create a pretty plot.
    """
    outName  = []
    outRange = []
    outValue = []
    outType  = []
    types = results[ options.rootName ].keys()
    for r in results:
        for t in types:
            if t != 'NGE' and t != 'NXE':
                for i in range(len(results[r][t]['rangeStart'])):
                    myRange = '"%s - %s"' %(results[r][t]['rangeStart'][i], results[r][t]['rangeEnd'][i])
                    myName = '"%s"' %r
                    myType = '"%s"' %t
                    outName.append(myName)
                    outType.append(myType)
                    outRange.append(myRange)
                    outValue.append(results[r][t]['values'][i])
    sys.stdout.write( 'outName = c(')
    for i in range(len(outName)):
        if i == len(outName) - 1:
            s = ')\n'
        else:
            s = ', '
        sys.stdout.write("%s%s" %(outName[i], s))
    sys.stdout.write( 'outValue = c(')
    for i in range(len(outValue)):
        if i == len(outValue) - 1:
            s = ')\n'
        else:
            s = ', '
        sys.stdout.write("%s%s" %(outValue[i], s))
    sys.stdout.write( 'outRange = c(')
    for i in range(len(outRange)):
        if i == len(outRange) - 1:
            s = ')\n'
        else:
            s = ', '
        sys.stdout.write("%s%s" %(outRange[i], s))
    sys.stdout.write( 'outType = c(')
    for i in range(len(outType)):
        if i == len(outType) - 1:
            s = ')\n'
        else:
            s = ', '
        sys.stdout.write("%s%s" %(outType[i], s))
    sys.stdout.write('outName = factor(outName, levels=c(\'%s\', ' % options.rootName)
    alphabetized = results.keys()
    alphabetized.sort(key=lambda x: x.lower())
    for i in range(len(alphabetized)):
        if alphabetized[i] != options.rootName:
            if i == len(alphabetized) - 1:
                s = '))\n'
            else:
                s = ', '
            sys.stdout.write("'%s'%s" %(alphabetized[i], s))
    print 'df = data.frame(Range=outRange, Name=outName, Value=outValue, Type=outType)'
    print 'require(ggplot2)'
    print '# qplot(x=Range, y=Value, data=df, geom=\'bar\', facets = Name ~ Type, color=Type, scales="free")'
    print 'mt = ggplot(df, aes(x=factor(Range, levels=factor(Range)), y=Value)) + geom_bar()'
    print '# mt + facet_grid(Type ~ Name, scales="free") + aes(fill=Type)'
    print '# mt + opts(axis.text.x = theme_text(size = 4, angle=45, hjust=1)) + facet_wrap(Name ~ Type, scales="free", ncol=4) + aes(fill=Type)'

def main():
    parser=OptionParser()
    initOptions( parser )
    ( options, args ) = parser.parse_args()
    checkOptions( options )
    nt = newickTreeParser( options.inputNewick, 0.0)
    if nt.iD == None:
        nt.iD= options.rootName
        
    leaves={}
    extractLeavesAndIntBranches( nt, options, leaves )
    leaves[ options.rootName ] = True
    results = parseStats( options, leaves )
    standardizeResults( options, results )
    printStats( options, results )
    printScript( options, results )
    printScriptGG( options, results )

if __name__ == "__main__":
    main()
