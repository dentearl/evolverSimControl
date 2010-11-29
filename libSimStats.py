def initOptions(parser):
    parser.add_option('-p', '--parentDir',dest='parentDir',
                      help='Parent directory.')
    parser.add_option('-d', '--childDir',dest='childDir',
                      help='cycle directory.')
    parser.add_option('-j', '--jobFile',dest='jobFile',
                      help='jobFile, passed in by jobTree.')
    
def checkOptions( options, usage ):
    import xml.etree.ElementTree as ET
    import os, sys
    if options.childDir == None:
        sys.stderr.write('%s: Error, specify dir.\n' % sys.argv[0])
        usage()
    if not os.path.isdir(options.childDir):
        sys.stderr.write('%s: Error, directory "%s" is not a directory!\n' % (sys.argv[0], options.dir))
        usage()
    if options.jobFile == None:
        sys.stderr.write('%s: Error, specify jobFile.\n' % sys.argv[0])
        usage()
    if not os.path.exists(options.jobFile):
        sys.stderr.write('%s: Error, jobFile does not exist.\n' % sys.argv[0])
        usage()
    if options.parentDir == None:
        sys.stderr.write('%s: Error, specify parent dir.\n' % sys.argv[0])
        usage()
    if not os.path.isdir(options.parentDir):
        sys.stderr.write('%s: Error, Parent dir "%s" not a directory!\n' % (sys.argv[0], options.parentDir))
        usage()
    options.childDir = os.path.abspath(options.childDir)
    (simDir, tail)=os.path.split(options.childDir)
    if not os.path.exists(os.path.join(simDir,'simulationInfo.xml')):
        sys.stderr.write('%s: Error, simulationInfo.xml does not exist.\n' % sys.argv[0])
        usage()
    xmlTree = ET.parse(os.path.join(simDir,'simulationInfo.xml'))
    jobElm=xmlTree.getroot()
    childrenElm = xmlTree.find('rootDir')
    options.rootDir=os.path.abspath(childrenElm.text)
