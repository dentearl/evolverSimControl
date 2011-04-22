def standardOptions(parser):
    parser.add_option('-c', '--child',dest='childDir',
                      help='Child directory.')

def standardOptionsCheck(options, usage):
    import os, sys
    if (options.childDir == None):
        sys.stderr.write('%s: Error, specify child dir.\n' % sys.argv[0])
        usage()
    options.childDir   = os.path.abspath(options.childDir)
    options.theChild   = os.path.basename(options.childDir)

