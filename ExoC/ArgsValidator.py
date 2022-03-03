import os,sys
from ExoC.corelib import *

def opt_validate_stat(optparser):
    """Validate arguments from a OptParser object.
    """

    opt = optparser.parse_args()
    if not os.path.isfile(opt.sam):
        optparser.print_help()
        Info('ERROR: path to the SAM file is not valid')
        sys.exit(1)
        
    if not opt.name:
        opt.name = "NA"
        Info("ExoC will use 'NA' as the prefix for all the ouput files")

    if not opt.outdir:
        opt.outdir = os.getcwd()
        opt.outdir = os.path.join(opt.outdir,'stat_out')
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)
    else:
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)


    Info("Argument List: ")
    Info("SAM file = " + opt.sam)
    Info("Sample name = " + opt.name)
    Info("Output directory = " + opt.outdir)

    return opt

def opt_validate_cnv(optparser):
    """
    Validate arguments from a OptParser object.
    """

    opt = optparser.parse_args()
    if not os.path.isfile(opt.case):
        optparser.print_help()
        Info('ERROR: path to the case validpairs file is not valid')
        sys.exit(1)
        
    if not os.path.isfile(opt.control):
        optparser.print_help()
        
        Info('ERROR: path to the control validpairs file is not valid')
        sys.exit(1)
        
    if not os.path.isfile(opt.cnv):
        optparser.print_help()
        Info('ERROR: path to the CNV bed file is not valid')
        sys.exit(1)
        
    if not opt.name:
        opt.name = "NA"
        Info("ExoC will use 'NA' as the prefix for all the ouput files")

    if not opt.outdir:
        opt.outdir = os.getcwd()
        opt.outdir = os.path.join(opt.outdir,'stat_out')
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)
    else:
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)


    Info("Argument List: ")
    Info("Case validpairs file = " + ', '.join(opt.case))
    Info("Control validpairs file = " + ', '.join(opt.control))
    Info("Segmented CNV bed file = " + opt.cnv)
    Info("Sample name = " + opt.name)
    Info("Output directory = " + opt.outdir)

    return opt

def opt_validate_trans(optparser):
    """
    Validate arguments from a OptParser object.
    """

    
    opt = optparser.parse_args()
    if not os.path.isfile(opt.case):
        optparser.print_help()
        Info('ERROR: path to the case validpairs file is not valid')
        sys.exit(1)
        
    if not os.path.isfile(opt.control):
        optparser.print_help()
        
        Info('ERROR: path to the control validpairs file is not valid')
        sys.exit(1)
        
    if not opt.binsize:
        opt.binsize = 10000
    else:
        if int(opt.binsize) <= 0:
            optparser.print_help()
            Info('ERROR: Binsize value is not valid')
            sys.exit(1)
        
    if not opt.name:
        opt.name = "NA"
        Info("ExoC will use 'NA' as the prefix for all the ouput files")

    if not opt.outdir:
        opt.outdir = os.getcwd()
        opt.outdir = os.path.join(opt.outdir,'stat_out')
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)
    else:
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)
            
    if not opt.thr_frame:
        opt.thr_frame = 0


    Info("Argument List: ")
    Info("Case validpairs file = " + opt.case)
    Info("Control validpairs file = " + opt.control)
    Info("Binsize = " + opt.binsize)
    Info("Sample name = " + opt.name)
    Info("Output directory = " + opt.outdir)

    return opt
