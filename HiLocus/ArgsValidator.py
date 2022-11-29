import os,sys
from HiLocus.corelib import *

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
        Info("HiLocus will use 'NA' as the prefix for all the ouput files")

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
        opt.binsize = "10000"
    else:
        if int(opt.binsize) <= 0:
            optparser.print_help()
            Info('ERROR: Binsize value is not valid')
            sys.exit(1)
        
    if not opt.name:
        opt.name = "NA"
        Info("HiLocus will use 'NA' as the prefix for all the ouput files")

    if not opt.outdir:
        opt.outdir = os.getcwd()
        opt.outdir = os.path.join(opt.outdir,'stat_out')
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)
    else:
        if not os.path.isdir(opt.outdir):
            os.mkdir(opt.outdir)
            
    if not opt.thr_intra:
        opt.thr_intra = "0.000001"
    else:
        if float(opt.thr_intra) <= 0 or float(opt.thr_intra) > 1:
            optparser.print_help()
            Info('ERROR: thr_intra value is not valid')
            sys.exit(1)
            
    if not opt.thr_inter:
        opt.thr_inter = "0.000001"
    else:
        if float(opt.thr_inter) <= 0 or float(opt.thr_inter) > 1:
            optparser.print_help()
            Info('ERROR: thr_inter value is not valid')
            sys.exit(1)

    Info("Argument List: ")
    Info("Case validpairs file = " + opt.case)
    Info("Control validpairs file = " + opt.control)
    Info("Binsize = " + opt.binsize)
    Info("Sample name = " + opt.name)
    Info("Output directory = " + opt.outdir)
    Info("Probability threshold for intrachromosomal translocations = " + opt.thr_intra)
    Info("Probability threshold for interchromosomal translocations = " + opt.thr_inter)


    return opt
