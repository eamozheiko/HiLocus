#!/usr/bin/env python


"""
ExoC - Translocation detection from Exome capture Hi-C

# short decription

@version: $1.0$
@author: Evgeniy Mozheiko
@contact: eamozheiko@gmail.com
"""

import sys, os, time
import re
from pkg_resources import resource_filename
from ExoC.ArgsValidator import *
from ExoC.corelib import *

def statrun(args):
    opts = opt_validate_stat(args)
    #command = "Rscript %s %s %s %s"%(args.sam, args.name, args.outdir)
    #print(command)
    #run_cmd(command)
    #print(os.path.abspath('.'))
    #print(opts)
    command = "Rscript %s %s"%(os.path.abspath('../ExoC/t.R'), os.path.abspath('..'))
    run_cmd(command)
    print(command)
    # Rscript
    
    #Info("Done! Find your CNV results from %s" %(args.outdir))
    return

def cnvrun(args):
    opts = opt_validate_cnv(args)
    
    ## calculate cnv pattern
    command = "Rscript cnv_pattern.R %s %s %s %s"%(os.path.abspath('../ExoC/t.R'), outputsubdir,matrixfile,tempbpoutputfile)
    print(command)
    run_cmd(command)
    
    Info("Done! Find your CNV results from %s; %s; and %s ;)"%(bicseqOut,outfig,outfig2))
    return

def transrun(args):
    opts = opt_validate_trans(args)
    
    ## inter

    command = "Rscript %s %s %s %s %s %s %s %s"%(os.path.abspath('../ExoC/inter_trans.R'), os.path.abspath('..'), opts.outdir, opts.case, opts.control, opts.binsize, opts.name, opts.thr_frame)
    print(command)
    #run_cmd(command)
    

    command = "Rscript %s %s %s %s %s %s %s %s"%(os.path.abspath('../ExoC/intra_trans.R'), os.path.abspath('..'), opts.outdir, opts.case, opts.control, opts.binsize, opts.name, opts.thr_frame)
    print(command)
    run_cmd(command)
    #run_cmd(command)
    
    ## intra
    #command = "Rscript inter_trans.R %s %s %s %s"%("../ExoC/intra_trans.R", opts.case, opts.control, opts.binsize, opts.name, opts.outdir)
    #print(command)
    #run_cmd(command)
    
    Info("Done! Find your translocation breakpoints file from %ssmall_inter_translocations_%s"%(opts.outdir, opts.name))
    return
