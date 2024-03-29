#!/usr/bin/env python


"""
HiLocus - Translocation detection from capture Hi-C

# short decription

@version: $1.0$
@author: Evgeniy Mozheiko
@contact: eamozheiko@gmail.com
"""

import sys, os, time
import re
from pkg_resources import resource_filename
from HiLocus.ArgsValidator import *
from HiLocus.corelib import *

def statrun(args):
    opts = opt_validate_stat(args)
    
    path_to_stat_script = resource_filename('HiLocus', 'scripts/sam_to_valid_pairs.sh')
    command = "%s %s"%("chmod u+x", path_to_stat_script)
    run_cmd(command)
    
    command = "%s %s %s %s"%(path_to_stat_script, os.path.abspath(opts.sam), os.path.abspath(opts.outdir), opts.name)
    print(command)
    run_cmd(command)
    
    Info("Done! Find your stat results from %s" %(opts.outdir))
    return


def transrun(args):
    opts = opt_validate_trans(args)
    
    ## prepare input
    input_format_case = os.path.splitext(opts.case)[-1]
    input_format_control = os.path.splitext(opts.control)[-1]
    
    ## intra
    path_to_intra_trans_script = resource_filename('HiLocus', 'scripts/intra_trans.R')
    command = "Rscript %s %s %s %s %s %s %s %s %s"%(path_to_intra_trans_script, os.path.abspath(opts.outdir), os.path.abspath(opts.case), os.path.abspath(opts.control), opts.binsize, opts.name, opts.thr_intra, input_format_case, input_format_control)
    print(command)
    run_cmd(command)
    
    ## inter
    path_to_inter_trans_script = resource_filename('HiLocus', 'scripts/inter_trans.R')
    command = "Rscript %s %s %s %s %s %s %s %s %s"%(path_to_inter_trans_script, os.path.abspath(opts.outdir), os.path.abspath(opts.case), os.path.abspath(opts.control), opts.binsize, opts.name, opts.thr_inter, input_format_case, input_format_control)
    print(command)
    run_cmd(command)
    
    Info("Done! Find your translocation results from %s"%(opts.outdir))
    return
