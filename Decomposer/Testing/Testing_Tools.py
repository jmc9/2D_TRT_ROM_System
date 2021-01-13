#==================================================================================================================================#
#
# Testing_Tools.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#
from Test_Functions import expf
from IO_Functions import output_f as outf
from DecompSrc_Handling import exec_decompsrc as ds_exec
from DecompProc_Handling import exec_decompproc as dp_exec

#==================================================================================================================================#
#
#==================================================================================================================================#
def Run_Test(exec_dir, decomp_dir, proc_dir, drop):

    #
    f, fgrids = expf([.1,.2])
    # print(f)

    #
    outfile = 'test.h5'
    outf(f, fgrids, outfile)

    #
    decomp_infile = 'input.inp'
    decomp_outfile = 'test_dmd.h5'
    decomp_logfile = 'DecompLog.log'

    proc_infile = 'proc_input.inp'
    proc_outdir = 'processed_test'
    proc_logfile = 'ProgLog.log'

    ds_exec(decomp_infile, decomp_outfile, outfile, decomp_logfile, drop, exec_dir, decomp_dir)
    dp_exec(proc_infile, proc_outdir, decomp_outfile, outfile, proc_logfile, drop, exec_dir, proc_dir)

    return
