#!/usr/bin/env python3
#==================================================================================================================================#
#
# Decomp_Processor.py is the main python script to process the outputs of Decomposer.c
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import getopt, sys #handles user inputs from command line
import numpy as np
import os

#netcdf tools
from netCDF4 import Dataset as ncds

#local tools
import ToolBox as tb
import Input_Routines as inp
import DMD_Routines as dmdr
import Process_Routines as pr

#==================================================================================================================================#
#
#==================================================================================================================================#
def Dcmp_Proc(infile,proc_dir,plotdir):

    tb.BigTitle() #output code title to command line

    #===========================================================================#
    #                                                                           #
    #                              READING INPUTS                               #
    #                                                                           #
    #===========================================================================#
    (proc_dir,plotdir,dset,dset_dat,data_names) = inp.Input(infile,proc_dir,plotdir)
    tb.dirset(proc_dir)
    plotdir = proc_dir+'/'+plotdir
    tb.dirset(plotdir)

    #===========================================================================#
    #                                                                           #
    #                       COLLECTING DECOMPOSITION DATA                       #
    #                                                                           #
    #===========================================================================#
    dset = ncds(dset,'r') #opening specified dataset
    Dcmp_Data = inp.Dcmp_Parms(dset,data_names)

    #===========================================================================#
    #                                                                           #
    #                         COLLECTING REFERENCE DATA                         #
    #                                                                           #
    #===========================================================================#
    if dset_dat != '':
        dset_dat = ncds(dset_dat,'r')
        Prob_Data = inp.Dat_Parms(dset_dat,Dcmp_Data)

        if (len(Prob_Data) != len(Dcmp_Data)):
            print('len(Prob_Data) != len(Dcmp_Data)')
            quit()
    else:
        Prob_Data = []

    #===========================================================================#
    #                                                                           #
    #                              PROCESSING DATA                              #
    #                                                                           #
    #===========================================================================#
    plt_modes = [0,1]
    recon = pr.Process_Data(dset,Dcmp_Data,Prob_Data,plt_modes,plotdir)
    err = pr.Error_Calc(recon,Prob_Data)
    pr.Output(recon, err, Prob_Data)
    os.system('cp proc_summary.h5 {}'.format(plotdir))

    dset.close() #closing dataset file
    print('success!')

#==================================================================================================================================#
#
#==================================================================================================================================#
def Options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['help','input=','output='])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        quit()

    (infile,proc_dir,plotdir) = Default_Opts()

    for o, a in opts:
        if o in ['-h','--help']:
            usage()
            quit()
        elif o in ['-i','--input']:
            infile = a
        elif o in ['-o','--output']:
            proc_dir = a

    return (infile,proc_dir,plotdir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Default_Opts():
    infile   = 'input.inp'
    proc_dir = ''
    plotdir  = ''

    return (infile,proc_dir,plotdir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def usage():
    tb.BigTitle()
    print()
    print('- Joseph M. Coale, jmcoale@ncsu.edu, josephcoale912@gmail.com ')
    print()
    print('This python3 code is used to process 2D TRT data that has been decomposed with Decomposer.c (Author: Joseph M. Coale)')
    print('Below are all available options and their usages')
    print()
    print('Options:')
    print('-h [help]: prints this message')
    print('-i [input]: used to specify an input file')
    print('-o [output]: used to specify the output directory')
    print()

#==================================================================================================================================#
# MAIN PROGRAM
# This portion of the code is run when this file is the main program
# i.e. when the call is made './TRT_processor.py'
#==================================================================================================================================#
if __name__ == '__main__':
    (infile,proc_dir,plotdir) = Options()
    Dcmp_Proc(infile,proc_dir,plotdir)
