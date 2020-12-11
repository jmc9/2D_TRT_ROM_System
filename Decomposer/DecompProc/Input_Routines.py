#!/usr/bin/env python3
#==================================================================================================================================#
#
# Input_Routines.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#local tools
import ToolBox as tb

#==================================================================================================================================#
#
#==================================================================================================================================#
def Input(infile,proc_dir_,plotdir_):

    tb.FileCheck(infile) #check if input file exists

    (proc_dir,plotdir,dset) = Default_Inps() #load default inputs

    #if proc_dir, plotdir previously specified - use values
    if not proc_dir_ == '': proc_dir = proc_dir_
    if not plotdir_ == '': plotdir = plotdir_

    file = open(infile,'r') #open input file

    #reading input file
    for line in file:
        if (line.strip()): #disregarding blank lines
            input = line.split() #delimiting input line by space

            if input[0] == 'proc_dir': proc_dir = input[1]

            elif input[0] == 'plotdir': plotdir = input[1]

            elif input[0] == 'dset': dset = input[1]

    Check_Inps(dset) #checking for invalid inputs

    return (proc_dir,plotdir,dset)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Default_Inps():
    proc_dir = 'processed_data'
    plotdir  = 'plots'
    dset     = 'dcmp_out.h5'

    return (proc_dir,plotdir,dset)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Check_Inps(dset):
    tb.FileCheck(dset)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Dcmp_Parms(dset):
    dcmp_type = getattr(dset,'dcmp_type')
    dcmp_data = getattr(dset,'dcmp_data')

    (xp,yp,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x) = Domain_Parms(dset)

    if dcmp_data == 'QDf':
        BClen = dset.dimensions['BClen'].size
        rank_BC = dset.dimensions['rank_BC'].size
    else:
        BClen = 0.
        rank_BC = 0.

    clen_avg = dset.dimensions['clen_avg'].size
    clen_edgV = dset.dimensions['clen_edgV'].size
    clen_edgH = dset.dimensions['clen_edgH'].size

    rank_avg = dset.dimensions['rank_avg'].size
    rank_edgV = dset.dimensions['rank_edgV'].size
    rank_edgH = dset.dimensions['rank_edgH'].size

    return (dcmp_type,dcmp_data,xp,yp,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x,
    BClen,clen_avg,clen_edgV,clen_edgH,rank_BC,rank_avg,rank_edgV,rank_edgH)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Domain_Parms(dset):
    Delx = dset['Delx'][:]
    Dely = dset['Dely'][:]
    Delt = dset['Delt']

    N_t = dset.dimensions['N_t'].size
    N_g = dset.dimensions['N_g'].size
    N_y = dset.dimensions['N_y'].size
    N_x = dset.dimensions['N_x'].size

    A=[]
    for j in range(len(Dely)):
        for i in range(len(Delx)):
            A.append(Dely[j]*Delx[i])

    Delt = 2e-2

    (xp,yp) = tb.Cell_Coords(Delx,Dely)
    tp = []
    for i in range(N_t): tp.append((i+1)*Delt)

    return (xp,yp,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Name_Gen(dcmp_data):
    if dcmp_data == 'QDf':
        names = ['BCg','fg_avg_xx','fg_edgV_xx','fg_avg_yy','fg_edgH_yy','fg_edgV_xy','fg_edgH_xy']

    return names
