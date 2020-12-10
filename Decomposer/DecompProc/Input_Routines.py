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

    return (dcmp_type,dcmp_data,xp,yp,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x)

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

    (xp,yp) = Cell_Coords(Delx,Dely)
    tp = []
    for i in range(N_t): tp.append((i+1)*Delt)

    return (xp,yp,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Cell_Coords(Delx,Dely):
    xp = [0.]
    for i in range(len(Delx)): xp.append(sum(Delx[:i]) + Delx[i])
    yp = [0.]
    for i in range(len(Dely)): yp.append(sum(Dely[:i]) + Dely[i])

    return (xp,yp)
