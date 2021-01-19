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
import numpy as np

from Classes import Data
from Classes import DMD
from Classes import Grid
from Classes import Dimension

#==================================================================================================================================#
#
#==================================================================================================================================#
def Input(infile,proc_dir_,plotdir_):

    tb.FileCheck(infile) #check if input file exists

    (proc_dir,plotdir,dset,dset_dat,dcmp_data) = Default_Inps() #load default inputs

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

            elif input[0] == 'dset_dat': dset_dat = input[1]

            elif input[0] == 'dcmp_data': dcmp_data = input[1:]

    Check_Inps(dset) #checking for invalid inputs

    dcmp_data = Data_Defaults(dcmp_data)

    return (proc_dir,plotdir,dset,dset_dat,dcmp_data)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Default_Inps():
    proc_dir  = 'processed_data'
    plotdir   = 'plots'
    dset      = 'dcmp_out.h5'
    dset_dat  = ''
    dcmp_data = 'QDf'

    return (proc_dir,plotdir,dset,dset_dat,dcmp_data)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Check_Inps(dset):
    tb.FileCheck(dset)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Dcmp_Parms(dset,data_names):
    dcmp_type = getattr(dset,'dcmp_type') #type of decomposition found in the datafile

    Dcmp_Data = [] #array to hold all decomposition data
    N_data = len(data_names) #number of input datasets

    #collecting information on each dataset
    p = 0
    for i in range(N_data):
        #checking if dataset exists in datafile
        try:
            ncdat = dset[data_names[i]]

        except:
                print("Error! {} either does not exist, or has no attribute block".format(data_names[i]))
                print("- Henceforth ignoring {} in visualization".format(data_names[i]))
                continue

        Dcmp_Data.append(Data(data_names[i])) #allocating element of Dcmp_Data as Data object

        Dcmp_Data[p].typeset(dcmp_type)

        Dcmp_Data[p].opt = [0]

        #collecting dimensions of dataset
        N_dims = getattr(ncdat, 'N_dims')
        Dcmp_Data[p].dims = []
        for j in range(N_dims):
            Dcmp_Data[p].dims.append( Dimension( name = getattr( ncdat, 'dim{}'.format(j) ) ) )
            Dcmp_Data[p].dims[j].len = dset.dimensions[Dcmp_Data[p].dims[j].name].size

        #collecting grids on which the dataset resides
        try:
            N_grids = getattr(ncdat, 'N_grids')

            Dcmp_Data[p].grids = []
            for j in range(N_grids):
                gname = getattr( ncdat, 'grid{}'.format(j) )
                grid = dset[gname][:]
                bnds = dset[gname].bnds
                Dcmp_Data[p].grids.append( Grid(bnds,grid) )

        except:
            Dcmp_Data[p].opt = [1]

        p+=1

    if not Dcmp_Data:
        print("Error! None of the data input for processing was found, terminating program")
        quit()

    return Dcmp_Data

#==================================================================================================================================#
#
#==================================================================================================================================#
def Dat_Parms(dset,Dcmp_Data):
    Prob_Data = [] #array to hold all data from problem
    N_data = len(Dcmp_Data) #number of input datasets

    #collecting information on each dataset
    for i in range(N_data):

        Prob_Data.append(Data(Dcmp_Data[i].name)) #allocating element of Dcmp_Data as Data object

        #checking if dataset exists in datafile
        try:
            ncdat = dset[Dcmp_Data[i].name]

        except:
            Prob_Data[i].opt = [-1]
            continue

        Prob_Data[i].opt = [0]

        #collecting dimensions of dataset
        dnames = ncdat.dimensions
        N_dims = len(dnames)
        Prob_Data[i].dims = []
        for j in range(N_dims):
            Prob_Data[i].dims.append( Dimension( name = dnames[j], len = dset.dimensions[dnames[j]].size ) )

        #
        dat = ncdat[:]
        dlen = 1
        for j in range(1,N_dims):
            dlen = dlen * Prob_Data[i].dims[j].len
        Prob_Data[i].dat = np.zeros([Prob_Data[i].dims[0].len, dlen])
        if N_dims == 4:
            for t in range(Prob_Data[i].dims[0].len):
                Prob_Data[i].dat[t] = tb.Flatten3D(dat[t])
        elif N_dims == 3:
            for t in range(Prob_Data[i].dims[0].len):
                Prob_Data[i].dat[t] = tb.Flatten2D(dat[t])
        elif N_dims == 2:
            for t in range(Prob_Data[i].dims[0].len):
                Prob_Data[i].dat[t] = dat[t]
        elif N_dims > 4:
            print('cannot flatted data of dimensionality greater than 4 right now!')
            quit()

        #collecting grids on which the dataset resides
        try:
            N_grids = getattr(ncdat, 'N_grids')

        except:
            Prob_Data[i].opt = [1]

        if Prob_Data[i].opt == [0]:
            Prob_Data[i].grids = []
            for j in range(N_grids):
                gname = getattr( ncdat, 'grid{}'.format(j) )
                grid = dset[gname][:]
                bnds = dset[gname].bnds
                Prob_Data[i].grids.append( Grid(bnds,grid) )

    if not Prob_Data:
        print("Error! None of the data input for processing was found, terminating program")
        quit()

    return Prob_Data


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
def Data_Defaults(dcmp_data):
    data_names = []

    for i in range(len(dcmp_data)):
        if dcmp_data[i] == 'QDf':
            data_names.append('BCg')
            data_names.append('fg_avg_xx')
            data_names.append('fg_edgV_xx')
            data_names.append('fg_avg_yy')
            data_names.append('fg_edgH_yy')
            data_names.append('fg_edgV_xy')
            data_names.append('fg_edgH_xy')
        elif dcmp_data[i] == "I":
            data_names.append("I_avg")
            data_names.append("I_edgV")
            data_names.append("I_edgH")
        else:
            data_names.append(dcmp_data[i])

    return data_names
