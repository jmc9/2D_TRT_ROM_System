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

    (proc_dir, plotdir, dset, dset_dat, dcmp_data, trunc_eps, trunc_opt) = Default_Inps() #load default inputs

    #if proc_dir, plotdir previously specified - use values
    if not proc_dir_ == '': proc_dir = proc_dir_
    if not plotdir_ == '': plotdir = plotdir_

    file = open(infile,'r') #open input file

    #reading input file
    for line in file:
        if (line.strip()): #disregarding blank lines
            input = line.split() #delimiting input line by space

            if input[0] ==   'proc_dir' : proc_dir  = input[1]

            elif input[0] == 'plotdir'  : plotdir   = input[1]

            elif input[0] == 'dset'     : dset      = input[1]

            elif input[0] == 'dset_dat' : dset_dat  = input[1]

            elif input[0] == 'dcmp_data': dcmp_data = input[1:]

            elif input[0] == 'trunc_eps': trunc_eps = float(input[1])

            elif input[0] == 'trunc_opt': trunc_opt = int(input[1])
    if ((trunc_opt in [1,2])and(trunc_eps == 0.)): trunc_eps = 1.

    Check_Inps(dset) #checking for invalid inputs

    dcmp_data = Data_Defaults(dcmp_data)

    return (proc_dir, plotdir, dset, dset_dat, dcmp_data, trunc_eps, trunc_opt)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Default_Inps():
    proc_dir  = 'processed_data'
    plotdir   = 'plots'
    dset      = 'dcmp_out.h5'
    dset_dat  = ''
    dcmp_data = 'QDf'
    trunc_eps = 0.
    trunc_opt = 3

    return (proc_dir, plotdir, dset, dset_dat, dcmp_data, trunc_eps, trunc_opt)

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
            try:
                stack_dat = getattr(ncdat, 'stack_dat')
                Dcmp_Data[p].opt = [1, stack_dat]

            except:
                Dcmp_Data[p].opt = [2]

        print('Collected {} Decomposition'.format(Dcmp_Data[p].name))
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

        if Dcmp_Data[i].opt[0] == 0:
            dnames = [Dcmp_Data[i].name]
        elif Dcmp_Data[i].opt[0] == 1:
            dnames = Dcmp_Data[i].opt[1].split(',')
        else:
            print('Input_Routines/Dat_Parms, unrecognized Dcmp_data.opt ({})'.format(Dcmp_Data[i].opt[0]))
            quit()

        try:
            ncdat = [dset[name] for name in dnames]
        except:
            print('no sir')
            Prob_Data[i].opt = [-1]
            continue

        Prob_Data[i].opt = Dcmp_Data[i].opt
        Prob_Data[i].dims = Dcmp_Data[i].dims

        dlen = 1
        for d in Prob_Data[i].dims[1:]: dlen *= d.len
        if Prob_Data[i].opt[0] == 0:
            dat = np.array(ncdat[0][:]).reshape(Prob_Data[i].dims[0].len, dlen)

        elif Prob_Data[i].opt[0] == 1:
            stack_dat = [np.array(ncd[:]) for ncd in ncdat]
            dat = tb.stack_data(stack_dat).reshape(Prob_Data[i].dims[0].len, dlen)

        Prob_Data[i].dat = dat

        if Prob_Data[i].opt[0] == 0:
            N_grids = getattr(ncdat[0], 'N_grids')
            Prob_Data[i].grids = []
            for j in range(N_grids):
                gname = getattr( ncdat[0], 'grid{}'.format(j) )
                grid = dset[gname][:]
                bnds = dset[gname].bnds
                Prob_Data[i].grids.append( Grid(bnds,grid) )

        elif Prob_Data[i].opt[0] == 1:
            grid1 = np.linspace(0., 1., Prob_Data[i].dims[ 0 ].len)
            bnds1 = [0., 1.]
            grid2 = np.linspace(0., 1., Prob_Data[i].dims[ len(Prob_Data[i].dims) - 1 ].len)
            bnds2 = [0., 1.]
            Prob_Data[i].grids = [ Grid(bnds1,grid1), Grid(bnds2,grid2) ]

        print('Collected {} Training Data'.format(Prob_Data[i].name))

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
