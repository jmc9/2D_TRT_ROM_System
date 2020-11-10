#!/usr/bin/env python3
#==================================================================================================================================#
#
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#netcdf tools
# from netCDF4 import Dataset as ncds

#==================================================================================================================================#
#
#==================================================================================================================================#
def input(infile):
    file = open(infile,'r')

    dsets=[]

    for line in file:
        if (line.strip()):
            input = line.split()

            if input[0] == 'ref_dataset': rset=input[1]
            elif input[0] == 'comp_datasets': dsets=input[1:]

    if dsets:
        if dsets == ['none']: dsets = [rset]
        else: dsets.insert(0,rset)
    else: dsets = [rset]

    file.close()

    return (dsets)

#==================================================================================================================================#
#
#==================================================================================================================================#
def inp_flags(dset):
    inp_flags = [dset.getncattr('I_out'),\
    dset.getncattr('D_out'),\
    dset.getncattr('fg_out'),\
    dset.getncattr('old_parms_out'),\
    dset.getncattr('HO_Eg_out'),\
    dset.getncattr('HO_E_out'),\
    dset.getncattr('HO_Fg_out'),\
    dset.getncattr('HO_F_out'),\
    dset.getncattr('Eg_out'),\
    dset.getncattr('MGQD_E_out'),\
    dset.getncattr('Fg_out'),\
    dset.getncattr('MGQD_F_out'),\
    dset.getncattr('E_out'),\
    dset.getncattr('F_out'),\
    dset.getncattr('its_out'),\
    dset.getncattr('conv_out'),\
    dset.getncattr('Src_out')]

    return inp_flags

#==================================================================================================================================#
#
#==================================================================================================================================#
def domain_parms(dset):
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

    (xp,yp) = cell_coords(Delx,Dely)
    tp = []
    for i in range(N_t): tp.append((i+1)*Delt)

    return (xp,yp,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x)

#==================================================================================================================================#
#
#==================================================================================================================================#
def cell_coords(Delx,Dely):
    xp = [0.]
    for i in range(len(Delx)): xp.append(sum(Delx[:i]) + Delx[i])
    yp = [0.]
    for i in range(len(Dely)): yp.append(sum(Dely[:i]) + Dely[i])

    return (xp,yp)
