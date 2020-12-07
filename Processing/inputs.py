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

    dsets = []
    tp_plt = []
    Tbound = [0.,0.]
    Ebound = [0.,0.]
    fg_avg_xx_bnd = []
    ntrend_tp = []
    dsnames = []
    trend_names = []

    for line in file:
        if (line.strip()):
            input = line.split()

            if input[0] == 'ref_dataset': rset=input[1]

            elif input[0] == 'comp_datasets': dsets=input[1:]

            elif input[0] == 'dsnames': dsnames=input[1:]

            elif input[0] == 'plt_times_init':
                for inp in input[1:]: tp_plt.append(int(inp))

            elif input[0] == 'trend_names': trend_names=input[1:]

            elif input[0] == 'plt_tfreq': plt_tfreq = int(input[1])

            elif input[0] == 'nrm_trend_times':
                for i in range(len(input)-1): ntrend_tp.append(int(input[i+1]))

            elif input[0] == 'Tbound':
                Tbound[0] = float(input[1])
                Tbound[1] = float(input[2])

            elif input[0] == 'Ebound':
                Ebound[0] = float(input[1])
                Ebound[1] = float(input[2])

            elif input[0] == 'fg_avg_xx_bnd':
                fg_avg_xx_bnd.append(int(input[1]))
                fg_avg_xx_bnd.append(float(input[2]))
                fg_avg_xx_bnd.append(float(input[3]))

    if dsets:
        if dsets == ['none']: dsets = [rset]
        else: dsets.insert(0,rset)
    else: dsets = [rset]

    if dsnames:
        if len(dsets)!=1 and len(dsnames)!=len(dsets)-1:
            print('ERROR! not every comp_dataset has been asigned a "dsname" (or too many dsnames have been specified)')
            quit()
        else:
            for i in range(len(dsnames)):
                dsnames[i] = dsnames[i].replace('_',' ')
    else:
        if len(dsets)!=1:
            print('ERROR! "dsnames" not specified in input')
            quit()

    if trend_names:
        if len(dsets)!=1 and len(trend_names)!=len(dsets)-1:
            print('ERROR! not every comp_dataset has been asigned a "trend_name" (or too many dsnames have been specified)')
            quit()
        else:
            for i in range(len(trend_names)):
                trend_names[i] = trend_names[i].replace('_',' ')
    else:
        if len(dsets)!=1:
            print('ERROR! "trend_names" not specified in input')
            quit()

    file.close()

    return (dsets,dsnames,trend_names,tp_plt,plt_tfreq,ntrend_tp,Tbound,Ebound,fg_avg_xx_bnd)

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
