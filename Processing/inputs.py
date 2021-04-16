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
    plt_times = []
    Tbound = [0.,0.]
    Ebound = [0.,0.]
    fg_avg_xx_bnd = []
    ntrend_tp = []
    dsnames = []
    fignames = []
    trend_names = []
    cs_times = []
    csx_times = []
    csy_times = []
    csx = [0, 1, 0]
    csy = [0, 1, 1]
    switch_ranks = 1
    switch_norms = 0
    switch_cs = 0
    switch_sol = 1

    for line in file:
        if (line.strip()):
            input = line.split()

            if input[0] == 'ref_dataset': rset=input[1]

            elif input[0] == 'comp_datasets': dsets=input[1:]

            elif input[0] == 'dsnames': dsnames=input[1:]

            elif input[0] == 'fignames': fignames=input[1:]

            elif input[0] == 'plt_times_init':
                for inp in input[1:]: plt_times.append(float(inp))

            elif input[0] == 'trend_names': trend_names=input[1:]

            elif input[0] == 'plt_tfreq': plt_tfreq = int(input[1])

            elif input[0] == 'nrm_trend_times':
                for i in range(len(input)-1): ntrend_tp.append(int(input[i+1]))

            elif input[0] == 'cs_times':
                for i in range(len(input)-1): cs_times.append(float(input[i+1]))
            elif input[0] == 'csx':
                csx[0] = int(input[1])
                csx[1] = int(input[2])
                csx[2] = int(input[3])
                for i in range(len(input)-4): csx_times.append(float(input[i+4]))
            elif input[0] == 'csy':
                csy[0] = int(input[1])
                csy[1] = int(input[2])
                csy[2] = int(input[3])
                for i in range(len(input)-4): csy_times.append(float(input[i+4]))

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

            elif input[0] == 'switch_ranks': switch_ranks = int(input[1])
            elif input[0] == 'switch_norms': switch_norms = int(input[1])
            elif input[0] == 'switch_cs':    switch_cs    = int(input[1])
            elif input[0] == 'switch_sol':   switch_sol   = int(input[1])
            elif input[0] == 'switch_errplots': switch_errplots = int(input[1])

    if dsets:
        if dsets == ['none']: dsets = [rset]
        else: dsets.insert(0,rset)
    else: dsets = [rset]

    if (len(csx_times)>0 or len(csx_times)>0):
        if len(csy_times)==0: csy_times = [t for t in csx_times]
        if len(csx_times)==0: csx_times = [t for t in csy_times]
        cs_times = [t for t in csx_times]
        for t in csy_times:
            if not t in cs_times:
                for i in range(len(cs_times)):
                    if t>cs_times[i]:
                        cs_times.insert(i+1,t)
                        break
    else:
        csy_times = [t for t in cs_times]
        csx_times = [t for t in cs_times]

    if dsnames:
        if len(dsets)!=1 and len(dsnames)!=len(dsets)-1:
            print('ERROR! not every comp_dataset has been asigned a "dsname" (or too many dsnames have been specified)')
            quit()
        else:
            for i in range(len(dsnames)):
                dsnames[i] = dsnames[i].replace('~',' ')
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

    if fignames:
        if len(dsets)!=1 and len(trend_names)!=len(dsets)-1:
            print('ERROR! not every comp_dataset has been asigned a "fignames" (or too many fignames have been specified)')
            quit()
    else:
        if len(dsets)!=1:
            print('ERROR! "fignames" not specified in input')
            quit()

    # if ntrend_tp: switch_norms=1

    file.close()

    return (dsets,dsnames,trend_names,fignames,plt_times,plt_tfreq,ntrend_tp,Tbound,Ebound,fg_avg_xx_bnd,
    cs_times,csx_times,csy_times,csx,csy,switch_ranks,switch_norms,switch_cs,switch_sol,switch_errplots)

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

    (xp,yp,xp_c,yp_c) = cell_coords(Delx,Dely)
    # tp = []
    # for i in range(N_t): tp.append((i+1)*Delt)
    tp=[t*10. for t in dset['tpts'][:]]

    return (xp,yp,xp_c,yp_c,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x)

#==================================================================================================================================#
#
#==================================================================================================================================#
def cell_coords(Delx,Dely):
    xp = [0.]
    for i in range(len(Delx)): xp.append(sum(Delx[:i]) + Delx[i])
    yp = [0.]
    for i in range(len(Dely)): yp.append(sum(Dely[:i]) + Dely[i])

    xp_c = []
    for i in range(len(Delx)): xp_c.append(sum(Delx[:i]) + Delx[i]/2)
    yp_c = []
    for i in range(len(Dely)): yp_c.append(sum(Dely[:i]) + Dely[i]/2)

    return (xp,yp,xp_c,yp_c)

#==================================================================================================================================#
#
#==================================================================================================================================#
def sample_infile():
    file = open('Sample.inp','w')

    file.write("* General notes:\n")
    file.write("* -> comments follow the '*' character\n")
    file.write("* -> all files passed to the program must have been generated with the '2D TRT ROM code'\n")
    file.write("* -> an example of each input argument is given below\n")
    file.write("* -> before each example is a brief description of the input argument\n\n")

    file.write("* ref_dataset:\n")
    file.write("* Arguments: 1\n")
    file.write("* Function:  specifies the reference file (i.e. reference or full-order solution)\n")
    file.write("* Required?: REQUIRED\n")
    file.write("ref_dataset fom.h5\n\n")

    file.write("* comp_datasets:\n")
    file.write("* Arguments: arbitrary number\n")
    file.write("* Function:  specifies the 'comparison' files (i.e. approximate or reduced-order solutions)\n")
    file.write("* Required?: REQUIRED\n")
    file.write("comp_datasets diff.h5 fld.h5 p13.h5 p1.h5\n\n")

    file.write("* dsnames:\n")
    file.write("* Arguments: same number of arguments as comp_datasets\n")
    file.write("* Function:  specifies the names of each data file given in comp_datasets, used in plots\n")
    file.write("* Required?: REQUIRED\n")
    file.write("dsnames Diffusion FLD $P_{1/3}$ $P_1$\n\n")

    file.write("* trend_names:\n")
    file.write("* Arguments: same number of arguments as comp_datasets\n")
    file.write("* Function:  similar to dsnames, but for specifying parameters in the legend of trend plots\n")
    file.write("* Required?: REQUIRED\n")
    file.write("trend_names Diffusion FLD $P_{1/3}$ $P_1$\n\n")

    file.write("* fignames:\n")
    file.write("* Arguments: one more argument than comp_datasets\n")
    file.write("* Function:  specifies abridged names of each data file given in comp_datasets (and the reference file), used to save figures with\n")
    file.write("* Required?: REQUIRED\n")
    file.write("fignames fom diff fld p13 p1\n\n")

    file.write("* plt_times_init:\n")
    file.write("* Arguments: 1\n")
    file.write("* Function:  specifies the time step to begin plotting solutions\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("plt_times_init 0\n\n")

    file.write("* plt_tfreq:\n")
    file.write("* Arguments: 1\n")
    file.write("* Function:  specifies time-step frequency of solution plotting (after initial time)\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("plt_tfreq 10\n\n")

    file.write("* Tbound:\n")
    file.write("* Arguments: 2\n")
    file.write("* Function:  specifies lower, upper bounds for plots of temperature (in eV)\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("Tbound 0. 1000.\n\n")

    file.write("* Ebound:\n")
    file.write("* Arguments: 2\n")
    file.write("* Function:  specifies lower, upper bounds for plots of radiation energy density (in ergs X 10^{13})\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("Ebound 0. 10.\n\n")

    file.write("* nrm_trend_times:\n")
    file.write("* Arguments: arbirtrary\n")
    file.write("* Function:  specifies time steps to include in trend plots\n")
    file.write("* Required?: OPTIONAL (but required if using trend plots)\n")
    file.write("nrm_trend_times 0 49 99 149 199 249 299\n\n")

    file.write("* csx:\n")
    file.write("* Arguments: 3 + arbitrary\n")
    file.write("* Function (first 3 args): flag (0 or 1) for plotting cross-sections along x-axis on the (bottom, middle, top) of the domain\n")
    file.write("* Function (4+ args): times (in ns) to plot x- cross sections at\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("csx 0 1 0 .01 .05 .1 .2\n\n")

    file.write("* csy:\n")
    file.write("* Arguments: 3 + arbitrary\n")
    file.write("* Function (first 3 args): flag (0 or 1) for plotting cross-sections along y-axis on the (left, middle, right) of the domain\n")
    file.write("* Function (4+ args): times (in ns) to plot y- cross sections at\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("csy 0 1 1 .05 .1 .2\n\n")

    file.write("* switch_ranks:\n")
    file.write("* Arguments: 1\n")
    file.write("* Function:  flag (0 or 1) to turn off/on respectively plotting of ranks\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("switch_ranks 0\n\n")

    file.write("* switch_norms:\n")
    file.write("* Arguments: 1\n")
    file.write("* Function:  flag (0 or 1) to turn off/on respectively plotting of error norms\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("switch_norms 1\n\n")

    file.write("* switch_cs:\n")
    file.write("* Arguments: 1\n")
    file.write("* Function:  flag (0 or 1) to turn off/on respectively plotting of cross_sections\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("switch_cs 1\n\n")

    file.write("* switch_sol:\n")
    file.write("* Arguments: 1\n")
    file.write("* Function:  flag (0 or 1) to turn off/on respectively plotting of solutions\n")
    file.write("* Required?: OPTIONAL\n")
    file.write("switch_sol 0\n\n")

    file.close()
