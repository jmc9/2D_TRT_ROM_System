#!/usr/bin/env python3
#==================================================================================================================================#
#
# Process_Routines.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import math
import cmath
import numpy as np

#netcdf tools
from netCDF4 import Dataset as ncds

#local tools
import DMD_Routines as dmdr
import ToolBox as tb

#==================================================================================================================================#
#
#==================================================================================================================================#
def Process_Data(dset,Dcmp_Data,Prob_Data,plt_modes,plotdir):
    recon = []

    N_data = len(Dcmp_Data)
    for i in range(N_data):
        dcmp = Dcmp_Data[i]
        recon.append([])

        dir = '{}/{}'.format(plotdir,dcmp.name)
        tb.dirset(dir)

        if (dcmp.type in ['DMD','DMDg']):
            (dcmp.rank, dcmp.clen) = dmdr.Read_Dims(dset,dcmp.name)
            dcmp.dat = dmdr.Read_DMD(dset,dcmp.name)

            dcmp = dmdr.DMD_Sort(dcmp)

            if dcmp.opt[0] == 0:
                dmdr.Plot_DMD(dcmp,dir,plt_modes)

                if len(Prob_Data) != 0:
                    if Prob_Data[i].opt[0] == 0:
                        coef = dmdr.Coef_Calc(dcmp,Prob_Data[i].dat[0])

                        len1 = len(Prob_Data[i].dat)
                        len2 = len(Prob_Data[i].dat[0])
                        edat = np.zeros([len1, len2],dtype=np.complex_)
                        for j in range(len1):
                            edat[j] = dmdr.Expand(dcmp, coef, j)
                        recon[i] = edat

            else:
                dmdr.Plot_DMD(dcmp,dir,plt_modes,evecs=False)

    return recon

#==================================================================================================================================#
#
#==================================================================================================================================#
def Error_Calc(recon, Prob_Data):
    err = []

    #if there is no data, abort - otherwise save number of data
    if len(Prob_Data) != 0:
        N_data = len(Prob_Data)
    else:
        return err

    for j in range(N_data):
        tlen = len(Prob_Data[j].dat)
        clen = len(Prob_Data[j].dat[0])
        err.append([]) #appending blank array first -> stays consistent with N_data and allows easy flag for uncalculated errors

        #straightforward error calculations for opt[0] = 0 (standard data array)
        if Prob_Data[j].opt[0] == 0:
            err_ = np.zeros((tlen,clen), dtype='complex') #temporary complex array to hold j-th data errors
            for t in range(tlen):
                for i in range(clen):
                    err_[t][i] = Prob_Data[j].dat[t][i] - recon[j][t][i]
            err[j] = err_ #now putting j-th calculated errors into err array for output

    return err

#==================================================================================================================================#
#
#==================================================================================================================================#
def Output(recon, err, Prob_Data):
    dset = ncds('proc_summary.h5','w') #opening specified dataset

    #if there is no data, abort - otherwise save number of data
    if len(Prob_Data) != 0:
        N_data = len(Prob_Data)
    else:
        return

    #need np arrays to isolate real/imag components later
    err = np.array(err)
    recon = np.array(recon)

    dims = []
    d = 0
    for i in range(N_data):
        dlocs = []
        if len(recon[i]) != 0:
            ndims = len(Prob_Data[i].dims)
            ngrids = len(Prob_Data[i].grids)
            dname = Prob_Data[i].name

            #--------------------------------------------------
            # Collecting all unique dimensions
            #--------------------------------------------------
            for j in range(ndims):
                if Prob_Data[i].dims[j] in dims:
                    dlocs.append(j)
                else:
                    dims.append(Prob_Data[i].dims[j])
                    nc_dim = dset.createDimension("dim_{}".format(d), Prob_Data[i].dims[j])
                    d += 1
                    dlocs.append(j)

            #--------------------------------------------------
            # Creating err, recon variables (in 2D)
            #--------------------------------------------------
            if ndims == 2:
                #-------Create err, recon nc variables---------
                nc_err_r = dset.createVariable("{}_Err_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1])))
                nc_err_i = dset.createVariable("{}_Err_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1])))

                nc_recon_r = dset.createVariable("{}_Recon_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1])))
                nc_recon_i = dset.createVariable("{}_Recon_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1])))

                nc_2err_r = dset.createVariable("{}_Err_2_Real".format(dname), "f8", ("dim_{}".format(dlocs[0])))
                nc_2err_i = dset.createVariable("{}_Err_2_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0])))

                #----Calculate 2-norm of err (with 2 grids)----
                # 2 grids + 2 dimensions means data is structured like f[time][space(x)]
                if ngrids == 2:
                    err2 = np.zeros(dims[dlocs[0]], dtype='complex')
                    for j in range(dims[dlocs[0]]):
                        for k in range(dims[dlocs[1]]):
                            err2[j] = err2[j] + err[i][j][k]**2
                        err2[j] = cmath.sqrt(err2[j])

                #---skipping 2-norm if there are not 2 grids---
                # with 2 dimensions, anything other than 2 grids either doesn't make sense or doesn't allow -
                # - for a physically interpretable 2-norm
                else:
                    err2 = []



            #--------------------------------------------------
            # Creating err, recon variables (in 3D)
            #--------------------------------------------------
            elif ndims == 3:
                #-------Create err, recon nc variables---------
                nc_err_r = dset.createVariable("{}_Err_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2])))
                nc_err_i = dset.createVariable("{}_Err_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2])))

                nc_recon_r = dset.createVariable("{}_Recon_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2])))
                nc_recon_i = dset.createVariable("{}_Recon_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2])))

                np.reshape(err[i], (dims[dlocs[0]], dims[dlocs[1]], dims[dlocs[2]]))
                np.reshape(recon[i], (dims[dlocs[0]], dims[dlocs[1]], dims[dlocs[2]]))

                #----Calculate 2-norm of err (with 3 grids)----
                # 3 grids + 3 dimensions means data is structured like f[time][space(y)][space(x)]
                if ngrids == 3:
                    nc_2err_r = dset.createVariable("{}_Err_2_Real".format(dname), "f8", ("dim_{}".format(dlocs[0])))
                    nc_2err_i = dset.createVariable("{}_Err_2_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0])))

                    err2 = np.zeros(dims[dlocs[0]], dtype='complex')
                    for j in range(dims[dlocs[0]]):
                        for y in range(dims[dlocs[1]]):
                            for x in range(dims[dlocs[2]]):
                                err2[j] = err2[j] + err[i][j][y][x]**2
                        err2[j] = cmath.sqrt(err2[j])

                #----Calculate 2-norm of err (with 2 grids)----
                # 2 grids + 3 dimensions means data is structured like f[time][group][space(x)]
                elif ngrids == 2:
                    nc_2err_r = dset.createVariable("{}_Err_2_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))
                    nc_2err_i = dset.createVariable("{}_Err_2_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))

                    err2 = np.zeros((dims[dlocs[0]], dims[dlocs[1]]), dtype='complex')
                    for j in range(dims[dlocs[0]]):
                        for k in range(dims[dlocs[1]]):
                            for x in range(dims[dlocs[2]]):
                                err2[j][k] = err2[j][k] + err[i][j][k][x]**2
                            err2[j][k] = cmath.sqrt(err2[j][k])

                #---skipping 2-norm if there are not 2 or 3 grids---
                # with 3 dimensions, anything other than 2 or 3 grids either doesn't make sense or doesn't allow -
                # - for a physically interpretable 2-norm
                else:
                    err2 = []

            #--------------------------------------------------
            # Creating err, recon variables (in 4D)
            #--------------------------------------------------
            elif ndims == 4:
                #-------Create err, recon nc variables---------
                nc_err_r = dset.createVariable("{}_Err_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2]), "dim_{}".format(dlocs[3])))
                nc_err_i = dset.createVariable("{}_Err_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2]), "dim_{}".format(dlocs[3])))

                nc_recon_r = dset.createVariable("{}_Recon_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2]), "dim_{}".format(dlocs[3])))
                nc_recon_i = dset.createVariable("{}_Recon_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1]), "dim_{}".format(dlocs[2]), "dim_{}".format(dlocs[3])))

                np.reshape(err[i], (dims[dlocs[0]], dims[dlocs[1]], dims[dlocs[2]], dims[dlocs[3]]))
                np.reshape(recon[i], (dims[dlocs[0]], dims[dlocs[1]], dims[dlocs[2]], dims[dlocs[3]]))

                #----Calculate 2-norm of err (with 3 grids)----
                # 4 grids + 4 dimensions means data is structured like f[time][space(z)][space(y)][space(x)]
                if ngrids == 4:
                    nc_2err_r = dset.createVariable("{}_Err_2_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))
                    nc_2err_i = dset.createVariable("{}_Err_2_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))

                    err2 = np.zeros((dims[dlocs[0]]), dtype='complex')
                    for j in range(dims[dlocs[0]]):
                        for z in range(dims[dlocs[1]]):
                            for y in range(dims[dlocs[2]]):
                                for x in range(dims[dlocs[3]]):
                                    err2[j] = err2[j] + err[i][j][z][y][x]**2
                        err2[j] = cmath.sqrt(err2[j])

                #----Calculate 2-norm of err (with 3 grids)----
                # 3 grids + 4 dimensions means data is structured like f[time][group][space(y)][space(x)]
                elif ngrids == 3:
                    nc_2err_r = dset.createVariable("{}_Err_2_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))
                    nc_2err_i = dset.createVariable("{}_Err_2_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))

                    err2 = np.zeros((dims[dlocs[0]], dims[dlocs[1]]), dtype='complex')
                    for j in range(dims[dlocs[0]]):
                        for k in range(dims[dlocs[1]]):
                            for y in range(dims[dlocs[2]]):
                                for x in range(dims[dlocs[3]]):
                                    err2[j][k] = err2[j][k] + err[i][j][k][y][x]**2
                            err2[j][k] = cmath.sqrt(err2[j][k])

                #----Calculate 2-norm of err (with 2 grids)----
                # 2 grids + 4 dimensions means data is structured like f[time][group_1][group_2][space(x)]
                if ngrids == 2:
                    nc_2err_r = dset.createVariable("{}_Err_2_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))
                    nc_2err_i = dset.createVariable("{}_Err_2_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                    "dim_{}".format(dlocs[1])))

                    err2 = np.zeros((dims[dlocs[0]], dims[dlocs[1]], dims[dlocs[2]]), dtype='complex')
                    for j in range(dims[dlocs[0]]):
                        for k in range(dims[dlocs[1]]):
                            for l in range(dims[dlocs[2]]):
                                for x in range(dims[dlocs[3]]):
                                    err2[j][k][l] = err2[j][k][l] + err[i][j][k][l][x]**2
                                err2[j][k][l] = cmath.sqrt(err2[j][k][l])

                #---skipping 2-norm if there are not 2 or 3 grids---
                # with 4 dimensions, anything other than 2, 3 or 4 grids either doesn't make sense or doesn't allow -
                # - for a physically interpretable 2-norm
                else:
                    err2 = []

            else:
                print('unsupported dimension count encountered in Process_Routines.Output')
                quit()

            #--------------------------------------------------
            # Writing err, recon data to dataset
            #--------------------------------------------------
            #writing absolute errors at each grid point between reference and recon data
            nc_err_r[:] = err.real[i]
            nc_err_i[:] = err.imag[i]
            #writing recon data over whole phase space
            nc_recon_r[:] = recon.real[i]
            nc_recon_i[:] = recon.imag[i]

            #creating, recording the maximal absolute error across all grid points
            nc_err_max = dset.createVariable("{}_Err_max".format(dname), "f8", ())
            nc_err_max[:] = np.amax(np.absolute(err))

            #if the error 2-norm was not calculated, then err2 = []
            if len(err2) != 0:
                nc_2err_r[:] = err2.real
                nc_2err_i[:] = err2.imag

                #creating, recording the maximal absolute error in 2-norm across all grid points
                nc_2err_max = dset.createVariable("{}_Err_2max".format(dname), "f8", ())
                nc_2err_max[:] = np.amax(np.absolute(err2))

    dset.close()
