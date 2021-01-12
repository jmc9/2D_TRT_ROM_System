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
                nc_err_r = dset.createVariable("{}_Err_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1])))
                nc_err_i = dset.createVariable("{}_Err_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                               "dim_{}".format(dlocs[1])))

                nc_recon_r = dset.createVariable("{}_Recon_Real".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1])))
                nc_recon_i = dset.createVariable("{}_Recon_Imag".format(dname), "f8", ("dim_{}".format(dlocs[0]),
                                                 "dim_{}".format(dlocs[1])))

            #--------------------------------------------------
            # Creating err, recon variables (in 3D)
            #--------------------------------------------------
            elif ndims == 3:
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

            #--------------------------------------------------
            # Creating err, recon variables (in 4D)
            #--------------------------------------------------
            elif ndims == 4:
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

            else:
                print('unsupported dimension count encountered in Process_Routines.Output')
                quit()

            #--------------------------------------------------
            # Writing err, recon data to dataset
            #--------------------------------------------------
            nc_err_r[:] = err.real[i]
            nc_err_i[:] = err.imag[i]
            nc_recon_r[:] = recon.real[i]
            nc_recon_i[:] = recon.imag[i]

    dset.close()
