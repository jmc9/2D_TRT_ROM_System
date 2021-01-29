#==================================================================================================================================#
#
# TRT_2D_Handling.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import numpy as np
import math

#netcdf tools
from netCDF4 import Dataset as ncds

#==================================================================================================================================#
#
#==================================================================================================================================#
def Data_Select(file_name, data_name, tstart, tend, group=0):
    eps=1e-12
    dset = ncds(file_name,'r') #open ncdf datafile/dataset

    ncdat = dset[data_name]

    # dat = np.array(ncdat[:])
    dat = ncdat[:]
    tp_ = dset[ getattr( ncdat, 'grid0' ) ][:]
    yp = dset[ getattr( ncdat, 'grid1' ) ][:]
    xp = dset[ getattr( ncdat, 'grid2' ) ][:]

    tp = [t for t in tp_ if ((t>=tstart-eps)and(t<=tend+eps))]

    ps = -1
    pe = -1
    for i in range(len(tp_)):
        if ((ps == -1) and (tp_[i] > tstart-eps)):
            ps = i
        if ((pe == -1) and (tp_[i] > tend-eps) and (tp_[i] < tend+eps)):
            pe = i
            break

    if group == 0:
        f = dat[ps:pe+1, :, :, :]
    else:
        f = dat[ps:pe+1, group-1, :, :]

    dset.close()

    return f, tp, yp, xp
