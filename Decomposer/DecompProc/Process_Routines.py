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
import copy
import os

#netcdf tools
from netCDF4 import Dataset as ncds

#local tools
import DMD_Routines as dmdr
import ToolBox as tb
import Plotters as pltr

#==================================================================================================================================#
#
#==================================================================================================================================#
def Process_Data(dset, Dcmp_Data, Prob_Data, plotdir, trunc_eps, trunc_opt):
    recon = []
    n_trunc = []

    N_data = len(Dcmp_Data)
    for i in range(N_data):
        dcmp = Dcmp_Data[i]
        recon.append([])

        dir = '{}/{}'.format(plotdir,dcmp.name)
        tb.dirset(dir)

        if (dcmp.type in ['DMD','DMDg']):
            (dcmp.rank, dcmp.clen) = dmdr.Read_Dims(dset,dcmp.name)
            dcmp.dat = dmdr.Read_DMD(dset,dcmp.name,dcmp.clen,dcmp.rank)

            if dcmp.opt[0] == 0:
                n_trunc.append( dmdr.DMD_Sort(dcmp, trunc_eps) )
                if dcmp.type == 'DMD': plt_modes = np.arange(dcmp.dat.N_modes)
                else: plt_modes = [np.arange(nmodes) for nmodes in dcmp.dat.N_modes]
                # dmdr.Plot_DMD(dcmp,dir,plt_modes)

                if len(Prob_Data) != 0:
                    if Prob_Data[i].opt[0] == 0:
                        coef = dmdr.Coef_Calc(dcmp, Prob_Data[i].dat[0], n_trunc[i], trunc_opt)

                        len1 = len(Prob_Data[i].dat)
                        len2 = len(Prob_Data[i].dat[0])
                        edat = np.zeros([len1, len2],dtype='complex')
                        for j in range(len1):
                            edat[j] = dmdr.Expand(dcmp, coef, j, n_trunc[i], trunc_opt)
                        recon[i] = copy.deepcopy(edat)

            else:
                dmdr.Plot_DMD(dcmp,dir,plt_modes,evecs=False)

    return recon, n_trunc

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
        err.append([]) #appending blank array first -> stays consistent with N_data and allows easy flag for uncalculated errors

        #straightforward error calculations for opt[0] = 0 (standard data array)
        if Prob_Data[j].opt[0] == 0:
            tlen = len(Prob_Data[j].dat)
            clen = len(Prob_Data[j].dat[0])

            err_ = np.zeros((tlen,clen)) #temporary array to hold j-th data errors
            for t in range(tlen):
                for i in range(clen):
                    err_[t][i] = abs(Prob_Data[j].dat[t][i] - recon[j][t][i])
            err[j] = err_ #now putting j-th calculated errors into err array for output

    return err

#==================================================================================================================================#
#
#==================================================================================================================================#
def errplot(proc_dir, xp, err2_abs, err2_rel, err_inf_abs, err_inf_rel, xlabel=''):
    drop = proc_dir+'/err_plots'
    tb.dirset(drop)

    plotdat = np.transpose(err2_abs)
    drop_ = drop+'/abs_2'; tb.dirset(drop_)
    errplot_single(drop_, plotdat, xp, 'Err_2_abs', xlabel)

    plotdat = np.transpose(err2_rel)
    drop_ = drop+'/rel_2'; tb.dirset(drop_)
    errplot_single(drop_, plotdat, xp, 'Err_2_rel', xlabel)

    plotdat = np.transpose(err_inf_abs)
    drop_ = drop+'/abs_inf'; tb.dirset(drop_)
    errplot_single(drop_, plotdat, xp, 'Err_inf_abs', xlabel)

    plotdat = np.transpose(err_inf_rel)
    drop_ = drop+'/rel_inf'; tb.dirset(drop_)
    errplot_single(drop_, plotdat, xp, 'Err_inf_rel', xlabel)

#==================================================================================================================================#
#
#==================================================================================================================================#
def errplot_single(drop, err, xp, name, xlabel=''):
    dims = np.shape(err)
    ndims = len(dims)

    if ndims == 1:
        pltr.lineplot(name, xp, err, drop, yscale='log', xlabel=xlabel)

    elif ndims == 2:
        legend = ['g{}'.format(g+1) for g in range(dims[0])]
        pltr.lineplot('{}_allg'.format(name), xp, err, drop, yscale='log', xlabel=xlabel, legend=legend)
        for g in range(dims[0]):
            pltr.lineplot('{}_g{}'.format(name,g+1), xp, err[g], drop, yscale='log', xlabel=xlabel)

    else:
        print('Process_Routines/errplot_single cannot handle ndims not in [1,2] right now')
        quit()


#==================================================================================================================================#
#
#==================================================================================================================================#
def Output(recon, err, Prob_Data, n_trunc, proc_dir):
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
    all_dims = []
    d = 0
    for i in range(N_data):
        dlocs = []
        if Prob_Data[i].opt[0] == 0:
            ndims = len(Prob_Data[i].dims)
            ngrids = len(Prob_Data[i].grids)
            dname = Prob_Data[i].name
            dims = Prob_Data[i].dims

            for j in range(ndims):
                if not (dims[j].name in all_dims):
                    nc_dim = dset.createDimension(dims[j].name, dims[j].len)
                    all_dims.append(dims[j].name)

            #--------------------------------------------------
            # Creating err, recon variables (in 2D)
            #--------------------------------------------------
            if ndims == 2:
                #-------Create err, recon nc variables---------
                nc_err = dset.createVariable("{}_Err".format(dname), "f8", (dims[0].name, dims[1].name) )

                nc_recon_r = dset.createVariable("{}_Recon_Real".format(dname), "f8", (dims[0].name, dims[1].name) )
                nc_recon_i = dset.createVariable("{}_Recon_Imag".format(dname), "f8", (dims[0].name, dims[1].name) )

                nc_2err_abs = dset.createVariable("{}_Err_2_abs".format(dname), "f8", (dims[0].name) )
                nc_2err_rel = dset.createVariable("{}_Err_2_rel".format(dname), "f8", (dims[0].name) )
                nc_ierr_abs = dset.createVariable("{}_Err_inf_abs".format(dname), "f8", (dims[0].name) )
                nc_ierr_rel = dset.createVariable("{}_Err_inf_rel".format(dname), "f8", (dims[0].name) )

                nc_trunc = dset.createVariable("{}_n_trunc".format(dname), "f8", ())

                #--Move err, recon data into temporary arrays--
                err_ = err[i]
                recon_ = recon[i]
                pdat_ = Prob_Data[i].dat

                #----Calculate 2-norm of err (with 2 grids)----
                # 2 grids + 2 dimensions means data is structured like f[time][space(x)]
                if ngrids == 2:
                    err2_abs = np.zeros(dims[0].len)
                    err2_rel = np.zeros(dims[0].len)
                    err_inf_abs = np.zeros(dims[0].len)
                    err_inf_rel = np.zeros(dims[0].len)
                    for j in range(dims[0].len):
                        err2_abs[j] = tb.norm(err_[j],2)
                        err_inf_abs[j] = tb.norm(err_[j],0)

                        err2_rel[j] = err2_abs[j] / tb.norm(pdat_[j],2)
                        err_inf_rel[j] = err_inf_abs[j] / tb.norm(pdat_[j],0)

                    errplot(proc_dir, Prob_Data[i].grids[0].crds, err2_abs, err2_rel, err_inf_abs, err_inf_rel)

                #---skipping 2-norm if there are not 2 grids---
                # with 2 dimensions, anything other than 2 grids either doesn't make sense or doesn't allow -
                # - for a physically interpretable 2-norm
                else:
                    err2_abs = []
                    err2_rel = []
                    err_inf_abs = []
                    err_inf_rel = []

            #--------------------------------------------------
            # Creating err, recon variables (in 3D)
            #--------------------------------------------------
            elif ndims == 3:
                #-------Create err, recon nc variables---------
                nc_err = dset.createVariable("{}_Err".format(dname), "f8", (dims[0].name, dims[1].name, dims[2].name) )

                nc_recon_r = dset.createVariable("{}_Recon_Real".format(dname), "f8", (dims[0].name, dims[1].name, dims[2].name) )
                nc_recon_i = dset.createVariable("{}_Recon_Imag".format(dname), "f8", (dims[0].name, dims[1].name, dims[2].name) )

                #--Move err, recon data into temporary arrays--
                # these temporary arrays also have the correct shape for output
                err_ = np.reshape(err[i], (dims[0].len, dims[1].len, dims[2].len))
                recon_ = np.reshape(recon[i], (dims[0].len, dims[1].len, dims[2].len))
                pdat_ = np.reshape(Prob_Data[i].dat, (dims[0].len, dims[1].len, dims[2].len))

                #----Calculate 2-norm of err (with 3 grids)----
                # 3 grids + 3 dimensions means data is structured like f[time][space(y)][space(x)]
                if ngrids == 3:
                    nc_2err_abs = dset.createVariable("{}_Err_2_abs".format(dname), "f8", (dims[0].name) )
                    nc_2err_rel = dset.createVariable("{}_Err_2_rel".format(dname), "f8", (dims[0].name) )
                    nc_ierr_abs = dset.createVariable("{}_Err_inf_abs".format(dname), "f8", (dims[0].name) )
                    nc_ierr_rel = dset.createVariable("{}_Err_inf_rel".format(dname), "f8", (dims[0].name) )

                    nc_trunc = dset.createVariable("{}_n_trunc".format(dname), "f8", ())

                    err2_abs = np.zeros(dims[0].len)
                    err2_rel = np.zeros(dims[0].len)
                    err_inf_abs = np.zeros(dims[0].len)
                    err_inf_rel = np.zeros(dims[0].len)
                    for j in range(dims[0].len):
                        err2_abs[j] = tb.norm(err_[j],2)
                        err_inf_abs[j] = tb.norm(err_[j],0)

                        err2_rel[j] = err2_abs[j] / tb.norm(pdat_[j],2)
                        err_inf_rel[j] = err_inf_abs[j] / tb.norm(pdat_[j],0)

                    errplot(proc_dir, Prob_Data[i].grids[0].crds, err2_abs, err2_rel, err_inf_abs, err_inf_rel)

                #----Calculate 2-norm of err (with 2 grids)----
                # 2 grids + 3 dimensions means data is structured like f[time][group][space(x)]
                elif ngrids == 2:
                    nc_2err_abs = dset.createVariable("{}_Err_2_abs".format(dname), "f8", (dims[0].name, dims[1].name) )
                    nc_2err_rel = dset.createVariable("{}_Err_2_rel".format(dname), "f8", (dims[0].name, dims[1].name) )
                    nc_ierr_abs = dset.createVariable("{}_Err_inf_abs".format(dname), "f8", (dims[0].name, dims[1].name) )
                    nc_ierr_rel = dset.createVariable("{}_Err_inf_rel".format(dname), "f8", (dims[0].name, dims[1].name) )

                    if isinstance(n_trunc[i],(list, np.ndarray)):
                        nc_trunc = dset.createVariable("{}_n_trunc".format(dname), "f8", (dims[1].name))
                    else:
                        nc_trunc = dset.createVariable("{}_n_trunc".format(dname), "f8", ())

                    err2_abs = np.zeros((dims[0].len, dims[1].len))
                    err2_rel = np.zeros((dims[0].len, dims[1].len))
                    err_inf_abs = np.zeros((dims[0].len, dims[1].len))
                    err_inf_rel = np.zeros((dims[0].len, dims[1].len))
                    for j in range(dims[0].len):
                        for k in range(dims[1].len):
                            err2_abs[j][k] = tb.norm(err_[j][k],2)
                            err_inf_abs[j][k] = tb.norm(err_[j][k],0)

                            err2_rel[j][k] = err2_abs[j][k] / tb.norm(pdat_[j][k],2)
                            err_inf_rel[j][k] = err_inf_abs[j][k] / tb.norm(pdat_[j][k],0)

                    errplot(proc_dir, Prob_Data[i].grids[0].crds, err2_abs, err2_rel, err_inf_abs, err_inf_rel)

                #---skipping 2-norm if there are not 2 or 3 grids---
                # with 3 dimensions, anything other than 2 or 3 grids either doesn't make sense or doesn't allow -
                # - for a physically interpretable 2-norm
                else:
                    err2_abs = []
                    err2_rel = []
                    err_inf_abs = []
                    err_inf_rel = []

            #--------------------------------------------------
            # Creating err, recon variables (in 4D)
            #--------------------------------------------------
            elif ndims == 4:
                #-------Create err, recon nc variables---------
                nc_err = dset.createVariable("{}_Err".format(dname), "f8", (dims[0].name, dims[1].name, dims[2].name, dims[3].name) )

                nc_recon_r = dset.createVariable("{}_Recon_Real".format(dname), "f8", (dims[0].name, dims[1].name, dims[2].name, dims[3].name) )
                nc_recon_i = dset.createVariable("{}_Recon_Imag".format(dname), "f8", (dims[0].name, dims[1].name, dims[2].name, dims[3].name) )

                #--Move err, recon data into temporary arrays--
                # these temporary arrays also have the correct shape for output
                err_ = np.reshape(err[i], (dims[0].len, dims[1].len, dims[2].len, dims[3].len))
                recon_ = np.reshape(recon[i], (dims[0].len, dims[1].len, dims[2].len, dims[3].len))
                pdat_ = np.reshape(Prob_Data[i].dat, (dims[0].len, dims[1].len, dims[2].len, dims[3].len))

                #----Calculate 2-norm of err (with 3 grids)----
                # 4 grids + 4 dimensions means data is structured like f[time][space(z)][space(y)][space(x)]
                if ngrids == 4:
                    nc_2err_abs = dset.createVariable("{}_Err_2_abs".format(dname), "f8", (dims[0].name) )
                    nc_2err_rel = dset.createVariable("{}_Err_2_rel".format(dname), "f8", (dims[0].name) )
                    nc_ierr_abs = dset.createVariable("{}_Err_inf_abs".format(dname), "f8", (dims[0].name) )
                    nc_ierr_rel = dset.createVariable("{}_Err_inf_rel".format(dname), "f8", (dims[0].name) )

                    nc_trunc = dset.createVariable("{}_n_trunc".format(dname), "f8", ())

                    err2_abs = np.zeros(dims[0].len)
                    err2_rel = np.zeros(dims[0].len)
                    err_inf_abs = np.zeros(dims[0].len)
                    err_inf_rel = np.zeros(dims[0].len)
                    for j in range(dims[0].len):
                        err2_abs[j] = tb.norm(err_[j],2)
                        err_inf_abs[j] = tb.norm(err_[j],0)

                        err2_rel[j] = err2_abs[j] / tb.norm(pdat_[j],2)
                        err_inf_rel[j] = err_inf_abs[j] / tb.norm(pdat_[j],0)

                    errplot(proc_dir, Prob_Data[i].grids[0].crds, err2_abs, err2_rel, err_inf_abs, err_inf_rel)

                #----Calculate 2-norm of err (with 3 grids)----
                # 3 grids + 4 dimensions means data is structured like f[time][group][space(y)][space(x)]
                elif ngrids == 3:
                    nc_2err_abs = dset.createVariable("{}_Err_2_abs".format(dname), "f8", (dims[0].name, dims[1].name) )
                    nc_2err_rel = dset.createVariable("{}_Err_2_rel".format(dname), "f8", (dims[0].name, dims[1].name) )
                    nc_ierr_abs = dset.createVariable("{}_Err_inf_abs".format(dname), "f8", (dims[0].name, dims[1].name) )
                    nc_ierr_rel = dset.createVariable("{}_Err_inf_rel".format(dname), "f8", (dims[0].name, dims[1].name) )

                    if isinstance(n_trunc[i],(list, np.ndarray)):
                        nc_trunc = dset.createVariable("{}_n_trunc".format(dname), "f8", (dims[1].name))
                    else:
                        nc_trunc = dset.createVariable("{}_n_trunc".format(dname), "f8", ())

                    err2_abs = np.zeros((dims[0].len, dims[1].len))
                    err2_rel = np.zeros((dims[0].len, dims[1].len))
                    err_inf_abs = np.zeros((dims[0].len, dims[1].len))
                    err_inf_rel = np.zeros((dims[0].len, dims[1].len))
                    for j in range(dims[0].len):
                        for k in range(dims[1].len):
                            err2_abs[j][k] = tb.norm(err_[j][k],2)
                            err_inf_abs[j][k] = tb.norm(err_[j][k],0)

                            err2_rel[j][k] = err2_abs[j][k] / tb.norm(pdat_[j][k],2)
                            err_inf_rel[j][k] = err_inf_abs[j][k] / tb.norm(pdat_[j][k],0)

                    errplot(proc_dir, Prob_Data[i].grids[0].crds, err2_abs, err2_rel, err_inf_abs, err_inf_rel)
                    # errplot(proc_dir, Prob_Data[i].grids[0].crds*10., err2_abs, err2_rel, err_inf_abs, err_inf_rel, xlabel='ns')

                #----Calculate 2-norm of err (with 2 grids)----
                # 2 grids + 4 dimensions means data is structured like f[time][group_1][group_2][space(x)]
                elif ngrids == 2:
                    nc_2err_abs = dset.createVariable("{}_Err_2_abs".format(dname), "f8", (dims[0].name, dims[1].name, dims[3].name) )
                    nc_2err_rel = dset.createVariable("{}_Err_2_rel".format(dname), "f8", (dims[0].name, dims[1].name, dims[3].name) )
                    nc_ierr_abs = dset.createVariable("{}_Err_inf_abs".format(dname), "f8", (dims[0].name, dims[1].name, dims[3].name) )
                    nc_ierr_rel = dset.createVariable("{}_Err_inf_rel".format(dname), "f8", (dims[0].name, dims[1].name, dims[3].name) )

                    err2_abs = np.zeros((dims[0].len, dims[1].len, dims[2].len))
                    err2_rel = np.zeros((dims[0].len, dims[1].len, dims[2].len))
                    err_inf_abs = np.zeros((dims[0].len, dims[1].len, dims[2].len))
                    err_inf_rel = np.zeros((dims[0].len, dims[1].len, dims[2].len))
                    for j in range(dims[0].len):
                        for k in range(dims[1].len):
                            for l in range(dims[2].len):
                                err2_abs[j][k][l] = tb.norm(err_[j][k][l],2)
                                err_inf_abs[j][k][l] = tb.norm(err_[j][k][l],0)

                                err2_rel[j][k][l] = err2_abs[j][k][l] / tb.norm(pdat_[j][k][l],2)
                                err_inf_rel[j][k][l] = err_inf_abs[j][k][l] / tb.norm(pdat_[j][k][l],0)

                #---skipping 2-norm if there are not 2 or 3 grids---
                # with 4 dimensions, anything other than 2, 3 or 4 grids either doesn't make sense or doesn't allow -
                # - for a physically interpretable 2-norm
                else:
                    err2_abs = []
                    err2_rel = []
                    err_inf_abs = []
                    err_inf_rel = []

            else:
                print('unsupported dimension count encountered in Process_Routines.Output')
                quit()

            #--------------------------------------------------
            # Writing err, recon data to dataset
            #--------------------------------------------------
            #writing absolute errors at each grid point between reference and recon data
            nc_err[:] = err_
            #writing recon data over whole phase space
            nc_recon_r[:] = recon_.real
            nc_recon_i[:] = recon_.imag

            nc_trunc[:] = n_trunc[i]

            #if the error 2-norm was not calculated, then err2 = []
            if (len(err2_abs) != 0):
                nc_2err_abs[:] = err2_abs
                nc_2err_rel[:] = err2_rel
                nc_ierr_abs[:] = err_inf_abs
                nc_ierr_rel[:] = err_inf_rel

                #creating, recording the maximal absolute error across all grid points
                nc_ierr_abs_max = dset.createVariable("{}_Err_inf_max_abs".format(dname), "f8", ())
                nc_ierr_abs_max[:] = np.amax(np.absolute(err_inf_abs))

                nc_ierr_rel_max = dset.createVariable("{}_Err_inf_max_rel".format(dname), "f8", ())
                nc_ierr_rel_max[:] = np.amax(np.absolute(err_inf_rel))

                #creating, recording the maximal absolute error in 2-norm across all grid points
                nc_2err_abs_max = dset.createVariable("{}_Err_2_max_abs".format(dname), "f8", ())
                nc_2err_abs_max[:] = np.amax(np.absolute(err2_abs))

                nc_2err_rel_max = dset.createVariable("{}_Err_2_max_rel".format(dname), "f8", ())
                nc_2err_rel_max[:] = np.amax(np.absolute(err2_rel))

    dset.close()
    os.system('cp proc_summary.h5 {}'.format(proc_dir))

    return
