#!/usr/bin/env python3
#==================================================================================================================================#
#
# DMD_Routines.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import numpy as np

#local tools
import Plotters as pltr
import ToolBox as tb

from Classes import DMD

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD(DMD_Data,dir,plt_modes=[0,1,2],evecs=True,evals=True):

    if evecs:
        evec_dir = dir+'/Eigenvectors'
        tb.dirset(evec_dir)
        Plot_DMD_evecs(DMD_Data,plt_modes,evec_dir)

    if evals:
        eval_dir = dir+'/Eigenvalues'
        tb.dirset(eval_dir)
        Plot_DMD_evals(DMD_Data,eval_dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD_evecs(DMD_Data,plt_modes,dir):
    N_grids = len(DMD_Data.grids)
    N_dims  = len(DMD_Data.dims)

    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        N_g = DMD_Data.dims[1]
        for g in range(N_g):
            for m in plt_modes:
                p = g * DMD_Data.clen * DMD_Data.rank + m * DMD_Data.clen
                w = DMD_Data.dat.evec.real[p : p + DMD_Data.clen]
                title = '{}_m{}_g{}'.format(DMD_Data.name,g+1,m+1)

                if (N_grids == 2): #
                    print('Plot_DMD_evecs doesn not support N_grids = 2 yet for groupwise decompositions')
                    quit()

                elif (N_grids == 3): #
                    N_x = DMD_Data.dims[N_dims-1]
                    N_y = DMD_Data.dims[N_dims-2]

                    xp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-2])

                    pltr.plot_heatmap2D(w, title, xp, yp, dir, N_y, N_x)

    else:
        for m in plt_modes:
            p = m * DMD_Data.clen
            w = DMD_Data.dat.evec.real[p : p + DMD_Data.clen]
            title = '{}_m{}'.format(DMD_Data.name,m+1)

            if (N_grids == N_dims): #no 'group' dimensions
                if (N_grids == 2): #simple 2D data
                    print('Plot_DMD_evecs doesn not support N_grids = 2 yet for groupwise decompositions')
                    quit()

                elif (N_grids == 3): #
                    N_x = DMD_Data.dims[N_dims-1]
                    N_y = DMD_Data.dims[N_dims-2]

                    xp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-2])

                    pltr.plot_heatmap2D(w, title, xp, yp, dir, N_y, N_x)


            elif (N_grids == N_dims - 1): #only 1 'group' dimension exists
                if (N_grids == 2): #
                    print('Plot_DMD_evecs doesn not support N_grids = 2 yet for groupwise decompositions')
                    quit()

                elif (N_grids == 3): #
                    N_x = DMD_Data.dims[N_dims-1]
                    N_y = DMD_Data.dims[N_dims-2]
                    N_z = DMD_Data.dims[N_dims-3]

                    xp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-2])

                    pltr.plot_heatmap3D(w, title, xp, yp, dir, N_z, N_y, N_x, zlab='g')

            else:
                print('Plot_DMD_evecs can only handle data with a single group dimension right now!')
                quit()

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD_evals(DMD_Data,dir):
    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        N_g = DMD_Data.dims[1]
        for g in range(N_g):
            p = g * DMD_Data.rank
            lr = DMD_Data.dat.eval.real[p : p + DMD_Data.dat.N_modes[g]]
            li = DMD_Data.dat.eval.imag[p : p + DMD_Data.dat.N_modes[g]]

            title = '{}_Evals_g{}'.format(DMD_Data.name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

            title = '{}_Evals_g{}_ucirc'.format(DMD_Data.name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlim=[-1.5,1.5],ylim=[-1.5,1.5],cr=1,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

    else:
        lr = DMD_Data.dat.eval.real[0:DMD_Data.dat.N_modes]
        li = DMD_Data.dat.eval.imag[0:DMD_Data.dat.N_modes]

        title = '{}_Evals'.format(DMD_Data.name)
        pltr.scatterplot([lr,li],title,dir,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

        title = '{}_Evals_ucirc'.format(DMD_Data.name)
        pltr.scatterplot([lr,li],title,dir,xlim=[-1.05,1.05],ylim=[-1.05,1.05],cr=1,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_DMD(dset,name):

    N_modes = getattr(dset,'N_modes_'+name)

    Evals = Read_Evals(dset, name, N_modes)
    Evecs = Read_Evecs(dset, name, N_modes)

    DMD_dat = DMD(N_modes, Evals, Evecs)

    return DMD_dat

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Dims(dset,name):
    dims = dset['W_real_'+name].dimensions
    ndims = len(dims)

    rank = dset.dimensions[dims[0]].size
    clen = dset.dimensions[dims[ndims-1]].size

    return (rank, clen)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Evals(dset,name,Nmodes):
    #reading eigenvalues from dataset
    L_real = dset['L_real_'+name][:] #real component
    L_imag = dset['L_imag_'+name][:] #imaginary component

    #checking whether evals are multi-d
    if hasattr(Nmodes,'__len__'):
        #if evals are multi-d, flattening them into 1D arrays
        L = [tb.Flatten2D(L_real), tb.Flatten2D(L_imag)]

    else:
        L = [L_real, L_imag]

    return L

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Evecs(dset,name,Nmodes):
    #reading eigenvectors from dataset
    W_real = dset['W_real_'+name][:] #real component
    W_imag = dset['W_imag_'+name][:] #imaginary component

    #flattening eigenvector arrays into 1D arrays
    if hasattr(Nmodes,'__len__'):
        W = [tb.Flatten3D(W_real), tb.Flatten3D(W_imag)]

    else:
        W = [tb.Flatten2D(W_real), tb.Flatten2D(W_imag)]

    return W
