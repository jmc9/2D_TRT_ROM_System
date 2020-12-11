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
from Classes import DMD_opts

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD(dset,name,modes,clen,rank,N_g,N_y,N_x,xp,yp,dir,evecs=True,evals=True):
    DMD_dat = Read_DMD(dset,name)

    if evecs:
        evec_dir = dir+'/Eigenvectors'
        tb.dirset(evec_dir)
        Plot_DMD_evecs(DMD_dat,name,evec_dir,clen,rank,modes,N_g,N_y,N_x,xp,yp)

    if evals:
        eval_dir = dir+'/Eigenvalues'
        tb.dirset(eval_dir)
        Plot_DMD_evals(DMD_dat,name,eval_dir,rank,N_g)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD_evecs(DMD_dat,name,dir,clen,rank,modes,N_g,N_y,N_x,xp,yp):
    if hasattr(DMD_dat.N_modes,'__len__'):
        for g in range(N_g):
            for m in modes:
                p = g*clen*rank + m*clen
                w = DMD_dat.evec.real[p:p+clen]
                title = '{}_m{}_g{}'.format(name,g+1,m+1)
                pltr.plot_heatmap2D(w,title,xp,yp,dir,N_y,N_x)
    else:
        for m in modes:
            p = m*clen
            w = DMD_dat.evec.real[p:p+clen]
            title = '{}_m{}'.format(name,m+1)
            pltr.plot_heatmap3D(w,title,xp,yp,dir,N_g,N_y,N_x,zlab='g')

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD_evals(DMD_dat,name,dir,rank,N_g):
    if hasattr(DMD_dat.N_modes,'__len__'):
        for g in range(N_g):
            p = g*rank
            lr = DMD_dat.eval.real[p:p+DMD_dat.N_modes[g]]
            li = DMD_dat.eval.imag[p:p+DMD_dat.N_modes[g]]

            title = '{}_Evals_g{}'.format(name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

            title = '{}_Evals_g{}_ucirc'.format(name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlim=[-1.5,1.5],ylim=[-1.5,1.5],cr=1,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

    else:
        lr = DMD_dat.eval.real[0:DMD_dat.N_modes]
        li = DMD_dat.eval.imag[0:DMD_dat.N_modes]

        title = '{}_Evals'.format(name)
        pltr.scatterplot([lr,li],title,dir,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

        title = '{}_Evals_ucirc'.format(name)
        pltr.scatterplot([lr,li],title,dir,xlim=[-1.5,1.5],ylim=[-1.5,1.5],cr=1,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_DMD(dset,name):

    N_modes = getattr(dset,'N_modes_'+name)

    Evals = Read_Evals(dset,name,N_modes)
    Evecs = Read_Evecs(dset,name,N_modes)

    DMD_dat = DMD(name,N_modes,Evals,Evecs)

    return DMD_dat


#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Nmodes(dset,names):
    N_modes = [];

    for name in names:
        N_modes.append(getattr(dset,'N_modes_'+name))

    return N_modes

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
