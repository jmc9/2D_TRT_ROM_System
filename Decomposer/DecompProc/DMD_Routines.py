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
import math
import cmath

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
                p = g * DMD_Data.rank + m
                lr = DMD_Data.dat.eval.real[p]; (lr_e, lr_c) = tb.pow10(lr)
                li = DMD_Data.dat.eval.imag[p]; (li_e, li_c) = tb.pow10(li)
                l2 = math.sqrt(lr**2 + li**2);  (l2_e, l2_c) = tb.pow10(l2)
                if (li < 0): s = '-'
                else: s = '+'
                title_r = 'Re($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)
                title_i = 'Im($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)

                p = g * DMD_Data.clen * DMD_Data.rank + m * DMD_Data.clen
                wr = DMD_Data.dat.evec.real[p : p + DMD_Data.clen]
                wi = DMD_Data.dat.evec.imag[p : p + DMD_Data.clen]
                pltname_r = '{}_m{}_g{}_re'.format(DMD_Data.name,m+1,g+1)
                pltname_i = '{}_m{}_g{}_im'.format(DMD_Data.name,m+1,g+1)

                if (N_grids == 2): #
                    print('Plot_DMD_evecs doesn not support N_grids = 2 yet for groupwise decompositions')
                    quit()

                elif (N_grids == 3): #
                    N_x = DMD_Data.dims[N_dims-1]
                    N_y = DMD_Data.dims[N_dims-2]

                    xp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-2])

                    pltr.plot_heatmap2D(wr, pltname_r, xp, yp, dir, N_y, N_x, title=title_r)
                    if li != 0.: pltr.plot_heatmap2D(wi, pltname_i, xp, yp, dir, N_y, N_x, title=title_i)

    else:
        for m in plt_modes:
            lr = DMD_Data.dat.eval.real[m]; (lr_e, lr_c) = tb.pow10(lr)
            li = DMD_Data.dat.eval.imag[m]; (li_e, li_c) = tb.pow10(li)
            l2 = math.sqrt(lr**2 + li**2);  (l2_e, l2_c) = tb.pow10(l2)
            if (li < 0): s = '-'
            else: s = '+'
            title_r = 'Re($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)
            title_i = 'Im($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)

            p = m * DMD_Data.clen
            wr = DMD_Data.dat.evec.real[p : p + DMD_Data.clen]
            wi = DMD_Data.dat.evec.imag[p : p + DMD_Data.clen]
            pltname_r = '{}_m{}_re'.format(DMD_Data.name,m+1)
            pltname_i = '{}_m{}_im'.format(DMD_Data.name,m+1)

            if (N_grids == N_dims): #no 'group' dimensions
                if (N_grids == 2): #simple 2D data
                    xp = DMD_Data.grids[N_grids-1].crds

                    pltr.lineplot(pltname_r, xp, wr, dir, title=title_r)
                    if li != 0.: pltr.lineplot(pltname_i, xp, wi, dir, title=title_i)

                elif (N_grids == 3): #
                    N_x = DMD_Data.dims[N_dims-1]
                    N_y = DMD_Data.dims[N_dims-2]

                    xp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-2])

                    pltr.plot_heatmap2D(wr, pltname_r, xp, yp, dir, N_y, N_x, title=title_r)
                    if li != 0.: pltr.plot_heatmap2D(wi, pltname_i, xp, yp, dir, N_y, N_x, title=title_i)


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

                    pltr.plot_heatmap3D(wr, pltname_r, xp, yp, dir, N_z, N_y, N_x, zlab='g', title=title_r)
                    if li != 0.: pltr.plot_heatmap3D(wi, pltname_i, xp, yp, dir, N_z, N_y, N_x, zlab='g', title=title_i)

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

#==================================================================================================================================#
#
#==================================================================================================================================#
def DMD_Sort(DMD_Data):

    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        N_g = len(DMD_Data.dat.N_modes)
        for g in range(N_g):

            evl_p = g * DMD_Data.rank
            lr = DMD_Data.dat.eval.real[evl_p : evl_p + DMD_Data.dat.N_modes[g]]
            li = DMD_Data.dat.eval.imag[evl_p : evl_p + DMD_Data.dat.N_modes[g]]

            evc_p = g * DMD_Data.clen * DMD_Data.rank
            wr = DMD_Data.dat.evec.real[evc_p : evc_p + DMD_Data.dat.N_modes[g] * DMD_Data.clen]
            wi = DMD_Data.dat.evec.imag[evc_p : evc_p + DMD_Data.dat.N_modes[g] * DMD_Data.clen]

            for j in range(DMD_Data.dat.N_modes[g]):
                max = math.sqrt( lr[j]**2 + li[j]**2 )
                loc = j
                for k in range(j+1,DMD_Data.dat.N_modes[g]):
                    val = math.sqrt( lr[k]**2 + li[k]**2 )
                    if (val > max):
                        max = val
                        loc = k

                if (loc != j):
                    evlr = lr[j]
                    evli = li[j]

                    p = j * DMD_Data.clen
                    evcr = wr[p : p + DMD_Data.clen]
                    evci = wi[p : p + DMD_Data.clen]

                    lr[j] = lr[loc]
                    li[j] = li[loc]

                    p2 = loc * DMD_Data.clen
                    wr[p : p + DMD_Data.clen] = wr[p2 : p2 + DMD_Data.clen]
                    wi[p : p + DMD_Data.clen] = wi[p2 : p2 + DMD_Data.clen]

                    lr[loc] = evlr
                    li[loc] = evli

                    wr[p2 : p2 + DMD_Data.clen] = evcr
                    wi[p2 : p2 + DMD_Data.clen] = evci

            DMD_Data.dat.eval.real[evl_p : evl_p + DMD_Data.dat.N_modes[g]] = lr
            DMD_Data.dat.eval.imag[evl_p : evl_p + DMD_Data.dat.N_modes[g]] = li
            DMD_Data.dat.evec.real[evc_p : evc_p + DMD_Data.dat.N_modes[g] * DMD_Data.clen] = wr
            DMD_Data.dat.evec.imag[evc_p : evc_p + DMD_Data.dat.N_modes[g] * DMD_Data.clen] = wi

    else:
        for j in range(DMD_Data.dat.N_modes):
            max = math.sqrt( DMD_Data.dat.eval.real[j]**2 + DMD_Data.dat.eval.imag[j]**2 )
            loc = j
            for k in range(j+1,DMD_Data.dat.N_modes):
                val = math.sqrt( DMD_Data.dat.eval.real[k]**2 + DMD_Data.dat.eval.imag[k]**2 )
                if (val > max):
                    max = val
                    loc = k

            if (loc != j):
                evlr = DMD_Data.dat.eval.real[j]
                evli = DMD_Data.dat.eval.imag[j]

                p = j * DMD_Data.clen
                evcr = DMD_Data.dat.evec.real[p : p + DMD_Data.clen]
                evci = DMD_Data.dat.evec.imag[p : p + DMD_Data.clen]

                DMD_Data.dat.eval.real[j] = DMD_Data.dat.eval.real[loc]
                DMD_Data.dat.eval.imag[j] = DMD_Data.dat.eval.imag[loc]

                p2 = loc * DMD_Data.clen
                DMD_Data.dat.evec.real[p : p + DMD_Data.clen] = DMD_Data.dat.evec.real[p2 : p2 + DMD_Data.clen]
                DMD_Data.dat.evec.imag[p : p + DMD_Data.clen] = DMD_Data.dat.evec.imag[p2 : p2 + DMD_Data.clen]

                DMD_Data.dat.eval.real[loc] = evlr
                DMD_Data.dat.eval.imag[loc] = evli

                DMD_Data.dat.evec.real[p2 : p2 + DMD_Data.clen] = evcr
                DMD_Data.dat.evec.imag[p2 : p2 + DMD_Data.clen] = evci

    return DMD_Data

#==================================================================================================================================#
#
#==================================================================================================================================#
def Coef_Calc(DMD_Data,Prob_Data):

    b = []
    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        quit()

    else:
        init = Prob_Data.dat[0][0 : DMD_Data.clen]
        p = 0
        for j in range(DMD_Data.dat.N_modes):
            prod = complex(0., 0.)
            nrm  = complex(0., 0.)
            for i in range(DMD_Data.clen):
                c1 = complex(init[i], 0.)
                c2 = complex(DMD_Data.dat.evec.real[p], DMD_Data.dat.evec.imag[p])

                prod = prod + c1*c2
                nrm = nrm + c2*c2

            b.append( prod / nrm )
            p = p + DMD_Data.clen

    return b
