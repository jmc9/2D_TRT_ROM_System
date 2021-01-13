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
                lr = DMD_Data.dat.eval[p].real; (lr_e, lr_c) = tb.pow10(lr)
                li = DMD_Data.dat.eval[p].imag; (li_e, li_c) = tb.pow10(li)
                l2 = math.sqrt(lr**2 + li**2);  (l2_e, l2_c) = tb.pow10(l2)
                if (li < 0): s = '-'
                else: s = '+'
                title_r = 'Re($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)
                title_i = 'Im($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)

                p = g * DMD_Data.clen * DMD_Data.rank + m * DMD_Data.clen
                wr = DMD_Data.dat.evec[p : p + DMD_Data.clen].real
                wi = DMD_Data.dat.evec[p : p + DMD_Data.clen].imag
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
            lr = DMD_Data.dat.eval[m].real; (lr_e, lr_c) = tb.pow10(lr)
            li = DMD_Data.dat.eval[m].imag; (li_e, li_c) = tb.pow10(li)
            l2 = math.sqrt(lr**2 + li**2);  (l2_e, l2_c) = tb.pow10(l2)
            if (li < 0): s = '-'
            else: s = '+'
            title_r = 'Re($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)
            title_i = 'Im($\omega$) \n  $\lambda$ = ${:.2f}\cdot 10^{:d} {} {:.2f}\cdot 10^{:d} i = | {:.2f}\cdot 10^{:d} |$'.format(lr_c,lr_e,s,abs(li_c),li_e,l2_c,l2_e)

            p = m * DMD_Data.clen
            wr = DMD_Data.dat.evec[p : p + DMD_Data.clen].real
            wi = DMD_Data.dat.evec[p : p + DMD_Data.clen].imag
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
            lr = DMD_Data.dat.eval[p : p + DMD_Data.dat.N_modes[g]].real
            li = DMD_Data.dat.eval[p : p + DMD_Data.dat.N_modes[g]].imag

            title = '{}_Evals_g{}'.format(DMD_Data.name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

            title = '{}_Evals_g{}_ucirc'.format(DMD_Data.name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlim=[-1.5,1.5],ylim=[-1.5,1.5],cr=1,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

            lr = DMD_Data.dat.expval[p : p + DMD_Data.dat.N_modes[g]].real
            li = DMD_Data.dat.expval[p : p + DMD_Data.dat.N_modes[g]].imag

            title = '{}_Expvals_g{}'.format(DMD_Data.name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlabel='Re($e\lambda$)',ylabel='Imag($e\lambda$)')

            title = '{}_Expvals_g{}_ucirc'.format(DMD_Data.name,g+1)
            pltr.scatterplot([lr,li],title,dir,xlim=[-1.5,1.5],ylim=[-1.5,1.5],cr=1,xlabel='Re($e\lambda$)',ylabel='Imag($e\lambda$)')

    else:
        lr = DMD_Data.dat.eval[0:DMD_Data.dat.N_modes].real
        li = DMD_Data.dat.eval[0:DMD_Data.dat.N_modes].imag

        title = '{}_Evals'.format(DMD_Data.name)
        pltr.scatterplot([lr,li],title,dir,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

        title = '{}_Evals_ucirc'.format(DMD_Data.name)
        pltr.scatterplot([lr,li],title,dir,xlim=[-1.05,1.05],ylim=[-1.05,1.05],cr=1,xlabel='Re($\lambda$)',ylabel='Imag($\lambda$)')

        lr = DMD_Data.dat.expval[0:DMD_Data.dat.N_modes].real
        li = DMD_Data.dat.expval[0:DMD_Data.dat.N_modes].imag

        title = '{}_Expvals'.format(DMD_Data.name)
        pltr.scatterplot([lr,li],title,dir,xlabel='Re($e\lambda$)',ylabel='Imag($e\lambda$)')

        title = '{}_Expvals_ucirc'.format(DMD_Data.name)
        pltr.scatterplot([lr,li],title,dir,xlim=[-1.05,1.05],ylim=[-1.05,1.05],cr=1,xlabel='Re($e\lambda$)',ylabel='Imag($e\lambda$)')

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_DMD(dset,name):

    N_modes = getattr(dset,'N_modes_'+name)

    Evals = Read_Evals(dset, name, N_modes)
    Evecs = Read_Evecs(dset, name, N_modes)
    Expvals = Read_Expvals(dset, name, N_modes)

    DMD_dat = DMD(N_modes, Evals, Evecs, Expvals)

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
def Read_Evals(dset,name,Nmodes,evn='L'):
    #reading eigenvalues from dataset
    L_real = dset['{}_real_'.format(evn)+name][:] #real component
    L_imag = dset['{}_imag_'.format(evn)+name][:] #imaginary component

    #checking whether evals are multi-d
    if hasattr(Nmodes,'__len__'):
        #if evals are multi-d, flattening them into 1D arrays
        L_real = tb.Flatten2D(L_real)
        L_imag = tb.Flatten2D(L_imag)

    L = np.zeros(Nmodes, dtype='complex')
    for i in range(Nmodes):
        L[i] = complex(L_real[i], L_imag[i])

    return L

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Expvals(dset,name,Nmodes):
    eL = Read_Evals(dset,name,Nmodes,'eL')

    return eL

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Evecs(dset,name,Nmodes):
    #reading eigenvectors from dataset
    W_real = dset['W_real_'+name][:] #real component
    W_imag = dset['W_imag_'+name][:] #imaginary component

    #flattening eigenvector arrays into 1D arrays
    if hasattr(Nmodes,'__len__'):
        W_real = tb.Flatten3D(W_real)
        W_imag = tb.Flatten3D(W_imag)

    else:
        W_real = tb.Flatten2D(W_real)
        W_imag = tb.Flatten2D(W_imag)

    W = np.zeros(len(W_real), dtype='complex')
    for i in range(len(W_real)):
        W[i] = complex(W_real[i], W_imag[i])

    return W

#==================================================================================================================================#
#
#==================================================================================================================================#
def DMD_Sort(DMD_Data):

    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        N_g = len(DMD_Data.dat.N_modes)
        for g in range(N_g):

            evl_p = g * DMD_Data.rank
            l = DMD_Data.dat.eval[evl_p : evl_p + DMD_Data.dat.N_modes[g]]
            el = DMD_Data.dat.expval[evl_p : evl_p + DMD_Data.dat.N_modes[g]]

            evc_p = g * DMD_Data.clen * DMD_Data.rank
            w = DMD_Data.dat.evec[evc_p : evc_p + DMD_Data.dat.N_modes[g] * DMD_Data.clen]

            for j in range(DMD_Data.dat.N_modes[g]):
                max = math.sqrt( l[j].real**2 + l[j].imag**2 )
                loc = j
                for k in range(j+1,DMD_Data.dat.N_modes[g]):
                    val = math.sqrt( l[k].real**2 + l[k].imag**2 )
                    if (val > max):
                        max = val
                        loc = k

                if (loc != j):
                    evl = l[j]
                    expvl = el[j]

                    p = j * DMD_Data.clen
                    evc = w[p : p + DMD_Data.clen]

                    l[j] = l[loc]
                    el[j] = el[loc]

                    p2 = loc * DMD_Data.clen
                    w[p : p + DMD_Data.clen] = w[p2 : p2 + DMD_Data.clen]

                    l[loc] = evl
                    el[loc] = expvl

                    w[p2 : p2 + DMD_Data.clen] = evc

            DMD_Data.dat.eval[evl_p : evl_p + DMD_Data.dat.N_modes[g]] = l
            DMD_Data.dat.expval[evl_p : evl_p + DMD_Data.dat.N_modes[g]] = el
            DMD_Data.dat.evec[evc_p : evc_p + DMD_Data.dat.N_modes[g] * DMD_Data.clen] = w

    else:
        for j in range(DMD_Data.dat.N_modes):
            max = math.sqrt( DMD_Data.dat.eval[j].real**2 + DMD_Data.dat.eval[j].imag**2 )
            loc = j
            for k in range(j+1,DMD_Data.dat.N_modes):
                val = math.sqrt( DMD_Data.dat.eval[k].real**2 + DMD_Data.dat.eval[k].imag**2 )
                if (val > max):
                    max = val
                    loc = k

            if (loc != j):
                evl = DMD_Data.dat.eval[j]
                expvl = DMD_Data.dat.expval[j]

                p = j * DMD_Data.clen
                evc = np.zeros(DMD_Data.clen, dtype='complex')
                for q in range(DMD_Data.clen): evc[q] = DMD_Data.dat.evec[p + q]

                DMD_Data.dat.eval[j] = DMD_Data.dat.eval[loc]
                DMD_Data.dat.expval[j] = DMD_Data.dat.expval[loc]

                p2 = loc * DMD_Data.clen
                DMD_Data.dat.evec[p : p + DMD_Data.clen] = DMD_Data.dat.evec[p2 : p2 + DMD_Data.clen]

                DMD_Data.dat.eval[loc] = evl
                DMD_Data.dat.expval[loc] = expvl

                DMD_Data.dat.evec[p2 : p2 + DMD_Data.clen] = evc

    return DMD_Data

#==================================================================================================================================#
#
#==================================================================================================================================#
def Coef_Calc(DMD_Data,init):

    b = []
    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        quit()

    else:
        p = 0
        w = np.reshape(DMD_Data.dat.evec, (DMD_Data.dat.N_modes, DMD_Data.clen))
        b = np.linalg.solve(w, init)

    return b

#==================================================================================================================================#
#
#==================================================================================================================================#
def Expand(DMD_Data, Coef, t):

    Expn = np.zeros(DMD_Data.clen)
    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        quit()

    else:
        p = 0
        for i in range(DMD_Data.dat.N_modes):
            w = DMD_Data.dat.evec[p : p + DMD_Data.clen]
            l = DMD_Data.dat.eval[i]
            Expn = Expn + Coef[i] * w * l**t

            p += DMD_Data.clen

    return Expn
