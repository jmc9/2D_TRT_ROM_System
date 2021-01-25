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
import copy

#local tools
import Plotters as pltr
import ToolBox as tb

from Classes import DMD

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD(DMD_Data,dir,plt_modes=[0,1,2],evecs=True,evals=True):

    dat_ = copy.deepcopy(DMD_Data) #have to copy the DMD data object like this to avoid copying to dat_ by reference
    DMD_Sort(dat_) #sorting eigenvalues/eigenvectors by magnitude of eigenvalues

    if evecs:
        evec_dir = dir+'/Eigenvectors'
        tb.dirset(evec_dir)
        Plot_DMD_evecs(dat_,plt_modes,evec_dir)

    if evals:
        eval_dir = dir+'/Eigenvalues'
        tb.dirset(eval_dir)
        Plot_DMD_evals(dat_,eval_dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_DMD_evecs(DMD_Data,plt_modes,dir):
    N_grids = len(DMD_Data.grids)
    N_dims  = len(DMD_Data.dims)

    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        N_g = DMD_Data.dims[1].len
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
                    N_x = DMD_Data.dims[N_dims-1].len
                    N_y = DMD_Data.dims[N_dims-2].len

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
                    N_x = DMD_Data.dims[N_dims-1].len
                    N_y = DMD_Data.dims[N_dims-2].len

                    xp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(DMD_Data.grids[N_grids-2])

                    pltr.plot_heatmap2D(wr, pltname_r, xp, yp, dir, N_y, N_x, title=title_r)
                    if li != 0.: pltr.plot_heatmap2D(wi, pltname_i, xp, yp, dir, N_y, N_x, title=title_i)


            elif (N_grids == N_dims - 1): #only 1 'group' dimension exists
                if (N_grids == 2): #
                    print('Plot_DMD_evecs doesn not support N_grids = 2 yet for groupwise decompositions')
                    quit()

                elif (N_grids == 3): #
                    N_x = DMD_Data.dims[N_dims-1].len
                    N_y = DMD_Data.dims[N_dims-2].len
                    N_z = DMD_Data.dims[N_dims-3].len

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
        N_g = DMD_Data.dims[1].len
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
def Read_DMD(dset,name,clen,rank):

    N_modes = getattr(dset,'N_modes_'+name)

    Evals = Read_Evals(dset, name, N_modes)
    Evecs = Read_Evecs(dset, name, N_modes, clen)
    Expvals = Read_Expvals(dset, name, N_modes)
    r_Evecs = Read_r_Evecs(dset, name, N_modes, rank)
    Uvecs = Read_Vecs(dset, name, N_modes, clen)

    DMD_dat = DMD(N_modes, Evals, Evecs, Expvals, Uvecs, r_Evecs)

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
        # L_real = tb.Flatten2D(L_real)
        # L_imag = tb.Flatten2D(L_imag)
        L = np.zeros(sum(Nmodes), dtype='complex')
        p1 = 0
        for n in range(len(Nmodes)):
            for i in range(Nmodes[n]):
                L[p1] = complex(L_real[n][i], L_imag[n][i])
                p1 += 1

    else:
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
def Read_Evecs(dset,name,Nmodes,vlen,evcn='W',tlen=0):
    # #reading eigenvectors from dataset
    W_real = Read_Vecs(dset,name,Nmodes,vlen,'{}_real'.format(evcn),tlen)
    W_imag = Read_Vecs(dset,name,Nmodes,vlen,'{}_imag'.format(evcn),tlen)

    W = np.zeros(len(W_real), dtype='complex')
    for i in range(len(W_real)):
        W[i] = complex(W_real[i], W_imag[i])

    return W

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_r_Evecs(dset,name,Nmodes,rank):
    rW = Read_Evecs(dset,name,Nmodes,rank,'Wtild',Nmodes)

    return rW

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Vecs(dset,name,Nmodes,vlen,uvcn='Umat',tlen=0):
    if tlen == 0: tlen = vlen

    #reading eigenvectors from dataset
    # U = np.array(dset['{}_'.format(uvcn)+name])
    U_ = np.array(dset['{}_'.format(uvcn)+name])

    #flattening eigenvector arrays into 1D arrays
    if hasattr(Nmodes,'__len__'):
        # U_ = tb.Flatten3D(U)
        U = np.zeros(sum(Nmodes) * tlen) #creating vector of true size

        p1 = 0
        for n in range(len(Nmodes)):
            for i in range(Nmodes[n]):
                for j in range(tlen):
                    U[p1] = U_[n][i][j]
                    p1 += 1


    else:
        # U_ = tb.Flatten2D(U)
        U = np.zeros(Nmodes * tlen) #creating vector of true size

        #copying flattened vector back into U, truncating columns after tlen elements
        p1 = 0
        p2 = 0
        k = vlen - tlen
        for i in range(Nmodes):
            for j in range(tlen):
                U[p1] = U_[i][j]
                # U[p1] = U_[p2]
                p1 += 1
            #     p2 += 1
            # p2 += k

    return U

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
                    evc = np.zeros(DMD_Data.clen, dtype='complex')
                    for q in range(DMD_Data.clen): evc[q] = DMD_Data.dat.evec[p + q]
                    # evc = w[p : p + DMD_Data.clen]

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
        nm = len(DMD_Data.dat.N_modes)
        for n in range(nm):
            p1 = n * DMD_Data.dat.N_modes[n] * DMD_Data.dat.N_modes[n]
            p2 = p1 + DMD_Data.dat.N_modes[n] * DMD_Data.dat.N_modes[n]
            w = np.transpose( np.reshape(DMD_Data.dat.r_evec[p1 : p2], (DMD_Data.dat.N_modes[n], DMD_Data.dat.N_modes[n]) ) )

            p1 = n * DMD_Data.dat.N_modes[n] * DMD_Data.clen
            p2 = p1 + DMD_Data.dat.N_modes[n] * DMD_Data.clen
            ut = np.reshape(DMD_Data.dat.uvec[p1 : p2], (DMD_Data.dat.N_modes[n], DMD_Data.clen))
            ut = np.array(ut, dtype='complex')

            p1 = n * DMD_Data.clen
            p2 = p1 + DMD_Data.clen
            init_ = np.array(init[p1 : p2], dtype='complex')

            b.append( np.linalg.solve(w, np.dot(ut, init_) ) )

    else:
        w = np.transpose( np.reshape(DMD_Data.dat.r_evec, (DMD_Data.dat.N_modes, DMD_Data.dat.N_modes) ) )
        ut = np.reshape(DMD_Data.dat.uvec, (DMD_Data.dat.N_modes, DMD_Data.clen))
        ut = np.array(ut, dtype='complex')
        init_ = np.array(init, dtype='complex')
        b = np.linalg.solve(w, np.dot(ut, init_) )

    return b

#==================================================================================================================================#
#
#==================================================================================================================================#
def Expand(DMD_Data, Coef, t):

    # Expn = np.zeros(DMD_Data.clen, dtype='complex')
    if hasattr(DMD_Data.dat.N_modes,'__len__'):
        # quit()
        nm = len(DMD_Data.dat.N_modes)
        Expn = np.zeros(nm * DMD_Data.clen, dtype='complex')

        p1 = 0
        p2 = 0
        p3 = 0
        for n in range(nm):
            for i in range(DMD_Data.dat.N_modes):
                l = DMD_Data.dat.eval[p1]
                for j in range(DMD_Data.clen):
                    Expn[p3] = Expn[p3] + Coef[p1] * DMD_Data.dat.evec[p3] * l**t
                    p2 += 1
                    p3 += 1
                p3 -= DMD_Data.clen
                p1 += 1

    else:
        Expn = np.zeros(DMD_Data.clen, dtype='complex')

        p = 0
        for i in range(DMD_Data.dat.N_modes):
            l = DMD_Data.dat.eval[i]
            for j in range(DMD_Data.clen):
                Expn[j] = Expn[j] + Coef[i] * DMD_Data.dat.evec[p] * l**t
                p += 1

    return Expn
