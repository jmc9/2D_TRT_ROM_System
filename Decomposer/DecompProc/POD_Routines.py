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

from Classes import POD

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Dims(dset, name):
    dims = dset['U_'+name].dimensions
    ndims = len(dims)

    rank = dset.dimensions[dims[0]].size
    clen = dset.dimensions[dims[ndims-1]].size

    return (rank, clen)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Trunc_SVD(S, eps):
    n = len(S)
    sum1 = 0.
    for i in range(n):
        sum1 += S[i]**2

    sum2 = 0.
    for i in reversed(range(n)):
        sum2 += S[i]**2
        if ( math.sqrt(sum2 / sum1) >= eps ):
            return i + 1

    return -1

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_S(dset, name, eps):
    S_full = dset['S_{}'.format(name)][:]
    shp = np.shape(S_full)
    n = len(shp)

    if n == 1:
        t_rank = Trunc_SVD(S_full, eps)
        S = S_full[0 : t_rank]
    elif n == 2:
        t_rank = [Trunc_SVD(S_full[g], eps) for g in range(shp[0])]
        S = [s for (Sg, r) in zip(S_full, t_rank) for s in Sg[0 : r]]
    else:
        print('POD_Routines/Read_S - len( np.shape(S) ) > 2 !')
        quit()

    return S, t_rank

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_U(dset, name, t_rank):
    U_full = dset['U_{}'.format(name)][:]
    shp = np.shape(U_full)
    n = len(shp)

    if n == 2:
        U = U_full[0 : t_rank]
        U = [ui for u in U for ui in u]
    elif n == 3:
        U = [ui for (Ug, r) in zip(U_full, t_rank) for u in Ug[0 : r] for ui in u]
    else:
        print('POD_Routines/Read_U - len( np.shape(U) ) not in [2, 3] !')
        quit()

    return U

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_V(dset, name, t_rank):
    Vt_full = dset['Vt_{}'.format(name)][:]
    shp = np.shape(Vt_full)
    n = len(shp)

    if n == 2:
        V_full = np.transpose(Vt_full)
        V = V_full[0 : t_rank]
        V = [vi for v in V for vi in v]
    elif n == 3:
        V_full = np.transpose(Vt_full, (0,2,1))
        V = [vi for (Vg, r) in zip(V_full, t_rank) for v in Vg[0 : r] for vi in v]
    else:
        print('POD_Routines/Read_Vt - len( np.shape(V) ) not in [2, 3] !')
        quit()

    return V

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_C(dset, name):
    C = dset['C_{}'.format(name)][:]
    shp = np.shape(C)
    n = len(shp)

    if n == 1:
        pass
    elif n == 2:
        C = [c for Cg in C for c in Cg]
    else:
        print('POD_Routines/Read_C - len( np.shape(C) ) > 2 !')
        quit()

    return C

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_POD(dset, name, clen, rank, eps):

    S, t_rank = Read_S(dset, name, eps)
    U = Read_U(dset, name, t_rank)
    V = Read_V(dset, name, t_rank)
    C = Read_C(dset, name)

    POD_dat = POD(U, S, V, C, t_rank)

    return POD_dat

#==================================================================================================================================#
#
#==================================================================================================================================#
def Expand(POD_Data, t):

    if hasattr(POD_Data.dat.t_rank,'__len__'):
        nr = len(POD_Data.dat.t_rank)
        Expn = np.zeros(nr * POD_Data.clen)
        Nt = int( len(POD_Data.dat.V) / sum(POD_Data.dat.t_rank) )

        p1 = 0
        p2 = 0
        p3 = 0
        p4 = t
        for n in range(nr):
            for i in range(POD_Data.dat.t_rank[n]):
                for j in range(POD_Data.clen):
                    Expn[p1] += POD_Data.dat.S[p3] * POD_Data.dat.U[p2] * POD_Data.dat.V[p4]
                    p1 += 1
                    p2 += 1
                p1 -= POD_Data.clen
                p3 += 1
                p4 += Nt
            p1 += POD_Data.clen

    else:
        Expn = np.zeros(POD_Data.clen)
        Nt = int( len(POD_Data.dat.V) / POD_Data.dat.t_rank )

        p1 = 0
        p2 = t
        for i in range(POD_Data.dat.t_rank):
            for j in range(POD_Data.clen):
                Expn[j] += POD_Data.dat.S[i] * POD_Data.dat.U[p1] * POD_Data.dat.V[p2]
                p1 += 1
            p2 += Nt

    Expn += POD_Data.dat.C

    return Expn

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_POD(POD_Data, dir, plt_modes, svals=True, uvecs=True, vvecs=True):
    if svals:
        sval_dir = dir+'/Svals'
        tb.dirset(sval_dir)
        Plot_Svals(POD_Data, sval_dir)

    if uvecs:
        uvec_dir = dir+'/Uvecs'
        tb.dirset(uvec_dir)
        Plot_Uvecs(POD_Data, uvec_dir, plt_modes)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_Svals(POD_Data, drop):
    if hasattr(POD_Data.dat.t_rank,'__len__'):
        p = 0
        for i in range(len(POD_Data.dat.t_rank)):
            xp = np.arange(POD_Data.dat.t_rank[i])+1
            pdata = POD_Data.dat.S[p : p + POD_Data.dat.t_rank[i]]

            pltr.lineplot('Svals_g{}'.format(i+1), xp, pdata, drop, yscale='log', line='k.', marker='.')

            p += POD_Data.dat.t_rank[i]

    else:
        xp = np.arange(POD_Data.dat.t_rank)+1
        pdata = POD_Data.dat.S
        pltr.lineplot('Svals', xp, pdata, drop, yscale='log', line='k.', marker='.')

#==================================================================================================================================#
#
#==================================================================================================================================#
def Plot_Uvecs(POD_Data, dir, plt_modes):
    N_grids = len(POD_Data.grids)
    N_dims  = len(POD_Data.dims)

    if hasattr(POD_Data.dat.t_rank,'__len__'):
        N_g = POD_Data.dims[1].len
        for g in range(N_g):
            for m in plt_modes[g]:
                p = sum(POD_Data.dat.t_rank[:g]) + m

                p = POD_Data.clen * sum(POD_Data.dat.t_rank[:g]) + m * POD_Data.clen
                u = POD_Data.dat.U[p : p + POD_Data.clen]
                pltname = '{}_m{}_g{}'.format(POD_Data.name,m+1,g+1)

                if (N_grids == 2): #
                    xp = POD_Data.grids[N_grids-1].crds

                    pltr.lineplot(pltname, xp, u, dir)

                elif (N_grids == 3): #
                    N_x = POD_Data.dims[N_dims-1].len
                    N_y = POD_Data.dims[N_dims-2].len

                    xp = tb.Grid_Edg_Gen(POD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(POD_Data.grids[N_grids-2])

                    pltr.plot_heatmap2D(u, pltname, xp, yp, dir, N_y, N_x)

    else:
        for m in plt_modes:
            p = m * POD_Data.clen
            u = POD_Data.dat.U[p : p + POD_Data.clen]
            pltname = '{}_m{}'.format(POD_Data.name,m+1)

            if (N_grids == N_dims): #no 'group' dimensions
                if (N_grids == 2): #simple 2D data
                    xp = POD_Data.grids[N_grids-1].crds

                    pltr.lineplot(pltname, xp, u, dir)

                elif (N_grids == 3): #
                    N_x = POD_Data.dims[N_dims-1].len
                    N_y = POD_Data.dims[N_dims-2].len

                    xp = tb.Grid_Edg_Gen(POD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(POD_Data.grids[N_grids-2])

                    pltr.plot_heatmap2D(u, pltname, xp, yp, dir, N_y, N_x)


            elif (N_grids == N_dims - 1): #only 1 'group' dimension exists
                if (N_grids == 2): #
                    xp = POD_Data.grids[N_grids-1].crds
                    N_g = POD_Data.dims[N_dims-2].len

                    pltr.lineplot_grouped(pltname, xp, N_g, u, dir)

                elif (N_grids == 3): #
                    N_x = POD_Data.dims[N_dims-1].len
                    N_y = POD_Data.dims[N_dims-2].len
                    N_z = POD_Data.dims[N_dims-3].len

                    xp = tb.Grid_Edg_Gen(POD_Data.grids[N_grids-1])
                    yp = tb.Grid_Edg_Gen(POD_Data.grids[N_grids-2])

                    pltr.plot_heatmap3D(u, pltname, xp, yp, dir, N_z, N_y, N_x, zlab='g')

            else:
                print('Plot_DMD_evecs can only handle data with a single group dimension right now!')
                quit()
