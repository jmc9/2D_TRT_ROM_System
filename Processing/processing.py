#!/usr/bin/env python3
#==================================================================================================================================#
#
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
import numpy as np
import math_tools as mt
import os

#==================================================================================================================================#
#
#==================================================================================================================================#
def cross_section(data,shape,axis):

    pa=0
    for a in axis:
        if a<0:break
        pa+=1

    n=len(shape)
    dim=shape[pa]

    p=0
    for i in range(pa):
        prod=axis[i]
        for j in range(i,pa): prod*=shape[j]
        p+=prod
    for i in range(pa+1,n):
        # print(i,axis[i])
        prod=axis[i]
        for j in range(i,pa): prod*=shape[j]
        p+=prod


    l=1
    if pa<(n-1):
        for s in shape[pa+1:]: l*=s

    z = data[p:(p+dim*l)]

    cs = [ x for x in z[::l] ]

    return cs

#==================================================================================================================================#
#
#==================================================================================================================================#
def errcalc(data):
    n = len(data[0])

    err=[]; ref=[]; rerr=[]
    for i in range(n):
        ref.append(data[0][i])
        err.append(data[0][i] - data[1][i])
        rerr.append(abs(err[i]/ref[i]))

    return err, rerr, ref

#==================================================================================================================================#
#
#==================================================================================================================================#
def data_compare(data,A):
    err, rerr, ref = errcalc(data)

    enorm_inf = mt.norm_inf(err)
    enorm2 = mt.norm2(err)
    enormL2 = mt.normL2(err,A)

    rnorm_inf = enorm_inf/mt.norm_inf(ref)
    rnorm2 = enorm2/mt.norm2(ref)
    rnormL2 = enormL2/mt.normL2(ref,A)

    return (enorm_inf,rnorm_inf,enorm2,rnorm2,enormL2,rnormL2)

#==================================================================================================================================#
#
#==================================================================================================================================#
def gdata_compare(data,Ag):
    n = len(data[0])

    err=[]; ref=[];
    for i in range(n):
        ref.append(data[0][i])
        err.append(data[0][i] - data[1][i])

    enorm_inf = mt.norm_inf(err)
    enorm2 = mt.norm2(err)
    enormL2 = mt.normL2(err,Ag)

    rnorm_inf = enorm_inf/mt.norm_inf(ref)
    rnorm2 = enorm2/mt.norm2(ref)
    rnormL2 = enormL2/mt.normL2(ref,Ag)

    return (enorm_inf,rnorm_inf,enorm2,rnorm2,enormL2,rnormL2)

#==================================================================================================================================#
#
#==================================================================================================================================#
def get_data(dset,inp_flags,t,proc_qdf_mg,proc_flux_mg,proc_edens_mg):
    Temp = dset['Temperature'][t][:]; Temp = dat_flatten(Temp)

    E_avg=[]; E_edgV=[]; E_edgH=[]
    if inp_flags[12] == 1:
        E_avg = dset['E_avg_Grey'][t][:]; E_avg = dat_flatten(E_avg)
        E_edgV = dset['E_edgV_Grey'][t][:]; E_edgV = dat_flatten(E_edgV)
        E_edgH = dset['E_edgH_Grey'][t][:]; E_edgH = dat_flatten(E_edgH)

    E_avg_MGQD=[]; E_edgV_MGQD=[]; E_edgH_MGQD=[]
    if inp_flags[9] == 1:
        E_avg_MGQD = dset['E_avg_MGQD'][t][:]; E_avg_MGQD = dat_flatten(E_avg_MGQD)
        E_edgV_MGQD = dset['E_edgV_MGQD'][t][:]; E_edgV_MGQD = dat_flatten(E_edgV_MGQD)
        E_edgH_MGQD = dset['E_edgH_MGQD'][t][:]; E_edgH_MGQD = dat_flatten(E_edgH_MGQD)

    E_avg_HO=[]; E_edgV_HO=[]; E_edgH_HO=[]
    if inp_flags[5] == 1:
        E_avg_HO = dset['E_avg_HO'][t][:]; E_avg_HO = dat_flatten(E_avg_HO)
        E_edgV_HO = dset['E_edgV_HO'][t][:]; E_edgV_HO = dat_flatten(E_edgV_HO)
        E_edgH_HO = dset['E_edgH_HO'][t][:]; E_edgH_HO = dat_flatten(E_edgH_HO)

    Fx_edgV=[]; Fy_edgH=[]
    if inp_flags[13] == 1:
        Fx_edgV = dset['Fx_edgV_Grey'][t][:]; Fx_edgV = dat_flatten(Fx_edgV)
        Fy_edgH = dset['Fy_edgH_Grey'][t][:]; Fy_edgH = dat_flatten(Fy_edgH)

    Fx_edgV_MGQD=[]; Fy_edgH_MGQD=[]
    if inp_flags[11] == 1:
        Fx_edgV_MGQD = dset['Fx_edgV_MGQD'][t][:]; Fx_edgV_MGQD = dat_flatten(Fx_edgV_MGQD)
        Fy_edgH_MGQD = dset['Fy_edgH_MGQD'][t][:]; Fy_edgH_MGQD = dat_flatten(Fy_edgH_MGQD)

    Fx_edgV_HO=[]; Fy_edgH_HO=[]
    if inp_flags[7] == 1:
        Fx_edgV_HO = dset['Fx_edgV_HO'][t][:]; Fx_edgV_HO = dat_flatten(Fx_edgV_HO)
        Fy_edgH_HO = dset['Fy_edgH_HO'][t][:]; Fy_edgH_HO = dat_flatten(Fy_edgH_HO)

    Eg_avg=[]; Eg_edgV=[]; Eg_edgH=[]
    if ((inp_flags[8] == 1) and (proc_edens_mg==1)):
        Eg_avg = dset['Eg_avg_MGQD'][t][:]; Eg_avg = gdat_flatten(Eg_avg)
        Eg_edgV = dset['Eg_edgV_MGQD'][t][:]; Eg_edgV = gdat_flatten(Eg_edgV)
        Eg_edgH = dset['Eg_edgH_MGQD'][t][:]; Eg_edgH = gdat_flatten(Eg_edgH)

    Eg_avg_HO=[]; Eg_edgV_HO=[]; Eg_edgH_HO=[]
    if ((inp_flags[4] == 1) and (proc_edens_mg==1)):
        Eg_avg_HO = dset['Eg_avg_HO'][t][:]; Eg_avg_HO = gdat_flatten(Eg_avg_HO)
        Eg_edgV_HO = dset['Eg_edgV_HO'][t][:]; Eg_edgV_HO = gdat_flatten(Eg_edgV_HO)
        Eg_edgH_HO = dset['Eg_edgH_HO'][t][:]; Eg_edgH_HO = gdat_flatten(Eg_edgH_HO)

    Fxg_edgV=[]; Fyg_edgH=[]
    if ((inp_flags[10] == 1) and (proc_flux_mg==1)):
        Fxg_edgV = dset['Fxg_edgV_MGQD'][t][:]; Fxg_edgV = gdat_flatten(Fxg_edgV)
        Fyg_edgH = dset['Fyg_edgH_MGQD'][t][:]; Fyg_edgH = gdat_flatten(Fyg_edgH)

    Fxg_edgV_HO=[]; Fyg_edgH_HO=[]
    if ((inp_flags[6] == 1) and (proc_flux_mg==1)):
        Fxg_edgV_HO = dset['Fxg_edgV_HO'][t][:]; Fxg_edgV_HO = gdat_flatten(Fxg_edgV_HO)
        Fyg_edgH_HO = dset['Fyg_edgH_HO'][t][:]; Fyg_edgH_HO = gdat_flatten(Fyg_edgH_HO)

    fg_avg_xx=[]; fg_avg_xy=[]; fg_avg_yy=[]
    fg_edgV_xx=[]; fg_edgV_xy=[]
    fg_edgH_yy=[]; fg_edgH_xy=[]
    if ((inp_flags[2] == 1) and (proc_qdf_mg==1)):
        fg_avg_xx = dset['fg_avg_xx'][t][:]; fg_avg_xx = gdat_flatten(fg_avg_xx)
        fg_avg_xy = dset['fg_avg_xy'][t][:]; fg_avg_xy = gdat_flatten(fg_avg_xy)
        fg_avg_yy = dset['fg_avg_yy'][t][:]; fg_avg_yy = gdat_flatten(fg_avg_yy)

        fg_edgV_xx = dset['fg_edgV_xx'][t][:]; fg_edgV_xx = gdat_flatten(fg_edgV_xx)
        fg_edgV_xy = dset['fg_edgV_xy'][t][:]; fg_edgV_xy = gdat_flatten(fg_edgV_xy)

        fg_edgH_yy = dset['fg_edgH_yy'][t][:]; fg_edgH_yy = gdat_flatten(fg_edgH_yy)
        fg_edgH_xy = dset['fg_edgH_xy'][t][:]; fg_edgH_xy = gdat_flatten(fg_edgH_xy)

    T_xWvSpd = dset['Temp_xWave_Speed'][t]
    E_xWvSpd = dset['E_xWave_Speed'][t]

    return (Temp,E_avg,E_edgV,E_edgH,E_avg_MGQD,E_edgV_MGQD,E_edgH_MGQD,E_avg_HO,E_edgV_HO,E_edgH_HO,Fx_edgV,\
    Fy_edgH,Fx_edgV_MGQD,Fy_edgH_MGQD,Fx_edgV_HO,Fy_edgH_HO,Eg_avg,Eg_edgV,Eg_edgH,Eg_avg_HO,Eg_edgV_HO,\
    Eg_edgH_HO,Fxg_edgV,Fyg_edgH,Fxg_edgV_HO,Fyg_edgH_HO,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,\
    fg_edgH_yy,fg_edgH_xy,T_xWvSpd,E_xWvSpd)

def gdat_flatten(gdat):
    gdat = [x for G_list in gdat for Y_list in G_list for x in Y_list]
    return gdat

def dat_flatten(dat):
    dat = [x for Y_list in dat for x in Y_list]
    return dat

#==================================================================================================================================#
#
#==================================================================================================================================#
def init_arr():
    Temp=[[],[]]
    E_avg=[[],[]]; E_edgV=[[],[]]; E_edgH=[[],[]]
    E_avg_MGQD=[[],[]]; E_edgV_MGQD=[[],[]]; E_edgH_MGQD=[[],[]]
    E_avg_HO=[[],[]]; E_edgV_HO=[[],[]]; E_edgH_HO=[[],[]]
    Fx_edgV=[[],[]]; Fy_edgH=[[],[]]
    Fx_edgV_MGQD=[[],[]]; Fy_edgH_MGQD=[[],[]]
    Fx_edgV_HO=[[],[]]; Fy_edgH_HO=[[],[]]
    Eg_avg=[[],[]]; Eg_edgV=[[],[]]; Eg_edgH=[[],[]]
    Eg_avg_HO=[[],[]]; Eg_edgV_HO=[[],[]]; Eg_edgH_HO=[[],[]]
    Fxg_edgV=[[],[]]; Fyg_edgH=[[],[]]
    Fxg_edgV_HO=[[],[]]; Fyg_edgH_HO=[[],[]]
    fg_avg_xx=[[],[]]; fg_avg_xy=[[],[]]; fg_avg_yy=[[],[]]
    fg_edgV_xx=[[],[]]; fg_edgV_xy=[[],[]]
    fg_edgH_yy=[[],[]]; fg_edgH_xy=[[],[]]

    return (Temp,E_avg,E_edgV,E_edgH,E_avg_MGQD,E_edgV_MGQD,E_edgH_MGQD,E_avg_HO,E_edgV_HO,E_edgH_HO,Fx_edgV,\
    Fy_edgH,Fx_edgV_MGQD,Fy_edgH_MGQD,Fx_edgV_HO,Fy_edgH_HO,Eg_avg,Eg_edgV,Eg_edgH,Eg_avg_HO,Eg_edgV_HO,\
    Eg_edgH_HO,Fxg_edgV,Fyg_edgH,Fxg_edgV_HO,Fyg_edgH_HO,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,\
    fg_edgH_yy,fg_edgH_xy)

#==================================================================================================================================#
#
#==================================================================================================================================#
def A_calcs(A,N_g,N_y,N_x):

    #----------A_EdgV----------
    #contains the 'areas' of quantities residing on vertical (y=const) cell faces

    A_EdgV = np.zeros([N_y*(N_x+1)])
    p1 = 0
    p2 = 0
    for j in range(N_y):
        A_EdgV[p1] = A[p2]/2.

        for i in range(N_x-1):
            p1 = p1 + 1
            A_EdgV[p1] = (A[p2]+A[p2+1])/2.
            p2 = p2 + 1

        p1 = p1 + 1
        A_EdgV[p1] = A[p2]/2.

    #----------A_EdgH----------
    #contains the 'areas' of quantities residing on horizontal (x=const) cell faces

    A_EdgH = np.zeros([(N_y+1)*N_x])
    p1 = 0
    p2 = 0
    for i in range(N_x):
        A_EdgH[p1] = A[p2]/2.
        p1 = p1 + 1
        p2 = p2 + 1
    p2 = 0
    for j in range(N_y-1):
        for i in range(N_x):
            A_EdgV[p1] = (A[p2]+A[p2+N_x])/2.
            p1 = p1 + 1
            p2 = p2 + 1
    p2 = p2 - N_x
    for i in range(N_x):
        A_EdgH[p1] = A[p2]/2.
        p1 = p1 + 1
        p2 = p2 + 1

    #----------Ag----------
    #contains the group/spatial weights of multigroup cell-averaged quantities

    Ag = np.zeros([N_g*N_y*N_x])
    p1 = 0
    for g in range(N_g):
        p2 = 0
        for i in range(N_y*N_x):
            Ag[p1] = A[p2]
            p1 = p1 + 1
            p2 = p2 + 1

    #----------Ag_EdgV----------
    #contains the group/spatial weights of multigroup quantities residing on vertical (y=const) cell faces

    Ag_EdgV = np.zeros([N_g*N_y*(N_x+1)])
    p1 = 0
    for g in range(N_g):
        p2 = 0
        for i in range(N_y*(N_x+1)):
            Ag_EdgV[p1] = A_EdgV[p2]
            p1 = p1 + 1
            p2 = p2 + 1

    #----------Ag_EdgH----------
    #contains the group/spatial weights of multigroup quantities residing on horizontal (x=const) cell faces

    Ag_EdgH = np.zeros([N_g*(N_y+1)*N_x])
    p1 = 0
    for g in range(N_g):
        p2 = 0
        for i in range((N_y+1)*N_x):
            Ag_EdgH[p1] = A_EdgH[p2]
            p1 = p1 + 1
            p2 = p2 + 1

    return (A_EdgV, A_EdgH, Ag, Ag_EdgV, Ag_EdgH)

#==================================================================================================================================#
#
#==================================================================================================================================#
def rank_collect(dsets,nsets,casename,N_g,proc_dir):

    rankfile = open('Ranks.txt','w')

    rrank_BCg = []
    rrank_fg_avg_xx = []
    rrank_fg_edgV_xx = []
    rrank_fg_avg_yy = []
    rrank_fg_edgH_yy = []
    rrank_fg_edgV_xy = []
    rrank_fg_edgH_xy = []
    dset_times = []

    for i in range(nsets-1):
        r_BCg = []
        r_avg_xx = []
        r_edgV_xx = []
        r_avg_yy = []
        r_edgH_yy = []
        r_edgV_xy = []
        r_edgH_xy = []

        rankfile.write('Dataset: '+casename[i+1]+'\n')
        rankfile.write('--------------------------------------------------\n')

        set = dsets[i+1]
        run_type = set.getncattr("run_type")
        if run_type == "mg_pod":
            gsum = set['PODgsum'][:]
            dset_times.append(set.getncattr("dset_times"))
        elif run_type == "mg_dmd":
            gsum = set['DMDgsum'][:]
            dset_times.append(set.getncattr("dset_times"))
        nbases = len(gsum)

        rg = np.full(nbases, 1, dtype=int)
        for j in range(nbases):
            if (gsum[j] == 1):
                r_BCg.append(    [ set['rrank_BCg'][j][0]        ])
                r_avg_xx.append( [ set['rrank_fg_avg_xx'][j][0]  ])
                r_edgV_xx.append([ set['rrank_fg_edgV_xx'][j][0] ])
                r_avg_yy.append( [ set['rrank_fg_avg_yy'][j][0]  ])
                r_edgH_yy.append([ set['rrank_fg_edgH_yy'][j][0] ])
                r_edgV_xy.append([ set['rrank_fg_edgV_xy'][j][0] ])
                r_edgH_xy.append([ set['rrank_fg_edgH_xy'][j][0] ])

            else:
                rg[j] = N_g
                r_BCg.append(    set['rrank_BCg'][j][:])
                r_avg_xx.append( set['rrank_fg_avg_xx'][j][:])
                r_edgV_xx.append(set['rrank_fg_edgV_xx'][j][:])
                r_avg_yy.append( set['rrank_fg_avg_yy'][j][:])
                r_edgH_yy.append(set['rrank_fg_edgH_yy'][j][:])
                r_edgV_xy.append(set['rrank_fg_edgV_xy'][j][:])
                r_edgH_xy.append(set['rrank_fg_edgH_xy'][j][:])

        rank_write(rankfile, 'BCg:', nbases, rg, r_BCg)
        rank_write(rankfile, 'fg_avg_xx:', nbases, rg, r_avg_xx)
        rank_write(rankfile, 'fg_edgV_xx:', nbases, rg, r_edgV_xx)
        rank_write(rankfile, 'fg_avg_yy:', nbases, rg, r_avg_yy)
        rank_write(rankfile, 'fg_edgH_yy:', nbases, rg, r_edgH_yy)
        rank_write(rankfile, 'fg_edgV_xy:', nbases, rg, r_edgV_xy)
        rank_write(rankfile, 'fg_edgH_xy:', nbases, rg, r_edgH_xy)

        rrank_BCg.append( r_BCg )
        rrank_fg_avg_xx.append(  r_avg_xx )
        rrank_fg_edgV_xx.append( r_edgV_xx )
        rrank_fg_avg_yy.append(  r_avg_yy )
        rrank_fg_edgH_yy.append( r_edgH_yy )
        rrank_fg_edgV_xy.append( r_edgV_xy )
        rrank_fg_edgH_xy.append( r_edgH_xy )

        rankfile.write('\n')
    rankfile.close()
    os.system('mv Ranks.txt '+proc_dir)

    time_intervals = [0.]
    all_dset_times = [ t for times in dset_times for t in times ]
    ntimes = len(all_dset_times)
    for i in range(ntimes):
        t = max(all_dset_times)
        for j in range(ntimes):
            if ((all_dset_times[j] > time_intervals[i])and(all_dset_times[j] < t)):
                t = all_dset_times[j]
        time_intervals.append(t)
        if t==max(all_dset_times): break

    ranks = []
    ntimes = len(time_intervals)
    for i in range(ntimes-1):
        r_BCg = []
        r_avg_xx = []
        r_edgV_xx = []
        r_avg_yy = []
        r_edgH_yy = []
        r_edgV_xy = []
        r_edgH_xy = []
        for j in range(nsets-1):
            k = find_tloc(dset_times[j], time_intervals[i], time_intervals[i+1])
            r_BCg.append(    rrank_BCg[j][k][0] )
            r_avg_xx.append( rrank_fg_avg_xx[j][k][0]  )
            r_edgV_xx.append(rrank_fg_edgV_xx[j][k][0] )
            r_avg_yy.append( rrank_fg_avg_yy[j][k][0]  )
            r_edgH_yy.append(rrank_fg_edgH_yy[j][k][0] )
            r_edgV_xy.append(rrank_fg_edgV_xy[j][k][0] )
            r_edgH_xy.append(rrank_fg_edgH_xy[j][k][0] )
        ranks.append( [ r_BCg, r_avg_xx, r_edgV_xx, r_avg_yy, r_edgH_yy, r_edgV_xy, r_edgH_xy ] )

    return (time_intervals, ranks, gsum)

#==================================================================================================================================#
#
#==================================================================================================================================#
def rank_write(rankfile, name, nbases, rg, rank):
    rankfile.write(name)
    rankfile.write('\n')
    for j in range(nbases):
        rankfile.write(' ({}) '.format(j+1))
        for g in range(rg[j]): rankfile.write('  {}'.format(rank[j][g]))
        rankfile.write('\n')

#==================================================================================================================================#
#
#==================================================================================================================================#
def find_tloc(times, t_low, t_up):
    loc = len(times)
    for i in range(loc-1):
        # print(t, times[i], times[i+1])
        if ((t_low>=times[i])and(t_up<=times[i+1])):
            loc = i
            # print(loc)
            break
    return loc
