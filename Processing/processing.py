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
def data_compare(data,A):
    n = len(data[0])

    err=[]; ref=[]
    for i in range(n):
        ref.append(data[0][i])
        err.append(data[0][i] - data[1][i])

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
def get_data(dset,inp_flags,t):
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
    if inp_flags[8] == 1:
        Eg_avg = dset['Eg_avg_MGQD'][t][:]; Eg_avg = gdat_flatten(Eg_avg)
        Eg_edgV = dset['Eg_edgV_MGQD'][t][:]; Eg_edgV = gdat_flatten(Eg_edgV)
        Eg_edgH = dset['Eg_edgH_MGQD'][t][:]; Eg_edgH = gdat_flatten(Eg_edgH)

    Eg_avg_HO=[]; Eg_edgV_HO=[]; Eg_edgH_HO=[]
    if inp_flags[4] == 1:
        Eg_avg_HO = dset['Eg_avg_HO'][t][:]; Eg_avg_HO = gdat_flatten(Eg_avg_HO)
        Eg_edgV_HO = dset['Eg_edgV_HO'][t][:]; Eg_edgV_HO = gdat_flatten(Eg_edgV_HO)
        Eg_edgH_HO = dset['Eg_edgH_HO'][t][:]; Eg_edgH_HO = gdat_flatten(Eg_edgH_HO)

    Fxg_edgV=[]; Fyg_edgH=[]
    if inp_flags[10] == 1:
        Fxg_edgV = dset['Fxg_edgV_MGQD'][t][:]; Fxg_edgV = gdat_flatten(Fxg_edgV)
        Fyg_edgH = dset['Fyg_edgH_MGQD'][t][:]; Fyg_edgH = gdat_flatten(Fyg_edgH)

    Fxg_edgV_HO=[]; Fyg_edgH_HO=[]
    if inp_flags[6] == 1:
        Fxg_edgV_HO = dset['Fxg_edgV_HO'][t][:]; Fxg_edgV_HO = gdat_flatten(Fxg_edgV_HO)
        Fyg_edgH_HO = dset['Fyg_edgH_HO'][t][:]; Fyg_edgH_HO = gdat_flatten(Fyg_edgH_HO)

    fg_avg_xx=[]; fg_avg_xy=[]; fg_avg_yy=[]
    fg_edgV_xx=[]; fg_edgV_xy=[]
    fg_edgH_yy=[]; fg_edgH_xy=[]
    if inp_flags[2] == 1:
        fg_avg_xx = dset['fg_avg_xx'][t][:]; fg_avg_xx = gdat_flatten(fg_avg_xx)
        fg_avg_xy = dset['fg_avg_xy'][t][:]; fg_avg_xy = gdat_flatten(fg_avg_xy)
        fg_avg_yy = dset['fg_avg_yy'][t][:]; fg_avg_yy = gdat_flatten(fg_avg_yy)

        fg_edgV_xx = dset['fg_edgV_xx'][t][:]; fg_edgV_xx = gdat_flatten(fg_edgV_xx)
        fg_edgV_xy = dset['fg_edgV_xy'][t][:]; fg_edgV_xy = gdat_flatten(fg_edgV_xy)

        fg_edgH_yy = dset['fg_edgH_yy'][t][:]; fg_edgH_yy = gdat_flatten(fg_edgH_yy)
        fg_edgH_xy = dset['fg_edgH_xy'][t][:]; fg_edgH_xy = gdat_flatten(fg_edgH_xy)

    return (Temp,E_avg,E_edgV,E_edgH,E_avg_MGQD,E_edgV_MGQD,E_edgH_MGQD,E_avg_HO,E_edgV_HO,E_edgH_HO,Fx_edgV,\
    Fy_edgH,Fx_edgV_MGQD,Fy_edgH_MGQD,Fx_edgV_HO,Fy_edgH_HO,Eg_avg,Eg_edgV,Eg_edgH,Eg_avg_HO,Eg_edgV_HO,\
    Eg_edgH_HO,Fxg_edgV,Fyg_edgH,Fxg_edgV_HO,Fyg_edgH_HO,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,\
    fg_edgH_yy,fg_edgH_xy)

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

    ranks=[]
    rrank_BCg = []
    rrank_fg_avg_xx = []
    rrank_fg_edgV_xx = []
    rrank_fg_avg_yy = []
    rrank_fg_edgH_yy = []
    rrank_fg_edgV_xy = []
    rrank_fg_edgH_xy = []

    for i in range(nsets-1):

        rankfile.write('Dataset: '+casename[i+1]+'\n')
        rankfile.write('--------------------------------------------------\n')

        set = dsets[i+1]
        gsum = set.getncattr("PODgsum")
        if (gsum == 1):
            rrank_BCg.append(set['rrank_BCg'][0])
            rrank_fg_avg_xx.append(set['rrank_fg_avg_xx'][0])
            rrank_fg_edgV_xx.append(set['rrank_fg_edgV_xx'][0])
            rrank_fg_avg_yy.append(set['rrank_fg_avg_yy'][0])
            rrank_fg_edgH_yy.append(set['rrank_fg_edgH_yy'][0])
            rrank_fg_edgV_xy.append(set['rrank_fg_edgV_xy'][0])
            rrank_fg_edgH_xy.append(set['rrank_fg_edgH_xy'][0])

            rankfile.write('BCg: {}\n'.format(rrank_BCg[i]))
            rankfile.write('fg_avg_xx: {}\n'.format(rrank_fg_avg_xx[i]))
            rankfile.write('fg_avg_yy: {}\n'.format(rrank_fg_avg_yy[i]))
            rankfile.write('fg_edgV_xx: {}\n'.format(rrank_fg_edgV_xx[i]))
            rankfile.write('fg_edgH_yy: {}\n'.format(rrank_fg_edgH_yy[i]))
            rankfile.write('fg_edgV_xy: {}\n'.format(rrank_fg_edgV_xy[i]))
            rankfile.write('fg_edgH_xy: {}\n'.format(rrank_fg_edgH_xy[i]))

        else:
            rrank_BCg = set['rrank_BCg'][:]
            rrank_fg_avg_xx = set['rrank_fg_avg_xx'][:]
            rrank_fg_edgV_xx = set['rrank_fg_edgV_xx'][:]
            rrank_fg_avg_yy = set['rrank_fg_avg_yy'][:]
            rrank_fg_edgH_yy = set['rrank_fg_edgH_yy'][:]
            rrank_fg_edgV_xy = set['rrank_fg_edgV_xy'][:]
            rrank_fg_edgH_xy = set['rrank_fg_edgH_xy'][:]

            rankfile.write('BCg:')
            for g in range(N_g): rankfile.write('  {}'.format(rrank_BCg[g]))
            rankfile.write('\n')

            rankfile.write('fg_avg_yy:')
            for g in range(N_g): rankfile.write('  {}'.format(rrank_fg_avg_yy[g]))
            rankfile.write('\n')

            rankfile.write('fg_avg_yy:')
            for g in range(N_g): rankfile.write('  {}'.format(rrank_fg_avg_yy[g]))
            rankfile.write('\n')

            rankfile.write('fg_edgV_xx:')
            for g in range(N_g): rankfile.write('  {}'.format(rrank_fg_edgV_xx[g]))
            rankfile.write('\n')

            rankfile.write('fg_edgH_yy:')
            for g in range(N_g): rankfile.write('  {}'.format(rrank_fg_edgH_yy[g]))
            rankfile.write('\n')

            rankfile.write('fg_edgV_xy:')
            for g in range(N_g): rankfile.write('  {}'.format(rrank_fg_edgV_xy[g]))
            rankfile.write('\n')

            rankfile.write('fg_edgH_xy:')
            for g in range(N_g): rankfile.write('  {}'.format(rrank_fg_edgH_xy[g]))
            rankfile.write('\n')

        rankfile.write('\n')
    rankfile.close()
    os.system('mv Ranks.txt '+proc_dir)

    if (gsum == 1):
        # rleg = ['$\mathcal{C}$','$\mathcal{F}_{a,xx}$','$\mathcal{F}_{v,xx}$','$\mathcal{F}_{a,yy}$','$\mathcal{F}_{h,yy}$','$\mathcal{F}_{v,xy}$','$\mathcal{F}_{h,xy}$']
        ranks = [rrank_BCg,rrank_fg_avg_xx,rrank_fg_edgV_xx,rrank_fg_avg_yy,rrank_fg_edgH_yy,rrank_fg_edgV_xy,rrank_fg_edgH_xy]
        # pltr.lineplot("Ranks",trend_names,ranks,proc_dir,rleg,yscale='linear',xlabel=r'POD error $\varepsilon$',ylabel='Dimensionality (k)',marker='D',legloc='upper left')

    return (ranks,gsum)
