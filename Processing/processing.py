#!/usr/bin/env python3
#==================================================================================================================================#
#
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
# import numpy as np
import math_tools as mt

#==================================================================================================================================#
#
#==================================================================================================================================#
def data_compare(data,A):
    yp = len(data[0]); xp = len(data[0][0])

    err=[]; ref=[]
    for j in range(yp):
        for i in range(xp):
            ref.append(data[0][j][i])
            err.append(data[0][j][i] - data[1][j][i])

    enorm_inf = mt.norm_inf(err)
    enorm2 = mt.norm2(err)
    enormL2 = mt.normL2(err,A)

    rnorm_inf = mt.norm_inf(err)/mt.norm_inf(ref)
    rnorm2 = mt.norm2(err)/mt.norm2(ref)
    rnormL2 = mt.normL2(err,A)/mt.normL2(ref,A)

    return (enorm_inf,rnorm_inf,enorm2,rnorm2,enormL2,rnormL2)

#==================================================================================================================================#
#
#==================================================================================================================================#
def get_data(dset,inp_flags,t):
    Temp = dset['Temperature'][t][:]

    E_avg=[]; E_edgV=[]; E_edgH=[]
    if inp_flags[12] == 1:
        E_avg = dset['E_avg_Grey'][t][:]
        E_edgV = dset['E_edgV_Grey'][t][:]
        E_edgH = dset['E_edgH_Grey'][t][:]

    E_avg_MGQD=[]; E_edgV_MGQD=[]; E_edgH_MGQD=[]
    if inp_flags[9] == 1:
        E_avg_MGQD = dset['E_avg_MGQD'][t][:]
        E_edgV_MGQD = dset['E_edgV_MGQD'][t][:]
        E_edgH_MGQD = dset['E_edgH_MGQD'][t][:]

    E_avg_HO=[]; E_edgV_HO=[]; E_edgH_HO=[]
    if inp_flags[5] == 1:
        E_avg_HO = dset['E_avg_HO'][t][:]
        E_edgV_HO = dset['E_edgV_HO'][t][:]
        E_edgH_HO = dset['E_edgH_HO'][t][:]

    Fx_edgV=[]; Fy_edgH=[]
    if inp_flags[13] == 1:
        Fx_edgV = dset['Fx_edgV_Grey'][t][:]
        Fy_edgH = dset['Fy_edgH_Grey'][t][:]

    Fx_edgV_MGQD=[]; Fy_edgH_MGQD=[]
    if inp_flags[11] == 1:
        Fx_edgV_MGQD = dset['Fx_edgV_MGQD'][t][:]
        Fy_edgH_MGQD = dset['Fy_edgH_MGQD'][t][:]

    Fx_edgV_HO=[]; Fy_edgH_HO=[]
    if inp_flags[7] == 1:
        Fx_edgV_HO = dset['Fx_edgV_HO'][t][:]
        Fy_edgH_HO = dset['Fy_edgH_HO'][t][:]

    Eg_avg=[]; Eg_edgV=[]; Eg_edgH=[]
    if inp_flags[8] == 1:
        Eg_avg = dset['Eg_avg_MGQD'][t][:]
        Eg_edgV = dset['Eg_edgV_MGQD'][t][:]
        Eg_edgH = dset['Eg_edgH_MGQD'][t][:]

    Eg_avg_HO=[]; Eg_edgV_HO=[]; Eg_edgH_HO=[]
    if inp_flags[4] == 1:
        Eg_avg_HO = dset['Eg_avg_HO'][t][:]
        Eg_edgV_HO = dset['Eg_edgV_HO'][t][:]
        Eg_edgH_HO = dset['Eg_edgH_HO'][t][:]

    Fxg_edgV=[]; Fyg_edgH=[]
    if inp_flags[10] == 1:
        Fxg_edgV = dset['Fxg_edgV_MGQD'][t][:]
        Fyg_edgH = dset['Fyg_edgH_MGQD'][t][:]

    Fxg_edgV_HO=[]; Fyg_edgH_HO=[]
    if inp_flags[6] == 1:
        Fxg_edgV_HO = dset['Fxg_edgV_HO'][t][:]
        Fyg_edgH_HO = dset['Fyg_edgH_HO'][t][:]

    fg_avg_xx=[]; fg_avg_xy=[]; fg_avg_yy=[]
    fg_edgV_xx=[]; fg_edgV_xy=[]
    fg_edgH_yy=[]; fg_edgH_xy=[]
    if inp_flags[2] == 1:
        fg_avg_xx = dset['fg_avg_xx'][t][:]
        fg_avg_xy = dset['fg_avg_xy'][t][:]
        fg_avg_yy = dset['fg_avg_yy'][t][:]

        fg_edgV_xx = dset['fg_edgV_xx'][t][:]
        fg_edgV_xy = dset['fg_edgV_xy'][t][:]

        fg_edgH_yy = dset['fg_edgH_yy'][t][:]
        fg_edgH_xy = dset['fg_edgH_xy'][t][:]

    return (Temp,E_avg,E_edgV,E_edgH,E_avg_MGQD,E_edgV_MGQD,E_edgH_MGQD,E_avg_HO,E_edgV_HO,E_edgH_HO,Fx_edgV,\
    Fy_edgH,Fx_edgV_MGQD,Fy_edgH_MGQD,Fx_edgV_HO,Fy_edgH_HO,Eg_avg,Eg_edgV,Eg_edgH,Eg_avg_HO,Eg_edgV_HO,\
    Eg_edgH_HO,Fxg_edgV,Fyg_edgH,Fxg_edgV_HO,Fyg_edgH_HO,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,\
    fg_edgH_yy,fg_edgH_xy)

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
