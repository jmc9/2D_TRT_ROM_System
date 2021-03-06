#!/usr/bin/env python3
#==================================================================================================================================#
#
# ToolBox.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import os
import math
import numpy as np
import copy

#==================================================================================================================================#
#
#==================================================================================================================================#
def FileCheck(file):
    if not os.path.isfile(file):
        print('Input file, '+file+', does not exist. aborting program')
        quit()

#==================================================================================================================================#
# function dirset
#
# Usage: dirset(dir)
#
# dirset 'sets up' a directory on the given path 'dir'
# - if the directory does not exist, it is created
# - if the directory does exist, its contents are deleted
#
# INPUTS:
#   dir - path to target directory
#==================================================================================================================================#
def dirset(dir, replace=True):
    if os.path.isdir(dir): #if directory exists, clear all files from inside
        if ((len(os.listdir(dir)) > 0)and(replace)): os.system('rm -r '+dir+'/*')
    else: #if directory does not exist, create it
        os.mkdir(dir)


#==================================================================================================================================#
#
#==================================================================================================================================#
def Flatten2D(arr):
    arr = [x for Y_list in arr for x in Y_list]
    return arr

def Flatten3D(arr):
    arr = [x for Z_list in arr for Y_list in Z_list for x in Y_list]
    return arr

#==================================================================================================================================#
#
#==================================================================================================================================#
def Cell_Coords(Delx,Dely):
    xp = [0.]
    for i in range(len(Delx)): xp.append(sum(Delx[:i]) + Delx[i])
    yp = [0.]
    for i in range(len(Dely)): yp.append(sum(Dely[:i]) + Dely[i])

    return (xp,yp)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Grid_Edg_Gen(Grid):
    pts = len(Grid.crds) - 1
    grid_edg = [Grid.bnds[0]]
    for i in range(pts):
        grid_edg.append( ( Grid.crds[i] + Grid.crds[i+1] ) / 2. )
    grid_edg.append(Grid.bnds[1])

    return grid_edg

#==================================================================================================================================#
#
#==================================================================================================================================#
def norm(arr: np.ndarray, type):
    n = 0

    if type in [2,'two']:
        arrfl = arr.flatten()
        n = 0.
        for x in arrfl:
            n = n + x**2
        n = math.sqrt(n)

    elif type in [0,'inf']:
        n = np.amax(np.absolute(arr))

    return n


#==================================================================================================================================#
#
#==================================================================================================================================#
def pow10(num):
    if num!= 0.:
        l = math.log10(abs(num))
        if l >= 0:
            exp = math.floor(l)
        else:
            exp = math.ceil(l)
        coef = num / 10.**exp

    else:
        exp = 0
        coef = 0.

    return (exp, coef)

#==================================================================================================================================#
#
#==================================================================================================================================#
def stack_data(stack_dat_):
    shp = np.asarray(np.shape(stack_dat_))
    ns = shp[0]
    sdim = shp[len(shp)-1]
    dlen = np.prod(shp[1:-1])
    shift = (ns-1)*sdim

    stack_dat = []
    for k in range(ns):
        stack_dat.append(stack_dat_[k].flatten())

    dat = np.zeros(dlen*ns*sdim)

    p_stk = 0
    p_dat = np.linspace(0, (ns-1)*sdim, ns, dtype=int)
    for i in range(dlen):
        for j in range(sdim):
            for k in range(ns):
                dat[p_dat[k]] = stack_dat[k][p_stk]
                p_dat[k] += 1
            p_stk += 1
        p_dat += shift

    return dat

#==================================================================================================================================#
#
#==================================================================================================================================#
def BigTitle():
    print('  _____                                                          ')
    print(' |  __ \\                                                        ')
    print(' | |  | | ___  ___ ___  _ __ ___  _ __                           ')
    print(' | |  | |/ _ \/ __/ _ \| \'_ ` _ \\| \'_ \\                      ')
    print(' | |__| |  __/ (_| (_) | | | | | | |_) |                         ')
    print(' |_____/ \\___|\\___\\___/|_| |_| |_| .__/                       ')
    print('              |  __ \\            | |                            ')
    print('              | |__) | __ ___   _|_|___  ___ ___  ___  _ __      ')
    print('              |  ___/ \'__/ _ \\ / __/ _ \\/ __/ __|/ _ \\| \'__|')
    print('              | |   | | | (_) | (_|  __/\\__ \\__ \\ (_) | |     ')
    print('              |_|   |_|  \\___/ \\___\\___||___/___/\\___/|_|    ')
