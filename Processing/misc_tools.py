#!/usr/bin/env python3
#==================================================================================================================================#
#
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import os

#==================================================================================================================================#
#
#==================================================================================================================================#
def plotdirs(rootdir):
    plotdirs = [rootdir+'/I_avg',\
    rootdir+'/fg_xx',\
    rootdir+'/fg_xy',\
    rootdir+'/fg_yy',\
    rootdir+'/D',\
    rootdir+'/E_avg_Grey',\
    rootdir+'/E_avg_MGQD',\
    rootdir+'/E_avg_HO',\
    rootdir+'/Eg_avg_MGQD',\
    rootdir+'/Eg_avg_HO',\
    rootdir+'/Temperature']

    pd_names = ['I_avg',\
    'fg_xx',\
    'fg_xy',\
    'fg_yy',\
    'D',\
    'E_avg_Grey',\
    'E_avg_MGQD',\
    'E_avg_HO',\
    'Eg_avg_MGQD',\
    'Eg_avg_HO',\
    'Temperature']

    return (plotdirs,pd_names)

#==================================================================================================================================#
#
#==================================================================================================================================#
def locate_pd(name,pd_names):
    loc = -1
    n = len(pd_names)
    for i in range(n):
        if pd_names[i] == name: loc = i; break

    return loc

#==================================================================================================================================#
#
#==================================================================================================================================#
def setup_plotdirs(plotdirs,pd_names,inp_flags):

    loc = locate_pd('Temperature',pd_names); dir = plotdirs[loc]; dirset(dir)

    if inp_flags[2] == 1:
        loc = locate_pd('fg_xx',pd_names); dir = plotdirs[loc]; dirset(dir)
        loc = locate_pd('fg_xy',pd_names); dir = plotdirs[loc]; dirset(dir)
        loc = locate_pd('fg_yy',pd_names); dir = plotdirs[loc]; dirset(dir)

    if inp_flags[12] == 1: loc = locate_pd('E_avg_Grey',pd_names); dir = plotdirs[loc]; dirset(dir)
    if inp_flags[9] == 1: loc = locate_pd('E_avg_MGQD',pd_names); dir = plotdirs[loc]; dirset(dir)
    if inp_flags[5] == 1: loc = locate_pd('E_avg_HO',pd_names); dir = plotdirs[loc]; dirset(dir)
    if inp_flags[8] == 1: loc = locate_pd('Eg_avg_MGQD',pd_names); dir = plotdirs[loc]; dirset(dir)
    if inp_flags[4] == 1: loc = locate_pd('Eg_avg_HO',pd_names); dir = plotdirs[loc]; dirset(dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def dirset(dir):
    if os.path.isdir(dir): #if directory exists, clear all files from inside
        if len(os.listdir(dir)) > 0: os.system('rm -r '+dir+'/*')
    else: #if directory does not exist, create it
        os.mkdir(dir)
