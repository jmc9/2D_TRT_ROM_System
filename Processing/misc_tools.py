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
# function plotdirs
#
# Usage: (plotdirs,pd_names) = plotdirs(rootdir)
#
# plotdirs generates the directory paths for placing plots of a TRT solution appended to some root directory path, along with -
# - the names of the data that will be held in each directory for search/sort purposes
#
# INPUTS:
#   rootdir (string) - the root directory path in which to place the plot directories
#
# OUTPUTS:
#   plotdirs - array of directory paths
#   pd_names - array containing the 'names' corresponding to the data to be held in each plot directory
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
# function locate_pd
#
# Usage: loc = locate_pd(name,pd_names)
#
# locate_pd gives the position of array pd_names that corresponds to the given string 'name'
#
# INPUTS:
#   pd_names - array containing the 'names' corresponding to the data to be held in each plot directory
#   name - the data name whose location in pd_names is to be found
#
# OUTPUTS:
#   loc - position of 'name' in array pd_names
#
# ERRORS:
#   if loc returns a value of -1, then 'name' was not contained in pd_names
#
#==================================================================================================================================#
def locate_pd(name,pd_names):
    loc = -1
    n = len(pd_names)
    for i in range(n):
        if pd_names[i] == name: loc = i; break

    return loc

#==================================================================================================================================#
# function setup_plotdirs
#
# Usage: setup_plotdirs(plotdirs,pd_names,inp_flags)
#
# setup_plotdirs creates directories for the plots of a specific dataset's contents
#
# INPUTS:
#   plotdirs - array of directory paths
#   pd_names - array containing the 'names' corresponding to the data to be held in each plot directory
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
def dirset(dir):
    if os.path.isdir(dir): #if directory exists, clear all files from inside
        if len(os.listdir(dir)) > 0: os.system('rm -r '+dir+'/*')
    else: #if directory does not exist, create it
        os.mkdir(dir)
