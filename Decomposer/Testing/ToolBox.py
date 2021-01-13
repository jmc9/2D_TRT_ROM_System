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
def dirset(dir, path):
    if os.path.isdir(dir): #if directory exists, clear all files from inside
        if len(os.listdir(dir)) > 0: os.system('rm -r '+dir+'/*')
    else: #if directory does not exist, create it
        os.mkdir(dir)

    return os.path.join(path, dir)
