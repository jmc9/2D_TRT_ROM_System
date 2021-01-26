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
def dirset(dir, path=''):
    if path != '':
        dir = os.path.join(path, dir)

    if os.path.isdir(dir): #if directory exists, clear all files from inside
        if len(os.listdir(dir)) > 0: os.system('rm -r '+dir+'/*')
    else: #if directory does not exist, create it
        os.mkdir(dir)

    return dir


#==================================================================================================================================#
#
#==================================================================================================================================#
def File_Find(file, key):
    for line in file:
        if (line.strip()): #disregarding blank lines
            input = line.split() #delimiting input line by space

            if input[0] == key:
                arg = input[1]
                file.seek(0)
                return arg

    print('{} not found'.format(key))
    quit()
