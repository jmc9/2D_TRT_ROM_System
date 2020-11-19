#!/usr/bin/env python3
#==================================================================================================================================#
#
# order_plts.py is a small python script to take plots output from TRT_processor.py and rename them in list format
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
import os
import getopt, sys

#==================================================================================================================================#
#
#==================================================================================================================================#
def order(times,vname,indir,outdir):
    check_dirs(indir,outdir)
    ipath = indir+'/'
    opath = outdir+'/'

    i = 0
    for t in times:
        file = ipath+vname+'_'+'{:.2e}'.format(t)+'.pdf'
        if os.path.isfile(file):
            i = i + 1
            new = opath+vname+'_'+'{:d}'.format(i)+'.pdf'
            os.system('cp '+file+' '+new)

#==================================================================================================================================#
#
#==================================================================================================================================#
def gen_times(start,dt,end):
    tp = round((end-start)/dt)

    times = []
    time = start
    for i in range(tp):
        time = time + dt
        times.append(time)

    return times

#==================================================================================================================================#
#
#==================================================================================================================================#
def check_dirs(indir,outdir):
    if not os.path.isdir(indir):
        print('Could not locate input directory ('+indir+'), terminating program')
        quit()

    dirset(outdir)

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

#==================================================================================================================================#
#
#==================================================================================================================================#
def options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:o:v:t:d:s:', ['help','input=','output=','vname='])
    except getopt.GetoptError as err:
        print(str(err))
        usage(); quit()

    (indir,outdir,vname,time,dt,start) = defaults()

    for o, a in opts:
        if o in ['-h','--help']:
            usage(); quit()
        elif o in ['-i','--input']:
            indir = a
        elif o in ['-o','--output']:
            outdir = a
        elif o in ['-v','--vname']:
            vname = a
        elif o in ['-t','--time']:
            time = float(a)
        elif o in ['-d','--dt']:
            dt = float(a)
        elif o in ['-s','--start']:
            start = float(a)

    return (indir,outdir,vname,time,dt,start)

#==================================================================================================================================#
#
#==================================================================================================================================#
def defaults():
    indir = 'Temperature'
    outdir = 'Temperature_ordered'
    vname = 'Temperature'
    time = 6.
    dt = 2e-2
    start = 0.

    return (indir,outdir,vname,time,dt,start)

#==================================================================================================================================#
#
#==================================================================================================================================#
def usage():
    print()
    print('- order_plts.py')
    print('- Joseph M. Coale, jmcoale@ncsu.edu')
    print()
    print('This python3 code is used to order plots generated by TRT_processor.py in list format (Author: Joseph M. Coale)')
    print('Below are all available options and their usages')
    print()
    print('Options:')
    print('-h [help]: prints this message')
    print('-i [input]: used to specify the directory containing plots to order')
    print('-o [output]: used to specify the output directory')
    print('-v [vname]: specify name of the variable to order (i.e. Temperature)')
    print('-t [time]: end time')
    print('-d [dt]: time step length')
    print('-s [start]: start time')
    print()

#==================================================================================================================================#
# MAIN PROGRAM
# This portion of the code is run when this file is the main program
# i.e. when the call is made './TRT_processor.py'
#==================================================================================================================================#
if __name__ == '__main__':
    (indir,outdir,vname,time,dt,start) = options()
    times = gen_times(start,dt,time)
    order(times,vname,indir,outdir)
