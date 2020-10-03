#!/usr/bin/env python3
#==================================================================================================================================#
#
# TRT_processor.py is the main python script to process the outputs of the 2D TRT ROM code
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

#netcdf tools
from netCDF4 import Dataset as ncds

#==================================================================================================================================#
#
#==================================================================================================================================#
def TRT_process():
    inpdir = 'output.h5'
    dset = ncds(inpdir,'r') #opening dataset of target TRT results

    plotdir = 'plots'
    dirset(plotdir)

    Delx = dset['Delx'][:]; Dely = dset['Dely'][:]
    Delt = dset['Delt']; N_t = dset.dimensions['N_t'].size

    Delt = 2e-3

    (xp,yp) = cell_coords(Delx,Dely)
    tp = []
    for i in range(N_t): tp.append((i+1)*Delt)

    plot_data(dset,'E_avg_Grey',xp,yp,tp,plotdir)

    dset.close() #closing dataset file

#==================================================================================================================================#
#
#==================================================================================================================================#
def cell_coords(Delx,Dely):
    xp = [0.]
    for i in range(len(Delx)): xp.append(sum(Delx[:i]) + Delx[i])
    yp = [0.]
    for i in range(len(Dely)): yp.append(sum(Dely[:i]) + Dely[i])

    return (xp,yp)

#==================================================================================================================================#
#
#==================================================================================================================================#
def plot_data(dset,name,xp,yp,tp,drop):
    N_t = len(tp)
    data = dset[name][:]
    for i in range(N_t): heatmap2d(name,tp[i],xp,yp,data[i],drop)

#==================================================================================================================================#
#
#==================================================================================================================================#
def dirset(dir):
    if os.path.isdir(dir): #if directory exists, clear all files from inside
        os.system('rm -r '+dir+'/*')
    else: #if directory does not exist, create it
        os.mkdir(dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def heatmap2d(name,time,xp,yp,arr: np.ndarray,drop):
    fig, ax = plt.subplots() #creating the figure
    plt.pcolormesh(xp,yp,arr,cmap='inferno') #plotting
    plt.colorbar() #creating color bar
    tickfreq = 1 #frequency of ticks on each axis
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the x-axis
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the y-axis
    plt.savefig(drop+'/'+name+'_'+'{:.2e}'.format(time)+'.pdf') #saving the figure
    plt.close(fig) #closing out figure

#==================================================================================================================================#
#
#==================================================================================================================================#
if __name__ == '__main__':
    TRT_process()
