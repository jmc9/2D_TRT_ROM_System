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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
# import matplotlib.colors as mcolors

#local tools
import misc_tools as msc

#==================================================================================================================================#
#
#==================================================================================================================================#
def plot_norms(name,norms,tp,casenames,drop,ylabel,xlabel,marker=''):
    lineplot(name+"_abs_inf_Norm",tp,norms[0],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker)
    lineplot(name+"_abs_2_Norm",tp,norms[2],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker)
    lineplot(name+"_abs_L2_Norm",tp,norms[4],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker)
    lineplot(name+"_rel_inf_Norm",tp,norms[1],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker)
    lineplot(name+"_rel_2_Norm",tp,norms[3],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker)
    lineplot(name+"_rel_L2_Norm",tp,norms[5],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker)

#==================================================================================================================================#
# function plot_results
#
# Usage: plot_results(dset,plotdir,xp,yp,time,N_g,N_y,N_x,inp_flags,plotdirs,pd_names,Temp,fg_avg_xx,fg_avg_xy,fg_avg_yy,E_avg,E_avg_MGQD,\
#        E_avg_HO,Eg_avg,Eg_avg_HO)
#
# plot_results generates plots of the data contained in a given dataset
#
#==================================================================================================================================#
def plot_results(dset,plotdir,xp,yp,time,N_g,N_y,N_x,inp_flags,plotdirs,pd_names,Temp,fg_avg_xx,fg_avg_xy,fg_avg_yy,E_avg,E_avg_MGQD,\
E_avg_HO,Eg_avg,Eg_avg_HO,bnds):

    loc = msc.locate_pd('Temperature',pd_names); dir = plotdirs[loc]
    plot_data(Temp,'Temperature',xp,yp,time,dir,N_y,N_x,'inferno','eV',bnds[0])

    if inp_flags[2] == 1:

        loc = msc.locate_pd('fg_xx',pd_names); dir = plotdirs[loc]
        plot_gdata(fg_avg_xx,'fg_avg_xx',xp,yp,time,dir,N_g,N_y,N_x,'inferno',bnds=bnds[4])

        loc = msc.locate_pd('fg_yy',pd_names); dir = plotdirs[loc]
        plot_gdata(fg_avg_yy,'fg_avg_yy',xp,yp,time,dir,N_g,N_y,N_x,'inferno')

        loc = msc.locate_pd('fg_xy',pd_names); dir = plotdirs[loc]
        plot_gdata(fg_avg_xy,'fg_avg_xy',xp,yp,time,dir,N_g,N_y,N_x,'inferno')

    if inp_flags[12] == 1: loc = msc.locate_pd('E_avg_Grey',pd_names); dir = plotdirs[loc]; plot_data(E_avg,'E_avg_Grey',xp,yp,time,dir,N_y,N_x,'inferno','erg$\cdot 10^{13}$',bnds[1])

    if inp_flags[9] == 1: loc = msc.locate_pd('E_avg_MGQD',pd_names); dir = plotdirs[loc]; plot_data(E_avg_MGQD,'E_avg_MGQD',xp,yp,time,dir,N_y,N_x,'inferno','erg$\cdot 10^{13}$',bnds[2])

    if inp_flags[8] == 1: loc = msc.locate_pd('Eg_avg_MGQD',pd_names); dir = plotdirs[loc]; plot_gdata(Eg_avg,'Eg_avg_MGQD',xp,yp,time,dir,N_g,N_y,N_x,'inferno')

    if inp_flags[5] == 1: loc = msc.locate_pd('E_avg_HO',pd_names); dir = plotdirs[loc]; plot_data(E_avg_HO,'E_avg_HO',xp,yp,time,dir,N_y,N_x,'inferno','erg$\cdot 10^{13}$',bnds[3])

    if inp_flags[4] == 1: loc = msc.locate_pd('Eg_avg_HO',pd_names); dir = plotdirs[loc]; plot_gdata(Eg_avg_HO,'Eg_avg_HO',xp,yp,time,dir,N_g,N_y,N_x,'inferno')

#==================================================================================================================================#
# function plot_data / plot_gdata
#
# Usage: plot_data(data,name,xp,yp,time,drop)
#        plot_gdata(data,name,xp,yp,time,drop)
#
# INPUTS:
#   data - array of data to plot
#   name - name of data to be plot
#   xp - array of cell-face x-coordinates
#   yp - array of cell-face y-coordinates
#   time - time instant for data
#   drop - path to directory where plot should be dropped at
#   N_g - number of frequency groups
#   N_y - number of y grid points
#   N_x - number of x grid points
#
#==================================================================================================================================#
def plot_data(data,name,xp,yp,time,drop,N_y,N_x,c,units=[],bnds=[]):
    dat = np.reshape(data,(N_y,N_x)) #reforming vector of data to 2D array
    heatmap2d(dat,name,xp,yp,time,drop,units,bnds,c) #plotting data over spatial grid

def plot_gdata(data,name,xp,yp,time,drop,N_g,N_y,N_x,c,units=[],bnds=[]):
    gdat = np.reshape(data,(N_g,N_y,N_x)) #reforming vector of data to N_g 2D arrays
    for g in range(N_g):
        if bnds: heatmap2d(gdat[g],name+'_g'+str(g+1),xp,yp,time,drop,units,bnds[g],c) #plotting each group data over spatial grid
        else: heatmap2d(gdat[g],name+'_g'+str(g+1),xp,yp,time,drop,units,c=c)

#==================================================================================================================================#
#
#==================================================================================================================================#
def heatmap2d(arr: np.ndarray,name,xp,yp,time,drop,units=[],bnds=[],c='inferno'):
    fig, ax = plt.subplots() #creating the figure

    if bnds and not bnds == [0.,0.]:
        plt.pcolormesh(xp,yp,arr,cmap=c,vmin=bnds[0],vmax=bnds[1]) #plotting
    else:
        plt.pcolormesh(xp,yp,arr,cmap=c) #plotting

    clb = plt.colorbar() #creating color bar
    if units:
        clb.ax.set_title(units)

    tickfreq = 1 #frequency of ticks on each axis
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the x-axis
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the y-axis
    plt.savefig(drop+'/'+name+'_'+'{:.2e}'.format(time)+'.pdf') #saving the figure
    plt.close(fig) #closing out figure

#==================================================================================================================================#
#
#==================================================================================================================================#
def lineplot(name,tp,arr: np.ndarray,drop,legend,yscale='linear',xlabel='',ylabel='',marker='',legloc='upper right'):
    if hasattr(arr[0],'__len__'):
        sets = len(arr)

        if isinstance(marker,str):
            mrk = marker
            marker = []
            for i in range(sets): marker.append(mrk)

        elif isinstance(marker,list):
            if len(marker)!=sets:
                print('Error in lineplot - marker is a list but does not have as many elements as there are datasets')
                quit()

        else:
            print('Error in lineplot - marker not a list or a single string while multiple datasets exist')
            quit()

    else:
        sets = 0

        if not isinstance(marker,str):
            print('Error in lineplot - marker not a string while only one dataset exists')
            quit()

    plt.rc('font', family='serif', size=10)
    plt.rc('xtick', labelsize='medium')
    plt.rc('ytick', labelsize='medium')

    fig, ax = plt.subplots(figsize=(7,5)) #creating the figure
    # plt.subplots_adjust(left=0.15,bottom=0.15,right=.85,top=.95)

    if sets == 0:
        ax.plot(tp,arr,marker=marker)
    else:
        for i in range(sets):
            leg = legend[i]
            ax.plot(tp,arr[i],label=leg,linewidth=.75,marker=marker[i])

    ax.legend()
    plt.legend(loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize='small')
    # plt.legend(bbox_to_anchor=(1.2, 1), loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize='small')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.tick_params(axis='y',direction='in',labelleft=True,left=True,right=True)
    ax.tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True)

    ax.set_yscale(yscale)
    plt.grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--')
    plt.grid(which = 'minor', axis = 'y', color = 'lightgrey', linestyle = ':', linewidth=.75)

    plt.savefig(drop+'/'+name+'.pdf') #saving the figure
    plt.close(fig) #closing out figure
