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

#==================================================================================================================================#
# function plot_heatmap2D / plot_heatmap3D
#
# Usage: plot_heatmap2D(data,name,xp,yp,drop,N_y,N_x)
#        plot_heatmap3D(data,name,xp,yp,drop,N_z,N_y,N_x)
#
# INPUTS:
#   data - array of data to plot
#   name - name of data to be plot
#   xp - array of cell-face x-coordinates
#   yp - array of cell-face y-coordinates
#   drop - path to directory where plot should be dropped at
#   N_z - number of z grid points
#   N_y - number of y grid points
#   N_x - number of x grid points
#
#==================================================================================================================================#
def plot_heatmap2D(data,name,xp,yp,drop,N_y,N_x,c='inferno',units=[],bnds=[],title=''):
    dat = np.reshape(data,(N_y,N_x)) #reforming vector of data to 2D array
    heatmap2d(dat,name,xp,yp,drop,units,bnds,c,title=title) #plotting data over spatial grid

def plot_heatmap3D(data,name,xp,yp,drop,N_z,N_y,N_x,c='inferno',units=[],bnds=[],zlab='',title=''):
    dat = np.reshape(data,(N_z,N_y,N_x)) #reforming vector of data to N_g 2D arrays
    for z in range(N_z):
        if bnds: heatmap2d(dat[z],name+'_'+zlab+str(z+1),xp,yp,drop,units,bnds[g],c,title=title) #plotting each group data over spatial grid
        else: heatmap2d(dat[z],name+'_'+zlab+str(z+1),xp,yp,drop,units,c=c,title=title)

#==================================================================================================================================#
#
#==================================================================================================================================#
def heatmap2d(arr: np.ndarray,name,xp,yp,drop,units=[],bnds=[],c='inferno',title=''):
    fig, ax = plt.subplots() #creating the figure

    if bnds and not bnds == [0.,0.]:
        plt.pcolormesh(xp,yp,arr,cmap=c,vmin=bnds[0],vmax=bnds[1]) #plotting
    else:
        plt.pcolormesh(xp,yp,arr,cmap=c) #plotting

    clb = plt.colorbar() #creating color bar
    if units:
        clb.ax.set_title(units)

    if (title != ''):
        plt.title(title)

    tickfreq = 1 #frequency of ticks on each axis
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the x-axis
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the y-axis
    plt.savefig(drop+'/'+name+'.pdf') #saving the figure
    plt.close(fig) #closing out figure

#==================================================================================================================================#
#
#==================================================================================================================================#
def lineplot_grouped(name,xp,N_g,arr: np.ndarray,drop,legend='',yscale='linear',xlabel='',ylabel='',marker='',legloc='upper right',title=''):
    N_x = len(xp)
    dat = np.reshape(arr,(N_g,N_x))
    for g in range(N_g):
        lineplot(name+'_'+str(g+1),xp,dat[g],drop,legend,yscale,xlabel,ylabel,marker,legloc,title)

#==================================================================================================================================#
#
#==================================================================================================================================#
def lineplot(name,tp,arr: np.ndarray,drop,legend='',yscale='linear',xlabel='',ylabel='',marker='',legloc='upper right',title=''):
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

    if isinstance(legend,(list,np.ndarray)):
        ax.legend()
        plt.legend(loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize='small')
        # plt.legend(bbox_to_anchor=(1.2, 1), loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize='small')

    if (title != ''):
        plt.title(title)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.tick_params(axis='y',direction='in',labelleft=True,left=True,right=True)
    ax.tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True)

    ax.set_yscale(yscale)
    plt.grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--')
    # plt.grid(which = 'minor', axis = 'y', color = 'lightgrey', linestyle = ':', linewidth=.75)

    plt.savefig(drop+'/'+name+'.pdf') #saving the figure
    plt.close(fig) #closing out figure

#==================================================================================================================================#
#
#==================================================================================================================================#
def scatterplot(arr: np.ndarray, name, drop, xlabel='', ylabel='', marker='.', c='black', legend=[], yscale='linear',
    xlim=[], ylim=[], cr=0, ygrid=True, xgrid=True):

    if hasattr(arr[0][0],'__len__'):
        sets = len(arr[0])
    else:
        sets = 0

    plt.rc('font', family='serif', size=10)
    plt.rc('xtick', labelsize='medium')
    plt.rc('ytick', labelsize='medium')

    fig, ax = plt.subplots(figsize=(7,5)) #creating the figure

    if sets == 0:
        ax.scatter(arr[0],arr[1],marker=marker,c=c)
    else:
        for i in range(sets):
            leg = legend[i]
            ax.scatter(arr[i][0],arr[i][1],label=leg,marker=marker[i],c=c[i])

        ax.legend()
        plt.legend(loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize='small')
        # plt.legend(bbox_to_anchor=(1.2, 1), loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize='small')

    if cr != 0:
        ax.set_aspect(1)
        circ = plt.Circle((0,0),radius=cr,fill=False)
        ax.add_artist(circ)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    if xlim: ax.set_xlim(xlim)
    if ylim: ax.set_ylim(ylim)

    ax.tick_params(axis='y',direction='in',labelleft=True,left=True,right=True)
    ax.tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True)

    ax.set_yscale(yscale)

    if ygrid:
        plt.grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--')
        plt.grid(which = 'minor', axis = 'y', color = 'lightgrey', linestyle = ':', linewidth=.75)
    if xgrid:
        plt.grid(which = 'major', axis = 'x', color = 'lightgrey', linestyle = '--')
        plt.grid(which = 'minor', axis = 'x', color = 'lightgrey', linestyle = ':', linewidth=.75)

    ax.set_axisbelow(True)

    plt.savefig(drop+'/'+name+'.pdf') #saving the figure
    plt.close(fig) #closing out figure
