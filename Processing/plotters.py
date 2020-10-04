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

#local tools
import misc_tools as msc

#==================================================================================================================================#
#
#==================================================================================================================================#
def plot_results(dset,plotdir,xp,yp,time,inp_flags,plotdirs,pd_names,Temp,fg_avg_xx,fg_avg_xy,fg_avg_yy,E_avg,E_avg_MGQD,\
E_avg_HO,Eg_avg,Eg_avg_HO):

    loc = msc.locate_pd('Temperature',pd_names); dir = plotdirs[loc]
    plot_data(Temp,'Temperature',xp,yp,time,dir)

    if inp_flags[2] == 1:
        loc = msc.locate_pd('fg_xx',pd_names); dir = plotdirs[loc]
        plot_data(fg_avg_xx,'fg_avg_xx',xp,yp,time,dir)

        loc = msc.locate_pd('fg_yy',pd_names); dir = plotdirs[loc]
        plot_data(fg_avg_xy,'fg_avg_yy',xp,yp,time,dir)

        loc = msc.locate_pd('fg_xy',pd_names); dir = plotdirs[loc]
        plot_data(fg_avg_yy,'fg_avg_xy',xp,yp,time,dir)

    if inp_flags[12] == 1: loc = msc.locate_pd('E_avg_Grey',pd_names); dir = plotdirs[loc]; plot_data(E_avg,'E_avg_Grey',xp,yp,time,dir)

    if inp_flags[9] == 1: loc = msc.locate_pd('E_avg_MGQD',pd_names); dir = plotdirs[loc]; plot_data(E_avg_MGQD,'E_avg_MGQD',xp,yp,time,dir)

    if inp_flags[8] == 1: loc = msc.locate_pd('Eg_avg_MGQD',pd_names); dir = plotdirs[loc]; plot_data(Eg_avg,'Eg_avg_MGQD',xp,yp,time,dir)

    if inp_flags[5] == 1: loc = msc.locate_pd('E_avg_HO',pd_names); dir = plotdirs[loc]; plot_data(E_avg_HO,'E_avg_HO',xp,yp,time,dir)

    if inp_flags[4] == 1: loc = msc.locate_pd('Eg_avg_HO',pd_names); dir = plotdirs[loc]; plot_data(Eg_avg_HO,'Eg_avg_HO',xp,yp,time,dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def plot_data(data,name,xp,yp,time,drop):
    if len(np.shape(data)) == 3:
        gp = len(data[0])
        for g in range(gp): heatmap2d(name+'_g'+str(g+1),time,xp,yp,data[g],drop)
    else:
        heatmap2d(name,time,xp,yp,data,drop)

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
def lineplot(name,tp,arr: np.ndarray,drop,legend):
    if hasattr(arr[0],'__len__'): sets = len(arr)
    else: sets = 1

    fig, ax = plt.subplots() #creating the figure

    if sets == 1:
        ax.plot(tp,arr)
    else:
        for i in range(sets):
            leg = legend[i]
            ax.plot(tp,arr[i],label=leg)
