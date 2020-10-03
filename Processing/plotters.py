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
def plot_results(dset,plotdir,xp,yp,tp,inp_flags):

    plot_data(dset,'Temperature',xp,yp,tp,plotdir)

    if inp_flags[2] == 1:
        dir = plotdir+'/fgxx'; msc.dirset(dir)
        plot_data(dset,'fg_avg_xx',xp,yp,tp,dir)

        dir = plotdir+'/fgyy'; msc.dirset(dir)
        plot_data(dset,'fg_avg_yy',xp,yp,tp,dir)

        dir = plotdir+'/fgxy'; msc.dirset(dir)
        plot_data(dset,'fg_avg_xy',xp,yp,tp,dir)

    if inp_flags[12] == 1: plot_data(dset,'E_avg_Grey',xp,yp,tp,plotdir)

    if inp_flags[9] == 1: plot_data(dset,'E_avg_MGQD',xp,yp,tp,plotdir)

    if inp_flags[8] == 1: dir = plotdir+'/Eg_MGQD'; msc.dirset(dir); plot_data(dset,'Eg_avg_MGQD',xp,yp,tp,dir)

    if inp_flags[5] == 1: plot_data(dset,'E_avg_HO',xp,yp,tp,plotdir)

    if inp_flags[4] == 1: dir = plotdir+'/Eg_HO'; msc.dirset(dir); plot_data(dset,'Eg_avg_HO',xp,yp,tp,dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def plot_data(dset,name,xp,yp,tp,drop):
    N_t = len(tp)
    data = dset[name][:]
    if len(np.shape(data)) == 4:
        gp = len(data[0])
        for i in range(N_t):
            for g in range(gp): heatmap2d(name+'_g'+str(g+1),tp[i],xp,yp,data[i][g],drop)
    else:
        for i in range(N_t): heatmap2d(name,tp[i],xp,yp,data[i],drop)

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
