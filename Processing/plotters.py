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
import math
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.gridspec import GridSpec
# import matplotlib.colors as mcolors

#local tools
import misc_tools as msc
from processing import cross_section as cs

#==================================================================================================================================#
#
#==================================================================================================================================#
def Subplot_TE_norms(Tnorms, Ttrends, Enorms, Etrends, tp, trend_times, dsnames, trend_names, drop, TE_bnd=[], fsize=10, yscale='log'):
    n = len(dsnames)
    n2 = len(trend_times)
    Elabel='$\\frac{{|| {v}_{{FOM}} - {v}_{{ROM}} ||_2}}{{|| {v}_{{FOM}} ||_2}}$'.format(v='E')
    Tlabel='$\\frac{{|| {v}_{{FOM}} - {v}_{{ROM}} ||_2}}{{|| {v}_{{FOM}} ||_2}}$'.format(v='T')

    # if ylsize==0: ylsize = fsize
    # if xlsize==0: xlsize = fsize
    # if legsize==0: legsize = 8 #fsize

    plt.rc('font', family='serif', size=fsize)
    plt.rc('mathtext', fontset='stix')
    plt.rc('xtick', labelsize='medium')
    plt.rc('ytick', labelsize='medium')

    fig, axs = plt.subplots(2, 2, sharey='row', sharex='col', figsize=(14,10)) #creating the figure

    # axs[0,0].plot(tp, Enorms[i], label=leg, linewidth=.75, marker=marker[i])
    for i in range(n):
        axs[0,0].plot(tp, Enorms[3][i], label=dsnames[i], linewidth=.75)
        axs[1,0].plot(tp, Tnorms[3][i], label=dsnames[i], linewidth=.75)

    for i in range(n2):
        axs[0,1].plot(trend_names, Etrends[3][i], label=trend_times[i], linewidth=.75, marker='D')
        axs[1,1].plot(trend_names, Ttrends[3][i], label=trend_times[i], linewidth=.75, marker='D')

    axs[0,0].tick_params(axis='y',direction='in',labelleft=True,left=True,right=True, which='both')
    axs[1,0].tick_params(axis='y',direction='in',labelleft=True,left=True,right=True, which='both')
    axs[0,1].tick_params(axis='y',direction='in',labelleft=True,left=True,right=True, which='both')
    axs[1,1].tick_params(axis='y',direction='in',labelleft=True,left=True,right=True, which='both')

    axs[0,0].tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True, which='both')
    axs[1,0].tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True, which='both')
    axs[0,1].tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True, which='both')
    axs[1,1].tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True, which='both')

    plt.setp(axs[0,0].get_xticklabels(), visible=False)
    plt.setp(axs[0,1].get_xticklabels(), visible=False)
    plt.setp(axs[0,1].get_yticklabels(), visible=False)
    plt.setp(axs[1,1].get_yticklabels(), visible=False)

    axs[0,0].set_yscale(yscale)
    axs[1,0].set_yscale(yscale)
    axs[0,1].set_yscale(yscale)
    axs[1,1].set_yscale(yscale)

    axs[0,0].set_ylabel(Elabel, fontsize=fsize+4, labelpad=0)
    axs[1,0].set_ylabel(Tlabel, fontsize=fsize+4, labelpad=0)
    axs[1,0].set_xlabel('Time (ns)', fontsize=fsize)
    axs[1,1].set_xlabel('Truncation Criteria $\\xi_{{rel}}$', fontsize=fsize+4)

    axs[0,0].grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--', linewidth=.75)
    axs[1,0].grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--', linewidth=.75)
    axs[0,1].grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--', linewidth=.75)
    axs[1,1].grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--', linewidth=.75)

    if TE_bnd:
        axs[0,0].set_ylim(TE_bnd[3])
        axs[1,0].set_ylim(TE_bnd[3])
        axs[0,1].set_ylim(TE_bnd[3])
        axs[1,1].set_ylim(TE_bnd[3])

    axs[0,0].legend(loc='upper right', borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize=fsize)
    axs[0,1].legend(loc='upper right', borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize=fsize)
    axs[1,0].legend(loc='upper right', borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize=fsize)
    axs[1,1].legend(loc='upper right', borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize=fsize)
    # plt.legend(loc='upper right', borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize=fsize)


    fig.tight_layout(pad=1.)
    axs[0,0].minorticks_off()
    axs[0,1].minorticks_off()
    axs[1,0].minorticks_off()
    axs[1,1].minorticks_off()

    plt.savefig(drop+'/TE_combined_norms.pdf') #saving the figure
    plt.close(fig) #closing out figure

#==================================================================================================================================#
#
#==================================================================================================================================#
def plot_cs(name,cs,xp,times,dsnames,fignames,drop,xlabel,ylabel,textp=[]):
    nd, nt, nx = np.shape(cs)

    # one dataset, all times
    for d in range(nd):
        lineplot(name+"_"+fignames[d], xp, cs[d], drop, times, ylabel=ylabel, xlabel=xlabel)

    # one time, all datasets
    for t in range(nt):
        dat=[]
        for d in range(nd):
            dat.append(cs[d][t])
        lineplot('{}_t={}sh'.format(name,round(times[t],8)), xp, dat, drop, dsnames, ylabel=ylabel, xlabel=xlabel)

    # all times, all datasets
    plt.rc('font', family='serif', size=10)
    plt.rc('xtick', labelsize='medium')
    plt.rc('ytick', labelsize='medium')

    fig, ax = plt.subplots(figsize=(7,5)) #creating the figure

    marker = ['s','x','+','o','v','^','*','D','.']
    # style  = ['-','--','-.',':','densely dashdotted','-','-','-','-']
    style  = [(0,()), (0,(5,5)), (0,(1,1)), (0,(3,5,1,5)), (0,(5,2,5,5,1,4)), (0, (5, 1)), (0, (3, 1, 1, 1))]
    color  = ['black','tab:orange','tab:blue','tab:green','tab:purple','tab:brown','tab:pink','tab:green','tab:cyan']
    for t in range(nt):
        for d in range(nd):
            if t==0:
                ax.plot(xp,cs[d][t],label=dsnames[d],linewidth=.75, markersize=7, marker=marker[d], color=color[d], ls=style[d], markerfacecolor='none')
            else:
                ax.plot(xp,cs[d][t],linewidth=.75, markersize=7, marker=marker[d], color=color[d], ls=style[d], markerfacecolor='none')
        if len(textp)>0:
            plt.annotate('t={}ns'.format(round(times[t],8)),(xp[textp[t]],cs[0][t][textp[t]]),textcoords='offset points',xytext=(5,10))
    ax.legend()
    plt.legend(loc='upper right', borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=2., fontsize='small')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax.tick_params(axis='y',direction='in',labelleft=True,left=True,right=True)
    ax.tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True)

    plt.grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--')
    plt.grid(which = 'minor', axis = 'y', color = 'lightgrey', linestyle = ':', linewidth=.75)

    plt.savefig(drop+'/'+'{}_ALLPLOTS'.format(name)+'.pdf') #saving the figure
    plt.close(fig) #closing out figure

#==================================================================================================================================#
#
#==================================================================================================================================#
def plot_norms(name,norms,tp,casenames,drop,xlabel,ylabel='',marker='',vname='',bnds=[],fsize=10,leg_col=1,xspace=[0.,0.]):
    normlabl = True if ylabel=='' else False

    if normlabl: ylabel='$|| {v}_{{FOM}} - {v}_{{ROM}} ||_\infty$'.format(v=vname)
    bnd = bnds[0] if len(bnds)>0 else []
    lineplot(name+"_abs_inf_Norm",tp,norms[0],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker,ylim=bnd,
    grids='major',minorticks=False,fsize=fsize,leg_col=leg_col,xspace=xspace)

    if normlabl: ylabel='$|| {v}_{{FOM}} - {v}_{{ROM}} ||_2$'.format(v=vname)
    bnd = bnds[2] if len(bnds)>0 else []
    lineplot(name+"_abs_2_Norm",tp,norms[2],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker,ylim=bnd,
    grids='major',minorticks=False,fsize=fsize,leg_col=leg_col,xspace=xspace)

    if normlabl: ylabel='$|| {v}_{{FOM}} - {v}_{{ROM}} ||_{{L_2}}$'.format(v=vname)
    bnd = bnds[4] if len(bnds)>0 else []
    lineplot(name+"_abs_L2_Norm",tp,norms[4],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker,ylim=bnd,
    grids='major',minorticks=False,fsize=fsize,leg_col=leg_col,xspace=xspace)

    if normlabl: ylabel='$\\frac{{|| {v}_{{FOM}} - {v}_{{ROM}} ||_\infty}}{{|| {v}_{{FOM}} ||_\infty}}$'.format(v=vname)
    bnd = bnds[1] if len(bnds)>0 else []
    lineplot(name+"_rel_inf_Norm",tp,norms[1],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker,ylsize=fsize+4,ylpad=0,ylim=bnd,
    grids='major',minorticks=False,fsize=fsize,leg_col=leg_col,xspace=xspace)

    if normlabl: ylabel='$\\frac{{|| {v}_{{FOM}} - {v}_{{ROM}} ||_2}}{{|| {v}_{{FOM}} ||_2}}$'.format(v=vname)
    bnd = bnds[3] if len(bnds)>0 else []
    lineplot(name+"_rel_2_Norm",tp,norms[3],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker,ylsize=fsize+4,ylpad=0,ylim=bnd,
    grids='major',minorticks=False,fsize=fsize,leg_col=leg_col,xspace=xspace)

    if normlabl: ylabel='$\\frac{{|| {v}_{{FOM}} - {v}_{{ROM}} ||_{{L_2}}}}{{|| {v}_{{FOM}} ||_{{L_2}}}}$'.format(v=vname)
    bnd = bnds[5] if len(bnds)>0 else []
    lineplot(name+"_rel_L2_Norm",tp,norms[5],drop,casenames,yscale='log',ylabel=ylabel,xlabel=xlabel,marker=marker,ylsize=fsize+4,ylpad=0,ylim=bnd,
    grids='major',minorticks=False,fsize=fsize,leg_col=leg_col,xspace=xspace)

#==================================================================================================================================#
# function plot_results
#
# Usage: plot_results(dset,plotdir,xp,yp,time,N_g,N_y,N_x,inp_flags,plotdirs,pd_names,Temp,fg_avg_xx,fg_avg_xy,fg_avg_yy,E_avg,E_avg_MGQD,\
#        E_avg_HO,Eg_avg,Eg_avg_HO)
#
# plot_results generates plots of the data contained in a given dataset
#
#==================================================================================================================================#
def plot_results(dset,plotdir,xp,yp,time,N_g,N_y,N_x,inp_flags,plotdirs,pd_names,Temp,fg_avg_xx,fg_avg_xy,fg_avg_yy,E_avg,
E_avg_MGQD,E_avg_HO,Eg_avg,Eg_avg_HO,bnds,plot_qdf,plot_E_HO):

    loc = msc.locate_pd('Temperature',pd_names); dir = plotdirs[loc]
    Temp_ = [T/1000. for T in Temp]
    plot_data(Temp_,'Temperature',xp,yp,time,dir,N_y,N_x,'inferno','KeV',bnds[0],xlabel='x (cm)',ylabel='y (cm)')

    if ((inp_flags[2] == 1)and(plot_qdf)):

        loc = msc.locate_pd('fg_xx',pd_names); dir = plotdirs[loc]
        plot_gdata(fg_avg_xx,'fg_avg_xx',xp,yp,time,dir,N_g,N_y,N_x,'inferno',bnds=bnds[4])

        loc = msc.locate_pd('fg_yy',pd_names); dir = plotdirs[loc]
        plot_gdata(fg_avg_yy,'fg_avg_yy',xp,yp,time,dir,N_g,N_y,N_x,'inferno')

        loc = msc.locate_pd('fg_xy',pd_names); dir = plotdirs[loc]
        plot_gdata(fg_avg_xy,'fg_avg_xy',xp,yp,time,dir,N_g,N_y,N_x,'inferno')

    if inp_flags[12] == 1:
        loc = msc.locate_pd('E_avg_Grey',pd_names); dir = plotdirs[loc]
        plot_data(E_avg,'E_avg_Grey',xp,yp,time,dir,N_y,N_x,'inferno','$10^{13}$erg',bnds[1],xlabel='x (cm)',ylabel='y (cm)')

    if plot_E_HO:
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
def plot_data(data,name,xp,yp,time,drop,N_y,N_x,c,units=[],bnds=[],xlabel='',ylabel=''):
    dat = np.reshape(data,(N_y,N_x)) #reforming vector of data to 2D array
    heatmap2d(dat,name,xp,yp,time,drop,units,bnds,c,xlabel=xlabel,ylabel=ylabel) #plotting data over spatial grid

def plot_gdata(data,name,xp,yp,time,drop,N_g,N_y,N_x,c,units=[],bnds=[],xlabel='',ylabel=''):
    gdat = np.reshape(data,(N_g,N_y,N_x)) #reforming vector of data to N_g 2D arrays
    for g in range(N_g):
        if bnds: #plotting each group data over spatial grid
            heatmap2d(gdat[g],name+'_g'+str(g+1),xp,yp,time,drop,units,bnds[g],c,xlabel=xlabel,ylabel=ylabel)
        else:
            heatmap2d(gdat[g],name+'_g'+str(g+1),xp,yp,time,drop,units,c=c,xlabel=xlabel,ylabel=ylabel)

#==================================================================================================================================#
#
#==================================================================================================================================#
def heatmap2d(arr: np.ndarray,name,xp,yp,time,drop,units=[],bnds=[],c='inferno',save_format=0,clbsize=20.,fsize=15.,xlabel='',
ylabel='',clb_tpad=12,clb_tloc='left'):
    plt.rc('font', family='serif', size=fsize)
    fig, ax = plt.subplots(figsize=(7,5)) #creating the figure

    if (len(bnds)>0) and (not np.all(np.isclose(0.,bnds))):
        plt.pcolormesh(xp,yp,arr,cmap=c,vmin=bnds[0],vmax=bnds[1]) #plotting
    else:
        plt.pcolormesh(xp,yp,arr,cmap=c) #plotting

    clb = plt.colorbar(format=ticker.FuncFormatter(fmt)) #creating color bar
    if units:
        clb.ax.set_title(units, pad=clb_tpad, loc=clb_tloc)
    clb.ax.tick_params(labelsize=clbsize)

    ax.set_xlabel(xlabel, fontsize=fsize)
    ax.set_ylabel(ylabel, fontsize=fsize)

    tickfreq = 1 #frequency of ticks on each axis
    ax.xaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the x-axis
    ax.yaxis.set_major_locator(ticker.MultipleLocator(tickfreq)) #setting frequency of ticks along the y-axis
    if save_format==1: plt.savefig(drop+'/'+name+'_'+'{:.2e}'.format(time)+'.pdf') #saving the figure
    else: plt.savefig(drop+'/'+name+'_'+'{:.0f}ns'.format(time)+'.pdf') #saving the figure
    plt.close(fig) #closing out figure

def fmt(x, pos):
    a, b = '{:.1e}'.format(x).split('e')
    b = int(b)
    return r'${} \times 10^{{{}}}$'.format(a, b)

#==================================================================================================================================#
#
#==================================================================================================================================#
def lineplot(name,tp,arr: np.ndarray,drop,legend,yscale='linear',xlabel='',ylabel='',marker='',legloc='upper right',fsize=10,
ylsize=0,xlsize=0,ylpad=4,xlpad=4,legsize=0,ylim=[],grids='both',minorticks=True,leg_col=1,xspace=[0.,0.],xtsize=0,ytsize=0,
leganchor=(-1,-1),ylabside='left',linewidth=.75,linestyle='-',color=''):
    if hasattr(arr[0],'__len__'):
        sets = len(arr)

        ###
        ###
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

        ###
        ###
        if isinstance(linestyle,str):
            ls = linestyle
            linestyle = []
            for i in range(sets): linestyle.append(ls)

        elif isinstance(linestyle,list):
            if len(linestyle)!=sets:
                print('Error in lineplot - linestyle is a list but does not have as many elements as there are datasets')
                quit()

        else:
            print('Error in lineplot - linestyle not a list or a single string while multiple datasets exist')
            quit()

    else:
        sets = 0

        if not isinstance(marker,str):
            print('Error in lineplot - marker not a string while only one dataset exists')
            quit()

    if ylsize==0: ylsize = fsize
    if xlsize==0: xlsize = fsize
    if ytsize==0: ytsize = fsize
    if xtsize==0: xtsize = fsize
    if legsize==0: legsize = fsize-2

    plt.rc('font', family='serif', size=fsize)
    plt.rc('mathtext', fontset='stix')
    plt.rc('xtick', labelsize=xtsize)
    plt.rc('ytick', labelsize=ytsize)

    fig, ax = plt.subplots(figsize=(7,5)) #creating the figure
    # plt.subplots_adjust(left=0.15,bottom=0.15,right=.85,top=.95)

    if sets == 0:
        if color!='':
            ax.plot(tp,arr,linewidth=linewidth,marker=marker,linestyle=linestyle,c=color)
        else:
            ax.plot(tp,arr,linewidth=linewidth,marker=marker,linestyle=linestyle)
    else:
        for i in range(sets):
            leg = legend[i]
            ax.plot(tp,arr[i],label=leg,linewidth=linewidth,marker=marker[i],linestyle=linestyle[i])

    if legend:
        ax.legend()
        if leganchor==(-1,-1):
            plt.legend(loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1.,
            fontsize=legsize, ncol=leg_col)
        else:
            plt.legend(loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1.,
            fontsize=legsize, ncol=leg_col, bbox_to_anchor=leganchor)
        # plt.legend(bbox_to_anchor=(1.2, 1), loc=legloc, borderaxespad=0., edgecolor='black', framealpha=1., fancybox=False, handlelength=1., fontsize='small')

    ax.set_xlabel(xlabel, fontsize=xlsize, labelpad=xlpad)
    ax.set_ylabel(ylabel, fontsize=ylsize, labelpad=ylpad)

    if ylabside=='left':
        lableft=True; labright=False
    else:
        lableft=False; labright=True
    ax.tick_params(axis='y',direction='in',labelleft=lableft,labelright=labright,left=True,right=True, which='both')
    ax.tick_params(axis='x',direction='in',labelbottom=True,bottom=True,top=True, which='both')
    ax.yaxis.set_label_position(ylabside)

    ax.set_yscale(yscale)
    if grids in ['both','major']: plt.grid(which = 'major', axis = 'y', color = 'lightgrey', linestyle = '--')
    if grids in ['both','minor']: plt.grid(which = 'minor', axis = 'y', color = 'lightgrey', linestyle = ':', linewidth=.75)
    if not minorticks: plt.minorticks_off()

    if ylim:
        plt.ylim(ylim)

        if yscale=='log':
            pow = math.floor(math.log10(ylim[0]))+1
            p = pow if (pow%2==0) else (pow+1)
            ytk=[10**p]
            ytk_=10**(p+2)
            while ytk_<=ylim[1]:
                ytk.append(ytk_)
                ytk_*=100.
            ax.set_yticks(ytk)
    ax.set_xlim(np.array(xspace)+ax.get_xlim())

    plt.tight_layout()

    plt.savefig(drop+'/'+name+'.pdf') #saving the figure
    plt.close(fig) #closing out figure
