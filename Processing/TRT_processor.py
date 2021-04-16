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
#generic python tools
import getopt, sys #handles user inputs from command line
import numpy as np
import os

#netcdf tools
from netCDF4 import Dataset as ncds

#local tools
import inputs as inp
import plotters as pltr
import misc_tools as msc
import processing as proc
import math_tools as mt

#==================================================================================================================================#
#
#==================================================================================================================================#
def TRT_process(infile,proc_dir,plotdir,plt_tfreq):

    if not os.path.isfile(infile):
        print('Input file, '+infile+', does not exist. aborting program')
        quit()

    msc.dirset(proc_dir) #creating 'master' directory to contain all processed data

    #creating directories to put plots of error norms
    norm_dir = proc_dir+'/Norms'; msc.dirset(norm_dir)
    T_norm_dir = norm_dir+'/Temp'; msc.dirset(T_norm_dir)
    E_avg_norm_dir = norm_dir+'/E_avg'; msc.dirset(E_avg_norm_dir)
    fg_avg_xx_norm_dir = norm_dir+'/fg_avg_xx'; msc.dirset(fg_avg_xx_norm_dir)
    fg_avg_yy_norm_dir = norm_dir+'/fg_avg_yy'; msc.dirset(fg_avg_yy_norm_dir)
    fg_EdgV_xx_norm_dir = norm_dir+'/fg_EdgV_xx'; msc.dirset(fg_EdgV_xx_norm_dir)
    fg_EdgV_xy_norm_dir = norm_dir+'/fg_EdgV_xy'; msc.dirset(fg_EdgV_xy_norm_dir)
    fg_EdgH_yy_norm_dir = norm_dir+'/fg_EdgH_yy'; msc.dirset(fg_EdgH_yy_norm_dir)
    fg_EdgH_xy_norm_dir = norm_dir+'/fg_EdgH_xy'; msc.dirset(fg_EdgH_xy_norm_dir)

    cs_dir = proc_dir+'/Cross_Sections'; msc.dirset(cs_dir)
    T_cs_dir = cs_dir+'/Temp'; msc.dirset(T_cs_dir)
    E_cs_dir = cs_dir+'/E_avg'; msc.dirset(E_cs_dir)

    eplt_dir = proc_dir+'/Error_Plots'; msc.dirset(eplt_dir)
    T_eplt_dir = eplt_dir+'/Temp'; msc.dirset(T_eplt_dir)
    E_eplt_dir = eplt_dir+'/E_avg'; msc.dirset(E_eplt_dir)

    (dsets,dsnames,trend_names,fignames,plt_times,plt_tfreq,ntrend_tp,Tbound,Ebound,fg_avg_xx_bnd,
    cs_times,csx_times,csy_times,csx,csy,switch_ranks,switch_norms,switch_cs,switch_sol,switch_errplots,plt_consistency) = inp.input(infile) #reading input file, returning datasets to process
    nsets = len(dsets) #finding the number of datasets

    if switch_sol == 1:
        plot_qdf=False
        plot_E_HO=False
    elif switch_sol == 2:
        plot_qdf=True
        plot_E_HO=False
    else:
        plot_qdf=True
        plot_E_HO=True

    #initializing arrays that will hold simulation data
    (Temp,E_avg,E_edgV,E_edgH,E_avg_MGQD,E_edgV_MGQD,E_edgH_MGQD,E_avg_HO,E_edgV_HO,E_edgH_HO,Fx_edgV,\
    Fy_edgH,Fx_edgV_MGQD,Fy_edgH_MGQD,Fx_edgV_HO,Fy_edgH_HO,Eg_avg,Eg_edgV,Eg_edgH,Eg_avg_HO,Eg_edgV_HO,\
    Eg_edgH_HO,Fxg_edgV,Fyg_edgH,Fxg_edgV_HO,Fyg_edgH_HO,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,\
    fg_edgH_yy,fg_edgH_xy) = proc.init_arr()

    #initializing arrays that hold inputs and directory information
    inp_flags=[]; casedirs=[]; casename=[];
    plotdirs = []; pd_names = []
    for i in range(nsets): plotdirs.append([]); pd_names.append([])

    for n in range(nsets):
        casename.append(dsets[n][:-3]) #obtaining case names from dataset file names (subtracting the '.h5' tail)
        casedirs.append(proc_dir+'/'+casename[n]); msc.dirset(casedirs[n]) #creating case directory
        dsets[n] = ncds(dsets[n],'r') #opening each specified dataset of TRT results
        inp_flags.append(inp.inp_flags(dsets[n])) #checking which variables are given in this dataset

        (plotdirs[n],pd_names[n]) = msc.plotdirs(casedirs[n]) #generating directory paths
        msc.setup_plotdirs(plotdirs[n],pd_names[n],inp_flags[n]) #creating directories
    if not dsnames: dsnames = casename[1:]

    (xp,yp,xp_c,yp_c,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x) = inp.domain_parms(dsets[0]) #getting problem domain parameters
    (A_EdgV, A_EdgH, Ag, Ag_EdgV, Ag_EdgH) = proc.A_calcs(A,N_g,N_y,N_x) #calculating areas/weights for different types of functions
    tpts_ = [ [t*10. for t in ds['tpts'][:] ] for ds in dsets ]

    ###
    #
    tpts=[]
    tn=np.zeros((N_t,nsets),dtype=int)
    pn=0
    p=np.zeros(nsets,dtype=int)
    for tm in tpts_[0]:
        if tm in tpts: continue
        check=True
        p[1:]=0
        pd=1
        for t2 in tpts_[1:]:
            p2=0
            check2=False
            for tm2 in t2:
                if (tm>=tm2-1e-10)and(tm<=tm2+1e-10):
                    check2=True
                    p[pd]=p2
                    break
                p2 += 1
            if not check2:
                check = False
                break
            pd+=1
        if check:
            tpts.append(tm)
            for n in range(nsets): tn[pn][n] = p[n]
            pn += 1
        p[0] += 1
    #
    ###
    N_t2 = pn

    pltbnds = [Tbound,Ebound,Ebound,Ebound]
    bnd = []
    for g in range(N_g):
        boo = True
        j = 0
        for i in fg_avg_xx_bnd[::3]:
            if g==(i-1):
                boo = False
                bnd.append([fg_avg_xx_bnd[j*3+1],fg_avg_xx_bnd[j*3+2]])
                break
            j = j + 1
        if boo: bnd.append([0.,0.])
    pltbnds.append(bnd)

    #finding the instants of time for which to visualize data
    n_tplts=round(N_t/plt_tfreq)
    for i in range(n_tplts):
        j = (i+1)*plt_tfreq*Delt
        if not np.any(np.isclose(tp[tn[t][n]],plt_times)): tp_plt.append(j)

    #if any comp_datasets have been specified, collecting the ranks of their expansions and plotting them
    #currently only working for gsum=1
    if ((nsets>1)and(switch_ranks==1)):
        (time_intervals, ranks, gsum) = proc.rank_collect(dsets,nsets,casename,N_g,proc_dir)
        # rleg = ['$\mathcal{C}$','$\mathcal{F}_{a,xx}$','$\mathcal{F}_{v,xx}$','$\mathcal{F}_{a,yy}$','$\mathcal{F}_{h,yy}$','$\mathcal{F}_{v,xy}$','$\mathcal{F}_{h,xy}$']
        rleg = ['$\mathbf{A}^c_k$','$\mathbf{A}^{f_{xx,c}}_k$','$\mathbf{A}^{f_{xx,v}}_k$','$\mathbf{A}^{f_{yy,c}}_k$','$\mathbf{A}^{f_{yy,h}}_k$','$\mathbf{A}^{f_{xy,v}}_k$','$\mathbf{A}^{f_{xy,h}}_k$']
        for t in range(len(time_intervals)-1):
            pltr.lineplot("Ranks_t{}".format(t+1),trend_names,ranks[t],proc_dir,rleg,yscale='linear',xlabel=r'Truncation Criteria $\xi_{rel}$',ylabel='Dimensionality (k)',marker='D',legloc='upper left',legsize=12)

    #initializing arrays of error norms
    if (switch_norms==1):
        Temp_Norms = np.zeros([6,nsets-1,N_t2])
        E_avg_Norms = np.zeros([6,nsets-1,N_t2])
        fg_avg_xx_Norms = np.zeros([6,nsets-1,N_t2])
        fg_avg_yy_Norms = np.zeros([6,nsets-1,N_t2])
        fg_edgV_xx_Norms = np.zeros([6,nsets-1,N_t2])
        fg_edgV_xy_Norms = np.zeros([6,nsets-1,N_t2])
        fg_edgH_yy_Norms = np.zeros([6,nsets-1,N_t2])
        fg_edgH_xy_Norms = np.zeros([6,nsets-1,N_t2])

        Temp_Norm_Trends = np.zeros([6,len(ntrend_tp),nsets-1])
        E_avg_Norm_Trends = np.zeros([6,len(ntrend_tp),nsets-1])
        ntrend_times = []

    if (switch_cs==1):
        Temp_xcs = np.zeros([3,nsets,len(csx_times),N_x])
        Temp_ycs = np.zeros([3,nsets,len(csy_times),N_y])
        E_xcs = np.zeros([3,nsets,len(csx_times),N_x])
        E_ycs = np.zeros([3,nsets,len(csy_times),N_y])

    t2 = 0
    csx_t = 0
    csy_t = 0
    print('beginning data collection and plotting...')
    for t in range(N_t2):
        T_errs=[]

        print('... simulation time  =  '+str(tp[tn[t][0]]))
        for n in range(nsets):

            if n == 0: m = 0
            else: m = 1

            if ( (switch_norms==1)
               or ((switch_cs==1)and(np.any(np.isclose(tp[tn[t][n]],cs_times))))
               or ( ((switch_sol==1)or(switch_errplots==1)) and ( np.any(np.isclose(tp[tn[t][n]],plt_times)) )) ):
                (Temp[m],E_avg[m],E_edgV[m],E_edgH[m],E_avg_MGQD[m],E_edgV_MGQD[m],E_edgH_MGQD[m],E_avg_HO[m],E_edgV_HO[m],E_edgH_HO[m],Fx_edgV[m],\
                Fy_edgH,Fx_edgV_MGQD[m],Fy_edgH_MGQD[m],Fx_edgV_HO[m],Fy_edgH_HO[m],Eg_avg[m],Eg_edgV[m],Eg_edgH[m],Eg_avg_HO[m],Eg_edgV_HO[m],\
                Eg_edgH_HO[m],Fxg_edgV[m],Fyg_edgH[m],Fxg_edgV_HO[m],Fyg_edgH_HO[m],fg_avg_xx[m],fg_avg_xy[m],fg_avg_yy[m],fg_edgV_xx[m],fg_edgV_xy[m],\
                fg_edgH_yy[m],fg_edgH_xy[m]) = proc.get_data(dsets[n],inp_flags[n],tn[t][n])

            #plotting the results for case n for desired instants of time (contained in the array plt_times)
            if ((switch_sol>0)and( np.any(np.isclose(tp[tn[t][n]],plt_times)) )):
                pltr.plot_results(dsets[n],casedirs[n],xp,yp,tp[tn[t][n]],N_g,N_y,N_x,inp_flags[n],plotdirs[n],pd_names[n],Temp[m],fg_avg_xx[m],\
                fg_avg_xy[m],fg_avg_yy[m],E_avg[m],E_avg_MGQD[m],E_avg_HO[m],Eg_avg[m],Eg_avg_HO[m],pltbnds,plot_qdf,plot_E_HO)

            #if plotting cross-sections of solutions, collecting them here
            if (switch_cs==1):
                if np.any(np.isclose(tp[tn[t][n]],csx_times)):
                    if (csx[0]==1): #'bottom' of domain along x-axis
                        Temp_xcs[0][n][csx_t] = proc.cross_section(Temp[m],(N_y,N_x),(0,-1))
                        E_xcs[0][n][csx_t] = proc.cross_section(E_avg[m],(N_y,N_x),(0,-1))
                    if (csx[1]==1): #'middle' of domain along x-axis
                        Temp_xcs[1][n][csx_t] = proc.cross_section(Temp[m],(N_y,N_x),(N_y//2,-1))
                        E_xcs[1][n][csx_t] = proc.cross_section(E_avg[m],(N_y,N_x),(N_y//2,-1))
                    if (csx[2]==1): #'top' of domain along x-axis
                        Temp_xcs[2][n][csx_t] = proc.cross_section(Temp[m],(N_y,N_x),(N_y-1,-1))
                        E_xcs[2][n][csx_t] = proc.cross_section(E_avg[m],(N_y,N_x),(N_y-1,-1))

                if np.any(np.isclose(tp[tn[t][n]],csy_times)):
                    if (csy[0]==1): #'left' of domain along y-axis
                        Temp_ycs[0][n][csy_t] = proc.cross_section(Temp[m],(N_y,N_x),(-1,0))
                        E_ycs[0][n][csy_t] = proc.cross_section(E_avg[m],(N_y,N_x),(-1,0))
                    if (csy[1]==1): #'middle' of domain along y-axis
                        Temp_ycs[1][n][csy_t] = proc.cross_section(Temp[m],(N_y,N_x),(-1,N_x//2))
                        E_ycs[1][n][csy_t] = proc.cross_section(E_avg[m],(N_y,N_x),(-1,N_x//2))
                    if (csy[2]==1): #'right' of domain along y-axis
                        Temp_ycs[2][n][csy_t] = proc.cross_section(Temp[m],(N_y,N_x),(-1,N_x-1))
                        E_ycs[2][n][csy_t] = proc.cross_section(E_avg[m],(N_y,N_x),(-1,N_x-1))

            if ( (switch_errplots==1) and (n>0) and np.any(np.isclose(tp[tn[t][n]],plt_times)) ):
                err, rerr, ref = proc.errcalc(Temp)
                pltr.plot_data(rerr,'T_err_{}'.format(fignames[n]),xp,yp,tp[tn[t][0]],T_eplt_dir,N_y,N_x,'inferno')

                err, rerr, ref = proc.errcalc(E_avg)
                pltr.plot_data(rerr,'E_err_{}'.format(fignames[n]),xp,yp,tp[tn[t][0]],E_eplt_dir,N_y,N_x,'inferno')

            #if looking at a non-reference dataset, calculate norms of error compared to reference dataset
            if ((switch_norms==1)and(n>0)):
                (Temp_Norms[0][n-1][t],Temp_Norms[1][n-1][t],Temp_Norms[2][n-1][t],Temp_Norms[3][n-1][t],Temp_Norms[4][n-1][t],Temp_Norms[5][n-1][t]) = proc.data_compare(Temp,A)
                (E_avg_Norms[0][n-1][t],E_avg_Norms[1][n-1][t],E_avg_Norms[2][n-1][t],E_avg_Norms[3][n-1][t],E_avg_Norms[4][n-1][t],E_avg_Norms[5][n-1][t]) = proc.data_compare(E_avg,A)
                (fg_avg_xx_Norms[0][n-1][t],fg_avg_xx_Norms[1][n-1][t],fg_avg_xx_Norms[2][n-1][t],fg_avg_xx_Norms[3][n-1][t],fg_avg_xx_Norms[4][n-1][t],fg_avg_xx_Norms[5][n-1][t]) = proc.gdata_compare(fg_avg_xx,Ag)
                (fg_avg_yy_Norms[0][n-1][t],fg_avg_yy_Norms[1][n-1][t],fg_avg_yy_Norms[2][n-1][t],fg_avg_yy_Norms[3][n-1][t],fg_avg_yy_Norms[4][n-1][t],fg_avg_yy_Norms[5][n-1][t]) = proc.gdata_compare(fg_avg_yy,Ag)
                (fg_edgV_xx_Norms[0][n-1][t],fg_edgV_xx_Norms[1][n-1][t],fg_edgV_xx_Norms[2][n-1][t],fg_edgV_xx_Norms[3][n-1][t],fg_edgV_xx_Norms[4][n-1][t],fg_edgV_xx_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgV_xx,Ag_EdgV)
                (fg_edgV_xy_Norms[0][n-1][t],fg_edgV_xy_Norms[1][n-1][t],fg_edgV_xy_Norms[2][n-1][t],fg_edgV_xy_Norms[3][n-1][t],fg_edgV_xy_Norms[4][n-1][t],fg_edgV_xy_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgV_xy,Ag_EdgV)
                (fg_edgH_yy_Norms[0][n-1][t],fg_edgH_yy_Norms[1][n-1][t],fg_edgH_yy_Norms[2][n-1][t],fg_edgH_yy_Norms[3][n-1][t],fg_edgH_yy_Norms[4][n-1][t],fg_edgH_yy_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgH_yy,Ag_EdgH)
                (fg_edgH_xy_Norms[0][n-1][t],fg_edgH_xy_Norms[1][n-1][t],fg_edgH_xy_Norms[2][n-1][t],fg_edgH_xy_Norms[3][n-1][t],fg_edgH_xy_Norms[4][n-1][t],fg_edgH_xy_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgH_xy,Ag_EdgH)

                if t in ntrend_tp:
                    for i in range(6):
                        Temp_Norm_Trends[i][t2][n-1] = Temp_Norms[i][n-1][t]
                        E_avg_Norm_Trends[i][t2][n-1] = E_avg_Norms[i][n-1][t]
        if ((switch_norms==1)and(t in ntrend_tp)):
            ntrend_times.append('{} ns'.format(round(tp[tn[t][0]],8)))
            t2 = t2 + 1
        if ((switch_cs==1)and(np.any(np.isclose(tp[tn[t][n]],csx_times)))): csx_t+=1
        if ((switch_cs==1)and(np.any(np.isclose(tp[tn[t][n]],csy_times)))): csy_t+=1

    #plotting error norms
    if ((switch_norms==1)and(nsets>0)):
        if (plt_consistency==1):
            Tnrm_max = mt.norm_maxes(Temp_Norms)
            Enrm_max = mt.norm_maxes(E_avg_Norms)
            TE_upbnd = [max(Tmax,Emax) for (Tmax,Emax) in zip(Tnrm_max,Enrm_max)]

            Tnrm_min = mt.norm_mins(Temp_Norms)
            Enrm_min = mt.norm_mins(E_avg_Norms)
            TE_lowbnd = [min(Tmin,Emin) for (Tmin,Emin) in zip(Tnrm_min,Enrm_min)]

            TE_bnd = [[low,up] for (low,up) in zip(TE_lowbnd,TE_upbnd)]
            TE_bnd = [[mt.logbnd_low(low,5.), mt.logbnd_up(up,5.)] for (low,up) in zip(TE_lowbnd,TE_upbnd)]
        else:
            TE_bnd = []

        pltr.plot_norms("T",Temp_Norms,tpts,dsnames,T_norm_dir,'Time (ns)',vname='T',bnds=TE_bnd)
        pltr.plot_norms("E_avg",E_avg_Norms,tpts,dsnames,E_avg_norm_dir,'Time (ns)',vname='E',bnds=TE_bnd)
        pltr.plot_norms("fg_avg_xx",fg_avg_xx_Norms,tpts,dsnames,fg_avg_xx_norm_dir,'Time (ns)',vname='fxx')
        pltr.plot_norms("fg_avg_yy",fg_avg_yy_Norms,tpts,dsnames,fg_avg_yy_norm_dir,'Time (ns)',vname='fyy')
        pltr.plot_norms("fg_edgV_xx",fg_edgV_xx_Norms,tpts,dsnames,fg_EdgV_xx_norm_dir,'Time (ns)',vname='fxx')
        pltr.plot_norms("fg_edgV_xy",fg_edgV_xy_Norms,tpts,dsnames,fg_EdgV_xy_norm_dir,'Time (ns)',vname='fxy')
        pltr.plot_norms("fg_edgH_yy",fg_edgH_yy_Norms,tpts,dsnames,fg_EdgH_yy_norm_dir,'Time (ns)',vname='fyy')
        pltr.plot_norms("fg_edgH_xy",fg_edgH_xy_Norms,tpts,dsnames,fg_EdgH_xy_norm_dir,'Time (ns)',vname='fxy')

        pltr.plot_norms("T_trends",Temp_Norm_Trends,trend_names,ntrend_times,T_norm_dir,r'Truncation Criteria $\xi_{rel}$',vname='T',marker='D',bnds=TE_bnd)
        pltr.plot_norms("E_trends",E_avg_Norm_Trends,trend_names,ntrend_times,E_avg_norm_dir,r'Truncation Criteria $\xi_{rel}$',vname='E',marker='D',bnds=TE_bnd)

    #plotting solution cross sections
    if ((switch_cs==1)and(cs_times)):
        dsnames.insert(0,'fom')
        if (csx[0]==1):
            pltr.plot_cs('Temp_xcs_bot',Temp_xcs[0],xp_c,csx_times,dsnames,fignames,T_cs_dir,'x-length (cm)','Temperature (eV)',textp=[0,3,5,7,5])
            pltr.plot_cs('E_xcs_bot',E_xcs[0],xp_c,cs_times,dsnames,fignames,E_cs_dir,'x-length (cm)','Radiation Energy Density (erg$\cdot 1\\times 10^{13}$)',textp=[0,3,5,7,5])
        if (csx[1]==1):
            pltr.plot_cs('Temp_xcs_mid',Temp_xcs[1],xp_c,csx_times,dsnames,fignames,T_cs_dir,'x-length (cm)','Temperature (eV)',textp=[0,3,5,7,5])
            pltr.plot_cs('E_xcs_mid',E_xcs[1],xp_c,cs_times,dsnames,fignames,E_cs_dir,'x-length (cm)','Radiation Energy Density (erg$\cdot 1\\times 10^{13}$)',textp=[0,3,5,7,5])
        if (csx[2]==1):
            pltr.plot_cs('Temp_xcs_top',Temp_xcs[2],xp_c,csx_times,dsnames,fignames,T_cs_dir,'x-length (cm)','Temperature (eV)',textp=[0,3,5,7,5])
            pltr.plot_cs('E_xcs_top',E_xcs[2],xp_c,cs_times,dsnames,fignames,E_cs_dir,'x-length (cm)','Radiation Energy Density (erg$\cdot 1\\times 10^{13}$)',textp=[0,3,5,7,5])

        if (csy[0]==1):
            pltr.plot_cs('Temp_ycs_left',Temp_ycs[0],yp_c,csy_times,dsnames,fignames,T_cs_dir,'y-length (cm)','Temperature (eV)',textp=np.full(len(cs_times),N_y//2))
            pltr.plot_cs('E_ycs_left',E_ycs[0],yp_c,csy_times,dsnames,fignames,E_cs_dir,'y-length (cm)','Radiation Energy Density (erg$\cdot 1\\times 10^{13}$)',textp=np.full(len(cs_times),N_y//2))
        if (csy[1]==1):
            pltr.plot_cs('Temp_ycs_mid',Temp_ycs[1],yp_c,csy_times,dsnames,fignames,T_cs_dir,'y-length (cm)','Temperature (eV)',textp=np.full(len(cs_times),N_y//2))
            pltr.plot_cs('E_ycs_mid',E_ycs[1],yp_c,csy_times,dsnames,fignames,E_cs_dir,'y-length (cm)','Radiation Energy Density (erg$\cdot 1\\times 10^{13}$)',textp=np.full(len(cs_times),N_y//2))
        if (csy[2]==1):
            pltr.plot_cs('Temp_ycs_right',Temp_ycs[2],yp_c,csy_times,dsnames,fignames,T_cs_dir,'y-length (cm)','Temperature (eV)',textp=np.full(len(cs_times),N_y//2))
            pltr.plot_cs('E_ycs_right',E_ycs[2],yp_c,csy_times,dsnames,fignames,E_cs_dir,'y-length (cm)','Radiation Energy Density (erg$\cdot 1\\times 10^{13}$)',textp=np.full(len(cs_times),N_y//2))


    for n in range(nsets): dsets[n].close() #closing dataset file

#==================================================================================================================================#
#
#==================================================================================================================================#
def options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hzi:o:', ['help','input=','output='])
    except getopt.GetoptError as err:
        print(str(err))
        usage(); quit()

    (infile,proc_dir,plotdir,plt_tfreq) = defaults()

    for o, a in opts:
        if o in ['-h','--help']:
            usage(); quit()
        elif o in ['-z']:
            inp.sample_infile(); quit()
        elif o in ['-i','--input']:
            infile = a
        elif o in ['-o','--output']:
            proc_dir = a

    return (infile,proc_dir,plotdir,plt_tfreq)

#==================================================================================================================================#
#
#==================================================================================================================================#
def defaults():
    infile = 'input.inp'
    proc_dir = 'processed_data'
    plotdir = 'plots'
    plt_tfreq = 1

    return (infile,proc_dir,plotdir,plt_tfreq)

#==================================================================================================================================#
#
#==================================================================================================================================#
def usage():
    print()
    print('- TRT_processor')
    print('- Joseph M. Coale, jmcoale@ncsu.edu')
    print()
    print('This python3 code is used to process and compare results generated by the 2D TRT ROM code (Author: Joseph M. Coale)')
    print()
    print('If you would like a sample input file with descriptions of each input option, please call this code again with option -z\n--> i.e. call ./TRT_processor.py -z')
    print()
    print('Below are all available options and their usages')
    print()
    print('Options:')
    print('-h [help]: prints this message')
    print('-i [input]: used to specify an input file')
    print('-o [output]: used to specify the output directory')
    print()

#==================================================================================================================================#
# MAIN PROGRAM
# This portion of the code is run when this file is the main program
# i.e. when the call is made './TRT_processor.py'
#==================================================================================================================================#
if __name__ == '__main__':
    (infile,proc_dir,plotdir,plt_tfreq) = options()
    TRT_process(infile,proc_dir,plotdir,plt_tfreq)
