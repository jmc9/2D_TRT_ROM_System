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

    dsets = inp.input(infile) #reading input file, returning datasets to process
    nsets = len(dsets) #finding the number of datasets

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

    (xp,yp,tp,Delx,Dely,Delt,A,N_t,N_g,N_y,N_x) = inp.domain_parms(dsets[0]) #getting problem domain parameters
    (A_EdgV, A_EdgH, Ag, Ag_EdgV, Ag_EdgH) = proc.A_calcs(A,N_g,N_y,N_x) #calculating areas/weights for different types of functions

    #finding the instants of time for which to visualize data
    tp_plt=[0]; n_tplts=round(N_t/plt_tfreq)
    for i in range(n_tplts): tp_plt.append((i+1)*plt_tfreq-1)
    tp_plt=[1000]

    #initializing arrays of error norms
    Temp_Norms = np.zeros([6,nsets-1,N_t])
    E_avg_Norms = np.zeros([6,nsets-1,N_t])
    fg_avg_xx_Norms = np.zeros([6,nsets-1,N_t])
    fg_avg_yy_Norms = np.zeros([6,nsets-1,N_t])
    fg_edgV_xx_Norms = np.zeros([6,nsets-1,N_t])
    fg_edgV_xy_Norms = np.zeros([6,nsets-1,N_t])
    fg_edgH_yy_Norms = np.zeros([6,nsets-1,N_t])
    fg_edgH_xy_Norms = np.zeros([6,nsets-1,N_t])

    print('beginning data collection and plotting...')
    for t in range(N_t):
        print('... simulation time  =  '+str(tp[t]))
        for n in range(nsets):

            if n == 0: m = 0
            else: m = 1

            #collecting data from dataset
            (Temp[m],E_avg[m],E_edgV[m],E_edgH[m],E_avg_MGQD[m],E_edgV_MGQD[m],E_edgH_MGQD[m],E_avg_HO[m],E_edgV_HO[m],E_edgH_HO[m],Fx_edgV[m],\
            Fy_edgH,Fx_edgV_MGQD[m],Fy_edgH_MGQD[m],Fx_edgV_HO[m],Fy_edgH_HO[m],Eg_avg[m],Eg_edgV[m],Eg_edgH[m],Eg_avg_HO[m],Eg_edgV_HO[m],\
            Eg_edgH_HO[m],Fxg_edgV[m],Fyg_edgH[m],Fxg_edgV_HO[m],Fyg_edgH_HO[m],fg_avg_xx[m],fg_avg_xy[m],fg_avg_yy[m],fg_edgV_xx[m],fg_edgV_xy[m],\
            fg_edgH_yy[m],fg_edgH_xy[m]) = proc.get_data(dsets[n],inp_flags[n],t)

            #plotting the results for case n for desired instants of time (contained in the array tp_plt)
            if (t in tp_plt):
                pltr.plot_results(dsets[n],casedirs[n],xp,yp,tp[t],N_g,N_y,N_x,inp_flags[n],plotdirs[n],pd_names[n],Temp[m],fg_avg_xx[m],\
                fg_avg_xy[m],fg_avg_yy[m],E_avg[m],E_avg_MGQD[m],E_avg_HO[m],Eg_avg[m],Eg_avg_HO[m])

            #if looking at a non-reference dataset, calculate norms of error compared to reference dataset
            if n>0:
                (Temp_Norms[0][n-1][t],Temp_Norms[1][n-1][t],Temp_Norms[2][n-1][t],Temp_Norms[3][n-1][t],Temp_Norms[4][n-1][t],Temp_Norms[5][n-1][t]) = proc.data_compare(Temp,A)
                (E_avg_Norms[0][n-1][t],E_avg_Norms[1][n-1][t],E_avg_Norms[2][n-1][t],E_avg_Norms[3][n-1][t],E_avg_Norms[4][n-1][t],E_avg_Norms[5][n-1][t]) = proc.data_compare(E_avg,A)
                (fg_avg_xx_Norms[0][n-1][t],fg_avg_xx_Norms[1][n-1][t],fg_avg_xx_Norms[2][n-1][t],fg_avg_xx_Norms[3][n-1][t],fg_avg_xx_Norms[4][n-1][t],fg_avg_xx_Norms[5][n-1][t]) = proc.gdata_compare(fg_avg_xx,Ag)
                (fg_avg_yy_Norms[0][n-1][t],fg_avg_yy_Norms[1][n-1][t],fg_avg_yy_Norms[2][n-1][t],fg_avg_yy_Norms[3][n-1][t],fg_avg_yy_Norms[4][n-1][t],fg_avg_yy_Norms[5][n-1][t]) = proc.gdata_compare(fg_avg_yy,Ag)
                (fg_edgV_xx_Norms[0][n-1][t],fg_edgV_xx_Norms[1][n-1][t],fg_edgV_xx_Norms[2][n-1][t],fg_edgV_xx_Norms[3][n-1][t],fg_edgV_xx_Norms[4][n-1][t],fg_edgV_xx_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgV_xx,Ag_EdgV)
                (fg_edgV_xy_Norms[0][n-1][t],fg_edgV_xy_Norms[1][n-1][t],fg_edgV_xy_Norms[2][n-1][t],fg_edgV_xy_Norms[3][n-1][t],fg_edgV_xy_Norms[4][n-1][t],fg_edgV_xy_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgV_xy,Ag_EdgV)
                (fg_edgH_yy_Norms[0][n-1][t],fg_edgH_yy_Norms[1][n-1][t],fg_edgH_yy_Norms[2][n-1][t],fg_edgH_yy_Norms[3][n-1][t],fg_edgH_yy_Norms[4][n-1][t],fg_edgH_yy_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgH_yy,Ag_EdgH)
                (fg_edgH_xy_Norms[0][n-1][t],fg_edgH_xy_Norms[1][n-1][t],fg_edgH_xy_Norms[2][n-1][t],fg_edgH_xy_Norms[3][n-1][t],fg_edgH_xy_Norms[4][n-1][t],fg_edgH_xy_Norms[5][n-1][t]) = proc.gdata_compare(fg_edgH_xy,Ag_EdgH)

    #plotting error norms
    pltr.plot_norms("T",Temp_Norms,tp,casename[1:],T_norm_dir)
    pltr.plot_norms("E_avg",E_avg_Norms,tp,casename[1:],E_avg_norm_dir)
    pltr.plot_norms("fg_avg_xx",fg_avg_xx_Norms,tp,casename[1:],fg_avg_xx_norm_dir)
    pltr.plot_norms("fg_avg_yy",fg_avg_yy_Norms,tp,casename[1:],fg_avg_yy_norm_dir)
    pltr.plot_norms("fg_edgV_xx",fg_edgV_xx_Norms,tp,casename[1:],fg_EdgV_xx_norm_dir)
    pltr.plot_norms("fg_edgV_xy",fg_edgV_xy_Norms,tp,casename[1:],fg_EdgV_xy_norm_dir)
    pltr.plot_norms("fg_edgH_yy",fg_edgH_yy_Norms,tp,casename[1:],fg_EdgH_yy_norm_dir)
    pltr.plot_norms("fg_edgH_xy",fg_edgH_xy_Norms,tp,casename[1:],fg_EdgH_xy_norm_dir)

    for n in range(nsets): dsets[n].close() #closing dataset file

#==================================================================================================================================#
#
#==================================================================================================================================#
def options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'hi:o:', ['help','input=','output='])
    except getopt.GetoptError as err:
        print(str(err))
        usage(); quit()

    (infile,proc_dir,plotdir,plt_tfreq) = defaults()

    for o, a in opts:
        if o in ['-h','--help']:
            usage(); quit()
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
    plt_tfreq = 10

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
