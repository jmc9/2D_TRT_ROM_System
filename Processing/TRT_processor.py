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
def TRT_process(infile,proc_dir,plotdir):

    msc.dirset(proc_dir) #creating 'master' directory to contain all processed data

    dsets = inp.input(infile) #reading input file, returning datasets to process
    nsets = len(dsets) #finding the number of datasets

    (Temp,E_avg,E_edgV,E_edgH,E_avg_MGQD,E_edgV_MGQD,E_edgH_MGQD,E_avg_HO,E_edgV_HO,E_edgH_HO,Fx_edgV,\
    Fy_edgH,Fx_edgV_MGQD,Fy_edgH_MGQD,Fx_edgV_HO,Fy_edgH_HO,Eg_avg,Eg_edgV,Eg_edgH,Eg_avg_HO,Eg_edgV_HO,\
    Eg_edgH_HO,Fxg_edgV,Fyg_edgH,Fxg_edgV_HO,Fyg_edgH_HO,fg_avg_xx,fg_avg_xy,fg_avg_yy,fg_edgV_xx,fg_edgV_xy,\
    fg_edgH_yy,fg_edgH_xy) = proc.init_arr()

    inp_flags=[]; casedirs=[]
    plotdirs = []; pd_names = []
    for i in range(nsets): plotdirs.append([]); pd_names.append([])

    for n in range(nsets):
        casename = dsets[n][:-3] #obtaining case names from dataset file names (subtracting the '.h5' tail)
        casedirs.append(proc_dir+'/'+casename); msc.dirset(casedirs[n]) #creating case directory
        dsets[n] = ncds(dsets[n],'r') #opening each specified dataset of TRT results
        inp_flags.append(inp.inp_flags(dsets[n])) #checking which variables are given in this dataset

        (plotdirs[n],pd_names[n]) = msc.plotdirs(casedirs[n])
        msc.setup_plotdirs(plotdirs[n],pd_names[n],inp_flags[n])

    (xp,yp,tp,Delx,Dely,Delt,A) = inp.domain_parms(dsets[0]) #getting problem domain parameters
    N_t = len(tp)

    Temp_Norms = np.zeros([6,nsets,N_t])
    E_avg_Norms = np.zeros([6,nsets,N_t])
    for t in range(N_t):
        for n in range(nsets):

            if n == 0: m = 0
            else: m = 1

            (Temp[m],E_avg[m],E_edgV[m],E_edgH[m],E_avg_MGQD[m],E_edgV_MGQD[m],E_edgH_MGQD[m],E_avg_HO[m],E_edgV_HO[m],E_edgH_HO[m],Fx_edgV[m],\
            Fy_edgH,Fx_edgV_MGQD[m],Fy_edgH_MGQD[m],Fx_edgV_HO[m],Fy_edgH_HO[m],Eg_avg[m],Eg_edgV[m],Eg_edgH[m],Eg_avg_HO[m],Eg_edgV_HO[m],\
            Eg_edgH_HO[m],Fxg_edgV[m],Fyg_edgH[m],Fxg_edgV_HO[m],Fyg_edgH_HO[m],fg_avg_xx[m],fg_avg_xy[m],fg_avg_yy[m],fg_edgV_xx[m],fg_edgV_xy[m],\
            fg_edgH_yy[m],fg_edgH_xy[m]) = proc.get_data(dsets[n],inp_flags[n],t)

            pltr.plot_results(dsets[n],casedirs[n],xp,yp,tp[t],inp_flags[n],plotdirs[n],pd_names[n],Temp[m],fg_avg_xx[m],\
            fg_avg_xy[m],fg_avg_yy[m],E_avg[m],E_avg_MGQD[m],E_avg_HO[m],Eg_avg[m],Eg_avg_HO[m]) #plotting the results for case n

            if n>0:
                (Temp_Norms[0][n][t],Temp_Norms[1][n][t],Temp_Norms[2][n][t],Temp_Norms[3][n][t],Temp_Norms[4][n][t],Temp_Norms[5][n][t]) = proc.data_compare(Temp,A)
                (E_avg_Norms[0][n][t],E_avg_Norms[1][n][t],E_avg_Norms[2][n][t],E_avg_Norms[3][n][t],E_avg_Norms[4][n][t],E_avg_Norms[5][n][t]) = proc.data_compare(E_avg,A)

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

    (infile,proc_dir,plotdir) = defaults()

    for o, a in opts:
        if o in ['-h','--help']:
            usage(); quit()
        elif o in ['-i','--input']:
            infile = a
        elif o in ['-o','--output']:
            proc_dir = a

    return (infile,proc_dir,plotdir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def defaults():
    infile = 'input.inp'
    proc_dir = 'processed_data'
    plotdir = 'plots'

    return (infile,proc_dir,plotdir)

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
    (infile,proc_dir,plotdir) = options()
    TRT_process(infile,proc_dir,plotdir)
