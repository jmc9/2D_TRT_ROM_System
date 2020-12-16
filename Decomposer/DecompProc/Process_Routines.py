#!/usr/bin/env python3
#==================================================================================================================================#
#
# Process_Routines.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import math

#local tools
import DMD_Routines as dmdr
import ToolBox as tb

#==================================================================================================================================#
#
#==================================================================================================================================#
def Process_Data(dset,Dcmp_Data,Prob_Data,plt_modes,plotdir):

    N_data = len(Dcmp_Data)
    for i in range(N_data):
        dcmp = Dcmp_Data[i]

        dir = '{}/{}'.format(plotdir,dcmp.name)
        tb.dirset(dir)

        if (dcmp.type in ['DMD','DMDg']):
            (dcmp.rank, dcmp.clen) = dmdr.Read_Dims(dset,dcmp.name)
            dcmp.dat = dmdr.Read_DMD(dset,dcmp.name)

            dcmp = dmdr.DMD_Sort(dcmp)

            if dcmp.opt[0] == 0:
                dmdr.Plot_DMD(dcmp,dir,plt_modes)

                if len(Prob_Data) != 0:
                    coef = dmdr.Coef_Calc(dcmp,Prob_Data[i])

            else:
                dmdr.Plot_DMD(dcmp,dir,plt_modes,evecs=False)
