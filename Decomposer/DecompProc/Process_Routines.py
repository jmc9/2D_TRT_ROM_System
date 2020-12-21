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
import numpy as np

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
                    if Prob_Data[i].opt[0] == 0:
                        coef = dmdr.Coef_Calc(dcmp,Prob_Data[i])

                        len1 = len(Prob_Data[i].dat)
                        len2 = len(Prob_Data[i].dat[0])
                        edat = np.zeros([len1, len2],dtype=np.complex_)
                        for j in range(len1):
                            time = Prob_Data[i].grids[0].crds[j]
                            edat[j] = dmdr.Expand(dcmp, coef, time)

            else:
                dmdr.Plot_DMD(dcmp,dir,plt_modes,evecs=False)
