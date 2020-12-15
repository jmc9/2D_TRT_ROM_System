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

#local tools
import DMD_Routines as dmdr
import ToolBox as tb

#==================================================================================================================================#
#
#==================================================================================================================================#
def Process_Data(dset,Dcmp_Data,plt_modes,plotdir):

    N_data = len(Dcmp_Data)
    for i in range(N_data):
        dcmp = Dcmp_Data[i]

        dir = '{}/{}'.format(plotdir,dcmp.name)
        tb.dirset(dir)

        if (Dcmp_Data[i].type in ['DMD','DMDg']):
            (dcmp.rank, dcmp.clen) = dmdr.Read_Dims(dset,dcmp.name)
            dcmp.dat = dmdr.Read_DMD(dset,dcmp.name)

            if dcmp.opt[0] == 0:
                dmdr.Plot_DMD(dcmp,dir,plt_modes)
            else:
                dmdr.Plot_DMD(dcmp,dir,plt_modes,evecs=False)
