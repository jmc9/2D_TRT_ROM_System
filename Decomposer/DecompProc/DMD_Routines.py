#!/usr/bin/env python3
#==================================================================================================================================#
#
# DMD_Routines.py is the main python script to process the outputs of Decomposer.c
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Nmodes(dset,dcmp_data):
    N_modes = [];

    if dcmp_data == 'QDf':
        N_modes.append(getattr(dset,'N_modes_Cg'))
        N_modes.append(getattr(dset,'N_modes_fg_avg_xx'))
        N_modes.append(getattr(dset,'N_modes_fg_edgV_xx'))
        N_modes.append(getattr(dset,'N_modes_fg_avg_yy'))
        N_modes.append(getattr(dset,'N_modes_fg_edgH_yy'))
        N_modes.append(getattr(dset,'N_modes_fg_edgV_xy'))
        N_modes.append(getattr(dset,'N_modes_fg_edgH_xy'))

    return N_modes
