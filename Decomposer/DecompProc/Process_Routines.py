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
def Process_DMD(dset,plotdir,clen_avg,clen_edgV,clen_edgH,BClen,rank_avg,rank_edgV,rank_edgH,rank_BC,N_g,N_y,N_x,xp,yp):
    dir = '{}/{}'.format(plotdir,'fg_avg_xx'); tb.dirset(dir)
    dmdr.Plot_DMD(dset,'fg_avg_xx',range(10),clen_avg,rank_avg,N_g,N_y,N_x,xp,yp,dir)

    dir = '{}/{}'.format(plotdir,'fg_avg_yy'); tb.dirset(dir)
    dmdr.Plot_DMD(dset,'fg_avg_yy',range(10),clen_avg,rank_avg,N_g,N_y,N_x,xp,yp,dir)

    dir = '{}/{}'.format(plotdir,'BCg'); tb.dirset(dir)
    dmdr.Plot_DMD(dset,'BCg',[],BClen,rank_BC,N_g,N_y,N_x,xp,yp,dir,evecs=False)

    dir = '{}/{}'.format(plotdir,'fg_edgV_xx'); tb.dirset(dir)
    dmdr.Plot_DMD(dset,'fg_edgV_xx',[],clen_edgV,rank_edgV,N_g,N_y,N_x,xp,yp,dir,evecs=False)

    dir = '{}/{}'.format(plotdir,'fg_edgV_xy'); tb.dirset(dir)
    dmdr.Plot_DMD(dset,'fg_edgV_xy',[],clen_edgV,rank_edgV,N_g,N_y,N_x,xp,yp,dir,evecs=False)

    dir = '{}/{}'.format(plotdir,'fg_edgH_yy'); tb.dirset(dir)
    dmdr.Plot_DMD(dset,'fg_edgH_yy',[],clen_edgH,rank_edgH,N_g,N_y,N_x,xp,yp,dir,evecs=False)

    dir = '{}/{}'.format(plotdir,'fg_edgH_xy'); tb.dirset(dir)
    dmdr.Plot_DMD(dset,'fg_edgH_xy',[],clen_edgH,rank_edgH,N_g,N_y,N_x,xp,yp,dir,evecs=False)
