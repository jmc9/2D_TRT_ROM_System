#==================================================================================================================================#
#
# Testing_Tools.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#
import numpy as np
from Test_Functions import expf
from Test_Functions import ak_nonlf
from Test_Functions import trt_1df
from Test_Functions import trt_2df
from IO_Functions import output_f as outf
from DecompSrc_Handling import exec_decompsrc as ds_exec
from DecompProc_Handling import exec_decompproc as dp_exec

#==================================================================================================================================#
#
#==================================================================================================================================#
class peripherals():

    path = ''
    inp  = ''
    out  = ''
    log  = ''

    def __init__(self, path, inp='', out='', log=''):
        self.load(path=path, inp=inp, out=out, log=log)
        return

    def load(self, inp='0', out='0', log='0', path='0'):
        if path != '0': self.path = path
        if inp  != '0': self.inp  = inp
        if out  != '0': self.out  = out
        if log  != '0': self.log  = log
        return

#==================================================================================================================================#
#
#==================================================================================================================================#
def Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, drop, alpha=[.1], mesh=[[],[],[]], shp=['cvec', 'none'], scale=[], shift=[], eps=1e-12):

    #
    tp = []
    xp = []
    yp = []
    if isinstance(mesh, list):
        ml = len(mesh)
        if ml > 0:
            if isinstance(mesh[0], (list, np.ndarray)):
                tp = mesh[0]
                if ml > 1: xp = mesh[1]
                if ml > 2: yp = mesh[2]
            else:
                tp = mesh
    else:
        quit()

    #
    xshp = 'cvec'
    yshp = 'none'
    if isinstance(shp, (list, np.ndarray)):
        sl = len(shp)
        xshp = shp[0]
        if sl > 1: yshp = shp[1]
    elif isinstance (shp, str):
        xshp = shp
    else:
        quit()

    #
    f, fgrids = expf(alpha, tp, xp, yp, xshp, yshp, scale, shift)

    #
    outfile = 'test.h5'
    outf(f, fgrids, outfile)

    #
    ds_exec(decomp_perphs, outfile, drop, exec_dir, eps)
    dp_exec(proc_perphs, decomp_perphs.out, outfile, drop, exec_dir)

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def Run_Test_nak(exec_dir, decomp_perphs, proc_perphs, drop, tp, xp):
    #
    f, fgrids = ak_nonlf(tp, xp)

    #
    outfile = 'test.h5'
    outf(f, fgrids, outfile)

    #
    ds_exec(decomp_perphs, outfile, drop, exec_dir)
    dp_exec(proc_perphs, decomp_perphs.out, outfile, drop, exec_dir)

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, drop, file_name, data_name, mg_opt='', dcmp_type='DMD'):
    #
    f, fgrids = trt_1df(file_name, data_name, mg_opt)

    #
    outfile = 'test.h5'
    outf(f, fgrids, outfile)

    #
    ds_exec(decomp_perphs, outfile, drop, exec_dir, dcmp_type=dcmp_type)
    dp_exec(proc_perphs, decomp_perphs.out, outfile, drop, exec_dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, drop, file_name, data_name, tstart=0., tend=.6, dcmp_type='DMD', group=0, trunc_eps=1., trunc_opt=2, eps=1e-12):
    #
    f, fgrids = trt_2df(file_name, data_name, tstart, tend, group)

    #
    outfile = 'test.h5'
    outf(f, fgrids, outfile)

    #
    ds_exec(decomp_perphs, outfile, drop, exec_dir, dcmp_type=dcmp_type, eps=eps)
    dp_exec(proc_perphs, decomp_perphs.out, outfile, drop, exec_dir, trunc_eps=trunc_eps, trunc_opt=trunc_opt)
