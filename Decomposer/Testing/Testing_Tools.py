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
from Test_Functions import expf
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
def Run_Test(exec_dir, decomp_perphs, proc_perphs, drop, alpha=[.1], mesh=[[],[],[]], shp=['cvec', 'none']):

    #
    tp = []
    xp = []
    yp = []
    if isinstance(mesh, list):
        ml = len(mesh)
        if ml > 0:
            if isinstance(mesh[0], list):
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
    if isinstance(shp, list):
        sl = len(shp)
        xshp = shp[0]
        if sl > 1: yshp = shp[1]
    elif isinstance (shp, str):
        xshp = shp
    else:
        quit()
        
    #
    f, fgrids = expf(alpha, tp, xp, yp, xshp, yshp)

    #
    outfile = 'test.h5'
    outf(f, fgrids, outfile)

    #
    ds_exec(decomp_perphs, outfile, drop, exec_dir)
    dp_exec(proc_perphs, decomp_perphs.out, outfile, drop, exec_dir)

    return
