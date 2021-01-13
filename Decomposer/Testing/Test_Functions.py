#==================================================================================================================================#
#
# Testing_Functions.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import numpy as np
import math

import Grid_Handling as gh

#==================================================================================================================================#
#
#==================================================================================================================================#
def expf(alpha, tp=[], xp=[], yp=[], xshp='cvec', yshp='none'):

    #checking inputs for errors
    err, errs = expf_incheck(alpha, tp, xp, yp, xshp, yshp)
    if err != 0: #errors were encountered -> print error messages and terminate function call
        x = [print(messg) for messg in errs]
        return 1

    #finding number of exponentials
    N_a = len(alpha)

    #setting up grids over each dimension
    grids = expf_gcalc(N_a, tp, xp, yp, yshp)

    #
    w = w_gen(grids.x.pts, xshp, N_a)

    fdim = tuple(grids.dlist_tr())
    flen = grids.all.len

    f = np.zeros(flen)
    f = expf_calc(f, w, alpha, grids)
    f = np.reshape(f, fdim)

    return f, grids

#==================================================================================================================================#
#
#==================================================================================================================================#
def expf_gcalc(N_a, tp, xp, yp, yshp):
    grids = gh.grids()

    #calculating grid over temporal dimension
    grids.t = gh.grid_calc(N_a * 2, tp)

    #calculating grid over spatial dimension (in x-direction)
    grids.x = gh.grid_calc(N_a, xp)

    #calculating grid over spatial dimension (in y-direction)
    grids.y = gh.grid_calc(N_a, yp, yshp)

    grids.recalc()

    return grids

#==================================================================================================================================#
#
#==================================================================================================================================#
def expf_calc(f, w, alpha, grids):
    N_w = round(len(w)/grids.xy.len)
    for t in range(grids.t.len):
        pw = 0
        pf = t * grids.xy.len
        for j in range(N_w):
            for i in range(grids.xy.len):
                f[pf] = f[pf] + w[pw] * math.exp(alpha[j] * grids.t.pts[t])
                pw += 1
                pf += 1
            pf -= grids.xy.len
        pw -= N_w * grids.xy.len

    return f

#==================================================================================================================================#
#
#==================================================================================================================================#
def expf_incheck(alpha, tp=[], xp=[], yp=[], xshp='cvec', yshp='none'):
    check = 0
    errs = ['**********ERROR(s) DECTECTED**********',
            'The following error(s) were encountered in Test_Functions.expf :']

    if not isinstance(alpha,list):
        check = 1
        errs.append('alpha is not a list')

    if not isinstance(tp,list):
        check = 1
        errs.append('tp is not a list')
    else:
        for t in tp:
            if t < 0.:
                check = 1
                errs.append('at least one element of tp is negative')
                break

    if not isinstance(yp,list):
        check = 1
        errs.append('yp is not a list')
    else:
        for y in yp:
            if y < 0.:
                check = 1
                errs.append('at least one element of yp is negative')
                break

    if not isinstance(xp,list):
        check = 1
        errs.append('xp is not a list')
    else:
        for x in xp:
            if x < 0.:
                check = 1
                errs.append('at least one element of xp is negative')
                break

    if not xshp in ['cvec']:
        check = 1
        errs.append('invalid xshp detected ({})'.format(xshp))

    if not yshp in ['cvec','none']:
        check = 1
        errs.append('invalid yshp detected ({})'.format(xshp))

    return check, errs

#==================================================================================================================================#
#
#==================================================================================================================================#
def w_gen(xp, xshp='cvec', N_w=0, yp=[], yshp='none'):

    N_x = len(xp)

    if yshp == 'none':
        N_y = 0
        N_xy = N_x
    else:
        if len(yp) == 0:
            yp = xp
            N_y = N_x
        else:
            N_y = len(yp)
        N_xy = N_y * N_x

    if N_w == 0:
        N_w = N_xy
    w = np.zeros((N_w * N_xy))

    if xshp == 'cvec':
        if yshp in ['none','cvec']:
            N_v = min([N_xy, N_w])
            p = 0
            for i in range(N_v):
                w[p] = 1.
                p += N_xy + 1

    return w
