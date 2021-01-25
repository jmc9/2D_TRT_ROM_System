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
def ak_nonlf(tp, xp):

    grids = gh.grids()
    grids.t = gh.grid_calc(gp=tp)
    grids.x = gh.grid_calc(gp=xp)

    f = np.zeros(grids.t.len * grids.x.len)
    p = 0
    for t in grids.t.pts:
        for x in grids.x.pts:
            f[p] = (1. - x) * math.cos( 3. * math.pi * t * ( x + 1. ) ) * math.exp( -( 1. + x ) * t )
            p += 1
    f = np.reshape(f, (grids.t.len, grids.x.len))

    return f, grids

#==================================================================================================================================#
#
#==================================================================================================================================#
def expf(alpha, tp=[], xp=[], yp=[], xshp='cvec', yshp='none', scale=[], shift=[]):

    #checking inputs for errors
    err, errs = expf_incheck(alpha, tp, xp, yp, xshp, yshp, scale, shift)
    if err != 0: #errors were encountered -> print error messages and terminate function call
        print()
        x = [print(messg) for messg in errs]
        print()
        return 1

    #finding number of exponentials
    N_a = len(alpha)

    #
    if isinstance(xshp,str):
        shp = xshp
        xshp = []
        for i in range(N_a): xshp.append(shp)
    if isinstance(yshp,str):
        shp = yshp
        yshp = []
        for i in range(N_a): yshp.append(shp)

    #setting up grids over each dimension
    grids = expf_gcalc(N_a, tp, xp, yp, yshp[0])

    #
    w = w_gen(grids.x.pts, xshp, N_a, grids.y.pts, yshp, scale, shift)
    # print(np.resize(w,(N_a,grids.x.len)))

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

    return f

#==================================================================================================================================#
#
#==================================================================================================================================#
def expf_incheck(alpha, tp=[], xp=[], yp=[], xshp='cvec', yshp='none', scale=[], shift=[]):
    check = 0
    errs = ['**********ERROR(s) DECTECTED**********',
            'The following error(s) were encountered in Test_Functions.expf :']

    #
    if not isinstance(alpha,(list, np.ndarray)):
        check = 1
        errs.append('alpha is not a list')

    #
    if not isinstance(tp,(list, np.ndarray)):
        check = 1
        errs.append('tp is not a list')
    else:
        for t in tp:
            if t < 0.:
                check = 1
                errs.append('at least one element of tp is negative')
                break

    #
    if not isinstance(yp,(list, np.ndarray)):
        check = 1
        errs.append('yp is not a list')
    else:
        for y in yp:
            if y < 0.:
                check = 1
                errs.append('at least one element of yp is negative')
                break

    #
    if not isinstance(xp,(list, np.ndarray)):
        check = 1
        errs.append('xp is not a list')
    else:
        for x in xp:
            if x < 0.:
                check = 1
                errs.append('at least one element of xp is negative')
                break

    #
    if isinstance(xshp,str):
        if not xshp in ['cvec','linear','lin','const']:
            check = 1
            errs.append('invalid xshp detected ({})'.format(xshp))
    elif isinstance(xshp,(list, np.ndarray)):
        if len(xshp) < len(alpha):
            check = 1
            errs.append('len(xshp) < len(alpha)')
        for shp in xshp:
            if not shp in ['cvec','linear','lin','const']:
                check = 1
                errs.append('invalid xshp detected ({})'.format(shp))
    else:
        check = 1
        errs.append('xshp is not a string nor a list')

    #
    if isinstance(yshp,str):
        if not yshp in ['none','cvec','linear','lin','const']:
            check = 1
            errs.append('invalid yshp detected ({})'.format(yshp))
    elif isinstance(yshp,(list, np.ndarray)):
        if len(yshp) < len(alpha):
            check = 1
            errs.append('len(yshp) < len(alpha)')
        for shp in yshp:
            if not shp in ['none','cvec','linear','lin','const']:
                check = 1
                errs.append('invalid yshp detected ({})'.format(shp))
        if 'none' in yshp:
            for shp in yshp:
                if shp != 'none':
                    check = 1
                    errs.append('yshp list cannot mix "none" with other shapes')
                    break
    else:
        check = 1
        errs.append('yshp is not a string nor a list')

    #
    if isinstance(scale,(list, np.ndarray)):
        if (len(scale) < len(alpha))and(len(scale) != 0):
            check = 1
            errs.append('len(scale) < len(alpha)')

    #
    if isinstance(shift,(list, np.ndarray)):
        if (len(shift) < len(alpha))and(len(shift) != 0):
            check = 1
            errs.append('len(shift) < len(alpha)')

    #
    return check, errs

#==================================================================================================================================#
#
#==================================================================================================================================#
def w_gen(xp, xshp=['cvec'], N_w=1, yp=[], yshp=['none'], scale=[], shift=[]):

    #
    if len(scale) == 0:
        scale = []
        for i in range(N_w):
            scale.append([1., 1.])
    else:
        scale = scale.tolist()
        for i in range(len(scale)):
            if not isinstance(scale[i],(list, np.ndarray)):
                scale[i] = [scale[i], scale[i]]

    #
    if len(shift) == 0:
        shift = []
        for i in range(N_w):
            shift.append(0.)

    #
    N_x = len(xp)
    N_y = len(yp)
    N_xy = (N_y * N_x) if (N_y > 0) else N_x
    w = np.zeros((N_w * N_xy))

    p = 0
    for i in range(N_w):
        if (yshp[i] != 'none')and(N_y == 0): yshp[i] = 'none'

        if yshp[i] == 'none':
            if xshp[i] == 'cvec':
                w[p + i] = scale[i][0] + shift[i]
                p += N_x

            elif xshp[i] in ['linear','lin','const']:
                if xshp[i] == 'const': xvals = np.full(N_x, scale[i][0])
                elif xshp[i] in ['linear','lin']: xvals = np.linspace(0., scale[i][0], N_x)

                for j in range(N_x):
                    w[p] = xvals[j] + shift[i]
                    p += 1

        elif yshp[i] == 'cvec':
            if xshp[i] == 'cvec':
                w[p + i] = scale[i][0] + shift[i]
                p += N_xy

            elif xshp[i] == 'const':
                p += i * N_x
                for j in range(N_x):
                    w[p] = scale[i][0] + shift[i]
                    p += 1
                p += (N_y - i - 1) * N_x

            elif xshp[i] in ['linear','lin','const']:
                if xshp[i] == 'const': xvals = np.full(N_x, scale[i][0])
                elif xshp[i] in ['linear','lin']: xvals = np.linspace(0., scale[i][0], N_x)

                p += i * N_x
                for j in range(N_x):
                    w[p] = xvals[j] + shift[i]
                    p += 1
                p += (N_y - i - 1) * N_x

        elif yshp[i] in ['linear','lin','const']:
            if yshp[i] == 'const': yvals = np.full(N_y, scale[i][1])
            elif yshp[i] in ['linear','lin']: yvals = np.linspace(0., scale[i][1], N_y)

            if xshp[i] == 'cvec':
                for j in range(N_y):
                    w[p + i] = yvals[j] + shift[i]
                    p += N_x

            elif xshp[i] in ['linear','lin','const']:
                if xshp[i] == 'const': xvals = np.full(N_x, scale[i][0])
                elif xshp[i] in ['linear','lin']: xvals = np.linspace(0., scale[i][0], N_x)

                for j in range(N_y):
                    for k in range(N_x):
                        w[p] = xvals[k] + yvals[j] + shift[i]
                        p += 1

    return w
