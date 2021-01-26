#==================================================================================================================================#
#
# Dimension_Handling.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import numpy as np

#==================================================================================================================================#
#
#==================================================================================================================================#
class grid():

    # len = 0
    # pts = []
    # bnd = []

    def __init__(self, pts=[], bnd=[]):
        self.pts = pts
        self.len = len(pts)
        if bnd:
            self.bnd = bnds
        else:
            if self.len >= 2:
                self.bnd = [pts[0], pts[np-1]]
            else:
                self.bnd = []

#==================================================================================================================================#
#
#==================================================================================================================================#
class grids():

    # t = grid([],[])
    # g = grid([],[])
    # y = grid([],[])
    # x = grid([],[])
    #
    # xy  = grid()
    # all = grid()

    def __init__(self):
        self.t = grid()
        self.g = grid()
        self.y = grid()
        self.x = grid()
        self.xy = grid()
        all = grid()
        return

    def get_bnd(self, name):
        if   name == 't': return self.t.bnd
        elif name == 'y': return self.y.bnd
        elif name == 'x': return self.x.bnd
        else: return 0

    def get_len(self, name):
        if   name == 't': return self.t.len
        elif name == 'g': return self.g.len
        elif name == 'y': return self.y.len
        elif name == 'x': return self.x.len
        else: return 0

    def get_pts(self, name):
        if   name == 't': return self.t.pts
        elif name == 'y': return self.y.pts
        elif name == 'x': return self.x.pts
        else: return 0

    #returns list of all basic dimension lengths
    def dlist(self):
        return [self.t.len, self.y.len, self.x.len]

    #returns list of all basic dimension names (in the same order as given by dlist)
    def nlist(self):
        return ['t', 'y', 'x']

    #returns list of grid points all basic grids
    def glist(self):
        return [self.t.pts, self.y.pts, self.x.pts]

    #same as dlist, but the returned list does not include any 0-length dimensions
    def dlist_tr(self):
        dl = []
        for d in self.dlist():
            if d != 0: dl.append(d)
        return dl

    #same as nlist, but the returned list does not include the names of any 0-length dimensions
    def nlist_tr(self):
        nl = []
        for (d, n) in zip(self.dlist(), self.nlist()):
            if d != 0: nl.append(n)
        return nl

    #same as glist, but the returned list does not include any 0-length grids
    def glist_tr(self):
        gl = []
        for (d, g) in zip(self.dlist(), self.glist()):
            if d != 0: gl.append(g)
        return gl

    #calculates the total dimensionality of a list of dimension lengths
    def dcalc(self, dlist):
        l = 1
        c = False
        for d in dlist:
            if d != 0: l = l * d
            if d == 1: c = True
        if ((l == 1)and(not c)):
            return 0
        else:
            return l

    #upon call, updates all non-basic dimension lengths
    def recalc(self):
        self.all.len = self.dcalc(self.dlist())
        self.xy.len = self.dcalc([self.y.len, self.x.len])

#==================================================================================================================================#
#
#==================================================================================================================================#
def grid_calc(n=10., gp=[], opt=''):
    g = grid()
    if opt == 'none':
        g.len = 0
        g.bnd = [0., 0.]
        g.pts = []

    else:
        if len(gp) == 0: #no grid input -> make simple grid based on number of expansion terms
            g.len = n
            g.bnd = [0., g.len*1. - 1.]
            g.pts = np.linspace(g.bnd[0], g.bnd[1], g.len)
        else: #grid input found -> collect dimension and bounds
            g.len = len(gp)
            g.bnd = [gp[0], gp[g.len - 1]]
            g.pts = gp

    return g
