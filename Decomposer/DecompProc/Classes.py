#!/usr/bin/env python3
#==================================================================================================================================#
#
# Classes.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#

#==================================================================================================================================#
#
#==================================================================================================================================#
class complex:
    def __init__(self,real,imag):
        self.real = real
        self.imag = imag

#==================================================================================================================================#
#
#==================================================================================================================================#
class Grid:
    def __init__(self,bnds,crds):
        self.bnds = bnds
        self.crds = crds

#==================================================================================================================================#
#
#==================================================================================================================================#
class DMD:
    def __init__(self,N_modes=0,eval=[0,0],evec=[0,0]):
        self.N_modes = N_modes
        self.eval = complex(eval[0],eval[1])
        self.evec = complex(evec[0],evec[1])

#==================================================================================================================================#
#
#==================================================================================================================================#
class POD:
    U  = []
    S  = []
    Vt = []

#==================================================================================================================================#
#
#==================================================================================================================================#
class Data:
    def __init__(self,name):
        self.name  = name
        self.dat   = []
        self.grids = Grid([],[])
        self.dims  = []
        self.clen  = 0.
        self.rank  = 0.
        self.opt   = [0]

    def typeset(self,type):
        self.type = type

        if (type in ['DMD','DMDg']):
            self.dat = DMD
        elif (type in ['POD','PODg']):
            self.dat = POD
        else:
            print('Error! unsupported decomposition type detected')
            quit()