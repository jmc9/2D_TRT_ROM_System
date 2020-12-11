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

class DMD:
    def __init__(self,name,N_modes,eval,evec):
        self.name = name
        self.N_modes = N_modes
        self.eval = complex(eval[0],eval[1])
        self.evec = complex(evec[0],eval[1])

class DMD_opts:
    def __init__():
        self.names = names
        self.modes = modes
        self.plt_evecs = plt_evecs
        self.plt_evecs = plt_evals
        self.clen = clen
        self.rank = rank
