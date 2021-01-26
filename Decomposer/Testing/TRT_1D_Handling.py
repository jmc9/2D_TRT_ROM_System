#==================================================================================================================================#
#
# 1D_TRT_Handling.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import numpy as np
import math

from ToolBox import File_Find as ffind

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Block(file, data_name):
    block_found = False
    for line in file:
        if (line.strip()): #disregarding blank lines
            input = line.split() #delimiting input line by space

            if len(input)>2:
                if input[2] == '{}-01'.format(data_name):
                    block_found = True
                    break

    if not block_found: return

    f = []
    for line in file:
        if (line.strip()): #disregarding blank lines
            input = line.split() #delimiting input line by space

            f.append(input[2:])

        else:
            file.seek(0)
            break

    f = np.array(f)
    f = f.astype(np.float)
    f = np.transpose(f)

    return f

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_Data(file, data_name, N_g, mg_opt):

    if N_g == 0:
        f = Read_Block(file, data_name)

    else:
        N_g = 3
        fg = Read_Block(file, '{}_g01'.format(data_name))
        (N_t, N_x) = np.shape(fg)

        f = np.zeros((N_t, N_g, N_x))
        for g in range(N_g):
            fg = Read_Block(file, '{}_g{:02d}'.format(data_name, g+1))

            for t in range(N_t):
                for x in range(N_x):
                    f[t][g][x] = fg[t][x]

    return f

#==================================================================================================================================#
#
#==================================================================================================================================#
def Read_TRT(file_name, data_name, mg_opt=''):
    file = open(file_name, 'r')

    N_t = int( ffind(file, 'time_steps') )
    tlen = float( ffind(file, 'tlen') )
    N_x = int( ffind(file, 'cells') )
    xlen = float( ffind(file, 'xlen') )

    if mg_opt == '':
        N_g = 0
    else:
        N_g = int( ffind(file, 'groups') )

    f = Read_Data(file, data_name, N_g, mg_opt)

    file.close()

    return f, N_t, N_x, tlen, xlen
