#!/usr/bin/env python3
#==================================================================================================================================#
#
# Decomp_Tester.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import getopt, sys #handles user inputs from command line
import os

import Testing_Tools as tt
import ToolBox as tb

#==================================================================================================================================#
#
#==================================================================================================================================#
def Exec():
    exec_dir, decomp_dir, proc_dir = Path_Gen()
    Options()

    Decomp_Tester(exec_dir, decomp_dir, proc_dir)

#==================================================================================================================================#
#
#==================================================================================================================================#
def Decomp_Tester(exec_dir, decomp_dir, proc_dir, testdir='testdir'):

    testdir = tb.dirset(testdir, exec_dir)

    decomp_perphs = tt.peripherals(decomp_dir)
    proc_perphs = tt.peripherals(proc_dir)

    decomp_perphs.load('input.inp', 'test_dmd.h5', 'DecompLog.log')
    proc_perphs.load('proc_input.inp', 'processed_test', 'ProcLog.log')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 1                                               #
    # Simple test on basic calculations                                                                  #
    # single small exponential growth term in 1D with single grid point in space                         #
    #----------------------------------------------------------------------------------------------------#
    test1_drop = tb.dirset('Test_1', testdir)
    tt.Run_Test(exec_dir, decomp_perphs, proc_perphs, test1_drop, [.1], [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 2                                               #
    # multiple (5) exponential growth terms in 1D                                                        #
    # each exponential is attached to a unit vector in space (i.e. w_1=[0 1 0 0 0])                      #
    #----------------------------------------------------------------------------------------------------#
    test2_drop = tb.dirset('Test_2', testdir)
    tt.Run_Test(exec_dir, decomp_perphs, proc_perphs, test2_drop, [.1, .2, .3, .4, .5], [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 3                                               #
    # an exact copy of Test 2, but with exponential decay instead og growth                              #
    #----------------------------------------------------------------------------------------------------#
    test3_drop = tb.dirset('Test_3', testdir)
    tt.Run_Test(exec_dir, decomp_perphs, proc_perphs, test3_drop, [-.1, -.2, -.3, -.4, -.5], [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 4                                               #
    # now combines exponential growth and decay in 1D                                                    #
    # each exponential is attached to a unit vector in space (i.e. w_1=[0 1 0 0])                        #
    #----------------------------------------------------------------------------------------------------#
    test4_drop = tb.dirset('Test_4', testdir)
    tt.Run_Test(exec_dir, decomp_perphs, proc_perphs, test4_drop, [.1, -.1, .2, -.2], [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 5                                               #
    # The same test as Test 4 but in 2D geometry                                                         #
    # instead of 4 grid points along the x-axis, --                                                      #
    # -- there are 2 grid points along each of the x- and y- axes                                        #
    # each exponential is attached to a unit vector in 2D space (i.e. w_(1,0)=[ [0 1] [0 0] ])           #
    #----------------------------------------------------------------------------------------------------#
    test5_drop = tb.dirset('Test_5', testdir)
    tt.Run_Test(exec_dir, decomp_perphs, proc_perphs, test5_drop, [.1, -.1, .2, -.2], [[],[0.,1.],[0.,1.]], ['cvec', 'cvec'])

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 6                                               #
    # 2D geometry test with eigenvectors that are linear in the y-direction                              #
    # The x- and y- dimensions are now also different (non-square domain)                                #
    # There exist more grid points in space then there are eigenvalues                                   #
    #----------------------------------------------------------------------------------------------------#
    test6_drop = tb.dirset('Test_6', testdir)
    tt.Run_Test(exec_dir, decomp_perphs, proc_perphs, test6_drop, [.1, .3, -.2, -.4], [[],[],[0.,1.,2.,3.,4.,5.,6.,7.,8.,9.]], ['cvec', 'lin'])

#==================================================================================================================================#
#
#==================================================================================================================================#
def Path_Gen():
    decomp_dir = '../DecompSrc'
    proc_dir = '../DecompProc'

    exec_dir = os.getcwd()
    decomp_dir = os.path.join(exec_dir, decomp_dir)
    proc_dir = os.path.join(exec_dir, proc_dir)

    return exec_dir, decomp_dir, proc_dir

#==================================================================================================================================#
#
#==================================================================================================================================#
def Options():
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'ho:', ['help','output='])
    except getopt.GetoptError as err:
        print(str(err))
        usage()
        quit()

    Default_Opts()

    for o, a in opts:
        if o in ['-h','--help']:
            usage()
            quit()

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def Default_Opts():

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def usage():
    print()
    print('Options:')
    print('-h [help]: prints this message')
    print()

#==================================================================================================================================#
# MAIN PROGRAM
# This portion of the code is run when this file is the main program
# i.e. when the call is made './Decomp_Tester.py'
#==================================================================================================================================#
if __name__ == '__main__':
    Exec()
