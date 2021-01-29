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
import numpy as np
import math

import Testing_Tools as tt
import ToolBox as tb
import Test_Functions as tf


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
    # Single exponential with one grid point in space                                                    #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_1', testdir)
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, test_drop, [.1], [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 2                                               #
    # multiple exponentials in 1D                                                                        #
    # each exponential is attached to a unit vector in space (i.e. w_1=[0 1 0 0 0])                      #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_2', testdir)

    #Test 2-a
    #all positive exponentials
    subtest_drop = tb.dirset('Test_2-a', test_drop)
    evals = [.1, .2, .3, .4, .5]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 2-b
    #all negative exponentials
    subtest_drop = tb.dirset('Test_2-b', test_drop)
    evals = [-.1, -.2, -.3, -.4, -.5]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 2-c
    #all negative exponentials
    subtest_drop = tb.dirset('Test_2-c', test_drop)
    evals = [.1, .2, .3, .4, .5, -.1, -.2, -.3, -.4, -.5]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 3                                               #
    # Multiple exponentials of the same value                                                            #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_3', testdir)

    #Test 3-a
    #5 eigenvalues attached to unit vectors in space
    subtest_drop = tb.dirset('Test_3-a', test_drop)
    evals = [.1, .1, .1, .1, .1]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 3-b
    #5 eigenvalues attached to constant vectors in space
    subtest_drop = tb.dirset('Test_3-b', test_drop)
    evals = [.1, .1, .1, .1, .1]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'const')

    #Test 3-c
    #5 eigenvalues attached to varying vectors in space
    subtest_drop = tb.dirset('Test_3-c', test_drop)
    evals = [.1, .1, .1, .1, .1]
    xshp = ['lin', 'cvec', 'const', 'cvec', 'lin']
    shift = [0., 0., 0., 0., -2.]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], [xshp,'none'], shift=shift)

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 4                                               #
    # Extreme eigenvalues                                                                                #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_4', testdir)

    #Test 4-a
    #single positive eigenvalue
    subtest_drop = tb.dirset('Test_4-a', test_drop)
    evals = [100.]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 4-b
    #single negative eigenvalue
    subtest_drop = tb.dirset('Test_4-b', test_drop)
    evals = [-100.]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 4-c
    #single positive tiny eigenvalue
    subtest_drop = tb.dirset('Test_4-c', test_drop)
    evals = [1e-10]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 4-d
    #single negative tiny eigenvalue
    subtest_drop = tb.dirset('Test_4-d', test_drop)
    evals = [-1e-10]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 4-e
    #combined extreme positive and negative eigenvalues
    subtest_drop = tb.dirset('Test_4-e', test_drop)
    evals = [100., -100.]
    tp = [0., .5, 1., 1.5]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [tp,[],[]], 'cvec', eps=-1.)

    #Test 4-f
    #combined tiny positive and negative eigenvalues
    subtest_drop = tb.dirset('Test_4-f', test_drop)
    evals = [1e-10, -1e-10]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec', eps=-1.)

    #Test 4-g
    #combined extreme positive eigenvalue with several reasonable values
    subtest_drop = tb.dirset('Test_4-g', test_drop)
    evals = [100., .1, .2, -.1, -.2]
    tp = np.linspace(0.,1.,10)
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [tp,[],[]], 'cvec', eps=-1.)

    #Test 4-h
    #combined extreme negative eigenvalue with several reasonable values
    subtest_drop = tb.dirset('Test_4-h', test_drop)
    evals = [-100., .1, .2, -.1, -.2]
    tp = np.linspace(0.,1.,10)
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [tp,[],[]], 'cvec')

    #Test 4-i
    #single positive tiny eigenvalue with several reasonable values
    subtest_drop = tb.dirset('Test_4-i', test_drop)
    evals = [1e-10, .1, .2, -.1, -.2]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 4-j
    #single negative tiny eigenvalue with several reasonable values
    subtest_drop = tb.dirset('Test_4-j', test_drop)
    evals = [-1e-10, .1, .2, -.1, -.2]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 5                                               #
    # Clustered eigenvalues                                                                              #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_5', testdir)

    #Test 5-a
    #single cluster of eigenvalues
    subtest_drop = tb.dirset('Test_5-a', test_drop)
    # evals = [.1, .1+1e-5, .1+2e-5]
    evals = [.1, .1+1e-4, .1+2e-4]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #Test 5-b
    #two clusters of eigenvalues
    subtest_drop = tb.dirset('Test_5-b', test_drop)
    evals = [.1, .1+1e-4, .1+2e-4, -.3, -.3-1e-4, -.3-2e-4]
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 6                                               #
    # Many distinct eigenvalues                                                                          #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_6', testdir)

    #Test 6-a
    subtest_drop = tb.dirset('Test_6-a', test_drop)
    evals = np.linspace(-5., 5., 50)
    tp = np.linspace(0., 5., 101)
    xp = np.linspace(0., 1., 50)
    tt.Run_Test_exp(exec_dir, decomp_perphs, proc_perphs, subtest_drop, evals, [tp, xp, []], 'cvec')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 7                                               #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_7', testdir)

    #Test 7-a
    subtest_drop = tb.dirset('Test_7-a', test_drop)
    tp = np.linspace(1., math.pi, 51)
    xp = np.linspace(-1., 1., 101)

    tt.Run_Test_nak(exec_dir, decomp_perphs, proc_perphs, subtest_drop, tp, xp)

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 8                                               #
    # A suite of tests that uses data from the 1D TRT code                                               #
    # both grey and multigroup QD factors in 1D are used                                                 #
    # several groups are done independently, followed by groupwise and full-phase space DMD              #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_8', testdir)

    #Grey QD factors
    subtest_drop = tb.dirset('Test_8-a', test_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subtest_drop, 'QD_factors.out', 'f')

    #Selecting single group at a time of multigroup QD factors
    subtest_drop = tb.dirset('Test_8-b', test_drop)
    #
    #Group-2 of the multigroup QD factors
    subsub_drop = tb.dirset('Test_8-b-1', subtest_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'MG_QD_factors.out', 'f_g02')
    #
    #Group-8 of the multigroup QD factors
    subsub_drop = tb.dirset('Test_8-b-2', subtest_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'MG_QD_factors.out', 'f_g08')
    #
    #Group-17 of the multigroup QD factors
    subsub_drop = tb.dirset('Test_8-b-3', subtest_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'MG_QD_factors.out', 'f_g17')

    #groupwise decomposition of multigroup QD factors
    subtest_drop = tb.dirset('Test_8-c', test_drop)
    #
    #full dataset
    subsub_drop = tb.dirset('Test_8-c-1', subtest_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'MG_QD_factors.out', 'f', 'yes', 'DMDg')
    #
    #only from 1-6 ns
    subsub_drop = tb.dirset('Test_8-c-2', subtest_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'MG_QD_factors_t1-6.out', 'f', 'yes', 'DMDg')

    #full phase-space decomposition of multigroup QD factors
    subtest_drop = tb.dirset('Test_8-d', test_drop)
    #
    #full dataset
    subsub_drop = tb.dirset('Test_8-d-1', subtest_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'MG_QD_factors.out', 'f', 'yes', 'DMD')
    #
    #only from 1-6 ns
    subsub_drop = tb.dirset('Test_8-d-2', subtest_drop)
    tt.Run_Test_1dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'MG_QD_factors_t1-6.out', 'f', 'yes', 'DMD')

    #----------------------------------------------------------------------------------------------------#
    #                                               Test 9                                               #
    # A suite of tests that uses data from the 2D TRT code                                               #
    # Specifically this test looks at cell-averaged f_xx                                                 #
    #----------------------------------------------------------------------------------------------------#
    test_drop = tb.dirset('Test_9', testdir)

    subtest_drop = tb.dirset('Test_9-a', test_drop)

    #Selecting single group at a time of multigroup QD factors
    subtest_drop = tb.dirset('Test_9-a', test_drop)
    #
    #Group-2 of the multigroup QD factors
    subsub_drop = tb.dirset('Test_9-a-1', subtest_drop)
    tt.Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'fom_2d.h5', 'fg_avg_xx', dcmp_type='DMD', group=2)
    #
    #Group-8 of the multigroup QD factors
    subsub_drop = tb.dirset('Test_9-a-2', subtest_drop)
    tt.Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'fom_2d.h5', 'fg_avg_xx', dcmp_type='DMD', group=8)
    #
    #Group-17 of the multigroup QD factors
    subsub_drop = tb.dirset('Test_9-a-3', subtest_drop)
    tt.Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'fom_2d.h5', 'fg_avg_xx', dcmp_type='DMD', group=17)

    #groupwise decomposition of multigroup QD factors
    subtest_drop = tb.dirset('Test_9-b', test_drop)
    #
    #full dataset
    subsub_drop = tb.dirset('Test_9-b-1', subtest_drop)
    tt.Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'fom_2d.h5', 'fg_avg_xx', dcmp_type='DMDg')
    #
    #only from 1-6 ns
    subsub_drop = tb.dirset('Test_9-b-2', subtest_drop)
    tt.Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'fom_2d.h5', 'fg_avg_xx', tstart=.1, tend=.6, dcmp_type='DMDg')

    #full phase-space decomposition of multigroup QD factors
    subtest_drop = tb.dirset('Test_9-c', test_drop)
    #
    #full dataset
    subsub_drop = tb.dirset('Test_9-c-1', subtest_drop)
    tt.Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'fom_2d.h5', 'fg_avg_xx', dcmp_type='DMD')
    #
    #only from 1-6 ns
    subsub_drop = tb.dirset('Test_9-c-2', subtest_drop)
    tt.Run_Test_2dtrt(exec_dir, decomp_perphs, proc_perphs, subsub_drop, 'fom_2d.h5', 'fg_avg_xx', tstart=.1, tend=.6, dcmp_type='DMD')

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
