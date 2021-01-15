#==================================================================================================================================#
#
# DecompSrc_Handling.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#generic python tools
import os
import shutil

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_decompsrc(decomp_perphs, f_file, drop, exec_dir):

    #make all preparations for execution of DecompSrc
    exec_prep(decomp_perphs, f_file)

    #execute DecompSrc
    call_decompsrc(decomp_perphs, exec_dir)

    #clean and organize directories post DecompSrc execution/termination
    exec_clean(decomp_perphs, f_file, drop)

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def write_inpf(decomp_perphs, f_file, f_name='f', dcmp_type='DMD', svd_eps=1e-14):
    file = open(decomp_perphs.inp, 'w')
    file.write('dataset {}\n'.format(f_file))
    file.write('outfile {}\n'.format(decomp_perphs.out))
    file.write('dcmp_data 1 {}\n'.format(f_name))
    file.write('dcmp_type {}\n'.format(dcmp_type))
    file.write('svd_eps {}\n'.format(svd_eps))
    file.write('Prb_spec 0\n')
    file.write('Disc_wts 0\n')
    file.close()

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_prep(decomp_perphs, f_file):
    #create input file for DecompSrc
    write_inpf(decomp_perphs, f_file)

    #move required files for DecompSrc operation into DecompSrc directory from Testing directory
    shutil.move(decomp_perphs.inp, os.path.join(decomp_perphs.path, decomp_perphs.inp))
    shutil.move(f_file, os.path.join(decomp_perphs.path, f_file))

#==================================================================================================================================#
#
#==================================================================================================================================#
def call_decompsrc(decomp_perphs, exec_dir):
    os.chdir(decomp_perphs.path) #move into DecompSrc directory
    os.system('./decomposer.exe > {}'.format(decomp_perphs.log)) #excecute DecompSrc
    os.chdir(exec_dir) #move back into Testing directory

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_clean(decomp_perphs, f_file, drop):
    #move DecomSrc input file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_perphs.path, decomp_perphs.inp), os.path.join(drop, decomp_perphs.inp) )

    #move test data file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_perphs.path, f_file), os.path.join(drop, f_file) )

    #move DecomSrc output file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_perphs.path, decomp_perphs.out), os.path.join(drop, decomp_perphs.out) )

    #copy DecomSrc log file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_perphs.path, decomp_perphs.log), os.path.join(drop, decomp_perphs.log) )
