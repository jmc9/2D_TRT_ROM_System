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
def exec_decompproc(proc_perphs, decomp_outfile, f_file, drop, exec_dir, trunc_eps=1., trunc_opt=2):

    #make all preparations for execution of DecompProc
    exec_prep(proc_perphs, decomp_outfile, f_file, drop, trunc_eps, trunc_opt)

    #execute DecompProc
    call_decompproc(proc_perphs, exec_dir)

    #clean and organize directories post DecompProc execution/termination
    exec_clean(proc_perphs, decomp_outfile, f_file, drop)

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def write_inpf(proc_perphs, decomp_outfile, f_file, f_name='f', trunc_eps=1., trunc_opt=2):
    file = open(proc_perphs.inp, 'w')
    file.write('dset {}\n'.format(decomp_outfile))
    file.write('dset_dat {}\n'.format(f_file))
    file.write('dcmp_data {}\n'.format(f_name))
    file.write('dcmp_data {}\n'.format(f_name))
    file.write('trunc_eps {}\n'.format(trunc_eps))
    file.write('trunc_opt {}\n'.format(trunc_opt))
    file.close()

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_prep(proc_perphs, decomp_outfile, f_file, drop, trunc_eps=1., trunc_opt=2):
    #create input file for DecompSrc
    write_inpf(proc_perphs, decomp_outfile, f_file, trunc_eps=trunc_eps, trunc_opt=trunc_opt)

    #copy required files for DecompSrc operation into DecompSrc directory
    shutil.move( proc_perphs.inp, os.path.join(proc_perphs.path, proc_perphs.inp) )
    shutil.copyfile( os.path.join(drop, decomp_outfile), os.path.join(proc_perphs.path, decomp_outfile) )
    shutil.copyfile( os.path.join(drop, f_file), os.path.join(proc_perphs.path, f_file) )

#==================================================================================================================================#
#
#==================================================================================================================================#
def call_decompproc(proc_perphs, exec_dir):
    os.chdir(proc_perphs.path) #move into DecompProc directory
    os.system('./Decomp_Processer.py -i {} -o {} > {}'.format(proc_perphs.inp, proc_perphs.out, proc_perphs.log)) #excecute DecompProc
    os.chdir(exec_dir) #move back into Testing directory

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_clean(proc_perphs, decomp_outfile, f_file, drop):
    #move DecomProc input file to test directory from DecomProc directory
    shutil.move( os.path.join(proc_perphs.path, proc_perphs.inp), os.path.join(drop, proc_perphs.inp) )

    #move DecomProc log file to test directory from DecomProc directory
    shutil.move( os.path.join(proc_perphs.path, proc_perphs.log), os.path.join(drop, proc_perphs.log) )

    #move DecomProc output sub-directory to test directory from DecomProc directory
    shutil.move( os.path.join(proc_perphs.path, proc_perphs.out), os.path.join(drop, proc_perphs.out) )

    #delete test data file from DecomProc directory
    os.remove( os.path.join(proc_perphs.path, f_file) )

    #delete DecomSrc output file from DecomProc directory
    os.remove( os.path.join(proc_perphs.path, decomp_outfile) )
