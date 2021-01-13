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
def exec_decompsrc(decomp_infile, decomp_outfile, f_file, decomp_logfile, drop, exec_dir, decomp_dir):

    #make all preparations for execution of DecompSrc
    exec_prep(decomp_infile, decomp_outfile, f_file, decomp_dir)

    #execute DecompSrc
    call_decompsrc(decomp_logfile, decomp_dir, exec_dir)

    #clean and organize directories post DecompSrc execution/termination
    exec_clean(decomp_infile, decomp_outfile, f_file, decomp_logfile, decomp_dir, drop)

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def write_inpf(decomp_infile, decomp_outfile, f_file, f_name='f', dcmp_type='DMD', svd_eps=-1.):
    file = open(decomp_infile, 'w')
    file.write('dataset {}\n'.format(f_file))
    file.write('outfile {}\n'.format(decomp_outfile))
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
def exec_prep(decomp_infile, decomp_outfile, f_file, decomp_dir):
    #create input file for DecompSrc
    write_inpf(decomp_infile, decomp_outfile, f_file)

    #move required files for DecompSrc operation into DecompSrc directory from Testing directory
    shutil.move(decomp_infile, os.path.join(decomp_dir, decomp_infile))
    shutil.move(f_file, os.path.join(decomp_dir, f_file))

#==================================================================================================================================#
#
#==================================================================================================================================#
def call_decompsrc(decomp_logfile, decomp_dir, exec_dir):
    os.chdir(decomp_dir) #move into DecompSrc directory
    os.system('./decomposer.exe > {}'.format(decomp_logfile)) #excecute DecompSrc
    os.chdir(exec_dir) #move back into Testing directory

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_clean(decomp_infile, decomp_outfile, f_file, decomp_logfile, decomp_dir, drop):
    #move DecomSrc input file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_dir, decomp_infile), os.path.join(drop, decomp_infile) )

    #move test data file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_dir, f_file), os.path.join(drop, f_file) )

    #move DecomSrc output file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_dir, decomp_outfile), os.path.join(drop, decomp_outfile) )

    #copy DecomSrc log file to test directory from DecomSrc directory
    shutil.move( os.path.join(decomp_dir, decomp_logfile), os.path.join(drop, decomp_logfile) )
