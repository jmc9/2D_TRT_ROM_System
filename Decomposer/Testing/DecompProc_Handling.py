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
def exec_decompproc(proc_infile, proc_outdir, decomp_outfile, f_file, proc_logfile, drop, exec_dir, proc_dir):

    #make all preparations for execution of DecompProc
    exec_prep(proc_infile, decomp_outfile, f_file, proc_dir, drop)

    #execute DecompProc
    call_decompproc(proc_infile, proc_outdir, proc_logfile, proc_dir, exec_dir)

    #clean and organize directories post DecompProc execution/termination
    exec_clean(proc_infile, proc_outdir, decomp_outfile, f_file, proc_logfile, proc_dir, drop)

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def write_inpf(proc_infile, decomp_outfile, f_file, f_name='f'):
    file = open(proc_infile, 'w')
    file.write('dset {}\n'.format(decomp_outfile))
    file.write('dset_dat {}\n'.format(f_file))
    file.write('dcmp_data {}\n'.format(f_name))
    file.close()

    return

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_prep(proc_infile, decomp_outfile, f_file, proc_dir, drop):
    #create input file for DecompSrc
    write_inpf(proc_infile, decomp_outfile, f_file)

    #copy required files for DecompSrc operation into DecompSrc directory
    shutil.move( proc_infile, os.path.join(proc_dir, proc_infile) )
    shutil.copyfile( os.path.join(drop, decomp_outfile), os.path.join(proc_dir, decomp_outfile) )
    shutil.copyfile( os.path.join(drop, f_file), os.path.join(proc_dir, f_file) )

#==================================================================================================================================#
#
#==================================================================================================================================#
def call_decompproc(proc_infile, proc_outdir, proc_logfile, proc_dir, exec_dir):
    os.chdir(proc_dir) #move into DecompProc directory
    os.system('./Decomp_Processer.py -i {} -o {} > {}'.format(proc_infile, proc_outdir, proc_logfile)) #excecute DecompProc
    os.chdir(exec_dir) #move back into Testing directory

#==================================================================================================================================#
#
#==================================================================================================================================#
def exec_clean(proc_infile, proc_outdir, decomp_outfile, f_file, proc_logfile, proc_dir, drop):
    #move DecomProc input file to test directory from DecomProc directory
    shutil.move( os.path.join(proc_dir,proc_infile), os.path.join(drop, proc_infile) )

    #move DecomProc log file to test directory from DecomProc directory
    shutil.move( os.path.join(proc_dir, proc_logfile), os.path.join(drop, proc_logfile) )

    #move DecomProc output sub-directory to test directory from DecomProc directory
    shutil.move( os.path.join(proc_dir, proc_outdir), os.path.join(drop, proc_outdir) )

    #delete test data file from DecomProc directory
    os.remove( os.path.join(proc_dir, f_file) )

    #delete DecomSrc output file from DecomProc directory
    os.remove( os.path.join(proc_dir, decomp_outfile) )
