# 2D_TRT_ROM_System

## SUMMARY
The 2D_TRT_ROM_System is designed to solve thermal radiative transfer (TRT) problems in 2D geometry with data-driven reduced-order models. 
There are currently four distinct codes that comprise the overall system:
  1) The 2D TRT ROM code, contained in the MainSrc directory
  2) The TRT_processor code, contained in the Processing directory
  3) The Decomposer code, contained in the Decomposer/DecompSrc directory
  4) The Decomp_Processor code, contained in the Decomposer/DecompProc directory

## EXTERNAL DEPENDENCIES
Several external packages are required to compile and run the code system, including:
  1) GNU compilers (gfortran & gcc) <br/>
  Here's the documentation via [wiki](https://gcc.gnu.org/), and [GNU website](https://gnu.org) <br/>
  Installation instructions:
      - On Linux based systems, the following commands will generally suffice:  <br/>
        `sudo apt install gcc gfortran`
      - The gcc compiler as well as other useful tools like make are included in packages like: <br/>
        `sudo apt install build-essential`
    
  2) Python3 & pip <br/>
  Installation instructions:
      - On Linux based systems, the following commands will generally suffice: <br/>
        `sudo apt install python3` <br/>
        `sudo apt install python3-pip`
  
  3) LAPACK <br/>
  Documentation is found on [Netlib's site](http://www.netlib.org/lapack/) <br/>
  Installation instructions:
      - On Linux based systems, the following commands will generally suffice: <br/>
        `sudo apt-get install libblas-dev liblapack-dev`
    
  4) LAPACKE <br/>
  Documentation is found on [Netlib's site](http://www.netlib.org/lapack/lapacke.html) <br/>
  Installation instructions:
      - On Linux based systems that already have LAPACK installed, the following commands will generally suffice: <br/>
        `sudo apt-get install liblapacke-dev`
    
  5) SPARSKIT <br/>
  Documentation and downloads are on [Yousef Saad's site](https://www-users.cs.umn.edu/~saad/software/SPARSKIT/index.html) <br/>
  Installation instructions:
      - Instructions are included in the readme of the download from the above link
      - Once the tar.gz file has been downloaded and unpacked, makefiles are present that use the gfortran compiler
  
  6) NetCDF <br/>
  Documentation is found on [unidata's site](https://www.unidata.ucar.edu/software/netcdf/) <br/>
  Installation instructions:
      - On Linux based systems that already have LAPACK installed, the following commands will generally suffice: <br/>
        `sudo apt-get install libnetcdf-dev libnetcdff-dev`
      - On Linux based systems with pip3 installed, the following commands will generally suffice for the netcdf4 API: <br/>
        `pip install netCDF4`
  
## CODE DESCRIPTIONS
### 1) 2D TRT ROM Code
Directory: MainSrc <br/>
Language: FORTRAN <br/>
Compilation: Makefile included <br/>
Basic Description: Solves TRT problems

### 2) TRT Processor
Directory: Processing <br/>
Language: Python <br/>
Basic Description: Calculates errors in reduced-order solutions compared to the full-order solution and generates plots of the solution(s) and error(s)

### 3) Decomposer
Directory: MainSrc <br/>
Language: C <br/>
Compilation: Makefile included <br/>
Basic Description: Applies decomposition techniques to data

### 4) Decomp Processor
Directory: MainSrc <br/>
Language: Python <br/>
Basic Description: Visualizes the decompositions produced by the Decomposer code
 
  
  
  
  
