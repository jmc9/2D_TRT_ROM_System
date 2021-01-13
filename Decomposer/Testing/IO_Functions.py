#==================================================================================================================================#
#
# IO_Functions.py
# Author: Joseph M. Coale
# jmcoale@ncsu.edu - josephcoale912@gmail.com
#
#==================================================================================================================================#
#importing packages
#==================================================================================================================================#
#netcdf tools
from netCDF4 import Dataset as ncds

#==================================================================================================================================#
#
#==================================================================================================================================#
def output_f(f, grids, outfile):

    dset = ncds(outfile,'w') #create ncdf datafile/dataset

    #create dimensions and grids in datafile
    names = grids.nlist_tr() #get names of non-zero dimensions
    dnames = ['N_{}'.format(n) for n in names] #format dimension names
    gnames = ['{}p'.format(n) for n in names]  #format grid names
    nc_dims = []
    nc_grids = []
    for (n, dn, gn) in zip(names, dnames, gnames):
        nc_dims.append( dset.createDimension( dn, grids.get_len(n) ) ) #create dimension named 'dn' with length = grids.get_len(n)
        nc_grids.append( dset.createVariable( gn, 'f8', dn ) ) #create grid variable named 'gn' with dimension 'dn'

    #write grid bounds and grid data to datafile
    for (g, n) in zip(nc_grids, names):
        g.bnds = grids.get_bnd(n)
        g[:] = grids.get_pts(n)

    #create data (f) variable in datafile, write in values
    nc_f = dset.createVariable('f', 'f8', tuple(dnames) )
    nc_f[:] = f

    #write number of grids and each grid name as attributes of 'f' in the datafile
    nc_f.N_grids = len(gnames)
    for i in range(len(gnames)):
        nc_f.temp = gnames[i]
        nc_f.renameAttribute('temp', 'grid{}'.format(i))

    dset.close() #close dataset

    return
