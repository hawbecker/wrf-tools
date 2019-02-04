'''
  What it does: This dictionary contains functions for reading
                WRF (or netCDF) data.

  Who made it: patrick.hawbecker@nrel.gov
  When: 5/11/18
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm    
from netCDF4 import Dataset as ncdf
import pickle
import subprocess


def getdims(wrff):
    try:
        nx = wrff.dimensions['west_east'].size
    except KeyError:
        print 'No x-dimension'; nx = []
    try:
        ny = wrff.dimensions['south_north'].size
    except KeyError:
        print 'No y-dimension'; ny = []
    try:
        nz = wrff.dimensions['bottom_top'].size
    except KeyError:
        print 'No z-dimension'; nz = []
    try:
        nt = wrff.dimensions['Time'].size
    except KeyError:
        print 'No t-dimension'; nt = []
    return nt,nz,ny,nx

def getavgheight(wrff):
    nt = wrff.dimensions['Time'].size
    try:
        nz = wrff.dimensions['bottom_top'].size
    except KeyError:
        print 'No z-dimension'; return []
    if nt == 1:
        ph  = wrff.variables['PH'][0,:,:,:]
        phb = wrff.variables['PHB'][0,:,:,:]
        hgt = wrff.variables['HGT'][0,:,:]
        z   = np.mean(np.mean(((ph+phb)/9.81) - hgt,axis=1),axis=1)
        zs  = (z[1:] + z[:-1])*0.5
    else:
        z = np.zeros((nt,nz+1))
        for tt in range(0,nt):
            ph  = wrff.variables['PH'][tt,:,:,:]
            phb = wrff.variables['PHB'][tt,:,:,:]
            hgt = wrff.variables['HGT'][tt,:,:]
            z[tt,:] = np.mean(np.mean(((ph+phb)/9.81) - hgt,axis=1),axis=1)
        zs  = (z[:,1:] + z[:,:-1])*0.5
    return z,zs

def getheight(wrff):
    try:
        nz = wrff.dimensions['bottom_top'].size
    except KeyError:
        print 'No z-dimension'; return []
    ph  = wrff.variables['PH'][0,:,:,:]
    phb = wrff.variables['PHB'][0,:,:,:]
    hgt = wrff.variables['HGT'][0,:,:]

    z   = ((ph+phb)/9.81) - hgt
    zs  = (z[1:,:,:] + z[:-1,:,:])*0.5
    return z,zs

def getheightloc(wrff,y,x):
    nt = wrff.dimensions['Time'].size
    try:
        nz = wrff.dimensions['bottom_top'].size
    except KeyError:
        print 'No z-dimension'; return []
    if nt == 1:
        ph  = wrff.variables['PH'][0,:,y,x]
        phb = wrff.variables['PHB'][0,:,y,x]
        hgt = wrff.variables['HGT'][0,y,x]
        z   = ((ph+phb)/9.81) - hgt
        zs  = (z[1:] + z[:-1])*0.5
    else:
        z = np.zeros((nt,nz+1))
        for tt in range(0,nt):
            ph  = wrff.variables['PH'][tt,:,y,x]
            phb = wrff.variables['PHB'][tt,:,y,x]
            hgt = wrff.variables['HGT'][tt,y,x]
            z[tt,:] = ((ph+phb)/9.81) - hgt
        zs  = (z[:,1:] + z[:,:-1])*0.5
    return z,zs

def getwrffiles(fdir,fstr,returnFileNames=True):
    #fdir = file directory; fstr = file string structure (e.g. 'wrfout')
    nwrffs = subprocess.check_output('cd %s && ls %s*' % (fdir,fstr), shell=True).split()
    nt = np.shape(nwrffs)[0]
    if returnFileNames==True:
        return nwrffs,nt
    else:
        return nt

def latlon(wrff):
    lat = wrff.variables['XLAT'][0,:,:]
    lon = wrff.variables['XLONG'][0,:,:]
    return lat,lon

#=====================================#
# - - - - - - TOWER STUFF - - - - - - #
def gettowers(fdir,tstr):
    f = open('%s%s' % (fdir,tstr))
    nt = sum(1 for line in f)-3; f.close()
    f = open('%s%s' % (fdir,tstr))
    f.readline(); f.readline(); f.readline()
    tname = []; tij = np.zeros((2,nt))
    for tt in range(0,nt):
        line = f.readline().split()
        tname.append(line[1])
        tij[0,tt] = line[2]; tij[1,tt] = line[3]
    return tname,tij

def twrlocij(twrf):
    twr = open(twrf,'r')
    header = twr.readline().replace('(',' ').replace(')',' ').replace(',',' ').split()
    twr.close()
    stni = int(header[6]) - 1
    stnj = int(header[7]) - 1
    return stni,stnj

def twrlocll(twrf):
    twr = open(twrf,'r')
    header = twr.readline().replace('(',' ').replace(')',' ').replace(',',' ').split()
    twr.close()
    stni = float(header[9])
    stnj = float(header[8])
    return stnj,stni

class tower():
    def __init__(self,fstr):
        self.fstr = fstr
        self.getvars()
        self.getdata()
    def getvars(self):
        varns = subprocess.check_output('ls %s*' % (self.fstr), shell=True).split() 
        nvars = np.shape(varns)[0]
        for vv in range(0,nvars):
            varns[vv] = varns[vv].replace(self.fstr,'')
        self.varns = varns
        self.nvars = nvars
    def getdata(self):
        for vv in range(0,self.nvars):
            if self.varns[vv] != 'TS':
                f = open('%s%s' % (self.fstr,self.varns[vv]))
                nt = sum(1 for line in f)-1; f.close()
                f = open('%s%s' % (self.fstr,self.varns[vv]))
                self.header = f.readline().split()
                for tt in np.arange(0,nt):
                    line = f.readline().split()
                    if tt == 0: 
                        nz = np.shape(line)[0]-1
                        var = np.zeros((nt,nz))
                        ttime = np.zeros((nt))
                    var[tt,:] = line[1:]
                    ttime[tt] = np.float(line[0])
                self.nt   = nt
                self.time = ttime
                self.nz = nz
                if self.varns[vv] == 'PH':
                    self.ph = var
                elif self.varns[vv] == 'QV':
                    self.qv = var
                elif self.varns[vv] == 'TH':
                    self.th = var
                elif self.varns[vv] == 'UU':
                    self.uu = var
                elif self.varns[vv] == 'VV':
                    self.vv = var
                elif self.varns[vv] == 'WW':
                    self.ww = var
            elif self.varns[vv] == 'TS':
                f = open('%s%s' % (self.fstr,self.varns[vv]))
                nt = sum(1 for line in f)-1; f.close()
                f = open('%s%s' % (self.fstr,self.varns[vv]))
                f.readline()
                for tt in np.arange(0,nt):
                    line = f.readline().split()
                    if tt == 0: 
                        nv = np.shape(line)[0]-2
                        var = np.zeros((nt,nv))
                    var[tt,:] = line[2:]
                self.ts = var

#=====================================#

def wrftimes2hours(wrff):
    nt = np.shape(wrff.variables['Times'][:])[0]
    if nt == 1:
        time = ''.join(wrff.variables['Times'][0])
        year = np.float(time[:4]);    month = np.float(time[5:7])
        day  = np.float(time[8:10]);  hour  = np.float(time[11:13])
        minu = np.float(time[14:16]); sec   = np.float(time[17:19])
        hours = hour + minu/60.0 + sec/(60.0*60.0)
    else:
        year = np.asarray([]); month = np.asarray([])
        day  = np.asarray([]); hour  = np.asarray([])
        minu = np.asarray([]); sec   = np.asarray([])
        for tt in np.arange(0,nt):
            time = ''.join(wrff.variables['Times'][tt])
            year  = np.append(year,np.float(time[:4]))
            month = np.append(month,np.float(time[5:7]))
            day   = np.append(day,np.float(time[8:10]))
            hour  = np.append(hour,np.float(time[11:13]))
            minu  = np.append(minu,np.float(time[14:16]))
            sec   = np.append(sec,np.float(time[17:19]))
        hours = hour + minu/60.0 + sec/(60.0*60.0)
    return [year,month,day,hours]

def latlon2ij(wrff,latoi,lonoi):
    lat,lon = latlon(wrff)
    dist    = ((lat-latoi)**2 + (lon-lonoi)**2)**0.5
    jj,ii   = np.where(dist==np.min(dist))
    return ii[0],jj[0]

def unstagger2d(var,ax):
    if ax == 0:
        varu = (var[:-1,:] + var[1:,:])/2.0
    if ax == 1:
        varu = (var[:,:-1] + var[:,1:])/2.0
    return varu

def unstagger3d(var,ax):
    if ax == 0:
        varu = (var[:-1,:,:] + var[1:,:,:])/2.0
    if ax == 1:
        varu = (var[:,:-1,:] + var[:,1:,:])/2.0
    if ax == 2:
        varu = (var[:,:,:-1] + var[:,:,1:])/2.0
    return varu

def unstagger4d(var,ax):
    if ax == 1:
        varu = (var[:,:-1,:,:] + var[:,1:,:,:])/2.0
    if ax == 2:
        varu = (var[:,:,:-1,:] + var[:,:,1:,:])/2.0
    if ax == 3:
        varu = (var[:,:,:,:-1] + var[:,:,:,1:])/2.0
    return varu

