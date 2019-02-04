'''
  What it does: Plots time-series data from WRF's tslist output 
		

  Who made it: patrick.hawbecker@nrel.gov 
  When: 4/20/18
'''
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as ncdf
from matplotlib import cm
import subprocess
from matplotlib.colors import Normalize as Normalize
import wrfdict as wrfdict

# - - - - - - - - - - - - #
stnlat = 45.638004;            stnlon = -120.642973
fdir = '/projects/wfip/WFIP2/improved/extracted/'
savedir = '/projects/wfip/'
dom  = 1
# - - - - - - - - - - - - #

wrfoutf = subprocess.check_output('ls %swrfout_d0%d_*00' % (fdir,dom),shell=True).split()
wrfoutf = wrfoutf
nt = np.shape(wrfoutf)[0]

# Initialize time-series variables
hfx   = np.zeros((nt))
pblh  = np.zeros((nt))
ustr  = np.zeros((nt))
u10   = np.zeros((nt))
v10   = np.zeros((nt))
T2    = np.zeros((nt))
TH2   = np.zeros((nt))
swdwn = np.zeros((nt))
psfc  = np.zeros((nt))
time  = np.zeros((nt))
wdate = np.zeros((nt))

# - - - WRF Data - - -
cc = 0
for ff in wrfoutf:
    print ff
    wrfout     = ncdf(ff)
    year,month,day,hour = wrfdict.wrftimes2hours(wrfout)
    wdate[cc] = year*10000 + month*100 + day
    time[cc]  = hour
    if cc == 0:
        poii, poij = wrfdict.latlon2ij(wrfout,stnlat,stnlon)
        z,zs = wrfdict.getheight(wrfout)
        nz = np.shape(zs)[0]
        height = zs[:,poij,poii]
        u = np.zeros((nt,nz))
        v = np.zeros((nt,nz))
        w = np.zeros((nt,nz))
        T = np.zeros((nt,nz))
        P = np.zeros((nt,nz))
        Q = np.zeros((nt,nz))
    hfx[cc]    = wrfout.variables['HFX'][0,poij,poii]
    pblh[cc]   = wrfout.variables['PBLH'][0,poij,poii]
    psfc[cc]   = wrfout.variables['PSFC'][0,poij,poii]
    ustr[cc]   = wrfout.variables['UST'][0,poij,poii]
    u10[cc]    = wrfout.variables['U10'][0,poij,poii]
    v10[cc]    = wrfout.variables['V10'][0,poij,poii]
    T2[cc]     = wrfout.variables['T2'][0,poij,poii]
    TH2[cc]    = wrfout.variables['TH2'][0,poij,poii]
    swdwn[cc]  = wrfout.variables['SWDOWN'][0,poij,poii]
    u[cc]      = wrfdict.unstagger2d(wrfout.variables['U'][0,:,poij,poii-1:poii+1],ax=1)[:,0]
    v[cc]      = wrfdict.unstagger2d(wrfout.variables['V'][0,:,poij-1:poij+1,poii],ax=1)[:,0]
    ws         = wrfout.variables['W'][0,:,poij,poii]
    w[cc]      = (ws[1:] + ws[:-1])*0.5
    T[cc]      = wrfout.variables['T'][0,:,poij,poii]+300.0
    P[cc]      = wrfout.variables['P'][0,:,poij,poii]+wrfout.variables['PB'][0,:,poij,poii]
    Q[cc]      = wrfout.variables['QVAPOR'][0,:,poij,poii]

    cc += 1
newncdf = ncdf('%sExtractedWRFdataForWFIP2.nc' % savedir,'w',format='NETCDF4_CLASSIC')
newncdf.createDimension('NZ',nz)
newncdf.createDimension('time',nt)
newncdf.location = '(%f, %f)' % (wrfout.variables['XLAT'][0,poij,poii], wrfout.variables['XLONG'][0,poij,poii])
newncdf.elevation = '%f m' % wrfout.variables['HGT'][0,poij,poii]
newncdf.description = \
        'Extracted by Patrick Hawbecker (patrick.hawbecker@nrel.gov) on Feb. 4, 2019 from wrfout files located at %s' % fdir

times    = newncdf.createVariable('Time',np.float64, ('time',))
dates    = newncdf.createVariable('Date',np.float64, ('time',))
hgts     = newncdf.createVariable('Height',np.float64, ('NZ',))
hfxo     = newncdf.createVariable('HFX',np.float64, ('time',))
pblho    = newncdf.createVariable('PBLH',np.float64, ('time',))
psfco    = newncdf.createVariable('PSFC',np.float64, ('time',))
ustro    = newncdf.createVariable('USTAR',np.float64, ('time',))
u10o     = newncdf.createVariable('U10',np.float64, ('time',))
v10o     = newncdf.createVariable('V10',np.float64, ('time',))
T2o      = newncdf.createVariable('T2',np.float64, ('time',))
TH2o     = newncdf.createVariable('TH2',np.float64, ('time',))
swdwno   = newncdf.createVariable('SWDOWN',np.float64, ('time',))
uout     = newncdf.createVariable('U',np.float64, ('time','NZ',))
vout     = newncdf.createVariable('V',np.float64, ('time','NZ',))
wout     = newncdf.createVariable('W',np.float64, ('time','NZ',))
Tout     = newncdf.createVariable('T',np.float64, ('time','NZ',))
Pout     = newncdf.createVariable('P',np.float64, ('time','NZ',))
Qout     = newncdf.createVariable('Q',np.float64, ('time','NZ',))
times[:]   = time
dates[:]   = wdate
hgts[:]    = height
hfxo[:]    = hfx
pblho[:]   = pblh 
psfco[:]   = psfc 
ustro[:]   = ustr
u10o[:]    = u10
v10o[:]    = v10
T2o[:]     = T2
TH2o[:]    = TH2
swdwno[:]  = swdwn
uout[:]    = u
vout[:]    = v
wout[:]    = w
Tout[:]    = T
Pout[:]    = P
Qout[:]    = Q
newncdf.close()
print savedir
#    wlat,wlon = wrfdict.latlon(wrfout)
#    print poii,poij
#    plt.pcolormesh(wlon,wlat,hgt,cmap=cm.terrain)
#    plt.scatter(stnlon,stnlat,c='r',marker='x')
#    plt.scatter(wlon[poij,poii],wlat[poij,poii],c='k',marker='d')
#    plt.show()

