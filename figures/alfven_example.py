"""
Example to plot alfven speed in with vector field or contours

"""

# required imports

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import scipy.interpolate
import glob

#
# constants
#

temp_equiv = 1.0437e4   # temp_equiv normalization from multifluid
rho_equiv = 0.2e6       # temp_equiv normalization from multifluid
N = 121                 # size of grids in multifluid, equal to 1 + 2*limit
skip = 21               # number of skipped grid points for quiver/vectors
mu0 = 4*np.pi*1.e-7     # vacuum permeability
pmass = 1.67e-27        # proton mass

#
# setup the figure parameters (figsize defaults in inches)
# for multiple figures, you will want to change the first
# two indices and the figure size to compensate for scaling, e.g.,
# for a 2 row, 3 column figure:
#
#       f,axes = plt.subplots(2, 3, figsize(13,8), sharex=True, sharey=True)
#

f,axes = plt.subplots(1,1,figsize=(6,5),sharex=True,sharey=True)

#
# read in data files
#

x,y,qdens,qtemp,hdens,htemp,odens,otemp,edens,etemp = np.genfromtxt(r'./data/denstemp.dat',unpack=True,skip_header=1)
x,y,qvx,qvy,qvz,ovx,ovy,ovz,hvx,hvy,hvz = np.genfromtxt(r'./data/flowdata.dat',unpack=True,skip_header=1)
x,y,bx,by,bz,curx,cury,curz = np.genfromtxt(r'./data/fielddata.dat',unpack=True,skip_header=1)

#
# calculate equatorial flow speed and alfven speed
# so we can plot the ratio
#

alfspd = 1.e-9*np.sqrt(bx**2+by**2+bz**2)/np.sqrt(mu0*pmass*rho_equiv*(10**qdens+18*10**hdens+32*10**odens))
flowspd = 1.e6*np.sqrt( (0.333*(qvx+hvx+ovx))**2+(0.333*(qvy+hvy+ovy))**2)

#
# create plot
#

ax1 = plt.subplot(111)
ax1.set(aspect='equal')
xi = np.linspace(x.min(),x.max(),N)
yi = np.linspace(y.min(),y.max(),N)
zi = scipy.interpolate.griddata((x,y),flowspd/alfspd,(xi[None,:],yi[:,None]),method='cubic')

# uncomment to plot vector field
#plt.quiver(x[::skip],y[::skip],0.333e6*(qvx[::skip]+hvx[::skip]+ovx[::skip]),0.333e6*(qvy[::skip]+hvy[::skip]+ovy[::skip]),scale=2.e6,headwidth=5,color='blue')

# uncomment to plot contours
plt.contour(xi,yi,zi,colors='blue',linewidths=0.2,levels=np.arange(0.,8.,1.))

plt.imshow(zi,origin='lower',vmin=0.0,vmax=8.,extent=[-96.,96.,-96.,96.],cmap='afmhot',interpolation='bicubic')
ax1.add_artist( plt.Circle((0,0), 2., color='k'))
ax1.set_title(r'Aniso Alfv$\'{e}$n Mach',fontsize=16)
plt.tick_params(axis='both', which='minor', labelsize=12,direction='in')
plt.tick_params(axis='both', which='major', labelsize=12,direction='in')
ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(18))
ax1.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(6))
ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(16))
ax1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(8))
plt.yticks(np.arange(-72,73,18))
plt.colorbar()

plt.tight_layout()
plt.savefig('./plot_alfven_example.eps',format='eps',dpi=100)
