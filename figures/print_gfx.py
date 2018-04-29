"""
Generates polished plots from data files, output during multifluid runs. Handles all cuts at once, plotting each cut side-by-side on consistently scaled axes.
Required arguments:
	run_name
	n_grids
	nplots (I0.3 format)
	limit
	ut

Example from terminal: python3 ./figures/print_gfx.py debug 5 001 60 8.4
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from mpl_toolkits.mplot3d import axes3d
import scipy.interpolate
import glob
import sys

# physical constants
mu0 = 4.e-7*np.pi		# vacuum permeability
pmass = 1.67e-27		# proton mass

# file constants
n_headlines = 6			# Number of header lines in data files
data_dir = "./figures/data/"	# Note this script must be called from multifluid directory, as expected for a sim executable
gfx_dir = "./figures/images/"
xtn = 'eps'

# plotting constants
colormap = 'afmhot'		# Choice of color scheme for figures. afmhot is an approximate blackbody spectrum for an iron bar.
skip = 21				# number of skipped grid points for quiver/vectors
reduct = 4. 			# Reduction factor for number of data points to interpolate. Higher number means lower resolution.
n_contours = 8			# Number of contours for each plot
con_color = 'blue'		# Contour color
planet_color = 'lime'	# Color of circle at the origin
r_units = '$(R_E)$'		# Units to display on distance axes
fig_size = (12,4)		# In inches
cbar_pos = [0.94,0.1,0.02,0.75]	# Bottom-left corner x,y, w,h for colorbar, in % of fig_size

# contour level definitions
alfmax = [1., 2., 6., 10., 10.]

run_name = sys.argv[1]
n_grids_str = sys.argv[2]
fname_end = '_t' + sys.argv[3] + '.dat'
fig_end = '_t' + sys.argv[3] + '.' + xtn

limit = int(sys.argv[4])
nx = limit*2 + 1
ny = limit*2 + 1
nz = limit + 1

ut = float(sys.argv[5])

fname_start = run_name + "_"
n_grids = int(n_grids_str)

# Create figure pane and axes with appropriate general labeling:
fig = plt.subplots(1,3,figsize=fig_size,sharex=True,sharey=True)
fig = fig[0]
fig.subplots_adjust(left=0.06, right=0.92, wspace=0.25)
ax1 = plt.subplot(131)
ax2 = plt.subplot(132)
ax3 = plt.subplot(133)
ax1.set(aspect='equal')
ax2.set(aspect='equal')
ax3.set(aspect='equal')
ax1.set_title('xy',fontsize=16)
ax2.set_title('xz',fontsize=16)
ax3.set_title('yz',fontsize=16)
ax1.set_xlabel('x'+r_units)
ax1.set_ylabel('y'+r_units)
ax2.set_xlabel('x'+r_units)
ax2.set_ylabel('z'+r_units)
ax3.set_xlabel('y'+r_units)
ax3.set_ylabel('z'+r_units)


#for box in range(4,n_grids+1):
box = 4
# construct filenames
figpath_alf = gfx_dir + fname_start + 'gfx_alfmach_' + str(box) + fig_end
#	xy_fpath_plas = data_dir + fname_start + 'plas_' + str(box) + 'xy' + fname_end
xy_fpath_flow = data_dir + fname_start + 'flow_' + str(box) + 'xy' + fname_end
#	xy_fpath_pres = data_dir + fname_start + 'pres_' + str(box) + 'xy' + fname_end
#	xy_fpath_bande = data_dir + fname_start + 'bande_' + str(box) + 'xy' + fname_end
xy_fpath_model = data_dir + fname_start + 'model_' + str(box) + 'xy' + fname_end
#	xz_fpath_plas = data_dir + fname_start + 'plas_' + str(box) + 'xz' + fname_end
xz_fpath_flow = data_dir + fname_start + 'flow_' + str(box) + 'xz' + fname_end
#	xz_fpath_pres = data_dir + fname_start + 'pres_' + str(box) + 'xz' + fname_end
#	xz_fpath_bande = data_dir + fname_start + 'bande_' + str(box) + 'xz' + fname_end
xz_fpath_model = data_dir + fname_start + 'model_' + str(box) + 'xz' + fname_end
#	yz_fpath_plas = data_dir + fname_start + 'plas_' + str(box) + 'yz' + fname_end
yz_fpath_flow = data_dir + fname_start + 'flow_' + str(box) + 'yz' + fname_end
#	yz_fpath_pres = data_dir + fname_start + 'pres_' + str(box) + 'yz' + fname_end
#	yz_fpath_bande = data_dir + fname_start + 'bande_' + str(box) + 'yz' + fname_end
#	yz_fpath_model = data_dir + fname_start + 'model_' + str(box) + 'yz' + fname_end

# xy data:
#	xy_qdens,xy_qtemp,xy_hdens,xy_htemp,xy_odens,xy_otemp,xy_edens,xy_etemp = np.genfromtxt(xy_fpath_plas,unpack=True,skip_header=n_headlines)
xy_qvx,xy_qvy,xy_qvz,xy_ovx,xy_ovy,xy_ovz,xy_hvx,xy_hvy,xy_hvz,xy_tvx,xy_tvy,xy_tvz,xy_alfmach = np.genfromtxt(xy_fpath_flow,unpack=True,skip_header=n_headlines)
#	xy_qpara,xy_qperp,xy_qcross,xy_hpara,xy_hperp,xy_hcross,xy_opara,xy_operp,xy_ocross,xy_epres = np.genfromtxt(xy_fpath_pres,unpack=True,skip_header=n_headlines)
#	xy_bx,xy_by,xy_bz,xy_bmag,xy_ex,xy_ey,xy_ez,xy_emag,xy_curx,xy_cury,xy_curz = np.genfromtxt(xy_fpath_bande,unpack=True,skip_header=n_headlines)

# xz data:
#	xz_qdens,xz_qtemp,xz_hdens,xz_htemp,xz_odens,xz_otemp,xz_edens,xz_etemp = np.genfromtxt(xz_fpath_plas,unpack=True,skip_header=n_headlines)
xz_qvx,xz_qvy,xz_qvz,xz_ovx,xz_ovy,xz_ovz,xz_hvx,xz_hvy,xz_hvz,xz_tvx,xz_tvy,xz_tvz,xz_alfmach = np.genfromtxt(xz_fpath_flow,unpack=True,skip_header=n_headlines)
#	xz_qpara,xz_qperp,xz_qcross,xz_hpara,xz_hperp,xz_hcross,xz_opara,xz_operp,xz_ocross,xz_epres = np.genfromtxt(xz_fpath_pres,unpack=True,skip_header=n_headlines)
#	xz_bx,xz_by,xz_bz,xz_bmag,xz_ex,xz_ey,xz_ez,xz_emag,xz_curx,xz_cury,xz_curz = np.genfromtxt(xz_fpath_bande,unpack=True,skip_header=n_headlines)

# yz data:
#	yz_qdens,yz_qtemp,yz_hdens,yz_htemp,yz_odens,yz_otemp,yz_edens,yz_etemp = np.genfromtxt(yz_fpath_plas,unpack=True,skip_header=n_headlines)
yz_qvx,yz_qvy,yz_qvz,yz_ovx,yz_ovy,yz_ovz,yz_hvx,yz_hvy,yz_hvz,yz_tvx,yz_tvy,yz_tvz,yz_alfmach = np.genfromtxt(yz_fpath_flow,unpack=True,skip_header=n_headlines)
#	yz_qpara,yz_qperp,yz_qcross,yz_hpara,yz_hperp,yz_hcross,yz_opara,yz_operp,yz_ocross,yz_epres = np.genfromtxt(yz_fpath_pres,unpack=True,skip_header=n_headlines)
#	yz_bx,yz_by,yz_bz,yz_bmag,yz_ex,yz_ey,yz_ez,yz_emag,yz_curx,yz_cury,yz_curz = np.genfromtxt(yz_fpath_bande,unpack=True,skip_header=n_headlines)

# General position data:
x,y,z,r,bdx,bdy,bdz,bdipx,bdipy,bdipz,bdipmag = np.genfromtxt(xy_fpath_model,unpack=True,skip_header=n_headlines)
xmin = x.min()
xmax = x.max()
xi = np.linspace(xmin,xmax,nx)
ymin = y.min()
ymax = y.max()
yi = np.linspace(ymin,ymax,ny)
# Now we need to import a different cut to get max/min z values:
x,y,z,r,bdx,bdy,bdz,bdipx,bdipy,bdipz,bdipmag = np.genfromtxt(xz_fpath_model,unpack=True,skip_header=n_headlines)	
zmin = z.min()
zmax = z.max()
zi = np.linspace(zmin,zmax,nz)
deltax = xmax - xmin
deltay = ymax - ymin
deltaz = zmax - zmin

# Set ticks for all plots:
ax1.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltax/4))
ax1.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(deltax/8))
ax1.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltay/4))
ax1.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(deltay/8))
ax2.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltax/4))
ax2.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(deltax/8))
ax2.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltaz/4))
ax2.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(deltaz/8))
ax3.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltay/4))
ax3.xaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(deltay/8))
ax3.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltaz/4))
ax3.yaxis.set_minor_locator(matplotlib.ticker.MultipleLocator(deltaz/8))

# Alfven speed:
plot_title = run_name + r' Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + str(ut)
cbar_title = '|v|/V_Alf'
cbmin = 0.0
cbmax = alfmax[box-1]
contour_levels = np.arange(cbmin,cbmax,cbmax/n_contours)
hi1 = np.reshape(xy_alfmach,[ny,nx])
hi2 = np.reshape(xz_alfmach,[nz,nx])
hi3 = np.reshape(yz_alfmach,[nz,ny])

# Select each subplot and draw a contour on it while selected:
conxy = ax1.imshow(hi1,origin='lower',vmin=cbmin,vmax=cbmax,extent=[xmin,xmax,ymin,ymax],cmap=colormap,interpolation='bicubic')
ax1.contour(xi,yi,hi1,colors=con_color,linewidths=0.15,levels=contour_levels)

ax2.imshow(hi2,origin='lower',vmin=cbmin,vmax=cbmax,extent=[xmin,xmax,zmin,zmax],cmap=colormap,interpolation='bicubic')
ax2.contour(xi,zi,hi2,colors=con_color,linewidths=0.15,levels=contour_levels)

ax3.imshow(hi3,origin='lower',vmin=cbmin,vmax=cbmax,extent=[ymin,ymax,zmin,zmax],cmap=colormap,interpolation='bicubic')
ax3.contour(yi,zi,hi3,colors=con_color,linewidths=0.15,levels=contour_levels)

plt.suptitle(plot_title,fontsize=20)
ax1.scatter( 0,0, s=4.*np.pi, c=planet_color, zorder=10)
ax2.scatter( 0,0, s=4.*np.pi, c=planet_color, zorder=10)
ax3.scatter( 0,0, s=4.*np.pi, c=planet_color, zorder=10)

cbar_ax = fig.add_axes(cbar_pos)
cbar = plt.colorbar(conxy, ax=(ax1,ax2,ax3), cax=cbar_ax)
cbar.ax.set_title(cbar_title,size=14)
fig.savefig(figpath_alf,format=xtn,dpi=200)
