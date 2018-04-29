"""
Generates diagnostic plots from data files, output during multifluid runs. Handles all cuts at once, plotting them all on common xyz axes.
Required arguments:
	run_name
	n_grids
	nplots (I0.3 format)
	limit
	ut

Example from terminal: python3 ./figures/diag_gfx.py debug 5 001 60 8.4
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
from matplotlib.transforms import Bbox
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
xtn = 'png'

# plotting constants
colormap = 'afmhot'		# Choice of color scheme for figures. afmhot is an approximate blackbody spectrum for an iron bar.
skip = 21				# number of skipped grid points for quiver/vectors
reduct = 4. 			# Reduction factor for number of data points to interpolate. Higher number means lower resolution.
n_contours = 8			# Number of contours for each plot
con_color = 'blue'		# Contour color
planet_color = 'lime'	# Color of circle at the origin
r_units = '$(R_E)$'		# Units to display on distance axes
fig_size = (10,7.5)		# In inches
cbar_pos = [0.9,0.2,0.024,0.4]	# Bottom-left corner x,y, w,h for colorbar, in % of fig_size


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

# Create figure pane and 3D axes with appropriate labeling:
fig = plt.figure(figsize=fig_size)
# Roughly figure dimensions times these numbers determines the bounding box size.
# Uncomment here to see crop box for adjusting purposes.
#ax_crop = fig.add_axes([0.25, 0.1, 0.675, 0.625])
ax = fig.gca(projection='3d')
ax.set_aspect('equal')	# Has no effect in matplotlib v2.2.2 due to a bug.
gp = 0.65	# Gray % to fill in axes backgrounds on 3D plots
ax.w_xaxis.set_pane_color((gp+0.14,gp+0.14,gp+0.14,gp+0.14))
ax.w_zaxis.set_pane_color((gp+0.07,gp+0.07,gp+0.07,gp+0.07))
ax.w_yaxis.set_pane_color((gp,gp,gp,gp))
ax.set_xlabel('x'+r_units+r', wind $\rightarrow$')
ax.set_ylabel('y'+r_units)
ax.set_zlabel('z'+r_units)


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

#	# yz data:
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

# Tile 1D arrays to form consistent 2D arrays for use in making contours:
dbl_x = np.repeat(xi, nz)
dbl_x = np.reshape(dbl_x,[nx,nz])
dbl_x = np.transpose(dbl_x)
dbl_y = np.repeat(yi, nz)
dbl_y = np.reshape(dbl_y,[nx,nz])
dbl_y = np.transpose(dbl_y)
dbl_z = np.repeat(zi, nx)
dbl_z = np.reshape(dbl_z,[nz,nx])

# Set axis display limits and ticks:
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, 3*zmax)
ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltax/8))
ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltay/8))
ax.zaxis.set_major_locator(matplotlib.ticker.MultipleLocator(deltaz/4))
#ax.view_init(elev=20., azim=-70)

# Alfven speed:
plot_title = run_name + r': Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + str(ut)
cbar_title = '|v|/V_Alf'
cbmin = 0.0
cbmax = alfmax[box-1]
hi1 = np.reshape(xy_alfmach,[ny,nx])
hi2 = np.reshape(xz_alfmach,[nz,nx])
hi3 = np.reshape(yz_alfmach,[nz,ny])

# Generate contours for all 3 cuts:
contour_levels = np.arange(cbmin,cbmax,cbmax/n_contours)
fine_lvls = np.arange(cbmin,cbmax,cbmax/100)
plt.title(plot_title,fontsize=20,y=0.72, bbox=dict(facecolor='white'))
conxy = ax.contourf(xi, yi, hi1, zdir='z', offset=zmin, cmap=colormap, vmin=0., vmax=10., levels=fine_lvls)
ax.scatter( 0,0,zmin+1, s=8*np.pi, c=planet_color)
conxz = ax.contourf(dbl_x, hi2, dbl_z, zdir='y', offset=ymax, cmap=colormap, vmin=0., vmax=10., levels=fine_lvls)
ax.scatter( 0,ymax-1,0, s=8*np.pi, c=planet_color)	# Not visible due to a bug in zorder with Axes3D as of 4/28/18, mpl v2.2.2
conyz = ax.contourf(hi3, dbl_y, dbl_z, zdir='x', offset=xmin, cmap=colormap, vmin=0., vmax=10., levels=fine_lvls)
ax.scatter( xmin+1,0,0, s=8*np.pi, c=planet_color)

# Add colorbar, crop and save figure:
cbar_ax = fig.add_axes(cbar_pos)
cbar = fig.colorbar(conxy, cax=cbar_ax)
cbar.ax.set_title(cbar_title,size=16)
cbar.ax.tick_params(labelsize=14)
crop_bbox = Bbox.from_bounds(2.45, 0.75, 7.27, 4.5)
fig.savefig(figpath_alf, format=xtn, dpi=200, bbox_inches=crop_bbox)
