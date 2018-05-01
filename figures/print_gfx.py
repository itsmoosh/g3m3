"""
Generates polished plots from data files, output during multifluid runs. Handles all cuts at once, plotting each cut side-by-side on consistently scaled axes.
Required arguments:
	run_name
	n_grids
	nplots (I0.3 format)
	limit
	r_inner (in units of planetary radii)
	ut
	diagnostic (whether to print 3D plots)
	update_gifs (whether to generate gifs)

Example from terminal: python3 ./figures/print_gfx.py debug 5 001 60 2.0 8.4 False
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import axes3d
from matplotlib.patches import Circle, PathPatch
from matplotlib import ticker
import scipy.interpolate
import glob
import sys
import os
import gfx_functions as gfx

t_start = os.times().elapsed

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
r_units = ' $(R_E)$'	# Units to display on distance axes
fig_size = (12,4)		# In inches
fig3d_size = (10,7.5)
cbar_pos = [0.94,0.1,0.02,0.75]	# Bottom-left corner x,y, w,h for colorbar, in % of fig_size
cbar3d_pos = [0.9,0.15,0.024,0.4]
show_vecs = True		# Whether to print vectors or contours (for vector quantities)
fig_dpi = 200

# Contour level definitions
# Note: These values are appropriate for Saturn, probably not for Europa.
alfmax = [1., 2., 6., 10., 10.]	# Mach number
#speeds = [200., 250., 450., 450., 450.]	# km/s
bmags = [800., 200., 100., 100., 100.]

run_name = sys.argv[1]
n_grids_str = sys.argv[2]
nplots_str = sys.argv[3]

limit = int(sys.argv[4])
nx = limit*2 + 1
ny = limit*2 + 1
nz = limit + 1

r_inner = float(sys.argv[5])

ut = float(sys.argv[6])

if(sys.argv[7] == 'True'):
	diagnostic = True
	cbar_pos = cbar3d_pos
	fig_size = fig3d_size
else:
	diagnostic = False
crop = diagnostic

if(sys.argv[8] == 'True'):
	update_gifs = True
else:
	update_gifs = False

fname_start = run_name + "_"
n_grids = int(n_grids_str)
fname_end = '_t' + nplots_str + '.dat'
fig_end = '_t' + nplots_str + '.' + xtn
plot_opt = (skip,n_contours,con_color,colormap,planet_color,r_inner)

for box in range(1,n_grids+1):
	# Construct figure file paths
	gfxp1 = gfx_dir + fname_start + 'gfx_'
	gfxp2 = str(box) + fig_end
	qty_list = ('alfmach', 'qvel', 'hvel', 'ovel', 'bmag')
	figpath_alf = gfxp1 + 'alfmach' + gfxp2
	figpath_qvel = gfxp1 + 'qvel' + gfxp2
	figpath_hvel = gfxp1 + 'hvel' + gfxp2
	figpath_ovel = gfxp1 + 'ovel' + gfxp2
	figpath_bnet = gfxp1 + 'bnet' + gfxp2

	# Construct data file paths
	datap1 = data_dir + fname_start
	dataxy2 = '_' + str(box) + 'xy' + fname_end
	dataxz2 = '_' + str(box) + 'xz' + fname_end
	datayz2 = '_' + str(box) + 'yz' + fname_end
	#	xy_fpath_plas = datap1 + 'plas' + dataxy2
	xy_fpath_flow = datap1 + 'flow' + dataxy2
	#	xy_fpath_pres = datap1 + 'pres' + dataxy2
	xy_fpath_bande = datap1 + 'bande' + dataxy2
	xy_fpath_model = datap1 + 'model' + dataxy2
	#	xz_fpath_plas = datap1 + 'plas' + dataxz2
	xz_fpath_flow = datap1 + 'flow' + dataxz2
	#	xz_fpath_pres = datap1 + 'pres' + dataxz2
	xz_fpath_bande = datap1 + 'bande' + dataxz2
	xz_fpath_model = datap1 + 'model' + dataxz2
	#	yz_fpath_plas = datap1 + 'plas' + datayz2
	yz_fpath_flow = datap1 + 'flow' + datayz2
	#	yz_fpath_pres = datap1 + 'pres' + datayz2
	yz_fpath_bande = datap1 + 'bande' + datayz2
	yz_fpath_model = datap1 + 'model' + datayz2

	# xy data:
	#	xy_qdens,xy_qtemp,xy_hdens,xy_htemp,xy_odens,xy_otemp,xy_edens,xy_etemp = np.genfromtxt(xy_fpath_plas,unpack=True,skip_header=n_headlines)
	xy_qvx,xy_qvy,xy_qvz,xy_ovx,xy_ovy,xy_ovz,xy_hvx,xy_hvy,xy_hvz,xy_tvx,xy_tvy,xy_tvz,xy_alfmach = np.genfromtxt(xy_fpath_flow,unpack=True,skip_header=n_headlines)
	#	xy_qpara,xy_qperp,xy_qcross,xy_hpara,xy_hperp,xy_hcross,xy_opara,xy_operp,xy_ocross,xy_epres = np.genfromtxt(xy_fpath_pres,unpack=True,skip_header=n_headlines)
	xy_bx,xy_by,xy_bz,xy_bmag,xy_ex,xy_ey,xy_ez,xy_emag,xy_curx,xy_cury,xy_curz = np.genfromtxt(xy_fpath_bande,unpack=True,skip_header=n_headlines)

	# xz data:
	#	xz_qdens,xz_qtemp,xz_hdens,xz_htemp,xz_odens,xz_otemp,xz_edens,xz_etemp = np.genfromtxt(xz_fpath_plas,unpack=True,skip_header=n_headlines)
	xz_qvx,xz_qvy,xz_qvz,xz_ovx,xz_ovy,xz_ovz,xz_hvx,xz_hvy,xz_hvz,xz_tvx,xz_tvy,xz_tvz,xz_alfmach = np.genfromtxt(xz_fpath_flow,unpack=True,skip_header=n_headlines)
	#	xz_qpara,xz_qperp,xz_qcross,xz_hpara,xz_hperp,xz_hcross,xz_opara,xz_operp,xz_ocross,xz_epres = np.genfromtxt(xz_fpath_pres,unpack=True,skip_header=n_headlines)
	xz_bx,xz_by,xz_bz,xz_bmag,xz_ex,xz_ey,xz_ez,xz_emag,xz_curx,xz_cury,xz_curz = np.genfromtxt(xz_fpath_bande,unpack=True,skip_header=n_headlines)

	# yz data:
	#	yz_qdens,yz_qtemp,yz_hdens,yz_htemp,yz_odens,yz_otemp,yz_edens,yz_etemp = np.genfromtxt(yz_fpath_plas,unpack=True,skip_header=n_headlines)
	yz_qvx,yz_qvy,yz_qvz,yz_ovx,yz_ovy,yz_ovz,yz_hvx,yz_hvy,yz_hvz,yz_tvx,yz_tvy,yz_tvz,yz_alfmach = np.genfromtxt(yz_fpath_flow,unpack=True,skip_header=n_headlines)
	#	yz_qpara,yz_qperp,yz_qcross,yz_hpara,yz_hperp,yz_hcross,yz_opara,yz_operp,yz_ocross,yz_epres = np.genfromtxt(yz_fpath_pres,unpack=True,skip_header=n_headlines)
	yz_bx,yz_by,yz_bz,yz_bmag,yz_ex,yz_ey,yz_ez,yz_emag,yz_curx,yz_cury,yz_curz = np.genfromtxt(yz_fpath_bande,unpack=True,skip_header=n_headlines)

	# Position data:
	xy_x,xy_y,xy_z,xy_r,xy_bdx,xy_bdy,xy_bdz,xy_bdipx,xy_bdipy,xy_bdipz,xy_bdipmag = np.genfromtxt(xy_fpath_model,unpack=True,skip_header=n_headlines)	
	xz_x,xz_y,xz_z,xz_r,xz_bdx,xz_bdy,xz_bdz,xz_bdipx,xz_bdipy,xz_bdipz,xz_bdipmag = np.genfromtxt(xz_fpath_model,unpack=True,skip_header=n_headlines)
	yz_x,yz_y,yz_z,yz_r,yz_bdx,yz_bdy,yz_bdz,yz_bdipx,yz_bdipy,yz_bdipz,yz_bdipmag = np.genfromtxt(yz_fpath_model,unpack=True,skip_header=n_headlines)

	xmin = xy_x.min()
	xmax = xy_x.max()
	xi = np.linspace(xmin,xmax,nx)
	ymin = xy_y.min()
	ymax = xy_y.max()
	yi = np.linspace(ymin,ymax,ny)
	# We need a different cut to get max/min z values:
	zmin = xz_z.min()
	zmax = xz_z.max()
	zi = np.linspace(zmin,zmax,nz)
	deltax = xmax - xmin
	deltay = ymax - ymin
	deltaz = zmax - zmin
	minmax = (xmin,xmax,ymin,ymax,zmin,zmax)
	deltas = (deltax,deltay,deltaz)
	pos = (xi,yi,zi)

	# Alfven speed:
	fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
	cbar_title = r'$\frac{|\vec{v}_\mathrm{net}|}{V_\mathrm{Alf}}$'
	cbmin = 0.0
	cbmax = alfmax[box-1]
	cbparams = (cbmin,cbmax,cbar_pos,cbar_title)
	hi1 = np.reshape(xy_alfmach,[ny,nx])
	hi2 = np.reshape(xz_alfmach,[nz,nx])
	hi3 = np.reshape(yz_alfmach,[nz,ny])
	if(show_vecs):
		if(diagnostic):
			values = (
				((xy_x,xy_y,xy_z),(xy_tvx,xy_tvy,xy_tvz)),
				((xz_x,xz_y,xz_z),(xz_tvx,xz_tvy,xz_tvz)),
				((yz_x,yz_y,yz_z),(yz_tvx,yz_tvy,yz_tvz)),
				hi1,hi2,hi3 )		
		else:
				# showstreams patterning (diagnostic = False only):
				vec1 = np.stack((xy_tvx,xy_tvy,xy_tvz),axis=1)
				vec2 = np.stack((xz_tvx,xz_tvy,xz_tvz),axis=1)
				vec3 = np.stack((yz_tvx,yz_tvy,yz_tvz),axis=1)
				vec1 = np.reshape(vec1,[ny,nx,3])
				vec2 = np.reshape(vec2,[nz,nx,3])
				vec3 = np.reshape(vec3,[nz,ny,3])
				values = ( vec1,vec2,vec3, hi1,hi2,hi3 )
		plot_title = run_name + r' $\vec{v}_\mathrm{net}$ streamlines and Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + str(ut)
	else:
		plot_title = run_name + r' Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + str(ut)
		values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,show_vecs,True)
	gfx.save_fig(fig,figpath_alf,xtn,fig_dpi,crop)


	# q velocity:
	fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
	plot_title = run_name + r' species q velocity, box ' + str(box) + ', ut = ' + str(ut)
	cbar_title = r'$|\vec{v}_q|$ (km/s)'
	xy_qspd = np.sqrt( xy_qvx**2 + xy_qvy**2 + xy_qvz**2 )
	xz_qspd = np.sqrt( xz_qvx**2 + xz_qvy**2 + xz_qvz**2 )
	yz_qspd = np.sqrt( yz_qvx**2 + yz_qvy**2 + yz_qvz**2 )
	cbmin = xy_qspd.min()
	cbmax = xy_qspd.max()
	cbparams = (cbmin,cbmax,cbar_pos,cbar_title)
	hi1 = np.reshape(xy_qspd,[ny,nx])
	hi2 = np.reshape(xz_qspd,[nz,nx])
	hi3 = np.reshape(yz_qspd,[nz,ny])
	if(show_vecs):
		values = (
			((xy_x,xy_y,xy_z),(xy_qvx,xy_qvy,xy_qvz)),
			((xz_x,xz_y,xz_z),(xz_qvx,xz_qvy,xz_qvz)),
			((yz_x,yz_y,yz_z),(yz_qvx,yz_qvy,yz_qvz)),
			hi1,hi2,hi3 )
		gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,True,False)
	else:
		values = (hi1,hi2,hi3)
		gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False)
	gfx.save_fig(fig,figpath_qvel,xtn,fig_dpi,crop)


	# h velocity:
	fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
	plot_title = run_name + r' species h velocity, box ' + str(box) + ', ut = ' + str(ut)
	cbar_title = r'$|\vec{v}_h|$ (km/s)'
	xy_hspd = np.sqrt( xy_hvx**2 + xy_hvy**2 + xy_hvz**2 )
	xz_hspd = np.sqrt( xz_hvx**2 + xz_hvy**2 + xz_hvz**2 )
	yz_hspd = np.sqrt( yz_hvx**2 + yz_hvy**2 + yz_hvz**2 )
	cbmin = xy_hspd.min()
	cbmax = xy_hspd.max()
	cbparams = (cbmin,cbmax,cbar_pos,cbar_title)
	hi1 = np.reshape(xy_hspd,[ny,nx])
	hi2 = np.reshape(xz_hspd,[nz,nx])
	hi3 = np.reshape(yz_hspd,[nz,ny])
	if(show_vecs):
		values = (
			((xy_x,xy_y,xy_z),(xy_hvx,xy_hvy,xy_hvz)),
			((xz_x,xz_y,xz_z),(xz_hvx,xz_hvy,xz_hvz)),
			((yz_x,yz_y,yz_z),(yz_hvx,yz_hvy,yz_hvz)),
			hi1,hi2,hi3 )
		gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,True,False)
	else:
		values = (hi1,hi2,hi3)
		gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False)
	gfx.save_fig(fig,figpath_hvel,xtn,fig_dpi,crop)


#	# o velocity:
	fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
	plot_title = run_name + r' species o velocity, box ' + str(box) + ', ut = ' + str(ut)
	cbar_title = r'$|\vec{v}_o|$ (km/s)'
	xy_ospd = np.sqrt( xy_ovx**2 + xy_ovy**2 + xy_ovz**2 )
	xz_ospd = np.sqrt( xz_ovx**2 + xz_ovy**2 + xz_ovz**2 )
	yz_ospd = np.sqrt( yz_ovx**2 + yz_ovy**2 + yz_ovz**2 )
	cbmin = xy_ospd.min()
	cbmax = xy_ospd.max()
	cbparams = (cbmin,cbmax,cbar_pos,cbar_title)
	hi1 = np.reshape(xy_ospd,[ny,nx])
	hi2 = np.reshape(xz_ospd,[nz,nx])
	hi3 = np.reshape(yz_ospd,[nz,ny])
	if(show_vecs):
		values = (
			((xy_x,xy_y,xy_z),(xy_ovx,xy_ovy,xy_ovz)),
			((xz_x,xz_y,xz_z),(xz_ovx,xz_ovy,xz_ovz)),
			((yz_x,yz_y,yz_z),(yz_ovx,yz_ovy,yz_ovz)),
			hi1,hi2,hi3 )
		gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,True,False)
	else:
		values = (hi1,hi2,hi3)
		gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False)
	gfx.save_fig(fig,figpath_ovel,xtn,fig_dpi,crop)

	# magnetic field:
	fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
	plot_title = run_name + r' net magnetic field, box ' + str(box) + ', ut = ' + str(ut)
	cbar_title = r'$|\vec{B}_\mathrm{net}|$ (nT)'
	cbmin = 0.0
	cbmax = bmags[box-1]
	cbparams = (cbmin,cbmax,cbar_pos,cbar_title)
	hi1 = np.reshape(xy_bmag,[ny,nx])
	hi2 = np.reshape(xz_bmag,[nz,nx])
	hi3 = np.reshape(yz_bmag,[nz,ny])
	if(diagnostic):
		values = (
			((xy_x,xy_y,xy_z),(xy_bx,xy_by,xy_bz)),
			((xz_x,xz_y,xz_z),(xz_bx,xz_by,xz_bz)),
			((yz_x,yz_y,yz_z),(yz_bx,yz_by,yz_bz)),
			hi1,hi2,hi3 )		
	else:
		# showstreams patterning (diagnostic = False only):
		vec1 = np.stack((xy_bx,xy_by,xy_bz),axis=1)
		vec2 = np.stack((xz_bx,xz_by,xz_bz),axis=1)
		vec3 = np.stack((yz_bx,yz_by,yz_bz),axis=1)
		vec1 = np.reshape(vec1,[ny,nx,3])
		vec2 = np.reshape(vec2,[nz,nx,3])
		vec3 = np.reshape(vec3,[nz,ny,3])
		values = ( vec1,vec2,vec3, hi1,hi2,hi3 )
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,True,True)
	gfx.save_fig(fig,figpath_bnet,xtn,fig_dpi,crop)


# Create gifs for the past 10 figures
if(update_gifs):
	if((xtn == 'png') or (xtn == 'gif')):
		gfx.upd_gifs(qty_list, run_name, nplots_str, gfx_dir, n_grids, xtn)
	else:
		print("gif conversion not supported for "+xtn+". gifs not updated.")

t_end = os.times().elapsed
t_elapsed = (t_end - t_start) / 60
print("Debug: script duration " + str(t_elapsed) + " min")
