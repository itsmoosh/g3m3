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

Example from terminal: python3 ./figures/print_gfx.py debug 5 001 60 2.0 8.4 False False
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
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

# physical constants
mu0 = 4.e-7*np.pi		# vacuum permeability
pmass = 1.67e-27		# proton mass

# file constants
n_headlines = 6			# Number of header lines in data files
data_dir = "./figures/data/"	# Note this script must be called from multifluid directory, as expected for a sim executable
gfx_dir = "./figures/images/"
python_dir = "./figures/"
xtn = 'png'

# plotting constants
colormap = 'afmhot'		# Choice of color scheme for figures. afmhot is an approximate blackbody spectrum for an iron bar.
skip = 21				# number of skipped grid points for quiver/vectors
reduct = 4. 			# Reduction factor for number of data points to interpolate. Higher number means lower resolution.
show_vecs = True		# Whether to print vectors or contours (for vector quantities)
showstreams = False		# Whether to print streamlines to show fields and flows
stack = True			# Whether to print multiple rows of quantities on the same plot
draw_boxes = False		# Whether to draw the smaller boxes
fig_dpi = 200
n_contours = 8			# Number of contours for each plot
con_color = 'blue'		# Contour color
planet_color = 'lime'	# Color of circle at the origin
r_units = ' $(R_E)$'	# Units to display on distance axes
fig_size = (12,4)		# In inches
fig3d_size = (10,7.5)
vsp = 0.1				# Vertical space between graphical elements (e.g. colorbars)
cbx = 0.92
cby = vsp/2 + 0.0
cbw = 0.02
cbh = 0.75
cbar_pos = [cbx,cby,cbw,cbh+0.07]	# Bottom-left corner x,y, w,h for colorbar, in % of fig_size
cbar_pos22 = [cbx,cby+vsp/4,cbw,cbh/2]	# Same as above, but for the bottom row of a 2-stack plot
cbar_pos21 = [cbx,cby+cbh/2+vsp*1.25,cbw,cbh/2]	# For stack 2/2
cbar_pos31 = [cbx,cby*0.9+2*cbh/3+1.1*vsp,cbw,cbh/3]	# For top row of 3-stack
cbar_pos32 = [cbx,cby*0.9+cbh/3+1.1*vsp/2,cbw,cbh/3]	# For middle row of 3-stack
cbar_pos33 = [cbx,cby*0.9,cbw,cbh/3]	# For bottom row of 3-stack
cbar3d_pos = [0.9,0.15,0.024,0.4]	# For 3d diagnostic plots

# Contour level definitions
# Note: These values are appropriate for Saturn, probably not for Europa.
alfmax = [1., 2., 6., 10., 10.]	# Mach number
speeds = [200., 250., 450., 450., 450.]	# km/s
bmags = [800., 200., 100., 100., 100.]	# nT
currents = [0.2, 0.1, 0.07, 0.05, 0.02]	# nA/m^2
press = [3., 2., 2., 2., 2.]	# +/- log units of nPa
qpres_mid = -12	# Midpoint in log (nPa) of q pressure colorbars
hpres_mid = -13
opres_mid = -12
epres_mid = -12
emags = [30., 15., 7.5, 4., 1.5]	# mV/m
temps = [16., 14., 12., 8., 6.]	# keV
densities = [4., 4., 4., 4., 4.]	# +/- log units of #/cc
qdens_mid = 0
hdens_mid = -2
odens_mid = -1

run_name = sys.argv[1]
n_grids_str = sys.argv[2]
nplots_str = sys.argv[3]

limit = int(sys.argv[4])
nx = limit*2 + 1
ny = limit*2 + 1
nz = limit + 1

r_inner = float(sys.argv[5])

ut = float(sys.argv[6])
ut_str = sys.argv[6]

if(sys.argv[7] == 'True'):
	diagnostic = True
	cbar_pos = cbar3d_pos
	fig_size = fig3d_size
	stack = False
else:
	diagnostic = False
crop = diagnostic

if(sys.argv[8] == 'True'):
	update_gifs = True
else:
	update_gifs = False

if(stack):
	rows = 3
	brows = 2
else:
	rows = 1
	brows = 1

if(showstreams):
	vecs_title = "streamlines"
else:
	vecs_title = "vectors"

fname_start = run_name + "_"
n_grids = int(n_grids_str)
fname_end = '_t' + nplots_str + '.dat'
fig_end = '_t' + nplots_str + '.' + xtn
plot_opt = (skip,n_contours,con_color,colormap,planet_color,r_inner)
boxes_to_draw = ()
box_draw = (False,((),))

for box in range(1,n_grids+1):
	# Construct figure file paths
	gfxp1 = gfx_dir + fname_start + 'gfx_'
	gfxp2 = str(box) + fig_end
	if(stack):
		qty_list = ('bandalf', 'flows', 'qpres', 'hpres', 'opres',
			'elec', 'temps', 'dens')
		figpath_bandalf = gfxp1 + 'bandalf' + gfxp2
		figpath_flows = gfxp1 + 'flows' + gfxp2
		figpath_qpres = gfxp1 + 'qpres' + gfxp2
		figpath_hpres = gfxp1 + 'hpres' + gfxp2
		figpath_opres = gfxp1 + 'opres' + gfxp2
		figpath_elec = gfxp1 + 'elec' + gfxp2
		figpath_temps = gfxp1 + 'temps' + gfxp2
		figpath_dens = gfxp1 + 'dens' + gfxp2
	else:
		qty_list = ('alfmach', 'qvel', 'hvel', 'ovel', 'bnet', 'curr',
			'qpara', 'qperp', 'qcross',
			'hpara', 'hperp', 'hcross',
			'opara', 'operp', 'ocross',
			'efld',  'epres', 'etemp',
			'qtemp', 'htemp', 'otemp',
			'qdens', 'hdens', 'odens')
		figpath_alf = gfxp1 + 'alfmach' + gfxp2
		figpath_qvel = gfxp1 + 'qvel' + gfxp2
		figpath_hvel = gfxp1 + 'hvel' + gfxp2
		figpath_ovel = gfxp1 + 'ovel' + gfxp2
		figpath_bnet = gfxp1 + 'bnet' + gfxp2
		figpath_curr = gfxp1 + 'curr' + gfxp2
		figpath_qpara = gfxp1 + 'qpara' + gfxp2
		figpath_qperp = gfxp1 + 'qperp' + gfxp2
		figpath_qcross = gfxp1 + 'qcross' + gfxp2
		figpath_hpara = gfxp1 + 'hpara' + gfxp2
		figpath_hperp = gfxp1 + 'hperp' + gfxp2
		figpath_hcross = gfxp1 + 'hcross' + gfxp2
		figpath_opara = gfxp1 + 'opara' + gfxp2
		figpath_operp = gfxp1 + 'operp' + gfxp2
		figpath_ocross = gfxp1 + 'ocross' + gfxp2
		figpath_efld = gfxp1 + 'efld' + gfxp2
		figpath_epres = gfxp1 + 'epres' + gfxp2
		figpath_etemp = gfxp1 + 'etemp' + gfxp2
		figpath_qtemp = gfxp1 + 'qtemp' + gfxp2
		figpath_htemp = gfxp1 + 'htemp' + gfxp2
		figpath_otemp = gfxp1 + 'otemp' + gfxp2
		figpath_qdens = gfxp1 + 'qdens' + gfxp2
		figpath_hdens = gfxp1 + 'hdens' + gfxp2
		figpath_odens = gfxp1 + 'odens' + gfxp2


	# Construct data file paths
	datap1 = data_dir + fname_start
	dataxy2 = '_' + str(box) + 'xy' + fname_end
	dataxz2 = '_' + str(box) + 'xz' + fname_end
	datayz2 = '_' + str(box) + 'yz' + fname_end
	xy_fpath_plas = datap1 + 'plas' + dataxy2
	xy_fpath_flow = datap1 + 'flow' + dataxy2
	xy_fpath_pres = datap1 + 'pres' + dataxy2
	xy_fpath_bande = datap1 + 'bande' + dataxy2
	xy_fpath_model = datap1 + 'model' + dataxy2
	xz_fpath_plas = datap1 + 'plas' + dataxz2
	xz_fpath_flow = datap1 + 'flow' + dataxz2
	xz_fpath_pres = datap1 + 'pres' + dataxz2
	xz_fpath_bande = datap1 + 'bande' + dataxz2
	xz_fpath_model = datap1 + 'model' + dataxz2
	yz_fpath_plas = datap1 + 'plas' + datayz2
	yz_fpath_flow = datap1 + 'flow' + datayz2
	yz_fpath_pres = datap1 + 'pres' + datayz2
	yz_fpath_bande = datap1 + 'bande' + datayz2
	yz_fpath_model = datap1 + 'model' + datayz2

	# xy data:
	xy_qdens,xy_qtemp,xy_hdens,xy_htemp,xy_odens,xy_otemp,xy_edens,xy_etemp = np.genfromtxt(xy_fpath_plas,unpack=True,skip_header=n_headlines)
	xy_qvx,xy_qvy,xy_qvz,xy_ovx,xy_ovy,xy_ovz,xy_hvx,xy_hvy,xy_hvz,xy_tvx,xy_tvy,xy_tvz,xy_alfmach = np.genfromtxt(xy_fpath_flow,unpack=True,skip_header=n_headlines)
	xy_qpara,xy_qperp,xy_qcross,xy_hpara,xy_hperp,xy_hcross,xy_opara,xy_operp,xy_ocross,xy_epres = np.genfromtxt(xy_fpath_pres,unpack=True,skip_header=n_headlines)
	xy_bx,xy_by,xy_bz,xy_bmag,xy_ex,xy_ey,xy_ez,xy_emag,xy_curx,xy_cury,xy_curz = np.genfromtxt(xy_fpath_bande,unpack=True,skip_header=n_headlines)

	# xz data:
	xz_qdens,xz_qtemp,xz_hdens,xz_htemp,xz_odens,xz_otemp,xz_edens,xz_etemp = np.genfromtxt(xz_fpath_plas,unpack=True,skip_header=n_headlines)
	xz_qvx,xz_qvy,xz_qvz,xz_ovx,xz_ovy,xz_ovz,xz_hvx,xz_hvy,xz_hvz,xz_tvx,xz_tvy,xz_tvz,xz_alfmach = np.genfromtxt(xz_fpath_flow,unpack=True,skip_header=n_headlines)
	xz_qpara,xz_qperp,xz_qcross,xz_hpara,xz_hperp,xz_hcross,xz_opara,xz_operp,xz_ocross,xz_epres = np.genfromtxt(xz_fpath_pres,unpack=True,skip_header=n_headlines)
	xz_bx,xz_by,xz_bz,xz_bmag,xz_ex,xz_ey,xz_ez,xz_emag,xz_curx,xz_cury,xz_curz = np.genfromtxt(xz_fpath_bande,unpack=True,skip_header=n_headlines)

	# yz data:
	yz_qdens,yz_qtemp,yz_hdens,yz_htemp,yz_odens,yz_otemp,yz_edens,yz_etemp = np.genfromtxt(yz_fpath_plas,unpack=True,skip_header=n_headlines)
	yz_qvx,yz_qvy,yz_qvz,yz_ovx,yz_ovy,yz_ovz,yz_hvx,yz_hvy,yz_hvz,yz_tvx,yz_tvy,yz_tvz,yz_alfmach = np.genfromtxt(yz_fpath_flow,unpack=True,skip_header=n_headlines)
	yz_qpara,yz_qperp,yz_qcross,yz_hpara,yz_hperp,yz_hcross,yz_opara,yz_operp,yz_ocross,yz_epres = np.genfromtxt(yz_fpath_pres,unpack=True,skip_header=n_headlines)
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
	this_cbp = cbar_pos
	stream_opts = (False,'')


	#	###########################################	#
	#	Alfven speed, mag. field, and current stack	#
	#	###########################################	#
	#
	# Alfven speed:
	fig,b_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	cbar_title = r'$\frac{|\vec{v}_\mathrm{net}|}{V_\mathrm{Alf}}$'
	cbmin = 0.0
	cbmax = alfmax[box-1]
#	cbmax = xy_alfmach.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_alfmach,[ny,nx])
	hi2 = np.reshape(xz_alfmach,[nz,nx])
	hi3 = np.reshape(yz_alfmach,[nz,ny])
	if(show_vecs):
		if(diagnostic or not showstreams):
			values = (
				((xy_x,xy_y,xy_z),(xy_tvx,xy_tvy,xy_tvz)),
				((xz_x,xz_y,xz_z),(xz_tvx,xz_tvy,xz_tvz)),
				((yz_x,yz_y,yz_z),(yz_tvx,yz_tvy,yz_tvz)),
				hi1,hi2,hi3 )		
		else:
			cbar_repos = True
			stream_title = r'$|\vec{v}_\mathrm{net}|$ (km/s)'
			stream_opts = (cbar_repos,stream_title)
			# showstreams patterning (diagnostic = False only):
			vec1 = np.stack((xy_tvx,xy_tvy,xy_tvz),axis=1)
			vec2 = np.stack((xz_tvx,xz_tvy,xz_tvz),axis=1)
			vec3 = np.stack((yz_tvx,yz_tvy,yz_tvz),axis=1)
			vec1 = np.reshape(vec1,[ny,nx,3])
			vec2 = np.reshape(vec2,[nz,nx,3])
			vec3 = np.reshape(vec3,[nz,ny,3])
			values = ( vec1,vec2,vec3, hi1,hi2,hi3 )
		plot_title = run_name + r' $\vec{v}_\mathrm{net}$ ' + vecs_title + r' and Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + ut_str
	else:
		plot_title = run_name + r' Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + ut_str
		values = (hi1,hi2,hi3)

	gfx.gen_plot(diagnostic,fig,b_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,show_vecs,showstreams,box_draw,stream_opts)
	if(not stack):
		gfx.save_fig(fig,figpath_alf,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
		axes = axes[0]
	else:
		axes = b_axes[1]
		this_cbp = cbar_pos32

	# magnetic field:
	plot_title = run_name + r' net magnetic field, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$|\vec{B}_\mathrm{net}|$ (nT)'
	cbmin = 0.0
	cbmax = bmags[box-1]
#	cbmax = xy_bmag.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_bmag,[ny,nx])
	hi2 = np.reshape(xz_bmag,[nz,nx])
	hi3 = np.reshape(yz_bmag,[nz,ny])
	if(diagnostic or not showstreams):
		values = (
			((xy_x,xy_y,xy_z),(xy_bx,xy_by,xy_bz)),
			((xz_x,xz_y,xz_z),(xz_bx,xz_by,xz_bz)),
			((yz_x,yz_y,yz_z),(yz_bx,yz_by,yz_bz)),
			hi1,hi2,hi3 )		
	else:	# Problems with many zeros in B field arrays are preventing the script finishing. MJS europai run 5/11/18
		cbar_repos = False
		stream_title = ''
		stream_opts = (cbar_repos,stream_title)
		# showstreams patterning (diagnostic = False only):
		vec1 = np.stack((xy_bx,xy_by,xy_bz),axis=1)
		vec2 = np.stack((xz_bx,xz_by,xz_bz),axis=1)
		vec3 = np.stack((yz_bx,yz_by,yz_bz),axis=1)
		vec1 = np.reshape(vec1,[ny,nx,3])
		vec2 = np.reshape(vec2,[nz,nx,3])
		vec3 = np.reshape(vec3,[nz,ny,3])
		values = ( vec1,vec2,vec3, hi1,hi2,hi3 )
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,True,showstreams,box_draw,stream_opts)
	if(not stack):
		gfx.save_fig(fig,figpath_bnet,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
		axes = axes[0]
	else:
		axes = b_axes[2]
		this_cbp = cbar_pos33

	# net currents:
	plot_title = run_name + r' current density, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$|\vec{J}|$ (nA/m$^2$)'
	cbmin = 0.0
	cbmax = currents[box-1]
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	xy_curmag = np.sqrt( xy_curx**2 + xy_cury**2 + xy_curz**2 )
	xz_curmag = np.sqrt( xz_curx**2 + xz_cury**2 + xz_curz**2 )
	yz_curmag = np.sqrt( yz_curx**2 + yz_cury**2 + yz_curz**2 )
#	cbmax = xy_curmag.max()
#	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_curmag,[ny,nx])
	hi2 = np.reshape(xz_curmag,[nz,nx])
	hi3 = np.reshape(yz_curmag,[nz,ny])
	if(diagnostic or not showstreams):
		values = (
			((xy_x,xy_y,xy_z),(xy_curx,xy_cury,xy_curz)),
			((xz_x,xz_y,xz_z),(xz_curx,xz_cury,xz_curz)),
			((yz_x,yz_y,yz_z),(yz_curx,yz_cury,yz_curz)),
			hi1,hi2,hi3 )		
	else:
		cbar_repos = False
		stream_title = ''
		stream_opts = (cbar_repos,stream_title)
		# showstreams patterning (diagnostic = False only):
		vec1 = np.stack((xy_curx,xy_cury,xy_curz),axis=1)
		vec2 = np.stack((xz_curx,xz_cury,xz_curz),axis=1)
		vec3 = np.stack((yz_curx,yz_cury,yz_curz),axis=1)
		vec1 = np.reshape(vec1,[ny,nx,3])
		vec2 = np.reshape(vec2,[nz,nx,3])
		vec3 = np.reshape(vec3,[nz,ny,3])
		values = ( vec1,vec2,vec3, hi1,hi2,hi3 )
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,True,showstreams,box_draw,stream_opts)
	if(stack):
		gfx.save_fig(fig,figpath_bandalf,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_curr,xtn,fig_dpi,crop)


	#	######################	#
	#	species velocity stack	#
	#	######################	#
	#
	# q velocity:
	fig,fl_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	plot_title = run_name + r' species q velocity, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$|\vec{v}_q|$ (km/s)'
	cbmin = 0.0
	cbmax = speeds[box-1]
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	xy_qspd = np.sqrt( xy_qvx**2 + xy_qvy**2 + xy_qvz**2 )
	xz_qspd = np.sqrt( xz_qvx**2 + xz_qvy**2 + xz_qvz**2 )
	yz_qspd = np.sqrt( yz_qvx**2 + yz_qvy**2 + yz_qvz**2 )
#	cbmax = xy_qspd.max()
#	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_qspd,[ny,nx])
	hi2 = np.reshape(xz_qspd,[nz,nx])
	hi3 = np.reshape(yz_qspd,[nz,ny])
	if(show_vecs):
		values = (
			((xy_x,xy_y,xy_z),(xy_qvx,xy_qvy,xy_qvz)),
			((xz_x,xz_y,xz_z),(xz_qvx,xz_qvy,xz_qvz)),
			((yz_x,yz_y,yz_z),(yz_qvx,yz_qvy,yz_qvz)),
			hi1,hi2,hi3 )
	else:
		values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,fl_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,show_vecs,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_qvel,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = fl_axes[1]
		this_cbp = cbar_pos32

	# h velocity:
	plot_title = run_name + r' species h velocity, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$|\vec{v}_h|$ (km/s)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	xy_hspd = np.sqrt( xy_hvx**2 + xy_hvy**2 + xy_hvz**2 )
	xz_hspd = np.sqrt( xz_hvx**2 + xz_hvy**2 + xz_hvz**2 )
	yz_hspd = np.sqrt( yz_hvx**2 + yz_hvy**2 + yz_hvz**2 )
	hi1 = np.reshape(xy_hspd,[ny,nx])
	hi2 = np.reshape(xz_hspd,[nz,nx])
	hi3 = np.reshape(yz_hspd,[nz,ny])
	if(show_vecs):
		values = (
			((xy_x,xy_y,xy_z),(xy_hvx,xy_hvy,xy_hvz)),
			((xz_x,xz_y,xz_z),(xz_hvx,xz_hvy,xz_hvz)),
			((yz_x,yz_y,yz_z),(yz_hvx,yz_hvy,yz_hvz)),
			hi1,hi2,hi3 )
	else:
		values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,show_vecs,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_hvel,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = fl_axes[2]
		this_cbp = cbar_pos33

	# o velocity:
	plot_title = run_name + r' species o velocity, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$|\vec{v}_o|$ (km/s)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	xy_ospd = np.sqrt( xy_ovx**2 + xy_ovy**2 + xy_ovz**2 )
	xz_ospd = np.sqrt( xz_ovx**2 + xz_ovy**2 + xz_ovz**2 )
	yz_ospd = np.sqrt( yz_ovx**2 + yz_ovy**2 + yz_ovz**2 )
	hi1 = np.reshape(xy_ospd,[ny,nx])
	hi2 = np.reshape(xz_ospd,[nz,nx])
	hi3 = np.reshape(yz_ospd,[nz,ny])
	if(show_vecs):
		values = (
			((xy_x,xy_y,xy_z),(xy_ovx,xy_ovy,xy_ovz)),
			((xz_x,xz_y,xz_z),(xz_ovx,xz_ovy,xz_ovz)),
			((yz_x,yz_y,yz_z),(yz_ovx,yz_ovy,yz_ovz)),
			hi1,hi2,hi3 )
	else:
		values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,show_vecs,False,box_draw)
	if(stack):
		gfx.save_fig(fig,figpath_flows,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_ovel,xtn,fig_dpi,crop)


	#	###########################	#
	#	q pressure components stack	#
	#	###########################	#
	#
	# q para:
	fig,qpr_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	plot_title = run_name + r' q parallel pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{q||}$ (nPa)'
	cbmin = qpres_mid - press[box-1]
	cbmax = qpres_mid + press[box-1]
#	cbmin = xy_qpara.min()
#	cbmax = xy_qpara.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_qpara,[ny,nx])
	hi2 = np.reshape(xz_qpara,[nz,nx])
	hi3 = np.reshape(yz_qpara,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,qpr_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_qpara,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = qpr_axes[1]
		this_cbp = cbar_pos32

	# q perp:
	plot_title = run_name + r' q perpendicular pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{q\perp}$ (nPa)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_qperp,[ny,nx])
	hi2 = np.reshape(xz_qperp,[nz,nx])
	hi3 = np.reshape(yz_qperp,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_qperp,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = qpr_axes[2]
		this_cbp = cbar_pos33

	# q cross:
	plot_title = run_name + r' q cross pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{q\times}$ (nPa)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_qcross,[ny,nx])
	hi2 = np.reshape(xz_qcross,[nz,nx])
	hi3 = np.reshape(yz_qcross,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(stack):
		gfx.save_fig(fig,figpath_qpres,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_qcross,xtn,fig_dpi,crop)


	#	###########################	#
	#	h pressure components stack	#
	#	###########################	#
	#
	# h para:
	fig,hpr_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	plot_title = run_name + r' h parallel pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{h||}$ (nPa)'
	cbmin = hpres_mid - press[box-1]
	cbmax = hpres_mid + press[box-1]
#	cbmin = xy_hpara.min()
#	cbmax = xy_hpara.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_hpara,[ny,nx])
	hi2 = np.reshape(xz_hpara,[nz,nx])
	hi3 = np.reshape(yz_hpara,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,hpr_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_hpara,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = hpr_axes[1]
		this_cbp = cbar_pos32

	# h perp:
	plot_title = run_name + r' h perpendicular pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{h\perp}$ (nPa)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_hperp,[ny,nx])
	hi2 = np.reshape(xz_hperp,[nz,nx])
	hi3 = np.reshape(yz_hperp,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_hperp,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = hpr_axes[2]
		this_cbp = cbar_pos33

	# h cross:
	plot_title = run_name + r' h cross pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{h\times}$ (nPa)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_hcross,[ny,nx])
	hi2 = np.reshape(xz_hcross,[nz,nx])
	hi3 = np.reshape(yz_hcross,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(stack):
		gfx.save_fig(fig,figpath_hpres,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_hcross,xtn,fig_dpi,crop)


	#	###########################	#
	#	o pressure components stack	#
	#	###########################	#
	#
	# o para:
	fig,opr_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	plot_title = run_name + r' o parallel pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{o||}$ (nPa)'
	cbmin = opres_mid - press[box-1]
	cbmax = opres_mid + press[box-1]
#	cbmin = xy_opara.min()
#	cbmax = xy_opara.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_opara,[ny,nx])
	hi2 = np.reshape(xz_opara,[nz,nx])
	hi3 = np.reshape(yz_opara,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,opr_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_opara,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = opr_axes[1]
		this_cbp = cbar_pos32

	# o perp:
	plot_title = run_name + r' o perpendicular pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{o\perp}$ (nPa)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_operp,[ny,nx])
	hi2 = np.reshape(xz_operp,[nz,nx])
	hi3 = np.reshape(yz_operp,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_operp,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = opr_axes[2]
		this_cbp = cbar_pos33

	# o cross:
	plot_title = run_name + r' o cross pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{o\times}$ (nPa)'
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_ocross,[ny,nx])
	hi2 = np.reshape(xz_ocross,[nz,nx])
	hi3 = np.reshape(yz_ocross,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(stack):
		gfx.save_fig(fig,figpath_opres,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_ocross,xtn,fig_dpi,crop)


	#	######################################	#
	#	Electric field, electron T and P stack	#
	#	######################################	#
	#
	# electric field:
	plot_title = run_name + r' net electric field, box ' + str(box) + ', ut = ' + ut_str
	fig,e_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	cbar_title = r'$|\vec{E}_\mathrm{net}|$ (mV/m)'
	cbmin = 0.0
	cbmax = emags[box-1]
#	cbmax = xy_emag.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_emag,[ny,nx])
	hi2 = np.reshape(xz_emag,[nz,nx])
	hi3 = np.reshape(yz_emag,[nz,ny])
	if(diagnostic or not showstreams):
		values = (
			((xy_x,xy_y,xy_z),(xy_ex,xy_ey,xy_ez)),
			((xz_x,xz_y,xz_z),(xz_ex,xz_ey,xz_ez)),
			((yz_x,yz_y,yz_z),(yz_ex,yz_ey,yz_ez)),
			hi1,hi2,hi3 )		
	else:
		cbar_repos = False
		stream_title = ''
		stream_opts = (cbar_repos,stream_title)
		# showstreams patterning (diagnostic = False only):
		vec1 = np.stack((xy_ex,xy_ey,xy_ez),axis=1)
		vec2 = np.stack((xz_ex,xz_ey,xz_ez),axis=1)
		vec3 = np.stack((yz_ex,yz_ey,yz_ez),axis=1)
		vec1 = np.reshape(vec1,[ny,nx,3])
		vec2 = np.reshape(vec2,[nz,nx,3])
		vec3 = np.reshape(vec3,[nz,ny,3])
		values = ( vec1,vec2,vec3, hi1,hi2,hi3 )
	gfx.gen_plot(diagnostic,fig,e_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,True,showstreams,box_draw,stream_opts)
	if(not stack):
		gfx.save_fig(fig,figpath_efld,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
		axes = axes[0]
	else:
		axes = e_axes[1]
		this_cbp = cbar_pos32

	# electron temperature:
	plot_title = run_name + r' electron temperature, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$T_{e^-}$ (eV)'
	cbmin = 0.0
	cbmax = temps[box-1]/4
#	cbmax = xy_etemp.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_etemp,[ny,nx])
	hi2 = np.reshape(xz_etemp,[nz,nx])
	hi3 = np.reshape(yz_etemp,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_etemp,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,1,r_units)
		axes = axes[0]
	else:
		axes = e_axes[2]
		this_cbp = cbar_pos33

	# e press:
	plot_title = run_name + r' electron pressure, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,P_{e^-}$ (nPa)'
	cbmin = epres_mid - press[box-1]
	cbmax = epres_mid + press[box-1]
#	cbmin = xy_epres.min()
#	cbmax = xy_epres.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_epres,[ny,nx])
	hi2 = np.reshape(xz_epres,[nz,nx])
	hi3 = np.reshape(yz_epres,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(stack):
		gfx.save_fig(fig,figpath_elec,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_epres,xtn,fig_dpi,crop)


	#	#########################	#
	#	species temperature stack	#
	#	#########################	#
	#
	# q temperature:
	fig,temp_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	plot_title = run_name + r' species q temperature, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$k_BT_q$ (keV)'
	cbmin = 0.0
	cbmax = temps[box-1]*0.75
#	cbmax = xy_qtemp.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_qtemp,[ny,nx])
	hi2 = np.reshape(xz_qtemp,[nz,nx])
	hi3 = np.reshape(yz_qtemp,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,temp_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_qtemp,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = temp_axes[1]
		this_cbp = cbar_pos32

	# h temperature:
	plot_title = run_name + r' species h temperature, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$k_BT_h$ (keV)'
	cbmax = temps[box-1]
#	cbmax = xy_htemp.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_htemp,[ny,nx])
	hi2 = np.reshape(xz_htemp,[nz,nx])
	hi3 = np.reshape(yz_htemp,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_htemp,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = temp_axes[2]
		this_cbp = cbar_pos33

	# o temperature:
	plot_title = run_name + r' species o temperature, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'$k_BT_o$ (keV)'
	cbmax = temps[box-1]*2
#	cbmax = xy_otemp.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_otemp,[ny,nx])
	hi2 = np.reshape(xz_otemp,[nz,nx])
	hi3 = np.reshape(yz_otemp,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(stack):
		gfx.save_fig(fig,figpath_temps,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_otemp,xtn,fig_dpi,crop)


	#	#####################	#
	#	species density stack	#
	#	#####################	#
	#
	# q density:
	fig,dens_axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
	if(stack): this_cbp = cbar_pos31
	plot_title = run_name + r' species q density, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,n_q$ (cm$^{-3}$)'
	cbmin = qdens_mid - densities[box-1]
	cbmax = qdens_mid + densities[box-1]
#	cbmin = xy_qdens.min()
#	cbmax = xy_qdens.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_qdens,[ny,nx])
	hi2 = np.reshape(xz_qdens,[nz,nx])
	hi3 = np.reshape(yz_qdens,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,dens_axes[0],minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_qdens,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = dens_axes[1]
		this_cbp = cbar_pos32

	# h density:
	plot_title = run_name + r' species h density, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,n_h$ (cm$^{-3}$)'
	cbmin = hdens_mid - densities[box-1]
	cbmax = hdens_mid + densities[box-1]
#	cbmin = xy_hdens.min()
#	cbmax = xy_hdens.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_hdens,[ny,nx])
	hi2 = np.reshape(xz_hdens,[nz,nx])
	hi3 = np.reshape(yz_hdens,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(not stack):
		gfx.save_fig(fig,figpath_hdens,xtn,fig_dpi,crop)
		fig,axes = gfx.new_fig(fig_size,deltas,diagnostic,rows,r_units)
		axes = axes[0]
	else:
		axes = dens_axes[2]
		this_cbp = cbar_pos33

	# o density:
	plot_title = run_name + r' species o density, box ' + str(box) + ', ut = ' + ut_str
	cbar_title = r'log$\,n_h$ (cm$^{-3}$)'
	cbmin = odens_mid - densities[box-1]*1.5
	cbmax = odens_mid + densities[box-1]*1.5
#	cbmin = xy_odens.min()
#	cbmax = xy_odens.max()
	cbparams = (cbmin,cbmax,this_cbp,cbar_title)
	hi1 = np.reshape(xy_odens,[ny,nx])
	hi2 = np.reshape(xz_odens,[nz,nx])
	hi3 = np.reshape(yz_odens,[nz,ny])
	values = (hi1,hi2,hi3)
	gfx.gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,False,False,box_draw)
	if(stack):
		gfx.save_fig(fig,figpath_dens,xtn,fig_dpi,crop)
	else:
		gfx.save_fig(fig,figpath_odens,xtn,fig_dpi,crop)

	boxes_to_draw = boxes_to_draw + ((xmin,ymin,zmin,deltax,deltay,deltaz),)
	box_draw = (draw_boxes, boxes_to_draw)	# Set box_draw to indicate smaller boxes, but only do so after box 1 is done.

os.system( python_dir+"zip_data.run " + data_dir+' '+run_name+' '+nplots_str )

# Create gifs for the past 15 figures
if(update_gifs):
	if((xtn == 'png') or (xtn == 'gif')):
		gfx.upd_gifs(qty_list, run_name, nplots_str, gfx_dir, n_grids, xtn)
	else:
		print("gif conversion not supported for "+xtn+". gifs not updated.")
