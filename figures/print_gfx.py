"""
Generates polished plots from data files, output during multifluid runs. Handles all cuts at once, plotting each cut side-by-side on consistently scaled axes.
Required arguments:
	run_name
	n_grids
	nplots (I0.3 format)
	limit
	r_inner (in units of planetary radii)
	ut
	update_gifs

Example from terminal: python3 ./figures/print_gfx.py debug 5 001 60 2.0 8.4 False
"""

# Generates and returns a new figure pane and labeled axes
def new_fig(fig_size,deltas,r_units):
	deltax,deltay,deltaz = deltas

	fig = plt.subplots(1,3,figsize=fig_size,sharex=True,sharey=True)
	fig = fig[0]
	fig.subplots_adjust(left=0.06, right=0.92, wspace=0.25)
	ax1 = plt.subplot(131)
	ax2 = plt.subplot(132)
	ax3 = plt.subplot(133)
	axes = (ax1,ax2,ax3)
	ax1.set(aspect='equal')
	ax2.set(aspect='equal')
	ax3.set(aspect='equal')
	ax1.set_title('xy',fontsize=16)
	ax2.set_title('xz',fontsize=16)
	ax3.set_title('yz',fontsize=16)
	ax1.set_xlabel('x'+r_units+r', wind $\rightarrow$')
	ax1.set_ylabel('y'+r_units)
	ax2.set_xlabel('x'+r_units+r', wind $\rightarrow$')
	ax2.set_ylabel('z'+r_units)
	ax3.set_xlabel('y'+r_units+r', wind $\bigodot$')
	ax3.set_ylabel('z'+r_units)

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

	return fig,axes
####	END def new_fig	################################################


# Prints contours onto subplots and saves the figure
def gen_plot(fig,axes,minmax,xi,yi,zi,values,vecs,streams,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath,xtn,fig_dpi):
	ax1,ax2,ax3 = axes
	xmin,xmax,ymin,ymax,zmin,zmax = minmax
	cbmin,cbmax = cblims

	if(vecs):
		hi1,hi2,hi3,vec1,vec2,vec3 = values
		if(streams):
			lin_thk1 = np.sqrt(vec1[:,:,0]**2 + vec1[:,:,1]**2)
			lin_thk1 = lin_thk1 * 5/lin_thk1.max() + 0.25
			lin_thk2 = np.sqrt(vec2[:,:,0]**2 + vec2[:,:,1]**2)
			lin_thk2 = lin_thk2 * 5/lin_thk2.max() + 0.25
			lin_thk3 = np.sqrt(vec3[:,:,0]**2 + vec3[:,:,1]**2)
			lin_thk3 = lin_thk3 * 5/lin_thk3.max() + 0.25
			strm1=ax1.streamplot(xi,yi,vec1[:,:,0],vec1[:,:,1],density=0.25,linewidth=lin_thk1,arrowstyle='-|>',color=con_color)
			strm2=ax2.streamplot(xi,zi,vec2[:,:,0],vec2[:,:,1],density=0.25,linewidth=lin_thk2,arrowstyle='-|>',color=con_color)
			strm3=ax3.streamplot(yi,zi,vec3[:,:,0],vec3[:,:,1],density=0.25,linewidth=lin_thk3,arrowstyle='-|>',color=con_color)
		else:
			vec_scale = 6.e0 * np.arctan(2.e-2*(ymax-ymin)) * cbmax
			quiv1 = ax1.quiver(xi[::skip],yi[::skip],vec1[::skip,::skip,0],vec1[::skip,::skip,1],scale=vec_scale,headwidth=5,color='blue')
			quiv2 = ax2.quiver(xi[::skip],zi[::skip],vec2[::skip,::skip,0],vec2[::skip,::skip,1],scale=vec_scale,headwidth=5,color='blue')
			quiv3 = ax3.quiver(yi[::skip],zi[::skip],vec3[::skip,::skip,0],vec3[::skip,::skip,1],scale=vec_scale,headwidth=5,color='blue')
	else:
		hi1,hi2,hi3 = values
		contour_levels = np.arange(cbmin,cbmax,cbmax/n_contours)
		cont1 = ax1.contour(xi,yi,hi1,colors=con_color,linewidths=0.3,levels=contour_levels)
		cont2 = ax2.contour(xi,zi,hi2,colors=con_color,linewidths=0.3,levels=contour_levels)
		cont3 = ax3.contour(yi,zi,hi3,colors=con_color,linewidths=0.3,levels=contour_levels)

	# For some reason, these filled contours flip the vertical axis, so flip it first:
	hi1 = np.flipud(hi1)
	hi2 = np.flipud(hi2)
	hi3 = np.flipud(hi3)
	conxy = ax1.imshow(hi1,vmin=cbmin,vmax=cbmax,extent=[xmin,xmax,ymin,ymax],cmap=colormap,interpolation='bicubic')
	conxz = ax2.imshow(hi2,vmin=cbmin,vmax=cbmax,extent=[xmin,xmax,zmin,zmax],cmap=colormap,interpolation='bicubic')
	conyz = ax3.imshow(hi3,vmin=cbmin,vmax=cbmax,extent=[ymin,ymax,zmin,zmax],cmap=colormap,interpolation='bicubic')

	plt.suptitle(plot_title,fontsize=20)
	ax1.add_patch(plt.Circle((0.0,0.0), radius=r_inner, color=planet_color, zorder=10))
	ax2.add_patch(plt.Circle((0.0,0.0), radius=r_inner, color=planet_color, zorder=10))
	ax3.add_patch(plt.Circle((0.0,0.0), radius=r_inner, color=planet_color, zorder=10))

	cbar_ax = fig.add_axes(cbar_pos)
	cbar = plt.colorbar(conxy, ax=(ax1,ax2,ax3), cax=cbar_ax)
	cbar.ax.set_title(cbar_title,size=14, pad=10)
	fig.savefig(figpath, format=xtn, dpi=fig_dpi)
	plt.close(fig)
####	END def gen_plot	############################################


# Combines the 10 most recently printed plots into an animated gif
def upd_gifs(qty_list, gfx_dir, run_name, n_grids, nplots_str, xtn):
	nplots = int(nplots_str)
	gif_cmd = "convert -loop 0 -delay 100 "
	n_skipped = 10
	if(nplots < n_skipped):
		print("ERROR: update_gifs set but too few plots are present.")
		quit()
	for qty in qty_list:
		for box in range(1,n_grids+1):
			fname_start = gfx_dir + run_name + '_'
			fig_mid = "gfx_" + qty + str(box) + '_t'
			prev_plots = str(nplots - n_skipped).zfill(3)
			nplots_str = str(nplots).zfill(3)
			rangestr = 't' + prev_plots + '-' + nplots_str

			flist = ()
			for tnum in range(nplots - n_skipped, nplots, 1):
				fig_path = fname_start + fig_mid + str(tnum+1).zfill(3) + '.' + xtn
				flist = flist + (fig_path,)
			gif_name = fname_start + "anim_" + qty + str(box) + '_' + rangestr + '.gif'
			make_gifs_cmd = gif_cmd + ' '.join(flist) + ' ' + gif_name
#			os.system(make_gifs_cmd)
	print("Updated gifs for " + ' '.join(qty_list) + ', ' + rangestr )
####	END def upd_gifs	############################################


########################################################################
########################################################################
########################################################################
#						PROGRAM STARTS HERE
########################################################################
########################################################################
########################################################################

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import ticker
import scipy.interpolate
import glob
import sys
import os

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
cbar_pos = [0.94,0.1,0.02,0.75]	# Bottom-left corner x,y, w,h for colorbar, in % of fig_size
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
	update_gifs = True
else:
	update_gifs = False

fname_start = run_name + "_"
n_grids = int(n_grids_str)
fname_end = '_t' + nplots_str + '.dat'
fig_end = '_t' + nplots_str + '.' + xtn


for box in range(1,n_grids+1):
	# Construct figure file paths
	gfxp1 = gfx_dir + fname_start + 'gfx_'
	gfxp2 = str(box) + fig_end
	qty_list = ('alfmach', 'qvel', 'hvel', 'ovel', 'bmag')
	figpath_alf = gfxp1 + 'alfmach' + gfxp2
	figpath_qvel = gfxp1 + 'qvel' + gfxp2
	figpath_hvel = gfxp1 + 'hvel' + gfxp2
	figpath_ovel = gfxp1 + 'ovel' + gfxp2
	figpath_bmag = gfxp1 + 'bmag' + gfxp2

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
	minmax = (xmin,xmax,ymin,ymax,zmin,zmax)
	deltas = (deltax,deltay,deltaz)


#	# Alfven speed:
#	fig,axes = new_fig(fig_size,deltas,r_units)
#	cbar_title = r'$\frac{|\vec{v}_\mathrm{net}|}{V_\mathrm{Alf}}$'
#	cbmin = 0.0
#	cbmax = alfmax[box-1]
#	cblims = (cbmin,cbmax)
#	hi1 = np.reshape(xy_alfmach,[ny,nx])
#	hi2 = np.reshape(xz_alfmach,[nz,nx])
#	hi3 = np.reshape(yz_alfmach,[nz,ny])
#	if(show_vecs):
#		plot_title = run_name + r' $\vec{v}_\mathrm{net}$ streamlines and Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + str(ut)
#		vec1 = np.stack((xy_tvx,xy_tvy),axis=1)
#		vec2 = np.stack((xz_tvx,xz_tvz),axis=1)
#		vec3 = np.stack((yz_tvy,yz_tvz),axis=1)
#		vec1 = np.reshape(vec1,[ny,nx,2])
#		vec2 = np.reshape(vec2,[nz,nx,2])
#		vec3 = np.reshape(vec3,[nz,ny,2])
#		values = (hi1,hi2,hi3,vec1,vec2,vec3)
#		gen_plot(fig,axes,minmax,xi,yi,zi,values,True,True,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_alf,xtn,fig_dpi)
#	else:
#		plot_title = run_name + r' Alfv$\'{e}$n Mach, box ' + str(box) + ', ut = ' + str(ut)
#		values = (hi1,hi2,hi3)
#		gen_plot(fig,axes,minmax,xi,yi,zi,values,False,False,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_alf,xtn,fig_dpi)


	# q velocity:
	fig,axes = new_fig(fig_size,deltas,r_units)
	plot_title = run_name + r' species q velocity, box ' + str(box) + ', ut = ' + str(ut)
	cbar_title = r'$|\vec{v}_q|$ (km/s)'
	xy_qspd = np.sqrt( xy_qvx**2 + xy_qvy**2 + xy_qvz**2 )
	xz_qspd = np.sqrt( xz_qvx**2 + xz_qvy**2 + xz_qvz**2 )
	yz_qspd = np.sqrt( yz_qvx**2 + yz_qvy**2 + yz_qvz**2 )
	cbmin = xz_qspd.min()
	cbmax = xz_qspd.max()
	cblims = (cbmin,cbmax)
	hi1 = np.reshape(xy_qspd,[ny,nx])
	hi2 = np.reshape(xz_qspd,[nz,nx])
	hi3 = np.reshape(yz_qspd,[nz,ny])
	if(show_vecs):
		vec1 = np.stack((xy_qvx,xy_qvy),axis=1)
		vec2 = np.stack((xz_qvx,xz_qvz),axis=1)
		vec3 = np.stack((yz_qvy,yz_qvz),axis=1)
		vec1 = np.reshape(vec1,[ny,nx,2])
		vec2 = np.reshape(vec2,[nz,nx,2])
		vec3 = np.reshape(vec3,[nz,ny,2])
		values = (hi1,hi2,hi3,vec1,vec2,vec3)
		gen_plot(fig,axes,minmax,xi,yi,zi,values,True,False,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_qvel,xtn,fig_dpi)
	else:
		values = (hi1,hi2,hi3)
		gen_plot(fig,axes,minmax,xi,yi,zi,values,False,False,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_qvel,xtn,fig_dpi)


	# h velocity:
	fig,axes = new_fig(fig_size,deltas,r_units)
	plot_title = run_name + r' species h velocity, box ' + str(box) + ', ut = ' + str(ut)
	cbar_title = r'$|\vec{v}_h|$ (km/s)'
	xy_hspd = np.sqrt( xy_hvx**2 + xy_hvy**2 + xy_hvz**2 )
	xz_hspd = np.sqrt( xz_hvx**2 + xz_hvy**2 + xz_hvz**2 )
	yz_hspd = np.sqrt( yz_hvx**2 + yz_hvy**2 + yz_hvz**2 )
	cbmin = xz_hspd.min()
	cbmax = xz_hspd.max()
	cblims = (cbmin,cbmax)
	hi1 = np.reshape(xy_hspd,[ny,nx])
	hi2 = np.reshape(xz_hspd,[nz,nx])
	hi3 = np.reshape(yz_hspd,[nz,ny])
	if(show_vecs):
		vec1 = np.stack((xy_hvx,xy_hvy),axis=1)
		vec2 = np.stack((xz_hvx,xz_hvz),axis=1)
		vec3 = np.stack((yz_hvy,yz_hvz),axis=1)
		vec1 = np.reshape(vec1,[ny,nx,2])
		vec2 = np.reshape(vec2,[nz,nx,2])
		vec3 = np.reshape(vec3,[nz,ny,2])
		values = (hi1,hi2,hi3,vec1,vec2,vec3)
		gen_plot(fig,axes,minmax,xi,yi,zi,values,True,False,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_hvel,xtn,fig_dpi)
	else:
		values = (hi1,hi2,hi3)
		gen_plot(fig,axes,minmax,xi,yi,zi,values,False,False,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_hvel,xtn,fig_dpi)


	# o velocity:
	fig,axes = new_fig(fig_size,deltas,r_units)
	plot_title = run_name + r' species o velocity, box ' + str(box) + ', ut = ' + str(ut)
	cbar_title = r'$|\vec{v}_o|$ (km/s)'
	xy_ospd = np.sqrt( xy_ovx**2 + xy_ovy**2 + xy_ovz**2 )
	xz_ospd = np.sqrt( xz_ovx**2 + xz_ovy**2 + xz_ovz**2 )
	yz_ospd = np.sqrt( yz_ovx**2 + yz_ovy**2 + yz_ovz**2 )
	cbmin = xz_ospd.min()
	cbmax = xz_ospd.max()
	cblims = (cbmin,cbmax)
	hi1 = np.reshape(xy_ospd,[ny,nx])
	hi2 = np.reshape(xz_ospd,[nz,nx])
	hi3 = np.reshape(yz_ospd,[nz,ny])
	if(show_vecs):
		vec1 = np.stack((xy_ovx,xy_ovy),axis=1)
		vec2 = np.stack((xz_ovx,xz_ovz),axis=1)
		vec3 = np.stack((yz_ovy,yz_ovz),axis=1)
		vec1 = np.reshape(vec1,[ny,nx,2])
		vec2 = np.reshape(vec2,[nz,nx,2])
		vec3 = np.reshape(vec3,[nz,ny,2])
		values = (hi1,hi2,hi3,vec1,vec2,vec3)
		gen_plot(fig,axes,minmax,xi,yi,zi,values,True,False,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_ovel,xtn,fig_dpi)
	else:
		values = (hi1,hi2,hi3)
		gen_plot(fig,axes,minmax,xi,yi,zi,values,False,False,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_ovel,xtn,fig_dpi)


#	# magnetic field:
#	fig,axes = new_fig(fig_size,deltas,r_units)
#	plot_title = run_name + r' net magnetic field, box ' + str(box) + ', ut = ' + str(ut)
#	cbar_title = r'$|\vec{B}_\mathrm{net}|$ (nT)'
#	cbmin = 0.0
#	cbmax = bmags[box-1]
#	cblims = (cbmin,cbmax)
#	hi1 = np.reshape(xy_bmag,[ny,nx])
#	hi2 = np.reshape(xz_bmag,[nz,nx])
#	hi3 = np.reshape(yz_bmag,[nz,ny])
#	vec1 = np.stack((xy_bx,xy_by),axis=1)
#	vec2 = np.stack((xz_bx,xz_bz),axis=1)
#	vec3 = np.stack((yz_by,yz_bz),axis=1)
#	vec1 = np.reshape(vec1,[ny,nx,2])
#	vec2 = np.reshape(vec2,[nz,nx,2])
#	vec3 = np.reshape(vec3,[nz,ny,2])
#	values = (hi1,hi2,hi3,vec1,vec2,vec3)
#	gen_plot(fig,axes,minmax,xi,yi,zi,values,True,True,cblims,n_contours,con_color,colormap,plot_title,planet_color,r_inner,cbar_pos,cbar_title,figpath_bmag,xtn,fig_dpi)


# Create gifs for the past 10 figures
if(update_gifs):
	if(xtn == 'eps'):
		print("gif conversion not supported for eps. gifs not updated.")
	else:
		upd_gifs(qty_list, gfx_dir, run_name, n_grids, nplots_str, xtn)

t_end = os.times().elapsed
t_elapsed = (t_end - t_start) / 60
print("Debug: script duration " + str(t_elapsed) + " min")
