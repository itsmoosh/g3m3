"""
Contains plotting functions used by print_gfx.py to generate plots of physical quantities from multifluid simulations.
"""

import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.transforms import Bbox
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.patches import Circle, PathPatch
from matplotlib import ticker
import scipy.interpolate

def new_fig(fig_size,deltas,diagnostic,rows=1,r_units=' $(R_E)$'):
	"""
	Generates and returns a new figure pane and labeled axes, for one quantity.
	"""
	deltax,deltay,deltaz = deltas

	if(diagnostic):
		fig = plt.figure(figsize=fig_size)
		ax = fig.gca(projection='3d')
		ax.set_aspect('equal')	# Has no effect in matplotlib v2.2.2 due to a bug.
		gp = 0.65	# Gray % to fill in axes backgrounds on 3D plots
		ax.w_xaxis.set_pane_color((gp+0.14,gp+0.14,gp+0.14,gp+0.14))
		ax.w_zaxis.set_pane_color((gp+0.07,gp+0.07,gp+0.07,gp+0.07))
		ax.w_yaxis.set_pane_color((gp,gp,gp,gp))
		ax.set_xlabel('x'+r_units+r', wind $\rightarrow$')
		ax.set_ylabel('y'+r_units)
		ax.set_zlabel('z'+r_units)
		ax.xaxis.set_major_locator(mpl.ticker.MultipleLocator(deltax/8))
		ax.yaxis.set_major_locator(mpl.ticker.MultipleLocator(deltay/8))
		ax.zaxis.set_major_locator(mpl.ticker.MultipleLocator(deltaz/4))

		return fig,ax
	else:
		fig = plt.subplots(rows,3,figsize=fig_size,sharex=True,sharey=True)
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
		ax1.xaxis.set_major_locator(mpl.ticker.MultipleLocator(deltax/4))
		ax1.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(deltax/8))
		ax1.yaxis.set_major_locator(mpl.ticker.MultipleLocator(deltay/4))
		ax1.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(deltay/8))
		ax2.xaxis.set_major_locator(mpl.ticker.MultipleLocator(deltax/4))
		ax2.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(deltax/8))
		ax2.yaxis.set_major_locator(mpl.ticker.MultipleLocator(deltaz/4))
		ax2.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(deltaz/8))
		ax3.xaxis.set_major_locator(mpl.ticker.MultipleLocator(deltay/4))
		ax3.xaxis.set_minor_locator(mpl.ticker.MultipleLocator(deltay/8))
		ax3.yaxis.set_major_locator(mpl.ticker.MultipleLocator(deltaz/4))
		ax3.yaxis.set_minor_locator(mpl.ticker.MultipleLocator(deltaz/8))

		return fig,axes
####	END def new_fig	################################################


def gen_plot(diagnostic,fig,axes,minmax,pos,values,cbparams,plot_title,plot_opt,vecs=True,showstreams=False):
	"""
	Prints contours onto subplots, with optional vectors or streamlines.
	"""
	xi,yi,zi = pos
	xmin,xmax,ymin,ymax,zmin,zmax = minmax
	cbmin,cbmax,cbar_pos,cbar_title = cbparams
	skip,n_contours,con_color,colormap,planet_color,r_inner = plot_opt
	contour_levels = np.arange(cbmin,cbmax,cbmax/n_contours)		
	t_start = os.times().elapsed
	nx = len(xi)
	ny = len(yi)
	nz = len(zi)
	if(diagnostic):
		ax = axes
		fine_lvls = np.arange(cbmin,cbmax,cbmax/100)
		ax.set_xlim(xmin, xmax)
		ax.set_ylim(ymin, ymax)
		ax.set_zlim(zmin, 3*zmax)
		plt.title(plot_title,fontsize=20,x=0.6,y=0.76, bbox=dict(facecolor='white'))

		if(vecs):
			(pos_xy,vec1), (pos_xz,vec2), (pos_yz,vec3), hi1,hi2,hi3 = values
			skip = int((skip/1.5)**2)
#			vec_scale = 1.e-1 * np.arctan(2.e-1*(ymax-ymin)) * cbmax
#			vec1 /= vec_scale
#			vec2 /= vec_scale
#			vec3 /= vec_scale
			x,y,z = pos_xy
			vecx,vecy,vecz = vec1
			normlen= 10*(ymax-ymin)/skip
			quiv1 = ax.quiver(x[::skip],y[::skip],z[::skip],vecx[::skip],vecy[::skip],vecz[::skip],length=normlen,normalize=True,color=con_color)
			x,y,z = pos_xz
			vecx,vecy,vecz = vec2
			quiv2 = ax.quiver(x[::skip],y[::skip],z[::skip],vecx[::skip],vecy[::skip],vecz[::skip],length=normlen,normalize=True,color=con_color)
			x,y,z = pos_yz
			vecx,vecy,vecz = vec3
			quiv3 = ax.quiver(x[::skip],y[::skip],z[::skip],vecx[::skip],vecy[::skip],vecz[::skip],length=normlen,normalize=True,color=con_color)
		else:
			hi1,hi2,hi3 = values

		dbl_x = np.repeat(xi, nz)
		dbl_x = np.reshape(dbl_x,[nx,nz])
		dbl_x = np.transpose(dbl_x)
		dbl_y = np.repeat(yi, nz)
		dbl_y = np.reshape(dbl_y,[nx,nz])
		dbl_y = np.transpose(dbl_y)
		dbl_z = np.repeat(zi, nx)
		dbl_z = np.reshape(dbl_z,[nz,nx])

		# Generate contours for this plot
		conxy = ax.contourf(xi, yi, hi1, zdir='z', offset=zmin, cmap=colormap, vmin=cbmin, vmax=cbmax, levels=fine_lvls)
		conxz = ax.contourf(dbl_x, hi2, dbl_z, zdir='y', offset=ymax, cmap=colormap, vmin=cbmin, vmax=cbmax, levels=fine_lvls)
		conyz = ax.contourf(hi3, dbl_y, dbl_z, zdir='x', offset=xmin, cmap=colormap, vmin=cbmin, vmax=cbmax, levels=fine_lvls)
		# Marking the planet location looks terrible in the current implementation of pathpatch_2d_to_3d, largely because of z-fighting bugs.
#		pxy = Circle((0.0,0.0), r_inner, color=planet_color)
#		pxz = Circle((0.0,0.0), r_inner, color=planet_color)
#		pyz = Circle((0.0,0.0), r_inner, color=planet_color)
#		ax.add_patch(pxy)
#		ax.add_patch(pxz)
#		ax.add_patch(pyz)
#		art3d.pathpatch_2d_to_3d(pyz, z=xmin+1, zdir='x')
#		art3d.pathpatch_2d_to_3d(pxz, z=ymax-1, zdir='y')
#		art3d.pathpatch_2d_to_3d(pxy, z=zmin+1, zdir='z')

		# Add colorbar, crop and save figure:
		cbar_ax = fig.add_axes(cbar_pos)
		cbar = fig.colorbar(conxy, cax=cbar_ax)
		cbar_adj = 10	# Spacing in pt, meant to adjust for taller text like fractions
		cbar.ax.set_title(cbar_title,size=16, pad=cbar_adj)
		cbar.ax.tick_params(labelsize=14)
		axes = (ax,cbar_ax)
		return fig,axes

	else:
		ax1,ax2,ax3 = axes

		if(vecs):
			vec1,vec2,vec3, hi1,hi2,hi3 = values

			if(showstreams):
				lin_thk = np.sqrt(vec1[:,:,0]**2 + vec1[:,:,1]**2)
				lin_thk = lin_thk * 5/lin_thk.max() + 0.25
				strm1=ax1.streamplot(xi,yi,vec1[:,:,0],vec1[:,:,1],density=0.25,linewidth=lin_thk,arrowstyle='-|>',color=con_color)
				lin_thk = np.sqrt(vec2[:,:,0]**2 + vec2[:,:,2]**2)
				lin_thk = lin_thk * 5/lin_thk.max() + 0.25
				strm2=ax2.streamplot(xi,zi,vec2[:,:,0],vec2[:,:,2],density=0.25,linewidth=lin_thk,arrowstyle='-|>',color=con_color)
				lin_thk = np.sqrt(vec3[:,:,1]**2 + vec3[:,:,2]**2)
				lin_thk = lin_thk * 5/lin_thk.max() + 0.25
				strm3=ax3.streamplot(yi,zi,vec3[:,:,1],vec3[:,:,2],density=0.25,linewidth=lin_thk,arrowstyle='-|>',color=con_color)
			else:
				vec_scale = 6.e0 * np.arctan(2.e-2*(ymax-ymin)) * cbmax
				skip = int((skip/1.5)**2)
				(pos_xy,vec1), (pos_xz,vec2), (pos_yz,vec3) = vec1,vec2,vec3
				x,y,z = pos_xy
				vecx,vecy,vecz = vec1
				quiv1 = ax1.quiver(x[::skip],y[::skip],vecx[::skip],vecy[::skip],scale=vec_scale,headwidth=5,color=con_color)
				x,y,z = pos_xz
				vecx,vecy,vecz = vec2
				quiv2 = ax2.quiver(x[::skip],z[::skip],vecx[::skip],vecz[::skip],scale=vec_scale,headwidth=5,color=con_color)
				x,y,z = pos_yz
				vecx,vecy,vecz = vec3
				quiv3 = ax3.quiver(y[::skip],z[::skip],vecy[::skip],vecz[::skip],scale=vec_scale,headwidth=5,color=con_color)
		else:
			hi1,hi2,hi3 = values
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
		axes_row = (ax1,ax2,ax3,cbar_ax)
		return fig,axes_row
####	END def gen_plot	############################################


def save_fig(fig,figpath,xtn='png',fig_dpi=200,crop=False):
	"""
	Saves and closes a figure.
	"""
	if(crop):
		# Roughly, figure dimensions times these numbers determines the bounding box size.
		# Uncomment here to see crop box for adjusting purposes.
		#ax_crop = fig.add_axes([0.25, 0.1, 0.675, 0.625])
		crop_bbox = Bbox.from_bounds(1.9, 0.75, 8.0, 4.6)
		fig.savefig(figpath, format=xtn, dpi=fig_dpi, bbox_inches=crop_bbox)
	else:
		fig.savefig(figpath, format=xtn, dpi=fig_dpi)
	plt.close(fig)
	return
####	END def save_fig	############################################


def upd_gifs(qty_list, run_name, nplots_str, gfx_dir='./figures/images/', n_grids=5, xtn='png'):
	"""
	Combines the 10 most recently printed plots into an animated gif, then removes the still images.
	"""
	nplots = int(nplots_str)
	gif_cmd = "convert -loop 0 -delay 100 "
	n_skipped = 10
	if(nplots < n_skipped):
		print("WARNING: update_gifs set but too few plots are present. gifs not updated.")
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
				if(os.path.isfile(fig_path)):
					flist = flist + (fig_path,)
				else:
					print("Didn't find expected file: "+fig_path+", skipping.")
			gif_name = fname_start + "anim_" + qty + str(box) + '_' + rangestr + '.gif'
			img_files = ' '.join(flist)
			make_gifs_cmd = gif_cmd + img_files + ' ' + gif_name
			rm_images_cmd = "rm " + img_files
			os.system(make_gifs_cmd)
			os.system(rm_images_cmd)
	print("Updated gifs for " + ' '.join(qty_list) + ', ' + rangestr )
	return
####	END def upd_gifs	############################################


