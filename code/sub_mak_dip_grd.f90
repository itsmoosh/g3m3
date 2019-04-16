!
!	Initializes static magnetic field along entire grid
!
subroutine mak_dip_grd(bx0,by0,bz0,nx,ny,nz,n_grids, &
	mbndry,box,ijzero,numzero,mzero,r_inner, &
	grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)

	common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
		sin_tilt,cos_tilt,b0
	
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)
	dimension bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids), &
	bz0(nx,ny,nz,n_grids)
	integer ijzero(mbndry,3,mzero),numzero(mbndry)
	integer box
	
	sin_rot=sin(rot_angle)
	cos_rot=cos(rot_angle)
	
	!	write(*,*)'In mak dip: box, r_inner, b0'
	!	write(*,*) box, r_inner, b0
	
	dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
	dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
	dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
	
	!$omp  parallel do
	do k=1,nz
		az=grd_zmin(box)+dz*(k-1)
		z1=(az-zdip)
		
		do  j=1,ny
			ay=grd_ymin(box)+dy*(j-1)
			y1=(ay-ydip)
			
			do  i=1,nx
				ax=grd_xmin(box)+dx*(i-1)
				x1=(ax-xdip)
				
				!	Determine magnetic dipole field
				!	Rotate coordinates for planet motion
				
				xr=x1*cos_rot+y1*sin_rot
				yr=-x1*sin_rot+y1*cos_rot
				
				!	Tilt space to dipole space
				
				xp=xr*cos_tilt-z1*sin_tilt
				yp=yr
				zp=xr*sin_tilt+z1*cos_tilt
				
				x2=xp**2
				y2=yp**2
				z2=zp**2
				ar=sqrt(x2+y2+z2)
				
				!	Cartesian equivalent
				
				bmag=-b0/ar**5
				dbx=-3.*bmag*xp*zp
				dby=-3.*bmag*yp*zp
				
				dbz=bmag*(x2+y2-2.*z2)
				
				!	Tilt b field back to coordinate space
				
				rbx=dbx*cos_tilt+dbz*sin_tilt
				rby=dby
				rbz=-dbx*sin_tilt+dbz*cos_tilt
				
				!	Rotate b
				
				abx=rbx*cos_rot-rby*sin_rot
				aby=rbx*sin_rot+rby*cos_rot
				abz=rbz
				
				if(ar.gt.r_inner-1.5) then
					bx0(i,j,k,box)=abx
					by0(i,j,k,box)=aby
					bz0(i,j,k,box)=abz
				else
					bx0(i,j,k,box)=0.
					by0(i,j,k,box)=0.
					bz0(i,j,k,box)=0.
				endif
				
			enddo
		enddo
	enddo
	!
	!
	!     boundary conditions
	!
	if(box.le.mbndry)then
	do n=1,numzero(box)
	!
	!        get coords of point
	!
	i=ijzero(box,1,n)
	j=ijzero(box,2,n)
	k=ijzero(box,3,n)
	!
	bx0(i,j,k,box)=0.
	by0(i,j,k,box)=0.
	bz0(i,j,k,box)=0.
	enddo
	endif
	!
	return
end subroutine mak_dip_grd
