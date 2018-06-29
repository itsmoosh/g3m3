!
!	This file contains three subroutines:
!	mak_dip_all -- DEPRECATED
!	mak_dip_grd
!	mak_dip_moon -- DEPRECATED
!
subroutine mak_dip_all(bx0,by0,bz0,nx,ny,nz,n_grids, &
    mbndry,ijzero,numzero,mzero,r_inner, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      initialize static magnetic field along entire grid
	!		DEPRECATED
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
        grd_ymin(n_grids),grd_ymax(n_grids), &
        grd_zmin(n_grids),grd_zmax(n_grids)
    dimension bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids), &
        bz0(nx,ny,nz,n_grids)
    integer ijzero(mbndry,3,mzero),numzero(mbndry)
    integer box
    !
    !
	write(*,*) "WARNING: DEPRECATED. May not work as expected."
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    !$omp  parallel do
    do box=1,n_grids
        dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
        dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
        dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
        !
        do k=1,nz
            az=grd_zmin(box)+dz*(k-1)
            z1=(az-zdip)
            !
            do  j=1,ny
                ay=grd_ymin(box)+dy*(j-1)
                y1=(ay-ydip)
                !
                do  i=1,nx
                    ax=grd_xmin(box)+dx*(i-1)
                    x1=(ax-xdip)
                    !
                    !       rotate coordinates for planet motion
                    !
                    xr=x1*cos_rot+y1*sin_rot
                    yr=-x1*sin_rot+y1*cos_rot
                    !
                    !
                    !       tilt space to dipole space
                    !
                    xp=xr*cos_tilt-z1*sin_tilt
                    yp=yr
                    zp=xr*sin_tilt+z1*cos_tilt
                    !
                    x2=xp**2
                    y2=yp**2
                    z2=zp**2
                    ar=sqrt(x2+y2+z2)
                    !
                    !        cartesian equivalent
                    !
                    bmag=-b0/ar**5
                    dbx=-3.*bmag*xp*zp
                    dby=-3.*bmag*yp*zp
    				!
                    dbz=bmag*(x2+y2-2.*z2)
                    !
                    !        tilt b field back to coordinate space
                    !
                    rbx=dbx*cos_tilt+dbz*sin_tilt
                    rby=dby
                    rbz=-dbx*sin_tilt+dbz*cos_tilt
                    !
                    !         rotate b
                    !
                    abx=rbx*cos_rot-rby*sin_rot
                    aby=rbx*sin_rot+rby*cos_rot
                    abz=rbz
                    !
                    if(ar.gt.r_inner-1.5) then
                        bx0(i,j,k,box)=abx
                        by0(i,j,k,box)=aby
                        bz0(i,j,k,box)=abz
                    else
                        bx0(i,j,k,box)=0.
                        by0(i,j,k,box)=0.
                        bz0(i,j,k,box)=0.
                    endif
                    !
                enddo
            enddo
        enddo
        !
    enddo
    !
    !     boundary conditions
    !
    do box=1,mbndry
        !$omp  parallel do
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
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine mak_dip_grd(b0,xdip,ydip,zdip,rot_angle,sin_tilt,cos_tilt,r_inner, &
	nx,ny,nz,n_grids,mbndry,box,ijzero,numzero,mzero,grid_minvals,grid_maxvals, &
	bx0,by0,bz0)
	!
	!	Initialize dipole magnetic field along entire grid
	!	Only called when tilting = .true.
	!
	implicit none
	!
	integer, intent(in) :: nx, ny, nz, n_grids, box, mbndry, mzero, &
		ijzero(mbndry,3,mzero), numzero(mbndry)
	real, intent(in) :: r_inner, grid_minvals(3,n_grids), grid_maxvals(3,n_grids), &
		rot_angle, xdip, ydip, zdip, sin_tilt, cos_tilt, b0
	
	real, intent(out) :: bx0(nx,ny,nz,n_grids), by0(nx,ny,nz,n_grids), bz0(nx,ny,nz,n_grids)

	integer i, j, k, n	! loop counters/indices
	real sin_rot, cos_rot	! sin and cos of meridian longitude
	real dx, dy, dz	! grid spacing
	real ax, ay, az	! r vector of current evaluation point
	real x1, y1, z1	! script-r location vector (position relative to dipole moment location)
	real xr, yr, zr	! time-rotated version of script-r vector
	real xp, yp, zp	! rotation around normal to the plan containing dipole axis and rotation axis (chosen to be y, making the meridion longitude match the dipole tilt direction in the xy-plane)
	real x2, y2, z2	! squares of rotated script-r
	real ar	! distance from dipole center
	real bmag	! pre-factor in front of dipole field formula, NOT magnitude of B
	real dbx, dby, dbz	! dipole field vector, relative to dipole axis as z
	real rbx, rby, rbz	! dipole field rotated back from tilting of axis
	real abx, aby, abz	! dipole field rotated back from body rotation
	
	!
	sin_rot=sin(rot_angle)
	cos_rot=cos(rot_angle)
	!
	!     write(*,*)'In mak dip: box, r_inner, b0'
	!		write(*,*) box, r_inner, b0
	!
	dx=(grid_maxvals(1,box)-grid_minvals(1,box))/(nx-1.)
	dy=(grid_maxvals(2,box)-grid_minvals(2,box))/(ny-1.)
	dz=(grid_maxvals(3,box)-grid_minvals(3,box))/(nz-1.)
	!
	!$omp  parallel do
	do k=1, nz
		az = grid_minvals(3,box) + dz*(k-1)
		z1 = ( az - zdip )
		!
		do  j=1, ny
			ay = grid_minvals(2,box) + dy*(j-1)
			y1 = ( ay - ydip )
			!
			do  i=1, nx
			ax = grid_minvals(1,box) + dx*(i-1)
			x1 = ( ax - xdip )
			!
			!	Rotate coordinates for planet motion
			!
			xr =  x1*cos_rot + y1*sin_rot
			yr = -x1*sin_rot + y1*cos_rot
			!
			!	Tilt space to dipole space
			!
			xp = xr*cos_tilt - z1*sin_tilt
			yp = yr
			zp = xr*sin_tilt + z1*cos_tilt
			!
			x2 = xp**2
			y2 = yp**2
			z2 = zp**2
			ar = sqrt( x2 + y2 + z2 )
			!
			!	Calculate dipole field in tilted space (dipole axis along z)
			!
			! b0 is mu_0 * m, where m is magnetic dipole moment
			bmag = -b0/ar**5

			dbx = -3.*bmag*xp*zp
			dby = -3.*bmag*yp*zp
			dbz = bmag*(x2+y2-2.*z2)
			!
			!	Tilt B back to coordinate space
			!
			rbx =  dbx*cos_tilt + dbz*sin_tilt
			rby =  dby
			rbz = -dbx*sin_tilt + dbz*cos_tilt
			!
			!	Rotate B for time passed
			!
			abx = rbx*cos_rot - rby*sin_rot
			aby = rbx*sin_rot + rby*cos_rot
			abz = rbz
			!
			!	Zero out B inside the body
			!	This seems redundant given the conditional for the interior boundary following the loop over grid points. Remove? MJS 6/29/18
			!
			if(ar.gt.r_inner-1.5) then
				bx0(i,j,k,box) = abx
				by0(i,j,k,box) = aby
				bz0(i,j,k,box) = abz
			else
				bx0(i,j,k,box) = 0.
				by0(i,j,k,box) = 0.
				bz0(i,j,k,box) = 0.
			endif
			!
			enddo
		enddo
	enddo
	!
	!
	!	Interior boundary conditions
	!
	if(box.le.mbndry)then
		do n=1, numzero(box)
		!
		!	Get coords of point
		!
		i = ijzero(box,1,n)
		j = ijzero(box,2,n)
		k = ijzero(box,3,n)
		!
		!	Zero B at grid points inside body interior
		!
		bx0(i,j,k,box) = 0.
		by0(i,j,k,box) = 0.
		bz0(i,j,k,box) = 0.
		enddo
	endif
	!
	return
end
!
!
!	****************************************
!
!
subroutine mak_dip_moon(bx0,by0,bz0,nx,ny,nz,n_grids, &
    mbndry,box,ijzero,numzero,mzero,r_inner, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      initialize static magnetic field along entire grid
	!		DEPRECATED
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
        grd_ymin(n_grids),grd_ymax(n_grids), &
        grd_zmin(n_grids),grd_zmax(n_grids)
    dimension bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids), &
        bz0(nx,ny,nz,n_grids)
    integer ijzero(mbndry,3,mzero),numzero(mbndry)
    integer box
    !
	write(*,*) "WARNING: DEPRECATED. May not work as expected."
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    !     write(*,*)'In mak dip moon: box, nx, ny, nz, n_grids, mbndry'
	!		write(*,*) box, nx, ny, nz, n_grids, mbndry
    !
    dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    !$omp  parallel do
    do k=1,nz
        az=grd_zmin(box)+dz*(k-1)
        z1=(az-zdip)
        !
        do  j=1,ny
            ay=grd_ymin(box)+dy*(j-1)
            y1=(ay-ydip)
            !
            do  i=1,nx
                ax=grd_xmin(box)+dx*(i-1)
                x1=(ax-xdip)
                !
                !        determine magnetic dipole field
                !
                !
                !        determine magnetic dipole field
                !
                !       rotate coordinates for planet motion
                !
                xr=x1*cos_rot+y1*sin_rot
                yr=-x1*sin_rot+y1*cos_rot
                !
                !
                !       tilt space to dipole space
                !
                xp=xr*cos_tilt-z1*sin_tilt
                yp=yr
                zp=xr*sin_tilt+z1*cos_tilt
                !
                x2=xp**2
                y2=yp**2
                z2=zp**2
                ar=sqrt(x2+y2+z2)+1.e-8
                !        write(*,*) 'x2, y2, z2, ar'
                !        write(*,*) x2, y2, z2, ar
                !
                !        cartesian equivalent
                !
                bmag=-b0/ar**5
                dbx=-3.*bmag*xp*zp
                dby=-3.*bmag*yp*zp
    
                dbz=bmag*(x2+y2-2.*z2)
                !
                !        tilt b field back to coordinate space
                !
                rbx=dbx*cos_tilt+dbz*sin_tilt
                rby=dby
                rbz=-dbx*sin_tilt+dbz*cos_tilt
                !
                !         rotate b
                !
                abx=rbx*cos_rot-rby*sin_rot
                aby=rbx*sin_rot+rby*cos_rot
                abz=rbz
                !
                if(ar.gt.r_inner-1.5) then
                    bx0(i,j,k,box)=abx
                    by0(i,j,k,box)=aby
                    bz0(i,j,k,box)=abz
                else
                    bx0(i,j,k,box)=0.
                    by0(i,j,k,box)=0.
                    bz0(i,j,k,box)=0.
                endif
                !
            enddo
        enddo
    enddo
    !
    !
    !
    return
end
