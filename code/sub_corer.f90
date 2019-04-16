!
!	This subroutine fills in the core of the coarser (big) grid in
!	position nb using data from the finer (small) grid in position ns
!
subroutine corer(rho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)

	integer box
	dimension rho(nx,ny,nz,n_grids)
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)

	!	Set limits for coarse grid

	mb=box+1
	ms=box
	sx=(grd_xmax(ms)-grd_xmin(ms))/(nx-1.)
	sy=(grd_ymax(ms)-grd_ymin(ms))/(ny-1.)
	sz=(grd_zmax(ms)-grd_zmin(ms))/(nz-1.)
	dx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
	dy=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
	dz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
	
	nskip=dx/sx
	
	ai=1.+(grd_xmin(ms)-grd_xmin(mb))/dx
	nbx1=ai+.5
	aj=1.+(grd_ymin(ms)-grd_ymin(mb))/dy
	nby1=aj+.5
	ak=1.+(grd_zmin(ms)-grd_zmin(mb))/dz
	nbz1=ak+.5
	
	!$omp  parallel do
	do  ks=1+nskip,nz-nskip,nskip
		kb=nbz1+(ks-1)/nskip
		
		do  js=1+nskip,ny-nskip,nskip
			jb=nby1+(js-1)/nskip
			
			do  is=1+nskip,nx-nskip,nskip
				ib=nbx1+(is-1)/nskip
			
				!	Smooth to make compatible with coarse grid
			
				rho(ib,jb,kb,mb)=(6.*rho(is,js,ks,ms) + &
					rho(is+1,js,ks,ms)+rho(is-1,js,ks,ms) + &
					rho(is,js+1,ks,ms)+rho(is,js-1,ks,ms) + &
					rho(is,js,ks+1,ms)+rho(is,js,ks-1,ms))/12.
			enddo
		enddo
	enddo
	
	return
end subroutine corer
