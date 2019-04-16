!
!	This subroutine applies boundary conditions to all grid types
!	of the system and at any irregular boundaries
!
subroutine bndry_corer( &
	qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
	hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
	orho,opresx,opresy,opresz,opx,opy,opz, &
	epres,bx,by,bz, &
	qpresxy,qpresxz,qpresyz, &
	hpresxy,hpresxz,hpresyz, &
	opresxy,opresxz,opresyz, &
	nx,ny,nz,n_grids,box, &
	grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)

	integer box
	dimension &
	bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
	qrho(nx,ny,nz,n_grids),qpx(nx,ny,nz,n_grids), &
	qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
	qpresx(nx,ny,nz,n_grids),qpresy(nx,ny,nz,n_grids), &
	qpresz(nx,ny,nz,n_grids), qpresxy(nx,ny,nz,n_grids), &
	qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
	
	hrho(nx,ny,nz,n_grids),hpx(nx,ny,nz,n_grids), &
	hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
	hpresx(nx,ny,nz,n_grids),hpresy(nx,ny,nz,n_grids), &
	hpresz(nx,ny,nz,n_grids), hpresxy(nx,ny,nz,n_grids), &
	hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
	
	orho(nx,ny,nz,n_grids),opx(nx,ny,nz,n_grids), &
	opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
	opresx(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
	opresy(nx,ny,nz,n_grids), opresxy(nx,ny,nz,n_grids), &
	opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
	
	epres(nx,ny,nz,n_grids)
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)
	
	call corer(qrho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpresx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpresy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpresz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpresxy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpresxz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpresyz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(qpz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	
	call corer(hrho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpresx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpresy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpresz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpresxy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpresxz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpresyz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(hpz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	
	call corer(orho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opresx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opresy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opresz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opresxy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opresxz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opresyz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(opz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(epres,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(bx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(by,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	call corer(bz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
	grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	
	return
end subroutine bndry_corer
