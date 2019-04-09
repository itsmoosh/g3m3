!
!	This subroutine applies boundary conditions around the edges
!	of the system and at any irregular boundaries
!
subroutine bndry_moon(qrho,qpresx,qpresy,qpresz, &
	qpresxy,qpresxz,qpresyz, &
	qpx,qpy,qpz,rmassq, &
	hrho,hpresx,hpresy,hpresz, &
	hpresxy,hpresxz,hpresyz, &
	hpx,hpy,hpz,rmassh, &
	orho,opresx,opresy,opresz, &
	opresxy,opresxz,opresyz, &
	opx,opy,opz,rmasso,epres,bx,by,bz, &
	box,nx,ny,nz,n_grids,parm_srf,parm_mid,parm_zero, &
	ijsrf,numsrf,ijmid,nummid,ijzero, &
	numzero,mbndry,msrf,mmid,mzero, &
	qden_moon,hden_moon,oden_moon,cs_moon,gamma,ti_te, &
	vx_moon,vy_moon,vz_moon, &
	grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
	isotropic)

	integer box
	dimension qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
	qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
	qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
	qpresyz(nx,ny,nz,n_grids), &
	qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
	hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
	hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
	hpresxy(nx,ny,nz,n_grids),hpresxz(nx,ny,nz,n_grids), &
	hpresyz(nx,ny,nz,n_grids), &
	hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
	orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
	opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
	opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
	opresyz(nx,ny,nz,n_grids), &
	opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
	bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
	epres(nx,ny,nz,n_grids)
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)
	
	integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
	ijzero(mbndry,3,mzero)
	integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
	real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
	parm_zero(mbndry,7,mzero)
	
	tempi=cs_moon**2/gamma
	alpha_m=6.
	rvx=vx_moon
	rvy=vy_moon
	rvz=vz_moon
	
	dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
	dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
	dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
	
	!	Set atmospheric boundary points - specify magnetic field profile
	
	do n=1,nummid(box)
	
	!	Get coords of point
	
	i=ijmid(box,1,n)
	j=ijmid(box,2,n)
	k=ijmid(box,3,n)
	
	qrho(i,j,k,box)=parm_mid(box,1,n)
	hrho(i,j,k,box)=parm_mid(box,2,n)
	orho(i,j,k,box)=parm_mid(box,3,n)
	
	qpresx(i,j,k,box)=parm_mid(box,4,n)
	qpresy(i,j,k,box)=qpresx(i,j,k,box)
	qpresz(i,j,k,box)=qpresx(i,j,k,box)*0.9

	qpresxy(i,j,k,box)=0.
	qpresxz(i,j,k,box)=0.
	qpresyz(i,j,k,box)=0.
	
	hpresx(i,j,k,box)=parm_mid(box,5,n)
	hpresy(i,j,k,box)=hpresx(i,j,k,box)
	hpresz(i,j,k,box)=hpresx(i,j,k,box)*0.9
	
	hpresxy(i,j,k,box)=0.
	hpresxz(i,j,k,box)=0.
	hpresyz(i,j,k,box)=0.
	
	opresx(i,j,k,box)=parm_mid(box,6,n)
	opresy(i,j,k,box)=opresx(i,j,k,box)
	opresz(i,j,k,box)=opresx(i,j,k,box)*0.9
	
	opresxy(i,j,k,box)=0.
	opresxz(i,j,k,box)=0.
	opresyz(i,j,k,box)=0.
	
	epres(i,j,k,box)=parm_mid(box,7,n)
	
	opx(i,j,k,box)=orho(i,j,k,box)*rvx
	opy(i,j,k,box)=orho(i,j,k,box)*rvy
	opz(i,j,k,box)=orho(i,j,k,box)*rvz
	hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
	hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
	hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
	qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
	qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
	qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
	
	enddo
	
	!	Reset interior points
	
	do n=1,numzero(box)
		i=ijzero(box,1,n)
		j=ijzero(box,2,n)
		k=ijzero(box,3,n)
		
		qrho(i,j,k,box)=parm_zero(box,1,n)
		hrho(i,j,k,box)=parm_zero(box,2,n)
		orho(i,j,k,box)=parm_zero(box,3,n)
		
		qpresx(i,j,k,box)=parm_zero(box,4,n)
		qpresy(i,j,k,box)=qpresx(i,j,k,box)
		qpresz(i,j,k,box)=qpresx(i,j,k,box)*0.9
		
		qpresxy(i,j,k,box)=0.
		qpresxz(i,j,k,box)=0.
		qpresyz(i,j,k,box)=0.
		
		hpresx(i,j,k,box)=parm_zero(box,5,n)
		hpresy(i,j,k,box)=hpresx(i,j,k,box)
		hpresz(i,j,k,box)=hpresx(i,j,k,box)*0.9
		
		hpresxy(i,j,k,box)=0.
		hpresxz(i,j,k,box)=0.
		hpresyz(i,j,k,box)=0.
		
		opresx(i,j,k,box)=parm_zero(box,6,n)
		opresy(i,j,k,box)=opresx(i,j,k,box)
		opresz(i,j,k,box)=opresx(i,j,k,box)*0.9
		
		opresxy(i,j,k,box)=0.
		opresxz(i,j,k,box)=0.
		opresyz(i,j,k,box)=0.
		
		epres(i,j,k,box)=parm_zero(box,7,n)
		
		opx(i,j,k,box)=orho(i,j,k,box)*rvx
		opy(i,j,k,box)=orho(i,j,k,box)*rvy
		opz(i,j,k,box)=orho(i,j,k,box)*rvz
		hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
		hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
		hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
		qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
		qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
		qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
	
	enddo

	return
end subroutine bndry_moon
