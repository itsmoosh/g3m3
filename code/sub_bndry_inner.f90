!
!	This subroutine applies boundary conditions around the edges
!	of the system and at any irregular boundaries
!
subroutine bndry_inner( &
	qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
	hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
	orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
	epres,bx,by,bz, &
	qpresxy,qpresxz,qpresyz, &
	hpresxy,hpresxz,hpresyz, &
	opresxy,opresxz,opresyz, &
	nx,ny,nz,n_grids,parm_srf,parm_mid,parm_zero, &
	ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
	mbndry,msrf,mmid,mzero, &
	erho,epress,alpha_e,ti_te,o_conc, &
	re_equiv,r_inner,sbx_wind,spress, &
	grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
	grd_zmin,grd_zmax)
	
	integer box
	dimension &
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
	
	epres(nx,ny,nz,n_grids), &
	
	bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids)
	
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)
	
	integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
	ijzero(mbndry,3,mzero)
	integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
	real  parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
	parm_zero(mbndry,7,mzero)
	
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	
	common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
	sin_tilt,cos_tilt,b0
	
	!	Set surface conditions around earth
	!	with interior temperature as constant
	
	!	write(grdpt_f,*)'inside bndry inner with',mbndry
	!	write(grdpt_f,*)  numsrf,nummid,numzero
	
	aheight=r_inner+2.
	
	box=mbndry
	dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
	dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
	dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
	r_rot2=r_rot**2
	
	!	Set ionospheric boundary points - plasma quantities defined
	
	do n=1,numsrf(box)
	
		!	Get coords of point
	
		i=ijsrf(box,1,n)
		j=ijsrf(box,2,n)
		k=ijsrf(box,3,n)
		
		!	Set rotational speeds
		
		ay=grd_ymin(box)+dy*(j-1)
		ax=grd_xmin(box)+dx*(i-1)
		rx=ax*re_equiv
		ry=ay*re_equiv
		vfrac=1.
		rvy=rx*v_rot*vfrac
		rvx=-ry*v_rot*vfrac
		rvz=0.
		
		qden=amax1(parm_srf(box,1,n),qrho(i,j,k,box))
		hden=amax1(parm_srf(box,2,n),hrho(i,j,k,box))
		oden=amax1(parm_srf(box,3,n),orho(i,j,k,box))
		
		qrho(i,j,k,box)=qden
		hrho(i,j,k,box)=hden
		orho(i,j,k,box)=oden
		
		qpresx(i,j,k,box)=parm_srf(box,4,n)
		hpresx(i,j,k,box)=parm_srf(box,5,n)
		opresx(i,j,k,box)=parm_srf(box,6,n)
		
		qpresy(i,j,k,box)=parm_srf(box,4,n)
		hpresy(i,j,k,box)=parm_srf(box,5,n)
		opresy(i,j,k,box)=parm_srf(box,6,n)
		
		qpresz(i,j,k,box)=parm_srf(box,4,n)
		hpresz(i,j,k,box)=parm_srf(box,5,n)
		opresz(i,j,k,box)=parm_srf(box,6,n)
		
		epres(i,j,k,box)=parm_srf(box,7,n)
		
		hpx(i,j,k,box)=hden*rvx
		hpy(i,j,k,box)=hden*rvy
		hpz(i,j,k,box)=hden*rvz
		opx(i,j,k,box)=oden*rvx
		opy(i,j,k,box)=oden*rvy
		opz(i,j,k,box)=oden*rvz
		qpx(i,j,k,box)=qden*rvx
		qpy(i,j,k,box)=qden*rvy
		qpz(i,j,k,box)=qden*rvz
		
		!	Isotropic conditions: set off axis pressures to zero
		
		qpresxy(i,j,k,box)=0.
		qpresxz(i,j,k,box)=0.
		qpresyz(i,j,k,box)=0.
		hpresxy(i,j,k,box)=0.
		hpresxz(i,j,k,box)=0.
		hpresyz(i,j,k,box)=0.
		opresxy(i,j,k,box)=0.
		opresxz(i,j,k,box)=0.
		opresyz(i,j,k,box)=0.
		
	enddo
	
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
		hpresx(i,j,k,box)=parm_mid(box,5,n)
		opresx(i,j,k,box)=parm_mid(box,6,n)
		
		qpresy(i,j,k,box)=parm_mid(box,4,n)
		hpresy(i,j,k,box)=parm_mid(box,5,n)
		opresy(i,j,k,box)=parm_mid(box,6,n)
		
		qpresz(i,j,k,box)=parm_mid(box,4,n)
		hpresz(i,j,k,box)=parm_mid(box,5,n)
		opresz(i,j,k,box)=parm_mid(box,6,n)
		
		epres(i,j,k,box)=parm_mid(box,7,n)
		
		ay=grd_ymin(box)+dy*(j-1)
		ax=grd_xmin(box)+dx*(i-1)
		rx=ax*re_equiv
		ry=ay*re_equiv
		vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
		vfrac=1.
		rvy=rx*v_rot*vfrac
		rvx=-ry*v_rot*vfrac
		rvz=0.
		
		opx(i,j,k,box)=orho(i,j,k,box)*rvx
		opy(i,j,k,box)=orho(i,j,k,box)*rvy
		opz(i,j,k,box)=orho(i,j,k,box)*rvz
		hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
		hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
		hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
		qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
		qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
		qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
		
		qpresxy(i,j,k,box)=0.
		qpresxz(i,j,k,box)=0.
		qpresyz(i,j,k,box)=0.
		hpresxy(i,j,k,box)=0.
		hpresxz(i,j,k,box)=0.
		hpresyz(i,j,k,box)=0.
		opresxy(i,j,k,box)=0.
		opresxz(i,j,k,box)=0.
		opresyz(i,j,k,box)=0.
		
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
		hpresx(i,j,k,box)=parm_zero(box,5,n)
		opresx(i,j,k,box)=parm_zero(box,6,n)
		
		qpresy(i,j,k,box)=parm_zero(box,4,n)
		hpresy(i,j,k,box)=parm_zero(box,5,n)
		opresy(i,j,k,box)=parm_zero(box,6,n)
		
		qpresz(i,j,k,box)=parm_zero(box,4,n)
		hpresz(i,j,k,box)=parm_zero(box,5,n)
		opresz(i,j,k,box)=parm_zero(box,6,n)
		
		epres(i,j,k,box)=parm_zero(box,7,n)
		
		ay=grd_ymin(box)+dy*(j-1)
		ax=grd_xmin(box)+dx*(i-1)
		rx=ax*re_equiv
		ry=ay*re_equiv
		vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
		vfrac=1.
		rvy=rx*v_rot*vfrac
		rvx=-ry*v_rot*vfrac
		rvz=0.
		
		opx(i,j,k,box)=orho(i,j,k,box)*rvx
		opy(i,j,k,box)=orho(i,j,k,box)*rvy
		opz(i,j,k,box)=orho(i,j,k,box)*rvz
		hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
		hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
		hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
		qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
		qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
		qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
		
		qpresxy(i,j,k,box)=0.
		qpresxz(i,j,k,box)=0.
		qpresyz(i,j,k,box)=0.
		hpresxy(i,j,k,box)=0.
		hpresxz(i,j,k,box)=0.
		hpresyz(i,j,k,box)=0.
		opresxy(i,j,k,box)=0.
		opresxz(i,j,k,box)=0.
		opresyz(i,j,k,box)=0.
		
		bx(i,j,k,box)=0.
		by(i,j,k,box)=0.
		bz(i,j,k,box)=0.
		
	enddo
	
	return
end subroutine bndry_inner
