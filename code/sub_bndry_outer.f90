!
!	This subroutine applies boundary conditions around the edges
!	of the system and at any irregular boundaries
!
subroutine bndry_outer( &
	qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
	hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
	orho,opresx,opresy,opresz,opx,opy,opz, &
	epres, &
	qpresxy,qpresxz,qpresyz, &
	hpresxy,hpresxz,hpresyz, &
	opresxy,opresxz,opresyz, &
	rmassq,rmassh,rmasso, &
	bx,by,bz,bx0,by0,bz0,nx,ny,nz,n_grids, &
	srho,rho_frac,o_conc,spress,spx,spy,spz, &
	sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)

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
	
	bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
	bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids),bz0(nx,ny,nz,n_grids)
	
	!	x boundary conditions: set oxygen in solar wind at 2.5%
	
	box=n_grids
	ofrac=rho_frac*o_conc
	frac_h=rmassh/rmassq
	frac_o=rmasso/rmassq
	
	i=1
	ii=2
	do k=1,nz
		do j=1,ny
			
			!	Wind boundary conditions at i=1
			
			qrho(i,j,k,box)=srho
			qpx(i,j,k,box)=spx
			qpy(i,j,k,box)=spy
			qpz(i,j,k,box)=spz
			qpresx(i,j,k,box)=0.5*spress
			qpresy(i,j,k,box)=qpresx(i,j,k,box)
			qpresz(i,j,k,box)=qpresx(i,j,k,box)
			
			hrho(i,j,k,box)=srho*rho_frac*frac_h
			hpx(i,j,k,box)=spx*rho_frac*frac_h
			hpy(i,j,k,box)=spy*rho_frac*frac_h
			hpz(i,j,k,box)=spz*rho_frac*frac_h
			hpresx(i,j,k,box)=0.5*spress*rho_frac
			hpresy(i,j,k,box)=hpresx(i,j,k,box)
			hpresz(i,j,k,box)=hpresx(i,j,k,box)
			
			orho(i,j,k,box)=srho*ofrac*frac_o
			opx(i,j,k,box)=spx*ofrac*frac_o
			opy(i,j,k,box)=spy*ofrac*frac_o
			opz(i,j,k,box)=spz*ofrac*frac_o
			opresx(i,j,k,box)=0.5*spress*ofrac
			opresy(i,j,k,box)=opresx(i,j,k,box)
			opresz(i,j,k,box)=opresx(i,j,k,box)
			
			epres(i,j,k,box)=qpresx(i,j,k,box)/ti_te
			
			!        iostropic conditions: set off axis pressures to zero
			
			qpresxy(i,j,k,box)=0.
			qpresxz(i,j,k,box)=0.
			qpresyz(i,j,k,box)=0.
			hpresxy(i,j,k,box)=0.
			hpresxz(i,j,k,box)=0.
			hpresyz(i,j,k,box)=0.
			opresxy(i,j,k,box)=0.
			opresxz(i,j,k,box)=0.
			opresyz(i,j,k,box)=0.
			
			!	Set bx of wind to zero, by to the wind values
			!	and subtract any geomagnetic field
			
			bx(i,j,k,box)=sbx_wind-bx0(i,j,k,box)
			by(i,j,k,box)=sby_wind-by0(i,j,k,box)
			bz(i,j,k,box)=sbz_wind-bz0(i,j,k,box)
			
		enddo
	enddo
	!
	!	x boundary conditions - back wall
	!
	nx1=nx-1
	do k=1,nz
		do j=1,ny
			
			qrho(nx,j,k,box)=qrho(nx1,j,k,box)
			qpx(nx,j,k,box)=abs(qpx(nx1,j,k,box))
			qpy(nx,j,k,box)=qpy(nx1,j,k,box)
			qpz(nx,j,k,box)=qpz(nx1,j,k,box)
			qpresx(nx,j,k,box)=qpresx(nx1,j,k,box)
			qpresy(nx,j,k,box)=qpresy(nx1,j,k,box)
			qpresz(nx,j,k,box)=qpresz(nx1,j,k,box)
			
			hrho(nx,j,k,box)=hrho(nx1,j,k,box)
			hpx(nx,j,k,box)=abs(hpx(nx1,j,k,box))
			hpy(nx,j,k,box)=hpy(nx1,j,k,box)
			hpz(nx,j,k,box)=hpz(nx1,j,k,box)
			hpresx(nx,j,k,box)=hpresx(nx1,j,k,box)
			hpresy(nx,j,k,box)=hpresy(nx1,j,k,box)
			hpresz(nx,j,k,box)=hpresz(nx1,j,k,box)
			
			orho(nx,j,k,box)=orho(nx1,j,k,box)
			opx(nx,j,k,box)=abs(opx(nx1,j,k,box))
			opy(nx,j,k,box)=opy(nx1,j,k,box)
			opz(nx,j,k,box)=opz(nx1,j,k,box)
			opresx(nx,j,k,box)=opresx(nx1,j,k,box)
			
			opresy(nx,j,k,box)=opresy(nx1,j,k,box)
			opresz(nx,j,k,box)=opresz(nx1,j,k,box)
			
			epres(nx,j,k,box)=epres(nx1,j,k,box)
			
			!	Continuous boundary condition
			
			qpresxy(nx,j,k,box)=qpresxy(nx1,j,k,box)
			qpresxz(nx,j,k,box)=qpresxz(nx1,j,k,box)
			qpresyz(nx,j,k,box)=qpresyz(nx1,j,k,box)
			hpresxy(nx,j,k,box)=hpresxy(nx1,j,k,box)
			hpresxz(nx,j,k,box)=hpresxz(nx1,j,k,box)
			hpresyz(nx,j,k,box)=hpresyz(nx1,j,k,box)
			opresxy(nx,j,k,box)=opresxy(nx1,j,k,box)
			opresxz(nx,j,k,box)=opresxz(nx1,j,k,box)
			opresyz(nx,j,k,box)=opresyz(nx1,j,k,box)
			
			bx(nx,j,k,box)=bx(nx1,j,k,box)
			by(nx,j,k,box)=by(nx1,j,k,box)
			bz(nx,j,k,box)=bz(nx1,j,k,box)
		enddo
	enddo
	
	!	y boundary conditions
	
	ny1=ny-1
	do k=1,nz
		do i=2,nx
			i1=i
			
			qrho(i,1,k,box)=qrho(i1,2,k,box)
			qpresx(i,1,k,box)=qpresx(i1,2,k,box)
			qpresy(i,1,k,box)=qpresy(i1,2,k,box)
			qpresz(i,1,k,box)=qpresz(i1,2,k,box)
			
			qpx(i,1,k,box)=qpx(i1,2,k,box)
			qpy(i,1,k,box)=-abs(qpy(i1,2,k,box))
			qpz(i,1,k,box)=qpz(i1,2,k,box)
			
			hrho(i,1,k,box)=hrho(i1,2,k,box)
			hpresx(i,1,k,box)=hpresx(i1,2,k,box)
			hpresy(i,1,k,box)=hpresy(i1,2,k,box)
			hpresz(i,1,k,box)=hpresz(i1,2,k,box)
			
			hpx(i,1,k,box)=hpx(i1,2,k,box)
			hpy(i,1,k,box)=-abs(hpy(i1,2,k,box))
			hpz(i,1,k,box)=hpz(i1,2,k,box)
			
			orho(i,1,k,box)=orho(i1,2,k,box)
			opresx(i,1,k,box)=opresx(i1,2,k,box)
			opresy(i,1,k,box)=opresy(i1,2,k,box)
			opresz(i,1,k,box)=opresz(i1,2,k,box)
			
			opx(i,1,k,box)=opx(i1,2,k,box)
			opy(i,1,k,box)=-abs(opy(i1,2,k,box))
			opz(i,1,k,box)=opz(i1,2,k,box)
			
			epres(i,1,k,box)=epres(i1,2,k,box)
			
			!	Continuous boundary condition
			
			qpresxy(i,1,k,box)=qpresxy(i1,2,k,box)
			qpresxz(i,1,k,box)=qpresxz(i1,2,k,box)
			qpresyz(i,1,k,box)=qpresyz(i1,2,k,box)
			hpresxy(i,1,k,box)=hpresxy(i1,2,k,box)
			hpresxz(i,1,k,box)=hpresxz(i1,2,k,box)
			hpresyz(i,1,k,box)=hpresyz(i1,2,k,box)
			opresxy(i,1,k,box)=opresxy(i1,2,k,box)
			opresxz(i,1,k,box)=opresxz(i1,2,k,box)
			opresyz(i,1,k,box)=opresyz(i1,2,k,box)
			
			by(i,1,k,box)=by(i1,2,k,box)
			bx(i,1,k,box)=bx(i1,2,k,box)
			bz(i,1,k,box)=bz(i1,2,k,box)
			
			qrho(i,ny,k,box)=qrho(i1,ny1,k,box)
			qpresx(i,ny,k,box)=qpresx(i1,ny1,k,box)
			qpresy(i,ny,k,box)=qpresy(i1,ny1,k,box)
			qpresz(i,ny,k,box)=qpresz(i1,ny1,k,box)
			
			qpx(i,ny,k,box)=qpx(i1,ny1,k,box)
			qpy(i,ny,k,box)=abs(qpy(i1,ny1,k,box))
			qpz(i,ny,k,box)=qpz(i1,ny1,k,box)
			
			hrho(i,ny,k,box)=hrho(i1,ny1,k,box)
			hpresx(i,ny,k,box)=hpresx(i1,ny1,k,box)
			hpresy(i,ny,k,box)=hpresy(i1,ny1,k,box)
			hpresz(i,ny,k,box)=hpresz(i1,ny1,k,box)
			hpx(i,ny,k,box)=hpx(i1,ny1,k,box)
			hpy(i,ny,k,box)=abs(hpy(i1,ny1,k,box))
			hpz(i,ny,k,box)=hpz(i1,ny1,k,box)
			
			orho(i,ny,k,box)=orho(i1,ny1,k,box)
			opresx(i,ny,k,box)=opresx(i1,ny1,k,box)
			opresy(i,ny,k,box)=opresy(i1,ny1,k,box)
			opresz(i,ny,k,box)=opresz(i1,ny1,k,box)
			
			opx(i,ny,k,box)=opx(i1,ny1,k,box)
			opy(i,ny,k,box)=abs(opy(i1,ny1,k,box))
			opz(i,ny,k,box)=opz(i1,ny1,k,box)
			
			epres(i,ny,k,box)=epres(i1,ny1,k,box)
			
			!	Continuous boundary condition
			
			qpresxy(i,ny,k,box)=qpresxy(i1,ny1,k,box)
			qpresxz(i,ny,k,box)=qpresxz(i1,ny1,k,box)
			qpresyz(i,ny,k,box)=qpresyz(i1,ny1,k,box)
			hpresxy(i,ny,k,box)=hpresxy(i1,ny1,k,box)
			hpresxz(i,ny,k,box)=hpresxz(i1,ny1,k,box)
			hpresyz(i,ny,k,box)=hpresyz(i1,ny1,k,box)
			opresxy(i,ny,k,box)=opresxy(i1,ny1,k,box)
			opresxz(i,ny,k,box)=opresxz(i1,ny1,k,box)
			opresyz(i,ny,k,box)=opresyz(i1,ny1,k,box)
			
			by(i,ny,k,box)=by(i1,ny1,k,box)
			bx(i,ny,k,box)=bx(i1,ny1,k,box)
			bz(i,ny,k,box)=bz(i1,ny1,k,box)
		enddo
	enddo
	
	!	z boundary conditions
	
	nz1=nz-1
	do j=1,ny
		do i=2,nx
			i1=i
			
			qrho(i,j,nz,box)=qrho(i1,j,nz1,box)
			qpresx(i,j,nz,box)=qpresx(i1,j,nz1,box)
			qpresy(i,j,nz,box)=qpresy(i1,j,nz1,box)
			qpresz(i,j,nz,box)=qpresz(i1,j,nz1,box)
			
			qpx(i,j,nz,box)=qpx(i1,j,nz1,box)
			qpy(i,j,nz,box)=qpy(i1,j,nz1,box)
			qpz(i,j,nz,box)=abs(qpz(i1,j,nz1,box))
			
			hrho(i,j,nz,box)=hrho(i1,j,nz1,box)
			hpresx(i,j,nz,box)=hpresx(i1,j,nz1,box)
			hpresy(i,j,nz,box)=hpresy(i1,j,nz1,box)
			hpresz(i,j,nz,box)=hpresz(i1,j,nz1,box)
			
			hpx(i,j,nz,box)=hpx(i1,j,nz1,box)
			hpy(i,j,nz,box)=hpy(i1,j,nz1,box)
			hpz(i,j,nz,box)=abs(hpz(i1,j,nz1,box))
			
			orho(i,j,nz,box)=orho(i1,j,nz1,box)
			opresx(i,j,nz,box)=opresx(i1,j,nz1,box)
			opresy(i,j,nz,box)=opresy(i1,j,nz1,box)
			opresz(i,j,nz,box)=opresz(i1,j,nz1,box)
			
			opx(i,j,nz,box)=opx(i1,j,nz1,box)
			opy(i,j,nz,box)=opy(i1,j,nz1,box)
			opz(i,j,nz,box)=abs(opz(i1,j,nz1,box))
			
			epres(i,j,nz,box)=epres(i1,j,nz1,box)
			
			!	Off-diagonal elements
			
			qpresxy(i,j,nz,box)=qpresxy(i1,j,nz1,box)
			qpresxz(i,j,nz,box)=qpresxz(i1,j,nz1,box)
			qpresyz(i,j,nz,box)=qpresyz(i1,j,nz1,box)
			hpresxy(i,j,nz,box)=hpresxy(i1,j,nz1,box)
			hpresxz(i,j,nz,box)=hpresxz(i1,j,nz1,box)
			hpresyz(i,j,nz,box)=hpresyz(i1,j,nz1,box)
			opresxy(i,j,nz,box)=opresxy(i1,j,nz1,box)
			opresxz(i,j,nz,box)=opresxz(i1,j,nz1,box)
			opresyz(i,j,nz,box)=opresyz(i1,j,nz1,box)
			
			bz(i,j,nz,box)=bz(i1,j,nz1,box)
			bx(i,j,nz,box)=bx(i1,j,nz1,box)
			by(i,j,nz,box)=by(i1,j,nz1,box)
		enddo
	enddo
	
	!	Symmetric boundary conditions at k=1 and k=2
	
	do j=1,ny
		do i=2,nx
		
		qrho(i,j,1,box)=qrho(i,j,2,box)
		qpx(i,j,1,box)=qpx(i,j,2,box)
		qpy(i,j,1,box)=qpy(i,j,2,box)
		qpz(i,j,1,box)=-abs(qpz(i,j,2,box))
		qpresx(i,j,1,box)=qpresx(i,j,2,box)
		qpresy(i,j,1,box)=qpresy(i,j,2,box)
		qpresz(i,j,1,box)=qpresz(i,j,2,box)
		
		hrho(i,j,1,box)=hrho(i,j,2,box)
		hpx(i,j,1,box)=hpx(i,j,2,box)
		hpy(i,j,1,box)=hpy(i,j,2,box)
		hpz(i,j,1,box)=-abs(hpz(i,j,2,box))
		hpresx(i,j,1,box)=hpresx(i,j,2,box)
		hpresy(i,j,1,box)=hpresy(i,j,2,box)
		hpresz(i,j,1,box)=hpresz(i,j,2,box)
		
		orho(i,j,1,box)=orho(i,j,2,box)
		opx(i,j,1,box)=opx(i,j,2,box)
		opy(i,j,1,box)=opy(i,j,2,box)
		opz(i,j,1,box)=-abs(opz(i,j,2,box))
		opresx(i,j,1,box)=opresx(i,j,2,box)
		opresy(i,j,1,box)=opresy(i,j,2,box)
		opresz(i,j,1,box)=opresz(i,j,2,box)
		
		epres(i,j,1,box)=epres(i,j,2,box)
		
		!	Off-diagonal elements
		
		qpresxy(i,j,1,box)=qpresxy(i,j,2,box)
		qpresxz(i,j,1,box)=qpresxz(i,j,2,box)
		qpresyz(i,j,1,box)=qpresyz(i,j,2,box)
		hpresxy(i,j,1,box)=hpresxy(i,j,2,box)
		hpresxz(i,j,1,box)=hpresxz(i,j,2,box)
		hpresyz(i,j,1,box)=hpresyz(i,j,2,box)
		opresxy(i,j,1,box)=opresxy(i,j,2,box)
		opresxz(i,j,1,box)=opresxz(i,j,2,box)
		opresyz(i,j,1,box)=opresyz(i,j,2,box)
		
		!	Set bx of wind to zero, by to the wind values
		!	and subtract any geomagnetic field
		!	Force bz to be positive to attain equilibrium
		
		bx(i,j,1,box)=bx(i,j,2,box)
		by(i,j,1,box)=by(i,j,2,box)
		bz(i,j,1,box)=bz(i,j,2,box)
		enddo
	enddo
	
	!	Corner lines at left-right back wall
	
	do k=2,nz-1

		qrho(nx,1,k,box)=(qrho(nx1,1,k,box)+qrho(nx,2,k,box))/2.
		qrho(nx,ny,k,box)=(qrho(nx1,ny,k,box)+qrho(nx,ny1,k,box))/2.
		hrho(nx,1,k,box)=(hrho(nx1,1,k,box)+hrho(nx,2,k,box))/2.
		hrho(nx,ny,k,box)=(hrho(nx1,ny,k,box)+hrho(nx,ny1,k,box))/2.
		orho(nx,1,k,box)=(orho(nx1,1,k,box)+orho(nx,2,k,box))/2.
		orho(nx,ny,k,box)=(orho(nx1,ny,k,box)+orho(nx,ny1,k,box))/2.
	
		qpresx(nx,1,k,box)=(qpresx(nx1,1,k,box) + qpresx(nx,2,k,box))/2.
		qpresx(nx,ny,k,box)=(qpresx(nx1,ny,k,box) + qpresx(nx,ny1,k,box))/2.
		qpresy(nx,1,k,box)=(qpresy(nx1,1,k,box) + qpresy(nx,2,k,box))/2.
		qpresy(nx,ny,k,box)=(qpresy(nx1,ny,k,box) + qpresy(nx,ny1,k,box))/2.
		qpresz(nx,1,k,box)=(qpresz(nx1,1,k,box) + qpresz(nx,2,k,box))/2.
		qpresz(nx,ny,k,box)=(qpresz(nx1,ny,k,box) + qpresz(nx,ny1,k,box))/2.
	
		hpresx(nx,1,k,box)=(hpresx(nx1,1,k,box) + hpresx(nx,2,k,box))/2.
		hpresx(nx,ny,k,box)=(hpresx(nx1,ny,k,box) + hpresx(nx,ny1,k,box))/2.
		hpresy(nx,1,k,box)=(hpresy(nx1,1,k,box) + hpresy(nx,2,k,box))/2.
		hpresy(nx,ny,k,box)=(hpresy(nx1,ny,k,box) + hpresy(nx,ny1,k,box))/2.
		hpresz(nx,1,k,box)=(hpresz(nx1,1,k,box) + hpresz(nx,2,k,box))/2.
		hpresz(nx,ny,k,box)=(hpresz(nx1,ny,k,box) + hpresz(nx,ny1,k,box))/2.
	
		opresx(nx,1,k,box)=(opresx(nx1,1,k,box) + opresx(nx,2,k,box))/2.
		opresx(nx,ny,k,box)=(opresx(nx1,ny,k,box) + opresx(nx,ny1,k,box))/2.
		opresy(nx,1,k,box)=(opresy(nx1,1,k,box) + opresy(nx,2,k,box))/2.
		opresy(nx,ny,k,box)=(opresy(nx1,ny,k,box) + opresy(nx,ny1,k,box))/2.
		opresz(nx,1,k,box)=(opresz(nx1,1,k,box) + opresz(nx,2,k,box))/2.
		opresz(nx,ny,k,box)=(opresz(nx1,ny,k,box) + opresz(nx,ny1,k,box))/2.
	
		epres(nx,1,k,box) = (epres(nx1,1,k,box) + epres(nx,2,k,box))/2.
		epres(nx,ny,k,box)=(epres(nx1,ny,k,box) + epres(nx,ny1,k,box))/2.
	
		!	Off-diagonal elements
	
		qpresxy(nx,1,k,box)=(qpresxy(nx1,1,k,box) &
		+qpresxy(nx,2,k,box))/2.
		qpresxy(nx,ny,k,box)=(qpresxy(nx1,ny,k,box) &
		+qpresxy(nx,ny1,k,box))/2.
		qpresxz(nx,1,k,box)=(qpresxz(nx1,1,k,box) &
		+qpresxz(nx,2,k,box))/2.
		qpresxz(nx,ny,k,box)=(qpresxz(nx1,ny,k,box) &
		+qpresxz(nx,ny1,k,box))/2.
		qpresyz(nx,1,k,box)=(qpresyz(nx1,1,k,box) &
		+qpresyz(nx,2,k,box))/2.
		qpresyz(nx,ny,k,box)=(qpresyz(nx1,ny,k,box) &
		+qpresyz(nx,ny1,k,box))/2.
		hpresxy(nx,1,k,box)=(hpresxy(nx1,1,k,box) &
		+hpresxy(nx,2,k,box))/2.
		hpresxy(nx,ny,k,box)=(hpresxy(nx1,ny,k,box) &
		+hpresxy(nx,ny1,k,box))/2.
		hpresxz(nx,1,k,box)=(hpresxz(nx1,1,k,box) &
		+hpresxz(nx,2,k,box))/2.
		hpresxz(nx,ny,k,box)=(hpresxz(nx1,ny,k,box) &
		+hpresxz(nx,ny1,k,box))/2.
		hpresyz(nx,1,k,box)=(hpresyz(nx1,1,k,box) &
		+hpresyz(nx,2,k,box))/2.
		hpresyz(nx,ny,k,box)=(hpresyz(nx1,ny,k,box) &
		+hpresyz(nx,ny1,k,box))/2.
		opresxy(nx,1,k,box)=(opresxy(nx1,1,k,box) &
		+opresxy(nx,2,k,box))/2.
		opresxy(nx,ny,k,box)=(opresxy(nx1,ny,k,box) &
		+opresxy(nx,ny1,k,box))/2.
		opresxz(nx,1,k,box)=(opresxz(nx1,1,k,box) &
		+opresxz(nx,2,k,box))/2.
		opresxz(nx,ny,k,box)=(opresxz(nx1,ny,k,box) &
		+opresxz(nx,ny1,k,box))/2.
		opresyz(nx,1,k,box)=(opresyz(nx1,1,k,box) &
		+opresyz(nx,2,k,box))/2.
		opresyz(nx,ny,k,box)=(opresyz(nx1,ny,k,box) &
		+opresyz(nx,ny1,k,box))/2.
	
		qpx(nx,1,k,box)=(qpx(nx1,1,k,box)+qpx(nx,2,k,box))/2.
		qpx(nx,ny,k,box)=(qpx(nx1,ny,k,box)+qpx(nx,ny1,k,box))/2.
	
		qpy(nx,1,k,box)=(qpy(nx1,1,k,box)+qpy(nx,2,k,box))/2.
		qpy(nx,ny,k,box)=(qpy(nx1,ny,k,box)+qpy(nx,ny1,k,box))/2.
	
		qpz(nx,1,k,box)=(qpz(nx1,1,k,box)+qpz(nx,2,k,box))/2.
		qpz(nx,ny,k,box)=(qpz(nx1,ny,k,box)+qpz(nx,ny1,k,box))/2.
	
		hpx(nx,1,k,box)=(hpx(nx1,1,k,box)+hpx(nx,2,k,box))/2.
		hpx(nx,ny,k,box)=(hpx(nx1,ny,k,box)+hpx(nx,ny1,k,box))/2.
	
		hpy(nx,1,k,box)=(hpy(nx1,1,k,box)+hpy(nx,2,k,box))/2.
		hpy(nx,ny,k,box)=(hpy(nx1,ny,k,box)+hpy(nx,ny1,k,box))/2.
	
		hpz(nx,1,k,box)=(hpz(nx1,1,k,box)+hpz(nx,2,k,box))/2.
		hpz(nx,ny,k,box)=(hpz(nx1,ny,k,box)+hpz(nx,ny1,k,box))/2.
	
		opx(nx,1,k,box)=(opx(nx1,1,k,box)+opx(nx,2,k,box))/2.
		opx(nx,ny,k,box)=(opx(nx1,ny,k,box)+opx(nx,ny1,k,box))/2.
	
		opy(nx,1,k,box)=(opy(nx1,1,k,box)+opy(nx,2,k,box))/2.
		opy(nx,ny,k,box)=(opy(nx1,ny,k,box)+opy(nx,ny1,k,box))/2.
	
		opz(nx,1,k,box)=(opz(nx1,1,k,box)+opz(nx,2,k,box))/2.
		opz(nx,ny,k,box)=(opz(nx1,ny,k,box)+opz(nx,ny1,k,box))/2.
	
		bx(nx,1,k,box)=bx(nx1,1,k,box)
		bx(nx,ny,k,box)=bx(nx1,ny,k,box)
	
		by(nx,1,k,box)=(by(nx1,1,k,box) + by(nx,2,k,box))/2.
		by(nx,ny,k,box)=(by(nx1,ny,k,box) + by(nx,ny1,k,box))/2.
	
		bz(nx,1,k,box)=(bz(nx1,1,k,box) + bz(nx,2,k,box))/2.
		bz(nx,ny,k,box)=(bz(nx1,ny,k,box) + bz(nx,ny1,k,box))/2.

	enddo
	
	!	Corner lines at top and bottom back wall
	
	do j=2,ny-1

		qrho(nx,j,1,box)=qrho(nx,j,2,box)
		qrho(nx,j,nz,box)=qrho(nx1,j,nz,box)
		hrho(nx,j,1,box)=hrho(nx,j,2,box)
		hrho(nx,j,nz,box)=hrho(nx1,j,nz,box)
		orho(nx,j,1,box)=orho(nx,j,2,box)
		orho(nx,j,nz,box)=orho(nx1,j,nz,box)
		
		qpresx(nx,j,1,box)=qpresx(nx,j,2,box)
		qpresx(nx,j,nz,box)=qpresx(nx1,j,nz,box)
		qpresy(nx,j,1,box)=qpresy(nx,j,2,box)
		qpresy(nx,j,nz,box)=qpresy(nx1,j,nz,box)
		qpresz(nx,j,1,box)=qpresz(nx,j,2,box)
		qpresz(nx,j,nz,box)=qpresz(nx1,j,nz,box)
		
		hpresx(nx,j,1,box)=hpresx(nx,j,2,box)
		hpresx(nx,j,nz,box)=hpresx(nx1,j,nz,box)
		hpresy(nx,j,1,box)=hpresy(nx,j,2,box)
		hpresy(nx,j,nz,box)=hpresy(nx1,j,nz,box)
		hpresz(nx,j,1,box)=hpresz(nx,j,2,box)
		hpresz(nx,j,nz,box)=hpresz(nx1,j,nz,box)
		
		opresx(nx,j,1,box)=opresx(nx,j,2,box)
		opresx(nx,j,nz,box)=opresx(nx1,j,nz,box)
		opresy(nx,j,1,box)=opresy(nx,j,2,box)
		opresy(nx,j,nz,box)=opresy(nx1,j,nz,box)
		opresz(nx,j,1,box)=opresz(nx,j,2,box)
		opresz(nx,j,nz,box)=opresz(nx1,j,nz,box)
		
		epres(nx,j,1,box)=epres(nx,j,2,box)
		
		!	Off-diagonal elements
		
		qpresxy(nx,j,1,box)=qpresxy(nx,j,2,box)
		qpresxy(nx,j,nz,box)=qpresxy(nx1,j,nz,box)
		qpresxz(nx,j,1,box)=qpresxz(nx,j,2,box)
		qpresxz(nx,j,nz,box)=qpresxz(nx1,j,nz,box)
		qpresyz(nx,j,1,box)=qpresyz(nx,j,2,box)
		qpresyz(nx,j,nz,box)=qpresyz(nx1,j,nz,box)
		hpresxy(nx,j,1,box)=hpresxy(nx,j,2,box)
		hpresxy(nx,j,nz,box)=hpresxy(nx1,j,nz,box)
		hpresxz(nx,j,1,box)=hpresxz(nx,j,2,box)
		hpresxz(nx,j,nz,box)=hpresxz(nx1,j,nz,box)
		hpresyz(nx,j,1,box)=hpresyz(nx,j,2,box)
		hpresyz(nx,j,nz,box)=hpresyz(nx1,j,nz,box)
		opresxy(nx,j,1,box)=opresxy(nx,j,2,box)
		opresxy(nx,j,nz,box)=opresxy(nx1,j,nz,box)
		opresxz(nx,j,1,box)=opresxz(nx,j,2,box)
		opresxz(nx,j,nz,box)=opresxz(nx1,j,nz,box)
		opresyz(nx,j,1,box)=opresyz(nx,j,2,box)
		opresyz(nx,j,nz,box)=opresyz(nx1,j,nz,box)
		
		qpx(nx,j,1,box)=qpx(nx,j,2,box)
		qpx(nx,j,nz,box)=qpx(nx1,j,nz,box)
		
		qpy(nx,j,1,box)=qpy(nx,j,2,box)
		qpy(nx,j,nz,box)=qpy(nx1,j,nz,box)
		
		qpz(nx,j,1,box)=qpz(nx,j,2,box)
		qpz(nx,j,nz,box)=qpz(nx1,j,nz,box)
		
		hpx(nx,j,1,box)=hpx(nx,j,2,box)
		hpx(nx,j,nz,box)=hpx(nx1,j,nz,box)
		
		hpy(nx,j,1,box)=hpy(nx,j,2,box)
		hpy(nx,j,nz,box)=hpy(nx1,j,nz,box)
		
		hpz(nx,j,1,box)=hpz(nx,j,2,box)
		hpz(nx,j,nz,box)=hpz(nx1,j,nz,box)
		
		opx(nx,j,1,box)=opx(nx,j,2,box)
		opx(nx,j,nz,box)=opx(nx1,j,nz,box)
		
		opy(nx,j,1,box)=opy(nx,j,2,box)
		opy(nx,j,nz,box)=opy(nx1,j,nz,box)
		
		opz(nx,j,1,box)=opz(nx,j,2,box)
		opz(nx,j,nz,box)=opz(nx1,j,nz,box)
		
		bx(nx,j,1,box)=bx(nx,j,2,box)
		bx(nx,j,nz,box)=bx(nx1,j,nz,box)
		
		by(nx,j,1,box)=by(nx,j,2,box)
		by(nx,j,nz,box)=by(nx1,j,nz,box)
		
		bz(nx,j,1,box)=bz(nx,j,2,box)
		bz(nx,j,nz,box)=bz(nx1,j,nz,box)
	
	enddo
	
	do i=2,nx

		qrho(i,1,1,box)=qrho(i,1,2,box)
		qrho(i,ny,1,box)=qrho(i,ny,2,box)
		qrho(i,1,nz,box)=qrho(i,1,nz1,box)
		qrho(i,ny,nz,box)=qrho(i,ny,nz1,box)
		
		hrho(i,1,1,box)=hrho(i,1,2,box)
		hrho(i,ny,1,box)=hrho(i,ny,2,box)
		hrho(i,1,nz,box)=hrho(i,1,nz1,box)
		hrho(i,ny,nz,box)=hrho(i,ny,nz1,box)
		
		orho(i,1,1,box)=orho(i,1,2,box)
		orho(i,ny,1,box)=orho(i,ny,2,box)
		orho(i,1,nz,box)=orho(i,1,nz1,box)
		orho(i,ny,nz,box)=orho(i,ny,nz1,box)
		
		qpresx(i,1,1,box)=qpresx(i,1,2,box)
		qpresx(i,ny,1,box)=qpresx(i,ny,2,box)
		qpresx(i,1,nz,box)=qpresx(i,1,nz1,box)
		qpresx(i,ny,nz,box)=qpresx(i,ny,nz1,box)
		qpresy(i,1,1,box)=qpresy(i,1,2,box)
		qpresy(i,ny,1,box)=qpresy(i,ny,2,box)
		qpresy(i,1,nz,box)=qpresy(i,1,nz1,box)
		qpresy(i,ny,nz,box)=qpresy(i,ny,nz1,box)
		qpresz(i,1,1,box)=qpresz(i,1,2,box)
		qpresz(i,ny,1,box)=qpresz(i,ny,2,box)
		qpresz(i,1,nz,box)=qpresz(i,1,nz1,box)
		qpresz(i,ny,nz,box)=qpresz(i,ny,nz1,box)
		
		hpresx(i,1,1,box)=hpresx(i,1,2,box)
		hpresx(i,ny,1,box)=hpresx(i,ny,2,box)
		hpresx(i,1,nz,box)=hpresx(i,1,nz1,box)
		hpresx(i,ny,nz,box)=hpresx(i,ny,nz1,box)
		hpresy(i,1,1,box)=hpresy(i,1,2,box)
		hpresy(i,ny,1,box)=hpresy(i,ny,2,box)
		hpresy(i,1,nz,box)=hpresy(i,1,nz1,box)
		hpresy(i,ny,nz,box)=hpresy(i,ny,nz1,box)
		hpresz(i,1,1,box)=hpresz(i,1,2,box)
		hpresz(i,ny,1,box)=hpresz(i,ny,2,box)
		hpresz(i,1,nz,box)=hpresz(i,1,nz1,box)
		hpresz(i,ny,nz,box)=hpresz(i,ny,nz1,box)
		
		opresx(i,1,1,box)=opresx(i,1,2,box)
		opresx(i,ny,1,box)=opresx(i,ny,2,box)
		opresx(i,1,nz,box)=opresx(i,1,nz1,box)
		opresx(i,ny,nz,box)=opresx(i,ny,nz1,box)
		opresy(i,1,1,box)=opresy(i,1,2,box)
		opresy(i,ny,1,box)=opresy(i,ny,2,box)
		opresy(i,1,nz,box)=opresy(i,1,nz1,box)
		opresy(i,ny,nz,box)=opresy(i,ny,nz1,box)
		opresz(i,1,1,box)=opresz(i,1,2,box)
		opresz(i,ny,1,box)=opresz(i,ny,2,box)
		opresz(i,1,nz,box)=opresz(i,1,nz1,box)
		opresz(i,ny,nz,box)=opresz(i,ny,nz1,box)
		
		epres(i,1,1,box)=epres(i,1,2,box)
		epres(i,ny,1,box)=epres(i,ny,2,box)
		epres(i,1,nz,box)=epres(i,1,nz1,box)
		epres(i,ny,nz,box)=epres(i,ny,nz1,box)
		
		!	Off-diagonal elements
		
		qpresxy(i,1,1,box)=qpresxy(i,1,2,box)
		qpresxy(i,ny,1,box)=qpresxy(i,ny,2,box)
		qpresxy(i,1,nz,box)=qpresxy(i,1,nz1,box)
		qpresxy(i,ny,nz,box)=qpresxy(i,ny,nz1,box)
		qpresxz(i,1,1,box)=qpresxz(i,1,2,box)
		qpresxz(i,ny,1,box)=qpresxz(i,ny,2,box)
		qpresxz(i,1,nz,box)=qpresxz(i,1,nz1,box)
		qpresxz(i,ny,nz,box)=qpresxz(i,ny,nz1,box)
		qpresyz(i,1,1,box)=qpresyz(i,1,2,box)
		qpresyz(i,ny,1,box)=qpresyz(i,ny,2,box)
		qpresyz(i,1,nz,box)=qpresyz(i,1,nz1,box)
		qpresyz(i,ny,nz,box)=qpresyz(i,ny,nz1,box)
		
		hpresxy(i,1,1,box)=hpresxy(i,1,2,box)
		hpresxy(i,ny,1,box)=hpresxy(i,ny,2,box)
		hpresxy(i,1,nz,box)=hpresxy(i,1,nz1,box)
		hpresxy(i,ny,nz,box)=hpresxy(i,ny,nz1,box)
		hpresxz(i,1,1,box)=hpresxz(i,1,2,box)
		hpresxz(i,ny,1,box)=hpresxz(i,ny,2,box)
		hpresxz(i,1,nz,box)=hpresxz(i,1,nz1,box)
		hpresxz(i,ny,nz,box)=hpresxz(i,ny,nz1,box)
		hpresyz(i,1,1,box)=hpresyz(i,1,2,box)
		hpresyz(i,ny,1,box)=hpresyz(i,ny,2,box)
		hpresyz(i,1,nz,box)=hpresyz(i,1,nz1,box)
		hpresyz(i,ny,nz,box)=hpresyz(i,ny,nz1,box)
		
		opresxy(i,1,1,box)=opresxy(i,1,2,box)
		opresxy(i,ny,1,box)=opresxy(i,ny,2,box)
		opresxy(i,1,nz,box)=opresxy(i,1,nz1,box)
		opresxy(i,ny,nz,box)=opresxy(i,ny,nz1,box)
		opresxz(i,1,1,box)=opresxz(i,1,2,box)
		opresxz(i,ny,1,box)=opresxz(i,ny,2,box)
		opresxz(i,1,nz,box)=opresxz(i,1,nz1,box)
		opresxz(i,ny,nz,box)=opresxz(i,ny,nz1,box)
		opresyz(i,1,1,box)=opresyz(i,1,2,box)
		opresyz(i,ny,1,box)=opresyz(i,ny,2,box)
		opresyz(i,1,nz,box)=opresyz(i,1,nz1,box)
		opresyz(i,ny,nz,box)=opresyz(i,ny,nz1,box)
		
		qpx(i,1,1,box)=qpx(i,1,2,box)
		qpx(i,ny,1,box)=qpx(i,ny,2,box)
		qpx(i,1,nz,box)=qpx(i,1,nz1,box)
		qpx(i,ny,nz,box)=qpx(i,ny,nz1,box)
		
		qpy(i,1,1,box)=qpy(i,1,2,box)
		qpy(i,ny,1,box)=qpy(i,ny,2,box)
		qpy(i,1,nz,box)=qpy(i,1,nz1,box)
		qpy(i,ny,nz,box)=qpy(i,ny,nz1,box)
		
		qpz(i,1,1,box)=qpz(i,1,2,box)
		qpz(i,ny,1,box)=qpz(i,ny,2,box)
		qpz(i,1,nz,box)=qpz(i,1,nz1,box)
		qpz(i,ny,nz,box)=qpz(i,ny,nz1,box)
		
		hpx(i,1,1,box)=hpx(i,1,2,box)
		hpx(i,ny,1,box)=hpx(i,ny,2,box)
		hpx(i,1,nz,box)=hpx(i,1,nz1,box)
		hpx(i,ny,nz,box)=hpx(i,ny,nz1,box)
		
		hpy(i,1,1,box)=hpy(i,1,2,box)
		hpy(i,ny,1,box)=hpy(i,ny,2,box)
		hpy(i,1,nz,box)=hpy(i,1,nz1,box)
		hpy(i,ny,nz,box)=hpy(i,ny,nz1,box)
		
		hpz(i,1,1,box)=hpz(i,1,2,box)
		hpz(i,ny,1,box)=hpz(i,ny,2,box)
		hpz(i,1,nz,box)=hpz(i,1,nz1,box)
		hpz(i,ny,nz,box)=hpz(i,ny,nz1,box)
		
		opx(i,1,1,box)=opx(i,1,2,box)
		opx(i,ny,1,box)=opx(i,ny,2,box)
		opx(i,1,nz,box)=opx(i,1,nz1,box)
		opx(i,ny,nz,box)=opx(i,ny,nz1,box)
		
		opy(i,1,1,box)=opy(i,1,2,box)
		opy(i,ny,1,box)=opy(i,ny,2,box)
		opy(i,1,nz,box)=opy(i,1,nz1,box)
		opy(i,ny,nz,box)=opy(i,ny,nz1,box)
		
		opz(i,1,1,box)=opz(i,1,2,box)
		opz(i,ny,1,box)=opz(i,ny,2,box)
		opz(i,1,nz,box)=opz(i,1,nz1,box)
		opz(i,ny,nz,box)=opz(i,ny,nz1,box)
		
		bx(i,1,1,box)=bx(i,1,2,box)
		bx(i,ny,1,box)=bx(i,ny,2,box)
		bx(i,1,nz,box)=bx(i,1,nz1,box)
		bx(i,ny,nz,box)=bx(i,ny,nz1,box)
		
		by(i,1,1,box)=by(i,1,2,box)
		by(i,ny,1,box)=by(i,ny,2,box)
		by(i,1,nz,box)=by(i,1,nz1,box)
		by(i,ny,nz,box)=by(i,ny,nz1,box)
		
		bz(i,1,1,box)=bz(i,1,2,box)
		bz(i,ny,1,box)=bz(i,ny,2,box)
		bz(i,1,nz,box)=bz(i,1,nz1,box)
		bz(i,ny,nz,box)=bz(i,ny,nz1,box)
	
	enddo
	
	return
end subroutine bndry_outer
