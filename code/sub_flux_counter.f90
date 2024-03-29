!
!	Calculates the mass flux being generated by the moon:
!	particles/cm-3/s
!
subroutine flux_counter(qpx,qpy,qpz,hpx,hpy,hpz, &
	opx,opy,opz,nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso, &
	grid_spacing,qflux_in,qflux_out,hflux_in,hflux_out, &
	oflux_in,oflux_out)

	integer box
	dimension &
	qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
	hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
	opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids)
	
	!	Initialize fluxes
	
	qflux_in=0.
	qflux_out=0.
	hflux_in=0.
	hflux_out=0.
	oflux_in=0.
	oflux_out=0.
	
	!	Front and back walls
	
	i1=2
	i2=nx-1
	
	qflux1=0.
	qflux2=0.
	hflux1=0.
	hflux2=0.
	oflux1=0.
	oflux2=0.
	
	do k=1,nz
		do j=1,ny
			qflux1=qflux1+qpx(i1,j,k,box)
			hflux1=hflux1+hpx(i1,j,k,box)
			oflux1=oflux1+opx(i1,j,k,box)
			
			qflux2=qflux2+qpx(i2,j,k,box)
			hflux2=hflux2+hpx(i2,j,k,box)
			oflux2=oflux2+opx(i2,j,k,box)
		enddo
	enddo

	if(hflux1.gt.0) then
		qflux_in=qflux_in+qflux1*grid_spacing/rmassq
		hflux_in=hflux_in+hflux1*grid_spacing/rmassh
		oflux_in=oflux_in+oflux1*grid_spacing/rmasso
	else
		qflux_out=qflux_out-qflux1*grid_spacing/rmassq
		hflux_out=hflux_out-hflux1*grid_spacing/rmassh
		oflux_out=oflux_out-oflux1*grid_spacing/rmasso
	endif

	if(hflux2.lt.0) then
		qflux_in=qflux_in-qflux2*grid_spacing/rmassq
		hflux_in=hflux_in-hflux2*grid_spacing/rmassh
		oflux_in=oflux_in-oflux2*grid_spacing/rmasso
	else
		qflux_out=qflux_out+qflux2*grid_spacing/rmassq
		hflux_out=hflux_out+hflux2*grid_spacing/rmassh
		oflux_out=oflux_out+oflux2*grid_spacing/rmasso
	endif
	
	!	Left and right walls
	
	j1=2
	j2=ny-1
	
	qflux1=0.
	qflux2=0.
	hflux1=0.
	hflux2=0.
	oflux1=0.
	oflux2=0.
	
	do k=1,nz
		do i=1,nx
			qflux1=qflux1+qpy(i,j1,k,box)
			hflux1=hflux1+hpy(i,j1,k,box)
			oflux1=oflux1+opy(i,j1,k,box)
			
			qflux2=qflux2+qpy(i,j2,k,box)
			hflux2=hflux2+hpy(i,j2,k,box)
			oflux2=oflux2+opy(i,j2,k,box)
		enddo
	enddo

	if(hflux1.gt.0) then
		qflux_in=qflux_in+qflux1*grid_spacing/rmassq
		hflux_in=hflux_in+hflux1*grid_spacing/rmassh
		oflux_in=oflux_in+oflux1*grid_spacing/rmasso
	else
		qflux_out=qflux_out-qflux1*grid_spacing/rmassq
		hflux_out=hflux_out-hflux1*grid_spacing/rmassh
		oflux_out=oflux_out-oflux1*grid_spacing/rmasso
	endif

	if(hflux2.lt.0) then
		qflux_in=qflux_in-qflux2*grid_spacing/rmassq
		hflux_in=hflux_in-hflux2*grid_spacing/rmassh
		oflux_in=oflux_in-oflux2*grid_spacing/rmasso
	else
		qflux_out=qflux_out+qflux2*grid_spacing/rmassq
		hflux_out=hflux_out+hflux2*grid_spacing/rmassh
		oflux_out=oflux_out+oflux2*grid_spacing/rmasso
	endif
	
	!	Front and back walls
	
	k1=2
	k2=nz-1
	
	qflux1=0.
	qflux2=0.
	hflux1=0.
	hflux2=0.
	oflux1=0.
	oflux2=0.
	
	do j=1,ny
		do i=1,nx
			qflux1=qflux1+qpz(i,j,k1,box)
			hflux1=hflux1+hpz(i,j,k1,box)
			oflux1=oflux1+opz(i,j,k1,box)
			
			qflux2=qflux2+qpz(i,j,k2,box)
			hflux2=hflux2+hpz(i,j,k2,box)
			oflux2=oflux2+opz(i,j,k2,box)
		enddo
	enddo

	if(hflux1.gt.0) then
		qflux_in=qflux_in+qflux1*grid_spacing/rmassq
		hflux_in=hflux_in+hflux1*grid_spacing/rmassh
		oflux_in=oflux_in+oflux1*grid_spacing/rmasso
	else
		qflux_out=qflux_out-qflux1*grid_spacing/rmassq
		hflux_out=hflux_out-hflux1*grid_spacing/rmassh
		oflux_out=oflux_out-oflux1*grid_spacing/rmasso
	endif

	if(hflux2.lt.0) then
		qflux_in=qflux_in-qflux2*grid_spacing/rmassq
		hflux_in=hflux_in-hflux2*grid_spacing/rmassh
		oflux_in=oflux_in-oflux2*grid_spacing/rmasso
	else
		qflux_out=qflux_out+qflux2*grid_spacing/rmassq
		hflux_out=hflux_out+hflux2*grid_spacing/rmassh
		oflux_out=oflux_out+oflux2*grid_spacing/rmasso
	endif
	
	return
end subroutine flux_counter
