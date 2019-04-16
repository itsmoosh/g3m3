!
!	Checks for minimum rho and negative pressure
!	and resets value if necessary
!
subroutine set_rho(qrho,qpresx,qpresy,qpresz, &
	qpresxy,qpresxz,qpresyz,rmassq, &
	hrho,hpresx,hpresy,hpresz, &
	hpresxy,hpresxz,hpresyz,rmassh, &
	orho,opresx,opresy,opresz, &
	opresxy,opresxz,opresyz,rmasso, &
	epres,nx,ny,nz,n_grids,box,o_conc)

	integer box
	dimension qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
	qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
	qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
	qpresyz(nx,ny,nz,n_grids), &
	hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
	hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
	hpresxy(nx,ny,nz,n_grids),hpresxz(nx,ny,nz,n_grids), &
	hpresyz(nx,ny,nz,n_grids), &
	orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
	opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
	opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
	opresyz(nx,ny,nz,n_grids), &
	epres(nx,ny,nz,n_grids)
	
	!	d_min is the minimum allowable density
	
	q_min=2.e-4	! min q absolute density
	h_min=1.e-4	! min h absolute density
	o_min=2.e-5	! min o absolute density
	dstep=20.
	
	!$omp  parallel do
	do  k=1,nz
		do  j=1,ny
			do  i=1,nx
		
				tden = abs( qrho(i,j,k,box) )/rmassq &
					+ abs( hrho(i,j,k,box) ) /rmassh &
					+ abs( orho(i,j,k,box) )/rmasso + 1.e-6
				rhdens = abs( hrho(i,j,k,box) ) /rmassh/tden
				rodens = abs( orho(i,j,k,box) ) /rmasso/tden
				
				if( (qrho(i,j,k,box).lt.q_min) &
					.or. (qpresx(i,j,k,box).lt.0.) &
					.or. (qpresy(i,j,k,box).lt.0.) &
					.or. (qpresz(i,j,k,box).lt.0.) ) then
				
					qrho(i,j,k,box)=amax1( abs(qrho(i,j,k,box)),1.05*q_min )
					apres = ( abs(qpresx(i,j,k,box)) &
						+ abs(qpresy(i,j,k,box)) &
						+ abs(qpresz(i,j,k,box)) )/3.
					qpresx(i,j,k,box)=apres
					qpresy(i,j,k,box)=apres
					qpresz(i,j,k,box)=apres
					qpresxy(i,j,k,box)=0.
					qpresxz(i,j,k,box)=0.
					qpresyz(i,j,k,box)=0.
				endif
				
				if( (hrho(i,j,k,box).lt.h_min) &
					.or. (hpresx(i,j,k,box).lt.0.) &
					.or. (hpresy(i,j,k,box).lt.0.) &
					.or. (hpresz(i,j,k,box).lt.0.) ) then

					hrho(i,j,k,box)=amax1( abs(hrho(i,j,k,box)),1.05*h_min )
					apres = ( abs(hpresx(i,j,k,box)) &
						+ abs(hpresy(i,j,k,box)) &
						+ abs(hpresz(i,j,k,box)) )/3.
					hpresx(i,j,k,box)=apres
					hpresy(i,j,k,box)=apres
					hpresz(i,j,k,box)=apres
					hpresxy(i,j,k,box)=0.
					hpresxz(i,j,k,box)=0.
					hpresyz(i,j,k,box)=0.
				endif
				
				if( (orho(i,j,k,box).lt.o_min) &
					.or. (opresx(i,j,k,box).lt.0.) &
					.or. (opresy(i,j,k,box).lt.0.) &
					.or.(opresz(i,j,k,box).lt.0.) ) then

					orho(i,j,k,box) = amax1( abs(orho(i,j,k,box)), &
						1.05*o_min )
					apres = ( abs(opresx(i,j,k,box)) &
						+ abs(opresy(i,j,k,box)) &
						+ abs(opresz(i,j,k,box)) )/3.
					opresx(i,j,k,box)=apres
					opresy(i,j,k,box)=apres
					opresz(i,j,k,box)=apres
					opresxy(i,j,k,box)=0.
					opresxz(i,j,k,box)=0.
					opresyz(i,j,k,box)=0.
				endif
				
				epres(i,j,k,box)=abs(epres(i,j,k,box))
			enddo
		enddo
	enddo
	
	return
end subroutine set_rho
