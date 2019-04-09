!
!	Evolves the electron pressure equation
!
subroutine push_elec(wrkepres,oldepres,epres,evx,evy,evz, &
	gamma,gamma1,nx,ny,nz,n_grids,box,delt,rx,ry,rz)

	integer box
	dimension epres(nx,ny,nz,n_grids),oldepres(nx,ny,nz,n_grids), &
	wrkepres(nx,ny,nz,n_grids)
	
	dimension evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
	
	!	Set physical grid spacing
	
	dxt=2.*rx
	dyt=2.*ry
	dzt=2.*rz
	
	!$omp  parallel do
	do k=2,nz-1
		km=k-1
		kk=k+1
		
		do j=2,ny-1
			jm=j-1
			jj=j+1
			
			do i=2,nx-1
				ii=i+1
				im=i-1
				egradp_x=(wrkepres(ii,j,k,box)-wrkepres(im,j,k,box))/dxt
				egradp_y=(wrkepres(i,jj,k,box)-wrkepres(i,jm,k,box))/dyt
				egradp_z=(wrkepres(i,j,kk,box)-wrkepres(i,j,km,box))/dzt
				
				!	Pressure equations
				
				epres(i,j,k,box)=oldepres(i,j,k,box)-delt*gamma* &
					( ( (wrkepres(ii,j,k,box)*evx(ii,j,k) &
					-wrkepres(im,j,k,box)*evx(im,j,k))/dxt ) + &
					( (wrkepres(i,jj,k,box)*evy(i,jj,k) &
					-wrkepres(i,jm,k,box)*evy(i,jm,k))/dyt ) + &
					( (wrkepres(i,j,kk,box)*evz(i,j,kk) &
					-wrkepres(i,j,km,box)*evz(i,j,km))/dzt ) )
				epres(i,j,k,box)=epres(i,j,k,box) &
					+ delt*gamma1*( &
					evx(i,j,k)*egradp_x+evy(i,j,k)*egradp_y &
					+evz(i,j,k)*egradp_z )
			enddo
		enddo
	enddo
	
	return
end subroutine push_elec
