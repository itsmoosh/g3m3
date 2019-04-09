!
!	Standard Runge-Kutta time step
!
subroutine push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
	efldx,efldy,efldz,nx,ny,nz,n_grids,box,delt, &
	rx,ry,rz)

	integer box
	dimension bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
	oldbx(nx,ny,nz,n_grids),oldby(nx,ny,nz,n_grids),oldbz(nx,ny,nz,n_grids), &
	efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz)
	
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
				
				!	Induction equation
				
				bx(i,j,k,box)=oldbx(i,j,k,box)-delt*( &
					( (efldz(i,jj,k)-efldz(i,jm,k))/dyt ) &
					- ( (efldy(i,j,kk)-efldy(i,j,km))/dzt ) )
				
				by(i,j,k,box)=oldby(i,j,k,box)+delt*( &
					( (efldz(ii,j,k)-efldz(im,j,k))/dxt ) &
					- ( (efldx(i,j,kk)-efldx(i,j,km))/dzt ) )
				
				bz(i,j,k,box)=oldbz(i,j,k,box)-delt*( &
					( (efldy(ii,j,k)-efldy(im,j,k))/dxt) &
					- ( (efldx(i,jj,k)-efldx(i,jm,k))/dyt ) )
			enddo
		enddo
	enddo

	return
end subroutine push_bfld
