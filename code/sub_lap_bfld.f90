!
!	Apply the Lapidus smoothing technique for field components
!
subroutine lap_bfld(bx,by,bz,wrkbx,wrkby,wrkbz,vx,vy,vz, &
	nx,ny,nz,n_grids,box,chirho,chipxyz,chierg,delt)

	integer box
	dimension bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
	wrkbx(nx,ny,nz,n_grids),wrkby(nx,ny,nz,n_grids), &
	wrkbz(nx,ny,nz,n_grids), &
	vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
	
	!	Equal timing irrespective of grid size
	
	deltz=delt
	delty=delt
	deltx=delt
	
	!$omp  parallel do
	do k=2,nz-1
		km=k-1
		kp=k+1

		do j=2,ny-1
			jm=j-1
			jp=j+1

			do i=2,nx-1
				im=i-1
				ip=i+1
				
				uxp1=abs(vx(ip,j,k)-vx(i,j,k))
				uxm1=abs(vx(i,j,k)-vx(im,j,k))
				
				uyp1=abs(vy(i,jp,k)-vy(i,j,k))
				uym1=abs(vy(i,j,k)-vy(i,jm,k))
				
				uzp1=abs(vz(i,j,kp)-vz(i,j,k))
				uzm1=abs(vz(i,j,k)-vz(i,j,km))
				
				wrkbx(i,j,k,box)=bx(i,j,k,box)+chierg*( &
					deltx*(uxp1*(bx(ip,j,k,box)-bx(i,j,k,box)) &
					-uxm1*(bx(i,j,k,box)-bx(im,j,k,box))) &
					+delty*(uyp1*(bx(i,jp,k,box)-bx(i,j,k,box)) &
					-uym1*(bx(i,j,k,box)-bx(i,jm,k,box))) &
					+deltz*(uzp1*(bx(i,j,kp,box)-bx(i,j,k,box)) &
					-uzm1*(bx(i,j,k,box)-bx(i,j,km,box))) )
				wrkby(i,j,k,box)=by(i,j,k,box)+chierg*( &
					deltx*(uxp1*(by(ip,j,k,box)-by(i,j,k,box)) &
					-uxm1*(by(i,j,k,box)-by(im,j,k,box))) &
					+delty*(uyp1*(by(i,jp,k,box)-by(i,j,k,box)) &
					-uym1*(by(i,j,k,box)-by(i,jm,k,box))) &
					+deltz*(uzp1*(by(i,j,kp,box)-by(i,j,k,box)) &
					-uzm1*(by(i,j,k,box)-by(i,j,km,box))) )
				wrkbz(i,j,k,box)=bz(i,j,k,box)+chierg*( &
					deltx*(uxp1*(bz(ip,j,k,box)-bz(i,j,k,box)) &
					-uxm1*(bz(i,j,k,box)-bz(im,j,k,box))) &
					+delty*(uyp1*(bz(i,jp,k,box)-bz(i,j,k,box)) &
					-uym1*(bz(i,j,k,box)-bz(i,jm,k,box))) &
					+deltz*(uzp1*(bz(i,j,kp,box)-bz(i,j,k,box)) &
					-uzm1*(bz(i,j,k,box)-bz(i,j,km,box))) )
			enddo
		enddo
	enddo
	
	return
end subroutine lap_bfld
