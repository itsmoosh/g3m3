!
!	Apply the Lapidus smoothing technique for electrons
!
subroutine lap_elec(ppres,wrkppres,vx,vy,vz, &
	nx,ny,nz,n_grids,box,chirho,chipxyz,chierg,delt)

	integer box
	dimension ppres(nx,ny,nz,n_grids),wrkppres(nx,ny,nz,n_grids), &
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
				
				wrkppres(i,j,k,box)=ppres(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(ppres(ip,j,k,box)-ppres(i,j,k,box)) &
					-uxm1*(ppres(i,j,k,box)-ppres(im,j,k,box))) &
					+delty*(uyp1*(ppres(i,jp,k,box)-ppres(i,j,k,box)) &
					-uym1*(ppres(i,j,k,box)-ppres(i,jm,k,box))) &
					+deltz*(uzp1*(ppres(i,j,kp,box)-ppres(i,j,k,box)) &
					-uzm1*(ppres(i,j,k,box)-ppres(i,j,km,box))) )
			enddo
		enddo
	enddo
	
	return
end subroutine lap_elec
