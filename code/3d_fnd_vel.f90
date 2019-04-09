!
!	Converts momentum into velocity for plotting purposes.
!
subroutine fnd_vel(px,py,pz,rho,vx,vy,vz,nx,ny,nz,n_grids,box)

	implicit none

	integer, intent(in) :: nx, ny, nz, n_grids, box
	real, intent(in) :: px(nx,ny,nz,n_grids),py(nx,ny,nz,n_grids), &
		pz(nx,ny,nz,n_grids),rho(nx,ny,nz,n_grids)
	real, intent(out) :: vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)

	integer i,j,k
	real arho

	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				arho=amax1(rho(i,j,k,box),0.0001)
				vx(i,j,k)=px(i,j,k,box)/arho
				vy(i,j,k)=py(i,j,k,box)/arho
				vz(i,j,k)=pz(i,j,k,box)/arho
			enddo
		enddo
	enddo

	return
end subroutine fnd_vel
