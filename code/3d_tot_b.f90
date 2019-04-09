!
!	Initializes static magnetic field along entire grid
!
subroutine tot_b(btot,bsx,bsy,bsz,nx,ny,nz)

	implicit none

	integer, intent(in) :: nx, ny, nz
	real, intent(in) :: bsx(nx,ny,nz), bsy(nx,ny,nz), bsz(nx,ny,nz)

	real, intent(out) :: btot(nx,ny,nz)

	integer i,j,k
	real atot
	
	do k=1,nz
		do j=1,ny
			do i=1,nx
				atot=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2+bsz(i,j,k)**2)
				btot(i,j,k)=amax1(atot,1.e-5)
			enddo
		enddo
	enddo

	return
end subroutine tot_b
