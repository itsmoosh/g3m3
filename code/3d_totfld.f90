!
!	Calculates the total magnetic field from the perturbed and
!	stationary magnetic field
!
subroutine totfld(bx,bx0,btx,nx,ny,nz,n_grids,box)

	implicit none

	integer, intent(in) :: nx, ny, nz, n_grids, box
	real, intent(in) :: bx(nx,ny,nz,n_grids), bx0(nx,ny,nz,n_grids)

	real, intent(out) :: btx(nx,ny,nz)

	integer i,j,k
	
	!$omp parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				btx(i,j,k)=bx0(i,j,k,box)+bx(i,j,k,box)
			enddo
		enddo
	enddo

	return
end subroutine totfld
