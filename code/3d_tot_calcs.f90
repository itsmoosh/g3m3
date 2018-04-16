!
!	This file contains two subroutines:
!	totfld
!	tot_b
!
subroutine totfld(bx,bx0,btx,nx,ny,nz,n_grids,box)
    !
    !     calculates the total magnetic field from the perturbed and
    !        stationary magnetic field
    !
    integer, intent(in) :: nx, ny, nz, n_grids, box
    real, intent(in) :: bx(nx,ny,nz,n_grids), bx0(nx,ny,nz,n_grids)

	real, intent(out) :: btx(nx,ny,nz)

	integer i,j,k
    !
    !$omp parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                btx(i,j,k)=bx0(i,j,k,box)+bx(i,j,k,box)
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
    !
    !      initialize static magnetic field along entire grid
    !
	integer, intent(in) :: nx, ny, nz
	real, intent(in) :: bsx(nx,ny,nz), bsy(nx,ny,nz), bsz(nx,ny,nz)

	real, intent(out) :: btot(nx,ny,nz)

	integer i,j,k
	real atot
    !
    do k=1,nz
        do j=1,ny
            do i=1,nx
                atot=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2+bsz(i,j,k)**2)
                btot(i,j,k)=amax1(atot,1.e-5)
            enddo
        enddo
    enddo
    !
    return
end
