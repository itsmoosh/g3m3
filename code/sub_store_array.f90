subroutine store_array(bx,bx0,nx,ny,nz,n_grids,box)
    !
    !     calculates the total magnetic field from the perturbed and
    !        stationary magnetic field
	! REDUNDANT AND SCHEDULED FOR TRASHING
    !
    integer box
    dimension bx(nx,ny,nz,n_grids),bx0(nx,ny,nz,n_grids)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                bx(i,j,k,box)=bx0(i,j,k,box)
            enddo
        enddo
    enddo
    !
    return
end subroutine store_array
