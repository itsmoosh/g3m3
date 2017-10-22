subroutine ssmooth(px,wrkpx,nx,ny,nz,n_grids,box,chipx)
    !
    !     applies straight diffusion
    !
    integer box
    dimension px(nx,ny,nz,n_grids),wrkpx(nx,ny,nz,n_grids)
    !
    !     write(*,*) 'ssmooth: nx, ny, nz, n_grids, box, chipx'
	!		write(*,*) nx, ny, nz, n_grids, box, chipx
    !
    !         step 1:   diffuse 3-d direction
    !
    !$omp  parallel do
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                px(i,j,k,box)=wrkpx(i,j,k,box)+chipx*( &
                wrkpx(i+1,j,k,box)+wrkpx(i-1,j,k,box) &
                + wrkpx(i,j+1,k,box)+wrkpx(i,j-1,k,box) &
                + wrkpx(i,j,k+1,box)+wrkpx(i,j,k-1,box) &
                -6.*wrkpx(i,j,k,box))
            enddo
        enddo
    enddo
    !
    !     write onto original array
    !$omp  parallel do
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                wrkpx(i,j,k,box)=px(i,j,k,box)
            enddo
        enddo
    enddo
    !
    return
end
