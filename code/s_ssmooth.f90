subroutine ssmooth(px,wrkpx,nx,ny,nz,ngrd,m,chipx)
    !
    !     applies straight diffusion
    !
    dimension px(nx,ny,nz,ngrd),wrkpx(nx,ny,nz,ngrd)
    !
    !     write(6,*)'ssmooth',nx,ny,nz,ngrd,m,chipx
    !
    !         step 1:   diffuse 3-d direction
    !
    !$omp  parallel do
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                wrkpx(i+1,j,k,m)+wrkpx(i-1,j,k,m) &
                + wrkpx(i,j+1,k,m)+wrkpx(i,j-1,k,m) &
                + wrkpx(i,j,k+1,m)+wrkpx(i,j,k-1,m) &
                -6.*wrkpx(i,j,k,m))
            enddo
        enddo
    enddo
    !
    !     write onto original array
    !$omp  parallel do
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                wrkpx(i,j,k,m)=px(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end
