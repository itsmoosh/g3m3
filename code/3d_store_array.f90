subroutine store_array(bx,bx0,nx,ny,nz,ngrd,m)
    !
    !     calculates the total magnetic field from the perturbed and
    !        stationary magnetic field
    !
    dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                bx(i,j,k,m)=bx0(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end
