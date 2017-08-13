!
!	This file contains two subroutines:
!	totfld
!	tot_b
!
subroutine totfld(bx,bx0,btx,nx,ny,nz,ngrd,m)
    !
    !     calculates the total magnetic field from the perturbed and
    !        stationary magnetic field
    !
    dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd),btx(nx,ny,nz)
    !
    !$omp parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                btx(i,j,k)=bx0(i,j,k,m)+bx(i,j,k,m)
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
    dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
        btot(nx,ny,nz)
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
