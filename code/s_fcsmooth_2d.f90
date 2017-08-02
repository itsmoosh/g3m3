subroutine fcsmooth_2d(px,oldpx,wrkpx,nx,ny,nz,ngrd,m,chipx, &
    delta,delta1,fn)
    !
    !     applies flux correction smoothing to field quantities
    !           oldpx is the qunatity at the start of the time step
    !           wrkpx is the unsmoothed qunatity
    !           px is the final product
    !     wrk arrays assumed to have dimension larger than nx,ny,nz
    !
    dimension px(nx,ny,nz,ngrd),oldpx(nx,ny,nz,ngrd), &
    wrkpx(nx,ny,nz,ngrd), &
    delta(nx,ny,nz),delta1(nx,ny,nz),fn(nx,ny,nz)
    !
    !         step 1:   diffuse initial in x direction
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                oldpx(i+1,j,k,m)+oldpx(i-1,j,k,m) &
                -2.*oldpx(i,j,k,m))
            enddo
        enddo
    enddo
    !
    !     set boundary points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            px(1,j,k,m)=wrkpx(1,j,k,m)
            px(nx,j,k,m)=wrkpx(nx,j,k,m)
        enddo
    enddo
    !
    !     calculate the difference  between x nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx-1
                i1=i+1
                delta(i,j,k)=wrkpx(i1,j,k,m)-wrkpx(i,j,k,m)   !d-non-diffuse n+1 sol'n
                delta1(i,j,k)=px(i1,j,k,m)-px(i,j,k,m)   !d1-diffused n+1 sol'n
            enddo
        enddo
    enddo
    !
    do k=1,nz
        do j=1,ny
            delta(nx,j,k)=delta(nx-1,j,k)
            delta1(nx,j,k)=delta1(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                s=sign(1.,delta(i,j,k))
                fn(i,j,k)=s*amax1(0.,amin1(s*delta1(i-1,j,k), &
                chipx*abs(delta(i,j,k)),s*delta1(i+1,j,k)))
            enddo
        enddo
    enddo
    !
    do k=1,nz
        do j=1,ny
            fn(1,j,k)=fn(2,j,k)
            fn(nx,j,k)=fn(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                i1=i-1
                wrkpx(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i1,j,k) ! all smoothing
                !        px(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i1,j,k)  ! x-smooth only
            enddo
        enddo
    enddo
    !
    !         step 2:   diffuse initial in y -dreiction
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=1,nx
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                oldpx(i,jp,k,m)+oldpx(i,jm,k,m) &
                -2.*oldpx(i,j,k,m))
            enddo
        enddo
    enddo
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            px(i,1,k,m)=wrkpx(i,1,k,m)
            px(i,ny,k,m)=wrkpx(i,ny,k,m)
        enddo
    enddo
    !
    !     calculate the difference  between x nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny-1
            j1=j+1
            do i=1,nx
                delta(i,j,k)=wrkpx(i,j1,k,m)-wrkpx(i,j,k,m)  !non-diffuse n+1 solution
                delta1(i,j,k)=px(i,j1,k,m)-px(i,j,k,m)   !diffused n+1 solution
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            delta(i,ny,k)=delta(i,ny-1,k)
            delta1(i,ny,k)=delta1(i,ny-1,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            do i=1,nx
                s=sign(1.,delta(i,j,k))
                fn(i,j,k)=s*amax1(0.,amin1(s*delta1(i,j-1,k), &
                chipx*abs(delta(i,j,k)),s*delta1(i,j+1,k)))
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            fn(i,1,k)=fn(i,2,k)
            fn(i,ny,k)=fn(i,ny-1,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            j1=j-1
            do i=1,nx
                !        wrkpx(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i,j1,k) ! for x-y-z smoothing
                px(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i,j1,k) ! for x-y smoothing
            enddo
        enddo
    enddo
    !
    return
end
