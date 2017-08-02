subroutine fcsmooth(px,oldpx,wrkpx,nx,ny,nz,ngrd,m,chipx, &
    tempx,tempy,tempz)
    !
    !     applies flux correction smoothing to field quantities
    !           oldpx is the qunatity at the start of the time step
    !           wrkpx is the unsmoothed qunatity
    !           px is the final product
    !     wrk arrays assumed to have dimension larger than nx,ny,nz
    !
    dimension px(nx,ny,nz,ngrd),oldpx(nx,ny,nz,ngrd), &
    wrkpx(nx,ny,nz,ngrd), &
    tempx(nx,ny,nz),tempy(nx,ny,nz),tempz(nx,ny,nz)
    !
    real,dimension(nx,ny,nz) :: deltax,delta1x,fnx, &
    deltay,delta1y,fny, &
    deltaz,delta1z,fnz
    !
    !      save initial unperturbed values in temp array
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                tempx(i,j,k)=wrkpx(i,j,k,m)
            enddo
        enddo
    enddo
    !
    !         step 1:   diffuse initial in all directions
    !$omp  parallel do
    do k=2,nz-1
        kp=k+1
        km=k-1
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=2,nx-1
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                oldpx(i+1,j,k,m)+oldpx(i-1,j,k,m) &
                +oldpx(i,jp,k,m)+oldpx(i,jm,k,m) &
                +oldpx(i,j,kp,m)+oldpx(i,j,km,m) &
                -6.*oldpx(i,j,k,m))
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
            tempx(1,j,k)=wrkpx(1,j,k,m)
            tempx(nx,j,k)=wrkpx(nx,j,k,m)
        enddo
    enddo
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            px(i,j,1,m)=wrkpx(i,j,1,m)
            px(i,j,nz,m)=wrkpx(i,j,nz,m)
            tempx(i,j,1)=wrkpx(i,j,1,m)
            tempx(i,j,nz)=wrkpx(i,j,nz,m)
        enddo
    enddo
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            px(i,1,k,m)=wrkpx(i,1,k,m)
            px(i,ny,k,m)=wrkpx(i,ny,k,m)
            tempx(i,1,k)=wrkpx(i,1,k,m)
            tempx(i,ny,k)=wrkpx(i,ny,k,m)
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
                deltax(i,j,k)=wrkpx(i1,j,k,m)-wrkpx(i,j,k,m)   !d-non-diffuse n+1 sol'n
                delta1x(i,j,k)=px(i1,j,k,m)-px(i,j,k,m)   !d1-diffused n+1 sol'n
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            deltax(nx,j,k)=deltax(nx-1,j,k)
            delta1x(nx,j,k)=delta1x(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                s=sign(1.,deltax(i,j,k))
                fnx(i,j,k)=s*amax1(0.,amin1(s*delta1x(i-1,j,k), &
                chipx*abs(deltax(i,j,k)),s*delta1x(i+1,j,k)))
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            fnx(1,j,k)=fnx(2,j,k)
            fnx(nx,j,k)=fnx(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                i1=i-1
                tempx(i,j,k)=px(i,j,k,m)-fnx(i,j,k)+fnx(i1,j,k) ! all smoothing
                !        px(i,j,k,m)=px(i,j,k,m)-fnx(i)+fnx(i1)  ! x-smooth only
            enddo
        enddo
    enddo
    !
    !
    !
    !     calculate the difference  between y nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny-1
            j1=j+1
            do i=1,nx
                deltay(i,j,k)=wrkpx(i,j1,k,m)-wrkpx(i,j,k,m)  !non-diffuse n+1 solution
                delta1y(i,j,k)=px(i,j1,k,m)-px(i,j,k,m)   !diffused n+1 solution
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            deltay(i,ny,j)=deltay(i,ny-1,k)
            delta1y(i,ny,j)=delta1y(i,ny-1,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            do i=1,nx
                s=sign(1.,deltay(i,j,k))
                fny(i,j,k)=s*amax1(0.,amin1(s*delta1y(i,j-1,k), &
                chipx*abs(deltay(i,j,k)),s*delta1y(i,j+1,k)))
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            fny(i,1,k)=fny(i,2,k)
            fny(i,ny,k)=fny(i,ny-1,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            j1=j-1
            do i=1,nx
                tempx(i,j,k)=tempx(i,j,k)-fny(i,j,k)+fny(i,j1,k) ! for x-y-z smoothing
                !        px(i,j,k,m)=tempx(i,j,k)-fny(j)+fny(j1) ! for x-y smoothing
            enddo
        enddo
    enddo
    !
    !     calculatethe difference  between z nearest grid points
    !
    !$omp  parallel do
    do k=1,nz-1
        do j=1,ny
            do i=1,nx
                deltaz(i,j,k)=wrkpx(i,j,k+1,m)-wrkpx(i,j,k,m)  !non-diffuse n+1 soln
                delta1z(i,j,k)=px(i,j,k+1,m)-px(i,j,k,m)   !diffused n+1 solution
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            deltaz(i,j,nz)=deltaz(i,j,nz-1)
            delta1z(i,j,nz)=delta1z(i,j,nz-1)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=2,nz-1
        do j=1,ny
            do i=1,nx
                s=sign(1.,deltaz(i,j,k))
                fnz(i,j,k)=s*amax1(0.,amin1(s*delta1z(i,j,k-1), &
                chipx*abs(deltaz(i,j,k)),s*delta1z(i,j,k+1)))
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            fnz(i,j,1)=fnz(i,j,2)
            fnz(i,j,nz)=fnz(i,j,nz-1)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=2,nz-1
        do j=1,ny
            do i=1,nx
                px(i,j,k,m)=tempx(i,j,k)-fnz(i,j,k)+fnz(i,j,k-1)
            enddo
        enddo
    enddo
    !
    return
end
