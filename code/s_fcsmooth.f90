!
!	This file contains two subroutines:
!	fcsmooth
!	fcsmooth_2d
!
subroutine fcsmooth(px,oldpx,wrkpx,nx,ny,nz,n_grids,box,chipx, &
    tempx,tempy,tempz)
    !
    !     applies flux correction smoothing to field quantities
    !           oldpx is the qunatity at the start of the time step
    !           wrkpx is the unsmoothed qunatity
    !           px is the final product
    !     wrk arrays assumed to have dimension larger than nx,ny,nz
    !
    integer box
    dimension px(nx,ny,nz,n_grids),oldpx(nx,ny,nz,n_grids), &
    wrkpx(nx,ny,nz,n_grids), &
    tempx(nx,ny,nz),tempy(nx,ny,nz),tempz(nx,ny,nz)
    !
    real,dimension(nx,ny,nz) :: deltax,delta1x,fnx, &
    deltay,delta1y,fny, &
    deltaz,delta1z,fnz
    !
    !      save initial unperturbed values in temp array
	!		MJS 4/3/19: What's the point of this? We just assign
	!			values on top of these in the next step, without
	!			using the temporary values in the interim.
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                tempx(i,j,k)=wrkpx(i,j,k,box)
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
                px(i,j,k,box)=wrkpx(i,j,k,box)+chipx*( &
                oldpx(i+1,j,k,box)+oldpx(i-1,j,k,box) &
                +oldpx(i,jp,k,box)+oldpx(i,jm,k,box) &
                +oldpx(i,j,kp,box)+oldpx(i,j,km,box) &
                -6.*oldpx(i,j,k,box))
            enddo
        enddo
    enddo
    !
    !     set boundary points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            px(1,j,k,box)=wrkpx(1,j,k,box)
            px(nx,j,k,box)=wrkpx(nx,j,k,box)
            tempx(1,j,k)=wrkpx(1,j,k,box)
            tempx(nx,j,k)=wrkpx(nx,j,k,box)
        enddo
    enddo
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            px(i,j,1,box)=wrkpx(i,j,1,box)
            px(i,j,nz,box)=wrkpx(i,j,nz,box)
            tempx(i,j,1)=wrkpx(i,j,1,box)
            tempx(i,j,nz)=wrkpx(i,j,nz,box)
        enddo
    enddo
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            px(i,1,k,box)=wrkpx(i,1,k,box)
            px(i,ny,k,box)=wrkpx(i,ny,k,box)
            tempx(i,1,k)=wrkpx(i,1,k,box)
            tempx(i,ny,k)=wrkpx(i,ny,k,box)
        enddo
    enddo
    !
    !     calculate the x difference between nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx-1
                i1=i+1
                deltax(i,j,k)=wrkpx(i1,j,k,box)-wrkpx(i,j,k,box)   !d-non-diffuse n+1 sol'n
                delta1x(i,j,k)=px(i1,j,k,box)-px(i,j,k,box)   !d1-diffused n+1 sol'n
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
                tempx(i,j,k)=px(i,j,k,box)-fnx(i,j,k)+fnx(i1,j,k) ! all smoothing
                !        px(i,j,k,box)=px(i,j,k,box)-fnx(i)+fnx(i1)  ! x-smooth only
            enddo
        enddo
    enddo
    !
    !     calculate the y difference between nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny-1
            j1=j+1
            do i=1,nx
                deltay(i,j,k)=wrkpx(i,j1,k,box)-wrkpx(i,j,k,box)  !non-diffuse n+1 solution
                delta1y(i,j,k)=px(i,j1,k,box)-px(i,j,k,box)   !diffused n+1 solution
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            deltay(i,ny,k)=deltay(i,ny-1,k)
            delta1y(i,ny,k)=delta1y(i,ny-1,k)
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
                !        px(i,j,k,box)=tempx(i,j,k)-fny(j)+fny(j1) ! for x-y smoothing
            enddo
        enddo
    enddo
    !
    !     calculate the z difference between nearest grid points
    !
    !$omp  parallel do
    do k=1,nz-1
        do j=1,ny
            do i=1,nx
                deltaz(i,j,k)=wrkpx(i,j,k+1,box)-wrkpx(i,j,k,box)  !non-diffuse n+1 soln
                delta1z(i,j,k)=px(i,j,k+1,box)-px(i,j,k,box)   !diffused n+1 solution
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
                px(i,j,k,box)=tempx(i,j,k)-fnz(i,j,k)+fnz(i,j,k-1)
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine fcsmooth_2d(px,oldpx,wrkpx,nx,ny,nz,n_grids,box,chipx, &
    delta,delta1,fn)
    !
    !     applies flux correction smoothing to field quantities
    !           oldpx is the quantity at the start of the time step
    !           wrkpx is the unsmoothed quantity
    !           px is the final product
    !     wrk arrays assumed to have dimension larger than nx,ny,nz
    !
    integer box
    dimension px(nx,ny,nz,n_grids),oldpx(nx,ny,nz,n_grids), &
    wrkpx(nx,ny,nz,n_grids), &
    delta(nx,ny,nz),delta1(nx,ny,nz),fn(nx,ny,nz)
    !
    !         step 1:   diffuse initial in x direction
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                px(i,j,k,box)=wrkpx(i,j,k,box)+chipx*( &
                oldpx(i+1,j,k,box)+oldpx(i-1,j,k,box) &
                -2.*oldpx(i,j,k,box))
            enddo
        enddo
    enddo
    !
    !     set boundary points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            px(1,j,k,box)=wrkpx(1,j,k,box)
            px(nx,j,k,box)=wrkpx(nx,j,k,box)
        enddo
    enddo
    !
    !     calculate the difference between nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx-1
                i1=i+1
                delta(i,j,k)=wrkpx(i1,j,k,box)-wrkpx(i,j,k,box)   !d-non-diffuse n+1 sol'n
                delta1(i,j,k)=px(i1,j,k,box)-px(i,j,k,box)   !d1-diffused n+1 sol'n
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
                wrkpx(i,j,k,box)=px(i,j,k,box)-fn(i,j,k)+fn(i1,j,k) ! all smoothing
                !        px(i,j,k,box)=px(i,j,k,box)-fn(i,j,k)+fn(i1,j,k)  ! x-smooth only
            enddo
        enddo
    enddo
    !
    !         step 2:   diffuse initial in y-direction
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=1,nx
                px(i,j,k,box)=wrkpx(i,j,k,box)+chipx*( &
                oldpx(i,jp,k,box)+oldpx(i,jm,k,box) &
                -2.*oldpx(i,j,k,box))
            enddo
        enddo
    enddo
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            px(i,1,k,box)=wrkpx(i,1,k,box)
            px(i,ny,k,box)=wrkpx(i,ny,k,box)
        enddo
    enddo
    !
    !     calculate the difference between nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny-1
            j1=j+1
            do i=1,nx
                delta(i,j,k)=wrkpx(i,j1,k,box)-wrkpx(i,j,k,box)  !non-diffuse n+1 solution
                delta1(i,j,k)=px(i,j1,k,box)-px(i,j,k,box)   !diffused n+1 solution
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
                !        wrkpx(i,j,k,box)=px(i,j,k,box)-fn(i,j,k)+fn(i,j1,k) ! for x-y-z smoothing
                px(i,j,k,box)=px(i,j,k,box)-fn(i,j,k)+fn(i,j1,k) ! for x-y smoothing
            enddo
        enddo
    enddo
    !
    return
end
