subroutine calcur(bx,by,bz,nx,ny,nz,n_grids,box,curx,cury,curz, &
    rx,ry,rz)
    !
    !     this calculates the current associated with the perturbed b field
    !          i.e. curz= dby/dx - dbx/dy
    !
    dimension bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz)
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kp=k+1
    	!
        do j=2,ny-1
            jm=j-1
            jp=j+1
    		!
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                curx(i,j,k)=(bz(i,jp,k,box)-bz(i,jm,k,box))/dyt &
                - (by(i,j,kp,box)-by(i,j,km,box))/dzt
                !
                cury(i,j,k)=(bx(i,j,kp,box)-bx(i,j,km,box))/dzt &
                - (bz(ip,j,k,box)-bz(im,j,k,box))/dxt
                !
                curz(i,j,k)=(by(ip,j,k,box)-by(im,j,k,box))/dxt &
                - (bx(i,jp,k,box)-bx(i,jm,k,box))/dyt
            enddo
        enddo
    enddo
    !
    !     following boundary conditions are set so there is no forward
    !      communication, i.e. same x required
    !      and symmetry between k=1 and k=2
    !
    !     set boundary regions -- bottom and top panels
    !
    nx1=nx-1
    ny1=ny-1
    nz1=nz-1
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do j=2,ny1
        do i=2,nx1
            !
            !       regular boundary conditions
            !
            curx(i,j,1)=curx(i,j,2)
            curx(i,j,nz)=curx(i,j,nz1)
            !
            cury(i,j,1)=cury(i,j,2)
            cury(i,j,nz)=cury(i,j,nz1)
            !
            curz(i,j,1)=curz(i,j,2)
            curz(i,j,nz)=curz(i,j,nz1)
        enddo
    enddo
    !
    !       set boundary regions -- front and back
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz1
        do j=2,ny1
            !
            curx(1,j,k)=curx(2,j,k)
            curx(nx,j,k)=curx(nx1,j,k)
            !
            cury(1,j,k)=cury(2,j,k)
            cury(nx,j,k)=cury(nx1,j,k)
            !
            curz(1,j,k)=curz(2,j,k)
            curz(nx,j,k)=curz(nx1,j,k)
        enddo
    enddo
    !
    !       set boundary regions -- left and right
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do j=2,nz1
        do i=2,nx1
            curx(i,1,j)=curx(i,2,j)
            curx(i,ny,j)=curx(i,ny1,j)
            !
            cury(i,1,j)=cury(i,2,j)
            cury(i,ny,j)=cury(i,ny1,j)
            !
            curz(i,1,j)=curz(i,2,j)
            curz(i,ny,j)=curz(i,ny1,j)
        enddo
    enddo
    !
    !     set corner lines
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do i=2,nx1
        curx(i,1,1)=curx(i,1,2)
        curx(i,1,nz)=curx(i,1,nz1)
        cury(i,1,1)=cury(i,1,2)
        cury(i,1,nz)=cury(i,1,nz1)
        curz(i,1,1)=curz(i,1,2)
        curz(i,1,nz)=curz(i,1,nz1)
        !
        curx(i,ny,1)=curx(i,ny,2)
        curx(i,ny,nz)=curx(i,ny1,nz)
        cury(i,ny,1)=cury(i,ny,2)
        cury(i,ny,nz)=cury(i,ny1,nz)
        curz(i,ny,1)=curz(i,ny,2)
        curz(i,ny,nz)=curz(i,ny1,nz)
    enddo
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do j=2,ny1
        curx(1,j,1)=curx(1,j,2)
        curx(1,j,nz)=curx(1,j,nz1)
        cury(1,j,1)=cury(1,j,2)
        cury(1,j,nz)=cury(1,j,nz1)
        curz(1,j,1)=curz(1,j,2)
        curz(1,j,nz)=curz(1,j,nz1)
        !
        curx(nx,j,1)=curx(nx,j,2)
        curx(nx,j,nz)=curx(nx,j,nz1)
        cury(nx,j,1)=cury(nx,j,2)
        cury(nx,j,nz)=cury(nx,j,nz1)
        curz(nx,j,1)=curz(nx,j,2)
        curz(nx,j,nz)=curz(nx,j,nz1)
    enddo
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz1
        curx(1,1,k)=curx(1,2,k)
        curx(nx,1,k)=curx(nx,2,k)
        cury(1,1,k)=cury(1,2,k)
        cury(nx,1,k)=cury(nx,2,k)
        curz(1,1,k)=curz(1,2,k)
        curz(nx,1,k)=curz(nx,2,k)
        !
        curx(1,ny,k)=curx(1,ny1,k)
        curx(nx,ny,k)=curx(nx,ny1,k)
        cury(1,ny,k)=cury(1,ny1,k)
        cury(nx,ny,k)=cury(nx,ny1,k)
        curz(1,ny,k)=curz(1,ny1,k)
        curz(nx,ny,k)=curz(nx,ny1,k)
        !
    enddo
    !
    curx(1,1,1)=curx(1,1,2)
    curx(1,ny,1)=curx(1,ny,2)
    curx(1,1,nz)=(curx(1,1,nz1)+curx(1,2,nz))/2.
    curx(1,ny,nz)=(curx(1,ny,nz1)+curx(1,ny1,nz))/2.
    curx(nx,1,1)=(curx(nx,1,2)+curx(nx,2,1)+curx(nx1,1,1))/3.
    curx(nx,ny,1)=(curx(nx,ny,2)+curx(nx,ny1,1)+curx(nx1,ny,1))/3.
    curx(nx,1,nz)=(curx(nx,1,nz1)+curx(nx,2,nz)+curx(nx1,1,nz))/3.
    curx(nx,ny,nz)=(curx(nx,ny,nz1)+curx(nx,ny1,nz) &
    +curx(nx1,ny,nz))/3.
    !
    cury(1,1,1)=cury(1,1,2)
    cury(1,ny,1)=cury(1,ny,2)
    cury(1,1,nz)=(cury(1,1,nz1)+cury(1,2,nz))/2.
    cury(1,ny,nz)=(cury(1,ny,nz1)+cury(1,ny1,nz))/2.
    cury(nx,1,1)=(cury(nx,1,2)+cury(nx,2,1)+cury(nx1,1,1))/3.
    cury(nx,ny,1)=(cury(nx,ny,2)+cury(nx,ny1,1)+cury(nx1,ny,1))/3.
    cury(nx,1,nz)=(cury(nx,1,nz1)+cury(nx,2,nz)+cury(nx1,1,nz))/3.
    cury(nx,ny,nz)=(cury(nx,ny,nz1)+cury(nx,ny1,nz) &
    +cury(nx1,ny,nz))/3.
    !
    curz(1,1,1)=curz(1,1,2)
    curz(1,ny,1)=curz(1,ny,2)
    curz(1,1,nz)=(curz(1,1,nz1)+curz(1,2,nz))/2.
    curz(1,ny,nz)=(curz(1,ny,nz1)+curz(1,ny1,nz))/2.
    curz(nx,1,1)=(curz(nx,1,2)+curz(nx,2,1)+curz(nx1,1,1))/3.
    curz(nx,ny,1)=(curz(nx,ny,2)+curz(nx,ny1,1)+curz(nx1,ny,1))/3.
    curz(nx,1,nz)=(curz(nx,1,nz1)+curz(nx,2,nz)+curz(nx1,1,nz))/3.
    curz(nx,ny,nz)=(curz(nx,ny,nz1)+curz(nx,ny1,nz) &
    +curz(nx1,ny,nz))/3.
    !
    return
end
