subroutine csmooth(curx,cury,curz,fldx,fldy,fldz, &
    nx,ny,nz)
    !
    !     force smoothing of surface currents to see if we can
    !             get better graphics
    !
    dimension curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    fldx(nx,ny,nz),fldy(nx,ny,nz),fldz(nx,ny,nz)
    !
    do k=1,nz
        do j=1,ny
            do i=1,nx
                fldx(i,j,k)=0.
                fldy(i,j,k)=0.
                fldz(i,j,k)=0.
            enddo
        enddo
    enddo
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                !
                fldx(i,j,k)= &
                ( curx(i+1,j,k)+curx(i-1,j,k) +curx(i,j+1,k)+curx(i,j-1,k) &
                +curx(i,j,k+1)+curx(i,j,k-1) +6.*curx(i,j,k))/12.
                fldy(i,j,k)= &
                ( cury(i+1,j,k)+cury(i-1,j,k) +cury(i,j+1,k)+cury(i,j-1,k) &
                +cury(i,j,k+1)+cury(i,j,k-1) +6.*cury(i,j,k))/12.
                fldz(i,j,k)= &
                ( curz(i+1,j,k)+curz(i-1,j,k) +curz(i,j+1,k)+curz(i,j-1,k) &
                +curz(i,j,k+1)+curz(i,j,k-1) +6.*curz(i,j,k))/12.
    			!
            enddo
        enddo
    enddo
    !
    return
end
