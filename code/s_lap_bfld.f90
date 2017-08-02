subroutine lap_bfld(bx,by,bz,wrkbx,wrkby,wrkbz,vx,vy,vz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt)
    !
    !     apply the ladipus smoothing technique for particular ion component
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
    wrkbz(nx,ny,nz,ngrd), &
    vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !     equal waiting irrespective of grid size
    !
    deltz=delt
    delty=delt
    deltx=delt
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kp=k+1
        do j=2,ny-1
            jm=j-1
            jp=j+1
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                uxp1=abs(vx(ip,j,k)-vx(i,j,k))
                uxm1=abs(vx(i,j,k)-vx(im,j,k))
                !
                uyp1=abs(vy(i,jp,k)-vy(i,j,k))
                uym1=abs(vy(i,j,k)-vy(i,jm,k))
                !
                uzp1=abs(vz(i,j,kp)-vz(i,j,k))
                uzm1=abs(vz(i,j,k)-vz(i,j,km))
                !
                wrkbx(i,j,k,m)=bx(i,j,k,m)+chierg*( &
                deltx*(uxp1*(bx(ip,j,k,m)-bx(i,j,k,m)) &
                -uxm1*(bx(i,j,k,m)-bx(im,j,k,m))) &
                +delty*(uyp1*(bx(i,jp,k,m)-bx(i,j,k,m)) &
                -uym1*(bx(i,j,k,m)-bx(i,jm,k,m))) &
                +deltz*(uzp1*(bx(i,j,kp,m)-bx(i,j,k,m)) &
                -uzm1*(bx(i,j,k,m)-bx(i,j,km,m))) )
                wrkby(i,j,k,m)=by(i,j,k,m)+chierg*( &
                deltx*(uxp1*(by(ip,j,k,m)-by(i,j,k,m)) &
                -uxm1*(by(i,j,k,m)-by(im,j,k,m))) &
                +delty*(uyp1*(by(i,jp,k,m)-by(i,j,k,m)) &
                -uym1*(by(i,j,k,m)-by(i,jm,k,m))) &
                +deltz*(uzp1*(by(i,j,kp,m)-by(i,j,k,m)) &
                -uzm1*(by(i,j,k,m)-by(i,j,km,m))) )
                wrkbz(i,j,k,m)=bz(i,j,k,m)+chierg*( &
                deltx*(uxp1*(bz(ip,j,k,m)-bz(i,j,k,m)) &
                -uxm1*(bz(i,j,k,m)-bz(im,j,k,m))) &
                +delty*(uyp1*(bz(i,jp,k,m)-bz(i,j,k,m)) &
                -uym1*(bz(i,j,k,m)-bz(i,jm,k,m))) &
                +deltz*(uzp1*(bz(i,j,kp,m)-bz(i,j,k,m)) &
                -uzm1*(bz(i,j,k,m)-bz(i,j,km,m))) )
            enddo
        enddo
    enddo
    !
    !
    return
end
