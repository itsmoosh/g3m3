subroutine lap_elec(ppres,wrkppres,vx,vy,vz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt)
    !
    !     apply the ladipus smoothing technique for particular ion component
    !
    dimension ppres(nx,ny,nz,ngrd),wrkppres(nx,ny,nz,ngrd), &
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
                wrkppres(i,j,k,m)=ppres(i,j,k,m)+chipxyz*( &
                deltx*(uxp1*(ppres(ip,j,k,m)-ppres(i,j,k,m)) &
                -uxm1*(ppres(i,j,k,m)-ppres(im,j,k,m))) &
                +delty*(uyp1*(ppres(i,jp,k,m)-ppres(i,j,k,m)) &
                -uym1*(ppres(i,j,k,m)-ppres(i,jm,k,m))) &
                +deltz*(uzp1*(ppres(i,j,kp,m)-ppres(i,j,k,m)) &
                -uzm1*(ppres(i,j,k,m)-ppres(i,j,km,m))) )
            enddo
        enddo
    enddo
    !
    return
end
