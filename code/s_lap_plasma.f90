subroutine lap_plasma(rho,px,py,pz, &
    presx,presy,presz,presxy,presxz,presyz, &
    wrkrho,wrkpx,wrkpy,wrkpz, &
    wrkpresx,wrkpresy,wrkpresz, &
    wrkpresxy,wrkpresxz,wrkpresyz, &
    vx,vy,vz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt,isotropic)
    !
    !     apply the ladipus smoothing technique for particular ion component
    !
    dimension rho(nx,ny,nz,ngrd),presx(nx,ny,nz,ngrd), &
    presy(nx,ny,nz,ngrd),presz(nx,ny,nz,ngrd), &
    px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd), &
    wrkrho(nx,ny,nz,ngrd),wrkpresx(nx,ny,nz,ngrd), &
    wrkpresy(nx,ny,nz,ngrd),wrkpresz(nx,ny,nz,ngrd), &
    wrkpx(nx,ny,nz,ngrd),wrkpy(nx,ny,nz,ngrd), &
    wrkpz(nx,ny,nz,ngrd), &
    vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    dimension presxy(nx,ny,nz,ngrd),presxz(nx,ny,nz,ngrd), &
    presyz(nx,ny,nz,ngrd),wrkpresxy(nx,ny,nz,ngrd), &
    wrkpresxz(nx,ny,nz,ngrd),wrkpresyz(nx,ny,nz,ngrd)
    !
    !
    common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9), &
    grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9), &
    rx,ry,rz,xdip,ydip,zdip,rearth,b0, &
    sin_tilt,cos_tilt
    !
    !     equal waiting irrespective of grid size
    !
    deltz=delt
    delty=delt
    deltx=delt
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do  k=2,nz-1
        km=k-1
        kp=k+1
        do  j=2,ny-1
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
                wrkrho(i,j,k,m)=rho(i,j,k,m)+chirho*( &
                    deltx*(uxp1*(rho(ip,j,k,m)-rho(i,j,k,m)) &
                    -uxm1*(rho(i,j,k,m)-rho(im,j,k,m))) &
                    +delty*(uyp1*(rho(i,jp,k,m)-rho(i,j,k,m)) &
                    -uym1*(rho(i,j,k,m)-rho(i,jm,k,m))) &
                    +deltz*(uzp1*(rho(i,j,kp,m)-rho(i,j,k,m)) &
                    -uzm1*(rho(i,j,k,m)-rho(i,j,km,m))) )
		!
                wrkpx(i,j,k,m)=px(i,j,k,m)+chipxyz*( &
                    +deltx*(uxp1*(px(ip,j,k,m)-px(i,j,k,m)) &
                    -uxm1*(px(i,j,k,m)-px(im,j,k,m))) &
                    +delty*(uyp1*(px(i,jp,k,m)-px(i,j,k,m)) &
                    -uym1*(px(i,j,k,m)-px(i,jm,k,m))) &
                    +deltz*(uzp1*(px(i,j,kp,m)-px(i,j,k,m)) &
                    -uzm1*(px(i,j,k,m)-px(i,j,km,m))) )
                !
                wrkpy(i,j,k,m)=py(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(py(ip,j,k,m)-py(i,j,k,m)) &
                    -uxm1*(py(i,j,k,m)-py(im,j,k,m))) &
                    +delty*(uyp1*(py(i,jp,k,m)-py(i,j,k,m)) &
                    -uym1*(py(i,j,k,m)-py(i,jm,k,m))) &
                    +deltz*(uzp1*(py(i,j,kp,m)-py(i,j,k,m)) &
                    -uzm1*(py(i,j,k,m)-py(i,j,km,m))) )
                !
                wrkpz(i,j,k,m)=pz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(pz(ip,j,k,m)-pz(i,j,k,m)) &
                    -uxm1*(pz(i,j,k,m)-pz(im,j,k,m))) &
                    +delty*(uyp1*(pz(i,jp,k,m)-pz(i,j,k,m)) &
                    -uym1*(pz(i,j,k,m)-pz(i,jm,k,m))) &
                    +deltz*(uzp1*(pz(i,j,kp,m)-pz(i,j,k,m)) &
                    -uzm1*(pz(i,j,k,m)-pz(i,j,km,m))) )

                !
                !      diagonal elements
                !
                wrkpresx(i,j,k,m)=presx(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presx(ip,j,k,m)-presx(i,j,k,m)) &
                    -uxm1*(presx(i,j,k,m)- presx(im,j,k,m))) &
                    +delty*(uyp1*(presx(i,jp,k,m)-presx(i,j,k,m)) &
                    -uym1*(presx(i,j,k,m)- presx(i,jm,k,m))) &
                    +deltz*(uzp1*(presx(i,j,kp,m)-presx(i,j,k,m)) &
                    -uzm1*(presx(i,j,k,m)- presx(i,j,km,m))) )
                !
                wrkpresy(i,j,k,m)=presy(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presy(ip,j,k,m)-presy(i,j,k,m)) &
                    -uxm1*(presy(i,j,k,m)- presy(im,j,k,m))) &
                    +delty*(uyp1*(presy(i,jp,k,m)-presy(i,j,k,m)) &
                    -uym1*(presy(i,j,k,m)- presy(i,jm,k,m))) &
                    +deltz*(uzp1*(presy(i,j,kp,m)-presy(i,j,k,m)) &
                    -uzm1*(presy(i,j,k,m)- presy(i,j,km,m))) )
                !
                wrkpresz(i,j,k,m)=presz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presz(ip,j,k,m)-presz(i,j,k,m)) &
                    -uxm1*(presz(i,j,k,m)- presz(im,j,k,m))) &
                    +delty*(uyp1*(presz(i,jp,k,m)-presz(i,j,k,m)) &
                    -uym1*(presz(i,j,k,m)- presz(i,jm,k,m))) &
                    +deltz*(uzp1*(presz(i,j,kp,m)-presz(i,j,k,m)) &
                    -uzm1*(presz(i,j,k,m)- presz(i,j,km,m))) )
                !
                !      off-diagonal elements
                !
                wrkpresxy(i,j,k,m)=presxy(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presxy(ip,j,k,m)-presxy(i,j,k,m)) &
                    -uxm1*(presxy(i,j,k,m)- presxy(im,j,k,m))) &
                    +delty*(uyp1*(presxy(i,jp,k,m)-presxy(i,j,k,m)) &
                    -uym1*(presxy(i,j,k,m)- presxy(i,jm,k,m))) &
                    +deltz*(uzp1*(presxy(i,j,kp,m)-presxy(i,j,k,m)) &
                    -uzm1*(presxy(i,j,k,m)- presxy(i,j,km,m))) )
                !
                wrkpresxz(i,j,k,m)=presxz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presxz(ip,j,k,m)-presxz(i,j,k,m)) &
                    -uxm1*(presxz(i,j,k,m)- presxz(im,j,k,m))) &
                    +delty*(uyp1*(presxz(i,jp,k,m)-presxz(i,j,k,m)) &
                    -uym1*(presxz(i,j,k,m)- presxz(i,jm,k,m))) &
                    +deltz*(uzp1*(presxz(i,j,kp,m)-presxz(i,j,k,m)) &
                    -uzm1*(presxz(i,j,k,m)- presxz(i,j,km,m))) )
                !
                wrkpresyz(i,j,k,m)=presyz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presyz(ip,j,k,m)-presyz(i,j,k,m)) &
                    -uxm1*(presyz(i,j,k,m)- presyz(im,j,k,m))) &
                    +delty*(uyp1*(presyz(i,jp,k,m)-presyz(i,j,k,m)) &
                    -uym1*(presyz(i,j,k,m)- presyz(i,jm,k,m))) &
                    +deltz*(uzp1*(presyz(i,j,kp,m)-presyz(i,j,k,m)) &
                    -uzm1*(presyz(i,j,k,m)- presyz(i,j,km,m))) )
            enddo
        enddo
    enddo
    !
    return
end
