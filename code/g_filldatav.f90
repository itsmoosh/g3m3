subroutine filldatav(stuff,vx,vy,vz,nx,ny,nz,ngrd,mm,m, &
    t,tt,tx,ty,tz,mx,my,mz,xmin,ymin,zmin, &
    delx,dely,delz,vm,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     fill in diagnostic data cubes
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension stuff(nx,ny,nz,ngrd),vx(nx,ny,nz), &
    vy(nx,ny,nz),vz(nx,ny,nz)
    dimension t(mx,my,mz),tt(mx,my,mz), &
    tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz)
    !
    !      load t stuff
    !
    ddx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    ddy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    ddz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    vm=0.00001
    do k=1,mz
        az=zmin+delz*(k-1)
        ak=1.+(az-grd_zmin(m))/ddz
        k1=ak
        k2=k1+1
        dz=ak-k1
        !
        do j=1,my
            ay=ymin+dely*(j-1)
            aj=1.+(ay-grd_ymin(m))/ddy
            j1=aj
            j2=j1+1
            dy=aj-j1
            !
            do i=1,mx
                ax=xmin+delx*(i-1)
                ai=1.+(ax-grd_xmin(m))/ddx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                t(i,j,k)=stuff(i1,j1,k1,mm)*(1.-dx)*(1.-dy)*(1.-dz) &
                +stuff(i1,j1,k2,mm)*(1.-dx)*(1.-dy)*(dz) &
                +stuff(i1,j2,k1,mm)*(1.-dx)*(dy)*(1.-dz) &
                +stuff(i1,j2,k2,mm)*(1.-dx)*(dy)*(dz) &
                +stuff(i2,j1,k1,mm)*(dx)*(1.-dy)*(1.-dz) &
                +stuff(i2,j1,k2,mm)*(dx)*(1.-dy)*(dz) &
                +stuff(i2,j2,k1,mm)*(dx)*(dy)*(1.-dz) &
                +stuff(i2,j2,k2,mm)*(dx)*(dy)*(dz)
                !
                radius=sqrt(ax**2+ay**2+az**2)
                tt(i,j,k)=radius
                !
                tx(i,j,k)=vx(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
                +vx(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
                +vx(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
                +vx(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
                +vx(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
                +vx(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
                +vx(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
                +vx(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                ty(i,j,k)=vy(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
                +vy(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
                +vy(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
                +vy(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
                +vy(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
                +vy(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
                +vy(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
                +vy(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                tz(i,j,k)=vz(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
                +vz(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
                +vz(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
                +vz(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
                +vz(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
                +vz(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
                +vz(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
                +vz(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                vm=amax1(vm,abs(tx(i,j,k)),abs(ty(i,j,k)),abs(tz(i,j,k)))
                !
            enddo
        enddo
    enddo
    vm=amax1(vm,0.00004)
    !
    return
end
