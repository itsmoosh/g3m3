!
!	This file contains 3 subroutines:
!	rungea
!	rungeb
!	rungem
!
subroutine rungea(bx,by,bz,press,rho,nx,ny,nz, &
    m,xi,yi,zi,xmin,xmax,ymin,ymax, &
    zmin,zmax,add_dip,ergies,rhod,maxpts,dir, &
    delx,dely,delz,rearth,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine uses a fourth order runge-kutta method
    !         to trace stream functions of current and magnetic fld
    !         xi,yi,zi, and find maximum temperature along them
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz), &
    press(nx,ny,nz),rho(nx,ny,nz)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    real xray(1000),yray(1000),zray(1000)
    !
    logical add_dip,roc
    !
    maxpts=1000
    !
    tstep=0.25*dir
    !
    xray(1)=xi
    yray(1)=yi
    zray(1)=zi
    !
    ts6=tstep/6.
    do n=2,maxpts
        !
        !           step 1
        !
        ax=xray(n-1)
        ay=yray(n-1)
        az=zray(n-1)
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx1,dy1,dz1,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 2
        !
        ax=xray(n-1)+tstep*0.5*dx1
        ay=yray(n-1)+tstep*0.5*dy1
        az=zray(n-1)+tstep*0.5*dz1
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx2,dy2,dz2,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 3
        !
        ax=xray(n-1)+tstep*0.5*dx2
        ay=yray(n-1)+tstep*0.5*dy2
        az=zray(n-1)+tstep*0.5*dz2
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx3,dy3,dz3,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 4
        !
        ax=xray(n-1)+tstep*dx3
        ay=yray(n-1)+tstep*dy3
        az=zray(n-1)+tstep*dz3
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx4,dy4,dz4,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       update new position and wavevector
        !
        xray(n)=xray(n-1)+ts6*(dx1+2.*dx2+2.*dx3+dx4)
        yray(n)=yray(n-1)+ts6*(dy1+2.*dy2+2.*dy3+dy4)
        zray(n)=zray(n-1)+ts6*(dz1+2.*dz2+2.*dz3+dz4)
        !
        !    find energy at the point
        !
        ax=xray(n)
        ay=yray(n)
        az=zray(n)
        !
        !          interpolate data to grid point
        !
        ak=1.+(az-grd_zmin(m))/delz
        k1=ak
        k2=k1+1
        dz=ak-k1
        !
        aj=1.+(ay-grd_ymin(m))/dely
        j1=aj
        j2=j1+1
        dy=aj-j1
        !
        ai=1.+(ax-grd_xmin(m))/delx
        i1=ai
        i2=i1+1
        dx=ai-i1
        !
        aerg=press(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
        +press(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
        +press(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
        +press(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
        +press(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
        +press(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
        +press(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
        +press(i2,j2,k2)*(dx)*(dy)*(dz)
        arho=rho(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
        +rho(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
        +rho(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
        +rho(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
        +rho(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
        +rho(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
        +rho(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
        +rho(i2,j2,k2)*(dx)*(dy)*(dz)
    
        !
        aergd=aerg/arho
        arho=amax1(arho,0.0001)
        !
        ergies=amax1(aergd,ergies)
        rhod=amin1(arho,rhod)
        !
        !       test to see if ray is within selected region
        !
        ar=sqrt((xray(n)-xdip)**2+(yray(n)-ydip)**2 &
        +(zray(n)-zdip)**2)
        radius=sqrt((xray(n)-xray(1))**2+(yray(n)-yray(1))**2 &
        +(zray(n)-zray(1))**2)
        if((xray(n).lt.xmin).or.(xray(n).gt.xmax).or. &
        (yray(n).lt.ymin).or.(yray(n).gt.ymax).or. &
        (zray(n).lt.zmin).or.(zray(n).gt.zmax).or. &
        (ar.lt.rearth+1.).or.(radius.gt.4.*rearth).or.(roc)) return
        !
    enddo
    !
    return
end
!
!
!	********************************************
!
!
subroutine rungeb(efldx,efldy,efldz,bx,by,bz,nx,ny,nz, &
    ngrd,m,rx,xi,yi,zi,xmin,xmax,ymin,ymax, &
    zmin,zmax,add_dip,potential,maxpts,dir,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine uses a fourth order runge-kutta method
    !         to trace  stream functions of current and magnetic fld
    !         xi,yi,zi, and integrates the field aligned potential drop
    !
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz)
    real xray(1000),yray(1000),zray(1000)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    logical add_dip,roc
    !
    rx2=rx*3.
    potential=0.
    maxpts=1000
    delx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dely=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    delz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    tstep=0.25*dir
    !
    xray(1)=xi
    yray(1)=yi
    zray(1)=zi
    !
    ts6=tstep/6.
    do n=2,maxpts
        !
        !           step 1
        !
        ax=xray(n-1)
        ay=yray(n-1)
        az=zray(n-1)
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx1,dy1,dz1,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 2
        !
        ax=xray(n-1)+tstep*0.5*dx1
        ay=yray(n-1)+tstep*0.5*dy1
        az=zray(n-1)+tstep*0.5*dz1
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx2,dy2,dz2,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 3
        !
        ax=xray(n-1)+tstep*0.5*dx2
        ay=yray(n-1)+tstep*0.5*dy2
        az=zray(n-1)+tstep*0.5*dz2
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx3,dy3,dz3,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 4
        !
        ax=xray(n-1)+tstep*dx3
        ay=yray(n-1)+tstep*dy3
        az=zray(n-1)+tstep*dz3
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx4,dy4,dz4,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       update new position and wavevector
        !
        xray(n)=xray(n-1)+ts6*(dx1+2.*dx2+2.*dx3+dx4)
        yray(n)=yray(n-1)+ts6*(dy1+2.*dy2+2.*dy3+dy4)
        zray(n)=zray(n-1)+ts6*(dz1+2.*dz2+2.*dz3+dz4)
        !
        !      define integral limits
        !
        xi=xray(n-1)
        yi=yray(n-1)
        zi=zray(n-1)
        !
        xf=xray(n)
        yf=yray(n)
        zf=zray(n)
        !
        !          interpolate data to grid point
        !
        ax=(xf+xi)/2.
        ay=(yf+yi)/2.
        az=(zf+zi)/2.
        !
        ak=1.+(az-grd_zmin(m))/delz
        k1=ak
        k2=k1+1
        dz=ak-k1
        ddz=1.-dz
        !
        aj=1.+(ay-grd_ymin(m))/dely
        j1=aj
        j2=j1+1
        dy=aj-j1
        ddy=1.-dy
        !
        ai=1.+(ax-grd_xmin(m))/delx
        i1=ai
        i2=i1+1
        dx=ai-i1
        ddx=1.-dx
        !
        aex=  efldx(i1,j1,k1)*ddx*ddy*ddz &
        +efldx(i1,j1,k2)*ddx*ddy*dz &
        +efldx(i1,j2,k1)*ddx*dy*ddz &
        +efldx(i1,j2,k2)*ddx*dy*dz &
        +efldx(i2,j1,k1)*dx*ddy*ddz &
        +efldx(i2,j1,k2)*dx*ddy*dz &
        +efldx(i2,j2,k1)*dx*dy*ddz &
        +efldx(i2,j2,k2)*dx*dy*dz
        aey=  efldy(i1,j1,k1)*ddx*ddy*ddz &
        +efldy(i1,j1,k2)*ddx*ddy*dz &
        +efldy(i1,j2,k1)*ddx*dy*ddz &
        +efldy(i1,j2,k2)*ddx*dy*dz &
        +efldy(i2,j1,k1)*dx*ddy*ddz &
        +efldy(i2,j1,k2)*dx*ddy*dz &
        +efldy(i2,j2,k1)*dx*dy*ddz &
        +efldy(i2,j2,k2)*dx*dy*dz
        aez=  efldz(i1,j1,k1)*ddx*ddy*ddz &
        +efldz(i1,j1,k2)*ddx*ddy*dz &
        +efldz(i1,j2,k1)*ddx*dy*ddz &
        +efldz(i1,j2,k2)*ddx*dy*dz &
        +efldz(i2,j1,k1)*dx*ddy*ddz &
        +efldz(i2,j1,k2)*dx*ddy*dz &
        +efldz(i2,j2,k1)*dx*dy*ddz &
        +efldz(i2,j2,k2)*dx*dy*dz
        !
        drx=xf-xi
        dry=yf-yi
        drz=zf-zi
        !
        potential=potential+drx*aex+dry*aey+drz*aez
        !
        !       test to see if ray is within selected region
        !
        ar=sqrt(xray(n)**2+yray(n)**2 +zray(n)**2)
    
        radius=sqrt((xray(n)-xray(1))**2+(yray(n)-yray(1))**2 &
        +(zray(n)-zray(1))**2)
        if((xray(n).lt.xmin).or.(xray(n).gt.xmax).or. &
        (yray(n).lt.ymin).or.(yray(n).gt.ymax).or. &
        (zray(n).lt.zmin).or.(zray(n).gt.zmax).or. &
        (ar.lt.rearth+rx2).or.(radius.gt.3.*rearth) &
        .or.(roc)) return
        !
    enddo
    !
    return
end
!
!
!	********************************************
!
!
subroutine rungem(bx,by,bz,nx,ny,nz,m,rx, &
    xi,yi,zi,xmin,xmax,ymin,ymax,zmin,zmax,rearth, &
    add_dip,xray,yray,zray,maxpts,npts,dir,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine uses a fourth order runge-kutta method
    !         to trace to trace stream functions of current and magnetic fld
    !         xi,yi,zi, until it hits boundary or maximum number
    !         points is reach
    !
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    dimension  xray(maxpts),yray(maxpts),zray(maxpts)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !
    logical add_dip,roc
    !
    !      maxpts=1000
    !
    tstep=0.1*dir
    !
    xray(1)=xi
    yray(1)=yi
    zray(1)=zi
    !
    ts6=tstep/6.
    do n=2,maxpts
        !
        !           step 1
        !
        ax=xray(n-1)
        ay=yray(n-1)
        az=zray(n-1)
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx1,dy1,dz1,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 2
        !
        ax=xray(n-1)+tstep*0.5*dx1
        ay=yray(n-1)+tstep*0.5*dy1
        az=zray(n-1)+tstep*0.5*dz1
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx2,dy2,dz2,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 3
        !
        ax=xray(n-1)+tstep*0.5*dx2
        ay=yray(n-1)+tstep*0.5*dy2
        az=zray(n-1)+tstep*0.5*dz2
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx3,dy3,dz3,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 4
        !
        ax=xray(n-1)+tstep*dx3
        ay=yray(n-1)+tstep*dy3
        az=zray(n-1)+tstep*dz3
        call intpol(bx,by,bz,nx,ny,nz,m,add_dip, &
        ax,ay,az,dx4,dy4,dz4,1.,roc,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       update new position and wavevector
        !
        npts=n
        xray(n)=xray(n-1)+ts6*(dx1+2.*dx2+2.*dx3+dx4)
        yray(n)=yray(n-1)+ts6*(dy1+2.*dy2+2.*dy3+dy4)
        zray(n)=zray(n-1)+ts6*(dz1+2.*dz2+2.*dz3+dz4)
        !
        !       test to see if ray is within selected region
        !
        ar=sqrt(xray(n)**2+yray(n) **2+zray(n)**2)
        if((xray(n).lt.xmin).or.(xray(n).gt.xmax).or. &
            (yray(n).lt.ymin).or.(yray(n).gt.ymax).or. &
            (zray(n).lt.zmin).or.(zray(n).gt.zmax).or. &
            (ar.lt.rearth+1.5*rx).or.(roc)) then
            !
            npts=npts-1
            return
        endif
    enddo
    !
    return
end
