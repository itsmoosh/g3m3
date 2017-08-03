subroutine trace(bx,by,bz,n1,n2,n3,m,add_dip, &
    xi,yi,zi,dir,np,xf,yf,zf,xx,yy,zz, &
    xmin,xmax,ymin,ymax,zmin,zmax,l,rearth,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !    3 field line tracing subroutines originally by tsyganenko.
    !    modified to trace field lines of bx,by,bz (nx,ny,nz)
    !
    !   traces field line from arbitrary point of space up to the earth
    !   surface or up to model limiting boundary.
    !-------------- input parameters:
    !   xi,yi,zi - gsm coords of initial point (in earth radii),
    !   dir - sign of tracing direction: if dir=1. then antiparallel to
    !     b vector (e.g. from northern to southern conjugate point),
    !     and if dir=-1. then parallel to b.
    !   np - upper estimate of number of steps along the field line
    !     (of the order of several hundreds).
    !-------------- output parameters:
    !   xf,yf,zf - gsm coords of final point
    !   xx,yy,zz - arrays (length np) containing coords of field line points
    !   l - actual number of field line points. if l exceeds np, tracing
    !     terminates, and a warning is displayed
    !
    !                   author: nikolai a. tsyganenko
    !                           institute of physics
    !                           st.-petersburg state university
    !                           stary petergof 198904
    !                           st.-petersburg
    !                           russia
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    dimension bx(n1,n2,n3),by(n1,n2,n3),bz(n1,n2,n3), &
    xx(np),yy(np),zz(np)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    logical roc,add_dip
    !
    err=0.05
    !     ds=0.5*dir
    x=xi
    y=yi
    z=zi
    ds=dir
    ds3=ds/3.
    l=0
    !
    call intpol(bx,by,bz,n1,n2,n3,m,add_dip, &
    x,y,z,r1,r2,r3,ds3,roc,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    if(roc) return
    !
    do l=1,np
        xx(l)=x
        yy(l)=y
        zz(l)=z
        !
        call step(bx,by,bz,n1,n2,n3,m,add_dip,x,y,z,ds,err,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        ar=sqrt((x-xdip)**2+(y-ydip)**2+(z-zdip)**2)
        if((x.lt.xmin).or.(x.gt.xmax).or. &
            (y.lt.ymin).or.(y.gt.ymax).or. &
            (z.lt.zmin).or.(z.gt.zmax).or. &
            (ar.lt.rearth+1.)) then
            xf=x
            yf=y
            zf=z
            xx(l)=x
            yy(l)=y
            zz(l)=z
            err=0.0005
            return
        endif
    	!
    enddo
    !
    !     not enough points
    !
    !     write(6,10)
    l=l-1
    return
end
