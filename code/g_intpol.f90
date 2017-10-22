subroutine intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
    ax,ay,az,r1,r2,r3,ds3,roc,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     interpolate the magnetic to from nearest grid points
    !     to the actual point ax,ay,az and then calculates the
    !     direction for the increments in r1,r2,r3
    !
    integer box
    dimension bx(n1,n2,n3),by(n1,n2,n3),bz(n1,n2,n3)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    logical roc,add_dip
    !
    !       altered to make plots below equator assuming symmetry
    !
    ddx=(grd_xmax(box)-grd_xmin(box))/(n1-1.)
    ddy=(grd_ymax(box)-grd_ymin(box))/(n2-1.)
    ddz=(grd_zmax(box)-grd_zmin(box))/(n3-1.)
    !
    ak=1.+(az-grd_zmin(box))/ddz
    k1=ak
    k2=k1+1
    dz=ak-k1
    !
    aj=1.+(ay-grd_ymin(box))/ddy
    j1=aj
    j2=j1+1
    dy=aj-j1
    !
    ai=1.+(ax-grd_xmin(box))/ddx
    i1=ai
    i2=i1+1
    dx=ai-i1
    !
    abx=0.
    aby=0.
    abz=0.
    !
    if(add_dip)call dipole(abx,aby,abz,ax,ay,az)
    !
    abx=abx+bx(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
    +bx(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
    +bx(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
    +bx(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
    +bx(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
    +bx(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
    +bx(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
    +bx(i2,j2,k2)*(dx)*(dy)*(dz)
    !
    aby=aby+by(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
    +by(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
    +by(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
    +by(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
    +by(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
    +by(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
    +by(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
    +by(i2,j2,k2)*(dx)*(dy)*(dz)
    !
    abz=abz+bz(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
    +bz(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
    +bz(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
    +bz(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
    +bz(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
    +bz(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
    +bz(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
    +bz(i2,j2,k2)*(dx)*(dy)*(dz)
    !
    sqb=sqrt(abx**2+aby**2+abz**2)
    if(sqb.lt.1.e-4) then
        roc=.true.
        return
    endif
    !
    roc=.false.
    b=ds3/sqb
    r1=abx*b
    r2=aby*b
    r3=abz*b
    !
    return
end
