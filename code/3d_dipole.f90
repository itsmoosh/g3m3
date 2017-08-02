subroutine dipole(abx,aby,abz,ax,ay,az)
    !
    !     calculates magnetic field for a dipole
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    x1=(ax-xdip)
    y1=(ay-ydip)
    z1=(az-zdip)
    !
    !       rotate coordinates for planet motion
    !
    xr=x1*cos_rot+y1*sin_rot
    yr=-x1*sin_rot+y1*cos_rot
    !
    !
    !       tilt space to dipole space
    !
    xp=xr*cos_tilt-z1*sin_tilt
    yp=yr
    zp=xr*sin_tilt+z1*cos_tilt
    !
    x2=xp**2
    y2=yp**2
    z2=zp**2
    ar=sqrt(x2+y2+z2)
    !
    !        cartesian equivalent
    !
    bmag=-b0/ar**5
    dbx=-3.*bmag*xp*zp
    dby=-3.*bmag*yp*zp
    dbz=bmag*(x2+y2-2.*z2)
    !
    !        tilt b field back to coordinate space
    !
    rbx=dbx*cos_tilt+dbz*sin_tilt
    rby=dby
    rbz=-dbx*sin_tilt+dbz*cos_tilt
    !
    !         rotate b
    !
    abx=rbx*cos_rot-rby*sin_rot
    aby=rbx*sin_rot+rby*cos_rot
    abz=rbz
    !
    return
end
