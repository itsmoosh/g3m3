subroutine mak_dip_moon(bx0,by0,bz0,nx,ny,nz,ngrd, &
    mbndry,m,ijzero,numzero,mzero,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      initialize static magnetic field along entire grid
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
        grd_ymin(ngrd),grd_ymax(ngrd), &
        grd_zmin(ngrd),grd_zmax(ngrd)
    dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd), &
        bz0(nx,ny,nz,ngrd)
    integer ijzero(mbndry,3,mzero),numzero(mbndry)
    !
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    !     write(6,*)'in mak dip moon',m,nx,ny,nz,ngrd,mbndry
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !$omp  parallel do
    do k=1,nz
        az=grd_zmin(m)+dz*(k-1)
        z1=(az-zdip)
        !
        do  j=1,ny
            ay=grd_ymin(m)+dy*(j-1)
            y1=(ay-ydip)
            !
            do  i=1,nx
                ax=grd_xmin(m)+dx*(i-1)
                x1=(ax-xdip)
                !
                !        determine magnetic dipole field
                !
                !
                !        determine magnetic dipole field
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
                ar=sqrt(x2+y2+z2)+1.e-8
                !        write(6,*)x2,y2,z2,ar
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
                if(ar.gt.rearth-1.5) then
                    bx0(i,j,k,m)=abx
                    by0(i,j,k,m)=aby
                    bz0(i,j,k,m)=abz
                else
                    bx0(i,j,k,m)=0.
                    by0(i,j,k,m)=0.
                    bz0(i,j,k,m)=0.
                endif
                !
            enddo
        enddo
    enddo
    !
    !
    !
    return
end
