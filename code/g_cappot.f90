subroutine cappot(chrg,pott,nx,ny,nz,ngrd,m,radstrt, &
    re_equiv,v_equiv,b_equiv,time,write_dat, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make potential drop
    !         over auroral oval, vx,vy,vz, assumed to
    !         be the inductive electric field
    !      answer is in kv: does both hemispheres at once to save work
    !      charge is the density over a sample region around
    !       the earth
    !
    common /space/charge(31,31,31),poten(31,31,31), &
    pot(61,61),clat(61,61)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension chrg(nx,ny,nz),pott(nx,ny,nz)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical write_dat
    !
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    !
    ilab=0
    !
    !      dimension for plotted array
    !
    mx=61
    my=61
    my2=my/2+1
    mx2=mx/2+1
    amx2=mx2-1.
    !
    msx=31
    msy=31
    msz=31
    ijk=(msx-1)/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !      re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.6108652
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      b_equiv=65.75
    !       b is in nt and v in km/s so total si conversion is 1e-6
    !       v_equiv=1e-3
    !       d is still in kilometers so the
    !        potential is in kilovolts
    !
    d_equiv=6371.*re_equiv
    scale=(1.e-6)*b_equiv*v_equiv*d_equiv
    !
    delx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dely=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    delz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    ncx=1.-grd_xmin(m)/delx
    ncy=1.-grd_ymin(m)/dely
    ncz=1.-grd_zmin(m)/delz
    !
    !     find the charge density equal to the div of vxb
    !
    do ks=1,msz
        k=(ks-ijk)+ncz
        do js=1,msy
            j=(js-ijk)+ncy
            do is=1,msx
                i=(is-ijk)+ncx
                charge(is,js,ks)=chrg(i,j,k)
                poten(is,js,ks)=pott(i,j,k)
            enddo
        enddo
    enddo
    !
    !      initialize work array
    !
    do mm=1,2
        iside=(-1)**(mm-1)
        do j=1,my
            do i=1,mx
                pot(i,j)=0.0
            enddo
        enddo
        !
        do j=1,my
            do i=1,mx
                !
                !         find equivalent latitude
                !
                dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
                dlat=amax1(dlat,0.0000001)
                alat=dlat*del_theta*theta_equiv
                if(alat.le.1.55) then
                    cost=cos(alat)
                    sint=sin(alat)
                    !
                    !         find equivalent longitude
                    !
                    cosp=(i-mx2)/dlat
                    sinp=(j-my2)/dlat
                    !
                    !          find position on grid
                    !
                    z1=iside*(radstrt*cost)
                    x1=radstrt*sint*cosp
                    !
                    !       dipole space to real space
                    !
                    ax=x1*cos_tilt+z1*sin_tilt+xdip
                    az=-x1*sin_tilt+z1*cos_tilt+zdip
                    ay=(radstrt*sint*sinp+ydip)
                    ar=sqrt((ax-xdip)**2+(ay-ydip)**2+(az-zdip)**2) &
                    +0.0000001
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
                    apot=chrg(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                    +chrg(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                    +chrg(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                    +chrg(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                    +chrg(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                    +chrg(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                    +chrg(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                    +chrg(i2,j2,k2)*(dx)*(dy)*(dz)
                    !
                    pot(i,j)=-apot
                endif
            enddo
        enddo
        !
        !
        call frame
        call gselnt(0)
        call gsplci(1)
        call agseti('frame.',2)
        call agseti('set.',-1)
        !
        if(iside.gt.0) title= 'nth cap pot'
        if(iside.le.0) title= 'sth cap pot'
        call wtstr(.5,.975,title,2,0,0)
        write(title,'(f7.3)')time
        title='ut = '//title
        call wtstr(.8,.975,title,2,0,0)
        !
        ncon=21
        fmax=pot(1,1)
        fmin=fmax
        do j=1,my
            do i=1,mx
                fmax=amax1(fmax,pot(i,j))
                fmin=amin1(fmin,pot(i,j))
            enddo
        enddo
        !
        if (fmin.ge.0.0)then
            cmin=1.05*fmin
        else
            cmin=0.95*fmin
        endif
        !
        if(fmax.ge.0.0) then
            cmax=0.95*fmax
        else
            cmax=1.05*fmax
        endif
        finc=(cmax-cmin)/float(ncon-1)
        write(title,'(f7.2)')fmin
        title='fmin'//title
        call wtstr(.1,.98,title,1,0,0)
        write(title,'(f7.2)')fmax
        title='fmax'//title
        call wtstr(.25,.98,title,1,0,0)
        call gsplci(1)
        call conrec(pot,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
        !
        !     make constant latitude circles
        !
        dd=degrees/float(mx2-1)
        fmax=0.
        do j=1,my
            do i=1,mx
                dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
                clat(i,j)=90.-dlat*dd
                fmax=amax1(fmax,clat(i,j))
            enddo
        enddo
        !
        !     ilab=1
        call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    enddo
    !
    !
    !      output data to bin file if necessary
    !
    if(write_dat)then
        write(9)time,scale,cos_tilt,sin_tilt,radstrt, &
        b_equiv,v_equiv,re_equiv,rx
        write(9)charge
    endif
    return
end
