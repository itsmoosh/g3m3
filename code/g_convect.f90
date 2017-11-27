subroutine convect(vx,vy,vz,nx,ny,nz,box,radstrt, &
    re_equiv,iside,time,label,write_dat,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make convection pattern
    !         over auroral oval
    !      sdata tdata included to allow passing of data between
    !         successive calls of aurora
	!
    common /space/sdata(61,61),tdata(61,61),work(31,31), &
    avx(31,31),avy(31,31)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    integer box
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz),spv(2)
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical write_dat
	!
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
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
    mx=31
    my=31
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latitude at radstrt
    !
    !       re_equiv=0.84
	!
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.6108652
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
            avx(i,j)=0.0
            avy(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
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
                ar=sqrt((ax-xdip)**2+(ay-ydip)**2+(az-zdip)**2)+0.0000001
                !
                !          interpolate data to grid point
                !
                ak=1.+(az-grd_zmin(box))/delz
                k1=ak
                k2=k1+1
                dz=ak-k1
                !
                aj=1.+(ay-grd_ymin(box))/dely
                j1=aj
                j2=j1+1
                dy=aj-j1
                !
                ai=1.+(ax-grd_xmin(box))/delx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                bvx=vx(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +vx(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +vx(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +vx(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +vx(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +vx(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +vx(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +vx(i2,j2,k2)*(dx)*(dy)*(dz)
                bvy=vy(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +vy(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +vy(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +vy(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +vy(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +vy(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +vy(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +vy(i2,j2,k2)*(dx)*(dy)*(dz)
                bvz=vz(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +vz(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +vz(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +vz(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +vz(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +vz(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +vz(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +vz(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                !          find tangential components
                !
                cost=iside*(az-zdip)/ar
                bvz=iside*bvz
                sint=sqrt((ax-xdip)**2+(ay-ydip)**2)/ar
                cosp=(ax-xdip)/(ar*sint+0.000001)
                sinp=(ay-ydip)/(ar*sint+0.000001)
                vr=bvx*sint*cosp+bvy*sint*sinp+bvz*cost
                vtheta=bvx*cost*cosp+bvy*cost*sinp-bvz*sint
                vphi=-bvx*sinp+bvy*cosp
                !
                !        restructure magnetic field without radial velocity component
                !
                vr=0.
                bvx=vr*sint*cosp+vtheta*cost*cosp-vphi*sinp
                bvy=vr*sint*sinp+vtheta*cost*sinp+vphi*cosp
                bvz=vr*cost-vtheta*sint
                avx(i,j)=vtheta*cosp-vphi*sinp
                avy(i,j)=vtheta*sinp+vphi*cosp
                !         avx(i,j)=bvx
                !         avy(i,j)=bvy
            endif
            !
        enddo
    enddo
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    !     set color table
    !
    call sflush
    call isorb
    !
    title=' : convect'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    xcmin=0.
    xcmax=0.
    ycmin=0.
    ycmax=0.
    do j=1,my
        do i=1,mx
            xcmav=amax1(avx(i,j),xcmav)
            xcmin=amin1(avx(i,j),xcmin)
            ycmav=amax1(avy(i,j),ycmav)
            ycmin=amin1(avy(i,j),ycmin)
        enddo
    enddo
    !
    if(abs(xcmin).gt.xcmax)xcmax=abs(xcmin)
    if(abs(ycmin).gt.ycmax)ycmax=abs(ycmin)
    cmax=amax1(xcmax,ycmax)
    if(cmax.gt.1.e-4)then
        title='maxspeed '
        call wtstr(.3,.96,title,1,0,0)
        write(title,'(f7.4)')cmax
        call wtstr(.4,.96,title,1,0,0)
        !      call ezvec(avx,avy,mx,my)
        !      call velvct(avx,mx,avy,mx,mx,my,0.1*cmax,0.,0,100,0,spv)
        call vvinit(avx,mx,avy,mx,dummy,0,mx,my,dummy,0)
        call vvsetr('vrl',0.20)
        call vvectr(avx,avy,dummy,dummy,dummy,dummy)
    endif
    !
    !     make constant latitude circles
    !
    dd=degrees/float(mx2-1)
    fmax=0.
    do j=1,my
        do i=1,mx
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            work(i,j)=90.-dlat*dd
            fmax=amax1(fmax,work(i,j))
        enddo
    enddo
    !
    !     ilab=1
    call conrec(work,mx,mx,my,60.,90.,10.,0,-1,-1012)
    !
    !      output data to data file if necessary
    !
    if(write_dat)then
        write(recdt_f,*) time,avx,avy
    endif
    return
end
