!
!	This file contains five subroutines:
!	aurora
!	aurora_bfld
!	auroras
!	aurora_cur
!	aurora_pot
!
subroutine aurora(stuff,nx,ny,nz,box,radstrt,re_equiv,iside, &
    time, save_dat,add_two,label,ncon,write_dat,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension stuff(nx,ny,nz)
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_two,write_dat
    !
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    call isoclrs
    ilab=0
    ioffm=1
    irecmj=1
    irectx=1
    irecmn=1
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    !      theta_range=0.698132  ! 40 degrees
    theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
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
                !
                !          interpolate data to grid point
                !
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
                acur=stuff(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +stuff(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +stuff(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +stuff(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +stuff(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +stuff(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +stuff(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +stuff(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                work(i,j)=acur
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
    title=' : aurora'
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
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
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
    write(title,'(f7.4)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    write(title,'(f7.4)')fmax
    title='fmax'//title
    call wtstr(.25,.98,title,1,0,0)
    call gsplci(1)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
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
    !
    !     save data if required
    !
    if(save_dat)then
        if(iside.gt.0) then
            do j=1,my
                do i=1,mx
                    sdata(i,j)=work(i,j)
                enddo
            enddo
        else
            do j=1,my
                do i=1,mx
                    tdata(i,j)=work(i,j)
                enddo
            enddo
        endif
    endif
    !
    !      output data to bin file if necessary
    !
    if(write_dat)then
        write(8)time,work
    endif
    !
    !
    if(.not.add_two)return
    !
    ilab=0
    if(iside.gt.0)then
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+sdata(i,j)
            enddo
        enddo
    else
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+tdata(i,j)
            enddo
        enddo
    endif
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : total aur'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
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
    write(title,'(f6.3)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    write(title,'(f6.3)')fmax
    title='fmax'//title
    call wtstr(.25,.98,title,1,0,0)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
    !
    !     ilab=1
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    return
end
!
!
!	****************************************
!
!
subroutine aurora_bfld(bx,by,bz,nx,ny,nz,box,rx, &
    xmin,xmax,ymin,ymax,zmin,zmax,iside, &
    add_dip,radstrt,re_equiv,r_inner,time,label,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !       this subroutine seeks to determine open field
    !             line positions
    !
    common /space/work(91,91),clat(91,91), &
    xray(1000),yray(1000),zray(1000)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical add_dip,roc
    !
    !      dimension for plotted array
    !
    maxpts=1000
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.6108652
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    rx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
        enddo
    enddo
    !
    do j=1,my
        do i=1,mx
            !
            closed=0.
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
                !
                dir=-rx
                call rungem(bx,by,bz,nx,ny,nz,box,rx, &
                ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax,r_inner, &
                add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                !
                !        find the radial distance of end point
                !
                ar=sqrt((xray(npts)-xdip)**2+(yray(npts)-ydip)**2 &
                +(zray(npts)-zdip)**2)
                if(ar.le.r_inner+2.*rx)closed=closed+1.
                !
                dir=rx
                call rungem(bx,by,bz,nx,ny,nz,box,rx, &
                ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax,r_inner, &
                add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                !
                !        find the radial distance of end point
                !
                ar=sqrt(xray(npts) **2+yray(npts) **2+zray(npts)**2)
                if(ar.le.r_inner+2.*rx)closed=closed+1.
                work(i,j)=closed
            endif
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
    title=' : separatrix'
    call wtstr(.3,.975,label,2,0,0)
    call wtstr(.5,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    call conrec(work,mx,mx,my,1.4,1.6,.1,0,-1,-1012)
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
    !
    return
end
!
!
!	****************************************
!
!
subroutine auroras(press,rho,bsx,bsy,bsz, &
    nx,ny,nz,n_grids,box,radstrt,re_equiv,r_inner,iside, &
    time,save_dat,add_dip,label,ncon,write_dat, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make a close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !    
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91),chot(91,91)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    dimension press(nx,ny,nz),rho(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz)
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_dip,write_dat
    !
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    !
    ilab=0
    add_dip=.false.
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.698132  ! 40 degrees
    !      theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
            chot(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    xmin=grd_xmin(box)+1.
    xmax=grd_xmax(box)-1.
    ymin=grd_ymin(box)+1.
    ymax=grd_ymax(box)-1.
    zmin=grd_zmin(box)+1.
    zmax=grd_zmax(box)-1.
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
                !
                !          interpolate data to grid point
                !
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
                !         trace fld line and find max energy along it
                !
                ergies=0.
                rhod=0.
                dir=-delx
                call rungea(bsx,bsy,bsz,press,rho,nx,ny,nz, &
                box,ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax, &
                add_dip,ergies,rhod,maxpts,dir, &
                delx,dely,delz,r_inner,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                dir=+delx
                call rungea(bsx,bsy,bsz,press,rho,nx,ny,nz, &
                box,ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax, &
                add_dip,ergies,rhod,maxpts,dir, &
                delx,dely,delz,r_inner,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    			!
                chot(i,j)=ergies
            endif
        enddo
    enddo
    !
    !
    !     initialize viewport and frame headings
    !
    !     plot intensity = current*energy
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,chot(i,j))
            fmin=amin1(fmin,chot(i,j))
        enddo
    enddo
    !
    call frame
    call gselnt(0)
    call gsplci(18)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : intensity'
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
    write(title,'(f7.4)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    title='fmax'//title
    write(title,'(f7.4)')fmax
    call wtstr(.25,.98,title,1,0,0)
    if(fmin*fmax.lt.0.)then
        cmax=amax1(fmax,abs(fmin))
        cmin=-cmax
    else
        cmax=fmax
        cmin=fmin
    endif
    clev=(cmax-cmin)/(ncon+1.)
    !     if(.not.start) call cnrccf(chot,mx,mx,my,cmin,
    !    +             cmax,clev,0,-1,-1012,2.5,0,1)
    irecmn=18
    call conrec(chot,mx,mx,my,cmin, &
    cmax,clev,0,-1,-1012)
    !
    !     plot constant latitude circles
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
    call conrec(clat,mx,mx,my,50.,90.,10.,0,-1,-1012)
    !
    !     save data if required
    !
    if(write_dat)then
        write(8)time,chot
    endif
    return
end
!
!
!	****************************************
!
!
subroutine aurora_cur(stuff,nx,ny,nz,box,radstrt,re_equiv, &
    r_inner,iside,time,save_dat,add_two,label,ncon,write_dat, &
    b_equiv,planet_rad,tot_cur,peak_cur,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91)
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension stuff(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_two,write_dat
    !
    !     for no line labeling set ilab to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    call isoclrs
    ilab=0
    ioffm=1
    irecmj=1
    irectx=1
    irecmn=1
    !
    !       conversion of current simulation units to ua/box^2 and to ma
    !                 map to earth's surface so a ~ r**3
    !
    scale_cur= 10. * b_equiv *((radstrt*re_equiv)**3)/ &
    (4.* 3.1416*re_equiv*planet_rad)
    !       write(*,*)b_equiv,radstrt,re_equiv,planet_rad
    !       write(*,*)' current density in',scale_cur,' ua/box^2'
    scale_amps=((planet_rad/1000.)**2)*scale_cur
    !       write(*,*)'totcur scale',scale_amps
    amp_up=0.
    amp_down=0.
    peak_cur_up=0.
    peak_cur_down=0.
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    !      theta_range=0.698132  ! 40 degrees
    theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    rads=degrees/57.3
    del_rads=rads/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
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
            alat=dlat*del_theta*theta_equiv  ! equiv theta at sample point
            !
            theta=dlat*del_rads  ! theta in ionosphere
            ascale=del_rads*sin(theta)/(dlat+1.)
            !
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
                !
                !          interpolate data to grid point
                !
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
                acur=stuff(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +stuff(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +stuff(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +stuff(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +stuff(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +stuff(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +stuff(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +stuff(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                acur=acur*scale_amps  !in physical unit ua/box**2
                work(i,j)=acur
                totcur=ascale*acur
                if(acur.gt.0)then
                    amp_up=amp_up+totcur
                    peak_cur_up=amax1(peak_cur_up,acur)
                else
                    amp_down=amp_down+totcur
                    peak_cur_down=amin1(peak_cur_down,acur)
                endif
            endif
            !
        enddo
    enddo
    !
    tot_cur=(amp_up-amp_down)/2.
    peak_cur=amax1(peak_cur_up,-peak_cur_down)
    write(6,*)'totcur',amp_up,amp_down,tot_cur
    write(6,*)'peak current',peak_cur_up,peak_cur_down,peak_cur
    !
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : aurora'
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
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
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
    call wtstr(.3,.98,title,1,0,0)
    !
    write(title,'(f7.3)')amp_up
    title='tot_up'//title
    call wtstr(.1,.96,title,1,0,0)
    write(title,'(f7.3)')amp_down
    title='tot_down'//title
    call wtstr(.45,.96,title,1,0,0)
    !
    call gsplci(1)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
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
    !
    !     save data if required
    !
    if(save_dat)then
        if(iside.gt.0) then
            do j=1,my
                do i=1,mx
                    sdata(i,j)=work(i,j)
                enddo
            enddo
        else
            do j=1,my
                do i=1,mx
                    tdata(i,j)=work(i,j)
                enddo
            enddo
        endif
    endif
    !
    !      output data to bin file if necessary
    !
    if(write_dat)then
        write(8)time,work
    endif
    !
    !
    if(.not.add_two)return
    !
    ilab=0
    if(iside.gt.0)then
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+sdata(i,j)
            enddo
        enddo
    else
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+tdata(i,j)
            enddo
        enddo
    endif
    !
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : total aur'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.85,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
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
    write(title,'(f6.3)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    write(title,'(f6.3)')fmax
    title='fmax'//title
    call wtstr(.25,.98,title,1,0,0)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
    !
    !     ilab=1
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    return
end
!
!
!	****************************************
!
!
subroutine aurora_pot(efldx,efldy,efldz,bsx,bsy,bsz, &
    nx,ny,nz,n_grids,box,radstrt,re_equiv,iside,r_inner, &
    time,save_dat,add_dip,label,ncon,write_dat, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make a close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91),chot(91,91)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_dip,write_dat
    !
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    !
    ilab=0
    add_dip=.false.
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    !      theta_range=0.698132  ! 40 degrees
    theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
            chot(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    xmin=grd_xmin(box)+1.
    xmax=grd_xmax(box)-1.
    ymin=grd_ymin(box)+1.
    ymax=grd_ymax(box)-1.
    zmin=grd_zmin(box)+1.
    zmax=grd_zmax(box)-1.
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
                !
                !         trace fld line and find max energy along it
                !
                dir=-delx
                call rungeb(efldx,efldy,efldz,bsx,bsy,bsz,nx,ny,nz, &
                n_grids,box,delx,ax,ay,az,xmin,xmax,ymin,ymax, &
                zmin,zmax,add_dip,pot1,maxpts,dir,r_inner, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                dir=+delx
                call rungeb(efldx,efldy,efldz,bsx,bsy,bsz,nx,ny,nz, &
                n_grids,box,delx,ax,ay,az,xmin,xmax,ymin,ymax, &
                zmin,zmax,add_dip,pot2,maxpts,dir,r_inner, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    			!
                chot(i,j)=iside*(pot1-pot2)
            endif
            !
        enddo
    enddo
    !
    !     initialize viewport and frame headings
    !
    !     plot intensity = current*energy
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,chot(i,j))
            fmin=amin1(fmin,chot(i,j))
        enddo
    enddo
    !
    call frame
    call gselnt(0)
    call gsplci(18)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : intensity'
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
    write(title,'(f7.4)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    title='fmax'//title
    write(title,'(f7.4)')fmax
    call wtstr(.25,.98,title,1,0,0)
    !     if(fmin*fmax.lt.0.)then
    !       cmax=amax1(fmax,abs(fmin))
    !       cmin=-cmax
    !     else
    !       cmax=0.5*fmax
    !       cmin=fmin
    !     endif
    cmax=0.15    !20 kev
    cmin=-0.15
    clev=(cmax-cmin)/(ncon+1.)
    !     if(.not.start) call cnrccf(chot,mx,mx,my,cmin,
    !    +             cmax,clev,0,-1,-1012,2.5,0,1)
    irecmn=18
    call conrec(chot,mx,mx,my,cmin, &
    cmax,clev,0,-1,-1012)
    !
    !     plot constant latitude circles
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
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    !
    !     save data if required
    !
    if(write_dat)then
        write(8)time,chot
    endif
    return
end
