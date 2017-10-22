subroutine contrace(stuff,bx,by,bz,nx,ny,nz,n_grids,box, &
    xmin,xmax,ymin,ymax,zmin,zmax,iside, &
    time,label,nlevs,ncon,add_dip, &
    radstrt,r_inner,nphi,theta1,theta2,ncuts, &
    t,tt,t3,t2,work,mx,my,mz,muvwp2,mz2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
    !        nlevs number of isosurfaces to be plotted	- max 4
    !        ncon  number of contours to be plotted		- max 14
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    real t(mx,my,mz),tt(mx,my,mz),t3(mx,my,mz),t2(mx,my,mz2), &
    work(muvwp2,muvwp2)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    real xray(1000),yray(1000),zray(1000)
    dimension stuff(nx,ny,nz,n_grids), &
    bx(nx,ny,nz),by(nx,ny,nz), &
    bz(nx,ny,nz)
    real eye(3),tlev(6),tcon(14)
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    logical add_dip,roc
    !
    my2=my/2+1
    z22=1.3*mz2
    maxpts=1000
    !
    !     set color table
    !
    call sflush
    call isoclrs
    !
    !      t is the stuff to be plotted covering the region
    !         xmin-xmax,ymin-ymax,zmin-zmax
    !
    !      set up evenly spaced gridding for t to be plotted
    !
    axmax=amin1(xmax,grd_xmax(box)-xdip-.00001)
    aymax=amin1(ymax,grd_ymax(box)-ydip-.00001)
    azmax=amin1(zmax,grd_zmax(box)-zdip-.00001)
    !
    delx=(axmax-xmin)/float(mx-1)
    dely=(aymax-ymin)/float(my-1)
    delz=(azmax-zmin)/float(mz-1)
    !
    ddx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    ddy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    ddz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    !      set decrements have even spacing
    !
    !      delx=amin1(delx,dely,delz)
    !      dely=amin1(delx,dely,delz)
    !      delz=amin1(delx,dely,delz)
    !
    axmax=xmin+delx*(mx-1)
    aymax=ymin+dely*(my-1)
    azmax=zmin+delz*(mz-1)
    !
    !      load t stuff
    !
    do k=1,mz
        az=iside*(zmin+delz*(k-1))
        ak=1.+(az-grd_zmin(box))/ddz
        k1=ak
        k2=k1+1
        dz=ak-k1
        !
        do j=1,my
            ay=iside*(ymin+dely*(j-1))
            aj=1.+(ay-grd_ymin(box))/ddy
            j1=aj
            j2=j1+1
            dy=aj-j1
            !
            !
            do i=1,mx
                ax=xmin+delx*(i-1)
                ai=1.+(ax-grd_xmin(box))/ddx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                t(i,j,k)=stuff(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
                +stuff(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
                +stuff(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
                +stuff(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
                +stuff(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
                +stuff(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
                +stuff(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
                +stuff(i2,j2,k2,box)*(dx)*(dy)*(dz)
                !
                radius=sqrt(ax**2+ay**2+az**2)
                tt(i,j,k)=radius
                !
            enddo
        enddo
    enddo
    !
    !      find max and minimum values in array t
    !
    tmi=t(1,1,1)
    tma=t(1,1,1)
    do k=1,mz
        do j=1,my
            do i=1,mx
                tmi=amin1(t(i,j,k),tmi)
                tma=amax1(t(i,j,k),tma)
            enddo
        enddo
    enddo
    !
    !      set iso levels for surface plots
    !
    dlev=(tma-tmi)/(nlevs+1.)
    do n=1,nlevs
        tlev(n)=tmi+dlev*n
    enddo
    !
    !     set higher resolution contouring
    !
    dcon=(tma-tmi)/(ncon+1.)
    ampl=(abs(tma)+abs(tmi))/2.
    write(magnif,'(1pe9.2)') ampl
    do n=1,ncon
        tcon(n)=tmi+dcon*n
        lind(n)=n+1
        write(llbs(n),'(f4.1)') tcon(n)/ampl
    enddo
    !
    !    set viewport size
    !
    vpl=3200
    vpr=32760
    vpb=3200
    vpt=32760
    !
    !     muvwp2=amax1(mx*1.,my*1.,mz*1.)+2
    !
    !     ***************************************
    !       half plane only start at zdip
    !
    !     zero subset array tt
    !
    
    do k=1,mz2
        do j=1,my
            do i=1,mx
                t2(i,j,k)=0.
            enddo
        enddo
    enddo
    !
    !     initialize eye position
    !
    eye(1)=mx/2.
    eye(2)=-iside*my2*5.
    eye(3)=mz2*3.5
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    call wtstr(.2,.975,label,2,0,0)
    if(iside.ge.0)then
        title='north-dusk'
    else
        title='south-dawn'
    endif
    call wtstr(.5,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='units '//magnif
    call wtstr(.3,.955,title,1,0,0)
    write(title,'(1pe9.2)')tmi
    title='min '//title
    call wtstr(.6,.955,title,1,0,0)
    write(title,'(1pe9.2)')tma
    title='max '//title
    call wtstr(.8,.955,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')zdip
    title='y axis'
    call wtstr(.93,.17,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.12,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.93,.9,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.85,title,1,0,0)
    !
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')zdip
    title='x axis'
    call wtstr(.18,.6,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.16,.55,title,1,0,0)
    !
    !          loading x-y plane - running from y = my/2 to  my
    !
    k=1
    do j=1,my
        do i=1,mx
            t2(i,j,k)=t(i,j,mz2)
        enddo
    enddo
    !
    !     draw contours in x-y plane
    !
    do n=1,ncon
        call gselnt(0)
        call sflush
        tisom=tcon(n)
        call gsplci(n+1)
        call isosrf(t2,mx,mx,my,my,mz2,eye,muvwp2,work,tisom,1, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !     place the earth and map it out
    !
    do k=1,mz2
        do j=1,my
            do i=1,mx
                t2(i,j,k)=tt(i,j,mz2+k-1)
            enddo
        enddo
    enddo
    call gsplci(2)
    tisom=r_inner
    call isosrf(t2,mx,mx,my,my,mz2,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !    draw axes lines
    !
    x1=1
    x2=mx
    y1=1
    y2=my
    z1=1
    z2=mz2
    !
    call set3(0.,1.,0.,1.,x1,x2,y1,y2,z1,z2,eye)
    call sflush
    call gsplci(1)
    call line3(x1,y1,z1,x2,y1,z1)
    call line3(x1,y1,z1,x1,y2,z1)
    call line3(x1,y1,z1,x1,y1,z2)
    call sflush
    !
    !     set up ray trace initial points
    !
    dphi=6.283185/float(nphi)
    deltheta=(theta2-theta1)/float(ncuts)
    if(iside.gt.0)then
        rzmin=zdip
        rzmax=azmax
    else
        rzmin=-azmax
        rzmax=zdip
    endif
    !
    do mm=1,ncuts
        theta=theta1+deltheta*(mm-1)
        sint=sin(theta)
        cost=cos(theta)
        ncol=15-2*(mm-1)
        mcol=5+2*(mm-1)
        !
        do n=1,nphi
            phi=dphi*(n-1)
            cosp=cos(phi)
            sinp=sin(phi)
            !
            x1=radstrt*sint*cosp
            yi=iside*(radstrt*sint*sinp+ydip)
            z1=iside*(radstrt*cost)
            !
            !       dipole space to real space
            !
            xi=x1*cos_tilt+z1*sin_tilt+xdip
            zi=-x1*sin_tilt+z1*cos_tilt+zdip
            !
            dir=-delx
            call rungem(bx,by,bz,nx,ny,nz,box,ddx, &
            xi,yi,zi,xmin,xmax,ymin,ymax,rzmin,rzmax,r_inner, &
            add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            !
            do nn=1,npts
                !
                !        convert spatial cordinate to grid cordinate
                !
                xray(nn)=1.+(xray(nn)-xmin)/delx
                yray(nn)=1.+(iside*yray(nn)-ymin)/dely
                zray(nn)=1.+iside*(zray(nn)-zdip)/delz
            enddo
            !
            !       make color choice and plot the grap
            !
            call sflush
            call gsplci(ncol)
            call curve3(xray,yray,zray,npts)
            !
            dir=delx
            call rungem(bx,by,bz,nx,ny,nz,box,ddx, &
            xi,yi,zi,xmin,xmax,ymin,ymax,rzmin,rzmax,r_inner, &
            add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            !
            do nn=1,npts
                !
                !        convert spatial cordinate to grid cordinate
                !
                xray(nn)=1.+(xray(nn)-xmin)/delx
                yray(nn)=1.+(iside*yray(nn)-ymin)/dely
                zray(nn)=1.+iside*(zray(nn)-zdip)/delz
            enddo
            !
            !       make color choice and plot the grap
            !
            call sflush
            call gsplci(mcol)
            call curve3(xray,yray,zray,npts)
            !
        enddo
    enddo
    !
    !     put color bar
    !
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    call gsfais(1)
    call lbseti('cbl',1)
    call pcsetr('cs',1.25)
    call lblbar(0,0.1,0.9,0.,.1,ncon,1.,.3,lind,0,llbs,ncon,1)
    !
    !     initialize eye position
    !
    eye(1)=mx/2.
    eye(2)=-iside*my*2.5
    eye(3)=mz2*3.5
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    call wtstr(.2,.975,label,2,0,0)
    if(iside.ge.0)then
        title='north-dusk'
    else
        title='south-dawn'
    endif
    call wtstr(.5,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='units '//magnif
    call wtstr(.3,.955,title,1,0,0)
    write(title,'(1pe9.2)')tmi
    title='min '//title
    call wtstr(.6,.955,title,1,0,0)
    write(title,'(1pe9.2)')tma
    title='max '//title
    call wtstr(.8,.955,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')zmin
    title='y axis'
    call wtstr(.98,.17,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.95,.12,title,1,0,0)
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.93,.9,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.85,title,1,0,0)
    !
    !     write(wd1,'(f4.0)')axmax
    !     write(wd2,'(f4.0)')aymin
    !     write(wd3,'(f4.0)')zmin
    !     title='x axis'
    !     call wtstr(.05,.6,title,1,0,0)
    !     title=wd1//','//wd2//','//wd3
    !     call wtstr(.03,.55,title,1,0,0)
    !
    call gsplci(2)
    tisom=r_inner
    call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !
    !     zero subset array tt
    !
    do k=1,mz
        do j=1,my
            do i=1,mx
                t3(i,j,k)=0.
            enddo
        enddo
    enddo
    !
    !     load data ------ right hand side
    !
    !          loading x-y plane - running from y = my/2 to  my
    !
    k=mz2
    do j=1,my
        do i=1,mx
            t3(i,j,k)=t(i,j,k)
        enddo
    enddo
    !
    !     draw contours in x-y plane
    !
    do n=1,ncon
        call gselnt(0)
        call sflush
        tisom=tcon(n)
        call gsplci(n+1)
        call isosrf(t3,mx,mx,my,my,mz,eye,muvwp2,work,tisom,1, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !    draw axes lines
    !
    x1=1
    x2=mx
    y1=1
    y2=my
    z1=1
    z2=mz
    !
    call set3(0.,1.,0.,1.,x1,x2,y1,y2,z1,z2,eye)
    call sflush
    call gsplci(1)
    call line3(x1,y1,z1,x2,y1,z1)
    call line3(x1,y1,z1,x1,y2,z1)
    call line3(x1,y1,z1,x1,y1,z2)
    !     call line3(x1,y1,z22,x2,y1,z22)
    !     call line3(x1,y1,z22,x1,y2,z22)
    !     call line3(x1,y1,z22,x1,y1,z2)
    !     call line3(x1,y2,z22,x2,y2,z22)
    !     call line3(x2,y1,z22,x2,y2,z22)
    call sflush
    !
    do mm=1,ncuts
        theta=theta1+deltheta*(mm-1)
        sint=sin(theta)
        cost=cos(theta)
        ncol=15-2*(mm-1)
        mcol=5+2*(mm-1)
        !
        do n=1,nphi
            phi=dphi*(n-1)
            cosp=cos(phi)
            sinp=sin(phi)
            !
            x1=radstrt*sint*cosp
            yi=iside*(radstrt*sint*sinp+ydip)
            z1=iside*(radstrt*cost)
            !
            !       dipole space to real space
            !
            xi=x1*cos_tilt+z1*sin_tilt+xdip
            zi=-x1*sin_tilt+z1*cos_tilt+zdip
            !
            dir=-delx
            call rungem(bx,by,bz,nx,ny,nz,box,ddx, &
            xi,yi,zi,xmin,xmax,ymin,ymax,-zmax,zmax,r_inner, &
            add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            !
            do nn=1,npts
                !
                !        convert spatial cordinate to grid cordinate
                !
                xray(nn)=1.+(xray(nn)-xmin)/delx
                yray(nn)=1.+(iside*yray(nn)-ymin)/dely
                zray(nn)=1.+(iside*zray(nn)-zmin)/delz
            enddo
            !
            !       make color choice and plot the grap
            !
            call sflush
            call gsplci(ncol)
            call curve3(xray,yray,zray,npts)
            !
            dir=+delx
            call rungem(bx,by,bz,nx,ny,nz,box,ddx, &
            xi,yi,zi,xmin,xmax,ymin,ymax,-zmax,zmax,r_inner, &
            add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            !
            do nn=1,npts
                !
                !        convert spatial cordinate to grid cordinate
                !
                xray(nn)=1.+(xray(nn)-xmin)/delx
                yray(nn)=1.+(iside*yray(nn)-ymin)/dely
                zray(nn)=1.+(iside*zray(nn)-zmin)/delz
            enddo
            !
            !       make color choice and plot the grap
            !
            call sflush
            call gsplci(mcol)
            call curve3(xray,yray,zray,npts)
            !
        enddo
    enddo
    !
    !     put color bar
    !
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    call gsfais(1)
    call lbseti('cbl',1)
    call pcsetr('cs',1.25)
    call lblbar(0,0.1,0.9,0.,.1,ncon,1.,.3,lind,0,llbs,ncon,1)
    return
end
