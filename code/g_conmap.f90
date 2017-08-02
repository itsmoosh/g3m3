subroutine conmap(stuff,vx,vy,vz,bx,by,bz,nx,ny,nz,ngrd,mm,m, &
    xcraft,ncraft,re_equiv,rearth, &
    xmin,xmax,ymin,ymax,zmin,zmax, &
    time,label,nlevs,ncon,strtch,add_dip,ivel,kht, &
    tx,ty,tz,t,tt,t2,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
    !        nlevs number of isossurfaces to be plotted - max 4
    !        ncon  number of contours to be plotted     - max 14
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    dimension t(mx,my,mz),tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz), &
    tt(mx,my,mz),t2(mx,my,mz2),work(muvwp2,muvwp2)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    real xray(1000),yray(1000),zray(1000), &
    xray1(1000),yray1(1000),zray1(1000)
    dimension stuff(nx,ny,nz,ngrd),vx(nx,ny,nz),xcraft(4,ncraft), &
    vy(nx,ny,nz),vz(nx,ny,nz), &
    bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    real eye(3),tlev(6),tcon(14)
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    logical add_dip
    !
    !      dimension for plotted array
    !
    mx=41
    my=41
    mz=21
    !
    mz2=mz/2+1
    my2=my/2+1
    maxpts=1000
    !
    !      effective step size between points - slight less than
    !         unity to ensure arrows do not go out of bounds
    !
    hi=.99
    hj=.99
    hk=.99
    !
    !     al = length of arrow in xy-plane
    !     stretch= enlargement of arrows
    !
    al=0.6*hi
    al2=0.2
    stretch=2.0*strtch
    vmin=0.00001
    !
    !      angle of arrow
    !
    beta=30./180.*3.1415926
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
    !
    axmax=amin1(xmax,grd_xmax(m)-.0001)
    aymax=amin1(ymax,grd_ymax(m)-.0001)
    azmax=amin1(zmax,grd_zmax(m)-.0001)
    !
    rx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    !
    delx=(axmax-xmin)/float(mx-1)
    dely=(aymax-ymin)/float(my-1)
    delz=(azmax-zmin)/float(mz-1)
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
    !      height=zmin+delz*0.7*(mz-1)
    height=delz*0.2*(mz-1)
    !
    !      load t stuff
    !
    call filldatav(stuff,vx,vy,vz,nx,ny,nz,ngrd,mm,m, &
    t,tt,tx,ty,tz,mx,my,mz,xmin,ymin,zmin, &
    delx,dely,delz,vm, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    do k=1,mz
        az=zmin+delz*(k-1)
        haz=az-height
        do j=1,my
            ay=ymin+dely*(j-1)
            do i=1,mx
                ax=xmin+delx*(i-1)
                radius=sqrt(ax**2+ay**2+haz**2)
                tt(i,j,k)=radius
                !
            enddo
        enddo
    enddo
    !
    !      find max and minimum values in array t
    !
    tmi=0.
    tma=0.
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
    !
    vpl=3200
    vpr=32760
    vpb=3200
    vpt=32760
    !
    !
    !     initialize eye position
    !
    eye(1)=mx/2.
    eye(2)=-my*2.5
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
    title=' : dawn dusk'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
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
    write(title,'(1pe9.2)')vm
    title='vm ='//title
    call wtstr(.7,.935,title,1,0,0)
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
    call gsplci(2)
    tisom=rearth
    call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !
    !     zero subset array tt
    !
    do k=1,mz
        do j=1,my
            do i=1,mx
                tt(i,j,k)=0.
            enddo
        enddo
    enddo
    !
    !     load data ------ right hand side
    !
    !          loading x-y plane - running from y = my/2 to  my
    !
    k=1
    kk=mz2
    do j=1,my
        do i=1,mx
            tt(i,j,k)=t(i,j,kk)
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
        call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,1, &
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
    z22=1.+mz/2+(height)/delz
    !
    call set3(0.,1.,0.,1.,x1,x2,y1,y2,z1,z2,eye)
    call sflush
    call gsplci(1)
    call line3(x1,y1,z1,x2,y1,z1)
    call line3(x1,y1,z1,x1,y2,z1)
    call line3(x1,y1,z1,x1,y1,z2)
    call line3(x1,y1,z22,x2,y1,z22)
    call line3(x1,y1,z22,x1,y2,z22)
    call line3(x1,y1,z22,x1,y1,z2)
    call line3(x1,y2,z22,x2,y2,z22)
    call line3(x2,y1,z22,x2,y2,z22)
    call sflush
    !
    !     draw points
    !
    ncol=3
    kk=2
    k=mz2+kht
    jskip=3
    iskip=4
    z1=kk*hk+0.05
    do j=1,my,jskip
        y1=j*hj+0.05
        do i=1,mx,iskip
            x1=i*hi+0.05
            !
            !        find relative size of arrow and skip out if too small
            !
            vmag=sqrt(tx(i,j,k)**2+ty(i,j,k)**2+tz(i,j,k)**2)
            call sflush
            call gsplci(2)
            if(vmag*strtch.ge.0.10*vm) then
                !
                !        arrow length
                !
                dx=iskip*hi*tx(i,j,k)/vm
                dy=jskip*hj*ty(i,j,k)/vm
                dz=iskip*hk*tz(i,j,k)/vm
                !
                dx=stretch*dx
                dy=stretch*dy
                dz=stretch*dz
                !
                x2=x1+dx
                y2=y1+dy
                z2=z1+dz
                !
                call sflush
                !
                call gsplci(1)
                !
                !        make simpler by not adding arrows in this panel
                !
                call line3(x1,y1,z1,x2,y2,z2)
                !
                !        test if you can fit an arrow head to the line
                !
                isize=0
                if((abs(dx).ge.al2).or.(abs(dy).ge.al2) &
                    .or.(abs(dz).ge.al2)) then
                    !
                    !      draw arrows - set color for arror heads
                    !
                    call sflush
                    if(((ivel.eq.1).and.(dx.ge.0.0)).or. &
                        ((ivel.eq.2).and.(dy.ge.0.0)).or. &
                        ((ivel.eq.3).and.(dz.ge.0.0)))then
                        call gsplci(15)
                    else
                        call gsplci(5)
                    endif
                    !
                    call arwxy(dx,dy,al,beta,adx1,ady1,adx2,ady2)
                    call arwxy(dy,dz,al,beta,ady1,adz1,ady2,adz2)
                    !
                    call line3(x2,y2,z2,x2+adx1,y2+ady1,z2+adz1)
                    call line3(x2,y2,z2,x2+adx2,y2+ady2,z2+adz2)
                    isize=1
                endif
            endif
            !
            !       draw magnetic field line
            !
            l=0
            dir=0.25*delx
            xi=xmin+delx*(i-1)
            yi=ymin+dely*(j-1)
            zi=zmin+delz*(k-1)
            call trace(bx,by,bz,nx,ny,nz,m,add_dip, &
            xi,yi,zi,dir,maxpts,xf,yf,zf, &
            xray,yray,zray,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l,rearth,ngrd, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            !
            rad1=sqrt(xi**2+yi**2+zi**2)-2.
            rad2=sqrt(xray(l)**2+yray(l)**2+zray(l)**2)-2.
            !
            do nn=1,l
                xray(nn)=1.+(xray(nn)-xmin)/delx
                yray(nn)=1.+(yray(nn)-ymin)/dely
                zray(nn)=1.+(zray(nn)+height-zmin)/delz
            enddo
            !
            !       call curve3(xray,yray,zray,l)
            !
            !       go the other direction just to make sure a ray plotted
            !
            l1=0
            dir=-0.25*delx
            call trace(bx,by,bz,nx,ny,nz,m,add_dip, &
            xi,yi,zi,dir,maxpts,xf1,yf1,zf1, &
            xray1,yray1,zray1,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l1,rearth,ngrd, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            !
            rad3=sqrt(xray1(l1)**2+yray1(l1)**2+zray1(l1)**2)-2.
            do nn=1,l1
                !
                !        convert spatial cordinate to grid cordinate
                !
                xray1(nn)=1.+(xray1(nn)-xmin)/delx
                yray1(nn)=1.+(yray1(nn)-ymin)/dely
                zray1(nn)=1.+(zray1(nn)+height-zmin)/delz
            enddo
            !
            call sflush
            if((xray(l).ge.mx-1.).and.(xray1(l1).ge.mx-1.))then
                call gsplci(12)
                call curve3(xray,yray,zray,l)
                call curve3(xray1,yray1,zray1,l1)
            else  if((rad2.le.rearth).and.(rad3.le.rearth))then
                call gsplci(5)
                call curve3(xray,yray,zray,l)
                call curve3(xray1,yray1,zray1,l1)
            else if((rad1.le.rearth).or.(rad2.le.rearth) &
                .or.(rad3.le.rearth))then
                call gsplci(7)
                call curve3(xray,yray,zray,l)
                call curve3(xray1,yray1,zray1,l1)
            else if(xray(1).le.3.)then
                call gsplci(15)
                call curve3(xray,yray,zray,l)
                call curve3(xray1,yray1,zray1,l1)
            endif
            !
            !
        enddo
    enddo
    !
    !      draw spacecraft - x-r coords and 3-d coords
    !
    do n=1,ncraft
        ax=xcraft(1,n)/re_equiv
        ay=xcraft(2,n)/re_equiv
        az=xcraft(3,n)/re_equiv
        aside=sign(1.,ay)
        ar=sqrt(ay**2+az**2)
        ar=aside*ar
        if( (ax.ge.xmin).and.(ax.le.xmax).and. &
            (ay.ge.ymin).and.(ay.le.ymax).and. &
            (az.ge.zmin).and.(az.le.zmax) )then
            x1=1.+(ax-0.5-xmin)/delx
            x2=1.+(ax+0.5-xmin)/delx
            y1=1.+(ay-0.5-ymin)/dely
            y2=1.+(ay+0.5-ymin)/dely
            z1=1.+(az+height-0.5-zmin)/delz
            z2=1.+(az+height+0.5-zmin)/delz
            r1=1.+(ar-0.5-ymin)/dely
            r2=1.+(ar+0.5-ymin)/dely
            !
            !        2-d projection
            !
            call sflush
            call gsplci(15)
            call line3(x1,r1,1.,x2,r1,1.)
            call line3(x1,r1,1.,x1,r2,1.)
            call line3(x2,r2,1.,x1,r2,1.)
            call line3(x2,r2,1.,x2,r1,1.)
            !
            !        3-d position
            !
            call sflush
            call gsplci(3)
            call line3(x1,y1,z1,x2,y1,z1)
            call line3(x1,y1,z1,x1,y2,z1)
            call line3(x1,y1,z1,x1,y1,z2)
            call line3(x2,y2,z1,x2,y1,z1)
            call line3(x2,y2,z1,x1,y2,z1)
            !
            call line3(x2,y2,z2,x2,y2,z1)
            call line3(x2,y2,z2,x1,y2,z2)
            call line3(x2,y2,z2,x2,y1,z2)
            call line3(x1,y1,z2,x1,y2,z2)
            call line3(x1,y1,z2,x2,y1,z2)
            !
            call line3(x1,y2,z1,x1,y2,z2)
            call line3(x2,y1,z1,x2,y1,z2)
        endif
        !
        dir=-0.25*delx
        call rungem(bx,by,bz,nx,ny,nz,m,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,rearth, &
        add_dip,xray,yray,zray,maxpts,npts,dir,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do nn=1,npts
            !
            !        convert spatial cordinate to grid cordinate
            !
            xray1(nn)=1.+(xray(nn)-xmin)/delx
            yray1(nn)=1.+(yray(nn)-ymin)/dely
            zray1(nn)=1.+(zray(nn)+height-zmin)/delz
        enddo
        !
        !       make color choice and plot the grap
        !           3d plot
        !
        call sflush
        call gsplci(1)
        call curve3(xray1,yray1,zray1,npts)
        !
        !          2d plot
        !
        !        side=(-1.)**n
        do nn=1,npts
            !
            !        convert spatial cordinate to grid cordinate
            !
            ar=sqrt((yray(nn)-ydip)**2+(zray(nn)-zdip)**2)
            ar=sign(ar,yray(nn))
            !        ar=ar*side
            xray1(nn)=1.+(xray(nn)-xmin)/delx
            yray1(nn)=1.+(ar-ymin)/dely
            zray1(nn)=1.
        enddo
        !        call sflush
        !        call gsplci(1)
        call curve3(xray1,yray1,zray1,npts)
        !
        !       go the other direction
        dir=+0.25*delx
        call rungem(bx,by,bz,nx,ny,nz,m,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,rearth, &
        add_dip,xray,yray,zray,maxpts,npts,dir,ngrd, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do nn=1,npts
            !
            !        convert spatial cordinate to grid cordinate
            !
            xray1(nn)=1.+(xray(nn)-xmin)/delx
            yray1(nn)=1.+(yray(nn)-ymin)/dely
            zray1(nn)=1.+(zray(nn)+height-zmin)/delz
        enddo
        !
        !       make color choice and plot the grap
        !
        !        call sflush
        !        call gsplci(9)
        call curve3(xray1,yray1,zray1,npts)
        !
        !          2d plot
        !
        !        side=(-1.)**n
        do nn=1,npts
            !
            !        convert spatial cordinate to grid cordinate
            !
            ar=sqrt((yray(nn)-ydip)**2+(zray(nn)-zdip)**2)
            ar=sign(ar,yray(nn))
            !        ar=sign(ar,xcraft(2,n))
            !        ar=ar*side
            xray1(nn)=1.+(xray(nn)-xmin)/delx
            yray1(nn)=1.+(ar-ymin)/dely
            zray1(nn)=1.
        enddo
        !        call sflush
        !        call gsplci(1)
        call curve3(xray1,yray1,zray1,npts)
        !
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
