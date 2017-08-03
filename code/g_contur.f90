subroutine contur(stuff,nx,ny,nz,ngrd,m,xmin,xmax, &
    ymin,ymax,zmin,zmax,time,label,nlevs,ncon,rearth, &
    t,tt,t2,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
    !        nlevs number of isosurfaces to be plotted	- max 4
    !        ncon number of contours to be plotted		- max 14
    !
    dimension t(mx,my,mz),tt(mx,my,mz2),t2(mx,my,mz2), &
    work(muvwp2,muvwp2)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension stuff(nx,ny,nz,ngrd)
    !
    real eye(3),tlev(6),tcon(14)
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    !
    !      dimension for plotted array
    !
    my2=my/2+1
    !
    !      t is the stuff to be plotted covering the region
    !         xmin-xmax,ymin-ymax,zmin-zmax
    !
    !      set up evenly spaced gridding for t to be plotted
    !
    axmax=amin1(xmax,grd_xmax(m)-.00001)
    aymax=amin1(ymax,grd_ymax(m)-.00001)
    azmax=amin1(zmax,grd_zmax(m)-.00001)
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
    ddx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    ddy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    ddz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !      load t stuff
    !
    do k=1,mz
        az=zmin+delz*(k-1)
        ak=1.+(az-grd_zmin(m))/ddz
        k1=ak
        k2=k1+1
        dz=ak-k1
        do j=1,my
            ay=ymin+dely*(j-1)
            aj=1.+(ay-grd_ymin(m))/ddy
            j1=aj
            j2=j1+1
            dy=aj-j1
            !
            do i=1,mx
                ax=xmin+delx*(i-1)
                ai=1.+(ax-grd_xmin(m))/ddx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                t(i,j,k)=stuff(i1,j1,k1,m)*(1.-dx)*(1.-dy)*(1.-dz) &
                +stuff(i1,j1,k2,m)*(1.-dx)*(1.-dy)*(dz) &
                +stuff(i1,j2,k1,m)*(1.-dx)*(dy)*(1.-dz) &
                +stuff(i1,j2,k2,m)*(1.-dx)*(dy)*(dz) &
                +stuff(i2,j1,k1,m)*(dx)*(1.-dy)*(1.-dz) &
                +stuff(i2,j1,k2,m)*(dx)*(1.-dy)*(dz) &
                +stuff(i2,j2,k1,m)*(dx)*(dy)*(1.-dz) &
                +stuff(i2,j2,k2,m)*(dx)*(dy)*(dz)
                !
                radius=sqrt(ax**2+ay**2+az**2)
                tt(i,j,k)=radius
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
    alim=0.5
    dlev=alim*(tma-tmi)/(nlevs+1.)
    do n=1,nlevs
        tlev(n)=tmi+dlev*n
    enddo
    !
    !     set higher resolution contouring
    !
    dcon=alim*(tma-tmi)/(ncon+1.)
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
    !     ****************************************
    !          make multiple surface plots
    !     ****************************************
    !
    !     a view from behind
    !     ------------------
    !
    !     initialize eye position
    !
    eye(1)=mx*5.
    eye(2)=my*4.
    eye(3)=mz*4.
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
    call isoclrs
    !
    title=' : back view'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
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
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')azmax
    title='x axis'
    call wtstr(.2,.25,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.2,.2,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.6,.92,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.6,.89,title,1,0,0)
    !
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')zmin
    title='y axis'
    call wtstr(.95,.87,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.92,.82,title,1,0,0)
    !
    !     draw desired isosurfaces
    !
    do n=1,nlevs
        call gselnt(0)
        call sflush
        tisom=tlev(n)
        ncol=n*(dlev/dcon)
        call gsplci(ncol)
        call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,3, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !     plot earth and grid references
    !
    call gsplci(2)
    tisom=rearth
    call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !     put color bar
    !
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    call gsfais(1)
    call lbseti('cbl',1)
    call pcsetr('cs',1.25)
    call lblbar(0,0.1,0.9,0.,.1,ncon,1.,.3,lind,0,llbs,ncon,1)
    !
    !     a view from in front
    !     --------------------
    !
    !     initialize eye position
    !
    eye(1)=-mx*5.
    eye(2)=my*4.
    eye(3)=mz*4.
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : front view'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
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
    call wtstr(.22,.3,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.2,.25,title,1,0,0)
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')azmax
    title='bck pos'
    call wtstr(.6,.9,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.6,.85,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.95,.8,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.92,.75,title,1,0,0)
    !
    !     draw desired isosurfaces
    !
    do n=1,nlevs
        call gselnt(0)
        call sflush
        tisom=tlev(n)
        ncol=n*(dlev/dcon)
        call gsplci(ncol)
        call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,3, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !     plot earth and grid references
    !
    call gsplci(2)
    tisom=rearth
    call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !     put color bar
    !
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    call gsfais(1)
    call lbseti('cbl',1)
    call pcsetr('cs',1.25)
    call lblbar(0,0.1,0.9,0.,.1,ncon,1.,.3,lind,0,llbs,ncon,1)
    !
    !     a view from the side and top
    !     ----------------------------
    !
    !     initialize eye position
    !
    eye(1)=mx/2.
    eye(2)=my*4.
    eye(3)=mz*4.
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : side view'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
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
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')zmin
    title='x axis'
    call wtstr(.22,.8,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.2,.75,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')zmin
    title='yaxis'
    call wtstr(.95,.3,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.92,.25,title,1,0,0)
    !
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.9,.92,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.88,title,1,0,0)
    !
    !     draw desired isosurfaces -- i hope
    !
    do n=1,nlevs
        call gselnt(0)
        call sflush
        tisom=tlev(n)
        ncol=n*(dlev/dcon)
        call gsplci(ncol)
        call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,3, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !     plot earth and grid references
    !
    call gsplci(2)
    tisom=rearth
    call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !     put color bar
    !
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    call gsfais(1)
    call lbseti('cbl',1)
    call pcsetr('cs',1.25)
    call lblbar(0,0.1,0.9,0.,.1,ncon,1.,.3,lind,0,llbs,ncon,1)
    !
    !     ***************************************
    !       make contour plots down the middle
    !
    my2=my/2+1
    aymin=(ymin+aymax)/2.
    !
    !     zero subset array tt
    !
     do k=1,mz
        do j=1,my2
            do i=1,mx
                t2(i,j,k)=0.
            enddo
        enddo
    enddo
    !
    !     initialize eye position
    !
    eye(1)=mx/2.
    eye(2)=-my2*5.
    eye(3)=mz*3.5
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : dusk noon'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
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
    call wtstr(.93,.17,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.12,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')aymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.93,.9,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.85,title,1,0,0)
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')aymin
    write(wd3,'(f4.0)')zmin
    title='x axis'
    call wtstr(.18,.6,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.16,.55,title,1,0,0)
    !
    !          loading x-y plane - running from y = my/2 to  my
    !
    k=1
    do j=1,my2
        do i=1,mx
            t2(i,j,k)=t(i,j,k)
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
        call isosrf(t2,mx,mx,my2,my2,mz,eye,muvwp2,work,tisom,1, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !     re-zero tt
    !
    k=1
    do j=1,my2
        do i=1,mx
            t2(i,j,k)=0.
        enddo
    enddo
    !
    !         loading x-z plane
    !
    j=my2
    do k=1,mz
        do i=1,mx
            t2(i,j,k)=t(i,j,k)
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
        call isosrf(t2,mx,mx,my2,my2,mz,eye,muvwp2,work,tisom,2, &
        vpl,vpr,vpb,vpt)
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
    eye(2)=-10.
    eye(3)=mz*8.
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
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')aymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.93,.9,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.85,title,1,0,0)
    !
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')aymin
    write(wd3,'(f4.0)')zmin
    title='x axis'
    call wtstr(.05,.6,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.03,.55,title,1,0,0)
    !
    !     load data ------ right hand side
    !
    !          loading x-y plane - running from y = my/2 to  my
    !
    k=1
    do j=1,my
        do i=1,mx
            tt(i,j,k)=t(i,j,k)
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
    !     put color bar
    !
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    call gsfais(1)
    call lbseti('cbl',1)
    call pcsetr('cs',1.25)
    call lblbar(0,0.1,0.9,0.,.1,ncon,1.,.3,lind,0,llbs,ncon,1)
    return
end
