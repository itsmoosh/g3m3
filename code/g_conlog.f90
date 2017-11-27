subroutine conlog(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
    xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
    time,label,nlevs,ncon,ivel,strtch,tlow,tlim, &
    tx,ty,tz,t,t2,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
    !      assumes log scales
    !        nlevs number of isossurfaces to be plotted - max 4
    !        ncon  number of contours to be plotted     - max 14
    !
    real stuff(nx,ny,nz,n_grids),vx(nx,ny,nz), &
    vy(nx,ny,nz),vz(nx,ny,nz)
    dimension t(mx,my,mz),tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz), &
    t2(mx,my,mz),work(muvwp2,muvwp2)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    real eye(3),tlev(6),tcon(18)
    character*4 llbs(18),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(18)
    integer box
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
    !
    !     skip parameters for field lines and arrows
    !
    jskip=(ny-1)/20
    iskip=(nx-1)/20+1
    !
    !      dimension for plotted array
    !
    my2=my/2+1
    !
    t=0.
    t2=0.
    tx=0.
    ty=0.
    tz=0.
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
    vm=0.0
    !
    !      angle of arrow
    !
    beta=30./180.*3.1415926
    !
    !     set color table
    !
    call sflush
    call isoclrs_hot
    !
    !      t is the stuff to be plotted covering the region
    !         xmin-xmax,ymin-ymax,zmin-zmax
    !
    !      set up evenly spaced gridding for t to be plotted
    !
    axmax=amin1(xmax,grd_xmax(box)-.0001)
    aymax=amin1(ymax,grd_ymax(box)-.0001)
    azmax=amin1(zmax,grd_zmax(box)-.0001)
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
    !      load t stuff
    !
    call filldatav(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
    t,t2,tx,ty,tz,mx,my,mz,xmin,ymin,zmin, &
    delx,dely,delz,vm, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
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
    dcon=(tlim-tlow)/(ncon-1.)
    !      ampl=tlim
    !      write(magnif,'(1pe9.2)') ampl
    do n=1,ncon
        tcon(n)=tlow+(n-1)*dcon
        lind(n)=n+1
        write(llbs(n),'(f4.1)') tcon(n)
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
    !     ***************************************
    !       make contour plots down the middle
    !
    my2=my/2+1
    aymin=(ymin+aymax)/2.
    !
    !     initialize eye position
    !
    !     eye(1)=-mx
    !     eye(2)=-my2*2.
    !     eye(3)=mz*3.0
    !
    eye(1)=mx/2
    eye(2)=-my2*2.
    eye(3)=mz2*3.5
    !
    !     initialize viewport and frame headings
    !
    !     zero subset array tt
    !
    do k=1,mz
        do j=1,my
            do i=1,mx
                t2(i,j,k)=0.
            enddo
        enddo
    enddo
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    call wtstr(.4,.975,label,2,0,0)
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
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')aymin
    write(wd3,'(f4.0)')zmin
    title='x axis'
    call wtstr(.18,.9,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.16,.85,title,1,0,0)
    !
    !          loading x-y plane - running from y = my/2 to  my
    !
    k=1
    kk=mz2
    do j=1,my
        do i=1,mx
            t2(i,j,k)=t(i,j,kk)
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
        call isosrf(t2,mx,mx,my,my,mz,eye,muvwp2,work,tisom,1, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !     re-zero tt
    !
    k=1
    do j=1,my
        do i=1,mx
            t2(i,j,k)=0.
        enddo
    enddo
    !
    !         loading x-z plane
    !
    j=my
    jj=my2
    do k=1,mz
        do i=1,mx
            t2(i,j,k)=t(i,jj,k)
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
        call isosrf(t2,mx,mx,my,my,mz,eye,muvwp2,work,tisom,2, &
        vpl,vpr,vpb,vpt)
    enddo
    !
    !     loading y-zplane
    !
    do k=1,mz
        do i=1,mx
            t2(i,j,k)=0.
        enddo
    enddo
    icut=1.+(xcut-xmin)/delx
    i=mx
    do k=1,mz
        do j=1,my
            t2(i,j,k)=t(icut,j,k)
        enddo
    enddo
    do n=1,ncon
        call gselnt(0)
        call sflush
        tisom=tcon(n)
        call gsplci(n+1)
        call isosrf(t2,mx,mx,my,my,mz,eye,muvwp2,work,tisom,4, &
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
    call set3(0.0,1.,0.0,1.,x1,x2,y1,y2,z1,z2,eye)
    call sflush
    call gsplci(1)
    call line3(x1,y1,z1,x2,y1,z1)
    !     call line3(x1,y1,z2,x2,y1,z2)
    call line3(x1,y1,z1,x1,y2,z1)
    !     call line3(x1,y1,z2,x1,y2,z2)
    !     call line3(x1,y1,z1,x1,y1,z2)
    call sflush
    !
    !     draw points - noon midnight sector
    !
    ncol=3
    j=my2
    jj=my
    !     iskip=2
    jcol=1
    do k=2,mz-1,2
        z1=k*hk+0.05
        y1=jj*hj-0.05
         do i=2,mx-1,iskip
            x1=i*hi+0.05
            !
            !        find relative size of arrow and skip out if too small
            !
            vmag=sqrt(tx(i,j,k)**2+ty(i,j,k)**2+tz(i,j,k)**2)
            if(vmag*strtch.ge.0.05*vm) then
                !
                !        arrow length
                !
                dx=iskip*hi*tx(i,j,k)/vm
                !        dy=iskip*hj*ty(i,j,k)/vm
                dy=0.
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
                call gsplci(jcol)
                call line3(x1,y1,z1,x2,y2,z2)
                !
                !        test if you can fit an arrow head to the line
                !
                if((abs(dx).ge.al2).or.(abs(dy).ge.al2) &
                    .or.(abs(dz).ge.al2)) then
                    !
                    !      draw arrows - set color for arror heads
                    !
                    if((ivel.eq.1).and.(dx.ge.0.0)) call gsplci(ncon)
                    if((ivel.eq.1).and.(dx.lt.0.0)) call gsplci(5)
                    if((ivel.eq.2).and.(dy.ge.0.0)) call gsplci(ncon)
                    if((ivel.eq.2).and.(dy.lt.0.0)) call gsplci(5)
                    if((ivel.eq.3).and.(dz.ge.0.0)) call gsplci(ncon)
                    if((ivel.eq.3).and.(dz.lt.0.0)) call gsplci(5)
                    !
                    call arwxy(dx,dy,al,beta,adx1,ady1,adx2,ady2)
                    call arwxy(dy,dz,al,beta,ady1,adz1,ady2,adz2)
                    !
                    call line3(x2,y2,z2,x2+adx1,y2+ady1,z2+adz1)
                    call line3(x2,y2,z2,x2+adx2,y2+ady2,z2+adz2)
                    !
                endif
            endif
        enddo
    enddo
    !
    !     draw points : equatorial
    !
    ncol=3
    k=mz2
    kk=1
    !     jskip=2
    !     iskip=2
    z1=kk*hk+0.05
    do j=2,my-1,jskip
        y1=j*hj+0.05
        do i=2,mx-1,iskip
            x1=i*hi+0.05
            !
            !        find relative size of arrow and skip out if too small
            !
            vmag=sqrt(tx(i,j,k)**2+ty(i,j,k)**2+tz(i,j,k)**2)
            if(vmag*strtch.ge.0.05*vm) then
                !
                !        arrow length
                !
                dx=iskip*hi*tx(i,j,k)/vm
                dy=iskip*hj*ty(i,j,k)/vm
                !        dz=iskip*hk*tz(i,j,k)/vm
                dz=0.
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
                call line3(x1,y1,z1,x2,y2,z2)
                !
                !        test if you can fit an arrow head to the line
                !
                if((abs(dx).ge.al2).or.(abs(dy).ge.al2) &
                    .or.(abs(dz).ge.al2)) then
                    !
                    !      draw arrows - set color for arror heads
                    !
                    call sflush
                    if(dx.ge.0.0) call gsplci(ncon)
                    if(dx.lt.0.0) call gsplci(5)
                    !
                    call arwxy(dx,dy,al,beta,adx1,ady1,adx2,ady2)
                    call arwxy(dy,dz,al,beta,ady1,adz1,ady2,adz2)
                    !
                    call line3(x2,y2,z2,x2+adx1,y2+ady1,z2+adz1)
                    call line3(x2,y2,z2,x2+adx2,y2+ady2,z2+adz2)
                    !
                endif
            endif
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
