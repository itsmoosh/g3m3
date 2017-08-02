subroutine vct3d(vx,vy,vz,nx,ny,nz,m,xmin,xmax,ymin, &
    ymax,zmin,zmax,time,label,iv,strtch,rearth, &
    t,tx,ty,tz,work,mx,my,mz,muvwp2,mz2,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     plots 3-d vector field - iv sets velocity direction for
    !           color changes in arrows
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    real t(mx,my,mz),tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz), &
    work(muvwp2,muvwp2)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
     dimension vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    real eye(3)
    integer vpl,vpr,vpb,vpt
    character*4 wd1,wd2,wd3
    character*12 label
    character*15 title
    !
    !      dimension for plotted array
     my2=my/2+1
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
    !     set viewport size
    !
    !     full screen would be vpl=1 (left) vpr=32760 (right)
    !                         vpb=1 (bottom) vpt=32760 (top)
    !     drawing lines over the top using set33 with
    !                xa,xb,ya,yb set to the percentage fraction
    !                of these numbers
    !
    vpl=3200
    vpr=32760
    vpb=3200
    vpt=32760
    !
    !      set up evenly spaced gridding for t to be plotted
    !
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
    delx=amin1(delx,dely,delz)
    dely=amin1(delx,dely,delz)
    delz=amin1(delx,dely,delz)
    !
    axmax=xmin+delx*(mx-1)
    aymax=ymin+dely*(my-1)
    azmax=zmin+delz*(mz-1)
    !
    !      note delx,dely,delz used to convert real space to a position
    !        on grid with
    !         x_grid = 1 +(x_real-xmin)/delx
    !                or
    !         x_real = xmin + (x_grid -1.)*delx
    !
    !     back view -------------
    !
    !     initialize eye position
    !
    eye(1)=mx*5.
    eye(2)=my*4.
    eye(3)=mz*4.
    !
    !     plot earth
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
    call wtstr(.4,.975,label,2,0,0)
    title='back view'
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    write(title,'(1pe9.2)')vm
    title='magn '//title
    call wtstr(.5,.95,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
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
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')zmin
    title='y axis'
    call wtstr(.95,.87,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.92,.82,title,1,0,0)
    !
    call gselnt(0)
    call sflush
    !
    !      load t stuff to position earth and find max vector
    !
    vm=0.
    ddx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    ddy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    ddz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    do  k=1,mz
        az=zmin+(k-1.)*delz
        do  j=1,my
            ay=ymin+(j-1.)*dely
            do  i=1,mx
                ax=xmin+(i-1.)*delx
                !
                radius=sqrt((ax-xdip)**2+(ay-ydip)**2 &
                +(az-zdip)**2)
                t(i,j,k)=radius
                !
            enddo
        enddo
    enddo
    !
    !     plot earth and grid references
    !
    call gsplci(2)
    tisom=rearth
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !    draw axes lines
    !
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
    call sflush
    !
    !     determine vector field
    !
    !
    vm=0
    do k=1,mz
        az=zmin+(k-1.)*delz
        ak=1.+(az-grd_zmin(m))/ddz
        k1=ak
        k2=k1+1
        dz=ak-k1
        !
        do j=1,my
            !
            ay=ymin+(j-1.)*dely
            aj=1.+(ay-grd_ymin(m))/ddy
            j1=aj
            j2=j1+1
            dy=aj-j1
            !
            do i=1,mx
                ax=xmin+(i-1.)*delx
                ai=1.+(ax-grd_xmin(m))/ddx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                tx(i,j,k)=vx(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
                +vx(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
                +vx(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
                +vx(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
                +vx(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
                +vx(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
                +vx(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
                +vx(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                ty(i,j,k)=vy(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
                +vy(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
                +vy(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
                +vy(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
                +vy(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
                +vy(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
                +vy(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
                +vy(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                tz(i,j,k)=vz(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
                +vz(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
                +vz(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
                +vz(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
                +vz(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
                +vz(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
                +vz(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
                +vz(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                vm=amax1(vm,abs(tx(i,j,k)),abs(ty(i,j,k)),abs(tz(i,j,k)))
                radius=sqrt((ax-xdip)**2+(ay-ydip)**2 &
                +(az-zdip)**2)
                 if(radius.lt.rearth+1.)then
                    tx(i,j,k)=0.
                    ty(i,j,k)=0.
                    tz(i,j,k)=0.
                endif
                !
            enddo
        enddo
    enddo
    !
    vm=amax1(vm,0.00004)
    !
    !     draw points
    !
    ncol=3
    ncol1=13
    jskip=4
    iskip=2
    do k=1,mz-5,2
        z1=k*hk+0.05
        do j=jskip/2,my,jskip
            y1=j*hj+0.05
            jcol=5+2*(j/jskip)
            do i=1,mx,iskip
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
                    dz=iskip*hk*tz(i,j,k)/vm
                    !
                    x2=x1+stretch*dx
                    y2=y1+stretch*dy
                    z2=z1+stretch*dz
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
                    call sflush
                    if(iv.eq.1)then
                        if(dx.ge.0.0)call gsplci(ncol)
                        if(dx.lt.0.0)call gsplci(ncol1)
                    else if(iv.eq.2)then
                        if(dy.ge.0.0)call gsplci(ncol)
                        if(dy.lt.0.0)call gsplci(ncol1)
                    else
                        if(dz.ge.0.0)call gsplci(ncol)
                        if(dz.lt.0.0)call gsplci(ncol1)
                    endif
                    call arwxy(dx,dy,al,beta,adx1,ady1,adx2,ady2)
                    call arwxy(dy,dz,al,beta,ady1,adz1,ady2,adz2)
                    !
                    call line3(x2,y2,z2,x2+adx1,y2+ady1,z2+adz1)
                    call line3(x2,y2,z2,x2+adx2,y2+ady2,z2+adz2)
                endif
            endif
            !
        enddo
    enddo
    enddo
    !
    !     front view -------------
    !
    !     initialize eye position
    !
    eye(1)=-mx*5.
    eye(2)=my*4.
    eye(3)=mz*4.
    !
    !     plot earth
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    call wtstr(.4,.975,label,2,0,0)
    title='front view'
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    write(title,'(1pe8.1)')vm
    title='magnit '//title
    call wtstr(.5,.95,title,1,0,0)
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
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')azmax
    title='z axis'
    call wtstr(.95,.8,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.92,.75,title,1,0,0)
    !
    call gselnt(0)
    call sflush
    !
    !     plot earth and grid references
    !
    call gsplci(2)
    tisom=rearth
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !    draw axes lines
    !
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
    call sflush
    !
    !
    !     draw points
    !
    ncol=3
    ncol1=13
    jskip=4
    iskip=2
    do k=1,mz-5,2
        z1=k*hk+0.05
        do j=jskip/2,my,jskip
            y1=j*hj+0.05
            jcol=5+2*(j/jskip)
            do i=1,mx,iskip
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
                    dz=iskip*hk*tz(i,j,k)/vm
                    !
                    x2=x1+stretch*dx
                    y2=y1+stretch*dy
                    z2=z1+stretch*dz
                    call sflush
                    call gsplci(jcol)
                    !
                    call line3(x1,y1,z1,x2,y2,z2)
                    !
                    !        test if you can fit an arrow head to the line
                    !
                    if((abs(dx).ge.al2).or.(abs(dy).ge.al2) &
                        .or.(abs(dz).ge.al2)) then
                        !
                        !      draw arrows
                        !
                        !
                        call sflush
                        if(iv.eq.1)then
                            if(dx.ge.0.0)call gsplci(ncol)
                            if(dx.lt.0.0)call gsplci(ncol1)
                        else if(iv.eq.2)then
                            if(dy.ge.0.0)call gsplci(ncol)
                            if(dy.lt.0.0)call gsplci(ncol1)
                        else
                            if(dz.ge.0.0)call gsplci(ncol)
                            if(dz.lt.0.0)call gsplci(ncol1)
                        endif
                        call arwxy(dx,dy,al,beta,adx1,ady1,adx2,ady2)
                        call arwxy(dy,dz,al,beta,ady1,adz1,ady2,adz2)
                        !
                        call line3(x2,y2,z2,x2+adx1,y2+ady1,z2+adz1)
                        call line3(x2,y2,z2,x2+adx2,y2+ady2,z2+adz2)
                    endif
                endif
                !
            enddo
        enddo
    enddo
    !
    !     top view
    !
    !
    !     initialize eye position
    !
    eye(1)=mx/6
    eye(2)=my*4.
    eye(3)=mz*4.
    !
    !     plot earth
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    call wtstr(.4,.975,label,2,0,0)
    title='top-dawn view'
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    write(title,'(1pe8.1)')vm
    title='magnit '//title
    call wtstr(.5,.95,title,1,0,0)
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
    call gselnt(0)
    call sflush
    !
    !     plot earth and grid references
    !
    call gsplci(2)
    tisom=rearth
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !    draw axes lines
    !
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
    call sflush
    !
    !
    !     draw points
    !
    ncol=3
    ncol1=13
    jskip=4
    iskip=2
    do k=1,mz-5,2
        z1=k*hk+0.05
        do j=jskip/2,my,jskip
            y1=j*hj+0.05
            jcol=5+2*(j/jskip)
            do i=1,mx,iskip
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
                    dz=iskip*hk*tz(i,j,k)/vm
                    !
                    x2=x1+stretch*dx
                    y2=y1+stretch*dy
                    z2=z1+stretch*dz
                    call sflush
                    call gsplci(jcol)
                    !
                    call line3(x1,y1,z1,x2,y2,z2)
                    !
                    !        test if you can fit an arrow head to the line
                    !
                    if((abs(dx).ge.al2).or.(abs(dy).ge.al2) &
                    .or.(abs(dz).ge.al2)) then
                        !
                        !      draw arrows
                        !
                        call sflush
                        !
                        if(iv.eq.1)then
                            if(dx.ge.0.0)call gsplci(ncol)
                            if(dx.lt.0.0)call gsplci(ncol1)
                        else if(iv.eq.2)then
                            if(dy.ge.0.0)call gsplci(ncol)
                            if(dy.lt.0.0)call gsplci(ncol1)
                        else
                            if(dz.ge.0.0)call gsplci(ncol)
                            if(dz.lt.0.0)call gsplci(ncol1)
                        endif
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
    enddo
    !
    !     top view -dusk side
    !
    !
    !     initialize eye position
    !
    eye(1)=mx/6
    eye(2)=-my*4.
    eye(3)=mz*4.
    !
    !     plot earth
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    call wtstr(.4,.975,label,2,0,0)
    title='top-dusk view'
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    write(title,'(1pe8.1)')vm
    title='magnit '//title
    call wtstr(.5,.95,title,1,0,0)
    !
    write(wd1,'(f4.0)')xmin
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')zmin
    title='y axis'
    call wtstr(.22,.8,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.2,.75,title,1,0,0)
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')ymin
    write(wd3,'(f4.0)')zmin
    title='xaxis'
    call wtstr(.95,.3,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.92,.25,title,1,0,0)
    !
    !
    write(wd1,'(f4.0)')axmax
    write(wd2,'(f4.0)')aymax
    write(wd3,'(f4.0)')azmax
    title='back top'
    call wtstr(.9,.92,title,1,0,0)
    title=wd1//','//wd2//','//wd3
    call wtstr(.9,.88,title,1,0,0)
    !
    call gselnt(0)
    call sflush
    !
    !     plot earth and grid references
    !
    call gsplci(2)
    tisom=rearth
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
    !
    !    draw axes lines
    !
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
    call sflush
    !
    !
    !     draw points
    !
    ncol=3
    ncol1=13
    jskip=4
    iskip=2
    do k=1,mz-5,2
        z1=k*hk+0.05
        do j=jskip/2,my,jskip
            y1=j*hj+0.05
            jcol=5+2*(j/jskip)
            do i=1,mx,iskip
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
                    dz=iskip*hk*tz(i,j,k)/vm
                    !
                    x2=x1+stretch*dx
                    y2=y1+stretch*dy
                    z2=z1+stretch*dz
                    call sflush
                    call  gsplci(jcol)
                    !
                    call line3(x1,y1,z1,x2,y2,z2)
                    !
                    !        test if you can fit an arrow head to the line
                    !
                    if((abs(dx).ge.al2).or.(abs(dy).ge.al2) &
                    .or.(abs(dz).ge.al2)) then
                        !
                        !      draw arrows
                        !
                        call sflush
                        !
                        if(iv.eq.1)then
                            if(dx.ge.0.0)call gsplci(ncol)
                            if(dx.lt.0.0)call gsplci(ncol1)
                        else if(iv.eq.2)then
                            if(dy.ge.0.0)call gsplci(ncol)
                            if(dy.lt.0.0)call gsplci(ncol1)
                        else
                            if(dz.ge.0.0)call gsplci(ncol)
                            if(dz.lt.0.0)call gsplci(ncol1)
                        endif
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
    enddo
    return
end
