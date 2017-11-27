!
!	This file contains two subroutines:
!	concraft
!	concraft_fix
!
subroutine concraft(stuff,cross,nx,ny,nz,box, &
    xcraft,ncraft,re_equiv,time,label,start,frac,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot cross-plane quantities at the
    !       spacecraft position
    !       percent - percentage level of max value for isosurface
    !        xcraft is the position of the spacecraft : assumed to
    !                  be in simulation units
    !        ncraft is the number of spacecraft to plot
    !
    !     dimension of work array show be ny-2,nz-2
    !          (do not plot edge arrays)
	!
    integer box
    logical start
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    real xrays(2),yrays(2)
    dimension xcraft(4,ncraft),stuff(nx,ny,nz),cross(ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    character*4 wd1,wd2,wd3
    character*12 label
    character*20 title
    !
    !      dimension for plotted array
    !
    my=ny
    mz=nz
    !
    !     set color table
    !
    call sflush
    call isoclrs
    !
    !      make a 2-d cross tail cut for either imp 8 and geotail
    !
    ncon=13
    xmin=grd_xmin(box)
    xmax=grd_xmax(box)
    ymin=grd_ymin(box)
    ymax=grd_ymax(box)
    zmin=grd_zmin(box)
    zmax=grd_zmax(box)
    delx=(xmax-xmin)/(nx-1.)
    dely=(ymax-ymin)/(ny-1.)
    delz=(zmax-zmin)/(nz-1.)
    !
    do n=ncraft,2,-1
        ax=xcraft(1,n)/re_equiv
        ay=xcraft(2,n)/re_equiv
        az=xcraft(3,n)/re_equiv
        ai=1.+(ax-xmin)/delx
        aj=1.+(ay-ymin)/dely
        ak=1.+(az-zmin)/delz
        if( (ax.ge.xmin).and.(ax.le.xmax).and. &
        (ay.ge.ymin).and.(ay.le.ymax).and. &
        (az.ge.zmin).and.(az.le.zmax) ) then
        !
        i1=ai
        i2=ai+1
        dx=ai-i1
        ddx=1.-dx
        fmin=0.
        fmax=0.
        do k=1,nz
            do j=1,ny
                cross(j,k)=stuff(i1,j,k)*ddx+ &
                stuff(i2,j,k)*dx
                fmin=amin1(fmin,cross(j,k))
                fmax=amax1(fmax,cross(j,k))
            enddo
        enddo
        !
        call frame
        call gselnt(0)
        call gsplci(1)
        call agseti('frame.',2)
        call agseti('set.',-1)
        !
        write(title,'(f7.3)')time
        title='ut = '//title
        call wtstr(.7,.975,title,2,0,0)
        call wtstr(.4,.975,label,2,0,0)
        write(title,'(f7.4)')fmin
        title='fmin'//title
        call wtstr(.15,.98,title,1,0,0)
        write(title,'(f7.4)')fmax
        call wtstr(.3,.98,title,1,0,0)
        title='fmax'//title
        if(fmin*fmax.lt.0.)then
            cmax=amax1(fmax,abs(fmin))
            cmin=-cmax
        else
            cmax=fmax/frac
            cmin=fmin/frac
        endif
        clev=(cmax-cmin)/(ncon+1.)
        if(.not.start) call cnrccf(cross,my,my,mz,cmin, &
            cmax,clev,0,-1,-1012,2.5,0,1)
            irecmn=18
            !     call conrec(cross,my,my,mz,cmin,
            !    +             cmax,clev,0,-1,-1012)
            !
            !     draw in position of spacecraft
            !
            call gsplci(16)
            nl=11
            n2=(nl-1)/2
            do nn=1,nl
                ns=(nn-n2)/2
                shift=ns*.025
                xrays(1)=aj-.5
                xrays(2)=aj+.5
                yrays(1)=ak+shift
                yrays(2)=ak+shift
                call curve(xrays,yrays,2)
                xrays(1)=aj+shift
                xrays(2)=aj+shift
                yrays(1)=ak+.5
                yrays(2)=ak-.5
                call curve(xrays,yrays,2)
            enddo
        endif
    enddo
        !
    return
end
!
!
!	******************************************
!
!
subroutine concraft_fix(stuff,cross,nx,ny,nz,box, &
    xcraft,ncraft, re_equiv,time,label,start,alo,ahi,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot cross-plane quantities at the
    !       spacecraft position
    !       percent - percentage level of max value for isosurface
    !        xcraft is the position of the spacecraft : assumed to
    !                  be in simulation units
    !        ncraft is the number of spacecraft to plot
    !
    !     dimension of work array show be ny-2,nz-2
    !          (do not plot edge arrays)
	!
    integer box
    logical start
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    real xrays(2),yrays(2)
    dimension xcraft(4,ncraft),stuff(nx,ny,nz),cross(ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    character*4 wd1,wd2,wd3
    character*12 label
    character*20 title
    !
    !      dimension for plotted array
    !
    my=ny
    mz=nz
    !
    !     set color table
    !
    call sflush
    call isoclrs
    !
    !      make a 2-d cross tail cut for either imp 8 and geotail
    !
    ncon=13
    xmin=grd_xmin(box)
    xmax=grd_xmax(box)
    ymin=grd_ymin(box)
    ymax=grd_ymax(box)
    zmin=grd_zmin(box)
    zmax=grd_zmax(box)
    delx=(xmax-xmin)/(nx-1.)
    dely=(ymax-ymin)/(ny-1.)
    delz=(zmax-zmin)/(nz-1.)
    !
    do n=ncraft,2,-1
        ax=xcraft(1,n)/re_equiv
        ay=xcraft(2,n)/re_equiv
        az=xcraft(3,n)/re_equiv
        ai=1.+(ax-xmin)/delx
        aj=1.+(ay-ymin)/dely
        ak=1.+(az-zmin)/delz
        if( (ax.ge.xmin).and.(ax.le.xmax).and. &
        (ay.ge.ymin).and.(ay.le.ymax).and. &
        (az.ge.zmin).and.(az.le.zmax) ) then
        !
        i1=ai
        i2=ai+1
        dx=ai-i1
        ddx=1.-dx
        fmin=0.
        fmax=0.
        do k=1,nz
            do j=1,ny
                cross(j,k)=stuff(i1,j,k)*ddx+ &
                stuff(i2,j,k)*dx
                fmin=amin1(fmin,cross(j,k))
                fmax=amax1(fmax,cross(j,k))
            enddo
        enddo
        !
        call frame
        call gselnt(0)
        call gsplci(1)
        call agseti('frame.',2)
        call agseti('set.',-1)
        !
        write(title,'(f7.3)')time
        title='ut = '//title
        call wtstr(.7,.975,title,2,0,0)
        call wtstr(.4,.975,label,2,0,0)
        write(title,'(f7.4)')fmin
        title='fmin'//title
        call wtstr(.15,.98,title,1,0,0)
        write(title,'(f7.4)')fmax
        call wtstr(.3,.98,title,1,0,0)
        title='fmax'//title
        !
        cmin=alo
        cmax=ahi
        clev=(cmax-cmin)/(ncon+1.)
        if(.not.start) call cnrccf(cross,my,my,mz,cmin, &
        cmax,clev,0,-1,-1012,2.5,0,1)
        irecmn=18
        !     call conrec(cross,my,my,mz,cmin,
        !    +             cmax,clev,0,-1,-1012)
        !
        !     draw in position of spacecraft
        !
        call gsplci(16)
        nl=11
        n2=(nl-1)/2
        do nn=1,nl
            ns=(nn-n2)/2
            shift=ns*.025
            xrays(1)=aj-.5
            xrays(2)=aj+.5
            yrays(1)=ak+shift
            yrays(2)=ak+shift
            call curve(xrays,yrays,2)
            xrays(1)=aj+shift
            xrays(2)=aj+shift
            yrays(1)=ak+.5
            yrays(2)=ak-.5
            call curve(xrays,yrays,2)
        enddo
    endif
    enddo
    !
    return
end
