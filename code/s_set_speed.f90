!
!	This file contains two subroutines:
!	set_speed
!	set_speed_agrd
!
subroutine set_speed( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz,opx,opy,opz, &
    epres,qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
    bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
    rho,presx,presy,presz,presxy,presxz,presyz,px,py,pz, &
    rmassq,rmassh,rmasso,nx,ny,nz,ngrd, &
    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
    vlim,alf_lim,o_conc,fastest, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !    checks for minimum rho and negative pressure
    !     and resets value if necessary
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    dimension rho(nx,ny,nz),presx(nx,ny,nz),presy(nx,ny,nz), &
    presz(nx,ny,nz),px(nx,ny,nz),py(nx,ny,nz),pz(nx,ny,nz), &
    presxy(nx,ny,nz),presxz(nx,ny,nz),presyz(nx,ny,nz)
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !     determine speeds on the code and changes accordingly
    !
    do m=1,ngrd
        fastest=0.
        !
        pxmax=0.
        pymax=0.
        pzmax=0.
        pmax=0.
        csmax=0.
        alfmax=0.
        !
        vqlim=vlim
        vhlim=1.4*vlim
        volim=2.*vlim
        !
        !     vqlim=0.5*sqrt(1.+m)
        !     vhlim=0.5*sqrt(1.+m)
        !     volim=0.5*sqrt(1.+m)
        !
        cqlim=1.0*vqlim
        chlim=2.0*vhlim
        colim=4.0*volim
        cslim=vlim
        !
        ani=0.99995
        !
        !       find bulk velocities: species 1
        !
        !$omp  parallel do
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    rho(i,j,k)=qrho(i,j,k,m)
                    px(i,j,k)=qpx(i,j,k,m)
                    py(i,j,k)=qpy(i,j,k,m)
                    pz(i,j,k)=qpz(i,j,k,m)
                    !
                    !      limit anisotropy
                    !
                    apres=(qpresx(i,j,k,m)+qpresy(i,j,k,m)+qpresz(i,j,k,m))/3.0
                    presx(i,j,k)=ani*qpresx(i,j,k,m)+(1.-ani)*apres
                    presy(i,j,k)=ani*qpresy(i,j,k,m)+(1.-ani)*apres
                    presz(i,j,k)=ani*qpresz(i,j,k,m)+(1.-ani)*apres
                    presxy(i,j,k)=ani*qpresxy(i,j,k,m)
                    presxz(i,j,k)=ani*qpresxz(i,j,k,m)
                    presyz(i,j,k)=ani*qpresyz(i,j,k,m)
                enddo
            enddo
        enddo
        !
        d_min=0.01
        !
        do k=2,nz-1
            kp=k+1
            km=k-1
            do j=2,ny-1
                jp=j+1
                jm=j-1
                do i=2,nx-1
                    ip=i+1
                    im=i-1
                    !
                    arho=qrho(i,j,k,m)+1.e-5
                    aqvx=abs(qpx(i,j,k,m)/arho)
                    aqvy=abs(qpy(i,j,k,m)/arho)
                    aqvz=abs(qpz(i,j,k,m)/arho)
                    vq=sqrt(aqvx**2+aqvy**2+aqvz**2)
                    cq=sqrt(0.33*(qpresx(i,j,k,m)+qpresy(i,j,k,m) &
                    +qpresz(i,j,k,m))/arho)
                    !
                    if ((vq.gt.vqlim).or.(cq.gt.cqlim))then
                        qrho(i,j,k,m)=(rho(ip,j,k)+rho(im,j,k)+rho(i,jp,k) &
                        +rho(i,jm,k)+rho(i,j,kp)+rho(i,j,km))/6.
                        qpx(i,j,k,m)=(px(ip,j,k)+px(im,j,k)+px(i,jp,k) &
                        +px(i,jm,k)+px(i,j,kp)+px(i,j,km))/6.
                        qpy(i,j,k,m)=(py(ip,j,k)+py(im,j,k)+py(i,jp,k) &
                        +py(i,jm,k)+py(i,j,kp)+py(i,j,km))/6.
                        qpz(i,j,k,m)=(pz(ip,j,k)+pz(im,j,k)+pz(i,jp,k) &
                        +pz(i,jm,k)+pz(i,j,kp)+pz(i,j,km))/6.
                        qpresx(i,j,k,m)=(presx(ip,j,k)+presx(im,j,k)+presx(i,jp,k) &
                        +presx(i,jm,k)+presx(i,j,kp)+presx(i,j,km))/6.
                        qpresy(i,j,k,m)=(presy(ip,j,k)+presy(im,j,k)+presy(i,jp,k) &
                            +presy(i,jm,k)+presy(i,j,kp)+presy(i,j,km))/6.
                        qpresz(i,j,k,m)=(presz(ip,j,k)+presz(im,j,k)+presz(i,jp,k) &
                            +presz(i,jm,k)+presz(i,j,kp)+presz(i,j,km))/6.
                        qpresxy(i,j,k,m)=(presxy(ip,j,k)+presxy(im,j,k)+presxy(i,jp,k) &
                            +presxy(i,jm,k)+presxy(i,j,kp)+presxy(i,j,km))/6.
                        qpresxz(i,j,k,m)=(presxz(ip,j,k)+presxz(im,j,k)+presxz(i,jp,k) &
                            +presxz(i,jm,k)+presxz(i,j,kp)+presxz(i,j,km))/6.
                        qpresyz(i,j,k,m)=(presyz(ip,j,k)+presyz(im,j,k)+presyz(i,jp,k) &
                            +presyz(i,jm,k)+presyz(i,j,kp)+presyz(i,j,km))/6.
                    endif
                    !
                    pxmax=amax1(pxmax,aqvx)
                    pymax=amax1(pymax,aqvy)
                    pzmax=amax1(pzmax,aqvz)
                    csmax=amax1(csmax,cq)
                enddo
            enddo
        enddo
        !
        !       find bulk velocities: species 2
        !
        !$omp  parallel do
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    rho(i,j,k)=hrho(i,j,k,m)
                    px(i,j,k)=hpx(i,j,k,m)
                    py(i,j,k)=hpy(i,j,k,m)
                    pz(i,j,k)=hpz(i,j,k,m)
                    !
                    !      limit anisotropy
                    !
                    apres=(hpresx(i,j,k,m)+hpresy(i,j,k,m)+hpresz(i,j,k,m))/3.0
                    presx(i,j,k)=ani*hpresx(i,j,k,m)+(1.-ani)*apres
                    presy(i,j,k)=ani*hpresy(i,j,k,m)+(1.-ani)*apres
                    presz(i,j,k)=ani*hpresz(i,j,k,m)+(1.-ani)*apres
                    presxy(i,j,k)=ani*hpresxy(i,j,k,m)
                    presxz(i,j,k)=ani*hpresxz(i,j,k,m)
                    presyz(i,j,k)=ani*hpresyz(i,j,k,m)
                enddo
            enddo
        enddo
        !
        d_min=0.001
        !
        do k=2,nz-1
            kp=k+1
            km=k-1
            do j=2,ny-1
                jp=j+1
                jm=j-1
                do i=2,nx-1
                    ip=i+1
                    im=i-1
                     arho=hrho(i,j,k,m)+1.e-5
                    ahvx=abs(hpx(i,j,k,m)/arho)
                    ahvy=abs(hpy(i,j,k,m)/arho)
                    ahvz=abs(hpz(i,j,k,m)/arho)
                    vh=sqrt(ahvx**2+ahvy**2+ahvz**2)
                    ch=sqrt(0.33*(hpresx(i,j,k,m)+hpresy(i,j,k,m) &
                    +hpresz(i,j,k,m))/arho)
                    !
                    if ((vh.gt.vhlim).or.(ch.gt.chlim))then
                        hrho(i,j,k,m)=(rho(ip,j,k)+rho(im,j,k)+rho(i,jp,k) &
                        +rho(i,jm,k)+rho(i,j,kp)+rho(i,j,km))/6.
                        hpx(i,j,k,m)=(px(ip,j,k)+px(im,j,k)+px(i,jp,k) &
                        +px(i,jm,k)+px(i,j,kp)+px(i,j,km))/6.
                        hpy(i,j,k,m)=(py(ip,j,k)+py(im,j,k)+py(i,jp,k) &
                        +py(i,jm,k)+py(i,j,kp)+py(i,j,km))/6.
                        hpz(i,j,k,m)=(pz(ip,j,k)+pz(im,j,k)+pz(i,jp,k) &
                        +pz(i,jm,k)+pz(i,j,kp)+pz(i,j,km))/6.
                        hpresx(i,j,k,m)=(presx(ip,j,k)+presx(im,j,k)+presx(i,jp,k) &
                        +presx(i,jm,k)+presx(i,j,kp)+presx(i,j,km))/6.
                        hpresy(i,j,k,m)=(presy(ip,j,k)+presy(im,j,k)+presy(i,jp,k) &
                            +presy(i,jm,k)+presy(i,j,kp)+presy(i,j,km))/6.
                        hpresz(i,j,k,m)=(presz(ip,j,k)+presz(im,j,k)+presz(i,jp,k) &
                            +presz(i,jm,k)+presz(i,j,kp)+presz(i,j,km))/6.
                        hpresxy(i,j,k,m)=(presxy(ip,j,k)+presxy(im,j,k)+presxy(i,jp,k) &
                            +presxy(i,jm,k)+presxy(i,j,kp)+presxy(i,j,km))/6.
                        hpresxz(i,j,k,m)=(presxz(ip,j,k)+presxz(im,j,k)+presxz(i,jp,k) &
                            +presxz(i,jm,k)+presxz(i,j,kp)+presxz(i,j,km))/6.
                        hpresyz(i,j,k,m)=(presyz(ip,j,k)+presyz(im,j,k)+presyz(i,jp,k) &
                            +presyz(i,jm,k)+presyz(i,j,kp)+presyz(i,j,km))/6.
                    endif
                    !
                    pxmax=amax1(pxmax,ahvx)
                    pymax=amax1(pymax,ahvy)
                    pzmax=amax1(pzmax,ahvz)
                    csmax=amax1(csmax,ch)
                    !
                enddo
            enddo
        enddo
        !
        !       find bulk velocities: species 3
        !
        !$omp  parallel do
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    rho(i,j,k)=orho(i,j,k,m)
                    px(i,j,k)=opx(i,j,k,m)
                    py(i,j,k)=opy(i,j,k,m)
                    pz(i,j,k)=opz(i,j,k,m)
                    !
                    !      limit anisotropy
                    !
                    apres=(opresx(i,j,k,m)+opresy(i,j,k,m)+opresz(i,j,k,m))/3.0
                    presx(i,j,k)=ani*opresx(i,j,k,m)+(1.-ani)*apres
                    presy(i,j,k)=ani*opresy(i,j,k,m)+(1.-ani)*apres
                    presz(i,j,k)=ani*opresz(i,j,k,m)+(1.-ani)*apres
                    presxy(i,j,k)=ani*opresxy(i,j,k,m)
                    presxz(i,j,k)=ani*opresxz(i,j,k,m)
                    presyz(i,j,k)=ani*opresyz(i,j,k,m)
                enddo
            enddo
        enddo
        !
        d_min=0.001
        !
        do k=2,nz-1
            kp=k+1
            km=k-1
            do j=2,ny-1
                jp=j+1
                jm=j-1
                do i=2,nx-1
                    ip=i+1
                    im=i-1
                    arho=orho(i,j,k,m)+1.e-5
                    aovx=abs(opx(i,j,k,m)/arho)
                    aovy=abs(opy(i,j,k,m)/arho)
                    aovz=abs(opz(i,j,k,m)/arho)
                    vo=sqrt(aovx**2+aovy**2+aovz**2)
                    co=sqrt(0.33*(opresx(i,j,k,m)+opresy(i,j,k,m) &
                    +opresz(i,j,k,m))/arho)
                    !
                    if ((vo.gt.volim).or.(co.gt.colim))then
                        orho(i,j,k,m)=(rho(ip,j,k)+rho(im,j,k)+rho(i,jp,k) &
                        +rho(i,jm,k)+rho(i,j,kp)+rho(i,j,km))/6.
                        opx(i,j,k,m)=(px(ip,j,k)+px(im,j,k)+px(i,jp,k) &
                        +px(i,jm,k)+px(i,j,kp)+px(i,j,km))/6.
                        opy(i,j,k,m)=(py(ip,j,k)+py(im,j,k)+py(i,jp,k) &
                        +py(i,jm,k)+py(i,j,kp)+py(i,j,km))/6.
                        opz(i,j,k,m)=(pz(ip,j,k)+pz(im,j,k)+pz(i,jp,k) &
                        +pz(i,jm,k)+pz(i,j,kp)+pz(i,j,km))/6.
                        opresx(i,j,k,m)=(presx(ip,j,k)+presx(im,j,k)+presx(i,jp,k) &
                        +presx(i,jm,k)+presx(i,j,kp)+presx(i,j,km))/6.
                        opresy(i,j,k,m)=(presy(ip,j,k)+presy(im,j,k)+presy(i,jp,k) &
                            +presy(i,jm,k)+presy(i,j,kp)+presy(i,j,km))/6.
                        opresz(i,j,k,m)=(presz(ip,j,k)+presz(im,j,k)+presz(i,jp,k) &
                            +presz(i,jm,k)+presz(i,j,kp)+presz(i,j,km))/6.
                        opresxy(i,j,k,m)=(presxy(ip,j,k)+presxy(im,j,k)+presxy(i,jp,k) &
                            +presxy(i,jm,k)+presxy(i,j,kp)+presxy(i,j,km))/6.
                        opresxz(i,j,k,m)=(presxz(ip,j,k)+presxz(im,j,k)+presxz(i,jp,k) &
                            +presxz(i,jm,k)+presxz(i,j,kp)+presxz(i,j,km))/6.
                        opresyz(i,j,k,m)=(presyz(ip,j,k)+presyz(im,j,k)+presyz(i,jp,k) &
                            +presyz(i,jm,k)+presyz(i,j,kp)+presyz(i,j,km))/6.
                    endif
                    !
                    pxmax=amax1(pxmax,aovx)
                    pymax=amax1(pymax,aovy)
                    pzmax=amax1(pzmax,aovz)
                    csmax=amax1(csmax,co)
                    !
                enddo
            enddo
        enddo
        !
        !     do electron pressure and alvfen speed
        !
        call qvset(0.,bsx,nx*ny*nz)
        call qvset(0.,bsy,nx*ny*nz)
        call qvset(0.,bsz,nx*ny*nz)
        !
        call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
        call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
        call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
        !
        !      find magnitude of b
        !
        call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
        !
        !$omp  parallel do
        do k=1,nz
             do j=1,ny
                do i=1,nx
                    px(i,j,k)=bx(i,j,k,m)
                    py(i,j,k)=by(i,j,k,m)
                    pz(i,j,k)=bz(i,j,k,m)
                enddo
            enddo
        enddo
        !
        do k=2,nz-1
            kp=k+1
            km=k-1
            do j=2,ny-1
                jp=j+1
                jm=j-1
                do i=2,nx-1
                    ip=i+1
                    im=i-1
                    !
                    !       electron pressure
                    !
                    arho=qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m)+1.e-5
                    cs=sqrt(epres(i,j,k,m)/arho)
                    if(cs.gt.cslim)then
                        epres(i,j,k,m)=epres(i,j,k,m)*cslim/cs
                    endif
                    !
                    !       find alfven speed
                    !
                    abfld=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                    +bsz(i,j,k)**2)
                    alf=abfld/sqrt(arho)
                    !
                    alfmax=amax1(alfmax,alf)
                    !
                enddo
            enddo
        enddo
        !
        pmax=sqrt(pxmax**2+pymax**2+pzmax**2)
        !
        write(6,195)m,csmax,alfmax,pxmax,pymax,pzmax
        195 format(1x,i2,5(1x,1pe12.5))
        !
        fastest=amax1(fastest,sqrt(pxmax**2+pymax**2+pzmax**2 &
        +csmax**2+alfmax**2))
        !
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine set_speed_agrd( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz,opx,opy,opz, &
    epres,qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
    bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
    rho,presx,presy,presz,presxy,presxz,presyz,px,py,pz, &
    rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
    pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
    vlim,alf_lim,o_conc,fastest,isotropic)
    !
    !	Checks for minimum rho and negative pressure
    !		and resets value if necessary
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    dimension rho(nx,ny,nz),presx(nx,ny,nz),presy(nx,ny,nz), &
    presz(nx,ny,nz),px(nx,ny,nz),py(nx,ny,nz),pz(nx,ny,nz), &
    presxy(nx,ny,nz),presxz(nx,ny,nz),presyz(nx,ny,nz)
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz)
    !
    common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9), &
    grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9), &
    rx,ry,rz,xdip,ydip,zdip,rearth,b0, &
    sin_tilt,cos_tilt
    !
    !	Determine speeds on the code and changes accordingly
    !
    fastest=0.
    !
    pxmax=0.
    pymax=0.
    pzmax=0.
    pmax=0.
    csmax=0.
    alfmax=0.
    !
    vqlim=vlim
    vhlim=1.4*vlim
    volim=2.*vlim
    !
    !     vqlim=0.5*sqrt(1.+m)
    !     vhlim=0.5*sqrt(1.+m)
    !     volim=0.5*sqrt(1.+m)
    !
    cqlim=1.0*vqlim
    chlim=1.0*vhlim
    colim=4.0*volim
    cslim=vlim
    !
    ani=0.9995
    !
    !       find bulk velocities: species 1
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rho(i,j,k)=qrho(i,j,k,m)
                px(i,j,k)=qpx(i,j,k,m)
                py(i,j,k)=qpy(i,j,k,m)
                pz(i,j,k)=qpz(i,j,k,m)
                !
                !      limit anisotropy
                !
                apres=(qpresx(i,j,k,m)+qpresy(i,j,k,m)+qpresz(i,j,k,m))/3.0
                presx(i,j,k)=ani*qpresx(i,j,k,m)+(1.-ani)*apres
                presy(i,j,k)=ani*qpresy(i,j,k,m)+(1.-ani)*apres
                presz(i,j,k)=ani*qpresz(i,j,k,m)+(1.-ani)*apres
                presxy(i,j,k)=ani*qpresxy(i,j,k,m)
                presxz(i,j,k)=ani*qpresxz(i,j,k,m)
                presyz(i,j,k)=ani*qpresyz(i,j,k,m)
            enddo
        enddo
    enddo
    !
    d_min=0.0333
    !
    do k=2,nz-1
        kp=k+1
        km=k-1
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=2,nx-1
                ip=i+1
                im=i-1
                !
                arho=qrho(i,j,k,m)+1.e-5
                aqvx=abs(qpx(i,j,k,m)/arho)
                aqvy=abs(qpy(i,j,k,m)/arho)
                aqvz=abs(qpz(i,j,k,m)/arho)
                vq=sqrt(aqvx**2+aqvy**2+aqvz**2)
                cq=sqrt(0.33*(qpresx(i,j,k,m)+qpresy(i,j,k,m) &
                +qpresz(i,j,k,m))/arho)
                !
                if ((vq.gt.vqlim).or.(cq.gt.cqlim))then
                    qrho(i,j,k,m)=(rho(ip,j,k)+rho(im,j,k)+rho(i,jp,k) &
                        +rho(i,jm,k)+rho(i,j,kp)+rho(i,j,km))/6.
                    qpx(i,j,k,m)=(px(ip,j,k)+px(im,j,k)+px(i,jp,k) &
                        +px(i,jm,k)+px(i,j,kp)+px(i,j,km))/6.
                    qpy(i,j,k,m)=(py(ip,j,k)+py(im,j,k)+py(i,jp,k) &
                        +py(i,jm,k)+py(i,j,kp)+py(i,j,km))/6.
                    qpz(i,j,k,m)=(pz(ip,j,k)+pz(im,j,k)+pz(i,jp,k) &
                        +pz(i,jm,k)+pz(i,j,kp)+pz(i,j,km))/6.
                    qpresx(i,j,k,m)=(presx(ip,j,k)+presx(im,j,k)+presx(i,jp,k) &
                        +presx(i,jm,k)+presx(i,j,kp)+presx(i,j,km))/6.
                    qpresy(i,j,k,m)=(presy(ip,j,k)+presy(im,j,k)+presy(i,jp,k) &
                        +presy(i,jm,k)+presy(i,j,kp)+presy(i,j,km))/6.
                    qpresz(i,j,k,m)=(presz(ip,j,k)+presz(im,j,k)+presz(i,jp,k) &
                        +presz(i,jm,k)+presz(i,j,kp)+presz(i,j,km))/6.
                    qpresxy(i,j,k,m)=(presxy(ip,j,k)+presxy(im,j,k)+presxy(i,jp,k) &
                        +presxy(i,jm,k)+presxy(i,j,kp)+presxy(i,j,km))/6.
                    qpresxz(i,j,k,m)=(presxz(ip,j,k)+presxz(im,j,k)+presxz(i,jp,k) &
                        +presxz(i,jm,k)+presxz(i,j,kp)+presxz(i,j,km))/6.
                    qpresyz(i,j,k,m)=(presyz(ip,j,k)+presyz(im,j,k)+presyz(i,jp,k) &
                        +presyz(i,jm,k)+presyz(i,j,kp)+presyz(i,j,km))/6.
                endif
                !
                pxmax=amax1(pxmax,aqvx)
                pymax=amax1(pymax,aqvy)
                pzmax=amax1(pzmax,aqvz)
                csmax=amax1(csmax,cq)
            enddo
        enddo
    enddo
    !
    !       find bulk velocities : species 2
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rho(i,j,k)=hrho(i,j,k,m)
                px(i,j,k)=hpx(i,j,k,m)
                py(i,j,k)=hpy(i,j,k,m)
                pz(i,j,k)=hpz(i,j,k,m)
                !
                !      limit anisotropy
                !
                apres=(hpresx(i,j,k,m)+hpresy(i,j,k,m)+hpresz(i,j,k,m))/3.0
                presx(i,j,k)=ani*hpresx(i,j,k,m)+(1.-ani)*apres
                presy(i,j,k)=ani*hpresy(i,j,k,m)+(1.-ani)*apres
                presz(i,j,k)=ani*hpresz(i,j,k,m)+(1.-ani)*apres
                presxy(i,j,k)=ani*hpresxy(i,j,k,m)
                presxz(i,j,k)=ani*hpresxz(i,j,k,m)
                presyz(i,j,k)=ani*hpresyz(i,j,k,m)
             enddo
        enddo
    enddo
    !
    d_min=0.001
    !
    do k=2,nz-1
        kp=k+1
        km=k-1
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=2,nx-1
                ip=i+1
                im=i-1
                 arho=hrho(i,j,k,m)+1.e-5
                ahvx=abs(hpx(i,j,k,m)/arho)
                ahvy=abs(hpy(i,j,k,m)/arho)
                ahvz=abs(hpz(i,j,k,m)/arho)
                vh=sqrt(ahvx**2+ahvy**2+ahvz**2)
                ch=sqrt(0.33*(hpresx(i,j,k,m)+hpresy(i,j,k,m) &
                +hpresz(i,j,k,m))/arho)
                !
                if ((vh.gt.vhlim).or.(ch.gt.chlim)) then
                    hrho(i,j,k,m)=(rho(ip,j,k)+rho(im,j,k)+rho(i,jp,k) &
                    +rho(i,jm,k)+rho(i,j,kp)+rho(i,j,km))/6.
                    hpx(i,j,k,m)=(px(ip,j,k)+px(im,j,k)+px(i,jp,k) &
                    +px(i,jm,k)+px(i,j,kp)+px(i,j,km))/6.
                    hpy(i,j,k,m)=(py(ip,j,k)+py(im,j,k)+py(i,jp,k) &
                    +py(i,jm,k)+py(i,j,kp)+py(i,j,km))/6.
                    hpz(i,j,k,m)=(pz(ip,j,k)+pz(im,j,k)+pz(i,jp,k) &
                    +pz(i,jm,k)+pz(i,j,kp)+pz(i,j,km))/6.
                    hpresx(i,j,k,m)=(presx(ip,j,k)+presx(im,j,k)+presx(i,jp,k) &
                    +presx(i,jm,k)+presx(i,j,kp)+presx(i,j,km))/6.
                    hpresy(i,j,k,m)=(presy(ip,j,k)+presy(im,j,k)+presy(i,jp,k) &
                        +presy(i,jm,k)+presy(i,j,kp)+presy(i,j,km))/6.
                    hpresz(i,j,k,m)=(presz(ip,j,k)+presz(im,j,k)+presz(i,jp,k) &
                        +presz(i,jm,k)+presz(i,j,kp)+presz(i,j,km))/6.
                    hpresxy(i,j,k,m)=(presxy(ip,j,k)+presxy(im,j,k)+presxy(i,jp,k) &
                        +presxy(i,jm,k)+presxy(i,j,kp)+presxy(i,j,km))/6.
                    hpresxz(i,j,k,m)=(presxz(ip,j,k)+presxz(im,j,k)+presxz(i,jp,k) &
                        +presxz(i,jm,k)+presxz(i,j,kp)+presxz(i,j,km))/6.
                    hpresyz(i,j,k,m)=(presyz(ip,j,k)+presyz(im,j,k)+presyz(i,jp,k) &
                        +presyz(i,jm,k)+presyz(i,j,kp)+presyz(i,j,km))/6.
                endif
                !
                pxmax=amax1(pxmax,ahvx)
                pymax=amax1(pymax,ahvy)
                pzmax=amax1(pzmax,ahvz)
                csmax=amax1(csmax,ch)
                !
            enddo
        enddo
    enddo
    !
    !       find bulk velocities: species 3
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                rho(i,j,k)=orho(i,j,k,m)
                px(i,j,k)=opx(i,j,k,m)
                py(i,j,k)=opy(i,j,k,m)
                pz(i,j,k)=opz(i,j,k,m)
                !
                !      limit anisotropy
                !
                apres=(opresx(i,j,k,m)+opresy(i,j,k,m)+opresz(i,j,k,m))/3.0
                presx(i,j,k)=ani*opresx(i,j,k,m)+(1.-ani)*apres
                presy(i,j,k)=ani*opresy(i,j,k,m)+(1.-ani)*apres
                presz(i,j,k)=ani*opresz(i,j,k,m)+(1.-ani)*apres
                presxy(i,j,k)=ani*opresxy(i,j,k,m)
                presxz(i,j,k)=ani*opresxz(i,j,k,m)
                presyz(i,j,k)=ani*opresyz(i,j,k,m)
            enddo
        enddo
    enddo
    !
    d_min=0.001
    !
    do k=2,nz-1
        kp=k+1
        km=k-1
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=2,nx-1
                ip=i+1
                im=i-1
                arho=orho(i,j,k,m)+1.e-5
                aovx=abs(opx(i,j,k,m)/arho)
                aovy=abs(opy(i,j,k,m)/arho)
                aovz=abs(opz(i,j,k,m)/arho)
                vo=sqrt(aovx**2+aovy**2+aovz**2)
                co=sqrt(0.33*(opresx(i,j,k,m)+opresy(i,j,k,m) &
                +opresz(i,j,k,m))/arho)
                !
                if ((vo.gt.volim).or.(co.gt.colim))then
                    orho(i,j,k,m)=(rho(ip,j,k)+rho(im,j,k)+rho(i,jp,k) &
                    +rho(i,jm,k)+rho(i,j,kp)+rho(i,j,km))/6.
                    opx(i,j,k,m)=(px(ip,j,k)+px(im,j,k)+px(i,jp,k) &
                    +px(i,jm,k)+px(i,j,kp)+px(i,j,km))/6.
                    opy(i,j,k,m)=(py(ip,j,k)+py(im,j,k)+py(i,jp,k) &
                    +py(i,jm,k)+py(i,j,kp)+py(i,j,km))/6.
                    opz(i,j,k,m)=(pz(ip,j,k)+pz(im,j,k)+pz(i,jp,k) &
                    +pz(i,jm,k)+pz(i,j,kp)+pz(i,j,km))/6.
                    opresx(i,j,k,m)=(presx(ip,j,k)+presx(im,j,k)+presx(i,jp,k) &
                    +presx(i,jm,k)+presx(i,j,kp)+presx(i,j,km))/6.
                    opresy(i,j,k,m)=(presy(ip,j,k)+presy(im,j,k)+presy(i,jp,k) &
                        +presy(i,jm,k)+presy(i,j,kp)+presy(i,j,km))/6.
                    opresz(i,j,k,m)=(presz(ip,j,k)+presz(im,j,k)+presz(i,jp,k) &
                        +presz(i,jm,k)+presz(i,j,kp)+presz(i,j,km))/6.
                    opresxy(i,j,k,m)=(presxy(ip,j,k)+presxy(im,j,k)+presxy(i,jp,k) &
                        +presxy(i,jm,k)+presxy(i,j,kp)+presxy(i,j,km))/6.
                    opresxz(i,j,k,m)=(presxz(ip,j,k)+presxz(im,j,k)+presxz(i,jp,k) &
                        +presxz(i,jm,k)+presxz(i,j,kp)+presxz(i,j,km))/6.
                    opresyz(i,j,k,m)=(presyz(ip,j,k)+presyz(im,j,k)+presyz(i,jp,k) &
                        +presyz(i,jm,k)+presyz(i,j,kp)+presyz(i,j,km))/6.
                endif
                !
                pxmax=amax1(pxmax,aovx)
                pymax=amax1(pymax,aovy)
                pzmax=amax1(pzmax,aovz)
                csmax=amax1(csmax,co)
                !
            enddo
        enddo
    enddo
    !
    !	Do electron pressure and alfven speed
    !
    call qvset(0.,bsx,nx*ny*nz)
    call qvset(0.,bsy,nx*ny*nz)
    call qvset(0.,bsz,nx*ny*nz)
    !
    call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
    call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
    call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
    !
    !	Find magnitude of b
    !
    call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
         do j=1,ny
            do i=1,nx
                px(i,j,k)=bx(i,j,k,m)
                py(i,j,k)=by(i,j,k,m)
                pz(i,j,k)=bz(i,j,k,m)
            enddo
        enddo
    enddo
    !
    do k=2,nz-1
        kp=k+1
        km=k-1
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=2,nx-1
                ip=i+1
                im=i-1
                !
                !	Electron pressure
				!
                arho=qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m)+1.e-5
                cs=sqrt(epres(i,j,k,m)/arho)
                if(cs.gt.cslim)then
                    epres(i,j,k,m)=epres(i,j,k,m)*cslim/cs
                endif
                !
                !	Find alfven speed
                !
                abfld=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                +bsz(i,j,k)**2)
                alf=abfld/sqrt(arho)
                 alfmax=amax1(alfmax,alf)
                !
            enddo
        enddo
    enddo
    !
    pmax=sqrt(pxmax**2+pymax**2+pzmax**2)
    !
    write(6,195)m,csmax,alfmax,pxmax,pymax,pzmax
    195 format(1x,i2,5(1x,1pe12.5))
    !
    fastest=amax1(fastest,sqrt(pxmax**2+pymax**2+pzmax**2 &
    +csmax**2+alfmax**2))
    !
    return
end
