    !
    !      **************************
    !
    subroutine push_elec(epres,oldepres,wrkepres,evx,evy,evz, &
    gamma,gamma1,nx,ny,nz,ngrd,m,delt,rx,ry,rz)
    !
    !      evolves the electron pressure equation
    !
    !
    dimension epres(nx,ny,nz,ngrd),oldepres(nx,ny,nz,ngrd), &
    wrkepres(nx,ny,nz,ngrd)
    !
    dimension evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
    !
    !
    !       set physical grid spacing
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    !
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kk=k+1
    
        do j=2,ny-1
            jm=j-1
            jj=j+1
    
            do i=2,nx-1
                ii=i+1
                im=i-1
                egradp_x=(wrkepres(ii,j,k,m)-wrkepres(im,j,k,m))/dxt
                egradp_y=(wrkepres(i,jj,k,m)-wrkepres(i,jm,k,m))/dyt
                egradp_z=(wrkepres(i,j,kk,m)-wrkepres(i,j,km,m))/dzt
                !
                !
                !       pressure equations:
                !
                epres(i,j,k,m)=oldepres(i,j,k,m)-delt*gamma* &
                ( ( (wrkepres(ii,j,k,m)*evx(ii,j,k) &
                -wrkepres(im,j,k,m)*evx(im,j,k))/dxt ) + &
                ( (wrkepres(i,jj,k,m)*evy(i,jj,k) &
                -wrkepres(i,jm,k,m)*evy(i,jm,k))/dyt ) + &
                ( (wrkepres(i,j,kk,m)*evz(i,j,kk) &
                -wrkepres(i,j,km,m)*evz(i,j,km))/dzt ) )
                epres(i,j,k,m)=epres(i,j,k,m) &
                + delt*gamma1*( &
                evx(i,j,k)*egradp_x+evy(i,j,k)*egradp_y &
                +evz(i,j,k)*egradp_z )
                !
            enddo
        enddo
    enddo
    !
    return
end

!
!      **************************
!
subroutine push_ion(qrho,qpresx,qpresy,qpresz, &
    qpx,qpy,qpz, &
    oldqrho,oldqpresx,oldqpresy,oldqpresz, &
    oldqpx,oldqpy,oldqpz, &
    wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
    wrkqpx,wrkqpy,wrkqpz, &
    qpresxy,qpresxz,qpresyz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    wrkqpresxy,wrkqpresxz,wrkqpresyz, &
    bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq, &
    qvx,qvy,qvz,tvx,tvy,tvz,gamma,gamma1, &
    nx,ny,nz,ngrd,m,delt,grav,re_equiv,reynolds, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,ani,isotropic)
    !
    !      standard runge-kutta push for ion equations
    !
    dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    qvx(nx,ny,nz),qvy(nx,ny,nz),qvz(nx,ny,nz), &
    tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
    btot(nx,ny,nz)
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    oldqrho(nx,ny,nz,ngrd),oldqpresx(nx,ny,nz,ngrd), &
    oldqpresy(nx,ny,nz,ngrd),oldqpresz(nx,ny,nz,ngrd), &
    oldqpresxy(nx,ny,nz,ngrd),oldqpresxz(nx,ny,nz,ngrd), &
    oldqpresyz(nx,ny,nz,ngrd), &
    oldqpx(nx,ny,nz,ngrd),oldqpy(nx,ny,nz,ngrd), &
    oldqpz(nx,ny,nz,ngrd), &
    wrkqrho(nx,ny,nz,ngrd),wrkqpresx(nx,ny,nz,ngrd), &
    wrkqpresy(nx,ny,nz,ngrd),wrkqpresz(nx,ny,nz,ngrd), &
    wrkqpresxy(nx,ny,nz,ngrd),wrkqpresxz(nx,ny,nz,ngrd), &
    wrkqpresyz(nx,ny,nz,ngrd), &
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    
    !
    !       maximum ion cyclotron frequency in normalized units
    !           that that can be rsolved
    bmax=0.20*rmassq/delt/reynolds
    !
    !        set distance scales
    !
    ddx=(grd_xmax(m)-grd_xmin(m))/(nx-1)
    ddy=(grd_ymax(m)-grd_ymin(m))/(ny-1)
    ddz=(grd_zmax(m)-grd_zmin(m))/(nz-1)
    !
    dxt=2.*ddx
    dyt=2.*ddy
    dzt=2.*ddz
    !
    d_min=1.e-6
    !
    !       maximum ion cyclotron frequency in normalized units
    !           that that can be rsolved
    rmass=rmassq
    bmax=0.20*rmassq/delt/reynolds
    !
    ! parallelizes loop rw, aug. 17, 2004
    !$omp  parallel do
    do k=1,nz
    
        do j=1,ny
    
            do i=1,nx
                aqrho=amax1(wrkqrho(i,j,k,m),d_min)
                !
                qvx(i,j,k)=wrkqpx(i,j,k,m)/aqrho
                qvy(i,j,k)=wrkqpy(i,j,k,m)/aqrho
                qvz(i,j,k)=wrkqpz(i,j,k,m)/aqrho
                !
            enddo          ! k loop
        enddo          ! j loop
    enddo          ! i loop
    !
    !     begin main loop for finding all terms at n+1/2 from i=2,nx-1
    !           and from j=2,ny-1 using boundary conditions
    
    ! parallelizes the loop bu, oct. 10 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kk=k+1
        az=grd_zmin(m)+ddz*(k-1)
    
        !        find estimates for the fluid at n+1/2
    
        do j=2,ny-1
            jm=j-1
            jj=j+1
            ay=grd_ymin(m)+ddy*(j-1)
            !
            do i=2,nx-1
                ax=grd_xmin(m)+ddx*(i-1)
                ii=i+1
                im=i-1
                !
                !       species 1
                !
                qrho(i,j,k,m)=oldqrho(i,j,k,m)-delt*( &
                ( (wrkqpx(ii,j,k,m)- wrkqpx(im,j,k,m))/dxt) &
                +( (wrkqpy(i,jj,k,m)- wrkqpy(i,jm,k,m))/dyt) &
                +( (wrkqpz(i,j,kk,m)- wrkqpz(i,j,km,m))/dzt) )
                !
                aqpx=oldqpx(i,j,k,m)-delt*( &
                ((wrkqpx(ii,j,k,m)*qvx(ii,j,k)- wrkqpx(im,j,k,m)*qvx(im,j,k)) &
                /dxt ) &
                +((wrkqpx(i,jj,k,m)*qvy(i,jj,k)- wrkqpx(i,jm,k,m)*qvy(i,jm,k)) &
                /dyt ) &
                +((wrkqpx(i,j,kk,m)*qvz(i,j,kk)- wrkqpx(i,j,km,m)*qvz(i,j,km)) &
                /dzt ) )
                !
                aqpy=oldqpy(i,j,k,m)-delt*( &
                ((wrkqpy(ii,j,k,m)*qvx(ii,j,k)- wrkqpy(im,j,k,m)*qvx(im,j,k)) &
                /dxt ) &
                +((wrkqpy(i,jj,k,m)*qvy(i,jj,k)- wrkqpy(i,jm,k,m)*qvy(i,jm,k)) &
                /dyt ) &
                +((wrkqpy(i,j,kk,m)*qvz(i,j,kk)- wrkqpy(i,j,km,m)*qvz(i,j,km)) &
                /dzt ) )
                !
                aqpz=oldqpz(i,j,k,m)-delt*( &
                ((wrkqpz(ii,j,k,m)*qvx(ii,j,k)- wrkqpz(im,j,k,m)*qvx(im,j,k)) &
                /dxt ) &
                +((wrkqpz(i,jj,k,m)*qvy(i,jj,k)- wrkqpz(i,jm,k,m)*qvy(i,jm,k)) &
                /dyt ) &
                +((wrkqpz(i,j,kk,m)*qvz(i,j,kk)- wrkqpz(i,j,km,m)*qvz(i,j,km)) &
                /dzt ) )
                !
                !       add increments from spatial derivatives
                !
                abx=bsx(i,j,k)
                aby=bsy(i,j,k)
                abz=bsz(i,j,k)
                bmag=btot(i,j,k)
                weight=amax1(1.,bmag/bmax)
                skin_depth=reynolds/weight
                !
                tvcrossb_x=(tvy(i,j,k)*abz-tvz(i,j,k)*aby)
                tvcrossb_y=-(tvx(i,j,k)*abz-tvz(i,j,k)*abx)
                tvcrossb_z=(tvx(i,j,k)*aby-tvy(i,j,k)*abx)
                !
                qvcrossb_x=(qvy(i,j,k)*abz-qvz(i,j,k)*aby)
                qvcrossb_y=-(qvx(i,j,k)*abz-qvz(i,j,k)*abx)
                qvcrossb_z=(qvx(i,j,k)*aby-qvy(i,j,k)*abx)
                
                qgradp_x=(wrkqpresx(ii,j,k,m)-wrkqpresx(im,j,k,m))/dxt &
                    +(wrkqpresxy(i,jj,k,m)-wrkqpresxy(i,jm,k,m))/dxt &
                    +(wrkqpresxz(i,j,kk,m)-wrkqpresxz(i,j,km,m))/dxt
                qgradp_y=(wrkqpresxy(ii,j,k,m)-wrkqpresxy(im,j,k,m))/dxt &
                    +(wrkqpresy(i,jj,k,m)-wrkqpresy(i,jm,k,m))/dyt &
                    +(wrkqpresyz(i,j,kk,m)-wrkqpresyz(i,j,km,m))/dxt
                qgradp_z=(wrkqpresxz(ii,j,k,m)-wrkqpresxz(im,j,k,m))/dxt &
                    +(wrkqpresyz(i,jj,k,m)-wrkqpresyz(i,jm,k,m))/dyt &
                    +(wrkqpresz(i,j,kk,m)-wrkqpresz(i,j,km,m))/dzt

                !
                !       using straight electric field approach
                !
                dele_x=(efldx(i,j,k)+tvcrossb_x)*reynolds
                dele_y=(efldy(i,j,k)+tvcrossb_y)*reynolds
                dele_z=(efldz(i,j,k)+tvcrossb_z)*reynolds
                !
                delv_x=(qvcrossb_x-tvcrossb_x)*skin_depth
                delv_y=(qvcrossb_y-tvcrossb_y)*skin_depth
                delv_z=(qvcrossb_z-tvcrossb_z)*skin_depth
                !
                qden=wrkqrho(i,j,k,m)/rmassq
                qpx(i,j,k,m)=aqpx+delt* &
                ( qden*(dele_x+delv_x) -qgradp_x )
                qpy(i,j,k,m)=aqpy+delt* &
                ( qden*(dele_y+delv_y) -qgradp_y )
                qpz(i,j,k,m)=aqpz+delt* &
                ( qden*(dele_z+delv_z) -qgradp_z )
                !
                !
                !       add in jxb and grad p forces
                !       qjcrossb_x=(cury(i,j,k)*abz-curz(i,j,k)*aby)
                !       qjcrossb_y=-(curx(i,j,k)*abz-curz(i,j,k)*abx)
                !       qjcrossb_z=(curx(i,j,k)*aby-cury(i,j,k)*abx)
                !
                !       qpx(i,j,k,m)=aqpx+delt*(qjcrossb_x-qgradp_x)
                !       qpy(i,j,k,m)=aqpy+delt*(qjcrossb_y-qgradp_y)
                !       qpz(i,j,k,m)=aqpz+delt*(qjcrossb_z-qgradp_z)
                !
                !       add gravity
                !
                radius=sqrt(ax**2+ay**2+az**2)+0.0000001
                g=delt*wrkqrho(i,j,k,m)*grav/(re_equiv*radius)**2
                qpx(i,j,k,m)=qpx(i,j,k,m)-g*ax/radius
                qpy(i,j,k,m)=qpy(i,j,k,m)-g*ay/radius
                qpz(i,j,k,m)=qpz(i,j,k,m)-g*az/radius
   
                !
                !       pressure equations: isotropic components
                !
                !
                qpresx(i,j,k,m)=oldqpresx(i,j,k,m)-delt* &
                    ( ( (wrkqpresx(ii,j,k,m)*qvx(ii,j,k) &
                    -wrkqpresx(im,j,k,m)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresx(i,jj,k,m)*qvy(i,jj,k) &
                    -wrkqpresx(i,jm,k,m)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresx(i,j,kk,m)*qvz(i,j,kk) &
                    -wrkqpresx(i,j,km,m)*qvz(i,j,km))/dzt ) ) &
                    -2.*delt*wrkqpresx(i,j,k,m) &
                    *(qvx(ii,j,k)-qvx(im,j,k))/dxt
                !
                qpresy(i,j,k,m)=oldqpresy(i,j,k,m)-delt* &
                    ( ( (wrkqpresy(ii,j,k,m)*qvx(ii,j,k) &
                    -wrkqpresy(im,j,k,m)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresy(i,jj,k,m)*qvy(i,jj,k) &
                    -wrkqpresy(i,jm,k,m)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresy(i,j,kk,m)*qvz(i,j,kk) &
                    -wrkqpresy(i,j,km,m)*qvz(i,j,km))/dzt ) ) &
                    -2.*delt*wrkqpresy(i,j,k,m) &
                    * (qvy(i,jj,k)-qvy(i,jm,k))/dyt
                !
                qpresz(i,j,k,m)=oldqpresz(i,j,k,m)-delt* &
                    ( ( (wrkqpresz(ii,j,k,m)*qvx(ii,j,k) &
                    -wrkqpresz(im,j,k,m)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresz(i,jj,k,m)*qvy(i,jj,k) &
                    -wrkqpresz(i,jm,k,m)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresz(i,j,kk,m)*qvz(i,j,kk) &
                    -wrkqpresz(i,j,km,m)*qvz(i,j,km))/dzt ) ) &
                    -2.*delt*wrkqpresz(i,j,k,m) &
                    *(qvz(i,j,kk)-qvz(i,j,km))/dzt
                !
                !       add in anisotropic components
                !
                abx=(6.*bsx(i,j,k)+bsx(i+1,j,k)+bsx(i-1,j,k) &
                    +bsx(i,j+1,k)+bsx(i,j-1,k) &
                    +bsx(i,j,k+1)+bsx(i,j,k-1))/12.
                    aby=(6.*bsy(i,j,k)+bsy(i+1,j,k)+bsy(i-1,j,k) &
                    +bsy(i,j+1,k)+bsy(i,j-1,k) &
                    +bsy(i,j,k+1)+bsy(i,j,k-1))/12.
                    abz=(6.*bsz(i,j,k)+bsz(i+1,j,k)+bsz(i-1,j,k) &
                    +bsz(i,j+1,k)+bsz(i,j-1,k) &
                    +bsz(i,j,k+1)+bsz(i,j,k-1))/12.
                bmag=sqrt(abx**2+aby**2+abz**2)
                !
                skin_factor=reynolds
                if(bmag.gt.bmax)skin_factor=reynolds*(bmax/bmag)
                !
                qpresx(i,j,k,m)=qpresx(i,j,k,m) &
                    -2.*delt*wrkqpresxy(i,j,k,m) &
                    *(qvx(i,jj,k)-qvx(i,jm,k))/dyt &
                    -2.*delt*wrkqpresxz(i,j,k,m) &
                    *(qvx(i,j,kk)-qvx(i,j,km))/dzt &
                    +2.*delt*skin_factor/rmass*( &
                    wrkqpresxy(i,j,k,m)*abz &
                    -wrkqpresxz(i,j,k,m)*aby )
        
                qpresy(i,j,k,m)=qpresy(i,j,k,m) &
                    -2.*delt*wrkqpresxy(i,j,k,m) &
                    * (qvy(ii,j,k)-qvy(im,j,k))/dxt &
                    -2.*delt*wrkqpresyz(i,j,k,m) &
                    * (qvy(i,j,kk)-qvy(i,j,km))/dzt &
                    +2.*delt*skin_factor/rmass*( &
                    wrkqpresyz(i,j,k,m)*abx &
                    -wrkqpresxy(i,j,k,m)*abz )
                !
                qpresz(i,j,k,m)=qpresz(i,j,k,m) &
                    -2.*delt*wrkqpresxz(i,j,k,m) &
                    *(qvz(ii,j,k)-qvz(im,j,k))/dxt &
                    -2.*delt*wrkqpresyz(i,j,k,m) &
                    *(qvz(i,jj,k)-qvz(i,jm,k))/dyt &
                    +2.*delt*skin_factor/rmass*( &
                    wrkqpresxz(i,j,k,m)*aby &
                    -wrkqpresyz(i,j,k,m)*abx )
                !
                !       pressure equations: offdiagonal elements
                !
                qpresxy(i,j,k,m)=oldqpresxy(i,j,k,m)-delt* &
                    ( ( (wrkqpresxy(ii,j,k,m)*qvx(ii,j,k) &
                    -wrkqpresxy(im,j,k,m)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresxy(i,jj,k,m)*qvy(i,jj,k) &
                    -wrkqpresxy(i,jm,k,m)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresxy(i,j,kk,m)*qvz(i,j,kk) &
                    -wrkqpresxy(i,j,km,m)*qvz(i,j,km))/dzt ) )

                qpresxy(i,j,k,m)=qpresxy(i,j,k,m) &
                    -delt*wrkqpresxy(i,j,k,m) &
                    *( (qvx(ii,j,k)-qvx(im,j,k))/dxt &
                    +(qvy(i,jj,k)-qvy(i,jm,k))/dyt ) &
                    -delt*wrkqpresy(i,j,k,m) &
                    *(qvx(i,jj,k)-qvx(i,jm,k))/dyt &
                    -delt*wrkqpresx(i,j,k,m) &
                    *(qvy(ii,j,k)-qvy(im,j,k))/dxt &
                    -delt*wrkqpresyz(i,j,k,m) &
                    *(qvx(i,j,kk)-qvx(i,j,km))/dzt &
                    -delt*wrkqpresxz(i,j,k,m) &
                    *(qvy(i,j,kk)-qvy(i,j,km))/dzt &
                    +delt*skin_factor/rmass*( &
                    abz*(wrkqpresy(i,j,k,m)-wrkqpresx(i,j,k,m)) &
                    -aby*wrkqpresyz(i,j,k,m)+abx*wrkqpresxz(i,j,k,m))
                !
                qpresxz(i,j,k,m)=oldqpresxz(i,j,k,m)-delt* &
                    ( ( (wrkqpresxz(ii,j,k,m)*qvx(ii,j,k) &
                    -wrkqpresxz(im,j,k,m)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresxz(i,jj,k,m)*qvy(i,jj,k) &
                    -wrkqpresxz(i,jm,k,m)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresxz(i,j,kk,m)*qvz(i,j,kk) &
                    -wrkqpresxz(i,j,km,m)*qvz(i,j,km))/dzt ) )
                qpresxz(i,j,k,m)=qpresxz(i,j,k,m) &
                    -delt*wrkqpresxz(i,j,k,m) &
                    *( (qvx(ii,j,k)-qvx(im,j,k))/dxt &
                    +(qvz(i,j,kk)-qvz(i,j,km))/dzt ) &
                    -delt*wrkqpresz(i,j,k,m) &
                    *(qvx(i,j,kk)-qvx(i,j,km))/dzt &
                    -delt*wrkqpresx(i,j,k,m) &
                    *(qvz(ii,j,k)-qvz(im,j,k))/dxt &
                    -delt*wrkqpresyz(i,j,k,m) &
                    *(qvx(i,jj,k)-qvx(i,jm,k))/dyt &
                    -delt*wrkqpresxy(i,j,k,m) &
                    *(qvz(i,jj,k)-qvz(i,jm,k))/dyt &
                    +delt*skin_factor/rmass*( &
                    aby*(wrkqpresx(i,j,k,m)-wrkqpresz(i,j,k,m)) &
                    -abx*wrkqpresxy(i,j,k,m)+abz*wrkqpresyz(i,j,k,m))
                !
                qpresyz(i,j,k,m)=oldqpresyz(i,j,k,m)-delt* &
                    ( ( (wrkqpresyz(ii,j,k,m)*qvx(ii,j,k) &
                    -wrkqpresyz(im,j,k,m)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresyz(i,jj,k,m)*qvy(i,jj,k) &
                    -wrkqpresyz(i,jm,k,m)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresyz(i,j,kk,m)*qvz(i,j,kk) &
                    -wrkqpresyz(i,j,km,m)*qvz(i,j,km))/dzt ) )
                qpresyz(i,j,k,m)=qpresyz(i,j,k,m) &
                    -delt*wrkqpresyz(i,j,k,m) &
                    *( (qvy(i,jj,k)-qvy(i,jm,k))/dyt &
                    +(qvz(i,j,kk)-qvz(i,j,km))/dzt ) &
                    -delt*wrkqpresy(i,j,k,m) &
                    *(qvz(i,jj,k)-qvz(i,jm,k))/dyt &
                    -delt*wrkqpresz(i,j,k,m) &
                    *(qvy(i,j,kk)-qvy(i,j,km))/dzt &
                    -delt*wrkqpresxy(i,j,k,m) &
                    *(qvz(ii,j,k)-qvz(im,j,k))/dxt &
                    -delt*wrkqpresxz(i,j,k,m) &
                    *(qvy(ii,j,k)-qvy(im,j,k))/dxt &
                    +delt*skin_factor/rmass*( &
                    abx*(wrkqpresz(i,j,k,m)-wrkqpresy(i,j,k,m)) &
                    -abz*wrkqpresxz(i,j,k,m)+aby*wrkqpresxy(i,j,k,m))
                !
                !
                !
            enddo             ! k loop
        enddo             ! j loop
        !
    enddo             ! i loop
    
    return
end

!
!      **************************
!
subroutine push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
    efldx,efldy,efldz,nx,ny,nz,ngrd,m,delt, &
    rx,ry,rz)
    !
    !      standard runge-kutta time step
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd),oldbz(nx,ny,nz,ngrd), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz)
    !
    !       set physical grid spacing
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    !
    ! parallelizes loop rw, aug. 17, 2004
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kk=k+1
        do j=2,ny-1
            jm=j-1
            jj=j+1
    
            do i=2,nx-1
                ii=i+1
                im=i-1
                !
                !       induction equation
                !
                bx(i,j,k,m)=oldbx(i,j,k,m)-delt*( &
                (  (efldz(i,jj,k)-efldz(i,jm,k))/dyt ) &
                -( (efldy(i,j,kk)-efldy(i,j,km))/dzt )  )
                !
                by(i,j,k,m)=oldby(i,j,k,m)+delt*( &
                (  (efldz(ii,j,k)-efldz(im,j,k))/dxt ) &
                -( (efldx(i,j,kk)-efldx(i,j,km))/dzt )  )
                !
                bz(i,j,k,m)=oldbz(i,j,k,m)-delt*( &
                (  (efldy(ii,j,k)-efldy(im,j,k))/dxt) &
                -( (efldx(i,jj,k)-efldx(i,jm,k))/dyt )     )
                !
            enddo             ! k loop
        enddo             ! j loop
        !
    enddo             ! i loop
    !
    return
end

!
!     ********************************************
!
subroutine calcur(bx,by,bz,nx,ny,nz,ngrd,m,curx,cury,curz, &
    rx,ry,rz)
    !
    !     this calculates the current associated with the perturbed b field
    !          i.e. curz= dby/dx - dbx/dy
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz)
    
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kp=k+1
    
        do j=2,ny-1
            jm=j-1
            jp=j+1
    
            do i=2,nx-1
                im=i-1
                ip=i+1
    
                !
                curx(i,j,k)=(bz(i,jp,k,m)-bz(i,jm,k,m))/dyt &
                - (by(i,j,kp,m)-by(i,j,km,m))/dzt
                !
                cury(i,j,k)=(bx(i,j,kp,m)-bx(i,j,km,m))/dzt &
                - (bz(ip,j,k,m)-bz(im,j,k,m))/dxt
                !
                curz(i,j,k)=(by(ip,j,k,m)-by(im,j,k,m))/dxt &
                - (bx(i,jp,k,m)-bx(i,jm,k,m))/dyt
            enddo
        enddo
    enddo
    !
    !     following boundary conditions are set so that no forward
    !      communication , i.e. same x required
    !      and symmetry between k=1 and k=2
    !
    !     set boundary regions - bottom and top panels
    !
    nx1=nx-1
    ny1=ny-1
    nz1=nz-1
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do j=2,ny1
        do i=2,nx1
            !
            !       regular boundary conditions
            !
            curx(i,j,1)=curx(i,j,2)
            curx(i,j,nz)=curx(i,j,nz1)
            !
            cury(i,j,1)=cury(i,j,2)
            cury(i,j,nz)=cury(i,j,nz1)
            !
            curz(i,j,1)=curz(i,j,2)
            curz(i,j,nz)=curz(i,j,nz1)
        enddo
    enddo
    !
    !       set boundary regions - front and back
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz1
        do j=2,ny1
            !
            curx(1,j,k)=curx(2,j,k)
            curx(nx,j,k)=curx(nx1,j,k)
            !
            cury(1,j,k)=cury(2,j,k)
            cury(nx,j,k)=cury(nx1,j,k)
            !
            curz(1,j,k)=curz(2,j,k)
            curz(nx,j,k)=curz(nx1,j,k)
        enddo
    enddo
    !
    !       set boundary regions - left and right
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do j=2,nz1
        do i=2,nx1
            curx(i,1,j)=curx(i,2,j)
            curx(i,ny,j)=curx(i,ny1,j)
            !
            cury(i,1,j)=cury(i,2,j)
            cury(i,ny,j)=cury(i,ny1,j)
            !
            curz(i,1,j)=curz(i,2,j)
            curz(i,ny,j)=curz(i,ny1,j)
        enddo
    enddo
    !
    !     set corner lines
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do i=2,nx1
        curx(i,1,1)=curx(i,1,2)
        curx(i,1,nz)=curx(i,1,nz1)
        cury(i,1,1)=cury(i,1,2)
        cury(i,1,nz)=cury(i,1,nz1)
        curz(i,1,1)=curz(i,1,2)
        curz(i,1,nz)=curz(i,1,nz1)
        !
        curx(i,ny,1)=curx(i,ny,2)
        curx(i,ny,nz)=curx(i,ny1,nz)
        cury(i,ny,1)=cury(i,ny,2)
        cury(i,ny,nz)=cury(i,ny1,nz)
        curz(i,ny,1)=curz(i,ny,2)
        curz(i,ny,nz)=curz(i,ny1,nz)
    enddo
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do j=2,ny1
        curx(1,j,1)=curx(1,j,2)
        curx(1,j,nz)=curx(1,j,nz1)
        cury(1,j,1)=cury(1,j,2)
        cury(1,j,nz)=cury(1,j,nz1)
        curz(1,j,1)=curz(1,j,2)
        curz(1,j,nz)=curz(1,j,nz1)
        !
        curx(nx,j,1)=curx(nx,j,2)
        curx(nx,j,nz)=curx(nx,j,nz1)
        cury(nx,j,1)=cury(nx,j,2)
        cury(nx,j,nz)=cury(nx,j,nz1)
        curz(nx,j,1)=curz(nx,j,2)
        curz(nx,j,nz)=curz(nx,j,nz1)
    enddo
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz1
        curx(1,1,k)=curx(1,2,k)
        curx(nx,1,k)=curx(nx,2,k)
        cury(1,1,k)=cury(1,2,k)
        cury(nx,1,k)=cury(nx,2,k)
        curz(1,1,k)=curz(1,2,k)
        curz(nx,1,k)=curz(nx,2,k)
        !
        curx(1,ny,k)=curx(1,ny1,k)
        curx(nx,ny,k)=curx(nx,ny1,k)
        cury(1,ny,k)=cury(1,ny1,k)
        cury(nx,ny,k)=cury(nx,ny1,k)
        curz(1,ny,k)=curz(1,ny1,k)
        curz(nx,ny,k)=curz(nx,ny1,k)
        !
    enddo
    !
    curx(1,1,1)=curx(1,1,2)
    curx(1,ny,1)=curx(1,ny,2)
    curx(1,1,nz)=(curx(1,1,nz1)+curx(1,2,nz))/2.
    curx(1,ny,nz)=(curx(1,ny,nz1)+curx(1,ny1,nz))/2.
    curx(nx,1,1)=(curx(nx,1,2)+curx(nx,2,1)+curx(nx1,1,1))/3.
    curx(nx,ny,1)=(curx(nx,ny,2)+curx(nx,ny1,1)+curx(nx1,ny,1))/3.
    curx(nx,1,nz)=(curx(nx,1,nz1)+curx(nx,2,nz)+curx(nx1,1,nz))/3.
    curx(nx,ny,nz)=(curx(nx,ny,nz1)+curx(nx,ny1,nz) &
    +curx(nx1,ny,nz))/3.
    !
    cury(1,1,1)=cury(1,1,2)
    cury(1,ny,1)=cury(1,ny,2)
    cury(1,1,nz)=(cury(1,1,nz1)+cury(1,2,nz))/2.
    cury(1,ny,nz)=(cury(1,ny,nz1)+cury(1,ny1,nz))/2.
    cury(nx,1,1)=(cury(nx,1,2)+cury(nx,2,1)+cury(nx1,1,1))/3.
    cury(nx,ny,1)=(cury(nx,ny,2)+cury(nx,ny1,1)+cury(nx1,ny,1))/3.
    cury(nx,1,nz)=(cury(nx,1,nz1)+cury(nx,2,nz)+cury(nx1,1,nz))/3.
    cury(nx,ny,nz)=(cury(nx,ny,nz1)+cury(nx,ny1,nz) &
    +cury(nx1,ny,nz))/3.
    !
    curz(1,1,1)=curz(1,1,2)
    curz(1,ny,1)=curz(1,ny,2)
    curz(1,1,nz)=(curz(1,1,nz1)+curz(1,2,nz))/2.
    curz(1,ny,nz)=(curz(1,ny,nz1)+curz(1,ny1,nz))/2.
    curz(nx,1,1)=(curz(nx,1,2)+curz(nx,2,1)+curz(nx1,1,1))/3.
    curz(nx,ny,1)=(curz(nx,ny,2)+curz(nx,ny1,1)+curz(nx1,ny,1))/3.
    curz(nx,1,nz)=(curz(nx,1,nz1)+curz(nx,2,nz)+curz(nx1,1,nz))/3.
    curz(nx,ny,nz)=(curz(nx,ny,nz1)+curz(nx,ny1,nz) &
    +curz(nx1,ny,nz))/3.
    !
    return
end

!
!     *****************************************
!
subroutine space_charge(charge,efldx,efldy,efldz, &
    bx0,by0,bz0,px,py,pz,rho,nx,ny,nz,ngrd,m,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      calculates the space-charge fields in the plasma
    !
    dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd), &
    bz0(nx,ny,nz,ngrd),charge(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd), &
    rho(nx,ny,nz,ngrd)
    !
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !
    rx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    ry=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    rz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    !
    d_min=0.001
    !
    !     take the divergence of e and correct
    !          for curl b_0 errors
    !
    do k=2,nz-1
        km=k-1
        kp=k+1
        do j=2,ny-1
            jm=j-1
            jp=j+1
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                !     div. e
                !
                charge(i,j,k)= &
                (efldx(ip,j,k)-efldx(im,j,k))/dxt &
                +(efldy(i,jp,k)-efldy(i,jm,k))/dyt &
                +(efldz(i,j,kp)-efldz(i,j,km))/dzt
                !
                !     correct for curl b_0 errors
                !
                arho=rho(i,j,k,m)
                arho=amax1(arho,d_min)
                !
                apx=px(i,j,k,m)/arho
                apy=py(i,j,k,m)/arho
                apz=pz(i,j,k,m)/arho
                !
                curlbx=( (bz0(i,jp,k,m)-bz0(i,jm,k,m))/dyt &
                -(by0(i,j,kp,m)-by0(i,j,km,m))/dzt )/arho
                curlby=( (bx0(i,j,kp,m)-bx0(i,j,km,m))/dzt &
                -(bz0(ip,j,k,m)-bz0(im,j,k,m))/dxt )/arho
                curlbz=( (by0(ip,j,k,m)-by0(im,j,k,m))/dxt &
                -(bx0(i,jp,k,m)-bx0(i,jm,k,m))/dyt )/arho
                !
                charge(i,j,k)=charge(i,j,k) &
                -(apx*curlbx+apy*curlby+apz*curlbz)
                !
            enddo
        enddo
    enddo
    !
    !      set flank bounday conditions
    !
    nx1=nx-1
    ny1=ny-1
    nz1=nz-1
    !
    do k=1,nz
        do j=1,ny
            charge(1,j,k)=charge(2,j,k)
            charge(nx,j,k)=charge(nx1,j,k)
        enddo
    enddo
    !
    do j=1,ny
        do i=1,nx
            charge(i,j,1)=charge(i,j,2)
            charge(i,j,nz)=charge(i,j,nz1)
        enddo
    enddo
    !
    do k=1,nz
        do i=1,nx
            charge(i,1,k)=charge(i,2,k)
            charge(i,ny,k)=charge(i,ny1,k)
        enddo
    enddo
    !
    !
    !
    return
end

!
!     *****************************************
!
subroutine bande(efldx,efldy,efldz,bsx,bsy,bsz, &
    curx,cury,curz,evx,evy,evz,btot, &
    epres,qrho,hrho,orho,rst,resist,reynolds, &
    nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
    ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
    rx,ry,rz)
    !
    !      calculates the surface magnetic field bs
    !      and the body electric field eb
    !
    dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz), &
    qrho(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd),btot(nx,ny,nz)
    dimension rst(nx,ny,nz,mbndry)
    integer ijmid(mbndry,3,mmid),ijzero(mbndry,3,mzero)
    integer nummid(mbndry),numzero(mbndry)
    !
    !
    !      ohm's law: ve = vi-j/ne
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                !
                avx=evx(i,j,k)
                avy=evy(i,j,k)
                avz=evz(i,j,k)
                !
                abx=bsx(i,j,k)
                aby=bsy(i,j,k)
                abz=bsz(i,j,k)
                !
                efldx(i,j,k)=-(avy*abz-avz*aby)
                efldy(i,j,k)=-(avz*abx-avx*abz)
                efldz(i,j,k)=-(avx*aby-avy*abx)
                !
            enddo
        enddo
    enddo
    !
    !     add in grad p term
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
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
                arho=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
                orho(i,j,k,m)/rmasso)*reynolds
                !
                efldx(i,j,k)=efldx(i,j,k)- &
                ((epres(ip,j,k,m)-epres(im,j,k,m))/dxt)/arho
                efldy(i,j,k)=efldy(i,j,k)- &
                ((epres(i,jp,k,m)-epres(i,jm,k,m))/dyt)/arho
                efldz(i,j,k)=efldz(i,j,k)- &
                ((epres(i,j,kp,m)-epres(i,j,km,m))/dzt)/arho
                !
            enddo
        enddo
    enddo
    !
    !      add in ionospheric resistance
    !
    if((m.le.mbndry).and.(resist.lt.5000.))then
        ! parallelizes loop rw, oct. 23, 2002
        !$omp  parallel do
        do n=1,nummid(m)
            i=ijmid(m,1,n)
            j=ijmid(m,2,n)
            k=ijmid(m,3,n)
            !
            eden=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
            orho(i,j,k,m)/rmasso)
            !
            efldx(i,j,k)=efldx(i,j,k)+rst(i,j,k,m)*curx(i,j,k)/eden
            efldy(i,j,k)=efldy(i,j,k)+rst(i,j,k,m)*cury(i,j,k)/eden
            efldz(i,j,k)=efldz(i,j,k)+rst(i,j,k,m)*curz(i,j,k)/eden
        enddo
        do n=1,numzero(m)
            i=ijzero(m,1,n)
            j=ijzero(m,2,n)
            k=ijzero(m,3,n)
            !
            eden=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
            orho(i,j,k,m)/rmasso)
            !
            efldx(i,j,k)=efldx(i,j,k)+rst(i,j,k,m)*curx(i,j,k)/eden
            efldy(i,j,k)=efldy(i,j,k)+rst(i,j,k,m)*cury(i,j,k)/eden
            efldz(i,j,k)=efldz(i,j,k)+rst(i,j,k,m)*curz(i,j,k)/eden
        enddo
    endif
    !
    !      boundary condtions for terrestrial surface currents
    !
    !     if(m.le.mbndry)then
    !      do n=1,nummid(m)
    !        i=ijmid(m,1,n)
    !        j=ijmid(m,2,n)
    !        k=ijmid(m,3,n)
    !        efldx(i,j,k)=0.
    !        efldy(i,j,k)=0.
    !        efldz(i,j,k)=0.
    !       enddo
    !      endif
    !
    !      set flank bounday conditions
    !
    nx1=nx-1
    ny1=ny-1
    nz1=nz-1
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            efldx(1,j,k)=efldx(2,j,k)
            efldx(nx,j,k)=efldx(nx1,j,k)
            efldy(1,j,k)=efldy(2,j,k)
            efldy(nx,j,k)=efldy(nx1,j,k)
            efldz(1,j,k)=efldz(2,j,k)
            efldz(nx,j,k)=efldz(nx1,j,k)
        enddo
    enddo
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            efldx(i,j,1)=efldx(i,j,2)
            efldx(i,j,nz)=efldx(i,j,nz1)
            efldy(i,j,1)=efldy(i,j,2)
            efldy(i,j,nz)=efldy(i,j,nz1)
            efldz(i,j,1)=efldz(i,j,2)
            efldz(i,j,nz)=efldz(i,j,nz1)
        enddo
    enddo
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            efldx(i,1,k)=efldx(i,2,k)
            efldx(i,ny,k)=efldx(i,ny1,k)
            efldy(i,1,k)=efldy(i,2,k)
            efldy(i,ny,k)=efldy(i,ny1,k)
            efldz(i,1,k)=efldz(i,2,k)
            efldz(i,ny,k)=efldz(i,ny1,k)
        enddo
    enddo
    !
    return
end

!
!     ***********************************************
!
subroutine qvset(a,v,n)
    !
    !     set an array v of total dimension n to value a
    !
    dimension v(1)
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do i=1,n
        v(i)=a
    enddo
    return
end

!
!     ********************************************
!
subroutine lap_test(qrho,qpres,qpx,qpy,qpz, &
    wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
    hrho,hpres,hpx,hpy,hpz, &
    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
    orho,opres,opx,opy,opz, &
    wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt, &
    rmassq,rmassh,rmasso)
    
    dimension qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd)
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpres(nx,ny,nz,ngrd), &
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd),wrkhrho(nx,ny,nz,ngrd), &
    wrkhpres(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd), &
    wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd), &
    wrkorho(nx,ny,nz,ngrd),wrkopres(nx,ny,nz,ngrd), &
    wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
    wrkopz(nx,ny,nz,ngrd)
    !
    write(6,*)'test lap entered'
    return
end

!
!     ********************************************
!
subroutine lap_plasma(rho,px,py,pz, &
    presx,presy,presz,presxy,presxz,presyz, &
    wrkrho,wrkpx,wrkpy,wrkpz, &
    wrkpresx,wrkpresy,wrkpresz, &
    wrkpresxy,wrkpresxz,wrkpresyz, &
    vx,vy,vz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt,isotropic)
    !
    !     apply the ladipus smoothing technique for particular ion component
    !
    dimension rho(nx,ny,nz,ngrd),presx(nx,ny,nz,ngrd), &
    presy(nx,ny,nz,ngrd),presz(nx,ny,nz,ngrd), &
    px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd), &
    wrkrho(nx,ny,nz,ngrd),wrkpresx(nx,ny,nz,ngrd), &
    wrkpresy(nx,ny,nz,ngrd),wrkpresz(nx,ny,nz,ngrd), &
    wrkpx(nx,ny,nz,ngrd),wrkpy(nx,ny,nz,ngrd), &
    wrkpz(nx,ny,nz,ngrd), &
    vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    dimension presxy(nx,ny,nz,ngrd),presxz(nx,ny,nz,ngrd), &
    presyz(nx,ny,nz,ngrd),wrkpresxy(nx,ny,nz,ngrd), &
    wrkpresxz(nx,ny,nz,ngrd),wrkpresyz(nx,ny,nz,ngrd)
    !
    
    common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9), &
    grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9), &
    rx,ry,rz,xdip,ydip,zdip,rearth,b0, &
    sin_tilt,cos_tilt
    !
    !     equal waiting irrespective of grid size
    !
    deltz=delt
    delty=delt
    deltx=delt
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do  k=2,nz-1
        km=k-1
        kp=k+1
        do  j=2,ny-1
            jm=j-1
            jp=j+1
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                uxp1=abs(vx(ip,j,k)-vx(i,j,k))
                uxm1=abs(vx(i,j,k)-vx(im,j,k))
                !
                uyp1=abs(vy(i,jp,k)-vy(i,j,k))
                uym1=abs(vy(i,j,k)-vy(i,jm,k))
                !
                uzp1=abs(vz(i,j,kp)-vz(i,j,k))
                uzm1=abs(vz(i,j,k)-vz(i,j,km))
                !
                wrkrho(i,j,k,m)=rho(i,j,k,m)+chirho*( &
                    deltx*(uxp1*(rho(ip,j,k,m)-rho(i,j,k,m)) &
                    -uxm1*(rho(i,j,k,m)-rho(im,j,k,m))) &
                    +delty*(uyp1*(rho(i,jp,k,m)-rho(i,j,k,m)) &
                    -uym1*(rho(i,j,k,m)-rho(i,jm,k,m))) &
                    +deltz*(uzp1*(rho(i,j,kp,m)-rho(i,j,k,m)) &
                    -uzm1*(rho(i,j,k,m)-rho(i,j,km,m))) )

                wrkpx(i,j,k,m)=px(i,j,k,m)+chipxyz*( &
                    +deltx*(uxp1*(px(ip,j,k,m)-px(i,j,k,m)) &
                    -uxm1*(px(i,j,k,m)-px(im,j,k,m))) &
                    +delty*(uyp1*(px(i,jp,k,m)-px(i,j,k,m)) &
                    -uym1*(px(i,j,k,m)-px(i,jm,k,m))) &
                    +deltz*(uzp1*(px(i,j,kp,m)-px(i,j,k,m)) &
                    -uzm1*(px(i,j,k,m)-px(i,j,km,m))) )
                !
                wrkpy(i,j,k,m)=py(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(py(ip,j,k,m)-py(i,j,k,m)) &
                    -uxm1*(py(i,j,k,m)-py(im,j,k,m))) &
                    +delty*(uyp1*(py(i,jp,k,m)-py(i,j,k,m)) &
                    -uym1*(py(i,j,k,m)-py(i,jm,k,m))) &
                    +deltz*(uzp1*(py(i,j,kp,m)-py(i,j,k,m)) &
                    -uzm1*(py(i,j,k,m)-py(i,j,km,m))) )
                !
                wrkpz(i,j,k,m)=pz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(pz(ip,j,k,m)-pz(i,j,k,m)) &
                    -uxm1*(pz(i,j,k,m)-pz(im,j,k,m))) &
                    +delty*(uyp1*(pz(i,jp,k,m)-pz(i,j,k,m)) &
                    -uym1*(pz(i,j,k,m)-pz(i,jm,k,m))) &
                    +deltz*(uzp1*(pz(i,j,kp,m)-pz(i,j,k,m)) &
                    -uzm1*(pz(i,j,k,m)-pz(i,j,km,m))) )

                !
                !      diagonal elements
                !
                wrkpresx(i,j,k,m)=presx(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presx(ip,j,k,m)-presx(i,j,k,m)) &
                    -uxm1*(presx(i,j,k,m)- presx(im,j,k,m))) &
                    +delty*(uyp1*(presx(i,jp,k,m)-presx(i,j,k,m)) &
                    -uym1*(presx(i,j,k,m)- presx(i,jm,k,m))) &
                    +deltz*(uzp1*(presx(i,j,kp,m)-presx(i,j,k,m)) &
                    -uzm1*(presx(i,j,k,m)- presx(i,j,km,m))) )
                !
                wrkpresy(i,j,k,m)=presy(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presy(ip,j,k,m)-presy(i,j,k,m)) &
                    -uxm1*(presy(i,j,k,m)- presy(im,j,k,m))) &
                    +delty*(uyp1*(presy(i,jp,k,m)-presy(i,j,k,m)) &
                    -uym1*(presy(i,j,k,m)- presy(i,jm,k,m))) &
                    +deltz*(uzp1*(presy(i,j,kp,m)-presy(i,j,k,m)) &
                    -uzm1*(presy(i,j,k,m)- presy(i,j,km,m))) )
                !
                wrkpresz(i,j,k,m)=presz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presz(ip,j,k,m)-presz(i,j,k,m)) &
                    -uxm1*(presz(i,j,k,m)- presz(im,j,k,m))) &
                    +delty*(uyp1*(presz(i,jp,k,m)-presz(i,j,k,m)) &
                    -uym1*(presz(i,j,k,m)- presz(i,jm,k,m))) &
                    +deltz*(uzp1*(presz(i,j,kp,m)-presz(i,j,k,m)) &
                    -uzm1*(presz(i,j,k,m)- presz(i,j,km,m))) )
                !
                !      off-diagonal elements
                !
                wrkpresxy(i,j,k,m)=presxy(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presxy(ip,j,k,m)-presxy(i,j,k,m)) &
                    -uxm1*(presxy(i,j,k,m)- presxy(im,j,k,m))) &
                    +delty*(uyp1*(presxy(i,jp,k,m)-presxy(i,j,k,m)) &
                    -uym1*(presxy(i,j,k,m)- presxy(i,jm,k,m))) &
                    +deltz*(uzp1*(presxy(i,j,kp,m)-presxy(i,j,k,m)) &
                    -uzm1*(presxy(i,j,k,m)- presxy(i,j,km,m))) )
                !
                wrkpresxz(i,j,k,m)=presxz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presxz(ip,j,k,m)-presxz(i,j,k,m)) &
                    -uxm1*(presxz(i,j,k,m)- presxz(im,j,k,m))) &
                    +delty*(uyp1*(presxz(i,jp,k,m)-presxz(i,j,k,m)) &
                    -uym1*(presxz(i,j,k,m)- presxz(i,jm,k,m))) &
                    +deltz*(uzp1*(presxz(i,j,kp,m)-presxz(i,j,k,m)) &
                    -uzm1*(presxz(i,j,k,m)- presxz(i,j,km,m))) )
                !
                wrkpresyz(i,j,k,m)=presyz(i,j,k,m)+chipxyz*( &
                    deltx*(uxp1*(presyz(ip,j,k,m)-presyz(i,j,k,m)) &
                    -uxm1*(presyz(i,j,k,m)- presyz(im,j,k,m))) &
                    +delty*(uyp1*(presyz(i,jp,k,m)-presyz(i,j,k,m)) &
                    -uym1*(presyz(i,j,k,m)- presyz(i,jm,k,m))) &
                    +deltz*(uzp1*(presyz(i,j,kp,m)-presyz(i,j,k,m)) &
                    -uzm1*(presyz(i,j,k,m)- presyz(i,j,km,m))) )
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine lap_elec(ppres,wrkppres,vx,vy,vz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt)
    !
    !     apply the ladipus smoothing technique for particular ion component
    !
    dimension ppres(nx,ny,nz,ngrd),wrkppres(nx,ny,nz,ngrd), &
    vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !     equal waiting irrespective of grid size
    !
    deltz=delt
    delty=delt
    deltx=delt
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kp=k+1
        do j=2,ny-1
            jm=j-1
            jp=j+1
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                uxp1=abs(vx(ip,j,k)-vx(i,j,k))
                uxm1=abs(vx(i,j,k)-vx(im,j,k))
                !
                uyp1=abs(vy(i,jp,k)-vy(i,j,k))
                uym1=abs(vy(i,j,k)-vy(i,jm,k))
                !
                uzp1=abs(vz(i,j,kp)-vz(i,j,k))
                uzm1=abs(vz(i,j,k)-vz(i,j,km))
                !
                wrkppres(i,j,k,m)=ppres(i,j,k,m)+chipxyz*( &
                deltx*(uxp1*(ppres(ip,j,k,m)-ppres(i,j,k,m)) &
                -uxm1*(ppres(i,j,k,m)-ppres(im,j,k,m))) &
                +delty*(uyp1*(ppres(i,jp,k,m)-ppres(i,j,k,m)) &
                -uym1*(ppres(i,j,k,m)-ppres(i,jm,k,m))) &
                +deltz*(uzp1*(ppres(i,j,kp,m)-ppres(i,j,k,m)) &
                -uzm1*(ppres(i,j,k,m)-ppres(i,j,km,m))) )
            enddo
        enddo
    enddo
    !
    !
    return
end

!
!     ********************************************
!
subroutine lap_bfld(bx,by,bz,wrkbx,wrkby,wrkbz,vx,vy,vz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt)
    !
    !     apply the ladipus smoothing technique for particular ion component
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
    wrkbz(nx,ny,nz,ngrd), &
    vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !     equal waiting irrespective of grid size
    !
    deltz=delt
    delty=delt
    deltx=delt
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kp=k+1
        do j=2,ny-1
            jm=j-1
            jp=j+1
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                uxp1=abs(vx(ip,j,k)-vx(i,j,k))
                uxm1=abs(vx(i,j,k)-vx(im,j,k))
                !
                uyp1=abs(vy(i,jp,k)-vy(i,j,k))
                uym1=abs(vy(i,j,k)-vy(i,jm,k))
                !
                uzp1=abs(vz(i,j,kp)-vz(i,j,k))
                uzm1=abs(vz(i,j,k)-vz(i,j,km))
                !
                wrkbx(i,j,k,m)=bx(i,j,k,m)+chierg*( &
                deltx*(uxp1*(bx(ip,j,k,m)-bx(i,j,k,m)) &
                -uxm1*(bx(i,j,k,m)-bx(im,j,k,m))) &
                +delty*(uyp1*(bx(i,jp,k,m)-bx(i,j,k,m)) &
                -uym1*(bx(i,j,k,m)-bx(i,jm,k,m))) &
                +deltz*(uzp1*(bx(i,j,kp,m)-bx(i,j,k,m)) &
                -uzm1*(bx(i,j,k,m)-bx(i,j,km,m))) )
                wrkby(i,j,k,m)=by(i,j,k,m)+chierg*( &
                deltx*(uxp1*(by(ip,j,k,m)-by(i,j,k,m)) &
                -uxm1*(by(i,j,k,m)-by(im,j,k,m))) &
                +delty*(uyp1*(by(i,jp,k,m)-by(i,j,k,m)) &
                -uym1*(by(i,j,k,m)-by(i,jm,k,m))) &
                +deltz*(uzp1*(by(i,j,kp,m)-by(i,j,k,m)) &
                -uzm1*(by(i,j,k,m)-by(i,j,km,m))) )
                wrkbz(i,j,k,m)=bz(i,j,k,m)+chierg*( &
                deltx*(uxp1*(bz(ip,j,k,m)-bz(i,j,k,m)) &
                -uxm1*(bz(i,j,k,m)-bz(im,j,k,m))) &
                +delty*(uyp1*(bz(i,jp,k,m)-bz(i,j,k,m)) &
                -uym1*(bz(i,j,k,m)-bz(i,jm,k,m))) &
                +deltz*(uzp1*(bz(i,j,kp,m)-bz(i,j,k,m)) &
                -uzm1*(bz(i,j,k,m)-bz(i,j,km,m))) )
            enddo
        enddo
    enddo
    !
    !
    return
end

!
!     ********************************************
!
subroutine flux_correct(qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
    wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
    wrkqpresxy,wrkqpresxz,wrkqpresyz, &
    wrkqpx,wrkqpy,wrkqpz, &
    oldqrho,oldqpresx,oldqpresy,oldqpresz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    oldqpx,oldqpy,oldqpz, &
    !
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
    wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
    wrkhpresxy,wrkhpresxz,wrkhpresyz, &
    wrkhpx,wrkhpy,wrkhpz, &
    oldhrho,oldhpresx,oldhpresy,oldhpresz, &
    oldhpresxy,oldhpresxz,oldhpresyz, &
    oldhpx,oldhpy,oldhpz, &
    !
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,opx,opy,opz, &
    wrkorho,wrkopresx,wrkopresy,wrkopresz, &
    wrkopresxy,wrkopresxz,wrkopresyz, &
    wrkopx,wrkopy,wrkopz, &
    oldorho,oldopresx,oldopresy,oldopresz, &
    oldopresxy,oldopresxz,oldopresyz, &
    oldopx,oldopy,oldopz, &
    !
    epres,wrkepres,oldepres, &
    bx,by,bz,wrkbx,wrkby,wrkbz, &
    oldbx,oldby,oldbz,vvx,vvy,vvz, &
    nx,ny,nz,ngrd,m,difrho,diferg,xspac, &
    isotropic)
    !
    !     apply flux correction smoothing technique
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd), &
    bz(nx,ny,nz,ngrd),epres(nx,ny,nz,ngrd), &
    vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz), &
    xspac(ngrd)
    !
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpx(nx,ny,nz,ngrd), &
    wrkqpy(nx,ny,nz,ngrd),wrkqpz(nx,ny,nz,ngrd), &
    wrkqpresx(nx,ny,nz,ngrd),wrkqpresy(nx,ny,nz,ngrd), &
    wrkqpresz(nx,ny,nz,ngrd),wrkqpresxy(nx,ny,nz,ngrd), &
    wrkqpresxz(nx,ny,nz,ngrd),wrkqpresyz(nx,ny,nz,ngrd), &
    !
    wrkhrho(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd), &
    wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd), &
    wrkhpresx(nx,ny,nz,ngrd),wrkhpresy(nx,ny,nz,ngrd), &
    wrkhpresz(nx,ny,nz,ngrd),wrkhpresxy(nx,ny,nz,ngrd), &
    wrkhpresxz(nx,ny,nz,ngrd),wrkhpresyz(nx,ny,nz,ngrd), &
    !
    wrkorho(nx,ny,nz,ngrd),wrkopx(nx,ny,nz,ngrd), &
    wrkopy(nx,ny,nz,ngrd),wrkopz(nx,ny,nz,ngrd), &
    wrkopresx(nx,ny,nz,ngrd),wrkopresy(nx,ny,nz,ngrd), &
    wrkopresz(nx,ny,nz,ngrd),wrkopresxy(nx,ny,nz,ngrd), &
    wrkopresxz(nx,ny,nz,ngrd),wrkopresyz(nx,ny,nz,ngrd), &
    !
    wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
    wrkbz(nx,ny,nz,ngrd),wrkepres(nx,ny,nz,ngrd)
    !
    !
    dimension oldqrho(nx,ny,nz,ngrd),oldqpx(nx,ny,nz,ngrd), &
    oldqpy(nx,ny,nz,ngrd),oldqpz(nx,ny,nz,ngrd), &
    oldqpresx(nx,ny,nz,ngrd),oldqpresy(nx,ny,nz,ngrd), &
    oldqpresz(nx,ny,nz,ngrd),oldqpresxy(nx,ny,nz,ngrd), &
    oldqpresxz(nx,ny,nz,ngrd),oldqpresyz(nx,ny,nz,ngrd), &
    !
    oldhrho(nx,ny,nz,ngrd),oldhpx(nx,ny,nz,ngrd), &
    oldhpy(nx,ny,nz,ngrd),oldhpz(nx,ny,nz,ngrd), &
    oldhpresx(nx,ny,nz,ngrd),oldhpresy(nx,ny,nz,ngrd), &
    oldhpresz(nx,ny,nz,ngrd),oldhpresxy(nx,ny,nz,ngrd), &
    oldhpresxz(nx,ny,nz,ngrd),oldhpresyz(nx,ny,nz,ngrd), &
    !
    oldorho(nx,ny,nz,ngrd),oldopx(nx,ny,nz,ngrd), &
    oldopy(nx,ny,nz,ngrd),oldopz(nx,ny,nz,ngrd), &
    oldopresx(nx,ny,nz,ngrd),oldopresy(nx,ny,nz,ngrd), &
    oldopresz(nx,ny,nz,ngrd),oldopresxy(nx,ny,nz,ngrd), &
    oldopresxz(nx,ny,nz,ngrd),oldopresyz(nx,ny,nz,ngrd), &
    !
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
    oldbz(nx,ny,nz,ngrd),oldepres(nx,ny,nz,ngrd)
    !
    !
    !      write(6,*)'in flux correct with',nx,ny,nz,ngrd,m,chirho,
    !    +              diferg,xspac
    !
    rx=xspac(m)
    ry=xspac(m)
    rz=xspac(m)
    chifcs=difrho
    call fcsmooth(qrho,oldqrho,wrkqrho,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)

    call fcsmooth(qpresx,oldqpresx,wrkqpresx,nx,ny,nz,ngrd,m, &
        chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresy,oldqpresy,wrkqpresy,nx,ny,nz,ngrd,m, &
        chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresz,oldqpresz,wrkqpresz,nx,ny,nz,ngrd,m, &
        chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresxy,oldqpresxy,wrkqpresxy,nx,ny,nz,ngrd, &
        m,chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresxz,oldqpresxz,wrkqpresxz,nx,ny,nz,ngrd, &
        m,chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresyz,oldqpresyz,wrkqpresyz,nx,ny,nz,ngrd, &
        m,chifcs,vvx,vvy,vvz)

    call fcsmooth(qpx,oldqpx,wrkqpx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(qpy,oldqpy,wrkqpy,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(qpz,oldqpz,wrkqpz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    call fcsmooth(hrho,oldhrho,wrkhrho,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)

    call fcsmooth(hpresx,oldhpresx,wrkhpresx,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresy,oldhpresy,wrkhpresy,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresz,oldhpresz,wrkhpresz,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresxy,oldhpresxy,wrkhpresxy,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresxz,oldhpresxz,wrkhpresxz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresyz,oldhpresyz,wrkhpresyz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)

    call fcsmooth(hpx,oldhpx,wrkhpx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(hpy,oldhpy,wrkhpy,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(hpz,oldhpz,wrkhpz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    call fcsmooth(orho,oldorho,wrkorho,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)

    call fcsmooth(opresx,oldopresx,wrkopresx,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(opresy,oldopresy,wrkopresy,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(opresz,oldopresz,wrkopresz,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(opresxy,oldopresxy,wrkopresxy,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(opresxz,oldopresxz,wrkopresxz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(opresyz,oldopresyz,wrkopresyz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)

    call fcsmooth(opx,oldopx,wrkopx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(opy,oldopy,wrkopy,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(opz,oldopz,wrkopz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    !
    call fcsmooth(epres,oldepres,wrkepres,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    chifcs=diferg
    call fcsmooth(bx,oldbx,wrkbx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(by,oldby,wrkby,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(bz,oldbz,wrkbz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    !
    return
end


!
!     ********************************************
!
    subroutine fcsmooth_2d(px,oldpx,wrkpx,nx,ny,nz,ngrd,m,chipx, &
    delta,delta1,fn)
    !
    !     applies flux correction smoothing to field quantities
    !           oldpx is the qunatity at the start of the time step
    !           wrkpx is the unsmoothed qunatity
    !           px is the final product
    !     wrk arrays assumed to have dimension larger than nx,ny,nz
    !
    dimension px(nx,ny,nz,ngrd),oldpx(nx,ny,nz,ngrd), &
    wrkpx(nx,ny,nz,ngrd), &
    delta(nx,ny,nz),delta1(nx,ny,nz),fn(nx,ny,nz)
    !
    !         step 1:   diffuse initial in x direction
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                oldpx(i+1,j,k,m)+oldpx(i-1,j,k,m) &
                -2.*oldpx(i,j,k,m))
            enddo
        enddo
    enddo
    !
    !     set boundary points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            px(1,j,k,m)=wrkpx(1,j,k,m)
            px(nx,j,k,m)=wrkpx(nx,j,k,m)
        enddo
    enddo
    !
    !     calculate the difference  between x nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx-1
                i1=i+1
                delta(i,j,k)=wrkpx(i1,j,k,m)-wrkpx(i,j,k,m)   !d-non-diffuse n+1 sol'n
                delta1(i,j,k)=px(i1,j,k,m)-px(i,j,k,m)   !d1-diffused n+1 sol'n
            enddo
        enddo
    enddo
    !
    do k=1,nz
        do j=1,ny
            delta(nx,j,k)=delta(nx-1,j,k)
            delta1(nx,j,k)=delta1(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                s=sign(1.,delta(i,j,k))
                fn(i,j,k)=s*amax1(0.,amin1(s*delta1(i-1,j,k), &
                chipx*abs(delta(i,j,k)),s*delta1(i+1,j,k)))
            enddo
        enddo
    enddo
    !
    do k=1,nz
        do j=1,ny
            fn(1,j,k)=fn(2,j,k)
            fn(nx,j,k)=fn(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                i1=i-1
                wrkpx(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i1,j,k) ! all smoothing
                !        px(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i1,j,k)  ! x-smooth only
            enddo
        enddo
    enddo
    !
    !         step 2:   diffuse initial in y -dreiction
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=1,nx
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                oldpx(i,jp,k,m)+oldpx(i,jm,k,m) &
                -2.*oldpx(i,j,k,m))
            enddo
        enddo
    enddo
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            px(i,1,k,m)=wrkpx(i,1,k,m)
            px(i,ny,k,m)=wrkpx(i,ny,k,m)
        enddo
    enddo
    !
    !     calculate the difference  between x nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny-1
            j1=j+1
            do i=1,nx
                delta(i,j,k)=wrkpx(i,j1,k,m)-wrkpx(i,j,k,m)  !non-diffuse n+1 solution
                delta1(i,j,k)=px(i,j1,k,m)-px(i,j,k,m)   !diffused n+1 solution
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            delta(i,ny,k)=delta(i,ny-1,k)
            delta1(i,ny,k)=delta1(i,ny-1,k)
        enddo
    enddo

    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            do i=1,nx
                s=sign(1.,delta(i,j,k))
                fn(i,j,k)=s*amax1(0.,amin1(s*delta1(i,j-1,k), &
                chipx*abs(delta(i,j,k)),s*delta1(i,j+1,k)))
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            fn(i,1,k)=fn(i,2,k)
            fn(i,ny,k)=fn(i,ny-1,k)
        enddo
    enddo

    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            j1=j-1
            do i=1,nx
                !        wrkpx(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i,j1,k) ! for x-y-z smoothing
                px(i,j,k,m)=px(i,j,k,m)-fn(i,j,k)+fn(i,j1,k) ! for x-y smoothing
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine fcsmooth(px,oldpx,wrkpx,nx,ny,nz,ngrd,m,chipx, &
    tempx,tempy,tempz)
    !
    !     applies flux correction smoothing to field quantities
    !           oldpx is the qunatity at the start of the time step
    !           wrkpx is the unsmoothed qunatity
    !           px is the final product
    !     wrk arrays assumed to have dimension larger than nx,ny,nz
    !
    dimension px(nx,ny,nz,ngrd),oldpx(nx,ny,nz,ngrd), &
    wrkpx(nx,ny,nz,ngrd), &
    tempx(nx,ny,nz),tempy(nx,ny,nz),tempz(nx,ny,nz)
    !
    real,dimension(nx,ny,nz) :: deltax,delta1x,fnx, &
    deltay,delta1y,fny, &
    deltaz,delta1z,fnz
    !
    !      save initial unperturbed values in temp array
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                tempx(i,j,k)=wrkpx(i,j,k,m)
            enddo
        enddo
    enddo
    !
    !         step 1:   diffuse initial in all directions
    !$omp  parallel do
    do k=2,nz-1
        kp=k+1
        km=k-1
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=2,nx-1
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                oldpx(i+1,j,k,m)+oldpx(i-1,j,k,m) &
                +oldpx(i,jp,k,m)+oldpx(i,jm,k,m) &
                +oldpx(i,j,kp,m)+oldpx(i,j,km,m) &
                -6.*oldpx(i,j,k,m))
            enddo
        enddo
    enddo
    !
    !     set boundary points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            px(1,j,k,m)=wrkpx(1,j,k,m)
            px(nx,j,k,m)=wrkpx(nx,j,k,m)
            tempx(1,j,k)=wrkpx(1,j,k,m)
            tempx(nx,j,k)=wrkpx(nx,j,k,m)
        enddo
    enddo
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            px(i,j,1,m)=wrkpx(i,j,1,m)
            px(i,j,nz,m)=wrkpx(i,j,nz,m)
            tempx(i,j,1)=wrkpx(i,j,1,m)
            tempx(i,j,nz)=wrkpx(i,j,nz,m)
        enddo
    enddo
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            px(i,1,k,m)=wrkpx(i,1,k,m)
            px(i,ny,k,m)=wrkpx(i,ny,k,m)
            tempx(i,1,k)=wrkpx(i,1,k,m)
            tempx(i,ny,k)=wrkpx(i,ny,k,m)
        enddo
    enddo
    !
    !     calculate the difference  between x nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx-1
                i1=i+1
                deltax(i,j,k)=wrkpx(i1,j,k,m)-wrkpx(i,j,k,m)   !d-non-diffuse n+1 sol'n
                delta1x(i,j,k)=px(i1,j,k,m)-px(i,j,k,m)   !d1-diffused n+1 sol'n
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            deltax(nx,j,k)=deltax(nx-1,j,k)
            delta1x(nx,j,k)=delta1x(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                s=sign(1.,deltax(i,j,k))
                fnx(i,j,k)=s*amax1(0.,amin1(s*delta1x(i-1,j,k), &
                chipx*abs(deltax(i,j,k)),s*delta1x(i+1,j,k)))
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            fnx(1,j,k)=fnx(2,j,k)
            fnx(nx,j,k)=fnx(nx-1,j,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=2,nx-1
                i1=i-1
                tempx(i,j,k)=px(i,j,k,m)-fnx(i,j,k)+fnx(i1,j,k) ! all smoothing
                !        px(i,j,k,m)=px(i,j,k,m)-fnx(i)+fnx(i1)  ! x-smooth only
            enddo
        enddo
    enddo
    !
    !
    !
    !     calculate the difference  between y nearest grid points
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny-1
            j1=j+1
            do i=1,nx
                deltay(i,j,k)=wrkpx(i,j1,k,m)-wrkpx(i,j,k,m)  !non-diffuse n+1 solution
                delta1y(i,j,k)=px(i,j1,k,m)-px(i,j,k,m)   !diffused n+1 solution
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            deltay(i,ny,k)=deltay(i,ny-1,k)
            delta1y(i,ny,k)=delta1y(i,ny-1,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            do i=1,nx
                s=sign(1.,deltay(i,j,k))
                fny(i,j,k)=s*amax1(0.,amin1(s*delta1y(i,j-1,k), &
                chipx*abs(deltay(i,j,k)),s*delta1y(i,j+1,k)))
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            fny(i,1,k)=fny(i,2,k)
            fny(i,ny,k)=fny(i,ny-1,k)
        enddo
    enddo
    !
    !$omp  parallel do
    do k=1,nz
        do j=2,ny-1
            j1=j-1
            do i=1,nx
                tempx(i,j,k)=tempx(i,j,k)-fny(i,j,k)+fny(i,j1,k) ! for x-y-z smoothing
                !        px(i,j,k,m)=tempx(i,j,k)-fny(j)+fny(j1) ! for x-y smoothing
            enddo
        enddo
    enddo
    !
    !     calculatethe difference  between z nearest grid points
    !
    !$omp  parallel do
    do k=1,nz-1
        do j=1,ny
            do i=1,nx
                deltaz(i,j,k)=wrkpx(i,j,k+1,m)-wrkpx(i,j,k,m)  !non-diffuse n+1 soln
                delta1z(i,j,k)=px(i,j,k+1,m)-px(i,j,k,m)   !diffused n+1 solution
            enddo
        enddo
    enddo
    !
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            deltaz(i,j,nz)=deltaz(i,j,nz-1)
            delta1z(i,j,nz)=delta1z(i,j,nz-1)
        enddo
    enddo
    
    !$omp  parallel do
    do k=2,nz-1
        do j=1,ny
            do i=1,nx
                s=sign(1.,deltaz(i,j,k))
                fnz(i,j,k)=s*amax1(0.,amin1(s*delta1z(i,j,k-1), &
                chipx*abs(deltaz(i,j,k)),s*delta1z(i,j,k+1)))
            enddo
        enddo
    enddo
    
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            fnz(i,j,1)=fnz(i,j,2)
            fnz(i,j,nz)=fnz(i,j,nz-1)
        enddo
    enddo
    
    !$omp  parallel do
    do k=2,nz-1
        do j=1,ny
            do i=1,nx
                px(i,j,k,m)=tempx(i,j,k)-fnz(i,j,k)+fnz(i,j,k-1)
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine diffuse(qrho,qpres,qpx,qpy,qpz, &
    hrho,hpres,hpx,hpy,hpz, &
    orho,opres,opx,opy,opz,epres,bx,by,bz, &
    wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,wrkorho, &
    wrkopres,wrkopx,wrkopy,wrkopz, &
    wrkepres,wrkbx,wrkby,wrkbz,nx,ny,nz,ngrd,m, &
    difrho,difpxyz,diferg,delt,xspac,yspac,zspac)
    !
    !     apply straight artifical diffusion :
    !
    
    dimension qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpres(nx,ny,nz,ngrd), &
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd),wrkhrho(nx,ny,nz,ngrd), &
    wrkhpres(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd), &
    wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd), &
    wrkorho(nx,ny,nz,ngrd),wrkopres(nx,ny,nz,ngrd), &
    wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
    wrkopz(nx,ny,nz,ngrd),wrkbx(nx,ny,nz,ngrd), &
    wrkby(nx,ny,nz,ngrd),wrkbz(nx,ny,nz,ngrd), &
    wrkepres(nx,ny,nz,ngrd)
    !
    difb=diferg*delt/xspac
    dife=difrho*delt/xspac
    difp=difpxyz*delt/xspac
    difr=difrho*delt/xspac
    difo=0.25*difr
    !
    ! parallelizes loop rw, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kp=k+1
        do j=2,ny-1
            jm=j-1
            jp=j+1
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                !       species 1
                !
                qrho(i,j,k,m)=wrkqrho(i,j,k,m) +difr*( &
                +wrkqrho(ip,j,k,m)+wrkqrho(im,j,k,m) &
                +wrkqrho(i,jp,k,m)+wrkqrho(i,jm,k,m) &
                +wrkqrho(i,j,kp,m)+wrkqrho(i,j,km,m) &
                -6.*wrkqrho(i,j,k,m) )
    
                qpx(i,j,k,m)=wrkqpx(i,j,k,m) +difp*( &
                +wrkqpx(ip,j,k,m)+wrkqpx(im,j,k,m) &
                +wrkqpx(i,jp,k,m)+wrkqpx(i,jm,k,m) &
                +wrkqpx(i,j,kp,m)+wrkqpx(i,j,km,m) &
                -6.*wrkqpx(i,j,k,m) )
                !
                qpy(i,j,k,m)=wrkqpy(i,j,k,m) +difp*( &
                +wrkqpy(ip,j,k,m)+wrkqpy(im,j,k,m) &
                +wrkqpy(i,jp,k,m)+wrkqpy(i,jm,k,m) &
                +wrkqpy(i,j,kp,m)+wrkqpy(i,j,km,m) &
                -6.*wrkqpy(i,j,k,m) )
                !
                qpz(i,j,k,m)=wrkqpz(i,j,k,m) +difp*( &
                +wrkqpz(ip,j,k,m)+wrkqpz(im,j,k,m) &
                +wrkqpz(i,jp,k,m)+wrkqpz(i,jm,k,m) &
                +wrkqpz(i,j,kp,m)+wrkqpz(i,j,km,m) &
                -6.*wrkqpz(i,j,k,m) )
                !
                qpres(i,j,k,m)=wrkqpres(i,j,k,m) +dife*( &
                + wrkqpres(ip,j,k,m)+wrkqpres(im,j,k,m) &
                + wrkqpres(i,jp,k,m)+wrkqpres(i,jm,k,m) &
                + wrkqpres(i,j,kp,m)+wrkqpres(i,j,km,m) &
                -6.*wrkqpres(i,j,k,m) )
                !
                !       species 2
                !
                hrho(i,j,k,m)=wrkhrho(i,j,k,m) +difr*( &
                + wrkhrho(ip,j,k,m)+wrkhrho(im,j,k,m) &
                + wrkhrho(i,jp,k,m)+wrkhrho(i,jm,k,m) &
                + wrkhrho(i,j,kp,m)+wrkhrho(i,j,km,m) &
                -6.*wrkhrho(i,j,k,m) )
    
                hpx(i,j,k,m)=wrkhpx(i,j,k,m) +difp*( &
                +wrkhpx(ip,j,k,m)+wrkhpx(im,j,k,m) &
                +wrkhpx(i,jp,k,m)+wrkhpx(i,jm,k,m) &
                +wrkhpx(i,j,kp,m)+wrkhpx(i,j,km,m) &
                -6.*wrkhpx(i,j,k,m) )
                !
                hpy(i,j,k,m)=wrkhpy(i,j,k,m) +difp*( &
                +wrkhpy(ip,j,k,m)+wrkhpy(im,j,k,m) &
                +wrkhpy(i,jp,k,m)+wrkhpy(i,jm,k,m) &
                +wrkhpy(i,j,kp,m)+wrkhpy(i,j,km,m) &
                -6.*wrkhpy(i,j,k,m) )
                !
                hpz(i,j,k,m)=wrkhpz(i,j,k,m) +difp*( &
                +wrkhpz(ip,j,k,m)+wrkhpz(im,j,k,m) &
                +wrkhpz(i,jp,k,m)+wrkhpz(i,jm,k,m) &
                +wrkhpz(i,j,kp,m)+wrkhpz(i,j,km,m) &
                -6.*wrkhpz(i,j,k,m) )
                !
                hpres(i,j,k,m)=wrkhpres(i,j,k,m) +dife*( &
                + wrkhpres(ip,j,k,m)+wrkhpres(im,j,k,m) &
                + wrkhpres(i,jp,k,m)+wrkhpres(i,jm,k,m) &
                + wrkhpres(i,j,kp,m)+wrkhpres(i,j,km,m) &
                -6.*wrkhpres(i,j,k,m) )
                !
                !       species 3
                !
                orho(i,j,k,m)=wrkorho(i,j,k,m) +difo*( &
                + wrkorho(ip,j,k,m)+wrkorho(im,j,k,m) &
                + wrkorho(i,jp,k,m)+wrkorho(i,jm,k,m) &
                + wrkorho(i,j,kp,m)+wrkorho(i,j,km,m) &
                -6.*wrkorho(i,j,k,m) )
    
                opx(i,j,k,m)=wrkopx(i,j,k,m) +difo*( &
                +wrkopx(ip,j,k,m)+wrkopx(im,j,k,m) &
                +wrkopx(i,jp,k,m)+wrkopx(i,jm,k,m) &
                +wrkopx(i,j,kp,m)+wrkopx(i,j,km,m) &
                -6.*wrkopx(i,j,k,m) )
                !
                opy(i,j,k,m)=wrkopy(i,j,k,m) +difo*( &
                +wrkopy(ip,j,k,m)+wrkopy(im,j,k,m) &
                +wrkopy(i,jp,k,m)+wrkopy(i,jm,k,m) &
                +wrkopy(i,j,kp,m)+wrkopy(i,j,km,m) &
                -6.*wrkopy(i,j,k,m) )
                !
                opz(i,j,k,m)=wrkopz(i,j,k,m) +difo*( &
                +wrkopz(ip,j,k,m)+wrkopz(im,j,k,m) &
                +wrkopz(i,jp,k,m)+wrkopz(i,jm,k,m) &
                +wrkopz(i,j,kp,m)+wrkopz(i,j,km,m) &
                -6.*wrkopz(i,j,k,m) )
                !
                opres(i,j,k,m)=wrkopres(i,j,k,m) +dife*( &
                + wrkopres(ip,j,k,m)+wrkopres(im,j,k,m) &
                + wrkopres(i,jp,k,m)+wrkopres(i,jm,k,m) &
                + wrkopres(i,j,kp,m)+wrkopres(i,j,km,m) &
                -6.*wrkopres(i,j,k,m) )
                !
                epres(i,j,k,m)=wrkepres(i,j,k,m) +dife*( &
                + wrkepres(ip,j,k,m)+wrkepres(im,j,k,m) &
                + wrkepres(i,jp,k,m)+wrkepres(i,jm,k,m) &
                + wrkepres(i,j,kp,m)+wrkepres(i,j,km,m) &
                -6.*wrkepres(i,j,k,m) )
                !
                bx(i,j,k,m)=wrkbx(i,j,k,m)+difb*( &
                + wrkbx(ip,j,k,m)+wrkbx(im,j,k,m) &
                + wrkbx(i,jp,k,m)+wrkbx(i,jm,k,m) &
                + wrkbx(i,j,kp,m)+wrkbx(i,j,km,m) &
                -6.*wrkbx(i,j,k,m) )
                !
                by(i,j,k,m)=wrkby(i,j,k,m)+difb*( &
                + wrkby(ip,j,k,m)+wrkby(im,j,k,m) &
                + wrkby(i,jp,k,m)+wrkby(i,jm,k,m) &
                + wrkby(i,j,kp,m)+wrkby(i,j,km,m) &
                -6.*wrkby(i,j,k,m) )
                !
                bz(i,j,k,m)=wrkbz(i,j,k,m)+difb*( &
                + wrkbz(ip,j,k,m)+wrkbz(im,j,k,m) &
                + wrkbz(i,jp,k,m)+wrkbz(i,jm,k,m) &
                + wrkbz(i,j,kp,m)+wrkbz(i,j,km,m) &
                -6.*wrkbz(i,j,k,m) )
            enddo
        enddo
    enddo
    !
    !     reset time index so that desired quanties lie at nt1 rather than nt2
    !
    !     ntnew=nt2
    !     nt2=nt1
    !     nt1=ntnew
    !
    return
end

!
!     ********************************************
!
subroutine csmooth(curx,cury,curz,fldx,fldy,fldz, &
    nx,ny,nz)
    !
    !     force smoothing of surface currents to see if we can
    !             get better graphics
    !
    dimension curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    fldx(nx,ny,nz),fldy(nx,ny,nz),fldz(nx,ny,nz)
    !
    do k=1,nz
        do j=1,ny
            do i=1,nx
                fldx(i,j,k)=0.
                fldy(i,j,k)=0.
                fldz(i,j,k)=0.
            enddo
        enddo
    enddo
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                !
                fldx(i,j,k)= &
                ( curx(i+1,j,k)+curx(i-1,j,k) +curx(i,j+1,k)+curx(i,j-1,k) &
                +curx(i,j,k+1)+curx(i,j,k-1) +6.*curx(i,j,k))/12.
                fldy(i,j,k)= &
                ( cury(i+1,j,k)+cury(i-1,j,k) +cury(i,j+1,k)+cury(i,j-1,k) &
                +cury(i,j,k+1)+cury(i,j,k-1) +6.*cury(i,j,k))/12.
                fldz(i,j,k)= &
                ( curz(i+1,j,k)+curz(i-1,j,k) +curz(i,j+1,k)+curz(i,j-1,k) &
                +curz(i,j,k+1)+curz(i,j,k-1) +6.*curz(i,j,k))/12.
    
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine ssmooth(px,wrkpx,nx,ny,nz,ngrd,m,chipx)
    !
    !     applies straight diffusion
    !
    dimension px(nx,ny,nz,ngrd),wrkpx(nx,ny,nz,ngrd)
    !
    !     write(6,*)'ssmooth',nx,ny,nz,ngrd,m,chipx
    !
    !         step 1:   diffuse 3-d direction
    !
    !$omp  parallel do
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                px(i,j,k,m)=wrkpx(i,j,k,m)+chipx*( &
                wrkpx(i+1,j,k,m)+wrkpx(i-1,j,k,m) &
                + wrkpx(i,j+1,k,m)+wrkpx(i,j-1,k,m) &
                + wrkpx(i,j,k+1,m)+wrkpx(i,j,k-1,m) &
                -6.*wrkpx(i,j,k,m))
            enddo
        enddo
    enddo
    !
    !     write onto original array
    !$omp  parallel do
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                wrkpx(i,j,k,m)=px(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end

!
!    ******************************************************
!
subroutine divb_cor(bx,by,bz,px,py,pz,rho,ppres, &
    nx,ny,nz,ngrd,m,srho,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this subroutine corrects for possible non-divergence of b
    !         it integrates from the wind boundary
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd), &
    bz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd),ppres(nx,ny,nz,ngrd), &
    px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !     set limit of x integration depending on when you
    !         hit magnetopause
    !
    r_max=4.5*srho
    r_min=0.15*srho
    rad_lim=3.*rearth
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    dxt=2.*dx
    dyt=2.*dy
    dzt=2.*dz
    !
    do k=2,nz-1
        ak=grd_zmin(m)+dz*(k-1)
        do j=2,ny-1
            aj=grd_ymin(m)+dy*(j-1)
            i=2
    
            ai=grd_xmin(m)+dx*(i-1)
            ar=sqrt(ai**2+ak**2+aj**2)
            arho=rho(i,j,k,m)
            if((arho.gt.r_min).or.(ar.gt.rad_lim))then
                bx(i,j,k,m)=bx(i-1,j,k,m) &
                -(dx/dyt)*((by(i-1,j+1,k,m)+by(i,j+1,k,m))/2. &
                -(by(i-1,j-1,k,m)+by(i,j-1,k,m))/2. ) &
                -(dx/dzt)*((bz(i-1,j,k+1,m)+bz(i,j,k+1,m))/2. &
                -(bz(i-1,j,k-1,m)+bz(i,j,k-1,m))/2. )
            endif
            do i=3,nx-1
                ai=grd_xmin(m)+dx*(i-1)
                ar=sqrt(ai**2+ak**2+aj**2)
                arho=rho(i,j,k,m)
                if((arho.lt.r_min).or.(ar.lt.rad_lim)) exit
                bx(i,j,k,m)=bx(i-2,j,k,m) &
                -(dxt/dyt)*(by(i-1,j+1,k,m)-by(i-1,j-1,k,m)) &
                -(dxt/dzt)*(bz(i-1,j,k+1,m)-bz(i-1,j,k-1,m))
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine bndry_inner( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
    orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
    epres,bx,by,bz, &
    qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz, &
    opresxy,opresxz,opresyz, &
    nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
    ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
    mbndry,msrf,mmid,mzero, &
    erho,epress,alpha_e,ti_te,o_conc, &
    re_equiv,rearth,sbx_wind,spress, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    
    !
    !     this routine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    dimension &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    ! 
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real  parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    !
    !     set surface conditions around earth
    !      withinterior temperature as constant
    !
    !     write(6,*)'inside bndry inner with',mbndry
    !     write(6,*)  numsrf,nummid,numzero
    !
    aheight=rearth+2.
    !
    m=1
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    r_rot2=r_rot**2
    !
    !     set ionospheric boundary points - plasma quantities defined
    !
    do n=1,numsrf(m)
        !
        !        get coords of point
        !
        i=ijsrf(m,1,n)
        j=ijsrf(m,2,n)
        k=ijsrf(m,3,n)
        !
        !         set rotational speeds
        !
        ay=grd_ymin(m)+dy*(j-1)
        ax=grd_xmin(m)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        !        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        !        qden=parm_srf(1,n)
        qden=amax1(parm_srf(m,1,n),qrho(i,j,k,m))
        hden=parm_srf(m,2,n)
        oden=parm_srf(m,3,n)
    
        qrho(i,j,k,m)=qden
        hrho(i,j,k,m)=hden
        orho(i,j,k,m)=oden
        !
        qpresx(i,j,k,m)=parm_srf(m,4,n)
        hpresx(i,j,k,m)=parm_srf(m,5,n)
        opresx(i,j,k,m)=parm_srf(m,6,n)
        !
        qpresy(i,j,k,m)=parm_srf(m,4,n)
        hpresy(i,j,k,m)=parm_srf(m,5,n)
        opresy(i,j,k,m)=parm_srf(m,6,n)
        !
        qpresz(i,j,k,m)=parm_srf(m,4,n)
        hpresz(i,j,k,m)=parm_srf(m,5,n)
        opresz(i,j,k,m)=parm_srf(m,6,n)
        !
        epres(i,j,k,m)=parm_srf(m,7,n)
        !
        hpx(i,j,k,m)=hden*rvx
        hpy(i,j,k,m)=hden*rvy
        hpz(i,j,k,m)=hden*rvz
        opx(i,j,k,m)=oden*rvx
        opy(i,j,k,m)=oden*rvy
        opz(i,j,k,m)=oden*rvz
        qpx(i,j,k,m)=qden*rvx
        qpy(i,j,k,m)=qden*rvy
        qpz(i,j,k,m)=qden*rvz
        !
        !        iostropic conditions: set off axis pressures to zero
        !
        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
    enddo
    !
    !     set atmospheric boundary points - specify magnetic field profile
    !
    do n=1,nummid(m)
        !
        !        get coords of point
        !
        i=ijmid(m,1,n)
        j=ijmid(m,2,n)
        k=ijmid(m,3,n)
        !
        qrho(i,j,k,m)=parm_mid(m,1,n)
        hrho(i,j,k,m)=parm_mid(m,2,n)
        orho(i,j,k,m)=parm_mid(m,3,n)
        !
        qpresx(i,j,k,m)=parm_mid(m,4,n)
        hpresx(i,j,k,m)=parm_mid(m,5,n)
        opresx(i,j,k,m)=parm_mid(m,6,n)
        !
        qpresy(i,j,k,m)=parm_mid(m,4,n)
        hpresy(i,j,k,m)=parm_mid(m,5,n)
        opresy(i,j,k,m)=parm_mid(m,6,n)
        !
        qpresz(i,j,k,m)=parm_mid(m,4,n)
        hpresz(i,j,k,m)=parm_mid(m,5,n)
        opresz(i,j,k,m)=parm_mid(m,6,n)
        !
        epres(i,j,k,m)=parm_mid(m,7,n)
        !
        ay=grd_ymin(m)+dy*(j-1)
        ax=grd_xmin(m)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        !        bx(i,j,k,m)=sbx_wind
        !        by(i,j,k,m)=0.
        !        bz(i,j,k,m)=0.
        !
    enddo
    !
    !
    !      reset interior points
    !
    do n=1,numzero(m)
        i=ijzero(m,1,n)
        j=ijzero(m,2,n)
        k=ijzero(m,3,n)
        !
        qrho(i,j,k,m)=parm_zero(m,1,n)
        hrho(i,j,k,m)=parm_zero(m,2,n)
        orho(i,j,k,m)=parm_zero(m,3,n)
        !
        qpresx(i,j,k,m)=parm_zero(m,4,n)
        hpresx(i,j,k,m)=parm_zero(m,5,n)
        opresx(i,j,k,m)=parm_zero(m,6,n)
        !
        qpresy(i,j,k,m)=parm_zero(m,4,n)
        hpresy(i,j,k,m)=parm_zero(m,5,n)
        opresy(i,j,k,m)=parm_zero(m,6,n)
        !
        qpresz(i,j,k,m)=parm_zero(m,4,n)
        hpresz(i,j,k,m)=parm_zero(m,5,n)
        opresz(i,j,k,m)=parm_zero(m,6,n)
        !
        epres(i,j,k,m)=parm_zero(m,7,n)
        !
        ay=grd_ymin(m)+dy*(j-1)
        ax=grd_xmin(m)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        bx(i,j,k,m)=0.
        by(i,j,k,m)=0.
        bz(i,j,k,m)=0.
        !
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine bndry_outer( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz,opx,opy,opz, &
    epres, &
    qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz, &
    opresxy,opresxz,opresyz, &
    rmassq,rmassh,rmasso, &
    bx,by,bz,bx0,by0,bz0,nx,ny,nz,ngrd, &
    srho,rho_frac,o_conc,spress,spx,spy,spz, &
    sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)
    !
    !     this routine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    !
    dimension &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    ! 
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
    !
    !     do x boundary conditions : set oxygen in solar wind
    !       at 2.5%
    !
    
    m=ngrd
    ofrac=rho_frac*o_conc
    frac_h=rmassh/rmassq
    frac_o=rmasso/rmassq
    !
    i=1
    ii=2
    do k=1,nz
        do j=1,ny
            !
            !       wind boundary conditions at i=1
            !
            qrho(i,j,k,m)=srho
            qpx(i,j,k,m)=spx
            qpy(i,j,k,m)=spy
            qpz(i,j,k,m)=spz
            qpresx(i,j,k,m)=0.5*spress
            qpresy(i,j,k,m)=qpresx(i,j,k,m)
            qpresz(i,j,k,m)=qpresx(i,j,k,m)
            !
            hrho(i,j,k,m)=srho*rho_frac*frac_h
            hpx(i,j,k,m)=spx*rho_frac*frac_h
            hpy(i,j,k,m)=spy*rho_frac*frac_h
            hpz(i,j,k,m)=spz*rho_frac*frac_h
            hpresx(i,j,k,m)=0.5*spress*rho_frac
            hpresy(i,j,k,m)=hpresx(i,j,k,m)
            hpresz(i,j,k,m)=hpresx(i,j,k,m)
            !
            orho(i,j,k,m)=srho*ofrac*frac_o
            opx(i,j,k,m)=spx*ofrac*frac_o
            opy(i,j,k,m)=spy*ofrac*frac_o
            opz(i,j,k,m)=spz*ofrac*frac_o
            opresx(i,j,k,m)=0.5*spress*ofrac
            opresy(i,j,k,m)=opresx(i,j,k,m)
            opresz(i,j,k,m)=opresx(i,j,k,m)
            !
            epres(i,j,k,m)=qpresx(i,j,k,m)/ti_te
            !
            !        iostropic conditions: set off axis pressures to zero
            !
            qpresxy(i,j,k,m)=0.
            qpresxz(i,j,k,m)=0.
            qpresyz(i,j,k,m)=0.
            hpresxy(i,j,k,m)=0.
            hpresxz(i,j,k,m)=0.
            hpresyz(i,j,k,m)=0.
            opresxy(i,j,k,m)=0.
            opresxz(i,j,k,m)=0.
            opresyz(i,j,k,m)=0.
            !
            !
            !       set bx of wind to zero, by to the wind values
            !          and subtract any geomagnetic field
            !
            bx(i,j,k,m)=sbx_wind-bx0(i,j,k,m)
            by(i,j,k,m)=sby_wind-by0(i,j,k,m)
            bz(i,j,k,m)=sbz_wind-bz0(i,j,k,m)
            !
            !       bx(i,j,k,m)=bx(2,j,k,m)
            !       by(i,j,k,m)=by(2,j,k,m)
            !       bz(i,j,k,m)=bz(2,j,k,m)
            !
        enddo
    enddo
    !
    !
    
    !     do x boundary conditions - back wall
    !
    nx1=nx-1
    do k=1,nz
        do j=1,ny
            !
            qrho(nx,j,k,m)=qrho(nx1,j,k,m)
            qpx(nx,j,k,m)=abs(qpx(nx1,j,k,m))
            qpy(nx,j,k,m)=qpy(nx1,j,k,m)
            qpz(nx,j,k,m)=qpz(nx1,j,k,m)
            qpresx(nx,j,k,m)=qpresx(nx1,j,k,m)
            qpresy(nx,j,k,m)=qpresy(nx1,j,k,m)
            qpresz(nx,j,k,m)=qpresz(nx1,j,k,m)
            !
            hrho(nx,j,k,m)=hrho(nx1,j,k,m)
            hpx(nx,j,k,m)=abs(hpx(nx1,j,k,m))
            hpy(nx,j,k,m)=hpy(nx1,j,k,m)
            hpz(nx,j,k,m)=hpz(nx1,j,k,m)
            hpresx(nx,j,k,m)=hpresx(nx1,j,k,m)
            hpresy(nx,j,k,m)=hpresy(nx1,j,k,m)
            hpresz(nx,j,k,m)=hpresz(nx1,j,k,m)
            !
            orho(nx,j,k,m)=orho(nx1,j,k,m)
            opx(nx,j,k,m)=abs(opx(nx1,j,k,m))
            opy(nx,j,k,m)=opy(nx1,j,k,m)
            opz(nx,j,k,m)=opz(nx1,j,k,m)
            opresx(nx,j,k,m)=opresx(nx1,j,k,m)

            opresy(nx,j,k,m)=opresy(nx1,j,k,m)
            opresz(nx,j,k,m)=opresz(nx1,j,k,m)
            !
            epres(nx,j,k,m)=epres(nx1,j,k,m)
            !
            !       continuous boundary condition
            !
            qpresxy(nx,j,k,m)=qpresxy(nx1,j,k,m)
            qpresxz(nx,j,k,m)=qpresxz(nx1,j,k,m)
            qpresyz(nx,j,k,m)=qpresyz(nx1,j,k,m)
            hpresxy(nx,j,k,m)=hpresxy(nx1,j,k,m)
            hpresxz(nx,j,k,m)=hpresxz(nx1,j,k,m)
            hpresyz(nx,j,k,m)=hpresyz(nx1,j,k,m)
            opresxy(nx,j,k,m)=opresxy(nx1,j,k,m)
            opresxz(nx,j,k,m)=opresxz(nx1,j,k,m)
            opresyz(nx,j,k,m)=opresyz(nx1,j,k,m)
            !
            bx(nx,j,k,m)=bx(nx1,j,k,m)
            by(nx,j,k,m)=by(nx1,j,k,m)
            bz(nx,j,k,m)=bz(nx1,j,k,m)
        enddo
    enddo
    !
    !     do y boundary conditions
    !
    ny1=ny-1
    do k=1,nz
        do i=2,nx
            i1=i
            !
            qrho(i,1,k,m)=qrho(i1,2,k,m)
            qpresx(i,1,k,m)=qpresx(i1,2,k,m)
            qpresy(i,1,k,m)=qpresy(i1,2,k,m)
            qpresz(i,1,k,m)=qpresz(i1,2,k,m)

            qpx(i,1,k,m)=qpx(i1,2,k,m)
            qpy(i,1,k,m)=-abs(qpy(i1,2,k,m))
            qpz(i,1,k,m)=qpz(i1,2,k,m)
            !
            hrho(i,1,k,m)=hrho(i1,2,k,m)
            hpresx(i,1,k,m)=hpresx(i1,2,k,m)
            hpresy(i,1,k,m)=hpresy(i1,2,k,m)
            hpresz(i,1,k,m)=hpresz(i1,2,k,m)

            hpx(i,1,k,m)=hpx(i1,2,k,m)
            hpy(i,1,k,m)=-abs(hpy(i1,2,k,m))
            hpz(i,1,k,m)=hpz(i1,2,k,m)
            !
            orho(i,1,k,m)=orho(i1,2,k,m)
            opresx(i,1,k,m)=opresx(i1,2,k,m)
            opresy(i,1,k,m)=opresy(i1,2,k,m)
            opresz(i,1,k,m)=opresz(i1,2,k,m)

            opx(i,1,k,m)=opx(i1,2,k,m)
            opy(i,1,k,m)=-abs(opy(i1,2,k,m))
            opz(i,1,k,m)=opz(i1,2,k,m)
            !
            epres(i,1,k,m)=epres(i1,2,k,m)
            !
            !       continuous boundary condition
            !
            qpresxy(i,1,k,m)=qpresxy(i1,2,k,m)
            qpresxz(i,1,k,m)=qpresxz(i1,2,k,m)
            qpresyz(i,1,k,m)=qpresyz(i1,2,k,m)
            hpresxy(i,1,k,m)=hpresxy(i1,2,k,m)
            hpresxz(i,1,k,m)=hpresxz(i1,2,k,m)
            hpresyz(i,1,k,m)=hpresyz(i1,2,k,m)
            opresxy(i,1,k,m)=opresxy(i1,2,k,m)
            opresxz(i,1,k,m)=opresxz(i1,2,k,m)
            opresyz(i,1,k,m)=opresyz(i1,2,k,m)
            !
            by(i,1,k,m)=by(i1,2,k,m)
            !    +             +(by0(i1,2,k,m)-by0(i,1,k,m))
            bx(i,1,k,m)=bx(i1,2,k,m)
            !    +              +(bx0(i1,2,k,m)-bx0(i,1,k,m))
            bz(i,1,k,m)=bz(i1,2,k,m)
            !    +              +(bz0(i1,2,k,m)-bz0(i,1,k,m))
            !
            qrho(i,ny,k,m)=qrho(i1,ny1,k,m)
            qpresx(i,ny,k,m)=qpresx(i1,ny1,k,m)
            qpresy(i,ny,k,m)=qpresy(i1,ny1,k,m)
            qpresz(i,ny,k,m)=qpresz(i1,ny1,k,m)

            qpx(i,ny,k,m)=qpx(i1,ny1,k,m)
            qpy(i,ny,k,m)=abs(qpy(i1,ny1,k,m))
            qpz(i,ny,k,m)=qpz(i1,ny1,k,m)
            !
            hrho(i,ny,k,m)=hrho(i1,ny1,k,m)
            hpresx(i,ny,k,m)=hpresx(i1,ny1,k,m)
            hpresy(i,ny,k,m)=hpresy(i1,ny1,k,m)
            hpresz(i,ny,k,m)=hpresz(i1,ny1,k,m)
            hpx(i,ny,k,m)=hpx(i1,ny1,k,m)
            hpy(i,ny,k,m)=abs(hpy(i1,ny1,k,m))
            hpz(i,ny,k,m)=hpz(i1,ny1,k,m)
            !
            orho(i,ny,k,m)=orho(i1,ny1,k,m)
            opresx(i,ny,k,m)=opresx(i1,ny1,k,m)
            opresy(i,ny,k,m)=opresy(i1,ny1,k,m)
            opresz(i,ny,k,m)=opresz(i1,ny1,k,m)

            opx(i,ny,k,m)=opx(i1,ny1,k,m)
            opy(i,ny,k,m)=abs(opy(i1,ny1,k,m))
            opz(i,ny,k,m)=opz(i1,ny1,k,m)
            !
            epres(i,ny,k,m)=epres(i1,ny1,k,m)
            !
            !       continuous boundary condition
            !
            !
            qpresxy(i,ny,k,m)=qpresxy(i1,ny1,k,m)
            qpresxz(i,ny,k,m)=qpresxz(i1,ny1,k,m)
            qpresyz(i,ny,k,m)=qpresyz(i1,ny1,k,m)
            hpresxy(i,ny,k,m)=hpresxy(i1,ny1,k,m)
            hpresxz(i,ny,k,m)=hpresxz(i1,ny1,k,m)
            hpresyz(i,ny,k,m)=hpresyz(i1,ny1,k,m)
            opresxy(i,ny,k,m)=opresxy(i1,ny1,k,m)
            opresxz(i,ny,k,m)=opresxz(i1,ny1,k,m)
            opresyz(i,ny,k,m)=opresyz(i1,ny1,k,m)

            !
            by(i,ny,k,m)=by(i1,ny1,k,m)
            !    +              +(by0(i1,ny1,k,m)-by0(i,ny,k,m))
            bx(i,ny,k,m)=bx(i1,ny1,k,m)
            !    +               +(bx0(i1,ny1,k,m)-bx0(i,ny,k,m))
            bz(i,ny,k,m)=bz(i1,ny1,k,m)
            !    +              +(bz0(i1,ny1,k,m)-bz0(i,ny,k,m))
        enddo
    enddo
    !
    !     do z boundary conditions
    !
    nz1=nz-1
    do j=1,ny
        do i=2,nx
            i1=i
            !
            qrho(i,j,nz,m)=qrho(i1,j,nz1,m)
            qpresx(i,j,nz,m)=qpresx(i1,j,nz1,m)
            qpresy(i,j,nz,m)=qpresy(i1,j,nz1,m)
            qpresz(i,j,nz,m)=qpresz(i1,j,nz1,m)

            qpx(i,j,nz,m)=qpx(i1,j,nz1,m)
            qpy(i,j,nz,m)=qpy(i1,j,nz1,m)
            qpz(i,j,nz,m)=abs(qpz(i1,j,nz1,m))
            !
            hrho(i,j,nz,m)=hrho(i1,j,nz1,m)
            hpresx(i,j,nz,m)=hpresx(i1,j,nz1,m)
            hpresy(i,j,nz,m)=hpresy(i1,j,nz1,m)
            hpresz(i,j,nz,m)=hpresz(i1,j,nz1,m)

            hpx(i,j,nz,m)=hpx(i1,j,nz1,m)
            hpy(i,j,nz,m)=hpy(i1,j,nz1,m)
            hpz(i,j,nz,m)=abs(hpz(i1,j,nz1,m))
            !
            orho(i,j,nz,m)=orho(i1,j,nz1,m)
            opresx(i,j,nz,m)=opresx(i1,j,nz1,m)
            opresy(i,j,nz,m)=opresy(i1,j,nz1,m)
            opresz(i,j,nz,m)=opresz(i1,j,nz1,m)

            opx(i,j,nz,m)=opx(i1,j,nz1,m)
            opy(i,j,nz,m)=opy(i1,j,nz1,m)
            opz(i,j,nz,m)=abs(opz(i1,j,nz1,m))
            !
            epres(i,j,nz,m)=epres(i1,j,nz1,m)
            !
            !       off diagnoal elements
            !
            qpresxy(i,j,nz,m)=qpresxy(i1,j,nz1,m)
            qpresxz(i,j,nz,m)=qpresxz(i1,j,nz1,m)
            qpresyz(i,j,nz,m)=qpresyz(i1,j,nz1,m)
            hpresxy(i,j,nz,m)=hpresxy(i1,j,nz1,m)
            hpresxz(i,j,nz,m)=hpresxz(i1,j,nz1,m)
            hpresyz(i,j,nz,m)=hpresyz(i1,j,nz1,m)
            opresxy(i,j,nz,m)=opresxy(i1,j,nz1,m)
            opresxz(i,j,nz,m)=opresxz(i1,j,nz1,m)
            opresyz(i,j,nz,m)=opresyz(i1,j,nz1,m)
            !
            !
            bz(i,j,nz,m)=bz(i1,j,nz1,m)
            !    +            +bz0(i1,j,nz1,m)-bz0(i,j,nz,m)
            bx(i,j,nz,m)=bx(i1,j,nz1,m)
            !    +            +bx0(i1,j,nz1,m)-bx0(i,j,nz,m)
            by(i,j,nz,m)=by(i1,j,nz1,m)
            !    +            +by0(i1,j,nz1,m)-by0(i,j,nz,m)
        enddo
    enddo
    !
    !       symmetry boundary conditions at k=1 and k=2
    !
    do j=1,ny
        do i=2,nx
            !
            qrho(i,j,1,m)=qrho(i,j,2,m)
            qpx(i,j,1,m)=qpx(i,j,2,m)
            qpy(i,j,1,m)=qpy(i,j,2,m)
            qpz(i,j,1,m)=-abs(qpz(i,j,2,m))
            qpresx(i,j,1,m)=qpresx(i,j,2,m)
            qpresy(i,j,1,m)=qpresy(i,j,2,m)
            qpresz(i,j,1,m)=qpresz(i,j,2,m)
            !
            hrho(i,j,1,m)=hrho(i,j,2,m)
            hpx(i,j,1,m)=hpx(i,j,2,m)
            hpy(i,j,1,m)=hpy(i,j,2,m)
            hpz(i,j,1,m)=-abs(hpz(i,j,2,m))
            hpresx(i,j,1,m)=hpresx(i,j,2,m)
            hpresy(i,j,1,m)=hpresy(i,j,2,m)
            hpresz(i,j,1,m)=hpresz(i,j,2,m)
            !
            orho(i,j,1,m)=orho(i,j,2,m)
            opx(i,j,1,m)=opx(i,j,2,m)
            opy(i,j,1,m)=opy(i,j,2,m)
            opz(i,j,1,m)=-abs(opz(i,j,2,m))
            opresx(i,j,1,m)=opresx(i,j,2,m)
            opresy(i,j,1,m)=opresy(i,j,2,m)
            opresz(i,j,1,m)=opresz(i,j,2,m)
            !
            epres(i,j,1,m)=epres(i,j,2,m)
            !
            !       off diagonal elements
            !
            qpresxy(i,j,1,m)=qpresxy(i,j,2,m)
            qpresxz(i,j,1,m)=qpresxz(i,j,2,m)
            qpresyz(i,j,1,m)=qpresyz(i,j,2,m)
            hpresxy(i,j,1,m)=hpresxy(i,j,2,m)
            hpresxz(i,j,1,m)=hpresxz(i,j,2,m)
            hpresyz(i,j,1,m)=hpresyz(i,j,2,m)
            opresxy(i,j,1,m)=opresxy(i,j,2,m)
            opresxz(i,j,1,m)=opresxz(i,j,2,m)
            opresyz(i,j,1,m)=opresyz(i,j,2,m)
            !
            !       set bx of wind to zero, by to the wind values
            !          and subtract any geomagnetic field
            !       force bz to be positive to attain equilibrium
            !
            bx(i,j,1,m)=bx(i,j,2,m)
            by(i,j,1,m)=by(i,j,2,m)
            bz(i,j,1,m)=bz(i,j,2,m)
        enddo
    enddo
    !
    !     do corner lines at left-right back wall
    !
    do k=2,nz-1
        qrho(nx,1,k,m)=(qrho(nx1,1,k,m)+qrho(nx,2,k,m))/2.
        qrho(nx,ny,k,m)=(qrho(nx1,ny,k,m)+qrho(nx,ny1,k,m))/2.
        hrho(nx,1,k,m)=(hrho(nx1,1,k,m)+hrho(nx,2,k,m))/2.
        hrho(nx,ny,k,m)=(hrho(nx1,ny,k,m)+hrho(nx,ny1,k,m))/2.
        orho(nx,1,k,m)=(orho(nx1,1,k,m)+orho(nx,2,k,m))/2.
        orho(nx,ny,k,m)=(orho(nx1,ny,k,m)+orho(nx,ny1,k,m))/2.
        !
        qpresx(nx,1,k,m)=(qpresx(nx1,1,k,m) &
        +qpresx(nx,2,k,m))/2.
        qpresx(nx,ny,k,m)=(qpresx(nx1,ny,k,m) &
        +qpresx(nx,ny1,k,m))/2.
        qpresy(nx,1,k,m)=(qpresy(nx1,1,k,m) &
            +qpresy(nx,2,k,m))/2.
        qpresy(nx,ny,k,m)=(qpresy(nx1,ny,k,m) &
            +qpresy(nx,ny1,k,m))/2.
        qpresz(nx,1,k,m)=(qpresz(nx1,1,k,m) &
            +qpresz(nx,2,k,m))/2.
        qpresz(nx,ny,k,m)=(qpresz(nx1,ny,k,m) &
            +qpresz(nx,ny1,k,m))/2.
        !
        hpresx(nx,1,k,m)=(hpresx(nx1,1,k,m) &
        +hpresx(nx,2,k,m))/2.
        hpresx(nx,ny,k,m)=(hpresx(nx1,ny,k,m) &
        +hpresx(nx,ny1,k,m))/2.
        hpresy(nx,1,k,m)=(hpresy(nx1,1,k,m) &
            +hpresy(nx,2,k,m))/2.
        hpresy(nx,ny,k,m)=(hpresy(nx1,ny,k,m) &
            +hpresy(nx,ny1,k,m))/2.
        hpresz(nx,1,k,m)=(hpresz(nx1,1,k,m) &
            +hpresz(nx,2,k,m))/2.
        hpresz(nx,ny,k,m)=(hpresz(nx1,ny,k,m) &
            +hpresz(nx,ny1,k,m))/2.
        !
        opresx(nx,1,k,m)=(opresx(nx1,1,k,m) &
        +opresx(nx,2,k,m))/2.
        opresx(nx,ny,k,m)=(opresx(nx1,ny,k,m) &
        +opresx(nx,ny1,k,m))/2.
        opresy(nx,1,k,m)=(opresy(nx1,1,k,m) &
            +opresy(nx,2,k,m))/2.
        opresy(nx,ny,k,m)=(opresy(nx1,ny,k,m) &
            +opresy(nx,ny1,k,m))/2.
        opresz(nx,1,k,m)=(opresz(nx1,1,k,m) &
            +opresz(nx,2,k,m))/2.
        opresz(nx,ny,k,m)=(opresz(nx1,ny,k,m) &
            +opresz(nx,ny1,k,m))/2.
        !
        epres(nx,1,k,m)=(epres(nx1,1,k,m) &
        +epres(nx,2,k,m))/2.
        epres(nx,ny,k,m)=(epres(nx1,ny,k,m) &
        +epres(nx,ny1,k,m))/2.
        !
        !       off diagonal elements
        !
        qpresxy(nx,1,k,m)=(qpresxy(nx1,1,k,m) &
            +qpresxy(nx,2,k,m))/2.
        qpresxy(nx,ny,k,m)=(qpresxy(nx1,ny,k,m) &
            +qpresxy(nx,ny1,k,m))/2.
        qpresxz(nx,1,k,m)=(qpresxz(nx1,1,k,m) &
            +qpresxz(nx,2,k,m))/2.
        qpresxz(nx,ny,k,m)=(qpresxz(nx1,ny,k,m) &
            +qpresxz(nx,ny1,k,m))/2.
        qpresyz(nx,1,k,m)=(qpresyz(nx1,1,k,m) &
            +qpresyz(nx,2,k,m))/2.
        qpresyz(nx,ny,k,m)=(qpresyz(nx1,ny,k,m) &
            +qpresyz(nx,ny1,k,m))/2.
        hpresxy(nx,1,k,m)=(hpresxy(nx1,1,k,m) &
            +hpresxy(nx,2,k,m))/2.
        hpresxy(nx,ny,k,m)=(hpresxy(nx1,ny,k,m) &
            +hpresxy(nx,ny1,k,m))/2.
        hpresxz(nx,1,k,m)=(hpresxz(nx1,1,k,m) &
            +hpresxz(nx,2,k,m))/2.
        hpresxz(nx,ny,k,m)=(hpresxz(nx1,ny,k,m) &
            +hpresxz(nx,ny1,k,m))/2.
        hpresyz(nx,1,k,m)=(hpresyz(nx1,1,k,m) &
            +hpresyz(nx,2,k,m))/2.
        hpresyz(nx,ny,k,m)=(hpresyz(nx1,ny,k,m) &
            +hpresyz(nx,ny1,k,m))/2.
        opresxy(nx,1,k,m)=(opresxy(nx1,1,k,m) &
            +opresxy(nx,2,k,m))/2.
        opresxy(nx,ny,k,m)=(opresxy(nx1,ny,k,m) &
            +opresxy(nx,ny1,k,m))/2.
        opresxz(nx,1,k,m)=(opresxz(nx1,1,k,m) &
            +opresxz(nx,2,k,m))/2.
        opresxz(nx,ny,k,m)=(opresxz(nx1,ny,k,m) &
            +opresxz(nx,ny1,k,m))/2.
        opresyz(nx,1,k,m)=(opresyz(nx1,1,k,m) &
            +opresyz(nx,2,k,m))/2.
        opresyz(nx,ny,k,m)=(opresyz(nx1,ny,k,m) &
            +opresyz(nx,ny1,k,m))/2.
        !
        qpx(nx,1,k,m)=(qpx(nx1,1,k,m)+qpx(nx,2,k,m))/2.
        qpx(nx,ny,k,m)=(qpx(nx1,ny,k,m)+qpx(nx,ny1,k,m))/2.
        !
        qpy(nx,1,k,m)=(qpy(nx1,1,k,m)+qpy(nx,2,k,m))/2.
        qpy(nx,ny,k,m)=(qpy(nx1,ny,k,m)+qpy(nx,ny1,k,m))/2.
        !
        qpz(nx,1,k,m)=(qpz(nx1,1,k,m)+qpz(nx,2,k,m))/2.
        qpz(nx,ny,k,m)=(qpz(nx1,ny,k,m)+qpz(nx,ny1,k,m))/2.
        !
        hpx(nx,1,k,m)=(hpx(nx1,1,k,m)+hpx(nx,2,k,m))/2.
        hpx(nx,ny,k,m)=(hpx(nx1,ny,k,m)+hpx(nx,ny1,k,m))/2.
        !
        hpy(nx,1,k,m)=(hpy(nx1,1,k,m)+hpy(nx,2,k,m))/2.
        hpy(nx,ny,k,m)=(hpy(nx1,ny,k,m)+hpy(nx,ny1,k,m))/2.
        !
        hpz(nx,1,k,m)=(hpz(nx1,1,k,m)+hpz(nx,2,k,m))/2.
        hpz(nx,ny,k,m)=(hpz(nx1,ny,k,m)+hpz(nx,ny1,k,m))/2.
        !
        opx(nx,1,k,m)=(opx(nx1,1,k,m)+opx(nx,2,k,m))/2.
        opx(nx,ny,k,m)=(opx(nx1,ny,k,m)+opx(nx,ny1,k,m))/2.
        !
        opy(nx,1,k,m)=(opy(nx1,1,k,m)+opy(nx,2,k,m))/2.
        opy(nx,ny,k,m)=(opy(nx1,ny,k,m)+opy(nx,ny1,k,m))/2.
        !
        opz(nx,1,k,m)=(opz(nx1,1,k,m)+opz(nx,2,k,m))/2.
        opz(nx,ny,k,m)=(opz(nx1,ny,k,m)+opz(nx,ny1,k,m))/2.
        !
        bx(nx,1,k,m)=bx(nx1,1,k,m)
        bx(nx,ny,k,m)=bx(nx1,ny,k,m)
        !
        by(nx,1,k,m)=(by(nx1,1,k,m) &
        +by(nx,2,k,m))/2.
        by(nx,ny,k,m)=(by(nx1,ny,k,m) &
        +by(nx,ny1,k,m))/2.
        !
        bz(nx,1,k,m)=(bz(nx1,1,k,m) &
        +bz(nx,2,k,m))/2.
        bz(nx,ny,k,m)=(bz(nx1,ny,k,m) &
        +bz(nx,ny1,k,m))/2.
        !
    enddo
    !
    !     do corner lines at up and bottom back wall
    !
    do j=2,ny-1
        qrho(nx,j,1,m)=qrho(nx,j,2,m)
        qrho(nx,j,nz,m)=qrho(nx1,j,nz,m)
        hrho(nx,j,1,m)=hrho(nx,j,2,m)
        hrho(nx,j,nz,m)=hrho(nx1,j,nz,m)
        orho(nx,j,1,m)=orho(nx,j,2,m)
        orho(nx,j,nz,m)=orho(nx1,j,nz,m)
        !
        qpresx(nx,j,1,m)=qpresx(nx,j,2,m)
        qpresx(nx,j,nz,m)=qpresx(nx1,j,nz,m)
        qpresy(nx,j,1,m)=qpresy(nx,j,2,m)
        qpresy(nx,j,nz,m)=qpresy(nx1,j,nz,m)
        qpresz(nx,j,1,m)=qpresz(nx,j,2,m)
        qpresz(nx,j,nz,m)=qpresz(nx1,j,nz,m)
        !
        hpresx(nx,j,1,m)=hpresx(nx,j,2,m)
        hpresx(nx,j,nz,m)=hpresx(nx1,j,nz,m)
        hpresy(nx,j,1,m)=hpresy(nx,j,2,m)
        hpresy(nx,j,nz,m)=hpresy(nx1,j,nz,m)
        hpresz(nx,j,1,m)=hpresz(nx,j,2,m)
        hpresz(nx,j,nz,m)=hpresz(nx1,j,nz,m)
        !
        opresx(nx,j,1,m)=opresx(nx,j,2,m)
        opresx(nx,j,nz,m)=opresx(nx1,j,nz,m)
        opresy(nx,j,1,m)=opresy(nx,j,2,m)
        opresy(nx,j,nz,m)=opresy(nx1,j,nz,m)
        opresz(nx,j,1,m)=opresz(nx,j,2,m)
        opresz(nx,j,nz,m)=opresz(nx1,j,nz,m)
        !
        epres(nx,j,1,m)=epres(nx,j,2,m)
        !
        !       off diagnonal elements
        !
        qpresxy(nx,j,1,m)=qpresxy(nx,j,2,m)
        qpresxy(nx,j,nz,m)=qpresxy(nx1,j,nz,m)
        qpresxz(nx,j,1,m)=qpresxz(nx,j,2,m)
        qpresxz(nx,j,nz,m)=qpresxz(nx1,j,nz,m)
        qpresyz(nx,j,1,m)=qpresyz(nx,j,2,m)
        qpresyz(nx,j,nz,m)=qpresyz(nx1,j,nz,m)
        hpresxy(nx,j,1,m)=hpresxy(nx,j,2,m)
        hpresxy(nx,j,nz,m)=hpresxy(nx1,j,nz,m)
        hpresxz(nx,j,1,m)=hpresxz(nx,j,2,m)
        hpresxz(nx,j,nz,m)=hpresxz(nx1,j,nz,m)
        hpresyz(nx,j,1,m)=hpresyz(nx,j,2,m)
        hpresyz(nx,j,nz,m)=hpresyz(nx1,j,nz,m)
        opresxy(nx,j,1,m)=opresxy(nx,j,2,m)
        opresxy(nx,j,nz,m)=opresxy(nx1,j,nz,m)
        opresxz(nx,j,1,m)=opresxz(nx,j,2,m)
        opresxz(nx,j,nz,m)=opresxz(nx1,j,nz,m)
        opresyz(nx,j,1,m)=opresyz(nx,j,2,m)
        opresyz(nx,j,nz,m)=opresyz(nx1,j,nz,m)
        !
        qpx(nx,j,1,m)=qpx(nx,j,2,m)
        qpx(nx,j,nz,m)=qpx(nx1,j,nz,m)
        !
        qpy(nx,j,1,m)=qpy(nx,j,2,m)
        qpy(nx,j,nz,m)=qpy(nx1,j,nz,m)
        !
        qpz(nx,j,1,m)=qpz(nx,j,2,m)
        qpz(nx,j,nz,m)=qpz(nx1,j,nz,m)
        !
        hpx(nx,j,1,m)=hpx(nx,j,2,m)
        hpx(nx,j,nz,m)=hpx(nx1,j,nz,m)
        !
        hpy(nx,j,1,m)=hpy(nx,j,2,m)
        hpy(nx,j,nz,m)=hpy(nx1,j,nz,m)
        !
        hpz(nx,j,1,m)=hpz(nx,j,2,m)
        hpz(nx,j,nz,m)=hpz(nx1,j,nz,m)
        !
        opx(nx,j,1,m)=opx(nx,j,2,m)
        opx(nx,j,nz,m)=opx(nx1,j,nz,m)
        !
        opy(nx,j,1,m)=opy(nx,j,2,m)
        opy(nx,j,nz,m)=opy(nx1,j,nz,m)
        !
        opz(nx,j,1,m)=opz(nx,j,2,m)
        opz(nx,j,nz,m)=opz(nx1,j,nz,m)
        !
        bx(nx,j,1,m)=bx(nx,j,2,m)
        bx(nx,j,nz,m)=bx(nx1,j,nz,m)
        !
        by(nx,j,1,m)=by(nx,j,2,m)
        by(nx,j,nz,m)=by(nx1,j,nz,m)
        !
        bz(nx,j,1,m)=bz(nx,j,2,m)
        bz(nx,j,nz,m)=bz(nx1,j,nz,m)
        !
    enddo
    !
    do i=2,nx
        qrho(i,1,1,m)=qrho(i,1,2,m)
        qrho(i,ny,1,m)=qrho(i,ny,2,m)
        qrho(i,1,nz,m)=qrho(i,1,nz1,m)
        qrho(i,ny,nz,m)=qrho(i,ny,nz1,m)
        !
        hrho(i,1,1,m)=hrho(i,1,2,m)
        hrho(i,ny,1,m)=hrho(i,ny,2,m)
        hrho(i,1,nz,m)=hrho(i,1,nz1,m)
        hrho(i,ny,nz,m)=hrho(i,ny,nz1,m)
        !
        orho(i,1,1,m)=orho(i,1,2,m)
        orho(i,ny,1,m)=orho(i,ny,2,m)
        orho(i,1,nz,m)=orho(i,1,nz1,m)
        orho(i,ny,nz,m)=orho(i,ny,nz1,m)
        !
        qpresx(i,1,1,m)=qpresx(i,1,2,m)
        qpresx(i,ny,1,m)=qpresx(i,ny,2,m)
        qpresx(i,1,nz,m)=qpresx(i,1,nz1,m)
        qpresx(i,ny,nz,m)=qpresx(i,ny,nz1,m)
        qpresy(i,1,1,m)=qpresy(i,1,2,m)
        qpresy(i,ny,1,m)=qpresy(i,ny,2,m)
        qpresy(i,1,nz,m)=qpresy(i,1,nz1,m)
        qpresy(i,ny,nz,m)=qpresy(i,ny,nz1,m)
        qpresz(i,1,1,m)=qpresz(i,1,2,m)
        qpresz(i,ny,1,m)=qpresz(i,ny,2,m)
        qpresz(i,1,nz,m)=qpresz(i,1,nz1,m)
        qpresz(i,ny,nz,m)=qpresz(i,ny,nz1,m)
        !
        hpresx(i,1,1,m)=hpresx(i,1,2,m)
        hpresx(i,ny,1,m)=hpresx(i,ny,2,m)
        hpresx(i,1,nz,m)=hpresx(i,1,nz1,m)
        hpresx(i,ny,nz,m)=hpresx(i,ny,nz1,m)
        hpresy(i,1,1,m)=hpresy(i,1,2,m)
        hpresy(i,ny,1,m)=hpresy(i,ny,2,m)
        hpresy(i,1,nz,m)=hpresy(i,1,nz1,m)
        hpresy(i,ny,nz,m)=hpresy(i,ny,nz1,m)
        hpresz(i,1,1,m)=hpresz(i,1,2,m)
        hpresz(i,ny,1,m)=hpresz(i,ny,2,m)
        hpresz(i,1,nz,m)=hpresz(i,1,nz1,m)
        hpresz(i,ny,nz,m)=hpresz(i,ny,nz1,m)
        !
        opresx(i,1,1,m)=opresx(i,1,2,m)
        opresx(i,ny,1,m)=opresx(i,ny,2,m)
        opresx(i,1,nz,m)=opresx(i,1,nz1,m)
        opresx(i,ny,nz,m)=opresx(i,ny,nz1,m)
        opresy(i,1,1,m)=opresy(i,1,2,m)
        opresy(i,ny,1,m)=opresy(i,ny,2,m)
        opresy(i,1,nz,m)=opresy(i,1,nz1,m)
        opresy(i,ny,nz,m)=opresy(i,ny,nz1,m)
        opresz(i,1,1,m)=opresz(i,1,2,m)
        opresz(i,ny,1,m)=opresz(i,ny,2,m)
        opresz(i,1,nz,m)=opresz(i,1,nz1,m)
        opresz(i,ny,nz,m)=opresz(i,ny,nz1,m)
        !
        epres(i,1,1,m)=epres(i,1,2,m)
        epres(i,ny,1,m)=epres(i,ny,2,m)
        epres(i,1,nz,m)=epres(i,1,nz1,m)
        epres(i,ny,nz,m)=epres(i,ny,nz1,m)
        !
        !       off diagonal elements
        !
        qpresxy(i,1,1,m)=qpresxy(i,1,2,m)
        qpresxy(i,ny,1,m)=qpresxy(i,ny,2,m)
        qpresxy(i,1,nz,m)=qpresxy(i,1,nz1,m)
        qpresxy(i,ny,nz,m)=qpresxy(i,ny,nz1,m)
        qpresxz(i,1,1,m)=qpresxz(i,1,2,m)
        qpresxz(i,ny,1,m)=qpresxz(i,ny,2,m)
        qpresxz(i,1,nz,m)=qpresxz(i,1,nz1,m)
        qpresxz(i,ny,nz,m)=qpresxz(i,ny,nz1,m)
        qpresyz(i,1,1,m)=qpresyz(i,1,2,m)
        qpresyz(i,ny,1,m)=qpresyz(i,ny,2,m)
        qpresyz(i,1,nz,m)=qpresyz(i,1,nz1,m)
        qpresyz(i,ny,nz,m)=qpresyz(i,ny,nz1,m)
        !
        hpresxy(i,1,1,m)=hpresxy(i,1,2,m)
        hpresxy(i,ny,1,m)=hpresxy(i,ny,2,m)
        hpresxy(i,1,nz,m)=hpresxy(i,1,nz1,m)
        hpresxy(i,ny,nz,m)=hpresxy(i,ny,nz1,m)
        hpresxz(i,1,1,m)=hpresxz(i,1,2,m)
        hpresxz(i,ny,1,m)=hpresxz(i,ny,2,m)
        hpresxz(i,1,nz,m)=hpresxz(i,1,nz1,m)
        hpresxz(i,ny,nz,m)=hpresxz(i,ny,nz1,m)
        hpresyz(i,1,1,m)=hpresyz(i,1,2,m)
        hpresyz(i,ny,1,m)=hpresyz(i,ny,2,m)
        hpresyz(i,1,nz,m)=hpresyz(i,1,nz1,m)
        hpresyz(i,ny,nz,m)=hpresyz(i,ny,nz1,m)
        !
        opresxy(i,1,1,m)=opresxy(i,1,2,m)
        opresxy(i,ny,1,m)=opresxy(i,ny,2,m)
        opresxy(i,1,nz,m)=opresxy(i,1,nz1,m)
        opresxy(i,ny,nz,m)=opresxy(i,ny,nz1,m)
        opresxz(i,1,1,m)=opresxz(i,1,2,m)
        opresxz(i,ny,1,m)=opresxz(i,ny,2,m)
        opresxz(i,1,nz,m)=opresxz(i,1,nz1,m)
        opresxz(i,ny,nz,m)=opresxz(i,ny,nz1,m)
        opresyz(i,1,1,m)=opresyz(i,1,2,m)
        opresyz(i,ny,1,m)=opresyz(i,ny,2,m)
        opresyz(i,1,nz,m)=opresyz(i,1,nz1,m)
        opresyz(i,ny,nz,m)=opresyz(i,ny,nz1,m)
        !
        qpx(i,1,1,m)=qpx(i,1,2,m)
        qpx(i,ny,1,m)=qpx(i,ny,2,m)
        qpx(i,1,nz,m)=qpx(i,1,nz1,m)
        qpx(i,ny,nz,m)=qpx(i,ny,nz1,m)
        !
        qpy(i,1,1,m)=qpy(i,1,2,m)
        qpy(i,ny,1,m)=qpy(i,ny,2,m)
        qpy(i,1,nz,m)=qpy(i,1,nz1,m)
        qpy(i,ny,nz,m)=qpy(i,ny,nz1,m)
        !
        qpz(i,1,1,m)=qpz(i,1,2,m)
        qpz(i,ny,1,m)=qpz(i,ny,2,m)
        qpz(i,1,nz,m)=qpz(i,1,nz1,m)
        qpz(i,ny,nz,m)=qpz(i,ny,nz1,m)
        !
        hpx(i,1,1,m)=hpx(i,1,2,m)
        hpx(i,ny,1,m)=hpx(i,ny,2,m)
        hpx(i,1,nz,m)=hpx(i,1,nz1,m)
        hpx(i,ny,nz,m)=hpx(i,ny,nz1,m)
        !
        hpy(i,1,1,m)=hpy(i,1,2,m)
        hpy(i,ny,1,m)=hpy(i,ny,2,m)
        hpy(i,1,nz,m)=hpy(i,1,nz1,m)
        hpy(i,ny,nz,m)=hpy(i,ny,nz1,m)
        !
        hpz(i,1,1,m)=hpz(i,1,2,m)
        hpz(i,ny,1,m)=hpz(i,ny,2,m)
        hpz(i,1,nz,m)=hpz(i,1,nz1,m)
        hpz(i,ny,nz,m)=hpz(i,ny,nz1,m)
        !
        opx(i,1,1,m)=opx(i,1,2,m)
        opx(i,ny,1,m)=opx(i,ny,2,m)
        opx(i,1,nz,m)=opx(i,1,nz1,m)
        opx(i,ny,nz,m)=opx(i,ny,nz1,m)
        !
        opy(i,1,1,m)=opy(i,1,2,m)
        opy(i,ny,1,m)=opy(i,ny,2,m)
        opy(i,1,nz,m)=opy(i,1,nz1,m)
        opy(i,ny,nz,m)=opy(i,ny,nz1,m)
        !
        opz(i,1,1,m)=opz(i,1,2,m)
        opz(i,ny,1,m)=opz(i,ny,2,m)
        opz(i,1,nz,m)=opz(i,1,nz1,m)
        opz(i,ny,nz,m)=opz(i,ny,nz1,m)
        !
        bx(i,1,1,m)=bx(i,1,2,m)
        bx(i,ny,1,m)=bx(i,ny,2,m)
        bx(i,1,nz,m)=bx(i,1,nz1,m)
        bx(i,ny,nz,m)=bx(i,ny,nz1,m)
        !
        by(i,1,1,m)=by(i,1,2,m)
        by(i,ny,1,m)=by(i,ny,2,m)
        by(i,1,nz,m)=by(i,1,nz1,m)
        by(i,ny,nz,m)=by(i,ny,nz1,m)
        !
        bz(i,1,1,m)=bz(i,1,2,m)
        bz(i,ny,1,m)=bz(i,ny,2,m)
        bz(i,1,nz,m)=bz(i,1,nz1,m)
        bz(i,ny,nz,m)=bz(i,ny,nz1,m)
        !
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine bndry_moon(qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz, &
    qpx,qpy,qpz,rmassq, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz, &
    hpx,hpy,hpz,rmassh, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz, &
    opx,opy,opz,rmasso,epres,bx,by,bz, &
    m,nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
    ijsrf,numsrf,ijmid,nummid,ijzero, &
    numzero,mbndry,msrf,mmid,mzero, &
    qden_moon,hden_moon,oden_moon,cs_moon,gamma,ti_te, &
    vx_moon,vy_moon,vz_moon, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    isotropic)
    
    !
    !     this routine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
    common /moon/xmoon,ymoon,zmoon,rmoon,b0_moon, &
    xdip_moon,ydip_moon,zdip_moon,offset
    !
    tempi=cs_moon**2/gamma
    alpha_m=6.
    rvx=vx_moon
    rvy=vy_moon
    rvz=vz_moon
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !     set atmospheric boundary points - specify magnetic field profile
    !
    do n=1,nummid(m)
        !
        !        get coords of point
        !
        i=ijmid(m,1,n)
        j=ijmid(m,2,n)
        k=ijmid(m,3,n)
        !
        qrho(i,j,k,m)=parm_mid(m,1,n)
        hrho(i,j,k,m)=parm_mid(m,2,n)
        orho(i,j,k,m)=parm_mid(m,3,n)
        !
        qpresx(i,j,k,m)=parm_mid(m,4,n)
        qpresy(i,j,k,m)=qpresx(i,j,k,m)
        qpresz(i,j,k,m)=qpresx(i,j,k,m)*0.9

        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        !
        hpresx(i,j,k,m)=parm_mid(m,5,n)
        hpresy(i,j,k,m)=hpresx(i,j,k,m)
        hpresz(i,j,k,m)=hpresx(i,j,k,m)*0.9

        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        !
        opresx(i,j,k,m)=parm_mid(m,6,n)
        opresy(i,j,k,m)=opresx(i,j,k,m)
        opresz(i,j,k,m)=opresx(i,j,k,m)*0.9

        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        epres(i,j,k,m)=parm_mid(m,7,n)
        !
        !        qrho(i,j,k,m)=amax1(0.25*qrho(i,j,k,m),parm_mid(m,1,n))
        !        hrho(i,j,k,m)=amax1(0.25*hrho(i,j,k,m),parm_mid(m,2,n))
        !        orho(i,j,k,m)=amax1(0.25*orho(i,j,k,m),parm_mid(m,3,n))
        !
        !        qpres(i,j,k,m)=amax1(0.25*qpres(i,j,k,m),parm_mid(m,4,n))
        !        hpres(i,j,k,m)=amax1(0.25*hpres(i,j,k,m),parm_mid(m,5,n))
        !        opres(i,j,k,m)=amax1(0.25*opres(i,j,k,m),parm_mid(m,6,n))
        !        epres(i,j,k,m)=amax1(0.25*epres(i,j,k,m),parm_mid(m,7,n))
        !
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
        !        bx(i,j,k,m)=sbx_wind
        !        by(i,j,k,m)=0.
        !        bz(i,j,k,m)=0.
        !
    enddo
    !
    !
    !      reset interior points
    !
    do n=1,numzero(m)
        i=ijzero(m,1,n)
        j=ijzero(m,2,n)
        k=ijzero(m,3,n)
        !
        qrho(i,j,k,m)=parm_zero(m,1,n)
        hrho(i,j,k,m)=parm_zero(m,2,n)
        orho(i,j,k,m)=parm_zero(m,3,n)
        !
        qpresx(i,j,k,m)=parm_zero(m,4,n)
        qpresy(i,j,k,m)=qpresx(i,j,k,m)
        qpresz(i,j,k,m)=qpresx(i,j,k,m)*0.9

        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        !
        hpresx(i,j,k,m)=parm_zero(m,5,n)
        hpresy(i,j,k,m)=hpresx(i,j,k,m)
        hpresz(i,j,k,m)=hpresx(i,j,k,m)*0.9

        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        !
        opresx(i,j,k,m)=parm_zero(m,6,n)
        opresy(i,j,k,m)=opresx(i,j,k,m)
        opresz(i,j,k,m)=opresx(i,j,k,m)*0.9

        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        epres(i,j,k,m)=parm_zero(m,7,n)
        !
        !
        !        qrho(i,j,k,m)=amax1(0.10*qrho(i,j,k,m),parm_zero(m,1,n))
        !        hrho(i,j,k,m)=amax1(0.10*hrho(i,j,k,m),parm_zero(m,2,n))
        !        orho(i,j,k,m)=amax1(0.10*orho(i,j,k,m),parm_zero(m,3,n))
        !
        !        qpres(i,j,k,m)=amax1(0.05*qpres(i,j,k,m),parm_zero(m,4,n))
        !        hpres(i,j,k,m)=amax1(0.05*hpres(i,j,k,m),parm_zero(m,5,n))
        !        opres(i,j,k,m)=amax1(0.05*opres(i,j,k,m),parm_zero(m,6,n))
        !        epres(i,j,k,m)=amax1(0.05*epres(i,j,k,m),parm_zero(m,7,n))
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
    
    enddo
    !
    return
end

!
!     ***********************************************
!
subroutine set_bndry_moon_ram(rmassq,rmassh,rmasso, &
    m,nx,ny,nz,ngrd, &
    parm_srf,parm_mid,parm_zero, &
    ijsrf,numsrf,ijmid,nummid,ijzero, &
    numzero,mbndry,msrf,mmid,mzero, &
    qden_moon,hden_moon,oden_moon, &
    vx_moon,vy_moon,tempi,gamma, &
    ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_m, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !      sputtered or ram produced ionosphere
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
    !      write(6,*)'bndry_moon values',
    !    +       qden_moon,hden_moon,oden_moon,tempi,gamma,
    !    +       ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_m
    !      write(6,*)mbndry,msrf,mmid,mzero
    !
    !      set scale lengths
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !      reset indices
    !
    numsrf(m)=0
    nummid(m)=0
    numzero(m)=0
    ntot=0
    !
    !      speed for ram induced ionosphere
    !
    vtot=sqrt(vx_moon**2+vy_moon**2)+1.e-11
    !
    !       set conditions at moon's new position
    !
    do k=1,nz
        az=grd_zmin(m)+dz*(k-1)
        zm=az-zmoon
        do j=1,ny
            ay=grd_ymin(m)+dy*(j-1)
            ym=ay-ymoon
            do i=1,nx
                ax=grd_xmin(m)+dx*(i-1)
                xm=ax-xmoon
                ar_moon=sqrt(xm**2+ym**2+zm**2)
                !
                ra_moon=((ar_moon+0.3*rmoon)/(1.3*rmoon))**(-alpha_m)
                ra_moon=amin1(1.,ra_moon)
                !
                rm=(vx_moon*xm+vy_moon*ym)/vtot
                !
                if(rm.lt.0.0)then
                    xscale=1.
                else
                    xr=sqrt(rm**2+(offset*rmoon)**2)-(offset*rmoon)
                    xscale=exp(-xr/(offset*rmoon))
                endif
                !
                qden=qden_moon*ra_moon*xscale
                hden=hden_moon*ra_moon*xscale
                oden=oden_moon*ra_moon*xscale
                !
                qpress=tempi*qden
                hpress=tempi*hden
                opress=tempi*oden
                epress=(qpress+hpress+opress)/ti_te_moon
    
                if(ar_moon.le.rmoon+0.6*dx)then
                    ntot=ntot+1
                endif
                !
                if(ar_moon.le.rmoon-1.5*dx)then
                    numzero(m)=numzero(m)+1
                    if(numzero(m).gt.mzero)then
                        write(6,*)'numzero too large',m,numzero(m),mzero
                        stop
                    endif
                    ijzero(m,1,numzero(m))=i
                    ijzero(m,2,numzero(m))=j
                    ijzero(m,3,numzero(m))=k
                    !
                    parm_zero(m,1,numzero(m))=qden*rmassq
                    parm_zero(m,2,numzero(m))=hden*rmassh
                    parm_zero(m,3,numzero(m))=oden*rmasso
                    parm_zero(m,4,numzero(m))=qpress
                    parm_zero(m,5,numzero(m))=hpress
                    parm_zero(m,6,numzero(m))=opress
                    parm_zero(m,7,numzero(m))=epress
                    !
                else  if(ar_moon.le.rmoon-0.5*dx) then
                    nummid(m)=nummid(m)+1
                    if(nummid(m).gt.mmid)then
                        write(6,*)'nummid too large',m,nummid(m),mmid
                        stop
                    endif
                    ijmid(m,1,nummid(m))=i
                    ijmid(m,2,nummid(m))=j
                    ijmid(m,3,nummid(m))=k
                    !
                    parm_mid(m,1,nummid(m))=qden*rmassq
                    parm_mid(m,2,nummid(m))=hden*rmassh
                    parm_mid(m,3,nummid(m))=oden*rmasso
                    parm_mid(m,4,nummid(m))=qpress
                    parm_mid(m,5,nummid(m))=hpress
                    parm_mid(m,6,nummid(m))=opress
                    parm_mid(m,7,nummid(m))=epress
                    !
                else if(ar_moon.le.rmoon+0.6*dx) then
                    numsrf(m)=numsrf(m)+1
                    if(numsrf(m).gt.msrf)then
                        write(6,*)'numsrf too large',m,numsrf(m),msrf
                        stop
                        !                endif
                        ijsrf(m,1,numsrf(m))=i
                        ijsrf(m,2,numsrf(m))=j
                        ijsrf(m,3,numsrf(m))=k
                        !
                        parm_srf(m,1,numsrf(m))=qden*rmassq
                        parm_srf(m,2,numsrf(m))=hden*rmassh
                        parm_srf(m,3,numsrf(m))=oden*rmasso
                        parm_srf(m,4,numsrf(m))=qpress
                        parm_srf(m,5,numsrf(m))=hpress
                        parm_srf(m,6,numsrf(m))=opress
                        parm_srf(m,7,numsrf(m))=epress
                        !
                    endif
                endif
                !
            enddo
        enddo
    enddo
    !
    !         write(6,*)'total pts',ntot
    !         write(6,*)'moon bndry_m pts',m,numsrf(m),nummid(m),numzero(m)
    !
    return
end

!
!     ***********************************************
!
subroutine set_bndry_moon(rmassq,rmassh,rmasso, &
    m,nx,ny,nz,ngrd, &
    parm_srf,parm_mid,parm_zero, &
    ijsrf,numsrf,ijmid,nummid,ijzero, &
    numzero,mbndry,msrf,mmid,mzero, &
    qden_moon,hden_moon,oden_moon,tempi,gamma, &
    ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_m, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !    photoionization ionosphere
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
    !      write(6,*)'bndry_moon values',
    !    +       qden_moon,hden_moon,oden_moon,tempi,gamma,
    !    +       ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_m
    !      write(6,*)mbndry,msrf,mmid,mzero
    !
    !      set scale lengths
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !      reset indices
    !
    numsrf(m)=0
    nummid(m)=0
    numzero(m)=0
    ntot=0
    !
    !
    !       set conditions at moon's new position
    !
    do k=1,nz
        az=grd_zmin(m)+dz*(k-1)
        zm=az-zmoon
        do j=1,ny
            ay=grd_ymin(m)+dy*(j-1)
            ym=ay-ymoon
            do i=1,nx
                ax=grd_xmin(m)+dx*(i-1)
                xm=ax-xmoon
                ar_moon=sqrt(xm**2+ym**2+zm**2)
                !
                ra_moon=((ar_moon+0.3*rmoon)/(1.3*rmoon))**(-alpha_m)
                ra_moon=amin1(1.,ra_moon)
                !
                !
                if(xm.lt.0)then
                    xscale=1.
                else
                    xr=sqrt(xm**2+(offset*rmoon)**2)-(offset*rmoon)
                    xscale=exp(-xr/(offset*rmoon))
                endif
                !
                qden=qden_moon*ra_moon*xscale
                hden=hden_moon*ra_moon*xscale
                oden=oden_moon*ra_moon*xscale
                !
                qpress=tempi*qden
                hpress=tempi*hden
                opress=tempi*oden
                epress=(qpress+hpress+opress)/ti_te_moon
    
                if(ar_moon.le.rmoon+0.6*dx)then
                    ntot=ntot+1
                endif
                !
                if(ar_moon.le.rmoon-1.5*dx)then
                    numzero(m)=numzero(m)+1
                    if(numzero(m).gt.mzero)then
                        write(6,*)'numzero too large',m,numzero(m),mzero
                        stop
                    endif
                    ijzero(m,1,numzero(m))=i
                    ijzero(m,2,numzero(m))=j
                    ijzero(m,3,numzero(m))=k
                    !
                    parm_zero(m,1,numzero(m))=qden*rmassq
                    parm_zero(m,2,numzero(m))=hden*rmassh
                    parm_zero(m,3,numzero(m))=oden*rmasso
                    parm_zero(m,4,numzero(m))=qpress
                    parm_zero(m,5,numzero(m))=hpress
                    parm_zero(m,6,numzero(m))=opress
                    parm_zero(m,7,numzero(m))=epress
                    !
                else  if(ar_moon.le.rmoon-0.5*dx) then
                    nummid(m)=nummid(m)+1
                    if(nummid(m).gt.mmid)then
                        write(6,*)'nummid too large',m,nummid(m),mmid
                        stop
                    endif
                    ijmid(m,1,nummid(m))=i
                    ijmid(m,2,nummid(m))=j
                    ijmid(m,3,nummid(m))=k
                    !
                    parm_mid(m,1,nummid(m))=qden*rmassq
                    parm_mid(m,2,nummid(m))=hden*rmassh
                    parm_mid(m,3,nummid(m))=oden*rmasso
                    parm_mid(m,4,nummid(m))=qpress
                    parm_mid(m,5,nummid(m))=hpress
                    parm_mid(m,6,nummid(m))=opress
                    parm_mid(m,7,nummid(m))=epress
                    !
                else if(ar_moon.le.rmoon+0.6*dx) then
                    numsrf(m)=numsrf(m)+1
                    if(numsrf(m).gt.msrf)then
                        write(6,*)'numsrf too large',m,numsrf(m),msrf
                        stop
                        !                endif
                        ijsrf(m,1,numsrf(m))=i
                        ijsrf(m,2,numsrf(m))=j
                        ijsrf(m,3,numsrf(m))=k
                        !
                        parm_srf(m,1,numsrf(m))=qden*rmassq
                        parm_srf(m,2,numsrf(m))=hden*rmassh
                        parm_srf(m,3,numsrf(m))=oden*rmasso
                        parm_srf(m,4,numsrf(m))=qpress
                        parm_srf(m,5,numsrf(m))=hpress
                        parm_srf(m,6,numsrf(m))=opress
                        parm_srf(m,7,numsrf(m))=epress
                        !
                    endif
                endif
                !
            enddo
        enddo
    enddo
    !
    !         write(6,*)'total pts',ntot
    !         write(6,*)'moon bndry_m pts',m,numsrf(m),nummid(m),numzero(m)
    !
    return
end

!
!     ********************************************
!
subroutine bndry_grd_core( &
    qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,opx,opy,opz, &
    epres,bx,by,bz, &
    nx,ny,nz,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    qrho_n,qpresx_n,qpresy_n,qpresz_n, &
    qpresxy_n,qpresxz_n,qpresyz_n, &
    qpx_n,qpy_n,qpz_n, &
    hrho_n,hpresx_n,hpresy_n,hpresz_n, &
    hpresxy_n,hpresxz_n,hpresyz_n, &
    hpx_n,hpy_n,hpz_n, &
    orho_n,opresx_n,opresy_n,opresz_n, &
    opresxy_n,opresxz_n,opresyz_n, &
    opx_n,opy_n,opz_n, &
    epres_n,bx_n,by_n,bz_n, &
    nx_n,ny_n,nz_n,ngrd_n,main_ngrd, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    
    !
    !     this routine applies boundary conditins to all grid types
    !     of the system and at any irregular boundaries
    !
    !
    dimension &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    ! 
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd)
    dimension qrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpy_n(nx_n,ny_n,nz_n,ngrd_n),qpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    hrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpx_n(nx_n,ny_n,nz_n,ngrd_n),hpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    orho_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opx_n(nx_n,ny_n,nz_n,ngrd_n),opy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    bx_n(nx_n,ny_n,nz_n,ngrd_n),by_n(nx_n,ny_n,nz_n,ngrd_n), &
    bz_n(nx_n,ny_n,nz_n,ngrd_n), &
    epres_n(nx_n,ny_n,nz_n,ngrd_n)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    call corer_grds(qrho,nx,ny,nz,ngrd,main_ngrd, &
    qrho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresx,nx,ny,nz,ngrd,main_ngrd, &
    qpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresy,nx,ny,nz,ngrd,main_ngrd, &
    qpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresz,nx,ny,nz,ngrd,main_ngrd, &
    qpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresxy,nx,ny,nz,ngrd,main_ngrd, &
    qpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresxz,nx,ny,nz,ngrd,main_ngrd, &
    qpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresyz,nx,ny,nz,ngrd,main_ngrd, &
    qpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpx,nx,ny,nz,ngrd,main_ngrd, &
    qpx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpy,nx,ny,nz,ngrd,main_ngrd, &
    qpy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpz,nx,ny,nz,ngrd,main_ngrd, &
    qpz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(hrho,nx,ny,nz,ngrd,main_ngrd, &
    hrho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresx,nx,ny,nz,ngrd,main_ngrd, &
    hpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresy,nx,ny,nz,ngrd,main_ngrd, &
    hpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresz,nx,ny,nz,ngrd,main_ngrd, &
    hpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresxy,nx,ny,nz,ngrd,main_ngrd, &
    hpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresxz,nx,ny,nz,ngrd,main_ngrd, &
    hpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresyz,nx,ny,nz,ngrd,main_ngrd, &
    hpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpx,nx,ny,nz,ngrd,main_ngrd, &
    hpx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpy,nx,ny,nz,ngrd,main_ngrd, &
    hpy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpz,nx,ny,nz,ngrd,main_ngrd, &
    hpz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(orho,nx,ny,nz,ngrd,main_ngrd, &
    orho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresx,nx,ny,nz,ngrd,main_ngrd, &
    opresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresy,nx,ny,nz,ngrd,main_ngrd, &
    opresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresz,nx,ny,nz,ngrd,main_ngrd, &
    opresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresxy,nx,ny,nz,ngrd,main_ngrd, &
    opresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresxz,nx,ny,nz,ngrd,main_ngrd, &
    opresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresyz,nx,ny,nz,ngrd,main_ngrd, &
    opresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opx,nx,ny,nz,ngrd,main_ngrd, &
    opx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opy,nx,ny,nz,ngrd,main_ngrd, &
    opy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opz,nx,ny,nz,ngrd,main_ngrd, &
    opz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(epres,nx,ny,nz,ngrd,main_ngrd, &
    epres_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(bx,nx,ny,nz,ngrd,main_ngrd, &
    bx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(by,nx,ny,nz,ngrd,main_ngrd, &
    by_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(bz,nx,ny,nz,ngrd,main_ngrd, &
    bz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    
    !
    return
end

!
!       ************************************************
!
subroutine corer_grds(rho,nx,ny,nz,ngrd,main_ngrd, &
    rho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    dimension rho(nx,ny,nz,ngrd),rho_n(nx_n,ny_n,nz_n,ngrd_n)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    !    fills in the core of the main grid in position mb
    !       using data from the finer (small) grid in position ms
    !
    mb=main_ngrd
    ms=ngrd_n
    !
    !     set limits for coarse grid
    !
    dx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
    dy=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
    dz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
    sx=(grd_xmax_n(ms)-grd_xmin_n(ms))/(nx_n-1.)
    sy=(grd_ymax_n(ms)-grd_ymin_n(ms))/(ny_n-1.)
    sz=(grd_zmax_n(ms)-grd_zmin_n(ms))/(nz_n-1.)
    nskip=dx/sx
    !
    ai=1.+(grd_xmin_n(ms)-grd_xmin(mb))/dx
    nbx1=ai
    aj=1.+(grd_ymin_n(ms)-grd_ymin(mb))/dy
    nby1=aj
    ak=1.+(grd_zmin_n(ms)-grd_zmin(mb))/dz
    nbz1=ak
    !
    !$omp  parallel do
    do ks=1+nskip,nz_n-nskip,nskip
        kb=nbz1+(ks-1)/nskip
        do js=1+nskip,ny_n-nskip,nskip
            jb=nby1+(js-1)/nskip
            do is=1+nskip,nx_n-nskip,nskip
                ib=nbx1+(is-1)/nskip
                !
                !         rho(ib,jb,kb,mb)=rho_n(is,js,ks,ms)
                !
                rho(ib,jb,kb,mb)=(6.*rho_n(is,js,ks,ms)+ &
                rho_n(is+1,js,ks,ms)+rho_n(is-1,js,ks,ms)+ &
                rho_n(is,js+1,ks,ms)+rho_n(is,js-1,ks,ms)+ &
                rho_n(is,js,ks+1,ms)+rho_n(is,js,ks-1,ms))/12.
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine bndry_corer( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz,opx,opy,opz, &
    epres,bx,by,bz, &
    qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz, &
    opresxy,opresxz,opresyz, &
    nx,ny,nz,ngrd,m, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
    !
    !     this routine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    !
    dimension &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    call corer(qrho,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresxy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresxz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresyz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call corer(hrho,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresxy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresxz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresyz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call corer(orho,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresxy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresxz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresyz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(epres,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(bx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(by,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(bz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    
    return
end

!
!       ************************************************
!
subroutine corer(rho,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      m is the coring subject grid index
    !
    dimension rho(nx,ny,nz,ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !    fills in the core of the coarser (big) grid in position nb
    !       using data from the finer (small) grid in position ns
    !
    !     set limits for coarse grid
    !
    mb=m+1
    ms=m
    sx=(grd_xmax(ms)-grd_xmin(ms))/(nx-1.)
    sy=(grd_ymax(ms)-grd_ymin(ms))/(ny-1.)
    sz=(grd_zmax(ms)-grd_zmin(ms))/(nz-1.)
    dx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
    dy=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
    dz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
    !
    nskip=dx/sx
    !
    ai=1.+(grd_xmin(ms)-grd_xmin(mb))/dx
    nbx1=ai+.5
    aj=1.+(grd_ymin(ms)-grd_ymin(mb))/dy
    nby1=aj+.5
    ak=1.+(grd_zmin(ms)-grd_zmin(mb))/dz
    nbz1=ak+.5
    !
    !$omp  parallel do
    do  ks=1+nskip,nz-nskip,nskip
        kb=nbz1+(ks-1)/nskip
        do  js=1+nskip,ny-nskip,nskip
            jb=nby1+(js-1)/nskip
            do  is=1+nskip,nx-nskip,nskip
                ib=nbx1+(is-1)/nskip
                !
                !       old version with no averaging
                !        rho(ib,jb,kb,mb)=rho(is,js,ks,ms)
                !
                !       smooth to make compatible with coarse grid
                !
                rho(ib,jb,kb,mb)=(6.*rho(is,js,ks,ms)+ &
                rho(is+1,js,ks,ms)+rho(is-1,js,ks,ms)+ &
                rho(is,js+1,ks,ms)+rho(is,js-1,ks,ms)+ &
                rho(is,js,ks+1,ms)+rho(is,js,ks-1,ms))/12.
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine bndry_flanks( &
    wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
    wrkqpresx,wrkqpresy,wrkqpresz, &
    wrkqpresxy,wrkqpresxz,wrkqpresyz, &
    wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
    wrkhpresx,wrkhpresy,wrkhpresz, &
    wrkhpresxy,wrkhpresxz,wrkhpresyz, &
    wrkorho,wrkopx,wrkopy,wrkopz, &
    wrkopresx,wrkopresy,wrkopresz, &
    wrkopresxy,wrkopresxz,wrkopresyz, &
    wrkepres, wrkbx,wrkby,wrkbz, &
    qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz, &
    hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz, &
    orho,opx,opy,opz,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz, &
    epres, bx,by,bz, &
    oldqrho,oldqpx,oldqpy,oldqpz, &
    oldqpresx,oldqpresy,oldqpresz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    oldhrho,oldhpx,oldhpy,oldhpz, &
    oldhpresx,oldhpresy,oldhpresz, &
    oldhpresxy,oldhpresxz,oldhpresyz, &
    oldorho,oldopx,oldopy,oldopz, &
    oldopresx,oldopresy,oldopresz, &
    oldopresxy,oldopresxz,oldopresyz, &
    oldepres, oldbx,oldby,oldbz,work, &
    nx,ny,nz,ngrd,m,t_old,t_new,t, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    
    !
    !     this routine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpresx(nx,ny,nz,ngrd), &
    wrkqpresy(nx,ny,nz,ngrd),wrkqpresz(nx,ny,nz,ngrd), &
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd), &
    wrkhrho(nx,ny,nz,ngrd),wrkhpresx(nx,ny,nz,ngrd), &
    wrkhpresy(nx,ny,nz,ngrd),wrkhpresz(nx,ny,nz,ngrd), &
    wrkhpx(nx,ny,nz,ngrd),wrkhpy(nx,ny,nz,ngrd), &
    wrkhpz(nx,ny,nz,ngrd), &
    wrkorho(nx,ny,nz,ngrd),wrkopresx(nx,ny,nz,ngrd), &
    wrkopresy(nx,ny,nz,ngrd),wrkopresz(nx,ny,nz,ngrd), &
    wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
    wrkopz(nx,ny,nz,ngrd), &
    wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
    wrkbz(nx,ny,nz,ngrd), &
    wrkepres(nx,ny,nz,ngrd)
    dimension oldqrho(nx,ny,nz,ngrd),oldqpresx(nx,ny,nz,ngrd), &
    oldqpresy(nx,ny,nz,ngrd),oldqpresz(nx,ny,nz,ngrd), &
    oldqpx(nx,ny,nz,ngrd),oldqpy(nx,ny,nz,ngrd), &
    oldqpz(nx,ny,nz,ngrd), &
    oldhrho(nx,ny,nz,ngrd),oldhpresx(nx,ny,nz,ngrd), &
    oldhpresy(nx,ny,nz,ngrd),oldhpresz(nx,ny,nz,ngrd), &
    oldhpx(nx,ny,nz,ngrd),oldhpy(nx,ny,nz,ngrd), &
    oldhpz(nx,ny,nz,ngrd), &
    oldorho(nx,ny,nz,ngrd),oldopresx(nx,ny,nz,ngrd), &
    oldopresy(nx,ny,nz,ngrd),oldopresz(nx,ny,nz,ngrd), &
    oldopx(nx,ny,nz,ngrd),oldopy(nx,ny,nz,ngrd), &
    oldopz(nx,ny,nz,ngrd), &
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
    oldbz(nx,ny,nz,ngrd), &
    oldepres(nx,ny,nz,ngrd)
    !
    dimension qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd),hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd)
    !
    dimension oldqpresxy(nx,ny,nz,ngrd),oldqpresxz(nx,ny,nz,ngrd), &
    oldqpresyz(nx,ny,nz,ngrd),oldhpresxy(nx,ny,nz,ngrd), &
    oldhpresxz(nx,ny,nz,ngrd),oldhpresyz(nx,ny,nz,ngrd), &
    oldopresxy(nx,ny,nz,ngrd),oldopresxz(nx,ny,nz,ngrd), &
    oldopresyz(nx,ny,nz,ngrd)
    dimension wrkqpresxy(nx,ny,nz,ngrd),wrkqpresxz(nx,ny,nz,ngrd), &
    wrkqpresyz(nx,ny,nz,ngrd),wrkhpresxy(nx,ny,nz,ngrd), &
    wrkhpresxz(nx,ny,nz,ngrd),wrkhpresyz(nx,ny,nz,ngrd), &
    wrkopresxy(nx,ny,nz,ngrd),wrkopresxz(nx,ny,nz,ngrd), &
    wrkopresyz(nx,ny,nz,ngrd)
    !
    dimension work(nx,ny,nz)
    !
    dimension t_old(ngrd),t_new(ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    call flanks(wrkqrho,qrho,oldqrho,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresx,qpresx,oldqpresx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresy,qpresy,oldqpresy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresz,qpresz,oldqpresz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresxy,qpresxy,oldqpresxy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresxz,qpresxz,oldqpresxz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresyz,qpresyz,oldqpresyz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpx,qpx,oldqpx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpy,qpy,oldqpy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpz,qpz,oldqpz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkhrho,hrho,oldhrho,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresx,hpresx,oldhpresx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresy,hpresy,oldhpresy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresz,hpresz,oldhpresz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresxy,hpresxy,oldhpresxy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresxz,hpresxz,oldhpresxz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresyz,hpresyz,oldhpresyz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpx,hpx,oldhpx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpy,hpy,oldhpy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpz,hpz,oldhpz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkorho,orho,oldorho,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresx,opresx,oldopresx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresy,opresy,oldopresy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresz,opresz,oldopresz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresxy,opresxy,oldopresxy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresxz,opresxz,oldopresxz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresyz,opresyz,oldopresyz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopx,opx,oldopx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopy,opy,oldopy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopz,opz,oldopz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkepres,epres,oldepres,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkbx,bx,oldbx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkby,by,oldby,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkbz,bz,oldbz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    return
end

!
!    ******************************************
!
subroutine flanks(wrkrho,rho,oldrho,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    dimension rho(nx,ny,nz,ngrd),oldrho(nx,ny,nz,ngrd), &
    wrkrho(nx,ny,nz,ngrd),work(nx,ny,nz)
    dimension t_new(ngrd),t_old(ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !     m is the suject array index
    !
    !      sets the outer boundaries of fine grid at position ns
    !      using data  from coarse grid at position nb
    !
    mb=m+1
    ms=m
    !
    t1=(t_new(mb)-t)/(t_new(mb)-t_old(mb))
    t2=1.-t1
    !
    delx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
    dely=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
    delz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
    !
    sx=(grd_xmax(ms)-grd_xmin(ms))/(nx-1.)
    sy=(grd_ymax(ms)-grd_ymin(ms))/(ny-1.)
    sz=(grd_zmax(ms)-grd_zmin(ms))/(nz-1.)
    
    !
    !     time interpolation
    !
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                work(i,j,k)=oldrho(i,j,k,mb)*t1+rho(i,j,k,mb)*t2
            enddo
        enddo
    enddo
    !
    !     interpolate onto grid
    !
    !     do  bottom panels
    !
    !$omp  parallel do
    do kk=2,1,-1
        ks=kk
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        frac=(kk-1)/2.
        frac1=1.-frac
        !
        do js=1,ny
            ay=grd_ymin(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            do is=1,nx
                ax=grd_xmin(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    !     do top panels
    !
    !$omp  parallel do
    do kk=2,1,-1
        ks=nz-kk+1
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        frac=(kk-1)/2.
        frac1=1.-frac
        !
        do js=1,ny
            ay=grd_ymin(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            do is=1,nx
                ax=grd_xmin(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
            enddo
        enddo
    enddo
    !
    !     do  front panels
    !
    !$omp  parallel do
    !
    do ks=1,nz
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        do jj=2,1,-1
            js=jj
            ay=grd_ymin(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            frac=(jj-1)/2.
            frac1=1.-frac
            do is=1,nx
                ax=grd_xmin(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    !     do  back panels
    !
    !$omp  parallel do
    !
    do ks=1,nz
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        do jj=2,1,-1
            js=ny-jj+1
            ay=grd_ymin(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            frac=(jj-1)/2.
            frac1=1.-frac
            do is=1,nx
                ax=grd_xmin(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    !     do  leftfront panels
    !
    !$omp  parallel do
    do ks=1,nz
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        do js=1,ny
            ay=grd_ymin(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            do ii=2,1,-1
                is=ii
                ax=grd_xmin(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                frac=(ii-1)/2.
                frac1=1.-frac
                !
                wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    do ks=1,nz
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        do js=1,ny
            ay=grd_ymin(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            do ii=2,1,-1
                is=nx-ii+1
                ax=grd_xmin(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                frac=(ii-1)/2.
                frac1=1.-frac
                !
                wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    return
end

!
!     ********************************************
!
subroutine bndry_grds( &
    qrho_n,qpresx_n,qpresy_n,qpresz_n, &
    qpresxy_n,qpresxz_n,qpresyz_n,qpx_n,qpy_n,qpz_n, &
    hrho_n,hpresx_n,hpresy_n,hpresz_n, &
    hpresxy_n,hpresxz_n,hpresyz_n,hpx_n,hpy_n,hpz_n, &
    orho_n,opresx_n,opresy_n,opresz_n, &
    opresxy_n,opresxz_n,opresyz_n,opx_n,opy_n,opz_n, &
    epres_n,bx_n,by_n,bz_n, &
    nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n,t, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n, &
    qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,opx,opy,opz, &
    epres,bx,by,bz, &
    oldqrho,oldqpresx,oldqpresy,oldqpresz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    oldqpx,oldqpy,oldqpz, &
    oldhrho,oldhpresx,oldhpresy,oldhpresz, &
    oldhpresxy,oldhpresxz,oldhpresyz, &
    oldhpx,oldhpy,oldhpz, &
    oldorho,oldopresx,oldopresy,oldopresz, &
    oldopresxy,oldopresxz,oldopresyz, &
    oldopx,oldopy,oldopz, &
    oldepres,oldbx,oldby,oldbz,work, &
    nx,ny,nz,ngrd,t_old,t_new, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    
    !
    !     this routine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    dimension qrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    hrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    orho_n(nx_n,ny_n,nz_n,ngrd_n), &
    opx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    bx_n(nx_n,ny_n,nz_n,ngrd_n),by_n(nx_n,ny_n,nz_n,ngrd_n), &
    bz_n(nx_n,ny_n,nz_n,ngrd_n),epres_n(nx_n,ny_n,nz_n,ngrd_n)
    !
    dimension qrho(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd),qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd),hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresy(nx,ny,nz,ngrd), &
    opresz(nx,ny,nz,ngrd),opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    dimension oldqrho(nx,ny,nz,ngrd), &
    oldqpresx(nx,ny,nz,ngrd),oldqpresy(nx,ny,nz,ngrd), &
    oldqpresz(nx,ny,nz,ngrd),oldqpresxy(nx,ny,nz,ngrd), &
    oldqpresxz(nx,ny,nz,ngrd),oldqpresyz(nx,ny,nz,ngrd), &
    oldqpx(nx,ny,nz,ngrd),oldqpy(nx,ny,nz,ngrd), &
    oldqpz(nx,ny,nz,ngrd), &
    !
    oldhrho(nx,ny,nz,ngrd), &
    oldhpresx(nx,ny,nz,ngrd),oldhpresy(nx,ny,nz,ngrd), &
    oldhpresz(nx,ny,nz,ngrd),oldhpresxy(nx,ny,nz,ngrd), &
    oldhpresxz(nx,ny,nz,ngrd),oldhpresyz(nx,ny,nz,ngrd), &
    oldhpx(nx,ny,nz,ngrd),oldhpy(nx,ny,nz,ngrd), &
    oldhpz(nx,ny,nz,ngrd), &
    !
    oldorho(nx,ny,nz,ngrd), &
    oldopresx(nx,ny,nz,ngrd),oldopresy(nx,ny,nz,ngrd), &
    oldopresz(nx,ny,nz,ngrd),oldopresxy(nx,ny,nz,ngrd), &
    oldopresxz(nx,ny,nz,ngrd),oldopresyz(nx,ny,nz,ngrd), &
    oldopx(nx,ny,nz,ngrd),oldopy(nx,ny,nz,ngrd), &
    oldopz(nx,ny,nz,ngrd), &
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
    oldbz(nx,ny,nz,ngrd),oldepres(nx,ny,nz,ngrd)
    !
    dimension work(nx,ny,nz)
    dimension t_old(ngrd),t_new(ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    call flanks_grds(qrho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qrho,oldqrho,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresx,oldqpresx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresy,oldqpresy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresz,oldqpresz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresxy,oldqpresxy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresxz,oldqpresxz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresyz,oldqpresyz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpx,oldqpx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpy,oldqpy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpz,oldqpz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(hrho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hrho,oldhrho,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresx,oldhpresx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresy,oldhpresy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresz,oldhpresz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresxy,oldhpresxy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresxz,oldhpresxz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresyz,oldhpresyz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpx,oldhpx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpy,oldhpy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpz,oldhpz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(orho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    orho,oldorho,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresx,oldopresx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresy,oldopresy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresz,oldopresz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresxy,oldopresxy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresxz,oldopresxz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresyz,oldopresyz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opx,oldopx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opy,oldopy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opz,oldopz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(epres_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    epres,oldepres,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(bx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    bx,oldbx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(by_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    by,oldby,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(bz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    bz,oldbz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    return
end

!
!    ******************************************
!
subroutine flanks_grds( &
    rho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd, &
    rho,oldrho,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    dimension rho(nx,ny,nz,ngrd),oldrho(nx,ny,nz,ngrd), &
    rho_n(nx_n,ny_n,nz_n,ngrd_n), &
    work(nx,ny,nz)
    dimension t_new(ngrd),t_old(ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    !     m is the suject array index
    !
    !      sets the outer boundaries of fine grid at position ns
    !      using data  from coarse grid at position nb
    !
    mb=main_ngrd
    !
    t1=(t_new(mb)-t)/(t_new(mb)-t_old(mb))
    t2=1.-t1
    !
    !     write(6,*)'flanks_grds',mb,t1,t2
    !
    delx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
    dely=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
    delz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
    !
    ms=ngrd_n
    !
    sx=(grd_xmax_n(ms)-grd_xmin_n(ms))/(nx_n-1.)
    sy=(grd_ymax_n(ms)-grd_ymin_n(ms))/(ny_n-1.)
    sz=(grd_zmax_n(ms)-grd_zmin_n(ms))/(nz_n-1.)
    !
    !     time interpolation
    !
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                work(i,j,k)=oldrho(i,j,k,mb)*t1+rho(i,j,k,mb)*t2
            enddo
        enddo
    enddo
    !
    !     do bottom panels
    !
    !$omp  parallel do
    do kk=2,1,-1
        ks=kk
        az=grd_zmin_n(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        frac=(kk-1)/2.
        frac1=1.-frac
        !
        do js=1,ny_n
            ay=grd_ymin_n(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            !
            do is=1,nx_n
                ax=grd_xmin_n(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                rho_n(is,js,ks,ms)=frac*rho_n(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    !     do top panels
    !
    !$omp  parallel do
    do kk=2,1,-1
        ks=nz_n-kk+1
        az=grd_zmin_n(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        frac=(kk-1)/2.
        frac1=1.-frac
        !
        do js=1,ny_n
            ay=grd_ymin_n(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            !
            do is=1,nx_n
                ax=grd_xmin_n(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                rho_n(is,js,ks,ms)=frac*rho_n(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    
    !
    !     do front panels
    !
    !$omp  parallel do
    do ks=1,nz_n
        az=grd_zmin_n(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        do jj=2,1,-1
            js=jj
            ay=grd_ymin_n(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            frac=(jj-1)/2.
            frac1=1.-frac
            do is=1,nx_n
                ax=grd_xmin_n(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                rho_n(is,js,ks,ms)=frac*rho_n(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    !     do back panels
    !
    !$omp  parallel do
    do ks=1,nz_n
        az=grd_zmin_n(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        !
        do jj=2,1,-1
            js=ny_n-jj+1
            ay=grd_ymin_n(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            frac=(jj-1)/2.
            frac1=1.-frac
            !
            do is=1,nx_n
                ax=grd_xmin_n(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                !
                rho_n(is,js,ks,ms)=frac*rho_n(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    
    !
    !     do left panels
    !
    !$omp  parallel do
    do ks=1,nz_n
        az=grd_zmin_n(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        !
        do js=1,ny_n
            ay=grd_ymin_n(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            !
            do ii=2,1,-1
                is=ii
                ax=grd_xmin_n(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                frac=(ii-1)/2.
                frac1=1.-frac
                !
                rho_n(is,js,ks,ms)=frac*rho_n(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
                !
            enddo
        enddo
    enddo
    !
    !     do right panel
    !
    !$omp  parallel do
    do ks=1,nz_n
        az=grd_zmin_n(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        !
        do js=1,ny_n
            ay=grd_ymin_n(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            !
            do ii=2,1,-1
                is=nx_n-ii+1
                ax=grd_xmin_n(ms)+sx*(is-1)
                ai=1.+(ax-grd_xmin(mb))/delx
                ib=ai
                ib1=ib+1
                dx=ai-ib
                ddx=1.-dx
                frac=(ii-1)/2.
                frac1=1.-frac
    
                rho_n(is,js,ks,ms)=frac*rho_n(is,js,ks,ms) &
                +frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
                +work(ib,jb,kb1)*ddx*ddy*dz &
                +work(ib,jb1,kb)*ddx*dy*ddz &
                +work(ib,jb1,kb1)*ddx*dy*dz &
                +work(ib1,jb,kb)*dx*ddy*ddz &
                +work(ib1,jb,kb1)*dx*ddy*dz &
                +work(ib1,jb1,kb)*dx*dy*ddz &
                +work(ib1,jb1,kb1)*dx*dy*dz)
                !
            enddo
        enddo
    enddo
    !
    return
end

!
!    ******************************************
!
subroutine flanks_synced(rho,nx,ny,nz,ngrd,m, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    dimension rho(nx,ny,nz,ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !
    !      sets the outer boundaries of fine grid at position ns
    !      using data  from coarse grid at position nb
    !
    mb=m+1
    nb=mb
    !
    ms=m
    ns=ms
    !
    delx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
    dely=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
    delz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
    !
    sx=(grd_xmax(ms)-grd_xmin(ms))/(nx-1.)
    sy=(grd_ymax(ms)-grd_ymin(ms))/(ny-1.)
    sz=(grd_zmax(ms)-grd_zmin(ms))/(nz-1.)
    !
    !
    ks_top=nz
    ks_bot=1
    js_top=ny
    js_bot=1
    is_top=nx
    is_bot=1
    !
    kb_top=1.+(grd_zmax(ms)-grd_zmin(mb))/delz
    kb_bot=1.+(grd_zmin(ms)-grd_zmin(mb))/delz
    jb_top=1.+(grd_ymax(ms)-grd_ymin(mb))/dely
    jb_bot=1.+(grd_ymin(ms)-grd_ymin(mb))/dely
    ib_top=1.+(grd_xmax(ms)-grd_xmin(mb))/delx
    ib_bot=1.+(grd_xmin(ms)-grd_xmin(mb))/delx
    !
    !     do top and bottom panels
    !
    !$omp  parallel do
    do  js=1,ny
        ay=grd_ymin(ms)+sy*(js-1)
        aj=1.+(ay-grd_ymin(mb))/dely
        jb=aj
        jb1=jb+1
        dy=aj-jb
        ddy=1.-dy
        do  is=1,nx
            ax=grd_xmin(ms)+sx*(is-1)
            ai=1.+(ax-grd_xmin(mb))/delx
            ib=ai
            ib1=ib+1
            dx=ai-ib
            ddx=1.-dx
            !
            rho(is,js,ks_top,ns)= &
            rho(ib,jb,kb_top,nb)*ddx*ddy &
            +rho(ib,jb1,kb_top,nb)*ddx*dy &
            +rho(ib1,jb,kb_top,nb)*dx*ddy &
            +rho(ib1,jb1,kb_top,nb)*dx*dy
            !
            rho(is,js,ks_bot,ns)= &
            rho(ib,jb,kb_bot,nb)*ddx*ddy &
            +rho(ib,jb1,kb_bot,nb)*ddx*dy &
            +rho(ib1,jb,kb_bot,nb)*dx*ddy &
            +rho(ib1,jb1,kb_bot,nb)*dx*dy
            !
        enddo
    enddo
    !
    !     do front  and back panels
    !
    !$omp  parallel do
    do  ks=1,nz
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        do is=1,nx
            ax=grd_xmin(ms)+sx*(is-1)
            ai=1.+(ax-grd_xmin(mb))/delx
            ib=ai
            ib1=ib+1
            dx=ai-ib
            ddx=1.-dx
            !
            rho(is,js_top,ks,ns)= &
            rho(ib,jb_top,kb,nb)*ddx*ddz &
            +rho(ib,jb_top,kb1,nb)*ddx*dz &
            +rho(ib1,jb_top,kb,nb)*dx*ddz &
            +rho(ib1,jb_top,kb1,nb)*dx*dz
            !
            rho(is,js_bot,ks,ns)= &
            rho(ib,jb_bot,kb,nb)*ddx*ddz &
            +rho(ib,jb_bot,kb1,nb)*ddx*dz &
            +rho(ib1,jb_bot,kb,nb)*dx*ddz &
            +rho(ib1,jb_bot,kb1,nb)*dx*dz
            !
        enddo
    enddo
    !
    !
    !     do left and right panels
    !
    !$omp  parallel do
    do  ks=1,nz
        az=grd_zmin(ms)+sz*(ks-1)
        ak=1.+(az-grd_zmin(mb))/delz
        kb=ak
        kb1=kb+1
        dz=ak-kb
        ddz=1.-dz
        do  js=1,ny
            ay=grd_ymin(ms)+sy*(js-1)
            aj=1.+(ay-grd_ymin(mb))/dely
            jb=aj
            jb1=jb+1
            dy=aj-jb
            ddy=1.-dy
            !
            rho(is_top,js,ks,ns)= &
            rho(ib_top,jb,kb,nb)*ddy*ddz &
            +rho(ib_top,jb1,kb,nb)*dy*ddz &
            +rho(ib_top,jb,kb1,nb)*ddy*dz &
            +rho(ib_top,jb1,kb1,nb)*dy*dz
            !
            rho(is_bot,js,ks,ns)= &
            rho(ib_bot,jb,kb,nb)*ddy*ddz &
            +rho(ib_bot,jb1,kb,nb)*dy*ddz &
            +rho(ib_bot,jb,kb1,nb)*ddy*dz &
            +rho(ib_bot,jb1,kb1,nb)*dy*dz
            !
        enddo
    enddo
    !
    return
end

!
!       ************************************************
!
subroutine refinement(rho,nx,ny,nz,ngrd,main_grid, &
    rho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    
    !
    dimension rho(nx,ny,nz,ngrd),rho_n(nx_n,ny_n,nz_n,ngrd_n)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    
    !
    !    fills in the core of the coarser (big) grid in position nb
    !       using data from the finer (small) grid in position ns
    !
    !     set limits for coarse grid
    !
    rho_n=0.
    m=main_grid
    m_n=ngrd_n
    rho_max=0.
    rho_n_max=0.
    
    rx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    ry=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    rz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    sx=(grd_xmax_n(m_n)-grd_xmin_n(m_n))/(nx_n-1.)
    sy=(grd_ymax_n(m_n)-grd_ymin_n(m_n))/(ny_n-1.)
    sz=(grd_zmax_n(m_n)-grd_zmin_n(m_n))/(nz_n-1.)
    !
    !    interpolate between main and sub grids
    !
    !
    !$omp  parallel do
    do k_n=1,nz_n
        az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
        !
        !     find position on main grid
        !
        ak=1.+(az_n-grd_zmin(m))/rz
        k=ak
        kk=k+1
        dz=ak-k
        ddz=1.-dz
        !
        do j_n=1,ny_n
            ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
            aj=1.+(ay_n-grd_ymin(m))/ry
            j=aj
            jj=j+1
            dy=aj-j
            ddy=1.-dy
            !
            do i_n=1,nx_n
                ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                ai=1.+(ax_n-grd_xmin(m))/rx
                i=ai
                ii=i+1
                dx=ai-i
                ddx=1.-dx
                !
                rho_n(i_n,j_n,k_n,m_n)= &
                rho(i,j,k,m)*ddx*ddy*ddz+ &
                rho(i,j,kk,m)*ddx*ddy*dz+ &
                rho(i,jj,k,m)*ddx*dy*ddz+ &
                rho(i,jj,kk,m)*ddx*dy*dz+ &
                rho(ii,j,k,m)*dx*ddy*ddz+ &
                rho(ii,j,kk,m)*dx*ddy*dz+ &
                rho(ii,jj,k,m)*dx*dy*ddz+ &
                rho(ii,jj,kk,m)*dx*dy*dz
                !
                !       rho_n(i_n,j_n,k_n,m_n)= rho(i,j,k,m)
                !
                !       write(7,*)i,j,k,m
                !          rho_max=amax1(rho_max,rho(i,j,k,m))
                !          rho_n_max=amax1(rho_n_max,rho_n(i_n,j_n,k_n,m_n))
            enddo
        enddo
    enddo
    !      write(6,*)'refine',rho_max,rho_n_max
    !
    !      interpolate onto subgrids
    !
    !
    do m_n=ngrd_n-1,1,-1
        m=m_n+1
        rx=(grd_xmax_n(m)-grd_xmin_n(m))/(nx_n-1.)
        ry=(grd_ymax_n(m)-grd_ymin_n(m))/(ny_n-1.)
        rz=(grd_zmax_n(m)-grd_zmin_n(m))/(nz_n-1.)
        sx=(grd_xmax_n(m_n)-grd_xmin_n(m_n))/(nx_n-1.)
        sy=(grd_ymax_n(m_n)-grd_ymin_n(m_n))/(ny_n-1.)
        sz=(grd_zmax_n(m_n)-grd_zmin_n(m_n))/(nz_n-1.)
        !
        !    interpolate between main and sub grids
        !
        do k_n=1,nz_n
            az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
            !
            !     find position on bigger grid
            !
            ak=1.+(az_n-grd_zmin_n(m))/rz
            k=ak
            kk=k+1
            dz=ak-k
            ddz=1.-dz
            !
            do j_n=1,ny_n
                ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
                aj=1.+(ay_n-grd_ymin_n(m))/ry
                j=aj
                jj=j+1
                dy=aj-j
                ddy=1.-dy
                !
                do i_n=1,nx_n
                    ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                    ai=1.+(ax_n-grd_xmin_n(m))/rx
                    i=ai
                    ii=i+1
                    dx=ai-i
                    ddx=1.-dx
                    !
                    rho_n(i_n,j_n,k_n,m_n)= &
                    rho_n(i,j,k,m)*ddx*ddy*ddz+ &
                    rho_n(i,j,kk,m)*ddx*ddy*dz+ &
                    rho_n(i,jj,k,m)*ddx*dy*ddz+ &
                    rho_n(i,jj,kk,m)*ddx*dy*dz+ &
                    rho_n(ii,j,k,m)*dx*ddy*ddz+ &
                    rho_n(ii,j,kk,m)*dx*ddy*dz+ &
                    rho_n(ii,jj,k,m)*dx*dy*ddz+ &
                    rho_n(ii,jj,kk,m)*dx*dy*dz
                    !
                enddo
            enddo
        enddo
    enddo
    !
    return
end

!
!       ************************************************
!
subroutine shift_grd(rho,nx,ny,nz,ngrd,main_grid, &
    rho_n,nx_n,ny_n,nz_n,ngrd_n,m_n,mx,my, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    
    !
    dimension rho(nx,ny,nz,ngrd),rho_n(nx_n,ny_n,nz_n,ngrd_n)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    
    !
    !     move grid so that it can keep up with an orbiting spacecraft
    !
    mx2=mx*2
    my2=my*2
    if (m_n.ne.ngrd_n)then
        !
        m=m_n+1
        rx=(grd_xmax_n(m)-grd_xmin_n(m))/(nx_n-1.)
        ry=(grd_ymax_n(m)-grd_ymin_n(m))/(ny_n-1.)
        rz=(grd_zmax_n(m)-grd_zmin_n(m))/(nz_n-1.)
        sx=(grd_xmax_n(m_n)-grd_xmin_n(m_n))/(nx_n-1.)
        sy=(grd_ymax_n(m_n)-grd_ymin_n(m_n))/(ny_n-1.)
        sz=(grd_zmax_n(m_n)-grd_zmin_n(m_n))/(nz_n-1.)
        !
        !      check for x-shift
        !
        if(mx.ne.0)then
            !
            !      set indices for grid shift
            !
    
            if(mx.gt.0)then      ! xshift
                n1=1
                n2=nx_n-mx2
                isign=1
                n3=n2+isign
                n4=nx_n
            else
                n1=nx_n
                n2=1-mx2
                isign=-1
                n3=n2+isign
                n4=1
            endif
            !
            !       shift array elements
            !
            !       write(6,*)'shifting xgrd',m_n,n1,n2,isign,n3,n4
            !
            !$omp  parallel do
            do k=1,nz_n
                do j=1,ny_n
                    do i=n1,n2,isign
                        ii=i+mx2
                        rho_n(i,j,k,m_n)=rho_n(ii,j,k,m_n)
                    enddo
                enddo
            enddo
            !
            !       set boundary conditions
            !
            !
            !    interpolate between main and sub grids
            !
            !$omp  parallel do
            do k_n=1,nz_n
                az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
                !
                !     find position on bigger grid
                !
                ak=1.+(az_n-grd_zmin_n(m))/rz
                k=ak
                kk=k+1
                dz=ak-k
                ddz=1.-dz
                !
                do j_n=1,ny_n
                    ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
                    aj=1.+(ay_n-grd_ymin_n(m))/ry
                    j=aj
                    jj=j+1
                    dy=aj-j
                    ddy=1.-dy
                    !
                    do i_n=n3,n4,isign
                        ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                        ai=1.+(ax_n-grd_xmin_n(m))/rx
                        i=ai
                        ii=i+1
                        dx=ai-i
                        ddx=1.-dx
                        !
                        rho_n(i_n,j_n,k_n,m_n)= &
                        rho_n(i,j,k,m)*ddx*ddy*ddz+ &
                        rho_n(i,j,kk,m)*ddx*ddy*dz+ &
                        rho_n(i,jj,k,m)*ddx*dy*ddz+ &
                        rho_n(i,jj,kk,m)*ddx*dy*dz+ &
                        rho_n(ii,j,k,m)*dx*ddy*ddz+ &
                        rho_n(ii,j,kk,m)*dx*ddy*dz+ &
                        rho_n(ii,jj,k,m)*dx*dy*ddz+ &
                        rho_n(ii,jj,kk,m)*dx*dy*dz
                        !
                    enddo
                enddo
            enddo
        endif    ! end xshift
        !
        !      check for y-shift
        !
        if(my.ne.0)then
            if(my.gt.0)then      ! yshift
                n1=1
                n2=ny_n-my2
                isign=1
                n3=n2+isign
                n4=ny_n
            else
                n1=ny_n
                n2=1-my2
                isign=-1
                n3=n2+isign
                n4=1
            endif
            !
            !     write(6,*)'y grid shft',m_n,n1,n2,isign,n3,n4
            !
            !       shift array elements
            !
            !$omp  parallel do
            do k=1,nz_n
                do j=n1,n2,isign
                    jj=j+my2
                    do i=1,nx_n
                        rho_n(i,j,k,m_n)=rho_n(i,jj,k,m_n)
                    enddo
                enddo
            enddo
            !
            !       set boundary conditions off bigger grid
            !
            !
            !$omp  parallel do
            do k_n=1,nz_n
                az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
                !
                !     find position on bigger grid
                !
                ak=1.+(az_n-grd_zmin_n(m))/rz
                k=ak
                kk=k+1
                dz=ak-k
                ddz=1.-dz
                !
                do j_n=n3,n4,isign
                    ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
                    aj=1.+(ay_n-grd_ymin_n(m))/ry
                    j=aj
                    jj=j+1
                    dy=aj-j
                    ddy=1.-dy
                    !
                    do i_n=1,nx_n
                        ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                        ai=1.+(ax_n-grd_xmin_n(m))/rx
                        i=ai
                        ii=i+1
                        dx=ai-i
                        ddx=1.-dx
                        !
                        rho_n(i_n,j_n,k_n,m_n)= &
                        rho_n(i,j,k,m)*ddx*ddy*ddz+ &
                        rho_n(i,j,kk,m)*ddx*ddy*dz+ &
                        rho_n(i,jj,k,m)*ddx*dy*ddz+ &
                        rho_n(i,jj,kk,m)*ddx*dy*dz+ &
                        rho_n(ii,j,k,m)*dx*ddy*ddz+ &
                        rho_n(ii,j,kk,m)*dx*ddy*dz+ &
                        rho_n(ii,jj,k,m)*dx*dy*ddz+ &
                        rho_n(ii,jj,kk,m)*dx*dy*dz
                        !
                    enddo
                enddo
            enddo
        endif    ! end yshift
        !
    else  ! ngrd_n shift
        !
        m=main_grid
    
        rx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
        ry=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
        rz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
        sx=(grd_xmax_n(m_n)-grd_xmin_n(m_n))/(nx_n-1.)
        sy=(grd_ymax_n(m_n)-grd_ymin_n(m_n))/(ny_n-1.)
        sz=(grd_zmax_n(m_n)-grd_zmin_n(m_n))/(nz_n-1.)
        !
        if(mx.ne.0)then
            if(mx.gt.0)then      ! xshift
                n1=1
                n2=nx_n-mx2
                isign=1
                n3=n2+isign
                n4=nx_n
            else
                n1=nx_n
                n2=1-mx2
                isign=-1
                n3=n2+isign
                n4=1
            endif
            !
            !       shift array elements
            !
            !       write(6,*)'shifting xgrd',m_n,n1,n2,isign,n3,n4
            !
            !$omp  parallel do
            do k=1,nz_n
                do j=1,ny_n
                    do i=n1,n2,isign
                        ii=i+mx2
                        rho_n(i,j,k,m_n)=rho_n(ii,j,k,m_n)
                    enddo
                enddo
            enddo
            !
            !       set boundary conditions of main grid
            !
            !
            !$omp  parallel do
            do k_n=1,nz_n
                az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
                !
                !       find position on main grid
                !
                ak=1.+(az_n-grd_zmin(m))/rz
                k=ak
                kk=k+1
                dz=ak-k
                ddz=1.-dz
                !
                do j_n=1,ny_n
                    ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
                    aj=1.+(ay_n-grd_ymin(m))/ry
                    j=aj
                    jj=j+1
                    dy=aj-j
                    ddy=1.-dy
                    !
                    do i_n=n3,n4,isign
                        ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                        ai=1.+(ax_n-grd_xmin(m))/rx
                        i=ai
                        ii=i+1
                        dx=ai-i
                        ddx=1.-dx
                        !
                        rho_n(i_n,j_n,k_n,m_n)= &
                        rho(i,j,k,m)*ddx*ddy*ddz+ &
                        rho(i,j,kk,m)*ddx*ddy*dz+ &
                        rho(i,jj,k,m)*ddx*dy*ddz+ &
                        rho(i,jj,kk,m)*ddx*dy*dz+ &
                        rho(ii,j,k,m)*dx*ddy*ddz+ &
                        rho(ii,j,kk,m)*dx*ddy*dz+ &
                        rho(ii,jj,k,m)*dx*dy*ddz+ &
                        rho(ii,jj,kk,m)*dx*dy*dz
    
                    enddo
                enddo
            enddo
        endif    ! xshift
        !
        !     do y-shift
        !
        if(my.ne.0)then
            if(my.gt.0)then      ! yshift
                n1=1
                n2=ny_n-my2
                isign=1
                n3=n2+isign
                n4=ny_n
            else
                n1=ny_n
                n2=1-my2
                isign=-1
                n3=n2+isign
                n4=1
            endif
            !
            !       shift array elements
            !
            !       write(6,*)'shifting ygrd',m_n,n1,n2,isign,n3,n4
            !
            !$omp  parallel do
            do k=1,nz_n
                do j=n1,n2,isign
                    jj=j+my2
                    do i=1,nx_n
                        rho_n(i,j,k,m_n)=rho_n(i,jj,k,m_n)
                    enddo
                enddo
            enddo
            !
            !       set boundary conditions of main grid
            !
            !
            !$omp  parallel do
            do k_n=1,nz_n
                az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
                !
                !       find position on main grid
                !
                ak=1.+(az_n-grd_zmin(m))/rz
                k=ak
                kk=k+1
                dz=ak-k
                ddz=1.-dz
                !
                do j_n=n3,n4,isign
                    ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
                    aj=1.+(ay_n-grd_ymin(m))/ry
                    j=aj
                    jj=j+1
                    dy=aj-j
                    ddy=1.-dy
                    !
                    do i_n=1,nx_n
                        ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                        ai=1.+(ax_n-grd_xmin(m))/rx
                        i=ai
                        ii=i+1
                        dx=ai-i
                        ddx=1.-dx
                        !
                        rho_n(i_n,j_n,k_n,m_n)= &
                        rho(i,j,k,m)*ddx*ddy*ddz+ &
                        rho(i,j,kk,m)*ddx*ddy*dz+ &
                        rho(i,jj,k,m)*ddx*dy*ddz+ &
                        rho(i,jj,kk,m)*ddx*dy*dz+ &
                        rho(ii,j,k,m)*dx*ddy*ddz+ &
                        rho(ii,j,kk,m)*dx*ddy*dz+ &
                        rho(ii,jj,k,m)*dx*dy*ddz+ &
                        rho(ii,jj,kk,m)*dx*dy*dz
    
                    enddo
                enddo
            enddo
        endif    !end y-shift
    
    endif    ! end ngrd_n shift
    !
    return
end

!
!     *********************************
!
subroutine set_rho(qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,rmassq, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,rmassh, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,rmasso, &
    epres,nx,ny,nz,ngrd,m,o_conc)
    !
    !    checks for minium rho and negative pressure
    !     and resets value if necessary
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    !     d_min is the minimum allowable density
    !
    q_min=2.e-4      ! min q aboslute density
    h_min=1.e-4      ! min h aboslute density
    o_min=2.e-5      ! min o aboslute density
    dstep=20.
    !
    !$omp  parallel do
    do  k=1,nz
        do  j=1,ny
            do  i=1,nx
                !
                tden=abs(qrho(i,j,k,m))/rmassq+abs(hrho(i,j,k,m))/rmassh &
                +abs(orho(i,j,k,m))/rmasso+1.e-6
                rhdens=abs(hrho(i,j,k,m))/rmassh/tden
                rodens=abs(orho(i,j,k,m))/rmasso/tden
                !
                if ((qrho(i,j,k,m).lt.q_min).or.(qpresx(i,j,k,m).lt.0.).or. &
                    (qpresy(i,j,k,m).lt.0.).or.(qpresz(i,j,k,m).lt.0.))then
                    !
                    qrho(i,j,k,m)=amax1(abs(qrho(i,j,k,m)),1.05*q_min)
                    apres=(abs(qpresx(i,j,k,m))+abs(qpresy(i,j,k,m))+ &
                    abs(qpresz(i,j,k,m)))/3.
                    qpresx(i,j,k,m)=apres
                    qpresy(i,j,k,m)=apres
                    qpresz(i,j,k,m)=apres
                    qpresxy(i,j,k,m)=0.
                    qpresxz(i,j,k,m)=0.
                    qpresyz(i,j,k,m)=0.
                endif
                !
                if ((hrho(i,j,k,m).lt.h_min).or.(hpresx(i,j,k,m).lt.0.).or. &
                    (hpresy(i,j,k,m).lt.0.).or.(hpresz(i,j,k,m).lt.0.))then
            
                    hrho(i,j,k,m)=amax1(abs(hrho(i,j,k,m)),1.05*h_min)
                    apres=(abs(hpresx(i,j,k,m))+abs(hpresy(i,j,k,m))+ &
                    abs(hpresz(i,j,k,m)))/3.
                    hpresx(i,j,k,m)=apres
                    hpresy(i,j,k,m)=apres
                    hpresz(i,j,k,m)=apres
                    hpresxy(i,j,k,m)=0.
                    hpresxz(i,j,k,m)=0.
                    hpresyz(i,j,k,m)=0.
                endif
                !
                if ((orho(i,j,k,m).lt.o_min).or.(opresx(i,j,k,m).lt.0.).or. &
                    (opresy(i,j,k,m).lt.0.).or.(opresz(i,j,k,m).lt.0.))then
                
                    orho(i,j,k,m)=amax1(abs(orho(i,j,k,m)),1.05*o_min)
                    apres=(abs(opresx(i,j,k,m))+abs(opresy(i,j,k,m))+ &
                    abs(opresz(i,j,k,m)))/3.
                    opresx(i,j,k,m)=apres
                    opresy(i,j,k,m)=apres
                    opresz(i,j,k,m)=apres
                    opresxy(i,j,k,m)=0.
                    opresxz(i,j,k,m)=0.
                    opresyz(i,j,k,m)=0.
                endif
                !
                epres(i,j,k,m)=abs(epres(i,j,k,m))
                !
            enddo
        enddo
    enddo
     
    return
end

!
!      **********************************************
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
    !    checks for minium rho and negative pressure
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
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz)
    !
    common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9), &
    grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9), &
    rx,ry,rz,xdip,ydip,zdip,rearth,b0, &
    sin_tilt,cos_tilt
    !
    !      determine speeds on the code and changes accordingly
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
    !
    !       find bulk velocities : species 1
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
    !       find bulk velocities : species 3
    !
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

!
!      **********************************************
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
    !    checks for minium rho and negative pressure
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
        !
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
        !
        !       find bulk velocities : species 1
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
        !       find bulk velocities : species 3
        !
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
