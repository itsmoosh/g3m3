!
!	This file contains three subroutines:
!	push_bfld
!	push_elec
!	push_ion
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
    ! parallelizes loop. RW, aug. 17, 2004
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
!
!	****************************************
!
!
subroutine push_elec(epres,oldepres,wrkepres,evx,evy,evz, &
    gamma,gamma1,nx,ny,nz,ngrd,m,delt,rx,ry,rz)
    !
    !      evolves the electron pressure equation
    !
    dimension epres(nx,ny,nz,ngrd),oldepres(nx,ny,nz,ngrd), &
    wrkepres(nx,ny,nz,ngrd)
    !
    dimension evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
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
    	!
        do j=2,ny-1
            jm=j-1
            jj=j+1
    		!
            do i=2,nx-1
                ii=i+1
                im=i-1
                egradp_x=(wrkepres(ii,j,k,m)-wrkepres(im,j,k,m))/dxt
                egradp_y=(wrkepres(i,jj,k,m)-wrkepres(i,jm,k,m))/dyt
                egradp_z=(wrkepres(i,j,kk,m)-wrkepres(i,j,km,m))/dzt
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
!
!	****************************************
!
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
    !           that can be resolved
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
    !           that can be resolved
	!
    rmass=rmassq
    bmax=0.20*rmassq/delt/reynolds
    !
    ! parallelizes loop. RW, aug. 17, 2004
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
    !
    ! parallelizes loop. BU, oct. 10 2002
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kk=k+1
        az=grd_zmin(m)+ddz*(k-1)
    	!
        !        find estimates for the fluid at n+1/2
    	!
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
                !
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
        		!
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
				!
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
    !
    if(ani.lt.1.00.and.(.not.isotropic))then
        !
        !$omp  parallel do
        do k=2,nz-1
            do j=2,ny-1
                do i=2,nx-1
                    !
                    apres=(qpresx(i,j,k,m)+qpresy(i,j,k,m)+qpresz(i,j,k,m))/3.0
    				!
                    qpresx(i,j,k,m)=ani*qpresx(i,j,k,m)+(1.-ani)*apres
                    qpresy(i,j,k,m)=ani*qpresy(i,j,k,m)+(1.-ani)*apres
                    qpresz(i,j,k,m)=ani*qpresz(i,j,k,m)+(1.-ani)*apres
                    qpresxy(i,j,k,m)=ani*qpresxy(i,j,k,m)
                    qpresxz(i,j,k,m)=ani*qpresxz(i,j,k,m)
                    qpresyz(i,j,k,m)=ani*qpresyz(i,j,k,m)
                    !
                enddo             ! k loop
            enddo             ! j loop
        enddo             ! i loop
        !
    endif
    return
end
