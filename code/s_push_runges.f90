!
!	This file contains three subroutines:
!	push_bfld
!	push_elec
!	push_ion
!
subroutine push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
    efldx,efldy,efldz,nx,ny,nz,n_grids,box,delt, &
    rx,ry,rz)
    !
    !      standard runge-kutta time step
    !
    integer box
    dimension bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    oldbx(nx,ny,nz,n_grids),oldby(nx,ny,nz,n_grids),oldbz(nx,ny,nz,n_grids), &
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
                bx(i,j,k,box)=oldbx(i,j,k,box)-delt*( &
                (  (efldz(i,jj,k)-efldz(i,jm,k))/dyt ) &
                -( (efldy(i,j,kk)-efldy(i,j,km))/dzt )  )
                !
                by(i,j,k,box)=oldby(i,j,k,box)+delt*( &
                (  (efldz(ii,j,k)-efldz(im,j,k))/dxt ) &
                -( (efldx(i,j,kk)-efldx(i,j,km))/dzt )  )
                !
                bz(i,j,k,box)=oldbz(i,j,k,box)-delt*( &
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
    gamma,gamma1,nx,ny,nz,n_grids,box,delt,rx,ry,rz)
    !
    !      evolves the electron pressure equation
    !
    integer box
    dimension epres(nx,ny,nz,n_grids),oldepres(nx,ny,nz,n_grids), &
    wrkepres(nx,ny,nz,n_grids)
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
                egradp_x=(wrkepres(ii,j,k,box)-wrkepres(im,j,k,box))/dxt
                egradp_y=(wrkepres(i,jj,k,box)-wrkepres(i,jm,k,box))/dyt
                egradp_z=(wrkepres(i,j,kk,box)-wrkepres(i,j,km,box))/dzt
                !
                !       pressure equations:
                !
                epres(i,j,k,box)=oldepres(i,j,k,box)-delt*gamma* &
                ( ( (wrkepres(ii,j,k,box)*evx(ii,j,k) &
                -wrkepres(im,j,k,box)*evx(im,j,k))/dxt ) + &
                ( (wrkepres(i,jj,k,box)*evy(i,jj,k) &
                -wrkepres(i,jm,k,box)*evy(i,jm,k))/dyt ) + &
                ( (wrkepres(i,j,kk,box)*evz(i,j,kk) &
                -wrkepres(i,j,km,box)*evz(i,j,km))/dzt ) )
                epres(i,j,k,box)=epres(i,j,k,box) &
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
    nx,ny,nz,n_grids,box,delt,grav,re_equiv,reynolds, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,isotropic)
    !
    !      standard runge-kutta push for ion equations
    !
    integer box
    dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    qvx(nx,ny,nz),qvy(nx,ny,nz),qvz(nx,ny,nz), &
    tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
    btot(nx,ny,nz)
    !
    dimension qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
    qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
    qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
    qpresyz(nx,ny,nz,n_grids), &
    qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    oldqrho(nx,ny,nz,n_grids),oldqpresx(nx,ny,nz,n_grids), &
    oldqpresy(nx,ny,nz,n_grids),oldqpresz(nx,ny,nz,n_grids), &
    oldqpresxy(nx,ny,nz,n_grids),oldqpresxz(nx,ny,nz,n_grids), &
    oldqpresyz(nx,ny,nz,n_grids), &
    oldqpx(nx,ny,nz,n_grids),oldqpy(nx,ny,nz,n_grids), &
    oldqpz(nx,ny,nz,n_grids), &
    wrkqrho(nx,ny,nz,n_grids),wrkqpresx(nx,ny,nz,n_grids), &
    wrkqpresy(nx,ny,nz,n_grids),wrkqpresz(nx,ny,nz,n_grids), &
    wrkqpresxy(nx,ny,nz,n_grids),wrkqpresxz(nx,ny,nz,n_grids), &
    wrkqpresyz(nx,ny,nz,n_grids), &
    wrkqpx(nx,ny,nz,n_grids),wrkqpy(nx,ny,nz,n_grids), &
    wrkqpz(nx,ny,nz,n_grids)
    !
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    !        set distance scales
    !
    ddx=(grd_xmax(box)-grd_xmin(box))/(nx-1)
    ddy=(grd_ymax(box)-grd_ymin(box))/(ny-1)
    ddz=(grd_zmax(box)-grd_zmin(box))/(nz-1)
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
                aqrho=amax1(wrkqrho(i,j,k,box),d_min)
                !
                qvx(i,j,k)=wrkqpx(i,j,k,box)/aqrho
                qvy(i,j,k)=wrkqpy(i,j,k,box)/aqrho
                qvz(i,j,k)=wrkqpz(i,j,k,box)/aqrho
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
        az=grd_zmin(box)+ddz*(k-1)
    	!
        !        find estimates for the fluid at n+1/2
    	!
        do j=2,ny-1
            jm=j-1
            jj=j+1
            ay=grd_ymin(box)+ddy*(j-1)
            !
            do i=2,nx-1
                ax=grd_xmin(box)+ddx*(i-1)
                ii=i+1
                im=i-1
                !
                !       species 1
                !
                qrho(i,j,k,box)=oldqrho(i,j,k,box)-delt*( &
                ( (wrkqpx(ii,j,k,box)- wrkqpx(im,j,k,box))/dxt) &
                +( (wrkqpy(i,jj,k,box)- wrkqpy(i,jm,k,box))/dyt) &
                +( (wrkqpz(i,j,kk,box)- wrkqpz(i,j,km,box))/dzt) )
                !
                aqpx=oldqpx(i,j,k,box)-delt*( &
                ((wrkqpx(ii,j,k,box)*qvx(ii,j,k)- wrkqpx(im,j,k,box)*qvx(im,j,k)) &
                /dxt ) &
                +((wrkqpx(i,jj,k,box)*qvy(i,jj,k)- wrkqpx(i,jm,k,box)*qvy(i,jm,k)) &
                /dyt ) &
                +((wrkqpx(i,j,kk,box)*qvz(i,j,kk)- wrkqpx(i,j,km,box)*qvz(i,j,km)) &
                /dzt ) )
                !
                aqpy=oldqpy(i,j,k,box)-delt*( &
                ((wrkqpy(ii,j,k,box)*qvx(ii,j,k)- wrkqpy(im,j,k,box)*qvx(im,j,k)) &
                /dxt ) &
                +((wrkqpy(i,jj,k,box)*qvy(i,jj,k)- wrkqpy(i,jm,k,box)*qvy(i,jm,k)) &
                /dyt ) &
                +((wrkqpy(i,j,kk,box)*qvz(i,j,kk)- wrkqpy(i,j,km,box)*qvz(i,j,km)) &
                /dzt ) )
                !
                aqpz=oldqpz(i,j,k,box)-delt*( &
                ((wrkqpz(ii,j,k,box)*qvx(ii,j,k)- wrkqpz(im,j,k,box)*qvx(im,j,k)) &
                /dxt ) &
                +((wrkqpz(i,jj,k,box)*qvy(i,jj,k)- wrkqpz(i,jm,k,box)*qvy(i,jm,k)) &
                /dyt ) &
                +((wrkqpz(i,j,kk,box)*qvz(i,j,kk)- wrkqpz(i,j,km,box)*qvz(i,j,km)) &
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
                qgradp_x=(wrkqpresx(ii,j,k,box)-wrkqpresx(im,j,k,box))/dxt &
                    +(wrkqpresxy(i,jj,k,box)-wrkqpresxy(i,jm,k,box))/dxt &
                    +(wrkqpresxz(i,j,kk,box)-wrkqpresxz(i,j,km,box))/dxt
                qgradp_y=(wrkqpresxy(ii,j,k,box)-wrkqpresxy(im,j,k,box))/dxt &
                    +(wrkqpresy(i,jj,k,box)-wrkqpresy(i,jm,k,box))/dyt &
                    +(wrkqpresyz(i,j,kk,box)-wrkqpresyz(i,j,km,box))/dxt
                qgradp_z=(wrkqpresxz(ii,j,k,box)-wrkqpresxz(im,j,k,box))/dxt &
                    +(wrkqpresyz(i,jj,k,box)-wrkqpresyz(i,jm,k,box))/dyt &
                    +(wrkqpresz(i,j,kk,box)-wrkqpresz(i,j,km,box))/dzt
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
                qden=wrkqrho(i,j,k,box)/rmassq
                qpx(i,j,k,box)=aqpx+delt* &
                ( qden*(dele_x+delv_x) -qgradp_x )
                qpy(i,j,k,box)=aqpy+delt* &
                ( qden*(dele_y+delv_y) -qgradp_y )
                qpz(i,j,k,box)=aqpz+delt* &
                ( qden*(dele_z+delv_z) -qgradp_z )
                !
                !       add in jxb and grad p forces
                !       qjcrossb_x=(cury(i,j,k)*abz-curz(i,j,k)*aby)
                !       qjcrossb_y=-(curx(i,j,k)*abz-curz(i,j,k)*abx)
                !       qjcrossb_z=(curx(i,j,k)*aby-cury(i,j,k)*abx)
                !
                !       qpx(i,j,k,box)=aqpx+delt*(qjcrossb_x-qgradp_x)
                !       qpy(i,j,k,box)=aqpy+delt*(qjcrossb_y-qgradp_y)
                !       qpz(i,j,k,box)=aqpz+delt*(qjcrossb_z-qgradp_z)
                !
                !       add gravity
                !
                radius=sqrt(ax**2+ay**2+az**2)+0.0000001
                g=delt*wrkqrho(i,j,k,box)*grav/(re_equiv*radius)**2
                qpx(i,j,k,box)=qpx(i,j,k,box)-g*ax/radius
                qpy(i,j,k,box)=qpy(i,j,k,box)-g*ay/radius
                qpz(i,j,k,box)=qpz(i,j,k,box)-g*az/radius
                !
                !       pressure equations: isotropic components
                !
                qpresx(i,j,k,box)=oldqpresx(i,j,k,box)-delt* &
                    ( ( (wrkqpresx(ii,j,k,box)*qvx(ii,j,k) &
                    -wrkqpresx(im,j,k,box)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresx(i,jj,k,box)*qvy(i,jj,k) &
                    -wrkqpresx(i,jm,k,box)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresx(i,j,kk,box)*qvz(i,j,kk) &
                    -wrkqpresx(i,j,km,box)*qvz(i,j,km))/dzt ) ) &
                    -2.*delt*wrkqpresx(i,j,k,box) &
                    *(qvx(ii,j,k)-qvx(im,j,k))/dxt
                !
                qpresy(i,j,k,box)=oldqpresy(i,j,k,box)-delt* &
                    ( ( (wrkqpresy(ii,j,k,box)*qvx(ii,j,k) &
                    -wrkqpresy(im,j,k,box)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresy(i,jj,k,box)*qvy(i,jj,k) &
                    -wrkqpresy(i,jm,k,box)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresy(i,j,kk,box)*qvz(i,j,kk) &
                    -wrkqpresy(i,j,km,box)*qvz(i,j,km))/dzt ) ) &
                    -2.*delt*wrkqpresy(i,j,k,box) &
                    * (qvy(i,jj,k)-qvy(i,jm,k))/dyt
                !
                qpresz(i,j,k,box)=oldqpresz(i,j,k,box)-delt* &
                    ( ( (wrkqpresz(ii,j,k,box)*qvx(ii,j,k) &
                    -wrkqpresz(im,j,k,box)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresz(i,jj,k,box)*qvy(i,jj,k) &
                    -wrkqpresz(i,jm,k,box)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresz(i,j,kk,box)*qvz(i,j,kk) &
                    -wrkqpresz(i,j,km,box)*qvz(i,j,km))/dzt ) ) &
                    -2.*delt*wrkqpresz(i,j,k,box) &
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
                qpresx(i,j,k,box)=qpresx(i,j,k,box) &
                    -2.*delt*wrkqpresxy(i,j,k,box) &
                    *(qvx(i,jj,k)-qvx(i,jm,k))/dyt &
                    -2.*delt*wrkqpresxz(i,j,k,box) &
                    *(qvx(i,j,kk)-qvx(i,j,km))/dzt &
                    +2.*delt*skin_factor/rmass*( &
                    wrkqpresxy(i,j,k,box)*abz &
                    -wrkqpresxz(i,j,k,box)*aby )
        		!
                qpresy(i,j,k,box)=qpresy(i,j,k,box) &
                    -2.*delt*wrkqpresxy(i,j,k,box) &
                    * (qvy(ii,j,k)-qvy(im,j,k))/dxt &
                    -2.*delt*wrkqpresyz(i,j,k,box) &
                    * (qvy(i,j,kk)-qvy(i,j,km))/dzt &
                    +2.*delt*skin_factor/rmass*( &
                    wrkqpresyz(i,j,k,box)*abx &
                    -wrkqpresxy(i,j,k,box)*abz )
                !
                qpresz(i,j,k,box)=qpresz(i,j,k,box) &
                    -2.*delt*wrkqpresxz(i,j,k,box) &
                    *(qvz(ii,j,k)-qvz(im,j,k))/dxt &
                    -2.*delt*wrkqpresyz(i,j,k,box) &
                    *(qvz(i,jj,k)-qvz(i,jm,k))/dyt &
                    +2.*delt*skin_factor/rmass*( &
                    wrkqpresxz(i,j,k,box)*aby &
                    -wrkqpresyz(i,j,k,box)*abx )
                !
                !       pressure equations: offdiagonal elements
                !
                qpresxy(i,j,k,box)=oldqpresxy(i,j,k,box)-delt* &
                    ( ( (wrkqpresxy(ii,j,k,box)*qvx(ii,j,k) &
                    -wrkqpresxy(im,j,k,box)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresxy(i,jj,k,box)*qvy(i,jj,k) &
                    -wrkqpresxy(i,jm,k,box)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresxy(i,j,kk,box)*qvz(i,j,kk) &
                    -wrkqpresxy(i,j,km,box)*qvz(i,j,km))/dzt ) )
				!
                qpresxy(i,j,k,box)=qpresxy(i,j,k,box) &
                    -delt*wrkqpresxy(i,j,k,box) &
                    *( (qvx(ii,j,k)-qvx(im,j,k))/dxt &
                    +(qvy(i,jj,k)-qvy(i,jm,k))/dyt ) &
                    -delt*wrkqpresy(i,j,k,box) &
                    *(qvx(i,jj,k)-qvx(i,jm,k))/dyt &
                    -delt*wrkqpresx(i,j,k,box) &
                    *(qvy(ii,j,k)-qvy(im,j,k))/dxt &
                    -delt*wrkqpresyz(i,j,k,box) &
                    *(qvx(i,j,kk)-qvx(i,j,km))/dzt &
                    -delt*wrkqpresxz(i,j,k,box) &
                    *(qvy(i,j,kk)-qvy(i,j,km))/dzt &
                    +delt*skin_factor/rmass*( &
                    abz*(wrkqpresy(i,j,k,box)-wrkqpresx(i,j,k,box)) &
                    -aby*wrkqpresyz(i,j,k,box)+abx*wrkqpresxz(i,j,k,box))
                !
                qpresxz(i,j,k,box)=oldqpresxz(i,j,k,box)-delt* &
                    ( ( (wrkqpresxz(ii,j,k,box)*qvx(ii,j,k) &
                    -wrkqpresxz(im,j,k,box)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresxz(i,jj,k,box)*qvy(i,jj,k) &
                    -wrkqpresxz(i,jm,k,box)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresxz(i,j,kk,box)*qvz(i,j,kk) &
                    -wrkqpresxz(i,j,km,box)*qvz(i,j,km))/dzt ) )
                qpresxz(i,j,k,box)=qpresxz(i,j,k,box) &
                    -delt*wrkqpresxz(i,j,k,box) &
                    *( (qvx(ii,j,k)-qvx(im,j,k))/dxt &
                    +(qvz(i,j,kk)-qvz(i,j,km))/dzt ) &
                    -delt*wrkqpresz(i,j,k,box) &
                    *(qvx(i,j,kk)-qvx(i,j,km))/dzt &
                    -delt*wrkqpresx(i,j,k,box) &
                    *(qvz(ii,j,k)-qvz(im,j,k))/dxt &
                    -delt*wrkqpresyz(i,j,k,box) &
                    *(qvx(i,jj,k)-qvx(i,jm,k))/dyt &
                    -delt*wrkqpresxy(i,j,k,box) &
                    *(qvz(i,jj,k)-qvz(i,jm,k))/dyt &
                    +delt*skin_factor/rmass*( &
                    aby*(wrkqpresx(i,j,k,box)-wrkqpresz(i,j,k,box)) &
                    -abx*wrkqpresxy(i,j,k,box)+abz*wrkqpresyz(i,j,k,box))
                !
                qpresyz(i,j,k,box)=oldqpresyz(i,j,k,box)-delt* &
                    ( ( (wrkqpresyz(ii,j,k,box)*qvx(ii,j,k) &
                    -wrkqpresyz(im,j,k,box)*qvx(im,j,k))/dxt ) + &
                    ( (wrkqpresyz(i,jj,k,box)*qvy(i,jj,k) &
                    -wrkqpresyz(i,jm,k,box)*qvy(i,jm,k))/dyt ) + &
                    ( (wrkqpresyz(i,j,kk,box)*qvz(i,j,kk) &
                    -wrkqpresyz(i,j,km,box)*qvz(i,j,km))/dzt ) )
                qpresyz(i,j,k,box)=qpresyz(i,j,k,box) &
                    -delt*wrkqpresyz(i,j,k,box) &
                    *( (qvy(i,jj,k)-qvy(i,jm,k))/dyt &
                    +(qvz(i,j,kk)-qvz(i,j,km))/dzt ) &
                    -delt*wrkqpresy(i,j,k,box) &
                    *(qvz(i,jj,k)-qvz(i,jm,k))/dyt &
                    -delt*wrkqpresz(i,j,k,box) &
                    *(qvy(i,j,kk)-qvy(i,j,km))/dzt &
                    -delt*wrkqpresxy(i,j,k,box) &
                    *(qvz(ii,j,k)-qvz(im,j,k))/dxt &
                    -delt*wrkqpresxz(i,j,k,box) &
                    *(qvy(ii,j,k)-qvy(im,j,k))/dxt &
                    +delt*skin_factor/rmass*( &
                    abx*(wrkqpresz(i,j,k,box)-wrkqpresy(i,j,k,box)) &
                    -abz*wrkqpresxz(i,j,k,box)+aby*wrkqpresxy(i,j,k,box))
                !
                !
                !
            enddo             ! k loop
        enddo             ! j loop
        !
    enddo             ! i loop
    !
    return
end
