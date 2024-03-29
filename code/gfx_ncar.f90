!
!	This file contains all subroutines related to ncar graphics.
!
subroutine visual( &
    qrho,qpresx,qpresy,qpresz,qpresxy, &
    qpresxz,qpresyz,qpx,qpy,qpz,rmassq, &
    hrho,hpresx,hpresy,hpresz,hpresxy, &
    hpresxz,hpresyz,hpx,hpy,hpz,rmassh, &
    orho,opresx,opresy,opresz,opresxy, &
    opresxz,opresyz,opx,opy,opz,rmasso, &
    epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
    curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
    tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
    nx,ny,nz,n_grids,xspac, &
    cross,along,flat,xcraft,ncraft,re_equiv, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,ut,b_equiv,ti_te,rho_equiv)
    !
	integer, parameter :: dp = kind(1.d0)
    integer box
	real(dp) ut
    real grd_xmin(n_grids),grd_xmax(n_grids),grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids),xspac(n_grids)
    real bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
    qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
    qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
    qpresyz(nx,ny,nz,n_grids), &
    opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
    opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
    opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
    opresyz(nx,ny,nz,n_grids), &
    hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
    hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
    hpresxy(nx,ny,nz,n_grids),hpresxz(nx,ny,nz,n_grids), &
    hpresyz(nx,ny,nz,n_grids), &
    epres(nx,ny,nz,n_grids)
    real bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids),bz0(nx,ny,nz,n_grids)
    real efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz)
    real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz), &
    tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2), &
    cross(ny,nz),along(nx,nz),flat(nx,ny)
    real xcraft(4,ncraft)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    character*4 wd1
    character*12 label
    !
    logical add_dip
    !
    do box=n_grids,1,-1
        rx=xspac(box)
        ry=xspac(box)
        rz=xspac(box)
        !
        ymin=grd_ymin(box)
        ymax=grd_ymax(box)
        zmin=grd_zmin(box)
        zmax=grd_zmax(box)
        xmin=grd_xmin(box)
        xmax=grd_xmax(box)
        !
        xcut=((xmin+xmax))/2.
        !
        add_dip=.false.
        !
		bsx(:,:,:) = bx0(:,:,:,box) + bx(:,:,:,box)
		bsy(:,:,:) = by0(:,:,:,box) + by(:,:,:,box)
		bsz(:,:,:) = bz0(:,:,:,box) + bz(:,:,:,box)
        !
		curx(:,:,:) = 0.
		cury(:,:,:) = 0.
		curz(:,:,:) = 0.
        !
        if(box.le.3)then
            preslim=20.0/float(box)
            po=preslim*0.33
        else
            preslim=20.0/(box-1.)
            po=preslim*0.33
        endif
        !
   		curx(:,:,:) = qpx(:,:,:,box) / amax1( qrho(:,:,:,box), smallbit )
		cury(:,:,:) = qpy(:,:,:,box) / amax1( qrho(:,:,:,box), smallbit )
		curz(:,:,:) = qpz(:,:,:,box) / amax1( qrho(:,:,:,box), smallbit )
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=qpresx(i,j,k,box)
                    presy=qpresy(i,j,k,box)
                    presz=qpresz(i,j,k,box)
                    presxy=qpresxy(i,j,k,box)
                    presxz=qpresxz(i,j,k,box)
                    presyz=qpresyz(i,j,k,box)
                    presmag=sqrt(presx**2+presy**2+presz**2+ &
                        2*(presxy**2)+2*(presyz**2)+2*(presxz**2))+1.e-11
                    apres=(presx+presy+presz+2*presxy+ &
                        2*presxz+2*presyz)/9.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    p_para=sqrt((presx*abx)**2+(presxy*aby)**2+(presxz*abz)**2 + &
                        (presxy*abx)**2+(presy*aby)**2+(presyz*abz)**2 + &
                        (presxz*abx)**2+(presyz*aby)**2+(presz*abz)**2)/(bmag)
					!
                    p_cross=sqrt((presx*vcrossb_x)**2+(presxy*vcrossb_y)**2+(presxz*vcrossb_z)**2+ &
                        (presxy*vcrossb_x)**2+(presy*vcrossb_y)**2+(presyz*vcrossb_z)**2+ &
                        (presxz*vcrossb_x)**2+(presyz*vcrossb_y)**2+(presz*vcrossb_z)**2)/(vcmag)
					!
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
    				!
                    tvx(i,j,k)=sqrt(apres)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=p_para/apres+0.0001
                    efldy(i,j,k)=p_cross/apres+0.0001
                    efldz(i,j,k)=p_perp/apres+0.0001
                enddo
            enddo
        enddo
        !
        wd1=''
        label=''
        write(wd1,'(i1)')box
        label='qpres '//wd1
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    		!
        label='qpara '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='qcross '//wd1
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='qperp '//wd1
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
    	!
		curx(:,:,:) = hpx(:,:,:,box) / amax1( hrho(:,:,:,box), smallbit )
		cury(:,:,:) = hpy(:,:,:,box) / amax1( hrho(:,:,:,box), smallbit )
		curz(:,:,:) = hpz(:,:,:,box) / amax1( hrho(:,:,:,box), smallbit )
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=hpresx(i,j,k,box)
                    presy=hpresy(i,j,k,box)
                    presz=hpresz(i,j,k,box)
                    presxy=hpresxy(i,j,k,box)
                    presxz=hpresxz(i,j,k,box)
                    presyz=hpresyz(i,j,k,box)
                    presmag=sqrt(presx**2+presy**2+presz**2+ &
                        2*(presxy**2)+2*(presxz**2)+2*(presyz**2))+1.e-11
                    apres=(presx+presy+presz+2*presxy+ &
                        2*presxz+2*presyz)/9.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    p_para=sqrt((presx*abx)**2+(presxy*aby)**2+(presxz*abz)**2 + &
                        (presxy*abx)**2+(presy*aby)**2+(presyz*abz)**2 + &
                        (presxz*abx)**2+(presyz*aby)**2+(presz*abz)**2)/(bmag)
					!
                    p_cross=sqrt((presx*vcrossb_x)**2+(presxy*vcrossb_y)**2+(presxz*vcrossb_z)**2+ &
                        (presxy*vcrossb_x)**2+(presy*vcrossb_y)**2+(presyz*vcrossb_z)**2+ &
                        (presxz*vcrossb_x)**2+(presyz*vcrossb_y)**2+(presz*vcrossb_z)**2)/(vcmag)

                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
    				!
                    tvx(i,j,k)=sqrt(apres)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=p_para/apres+0.0001
                    efldy(i,j,k)=p_cross/apres+0.0001
                    efldz(i,j,k)=p_perp/apres+0.0001
                enddo
            enddo
        enddo
        !
        label='hpres '//wd1
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        label='h_para '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)

        label='h_cross '//wd1
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)

        label='h_perp '//wd1
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
    	!
		curx(:,:,:) = opx(:,:,:,box) / amax1( orho(:,:,:,box), smallbit )
		cury(:,:,:) = opy(:,:,:,box) / amax1( orho(:,:,:,box), smallbit )
		curz(:,:,:) = opz(:,:,:,box) / amax1( orho(:,:,:,box), smallbit )
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=opresx(i,j,k,box)
                    presy=opresy(i,j,k,box)
                    presz=opresz(i,j,k,box)
                    presxy=opresxy(i,j,k,box)
                    presxz=opresxz(i,j,k,box)
                    presyz=opresyz(i,j,k,box)
                    presmag=sqrt(presx**2+presy**2+presz**2+ &
                        2*(presxy**2)+2*(presxz**2)+2*(presyz**2))+1.e-11
                    apres=(presx+presy+presz+2*presxy+ &
                        2*presxz+2*presyz)/9.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    p_para=sqrt((presx*abx)**2+(presxy*aby)**2+(presxz*abz)**2 + &
                        (presxy*abx)**2+(presy*aby)**2+(presyz*abz)**2 + &
                        (presxz*abx)**2+(presyz*aby)**2+(presz*abz)**2)/(bmag)
					!	
                    p_cross=sqrt((presx*vcrossb_x)**2+(presxy*vcrossb_y)**2+(presxz*vcrossb_z)**2+ &
                        (presxy*vcrossb_x)**2+(presy*vcrossb_y)**2+(presyz*vcrossb_z)**2+ &
                        (presxz*vcrossb_x)**2+(presyz*vcrossb_y)**2+(presz*vcrossb_z)**2)/(vcmag)
					!
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
    				!
                    tvx(i,j,k)=sqrt(apres)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=p_para/apres+0.0001
                    efldy(i,j,k)=p_cross/apres+0.0001
                    efldz(i,j,k)=p_perp/apres+0.0001
					!
                enddo
            enddo
        enddo
        !
        label='opres '//wd1
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='o_para '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='o_cross '//wd1
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='o_perp '//wd1
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=sqrt(epres(i,j,k,box))
                enddo
            enddo
        enddo
        !
        label='epres '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=abs(qpresx(i,j,k,box)+hpresx(i,j,k,box) &
                        +opresx(i,j,k,box))/3.
                    efldy(i,j,k)=abs(qpresy(i,j,k,box)+hpresy(i,j,k,box) &
                        +opresy(i,j,k,box))/3.
                    efldz(i,j,k)=abs(qpresz(i,j,k,box)+hpresz(i,j,k,box) &
                        +opresz(i,j,k,box))/3.
                    tvx(i,j,k)=sqrt(efldx(i,j,k)+efldy(i,j,k)+efldz(i,j,k) &
                        +epres(i,j,k,box))
                enddo
            enddo
        enddo
		!
        label='tpres '//wd1
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim*3., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !      trace velocity streams
        !
        call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,n_grids,box, &
            rmassq,rmassh,rmasso)
        !
        label='pres-vel '//wd1
        call conflow(tvx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,11,1,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !     find total magnetic field
        !
		bsx(:,:,:) = bx0(:,:,:,box) + bx(:,:,:,box)
		bsy(:,:,:) = by0(:,:,:,box) + by(:,:,:,box)
		bsz(:,:,:) = bz0(:,:,:,box) + bz(:,:,:,box)
        !
		curx(:,:,:) = 0.
		cury(:,:,:) = 0.
		curz(:,:,:) = 0.
        !
    	!
        label='box '//wd1
        call concross(tvx,curx,cury,curz,bsx,bsy,bsz, &
            nx,ny,nz,1,1,box,xcraft,ncraft,re_equiv,r_inner, &
            xmin,xmax,ymin,ymax,zmin,zmax, &
            ut,label,3,11,2.0,add_dip,1,-2,start, &
            tx,ty,tz,tt,tg1,tg2,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    		!
        call contop(tvx,curx,cury,curz,bsx,bsy,bsz, &
            nx,ny,nz,1,1,box,xcraft,ncraft,re_equiv,r_inner, &
            xmin,xmax,ymin,ymax,zmin,zmax, &
            ut,label,3,11,2.0,add_dip,1,0,start, &
            tx,ty,tz,tt,tg1,tg2,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       magnitude of b
        !
        blim=0.
    	!
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldx(i,j,k)=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                    +bsz(i,j,k)**2)
                    efldx(i,j,k)=alog10(b_equiv*efldx(i,j,k)+1.e-10)
                    blim=amax1(blim,efldx(i,j,k))
                enddo
            enddo
        enddo
        !
        label='bmag '//wd1
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,0.1,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        blim=0.001
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldx(i,j,k)=bsz(i,j,k)
                    if (efldx(i,j,k).gt.blim)efldx(i,j,k)=blim
                    if (efldx(i,j,k).lt.-blim)efldx(i,j,k)=-blim
                    efldx(i,j,k)=efldx(i,j,k)+1.001*blim
                enddo
            enddo
        enddo
		!
        label='bz '//wd1 
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,12,1,2.0,2.*blim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !      calculate alfven mach number
        !
        call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,n_grids,box, &
            rmassq,rmassh,rmasso)
        !
        dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
        dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
        dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    aden=qrho(i,j,k,box)+hrho(i,j,k,box)+orho(i,j,k,box)+1.e-8
                    efldx(i,j,k)=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                        +bsz(i,j,k)**2)
                    alf_spd=(efldx(i,j,k)+1.e-2)/sqrt(aden)
                    apres=tvx(i,j,k)+1.e-10
                    cs_spd=sqrt(apres/aden)
                    aspd=sqrt(curx(i,j,k)**2+cury(i,j,k)**2 &
                        +curz(i,j,k)**2)
                    efldy(i,j,k)=aspd/alf_spd
                    efldz(i,j,k)=aspd/cs_spd
                    !
                    !       calculate rotation velocity
                    !
                    ay=grd_ymin(box)+dy*(j-1)-ydip
                    ax=grd_xmin(box)+dx*(i-1)-xdip
                    ar=sqrt(ax**2+ay**2)+1.e-8
    				!
                    rspd=abs( (curx(i,j,k)*ay-cury(i,j,k)*ax) / ar)
                    radius=ar*re_equiv
                    rot_spd=radius*v_rot+.001
                    efldx(i,j,k)=rspd/rot_spd
                enddo
            enddo
        enddo
        !
        !
		curx(:,:,:) = 0.
		cury(:,:,:) = 0.
		curz(:,:,:) = 0.
        !
        label='alf_mach '//wd1 
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        label='cs_mach '//wd1 
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        label='rot_mach '//wd1 
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !      plot individual temperatures
        !
        tempx=0.5
        temph=2.*tempx
        tempo=4.*tempx
        if(box.gt.1) then
            rho_lim=10.
        else
            rho_lim=5.0
        endif
        rfrac=0.5
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    qden=qrho(i,j,k,box)/rmassq
                    if(qden.gt.0.001)then
                        apres=(qpresx(i,j,k,box)+qpresy(i,j,k,box) &
                            +qpresz(i,j,k,box))/3.
                        bsx(i,j,k)=amin1(tempx,sqrt(apres/qden))
                    else
                        bsx(i,j,k)=0.
                    endif
                    hden=hrho(i,j,k,box)/rmassh
                    if(hden.gt.0.0005)then
                        apres=(hpresx(i,j,k,box)+hpresy(i,j,k,box) &
                            +hpresz(i,j,k,box))/3.
                        bsy(i,j,k)=amin1(temph,sqrt(apres/hden))
                    else
                        bsy(i,j,k)=0.
                    endif
                    !
                    oden=orho(i,j,k,box)/rmasso
                    if(oden.gt.0.00002)then
                        apres=(opresx(i,j,k,box)+opresy(i,j,k,box) &
                            +opresz(i,j,k,box))/3.
                        bsz(i,j,k)=amin1(tempo,sqrt(apres/oden))
                    else
                        bsz(i,j,k)=0.
                    endif
                    !
                    efldx(i,j,k)=sqrt(epres(i,j,k,box)/(qden+hden+oden))
                    !
                enddo
            enddo
        enddo
        !
        label='q temp '//wd1 
        call conhot(bsx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,tempx, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='h temp '//wd1 
        call conhot(bsy,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,temph, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='o temp '//wd1 
        call conhot(bsz,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,tempo, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='e temp '//wd1 
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,tempx/sqrt(ti_te), &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : solar wind
        !
		curx(:,:,:) = 0.
		cury(:,:,:) = 0.
		curz(:,:,:) = 0.
        !       call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,n_grids,box)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=(qrho(i,j,k,box)/rmassq)
                    efldx(i,j,k)=alog10(efldx(i,j,k)*rho_equiv)+6. ! per box**3
                enddo
            enddo
        enddo
        !
        if(box.le.3)then
            plot_min=5.
            plot_max=8.5
        else
            plot_min=4.-0.5*(box-3)
            plot_max=7.5-0.5*(box-3)
        endif
        !
        label='q den '//wd1 
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : ionospheric h
        !
        !      call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,n_grids,box)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldy(i,j,k)=hrho(i,j,k,box)/rmassh
                    efldy(i,j,k)=alog10(efldy(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        label='h den '//wd1 
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,box,&
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows: ionospheric o
        !
        !       call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,n_grids,box)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldz(i,j,k)=orho(i,j,k,box)/rmasso
                    efldz(i,j,k)=alog10(efldz(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        label='o den '//wd1 
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows: total
        !
        !      call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
        !    +       opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,n_grids,box,
        !    +      rmassq,rmassh,rmasso)
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,box)/rmasso+ &
                        hrho(i,j,k,box)/rmassh+qrho(i,j,k,box)/rmassq)
                    efldx(i,j,k)=alog10(tden*rho_equiv)+6.
                    if(efldx(i,j,k).lt.0.003)then
                        curx(i,j,k)=0.
                        cury(i,j,k)=0.
                        curz(i,j,k)=0.
                    endif
                enddo
            enddo
        enddo
        !
        label='t den '//wd1 
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,box)/rmasso+ &
                    hrho(i,j,k,box)/rmassh+qrho(i,j,k,box)/rmassq)
                    efldy(i,j,k)=(hrho(i,j,k,box)/rmassh)/(tden+1.e-6)
                    efldz(i,j,k)=(orho(i,j,k,box)/rmasso)/(tden+1.e-6)
                    efldx(i,j,k)=(qrho(i,j,k,box)/rmassq)/(tden+1.e-6)
                enddo
            enddo
        enddo
        !
        label='rqdens x2 '//wd1 
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,0.5, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='rhdens '//wd1 
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,1.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='rodens x2 '//wd1 
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,box, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,0.5, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine arwxy(dx,dy,al,beta,adx1,ady1,adx2,ady2)
    !
    !     calculate relative positions of two points to draw arrows
    !        input: dx,dy  -- vector
    !           al - length  of arrow
    !           beta - angle of arrow
    !     output: (adx1,ady1), (adx2,ady2)  -- two relative positions
    !
    data pih,pi,pi32,pi2/1.57079,3.14159,4.71239,6.28318/
    !
    !       determine angle of the vector
    !
    if(abs(dx).le.0.00001) then
        theta=pih
        if(dy.lt.0.)  theta=pi32
    else if (abs(dy).lt.0.00001) then
        theta=0.
        if(dx.lt.0.) theta=pi
    else
        theta=atan(abs(dy/dx))
        if(dx.lt.0..and.dy.gt.0.) then
            theta=pi-theta
        else if(dx.lt.0..and.dy.lt.0.) then
            theta=pi+theta
        else if(dx.gt.0..and.dy.lt.0.) then
            theta=pi2-theta
        endif
    endif
    !
    !      calculate the relative position of two points
    !
    alfa1=theta-beta
    alfa2=pih-theta-beta
    adx1=-al*cos(alfa1)
    ady1=-al*sin(alfa1)
    adx2=-al*sin(alfa2)
    ady2=-al*cos(alfa2)
    !
    return
end
!
!
!	****************************************
!
!
subroutine aurora(stuff,nx,ny,nz,box,radstrt,re_equiv,iside, &
    time, save_dat,add_two,label,ncon,write_dat,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !
    integer box
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension stuff(nx,ny,nz)
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_two,write_dat
    !
    !     for no line labeling set ilab  to zero
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    call isoclrs
    ilab=0
    ioffm=1
    irecmj=1
    irectx=1
    irecmn=1
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    !      theta_range=0.698132  ! 40 degrees
    theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    do j=1,my
        do i=1,mx
            !
            !         find equivalent latitude
            !
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            dlat=amax1(dlat,0.0000001)
            alat=dlat*del_theta*theta_equiv
            if(alat.le.1.55) then
                cost=cos(alat)
                sint=sin(alat)
                !
                !         find equivalent longitude
                !
                cosp=(i-mx2)/dlat
                sinp=(j-my2)/dlat
                !
                !          find position on grid
                !
                z1=iside*(radstrt*cost)
                x1=radstrt*sint*cosp
                !
                !       dipole space to real space
                !
                ax=x1*cos_tilt+z1*sin_tilt+xdip
                az=-x1*sin_tilt+z1*cos_tilt+zdip
                ay=(radstrt*sint*sinp+ydip)
                !
                !          interpolate data to grid point
                !
                !
                ak=1.+(az-grd_zmin(box))/delz
                k1=ak
                k2=k1+1
                dz=ak-k1
                !
                aj=1.+(ay-grd_ymin(box))/dely
                j1=aj
                j2=j1+1
                dy=aj-j1
                !
                ai=1.+(ax-grd_xmin(box))/delx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                acur=stuff(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +stuff(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +stuff(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +stuff(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +stuff(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +stuff(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +stuff(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +stuff(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                work(i,j)=acur
            endif
            !
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
    title=' : aurora'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
        enddo
    enddo
    !
    if (fmin.ge.0.0)then
        cmin=1.05*fmin
    else
        cmin=0.95*fmin
    endif
    !
    if(fmax.ge.0.0) then
        cmax=0.95*fmax
    else
        cmax=1.05*fmax
    endif
    finc=(cmax-cmin)/float(ncon-1)
    write(title,'(f7.4)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    write(title,'(f7.4)')fmax
    title='fmax'//title
    call wtstr(.25,.98,title,1,0,0)
    call gsplci(1)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
    !
    !     make constant latitude circles
    !
    dd=degrees/float(mx2-1)
    fmax=0.
    do j=1,my
        do i=1,mx
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            clat(i,j)=90.-dlat*dd
            fmax=amax1(fmax,clat(i,j))
        enddo
    enddo
    !
    !     ilab=1
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    !
    !     save data if required
    !
    if(save_dat)then
        if(iside.gt.0) then
            do j=1,my
                do i=1,mx
                    sdata(i,j)=work(i,j)
                enddo
            enddo
        else
            do j=1,my
                do i=1,mx
                    tdata(i,j)=work(i,j)
                enddo
            enddo
        endif
    endif
    !
    !      output data to dat file if necessary
    !
    if(write_dat)then
        write(recdt_f,*) time, work
    endif
    !
    !
    if(.not.add_two)return
    !
    ilab=0
    if(iside.gt.0)then
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+sdata(i,j)
            enddo
        enddo
    else
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+tdata(i,j)
            enddo
        enddo
    endif
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : total aur'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
        enddo
    enddo
    !
    if (fmin.ge.0.0)then
        cmin=1.05*fmin
    else
        cmin=0.95*fmin
    endif
    !
    if(fmax.ge.0.0) then
        cmax=0.95*fmax
    else
        cmax=1.05*fmax
    endif
    finc=(cmax-cmin)/float(ncon-1)
    write(title,'(f6.3)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    write(title,'(f6.3)')fmax
    title='fmax'//title
    call wtstr(.25,.98,title,1,0,0)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
    !
    !     ilab=1
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    return
end
!
!
!	****************************************
!
!
subroutine aurora_bfld(bx,by,bz,nx,ny,nz,box,rx, &
    xmin,xmax,ymin,ymax,zmin,zmax,iside, &
    add_dip,radstrt,re_equiv,r_inner,time,label,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !       this subroutine seeks to determine open field
    !             line positions
    !
    common /space/work(91,91),clat(91,91), &
    xray(1000),yray(1000),zray(1000)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical add_dip,roc
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !      dimension for plotted array
    !
    maxpts=1000
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.6108652
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    rx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
        enddo
    enddo
    !
    do j=1,my
        do i=1,mx
            !
            closed=0.
            !
            !         find equivalent latitude
            !
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            dlat=amax1(dlat,0.0000001)
            alat=dlat*del_theta*theta_equiv
            if(alat.le.1.55) then
                cost=cos(alat)
                sint=sin(alat)
                !
                !         find equivalent longitude
                !
                cosp=(i-mx2)/dlat
                sinp=(j-my2)/dlat
                !
                !          find position on grid
                !
                z1=iside*(radstrt*cost)
                x1=radstrt*sint*cosp
                !
                !       dipole space to real space
                !
                ax=x1*cos_tilt+z1*sin_tilt+xdip
                az=-x1*sin_tilt+z1*cos_tilt+zdip
                ay=(radstrt*sint*sinp+ydip)
                !
                dir=-rx
                call rungem(bx,by,bz,nx,ny,nz,box,rx, &
                ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax,r_inner, &
                add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                !
                !        find the radial distance of end point
                !
                ar=sqrt((xray(npts)-xdip)**2+(yray(npts)-ydip)**2 &
                +(zray(npts)-zdip)**2)
                if(ar.le.r_inner+2.*rx)closed=closed+1.
                !
                dir=rx
                call rungem(bx,by,bz,nx,ny,nz,box,rx, &
                ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax,r_inner, &
                add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                !
                !        find the radial distance of end point
                !
                ar=sqrt(xray(npts) **2+yray(npts) **2+zray(npts)**2)
                if(ar.le.r_inner+2.*rx)closed=closed+1.
                work(i,j)=closed
            endif
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
    title=' : separatrix'
    call wtstr(.3,.975,label,2,0,0)
    call wtstr(.5,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    call conrec(work,mx,mx,my,1.4,1.6,.1,0,-1,-1012)
    !
    !     make constant latitude circles
    !
    dd=degrees/float(mx2-1)
    fmax=0.
    do j=1,my
        do i=1,mx
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            clat(i,j)=90.-dlat*dd
            fmax=amax1(fmax,clat(i,j))
        enddo
    enddo
    !
    !     ilab=1
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    !
    return
end
!
!
!	****************************************
!
!
subroutine auroras(press,rho,bsx,bsy,bsz, &
    nx,ny,nz,n_grids,box,radstrt,re_equiv,r_inner,iside, &
    time,save_dat,add_dip,label,ncon,write_dat, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make a close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !    
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91),chot(91,91)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    dimension press(nx,ny,nz),rho(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz)
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_dip,write_dat
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    !
    ilab=0
    add_dip=.false.
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.698132  ! 40 degrees
    !      theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
            chot(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    xmin=grd_xmin(box)+1.
    xmax=grd_xmax(box)-1.
    ymin=grd_ymin(box)+1.
    ymax=grd_ymax(box)-1.
    zmin=grd_zmin(box)+1.
    zmax=grd_zmax(box)-1.
    !
    do j=1,my
        do i=1,mx
            !
            !         find equivalent latitude
            !
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            dlat=amax1(dlat,0.0000001)
            alat=dlat*del_theta*theta_equiv
            if(alat.le.1.55) then
                cost=cos(alat)
                sint=sin(alat)
                !
                !         find equivalent longitude
                !
                cosp=(i-mx2)/dlat
                sinp=(j-my2)/dlat
                !
                !          find position on grid
                !
                z1=iside*(radstrt*cost)
                x1=radstrt*sint*cosp
                !
                !       dipole space to real space
                !
                ax=x1*cos_tilt+z1*sin_tilt+xdip
                az=-x1*sin_tilt+z1*cos_tilt+zdip
                ay=(radstrt*sint*sinp+ydip)
                !
                !          interpolate data to grid point
                !
                !
                ak=1.+(az-grd_zmin(box))/delz
                k1=ak
                k2=k1+1
                dz=ak-k1
                !
                aj=1.+(ay-grd_ymin(box))/dely
                j1=aj
                j2=j1+1
                dy=aj-j1
                !
                ai=1.+(ax-grd_xmin(box))/delx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                !         trace fld line and find max energy along it
                !
                ergies=0.
                rhod=0.
                dir=-delx
                call rungea(bsx,bsy,bsz,press,rho,nx,ny,nz, &
                box,ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax, &
                add_dip,ergies,rhod,maxpts,dir, &
                delx,dely,delz,r_inner,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                dir=+delx
                call rungea(bsx,bsy,bsz,press,rho,nx,ny,nz, &
                box,ax,ay,az,xmin,xmax,ymin,ymax,zmin,zmax, &
                add_dip,ergies,rhod,maxpts,dir, &
                delx,dely,delz,r_inner,n_grids, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    			!
                chot(i,j)=ergies
            endif
        enddo
    enddo
    !
    !
    !     initialize viewport and frame headings
    !
    !     plot intensity = current*energy
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,chot(i,j))
            fmin=amin1(fmin,chot(i,j))
        enddo
    enddo
    !
    call frame
    call gselnt(0)
    call gsplci(18)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : intensity'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    write(title,'(f7.4)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    title='fmax'//title
    write(title,'(f7.4)')fmax
    call wtstr(.25,.98,title,1,0,0)
    if(fmin*fmax.lt.0.)then
        cmax=amax1(fmax,abs(fmin))
        cmin=-cmax
    else
        cmax=fmax
        cmin=fmin
    endif
    clev=(cmax-cmin)/(ncon+1.)
    !     if(.not.start) call cnrccf(chot,mx,mx,my,cmin,
    !    +             cmax,clev,0,-1,-1012,2.5,0,1)
    irecmn=18
    call conrec(chot,mx,mx,my,cmin, &
    cmax,clev,0,-1,-1012)
    !
    !     plot constant latitude circles
    !
    dd=degrees/float(mx2-1)
    fmax=0.
    do j=1,my
        do i=1,mx
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            clat(i,j)=90.-dlat*dd
            fmax=amax1(fmax,clat(i,j))
        enddo
    enddo
    call conrec(clat,mx,mx,my,50.,90.,10.,0,-1,-1012)
    !
    !     save data if required
    !
    if(write_dat)then
        write(recdt_f,*) time, chot
    endif
    return
end
!
!
!	****************************************
!
!
subroutine aurora_cur(stuff,nx,ny,nz,box,radstrt,re_equiv, &
    r_inner,iside,time,save_dat,add_two,label,ncon,write_dat, &
    b_equiv,planet_rad,tot_cur,peak_cur,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91)
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension stuff(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_two,write_dat
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !     for no line labeling set ilab to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    call isoclrs
    ilab=0
    ioffm=1
    irecmj=1
    irectx=1
    irecmn=1
    !
    !       conversion of current simulation units to ua/box^2 and to ma
    !                 map to earth's surface so a ~ r**3
    !
    scale_cur= 10. * b_equiv *((radstrt*re_equiv)**3)/ &
    (4.* 3.1416*re_equiv*planet_rad)
    !       write(*,*)b_equiv,radstrt,re_equiv,planet_rad
    !       write(*,*)' current density in',scale_cur,' ua/box^2'
    scale_amps=((planet_rad/1000.)**2)*scale_cur
    !       write(*,*)'totcur scale',scale_amps
    amp_up=0.
    amp_down=0.
    peak_cur_up=0.
    peak_cur_down=0.
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    !      theta_range=0.698132  ! 40 degrees
    theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    rads=degrees/57.3
    del_rads=rads/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    do j=1,my
        do i=1,mx
            !
            !         find equivalent latitude
            !
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            dlat=amax1(dlat,0.0000001)
            alat=dlat*del_theta*theta_equiv  ! equiv theta at sample point
            !
            theta=dlat*del_rads  ! theta in ionosphere
            ascale=del_rads*sin(theta)/(dlat+1.)
            !
            if(alat.le.1.55) then
                cost=cos(alat)
                sint=sin(alat)
                !
                !         find equivalent longitude
                !
                cosp=(i-mx2)/dlat
                sinp=(j-my2)/dlat
                !
                !          find position on grid
                !
                z1=iside*(radstrt*cost)
                x1=radstrt*sint*cosp
                !
                !       dipole space to real space
                !
                ax=x1*cos_tilt+z1*sin_tilt+xdip
                az=-x1*sin_tilt+z1*cos_tilt+zdip
                ay=(radstrt*sint*sinp+ydip)
                !
                !          interpolate data to grid point
                !
                !
                ak=1.+(az-grd_zmin(box))/delz
                k1=ak
                k2=k1+1
                dz=ak-k1
                !
                aj=1.+(ay-grd_ymin(box))/dely
                j1=aj
                j2=j1+1
                dy=aj-j1
                !
                ai=1.+(ax-grd_xmin(box))/delx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                acur=stuff(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +stuff(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +stuff(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +stuff(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +stuff(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +stuff(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +stuff(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +stuff(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                acur=acur*scale_amps  !in physical unit ua/box**2
                work(i,j)=acur
                totcur=ascale*acur
                if(acur.gt.0)then
                    amp_up=amp_up+totcur
                    peak_cur_up=amax1(peak_cur_up,acur)
                else
                    amp_down=amp_down+totcur
                    peak_cur_down=amin1(peak_cur_down,acur)
                endif
            endif
            !
        enddo
    enddo
    !
    tot_cur=(amp_up-amp_down)/2.
    peak_cur=amax1(peak_cur_up,-peak_cur_down)
    write(fluxs_f,*) 'totcur: ', amp_up, amp_down, tot_cur
    write(fluxs_f,*) 'Peak current: ', peak_cur_up, peak_cur_down, peak_cur
    !
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : aurora'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
        enddo
    enddo
    !
    if (fmin.ge.0.0)then
        cmin=1.05*fmin
    else
        cmin=0.95*fmin
    endif
    !
    if(fmax.ge.0.0) then
        cmax=0.95*fmax
    else
        cmax=1.05*fmax
    endif
    finc=(cmax-cmin)/float(ncon-1)
    write(title,'(f7.2)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    write(title,'(f7.2)')fmax
    title='fmax'//title
    call wtstr(.3,.98,title,1,0,0)
    !
    write(title,'(f7.3)')amp_up
    title='tot_up'//title
    call wtstr(.1,.96,title,1,0,0)
    write(title,'(f7.3)')amp_down
    title='tot_down'//title
    call wtstr(.45,.96,title,1,0,0)
    !
    call gsplci(1)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
    !
    !     make constant latitude circles
    !
    dd=degrees/float(mx2-1)
    fmax=0.
    do j=1,my
        do i=1,mx
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            clat(i,j)=90.-dlat*dd
            fmax=amax1(fmax,clat(i,j))
        enddo
    enddo
    !
    !     ilab=1
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    !
    !     save data if required
    !
    if(save_dat)then
        if(iside.gt.0) then
            do j=1,my
                do i=1,mx
                    sdata(i,j)=work(i,j)
                enddo
            enddo
        else
            do j=1,my
                do i=1,mx
                    tdata(i,j)=work(i,j)
                enddo
            enddo
        endif
    endif
    !
    !      output data to log file if necessary
    !
    if(write_dat)then
        write(recdt_f,*) time, work
    endif
    !
    !
    if(.not.add_two)return
    !
    ilab=0
    if(iside.gt.0)then
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+sdata(i,j)
            enddo
        enddo
    else
        do j=1,my
            do i=1,mx
                work(i,j)=work(i,j)+tdata(i,j)
            enddo
        enddo
    endif
    !
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : total aur'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='t = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.85,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,work(i,j))
            fmin=amin1(fmin,work(i,j))
        enddo
    enddo
    !
    if (fmin.ge.0.0)then
        cmin=1.05*fmin
    else
        cmin=0.95*fmin
    endif
    !
    if(fmax.ge.0.0) then
        cmax=0.95*fmax
    else
        cmax=1.05*fmax
    endif
    finc=(cmax-cmin)/float(ncon-1)
    write(title,'(f6.3)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    write(title,'(f6.3)')fmax
    title='fmax'//title
    call wtstr(.25,.98,title,1,0,0)
    call conrec(work,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
    !
    !     ilab=1
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    return
end
!
!
!	****************************************
!
!
subroutine aurora_pot(efldx,efldy,efldz,bsx,bsy,bsz, &
    nx,ny,nz,n_grids,box,radstrt,re_equiv,iside,r_inner, &
    time,save_dat,add_dip,label,ncon,write_dat, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make a close up image
    !      of stuff near the auroral regions at a fixed distance radstrt
    !       re_equiv converts grid units to r_e
    !        add_two adds two current densities to produce
    !                 a total auroral map
    !
    common /space/sdata(91,91),tdata(91,91), &
    work(91,91),clat(91,91),chot(91,91)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    dimension efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical save_dat,add_dip,write_dat
    !
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    !
    ilab=0
    add_dip=.false.
    !
    !      dimension for plotted array
    !
    mx=91
    my=91
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !       re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    !      theta_range=0.698132  ! 40 degrees
    theta_range=0.6108652  ! 35 degrees
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
            chot(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    xmin=grd_xmin(box)+1.
    xmax=grd_xmax(box)-1.
    ymin=grd_ymin(box)+1.
    ymax=grd_ymax(box)-1.
    zmin=grd_zmin(box)+1.
    zmax=grd_zmax(box)-1.
    !
    do j=1,my
        do i=1,mx
            !
            !         find equivalent latitude
            !
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            dlat=amax1(dlat,0.0000001)
            alat=dlat*del_theta*theta_equiv
            if(alat.le.1.55) then
                cost=cos(alat)
                sint=sin(alat)
                !
                !         find equivalent longitude
                !
                cosp=(i-mx2)/dlat
                sinp=(j-my2)/dlat
                !
                !          find position on grid
                !
                z1=iside*(radstrt*cost)
                x1=radstrt*sint*cosp
                !
                !       dipole space to real space
                !
                ax=x1*cos_tilt+z1*sin_tilt+xdip
                az=-x1*sin_tilt+z1*cos_tilt+zdip
                ay=(radstrt*sint*sinp+ydip)
                !
                !         trace fld line and find max energy along it
                !
                dir=-delx
                call rungeb(efldx,efldy,efldz,bsx,bsy,bsz,nx,ny,nz, &
                n_grids,box,delx,ax,ay,az,xmin,xmax,ymin,ymax, &
                zmin,zmax,add_dip,pot1,maxpts,dir,r_inner, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                dir=+delx
                call rungeb(efldx,efldy,efldz,bsx,bsy,bsz,nx,ny,nz, &
                n_grids,box,delx,ax,ay,az,xmin,xmax,ymin,ymax, &
                zmin,zmax,add_dip,pot2,maxpts,dir,r_inner, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    			!
                chot(i,j)=iside*(pot1-pot2)
            endif
            !
        enddo
    enddo
    !
    !     initialize viewport and frame headings
    !
    !     plot intensity = current*energy
    !
    fmax=work(1,1)
    fmin=fmax
    do j=1,my
        do i=1,mx
            fmax=amax1(fmax,chot(i,j))
            fmin=amin1(fmin,chot(i,j))
        enddo
    enddo
    !
    call frame
    call gselnt(0)
    call gsplci(18)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    title=' : intensity'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    write(title,'(f7.4)')fmin
    title='fmin'//title
    call wtstr(.1,.98,title,1,0,0)
    title='fmax'//title
    write(title,'(f7.4)')fmax
    call wtstr(.25,.98,title,1,0,0)
    !     if(fmin*fmax.lt.0.)then
    !       cmax=amax1(fmax,abs(fmin))
    !       cmin=-cmax
    !     else
    !       cmax=0.5*fmax
    !       cmin=fmin
    !     endif
    cmax=0.15    !20 kev
    cmin=-0.15
    clev=(cmax-cmin)/(ncon+1.)
    !     if(.not.start) call cnrccf(chot,mx,mx,my,cmin,
    !    +             cmax,clev,0,-1,-1012,2.5,0,1)
    irecmn=18
    call conrec(chot,mx,mx,my,cmin, &
    cmax,clev,0,-1,-1012)
    !
    !     plot constant latitude circles
    !
    dd=degrees/float(mx2-1)
    fmax=0.
    do j=1,my
        do i=1,mx
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            clat(i,j)=90.-dlat*dd
            fmax=amax1(fmax,clat(i,j))
        enddo
    enddo
    call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    !
    !     save data if required
    !
    if(write_dat)then
        write(recdt_f,*) time, chot
    endif
    return
end
!
!
!	****************************************
!
!
subroutine cappot(chrg,pott,nx,ny,nz,n_grids,box,radstrt, &
    re_equiv,v_equiv,b_equiv,time,write_dat, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make potential drop
    !         over auroral oval, vx,vy,vz, assumed to
    !         be the inductive electric field
    !      answer is in kv: does both hemispheres at once to save work
    !      charge is the density over a sample region around
    !       the earth
    !
    common /space/charge(31,31,31),poten(31,31,31), &
    pot(61,61),clat(61,61)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    integer box
    dimension chrg(nx,ny,nz),pott(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical write_dat
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    !
    ilab=0
    !
    !      dimension for plotted array
    !
    mx=61
    my=61
    my2=my/2+1
    mx2=mx/2+1
    amx2=mx2-1.
    !
    msx=31
    msy=31
    msz=31
    ijk=(msx-1)/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latiude at radstrt
    !
    !      re_equiv=0.84
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.6108652
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      b_equiv=65.75
    !       b is in nt and v in km/s so total si conversion is 1e-6
    !       v_equiv=1e-3
    !       d is still in kilometers so the
    !        potential is in kilovolts
    !
    d_equiv=6371.*re_equiv
    scale=(1.e-6)*b_equiv*v_equiv*d_equiv
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    ncx=1.-grd_xmin(box)/delx
    ncy=1.-grd_ymin(box)/dely
    ncz=1.-grd_zmin(box)/delz
    !
    !     find the charge density equal to the div of vxb
    !
    do ks=1,msz
        k=(ks-ijk)+ncz
        do js=1,msy
            j=(js-ijk)+ncy
            do is=1,msx
                i=(is-ijk)+ncx
                charge(is,js,ks)=chrg(i,j,k)
                poten(is,js,ks)=pott(i,j,k)
            enddo
        enddo
    enddo
    !
    !      initialize work array
    !
    do mm=1,2
        iside=(-1)**(mm-1)
        do j=1,my
            do i=1,mx
                pot(i,j)=0.0
            enddo
        enddo
        !
        do j=1,my
            do i=1,mx
                !
                !         find equivalent latitude
                !
                dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
                dlat=amax1(dlat,0.0000001)
                alat=dlat*del_theta*theta_equiv
                if(alat.le.1.55) then
                    cost=cos(alat)
                    sint=sin(alat)
                    !
                    !         find equivalent longitude
                    !
                    cosp=(i-mx2)/dlat
                    sinp=(j-my2)/dlat
                    !
                    !          find position on grid
                    !
                    z1=iside*(radstrt*cost)
                    x1=radstrt*sint*cosp
                    !
                    !       dipole space to real space
                    !
                    ax=x1*cos_tilt+z1*sin_tilt+xdip
                    az=-x1*sin_tilt+z1*cos_tilt+zdip
                    ay=(radstrt*sint*sinp+ydip)
                    ar=sqrt((ax-xdip)**2+(ay-ydip)**2+(az-zdip)**2) &
                    +0.0000001
                    !
                    !          interpolate data to grid point
                    !
                    ak=1.+(az-grd_zmin(box))/delz
                    k1=ak
                    k2=k1+1
                    dz=ak-k1
                    !
                    aj=1.+(ay-grd_ymin(box))/dely
                    j1=aj
                    j2=j1+1
                    dy=aj-j1
                    !
                    ai=1.+(ax-grd_xmin(box))/delx
                    i1=ai
                    i2=i1+1
                    dx=ai-i1
                    !
                    apot=chrg(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                    +chrg(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                    +chrg(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                    +chrg(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                    +chrg(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                    +chrg(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                    +chrg(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                    +chrg(i2,j2,k2)*(dx)*(dy)*(dz)
                    !
                    pot(i,j)=-apot
                endif
            enddo
        enddo
        !
        !
        call frame
        call gselnt(0)
        call gsplci(1)
        call agseti('frame.',2)
        call agseti('set.',-1)
        !
        if(iside.gt.0) title= 'nth cap pot'
        if(iside.le.0) title= 'sth cap pot'
        call wtstr(.5,.975,title,2,0,0)
        write(title,'(f7.3)')time
        title='ut = '//title
        call wtstr(.8,.975,title,2,0,0)
        !
        ncon=21
        fmax=pot(1,1)
        fmin=fmax
        do j=1,my
            do i=1,mx
                fmax=amax1(fmax,pot(i,j))
                fmin=amin1(fmin,pot(i,j))
            enddo
        enddo
        !
        if (fmin.ge.0.0)then
            cmin=1.05*fmin
        else
            cmin=0.95*fmin
        endif
        !
        if(fmax.ge.0.0) then
            cmax=0.95*fmax
        else
            cmax=1.05*fmax
        endif
        finc=(cmax-cmin)/float(ncon-1)
        write(title,'(f7.2)')fmin
        title='fmin'//title
        call wtstr(.1,.98,title,1,0,0)
        write(title,'(f7.2)')fmax
        title='fmax'//title
        call wtstr(.25,.98,title,1,0,0)
        call gsplci(1)
        call conrec(pot,mx,mx,my,cmin,cmax,finc,0,-1,-1012)
        !
        !     make constant latitude circles
        !
        dd=degrees/float(mx2-1)
        fmax=0.
        do j=1,my
            do i=1,mx
                dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
                clat(i,j)=90.-dlat*dd
                fmax=amax1(fmax,clat(i,j))
            enddo
        enddo
        !
        !     ilab=1
        call conrec(clat,mx,mx,my,60.,90.,10.,0,-1,-1012)
    enddo
    !
    !
    !      output data to data file if necessary
    !
    if(write_dat)then
		write(recdt_f,*) time,scale,cos_tilt,sin_tilt,radstrt,b_equiv,v_equiv,re_equiv,rx,charge
    endif
    return
end
!
!
!	****************************************
!
!
subroutine cnrccf(zdat,kzdt,mzdt,nzdt,flow,fhgh,finc,nset,nhgh, &
    ndsh,smooth,lines,icolfg)
    !
    !  This is a version of the old conrec which does color-filled contours.  this
    !  is a modified cpcnrc which was written by dave kennison.  it requires
    !  the conpack library.   the first ten arguments are identical to the
    !  old conrec call.
    !
    !  smooth - zero for no smoothing; small values (e.g. .001) yield
    !     approximately cubic splines; large values (e.g. 50.) yield nearly
    !     polygonal curves; suggested starting value is 2.5; if negative,
    !     smoothing is done before coordinate mapping
    !
    !  lines - if nonzero, contour lines are drawn in the foreground color
    !
    !  icolfg - if zero, a default color table is used; if nonzero, user
    !     must set color indices 2 to ncl+2, where ncl = number of contour
    !     levels.  the number of colors required is always one more than the
    !     number of contour levels.  color index numbering for associated
    !     contour fill colors must begin with 2 since color indices 0 and 1
    !     are reserved for background and foreground, respectively.  these
    !     may be set by the user regardless of the value of icolfg.
    !     choosing the default color table will reset indices 2 to maxcol.
    !     a subset of these will actually be used.  see routine dfclrs for
    !     the value of maxcol.
    !
    integer box
    dimension zdat(kzdt,*)
    !
    ! define some needed dimensions.
    !
    parameter (lama=100000,lrwk=2000,liwk=1000,ncra=2000,ngrps=10, &
    locv=10)
    !
    ! define required workspace arrays.
    !
    dimension rwrk(lrwk),iwrk(liwk),iama(lama),xcra(ncra),ycra(ncra), &
    iaia(ngrps),igia(ngrps),rec(4)
    !
    ! logical to fix spaghetti code - mat
    !
    logical lexit
    !
    ! define a character variable to use for point-value labelling.
    !
    character*(locv) croz
    character*100 ilts
    common/mapawm/klat,icmp,icon
    !
    ! declare the contour-line drawing routine.
    !
    external colram
        lexit = .false.
    !
    !  ..set gks internal parameters
    !
    call gqclip(ier,iclip,rec)
    call gqfais(ier,ifais)
    call gsclip(0)
    call gsfais(1)
    !
    !  ..set the tension on the two dimensional smoother
    !
    call cpsetr('t2d',smooth)
    !
    ! arrange for the selection of contour levels as desired by the user.
    !
    if (finc.lt.0.) then
        call cpseti('cls - contour level selector',max(1,int(-finc)))
        call cpsetr('cis - contour interval specifier',0.)
        elseif (finc.eq.0.) then
        call cpseti('cls - contour level selector',16)
        call cpsetr('cis - contour interval specifier',0.)
    else
        call cpseti('cls - contour level selector',1)
        call cpsetr('cis - contour interval specifier',finc)
        if (flow.lt.fhgh) then
            call cpsetr('cmn - contour minimum',flow)
            call cpsetr('cmx - contour maximum',fhgh)
        endif
    endif
    !
    ! set up the desired mapping of output onto the plotter frame.
    !
    if (nset.lt.0) then
        call cpseti('set - do-set-call flag',1)
        call getset(xvpl,xvpr,yvpb,yvpt,xwdl,xwdr,ywdb,ywdt,lnlg)
        call cpsetr('vpl - viewport left edge',xvpl)
        call cpsetr('vpr - viewport right edge',xvpr)
        call cpsetr('vpb - viewport bottom edge',yvpb)
        call cpsetr('vpt - viewport top edge',yvpt)
        call cpseti('vps - viewport shape',0)
        elseif (nset.eq.0) then
        call cpseti('set - do-set-call flag',1)
        call cpsetr('vpl - viewport left edge',.05)
        call cpsetr('vpr - viewport right edge',.95)
        call cpsetr('vpb - viewport bottom edge',.05)
        call cpsetr('vpt - viewport top edge',.95)
        call cpseti('vps - viewport shape',4)
    else
        call cpseti('set - do-set-call flag',0)
    endif
    !
    ! decide what dash pattern to use.
    !
    if (lines.ne.0) then
        idsh=abs(ndsh)
        if (idsh.eq.0.or.idsh.eq.1.or.idsh.eq.1023) then
            idsh=ior(ishift(32767,1),1)
        else
            idsh=ior(ishift(idsh,6),iand(ishift(idsh,-4),63))
        endif
    endif
    !
    ! decide whether to label highs and lows or not.
    !
    if (nhgh.eq.0) then
        call cpsetc('hlt - high/low label text', &
        'h:b:$zdv$:e:''l:b:$zdv$:e:')
    else
        call cpsetc('hlt - high/low label text',' ')
    endif
    !
    ! initialize conpack and give it all array dimensions.
    !
    call cprect(zdat,kzdt,mzdt,nzdt,rwrk,lrwk,iwrk,liwk)
    !
    ! pick contour levels.
    !
    call cppkcl(zdat,rwrk,iwrk)
    !
    ! retrieve the contour levels selected, one at a time.  discard levels
    ! which are outside the range, if any, specified by the user-supplied
    ! values of flow and fhgh, and move the parameters for all remaining
    ! levels to the beginning of the parameter arrays.  set dash patterns
    ! for all levels.  the value of 'ciu' must be saved for possible
    ! restoration if it gets clobbered as a side effect of setting contour
    ! level 1.
    !
    call cpgetr('ciu - contour interval used',cinu)
    call cpgeti('ncl - number of contour levels',nclo)
    ncln=0
    do iclo=1,nclo
        call cpseti('pai - parameter array index',iclo)
        call cpgetr('clv - contour level',clev)
        if (flow.ge.fhgh.or.(clev.ge.flow.and.clev.le.fhgh)) then
            ncln=ncln+1
            if (ncln.ne.iclo) then
                call cpgeti('clu - contour level use flag',iclu)
                call cpseti('pai - parameter array index',ncln)
                call cpsetr('clv - contour level',clev)
                call cpseti('clu - contour level use flag',iclu)
                call cpseti('aia - area identifier above level',ncln+1)
                call cpseti('aib - area identifier below level',ncln)
                call cpseti('clc - contour line color index',-1)
                call cpsetc('cld - contour line dash pattern', &
                '$$$$$$$$$$$$$$$$')
                call cpseti('cll - contour line line width',-1)
                call cpseti('llc - line label color index',-1)
                call cpsetc('llt - line label text',' ')
            endif
        endif
        if (ndsh.gt.0.or.(ndsh.lt.0..and.clev.lt.0.)) &
        call cpseti('cld - contour line dash pattern',idsh)
    enddo
    !
    ! log an error if no contour levels were within the user's bounds.
    !
    if (ncln.eq.0) then
        call seter('cnrccf - no contour levels in specified range', &
        1,2)
        return
    endif
    !
    ! if the number of contour levels decreased, reset parameters affected.
    !
    if (ncln.lt.nclo) then
        call cpseti('ncl - number of contour levels',ncln)
        call cpsetr('ciu - contour interval used',cinu)
    endif
    !
    !  ..default color table
    !
    if (icolfg.eq.0) then
        ncol=ncln+1
        call dfclrs(ncol)
    endif
    !
    !  ..color filled contour
    !
    call arinam(iama,lama)
    call cpclam(zdat,rwrk,iwrk,iama)
    call arscam(iama,xcra,ycra,ncra,iaia,igia,ngrps,colram)
    !
    !  ..map
    !
    !     call getusv('ii',icolor)
    !     call setusv('ii',icmp)
    !     call mapgrd
    !     call maplbl
    !     call setusv('lw',2000)
    !     call maplot
    !     call setusv('lw',1000)
    !     call maplbm
    !     call setusv('ii',icolor)
    !
    ! if requested, put out a simple background.
    !
    if (nset.eq.0) call cpback(zdat,rwrk,iwrk)
    !
    ! see how the user has chosen to position contour levels.
    !
    call cpgeti('llp - line label positioning flag',llpf)
    !
    ! draw the contour lines, masking them if necessary.
    !
    if (lines.ne.0) then
        if (llpf.le.1) then
            call gqplci(ier,icolor)
            call gsplci(0)
            call cpcldr(zdat,rwrk,iwrk)
            call gsplci(icolor)
        else
            call arinam(iama,lama)
            call cplbam(zdat,rwrk,iwrk,iama)
            call cpcldm(zdat,rwrk,iwrk,iama,cpdrpl)
        endif
    endif
    !
    ! plot labels.
    !
    call cpgetc('ilt - informational label text',ilts)
    !     print*,'informational label is ',ilts
    call cplbdr(zdat,rwrk,iwrk)
    !
    ! if requested, label every point on the grid.
    !
    if (nhgh.gt.0) then
        call getset(xvpl,xvpr,yvpb,yvpt,xwdl,xwdr,ywdb,ywdt,lnlg)
        call cpgetr('cwm - character width multiplier',chwm)
        call cpgetr('hla - high/low label angle',angd)
        call cpgetr('hls - high/low label size',size)
        call cpgeti('map - mapping flag',imap)
        call cpgetr('orv - out-of-range value',orva)
        call cpgetr('spv - special value',spva)
        call cpgetr('xc1 - x coordinate at i = 1',xca1)
        call cpgetr('xcm - x coordinate at i = m',xcam)
        call cpgetr('yc1 - y coordinate at j = 1',yca1)
        call cpgetr('ycn - y coordinate at j = n',ycan)
        size=(xvpr-xvpl)*chwm*size
        if (xca1.eq.xcam) then
            xca1=1.
            xcam=real(mzdt)
        endif
        if (yca1.eq.ycan) then
            yca1=1.
            ycan=real(nzdt)
        endif
        do j=1,nzdt
            ypos=yca1+real(j-1)*(ycan-yca1)/real(nzdt-1)
            do i=1,mzdt
                xpos=xca1+real(i-1)*(xcam-xca1)/real(mzdt-1)
                if (spva.eq.0..or.zdat(i,j).ne.spva) then
                    call cpsetr('zdv - z data value',zdat(i,j))
                    call cpgetc('zdv - z data value',croz)
                    do k=locv,2,-1
                        if (croz(k:k).ne.' ') then
                            lcrz=k
                            lexit = .true.
                            exit
                        endif
                    enddo
                    if (.not.lexit) lcrz=1
                    if (imap.eq.0) then
                        call plchhq(xpos,ypos,croz(1:lcrz),size,angd,0.)
                    else
                        call cpmpxy(imap,xpos,ypos,xmpd,ympd)
                        if (orva.eq.0..or.xmpd.ne.orva) &
                        call plchhq(xmpd,ympd,croz(1:lcrz),size,angd,0.)
                    endif
                endif
            enddo
        enddo
    endif
    !
    !  ..reset gks internal parameters
    !
    call gsclip(iclip)
    call gsfais(ifais)
    return
end
!
!
!	****************************************
!
!
subroutine cpclrs
    !
    !      8 basic colors, 4 extras for curve tracing, 15 more for isoplots
    !
    real rgbv(3,12)
	!
	rgbv(:,1) = (/ 0.,0.,0. /)     !black 
   	rgbv(:,2) = (/ 1.,0.,0. /)    !red 
	rgbv(:,3) = (/ 0.,1.,0. /)   !green 
	rgbv(:,4) = (/ 0.,0.,1. /)   !blue 
	rgbv(:,5) = (/ 0.,1.,1. /)   !cyan 
	rgbv(:,6) = (/ 1.,0.,1. /)   !magenta 
	rgbv(:,7) = (/ 1.,1.,0. /)   !yellow 
	rgbv(:,8) = (/ 0.5, 1.0, 0.5 /)    !bright green 
	rgbv(:,9) = (/ 1.0, 0.5, 1.0 /)    !bright pink 
	rgbv(:,10) = (/ 0.3, 1.0, 1.0 /)    !bright blue 
	rgbv(:,11) = (/ 1.0, 1.0, 0.3 /)    !bright yellow 
	rgbv(:,12) = (/ 1.,1.,1. /)     !white
    !
    call gscr(1,0,1.,1.,1.)
    !
    do i = 1, 12
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine isoclrs
    real rgbv(3,16)
	!
	rgbv(:,1) = (/ 0.,0.,0. /)
	rgbv(:,2) = (/ .7,.7,.7 /)
	rgbv(:,3) = (/ .75,.5,1. /)
	rgbv(:,4) = (/ .5,0.,1. /)
	rgbv(:,5) = (/ 0.,0.,1. /)
	rgbv(:,6) = (/ 0.,.5,1. /)
	rgbv(:,7) = (/ 0.,1.,1. /)
	rgbv(:,8) = (/ 0.,1.,.6 /)
	rgbv(:,9) = (/ 0.,.85,0. /)
	rgbv(:,10) = (/ .7,1.,0. /)
	rgbv(:,11) = (/ 1.,1.,0. /)
	rgbv(:,12) = (/ 1.,.75,0. /)
	rgbv(:,13) = (/ 1.,.38,.38 /)
	rgbv(:,14) = (/ 1.,0.,.75 /)
	rgbv(:,15) = (/ 1.,0.,0. /)
	rgbv(:,16) = (/ 1.,1.,1. /)
    call gscr(1,0,1.,1.,1.)
    do i=1,16
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    return
end
!
!
!	****************************************
!
!
subroutine isoclrs_hot
    !
    !     arrays to plot colorbar
    !
    integer lind(20)
    real tcon(20)
    character*4 llbs(20)
    !
    real rgbv(3,21)
	!
	rgbv(:,1) = (/ 0.,0.,0. /)
	rgbv(:,2) = (/ .7,.7,.7 /)
	rgbv(:,3) = (/ .75,.5,1. /)
	rgbv(:,4) = (/ .5,0.,1. /)
	rgbv(:,5) = (/ .25,0.,1. /)
	rgbv(:,6) = (/ 0.,0.,1. /)
	rgbv(:,7) = (/ 0.,.25,1. /)
	rgbv(:,8) = (/ 0.,.5,1. /)
	rgbv(:,9) = (/ 0.,1.,1. /)
	rgbv(:,10) = (/ 0.,1.,.6 /)
	rgbv(:,11) = (/ 0.,.85,0. /)
	rgbv(:,12) = (/ 0.3,.85,0. /)
	rgbv(:,13) = (/ .7,1.,0. /)
	rgbv(:,14) = (/ 1.,1.,0. /)
	rgbv(:,15) = (/ 1.,.75,0. /)
	rgbv(:,16) = (/ 1.,.38,.38 /)
	rgbv(:,17) = (/ 1.,0.2,.2 /)
	rgbv(:,18) = (/ 1.,0.1,0.1 /)
	rgbv(:,19) = (/ 1.,0.,0. /)
	rgbv(:,20) = (/ 0.85,0.,0. /)
	rgbv(:,21) = (/ 1.,1.,1. /)
    call gscr(1,0,1.,1.,1.)
    do i=1,21
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    return
end
!
!
!	****************************************
!
!
subroutine isorb
	!
	!	creates a color map from the following RGB vectors:
	!
    real rgbv(3,16)
	!
	rgbv(:,1) = (/ 0.,0.,0. /)
	rgbv(:,2) = (/ .7,.7,.7 /)
	rgbv(:,3) = (/ .75,.5,1. /)
	rgbv(:,4) = (/ .5,0.,1. /)
	rgbv(:,5) = (/ 0.,0.,1. /)
	rgbv(:,6) = (/ 0.,.5,1. /)
	rgbv(:,7) = (/ 0.,1.,1. /)
	rgbv(:,8) = (/ 0.,1.,.6 /)
	rgbv(:,9) = (/ 1.,.38,.38 /)
	rgbv(:,10) = (/ 1.,0.,.38 /)
	rgbv(:,11) = (/ 1.,0.,0. /)
	rgbv(:,12) = (/ 0.,0.65,0. /)
	rgbv(:,13) = (/ 0.,.85,0. /)
	rgbv(:,14) = (/ .7,1.,0. /)
	rgbv(:,15) = (/ 1.,1.,0. /)
	rgbv(:,16) = (/ 1.,.75,0. /)
    call gscr(1,0,1.,1.,1.)
    do  i=1,16
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    return
end
!
!
!	****************************************
!
!
subroutine dfclrs(ncol)
    !
    !  reads from the default input ($in on the cray) a color table of
    !  maxcol colors in free format 8-bit integer rgb values (i.e. 0 to 255)
    !  and selects ncol colors evenly distributed among the maxcol colors;
    !  these are bound to color indices 2 to ncol+1
    !
    integer, parameter :: maxcol=250
    dimension irgbv(3,maxcol)
    data icall/0/
    save icall,irgbv
    if (ncol.le.0.or.ncol.gt.maxcol) then
        print*,'stop in dfclrs - input parameter ncol out of range'
        print*,'ncol = ',ncol
        stop
    endif
    if (icall.ne.1) then
        icall=1
        read(9,*) irgbv
    endif
    ainc=real(maxcol-1)/real(ncol-1)
    do jj=2,ncol+1
        j=nint(ainc*(jj-2))+1
        r=irgbv(1,j)/255.
        g=irgbv(2,j)/255.
        b=irgbv(3,j)/255.
        call gscr(1,jj,r,g,b)
    enddo
     return
end
!
!
!	****************************************
!
!
subroutine colram(xcra,ycra,ncra,iaia,igia,naia)
    dimension xcra(*),ycra(*),iaia(*),igia(*)
    if (iaia(1).le.0) return
    !
    !  ..area numbering starts with 1, but color index numbering starts
    !    with 2 because 0 and 1 are reserved for background and foreground,
    !    respectively
    !
    call gsfaci(iaia(1))
    call gfa(ncra-1,xcra,ycra)
    return
end
!
!
!	****************************************
!
!
subroutine conalong(stuff,cross,nx,ny,nz,box, &
    xcraft,ncraft,re_equiv,time,label,start,frac,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot x-z quantities at the
    !       spacecraft position
    !       percent - percentage level of max value for isosurface
    !        xcraft is the position of the spacecraft: assumed to
    !                  be in simulation units
    !        ncraft is the number of spacecraft to plot
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
    dimension xcraft(4,ncraft),stuff(nx,ny,nz),cross(nx,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    character*4 wd1,wd2,wd3
    character*12 label
    character*20 title
    !
    !      dimension for plotted array
    !
    my=nx
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
            j1=aj
            j2=aj+1
            dy=aj-j1
            ddy=1.-dy
            fmin=0.
            fmax=0.
            do k=1,nz
                do i=1,nx
                    cross(i,k)=stuff(i,j1,k)*ddy+ &
                    stuff(i,j2,k)*dy
                    fmin=amin1(fmin,cross(i,k))
                    fmax=amax1(fmax,cross(i,k))
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
                cmax=fmax
                cmin=fmin
            endif
            cmin=cmin/frac
            cmax=cmax/frac
            clev=(cmax-cmin)/(ncon+1.)
            if(.not.start) call cnrccf(cross,my,my,mz,cmin, &
            cmax,clev,0,-1,-1012,2.5,0,1)
            !     irecmn=18
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
                xrays(1)=ai-.5
                xrays(2)=ai+.5
                yrays(1)=ak+shift
                yrays(2)=ak+shift
                call curve(xrays,yrays,2)
                xrays(1)=ai+shift
                xrays(2)=ai+shift
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
!	****************************************
!
!
subroutine conalong_fix(stuff,cross,nx,ny,nz,box, &
    xcraft,ncraft,re_equiv,time,label,start,alo,ahi,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot x-z quantities at the
    !       spacecraft position
    !       percent - percentage level of max value for isosurface
    !        xcraft is the position of the spacecraft : assumed to
    !                  be in simulation units
    !        ncraft is the number of spacecraft to plot
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
    dimension xcraft(4,ncraft),stuff(nx,ny,nz),cross(nx,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    character*4 wd1,wd2,wd3
    character*12 label
    character*20 title
    !
    !      dimension for plotted array
    !
    my=nx
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
        j1=aj
        j2=aj+1
        dy=aj-j1
        ddy=1.-dy
        fmin=0.
        fmax=0.
        do k=1,nz
            do i=1,nx
                cross(i,k)=stuff(i,j1,k)*ddy+ &
                stuff(i,j2,k)*dy
                fmin=amin1(fmin,cross(i,k))
                fmax=amax1(fmax,cross(i,k))
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
        !
        !     draw in position of spacecraft
        !
        call gsplci(16)
        nl=11
        n2=(nl-1)/2
        do nn=1,nl
            ns=(nn-n2)/2
            shift=ns*.025
            xrays(1)=ai-.5
            xrays(2)=ai+.5
            yrays(1)=ak+shift
            yrays(2)=ak+shift
            call curve(xrays,yrays,2)
            xrays(1)=ai+shift
            xrays(2)=ai+shift
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
!	****************************************
!
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
!	****************************************
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
!
!
!	****************************************
!
!
subroutine concross(stuff,vx,vy,vz,bx,by,bz,nx,ny,nz, &
    n_grids,mm,box,xcraft,ncraft,re_equiv,r_inner, &
    xmin,xmax,ymin,ymax,zmin,zmax,time, &
    label,nlevs,ncon,strtch,add_dip,ivel,kht,start, &
    tx,ty,tz,tt,t,t2,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
    !        nlevs number of isossurfaces to be plotted - max 4
    !        ncon  number of contours to be plotted     - max 14
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    dimension t(mx,my,mz),tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz), &
    tt(mx,my,mz),t2(mx,my,mz2),work(muvwp2,muvwp2)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    real xray(1000),yray(1000),zray(1000), &
    xray1(1000),yray1(1000),zray1(1000)
    dimension stuff(nx,ny,nz,n_grids),vx(nx,ny,nz),xcraft(4,ncraft), &
    vy(nx,ny,nz),vz(nx,ny,nz), &
    bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    real xrays(2),yrays(2)
    real eye(3),tlev(6),tcon(14)
    character*4 llbs(14)
    character*5 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    integer box
    logical add_dip,start
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
    !
    ilab=0
    ioffm=10
    !
    !     skip parameters for field lines and arrows
    !
    jskip=(ny-1)/20
    iskip=(nx-1)/30+1
    !
    jskip=4
    iskip=4
    !
    !      dimension for 3-d plotted array
    !
    my2=my/2+1
    !
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
    axmax=amin1(xmax,grd_xmax(box)-.0001)
    aymax=amin1(ymax,grd_ymax(box)-.0001)
    azmax=amin1(zmax,grd_zmax(box)-.0001)
    rx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
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
    call filldatav(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
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
    tma=0.8*tma
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
    !     initialize eye position
    !
    eye(1)=mx/2.
    eye(2)=-my*2.5
    eye(3)=mz2*4.0
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    !     title=' : dawn dusk'
    call wtstr(.4,.975,label,2,0,0)
    !     call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
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
    ijump=-1
    ncol=3
    kk=2
    k=mz2+kht
    !     jskip=4
    !     iskip=8
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
                    !
                    !       draw magnetic field line
                    !
                endif
            endif
            l=0
            dir=0.25*delx
            xi=xmin+delx*(i-1)
            yi=ymin+dely*(j-1)
            zi=zmin+delz*(k-1)
            call trace(bx,by,bz,nx,ny,nz,box,add_dip, &
            xi,yi,zi,dir,maxpts,xf,yf,zf, &
            xray,yray,zray,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l,r_inner,n_grids, &
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
            call trace(bx,by,bz,nx,ny,nz,box,add_dip, &
            xi,yi,zi,dir,maxpts,xf1,yf1,zf1, &
            xray1,yray1,zray1,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l1,r_inner,n_grids, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            !
            rad3=sqrt(xray1(l1)**2+yray1(l1) **2+zray1(l1)**2)-2.
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
            else  if((rad2.le.r_inner).and.(rad3.le.r_inner))then
                call gsplci(5)
                call curve3(xray,yray,zray,l)
                call curve3(xray1,yray1,zray1,l1)
            else if((rad1.le.r_inner).or.(rad2.le.r_inner) &
                .or.(rad3.le.r_inner))then
                ijump=-ijump
                if(ijump.ge.0) then 
                    call gsplci(7)
                    call curve3(xray,yray,zray,l)
                    call curve3(xray1,yray1,zray1,l1)
                endif
            else
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
    !       wspace set writing distance to mlt and lat
    wspace=0.015
    nspace=0
    size=0.25*delx
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
            x1=1.+(ax-size-xmin)/delx
            x2=1.+(ax+size-xmin)/delx
            y1=1.+(ay-size-ymin)/dely
            y2=1.+(ay+size-ymin)/dely
            z1=1.+(az+height-size-zmin)/delz
            z2=1.+(az+height+size-zmin)/delz
            r1=1.+(ar-size-ymin)/dely
            r2=1.+(ar+size-ymin)/dely
            !
            !        2-d projection
            !
            call sflush
            call gsplci(9)
            call line3(x1,y1,1.,x2,y1,1.)
            call line3(x1,y1,1.,x1,y2,1.)
            call line3(x2,y2,1.,x1,y2,1.)
            call line3(x2,y2,1.,x2,y1,1.)
            !
            !        3-d position
            !
            call sflush
            call gsplci(9)
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
        call rungem(bx,by,bz,nx,ny,nz,box,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,r_inner, &
        add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        test for connection to earth and write out position
        !
        call gsplci(1)
        sx=(xray(npts)-xdip)*re_equiv
        sy=(yray(npts)-ydip)*re_equiv
        sz=(zray(npts)-zdip)*re_equiv
        rtot=sqrt(sx**2+sy**2+sz**2)
        if((rtot/re_equiv.le.r_inner+2.).and. &
            (rtot/re_equiv.gt.0.5*r_inner))then
            nspace=nspace+1
            !
            !        find polar coordinate : first rotate by dipole tilt
            !
            sx1=sx*cos_tilt-sz*sin_tilt
            sz1=sx*sin_tilt+sz*cos_tilt
            sy1=sy
            rtot=sqrt(sx1**2+sy1**2+sz1**2)
            if(sz1.ge.0)then
                colat=acos(sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=90.-colat
            else
                colat=acos(-sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=-(90.-colat)
            endif
            !
            !          find local time
            !
            if(sy1.ge.0.)then
                alt=atan2(sy1,sx1)
                alt=12.*alt/3.1416
            else
                alt=atan2(-sy1,sx1)
                alt=12.*(2.-alt/3.1416)
            endif
            arad=1.
            !        write(wd1,'(f5.1)')arad
            !        write(wd2,'(f5.1)')alat
            !        write(wd3,'(f5.1)')alt
            !        title=wd1//','//wd2//','//wd3//' '
            !        call wtstr(.4,.9-wspace*(nspace-1),title,1,0,0)
            !        write(title,'(i2)')n
            !        title='craft lat,mlt'//title
            !        call wtstr(.20,.9-wspace*(nspace-1),title,1,0,0)
        endif
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
            yray1(nn)=1.+(yray(nn)-ymin)/dely
            !        yray1(nn)=1.+(ar-ymin)/dely
            zray1(nn)=1.
        enddo
        !        call sflush
        !        call gsplci(1)
        call curve3(xray1,yray1,zray1,npts)
        !
        !       go the other direction
        dir=+0.25*delx
        call rungem(bx,by,bz,nx,ny,nz,box,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,r_inner, &
        add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        sx=(xray(npts)-xdip)*re_equiv
        sy=(yray(npts)-ydip)*re_equiv
        sz=(zray(npts)-zdip)*re_equiv
        rtot=sqrt(sx**2+sy**2+sz**2)
        if((rtot/re_equiv.le.r_inner+2.).and. &
            (rtot/re_equiv.gt.0.5*r_inner))then
            nspace=nspace+1
            !
            !        find polar coordinate : first rotate by dipole tilt
            !
            sx1=sx*cos_tilt-sz*sin_tilt
            sz1=sx*sin_tilt+sz*cos_tilt
            !
            sy1=sy
            rtot=sqrt(sx1**2+sy1**2+sz1**2)
            if(sz1.ge.0)then
                colat=acos(sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=90.-colat
            else
                colat=acos(-sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=-(90.-colat)
            endif
            !
            !          find local time
            !
            if(sy1.gt.0.)then
                alt=atan2(sy1,sx1)
                alt=12.*alt/3.1416
            else
                alt=atan2(-sy1,sx1)
                alt=12.*(2.-alt/3.1416)
            endif
            arad=1.
            !        write(wd1,'(f5.1)')arad
            !        write(wd2,'(f5.1)')alat
            !        write(wd3,'(f5.1)')alt
            !        title=wd1//','//wd2//','//wd3//' '
            !        call wtstr(.4,.9-wspace*(nspace-1),title,1,0,0)
            !        write(title,'(i2)')n
            !        title='craft lat,mlt'//title
            !        call wtstr(.20,.9-wspace*(nspace-1),title,1,0,0)
        endif
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
            yray1(nn)=1.+(yray(nn)-ymin)/dely
            !        yray1(nn)=1.+(ar-ymin)/dely
            zray1(nn)=1.
        enddo
        !        call sflush
        !        call gsplci(1)
        call curve3(xray1,yray1,zray1,npts)
        !
    enddo
    return
end
!
!
!	****************************************
!
!
subroutine conflow(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
    xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
    time,label,nlevs,ncon,ivel,strtch, &
    tx,ty,tz,t,t2,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
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
     real eye(3),tlev(6),tcon(14)
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    integer box
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
    !
    !     skip parameters for field lines and arrows
    !
    !      jskip=(ny-1)/40
    !      iskip=(nx-1)/40+1
    jskip=3
    iskip=3
    !
    !      dimension for plotted array
    !
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
    vm=0.0
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
    axmax=amin1(xmax,grd_xmax(box)-.00001)
    aymax=amin1(ymax,grd_ymax(box)-.00001)
    azmax=amin1(zmax,grd_zmax(box)-.00001)
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
    delx,dely,delz,vm,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
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
    ampl=(abs(tma)+abs(tmi))/2.+1.e-6
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
    !     ***************************************
    !       make contour plots down the middle
    !
    my2=my/2+1
    aymin=(ymin+aymax)/2.
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
                    if((ivel.eq.1).and.(dx.ge.0.0)) call gsplci(15)
                    if((ivel.eq.1).and.(dx.lt.0.0)) call gsplci(5)
                    if((ivel.eq.2).and.(dy.ge.0.0)) call gsplci(15)
                    if((ivel.eq.2).and.(dy.lt.0.0)) call gsplci(5)
                    if((ivel.eq.3).and.(dz.ge.0.0)) call gsplci(15)
                    if((ivel.eq.3).and.(dz.lt.0.0)) call gsplci(5)
                    !
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
    !
    !     draw points: equatorial
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
                    if(dx.ge.0.0) call gsplci(15)
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
!
!
!	****************************************
!
!
subroutine conhot(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
    xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
    time,label,nlevs,ncon,ivel,strtch,tlim, &
    tx,ty,tz,t,t2,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
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
    !      jskip=(ny-1)/40
    !      iskip=(nx-1)/40
    jskip=3
    iskip=3
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
    !      effective step size between points - slightly less than
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
    dcon=(tlim)/(ncon+1.)
    ampl=(abs(tlim))/2.
    write(magnif,'(1pe9.2)') ampl
    do n=1,ncon
        tcon(n)=dcon*n
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
    !     loading y-z plane
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
    !     draw points equatorial region
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
!
!
!	****************************************
!
!
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
!
!
!	****************************************
!
!
subroutine conmap(stuff,vx,vy,vz,bx,by,bz,nx,ny,nz,n_grids,mm,box, &
    xcraft,ncraft,re_equiv,r_inner, &
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
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    real xray(1000),yray(1000),zray(1000), &
    xray1(1000),yray1(1000),zray1(1000)
    dimension stuff(nx,ny,nz,n_grids),vx(nx,ny,nz),xcraft(4,ncraft), &
    vy(nx,ny,nz),vz(nx,ny,nz), &
    bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    real eye(3),tlev(6),tcon(14)
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    integer box
    logical add_dip
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
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
    axmax=amin1(xmax,grd_xmax(box)-.0001)
    aymax=amin1(ymax,grd_ymax(box)-.0001)
    azmax=amin1(zmax,grd_zmax(box)-.0001)
    !
    rx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
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
    call filldatav(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
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
    tisom=r_inner
    call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
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
            call trace(bx,by,bz,nx,ny,nz,box,add_dip, &
            xi,yi,zi,dir,maxpts,xf,yf,zf, &
            xray,yray,zray,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l,r_inner,n_grids, &
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
            call trace(bx,by,bz,nx,ny,nz,box,add_dip, &
            xi,yi,zi,dir,maxpts,xf1,yf1,zf1, &
            xray1,yray1,zray1,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l1,r_inner,n_grids, &
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
            else  if((rad2.le.r_inner).and.(rad3.le.r_inner))then
                call gsplci(5)
                call curve3(xray,yray,zray,l)
                call curve3(xray1,yray1,zray1,l1)
            else if((rad1.le.r_inner).or.(rad2.le.r_inner) &
                .or.(rad3.le.r_inner))then
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
        call rungem(bx,by,bz,nx,ny,nz,box,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,r_inner, &
        add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
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
        call rungem(bx,by,bz,nx,ny,nz,box,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,r_inner, &
        add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
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
!
!
!	****************************************
!
!
subroutine conplane(stuff,cross,nx,ny,nz,box, &
    xcraft,ncraft,re_equiv,time,label,start,frac,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot x-y quantities at the
    !       spacecraft position
    !       percent - percentage level of max value for isosurface
    !        xcraft is the position of the spacecraft: assumed to
    !                  be in simulation units
    !        ncraft is the number of spacecraft to plot
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
    common /recint/ irecmj,irecmn,irectx
    !
    real xrays(2),yrays(2)
    dimension xcraft(4,ncraft),stuff(nx,ny,nz),cross(nx,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    character*4 wd1,wd2,wd3
    character*12 label
    character*20 title
    !
    !      dimension for plotted array
    !
    my=nx
    mz=ny
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
            k1=ak
            k2=ak+1
            dz=ak-k1
            ddz=1.-dz
            fmin=0.
            fmax=0.
            do j=1,ny
                do i=1,nx
                    cross(i,j)=stuff(i,j,k1)*ddz+ &
                    stuff(i,j,k2)*dz
                    fmin=amin1(fmin,cross(i,j))
                    fmax=amax1(fmax,cross(i,j))
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
                cmax=fmax
                cmin=fmin
            endif
            cmin=cmin/frac
            cmax=cmax/frac
            clev=(cmax-cmin)/(ncon+1.)
            if(.not.start) call cnrccf(cross,my,my,mz,cmin, &
                cmax,clev,0,-1,-1012,2.5,0,1)
            !     irecmn=18
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
                xrays(1)=ai-.5
                xrays(2)=ai+.5
                yrays(1)=aj+shift
                yrays(2)=aj+shift
                call curve(xrays,yrays,2)
                xrays(1)=ai+shift
                xrays(2)=ai+shift
                yrays(1)=aj+.5
                yrays(2)=aj-.5
                call curve(xrays,yrays,2)
            enddo
        endif
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine conplane_fix(stuff,cross,nx,ny,nz,box, &
    xcraft,ncraft,re_equiv,time,label,start,alo,ahi,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot x-y quantities at the
    !       spacecraft position
    !       percent - percentage level of max value for isosurface
    !        xcraft is the position of the spacecraft: assumed to
    !                  be in simulation units
    !        ncraft is the number of spacecraft to plot
    !
    logical start
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj,irecmn,irectx
    !
    real xrays(2),yrays(2)
    dimension xcraft(4,ncraft),stuff(nx,ny,nz),cross(nx,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    character*4 wd1,wd2,wd3
    character*12 label
    character*20 title
    !
    !      dimension for plotted array
    !
    my=nx
    mz=ny
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
        k1=ak
        k2=ak+1
        dz=ak-k1
        ddz=1.-dz
        fmin=0.
        fmax=0.
        do j=1,ny
            do i=1,nx
                cross(i,j)=stuff(i,j,k1)*ddz+ &
                stuff(i,j,k2)*dz
                fmin=amin1(fmin,cross(i,j))
                fmax=amax1(fmax,cross(i,j))
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
        !     irecmn=18
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
            xrays(1)=ai-.5
            xrays(2)=ai+.5
            yrays(1)=aj+shift
            yrays(2)=aj+shift
            call curve(xrays,yrays,2)
            xrays(1)=ai+shift
            xrays(2)=ai+shift
            yrays(1)=aj+.5
            yrays(2)=aj-.5
            call curve(xrays,yrays,2)
        enddo
    endif
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine contop(stuff,vx,vy,vz,bx,by,bz,nx,ny,nz, &
    n_grids,mm,box,xcraft,ncraft,re_equiv,r_inner, &
    xmin,xmax,ymin,ymax,zmin,zmax,time, &
    label,nlevs,ncon,strtch,add_dip,ivel,kht,start, &
    tx,ty,tz,tt,t,t2,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      subroutine to plot nlevs isosurfaces and then plot contours along
    !        x-y and x-z planes
    !        nlevs number of isosurfaces to be plotted	- max 4
    !        ncon number of contours to be plotted		- max 14
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    common /recint/ irecmj     ,irecmn     ,irectx
    !
    dimension t(mx,my,mz),tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz), &
    tt(mx,my,mz),t2(mx,my,mz2),work(muvwp2,muvwp2)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    real xray(1000),yray(1000),zray(1000), &
    xray1(1000),yray1(1000),zray1(1000)
    dimension stuff(nx,ny,nz,n_grids),vx(nx,ny,nz),xcraft(4,ncraft), &
    vy(nx,ny,nz),vz(nx,ny,nz), &
    bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    real xrays(2),yrays(2)
    real eye(3),tlev(6),tcon(14)
    character*4 llbs(14)
    character*5 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    integer box
    logical add_dip,start
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
    !
    ilab=0
    ioffm=10
    !
    !     skip parameters for field lines and arrows
    !
    !      jskip=(ny-1)/20
    !      iskip=(nx-1)/30+1
    jskip=4
    iskip=4
    !
    !      dimension for 3-d plotted array
    !
    my2=my/2+1
    !
    maxpts=1000
    rx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    !
    !      effective step size between points - slightly less than
    !         unity to ensure arrows do not go out of bounds
    !
    hi=.99
    hj=.99
    hk=.99
    !
    !     al = length of arrow in xy-plane
    !     stretch = enlargement of arrows
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
    axmax=amin1(xmax,grd_xmax(box)-.00001)
    aymax=amin1(ymax,grd_ymax(box)-.00001)
    azmax=amin1(zmax,grd_zmax(box)-.00001)
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
    call filldatav(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
    t,tt,tx,ty,tz,mx,my,mz,xmin,ymin,zmin,delx,dely,delz,vm, &
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
    tma=0.8*tma
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
    !     initialize eye position
    !
    eye(1)=mx/2.
    eye(2)=my/2.
    eye(3)=mz*10.
    !
    !     initialize viewport and frame headings
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    !
    !     title=' : dawn dusk'
    call wtstr(.4,.975,label,2,0,0)
    !     call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    call gsplci(2)
    tisom=r_inner
    call isosrf(tt,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
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
    ijump=-1
    ncol=3
    kk=2
    k=mz2+kht
    !     jskip=2
    !     iskip=4
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
                    !
                    !       draw magnetic field line
                    !
                endif
            endif
            l=0
            dir=0.25*delx
            xi=xmin+delx*(i-1)
            yi=ymin+dely*(j-1)
            zi=zmin+delz*(k-1)
            call trace(bx,by,bz,nx,ny,nz,box,add_dip, &
            xi,yi,zi,dir,maxpts,xf,yf,zf, &
            xray,yray,zray,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l,r_inner,n_grids, &
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
            call trace(bx,by,bz,nx,ny,nz,box,add_dip, &
            xi,yi,zi,dir,maxpts,xf1,yf1,zf1, &
            xray1,yray1,zray1,xmin,axmax,ymin,aymax, &
            zmin+height/2.,azmax-height/2.,l1,r_inner,n_grids, &
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
            else  if((rad2.le.r_inner).and.(rad3.le.r_inner))then
                call gsplci(5)
                call curve3(xray,yray,zray,l)
                call curve3(xray1,yray1,zray1,l1)
            else if((rad1.le.r_inner).or.(rad2.le.r_inner) &
                .or.(rad3.le.r_inner))then
                ijump=-ijump
                if(ijump.ge.0)then 
                    call gsplci(7)
                    call curve3(xray,yray,zray,l)
                    call curve3(xray1,yray1,zray1,l1)
                endif
            else
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
    !       wspace set writing distance to mlt and lat
	!
    wspace=0.015
    nspace=0
    size=0.5*delx
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
            x1=1.+(ax-size-xmin)/delx
            x2=1.+(ax+size-xmin)/delx
            y1=1.+(ay-size-ymin)/dely
            y2=1.+(ay+size-ymin)/dely
            z1=1.+(az+height-size-zmin)/delz
            z2=1.+(az+height+size-zmin)/delz
            r1=1.+(ar-size-ymin)/dely
            r2=1.+(ar+size-ymin)/dely
            !
            !        2-d projection
            !
            !       call sflush
            !       call gsplci(15)
            !        call line3(x1,y1,1.,x2,y1,1.)
            !        call line3(x1,y1,1.,x1,y2,1.)
            !        call line3(x2,y2,1.,x1,y2,1.)
            !        call line3(x2,y2,1.,x2,y1,1.)
            !
            !        3-d position
            !
            call sflush
            call gsplci(9)
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
        dir=-0.25*rx
        call rungem(bx,by,bz,nx,ny,nz,box,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,r_inner, &
        add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        test for connection to earth and write out position
        !
        call gsplci(1)
        sx=(xray(npts)-xdip)*re_equiv
        sy=(yray(npts)-ydip)*re_equiv
        sz=(zray(npts)-zdip)*re_equiv
        rtot=sqrt(sx**2+sy**2+sz**2)
        if((rtot/re_equiv.le.r_inner+2.).and. &
            (rtot/re_equiv.gt.0.5*r_inner))then
            nspace=nspace+1
            !
            !        find polar coordinate : first rotate by dipole tilt
            !
            sx1=sx*cos_tilt-sz*sin_tilt
            sz1=sx*sin_tilt+sz*cos_tilt
            sy1=sy
            rtot=sqrt(sx1**2+sy1**2+sz1**2)
            if(sz1.ge.0)then
                colat=acos(sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=90.-colat
            else
                colat=acos(-sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=-(90.-colat)
            endif
            !
            !          find local time
            !
            if(sy1.ge.0.)then
                alt=atan2(sy1,sx1)
                alt=12.*alt/3.1416
            else
                alt=atan2(-sy1,sx1)
                alt=12.*(2.-alt/3.1416)
            endif
            arad=1.
        endif
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
        !
        !       go the other direction
        dir=+0.25*rx
        call rungem(bx,by,bz,nx,ny,nz,box,rx, &
        ax,ay,az,xmin,xmax,ymin,ymax, &
        zmin+height/2.,zmax-height/2.,r_inner, &
        add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        sx=(xray(npts)-xdip)*re_equiv
        sy=(yray(npts)-ydip)*re_equiv
        sz=(zray(npts)-zdip)*re_equiv
        rtot=sqrt(sx**2+sy**2+sz**2)
        if((rtot/re_equiv.le.r_inner+2.).and. &
            (rtot/re_equiv.gt.0.5*r_inner))then
            nspace=nspace+1
            !
            !        find polar coordinate : first rotate by dipole tilt
            !
            sx1=sx*cos_tilt-sz*sin_tilt
            sz1=sx*sin_tilt+sz*cos_tilt
            !
            sy1=sy
            rtot=sqrt(sx1**2+sy1**2+sz1**2)
            if(sz1.ge.0)then
                colat=acos(sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=90.-colat
            else
                colat=acos(-sz1/rtot)
                colat=57.296*colat/(sqrt(rtot))
                alat=-(90.-colat)
            endif
            !
            !          find local time
            !
            if(sy1.gt.0.)then
                alt=atan2(sy1,sx1)
                alt=12.*alt/3.1416
            else
                alt=atan2(-sy1,sx1)
                alt=12.*(2.-alt/3.1416)
            endif
            arad=1.
        endif
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
    enddo
    !
    return
end
!
!
!	****************************************
!
!
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
    integer box
    logical add_dip,roc
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
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
!
!
!	****************************************
!
!
subroutine contur(stuff,nx,ny,nz,n_grids,box,xmin,xmax, &
    ymin,ymax,zmin,zmax,time,label,nlevs,ncon,r_inner, &
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
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension stuff(nx,ny,nz,n_grids)
    !
    real eye(3),tlev(6),tcon(14)
    character*4 llbs(14),wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    integer vpl,vpr,vpb,vpt,lind(14)
    integer box
	integer, parameter :: dp = kind(1.d0)
	real(dp) time
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
    axmax=amin1(xmax,grd_xmax(box)-.00001)
    aymax=amin1(ymax,grd_ymax(box)-.00001)
    azmax=amin1(zmax,grd_zmax(box)-.00001)
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
    ddx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    ddy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    ddz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    !      load t stuff
    !
    do k=1,mz
        az=zmin+delz*(k-1)
        ak=1.+(az-grd_zmin(box))/ddz
        k1=ak
        k2=k1+1
        dz=ak-k1
        do j=1,my
            ay=ymin+dely*(j-1)
            aj=1.+(ay-grd_ymin(box))/ddy
            j1=aj
            j2=j1+1
            dy=aj-j1
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
    tisom=r_inner
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
    tisom=r_inner
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
    tisom=r_inner
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
!
!
!	****************************************
!
!
subroutine convect(vx,vy,vz,nx,ny,nz,box,radstrt, &
    re_equiv,iside,time,label,write_dat,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      this subroutine will try to make convection pattern
    !         over auroral oval
    !      sdata tdata included to allow passing of data between
    !         successive calls of aurora
	!
    common /space/sdata(61,61),tdata(61,61),work(31,31), &
    avx(31,31),avy(31,31)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    integer box
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz),spv(2)
    character*4 wd1,wd2,wd3
    character*12 label,magnif
    character*20 title
    logical write_dat
	!
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !     for no line labeling set ilab  to zero
    !
    common /conre4/ isizel,isizem,isizep,nrep, &
    ncrt,ilab,nulbll,ioffd, & 
    ext,ioffm,isolid ,nla, &
    nlm,xlt,ybt,side 
    !
    ilab=0
    !
    !      dimension for plotted array
    !
    mx=31
    my=31
    my2=my/2+1
    mx2=mx/2+1
    !
    !     re_equiv is the physical distance of a grid unit in re's
    !      theta_range is the latitudes for the polar plot and is set
    !                    at 35 degrees
    !       theta_equiv converts latitude at the earth surface to
    !                   latitude at radstrt
    !
    !       re_equiv=0.84
	!
    theta_equiv=sqrt(re_equiv*radstrt)
    theta_range=0.6108652
    degrees=57.3*theta_range
    del_theta=theta_range/float(mx2-1)
    !
    !      initialize work array
    !
    do j=1,my
        do i=1,mx
            work(i,j)=0.0
            avx(i,j)=0.0
            avy(i,j)=0.0
        enddo
    enddo
    !
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    do j=1,my
        do i=1,mx
            !
            !         find equivalent latitude
            !
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            dlat=amax1(dlat,0.0000001)
            alat=dlat*del_theta*theta_equiv
            if(alat.le.1.55) then
                cost=cos(alat)
                sint=sin(alat)
                !
                !         find equivalent longitude
                !
                cosp=(i-mx2)/dlat
                sinp=(j-my2)/dlat
                !
                !          find position on grid
                !
                z1=iside*(radstrt*cost)
                x1=radstrt*sint*cosp
                !
                !       dipole space to real space
                !
                ax=x1*cos_tilt+z1*sin_tilt+xdip
                az=-x1*sin_tilt+z1*cos_tilt+zdip
                ay=(radstrt*sint*sinp+ydip)
                ar=sqrt((ax-xdip)**2+(ay-ydip)**2+(az-zdip)**2)+0.0000001
                !
                !          interpolate data to grid point
                !
                ak=1.+(az-grd_zmin(box))/delz
                k1=ak
                k2=k1+1
                dz=ak-k1
                !
                aj=1.+(ay-grd_ymin(box))/dely
                j1=aj
                j2=j1+1
                dy=aj-j1
                !
                ai=1.+(ax-grd_xmin(box))/delx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                bvx=vx(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +vx(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +vx(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +vx(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +vx(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +vx(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +vx(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +vx(i2,j2,k2)*(dx)*(dy)*(dz)
                bvy=vy(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +vy(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +vy(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +vy(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +vy(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +vy(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +vy(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +vy(i2,j2,k2)*(dx)*(dy)*(dz)
                bvz=vz(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
                +vz(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
                +vz(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
                +vz(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
                +vz(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
                +vz(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
                +vz(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
                +vz(i2,j2,k2)*(dx)*(dy)*(dz)
                !
                !          find tangential components
                !
                cost=iside*(az-zdip)/ar
                bvz=iside*bvz
                sint=sqrt((ax-xdip)**2+(ay-ydip)**2)/ar
                cosp=(ax-xdip)/(ar*sint+0.000001)
                sinp=(ay-ydip)/(ar*sint+0.000001)
                vr=bvx*sint*cosp+bvy*sint*sinp+bvz*cost
                vtheta=bvx*cost*cosp+bvy*cost*sinp-bvz*sint
                vphi=-bvx*sinp+bvy*cosp
                !
                !        restructure magnetic field without radial velocity component
                !
                vr=0.
                bvx=vr*sint*cosp+vtheta*cost*cosp-vphi*sinp
                bvy=vr*sint*sinp+vtheta*cost*sinp+vphi*cosp
                bvz=vr*cost-vtheta*sint
                avx(i,j)=vtheta*cosp-vphi*sinp
                avy(i,j)=vtheta*sinp+vphi*cosp
                !         avx(i,j)=bvx
                !         avy(i,j)=bvy
            endif
            !
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
    !     set color table
    !
    call sflush
    call isorb
    !
    title=' : convect'
    call wtstr(.4,.975,label,2,0,0)
    call wtstr(.6,.975,title,2,0,0)
    write(title,'(f7.3)')time
    title='ut = '//title
    call wtstr(.8,.975,title,2,0,0)
    !
    title='lat. '
    call wtstr(.15,.96,title,1,0,0)
    write(title,'(f7.3)')degrees
    call wtstr(.2,.96,title,1,0,0)
    !
    xcmin=0.
    xcmax=0.
    ycmin=0.
    ycmax=0.
    do j=1,my
        do i=1,mx
            xcmav=amax1(avx(i,j),xcmav)
            xcmin=amin1(avx(i,j),xcmin)
            ycmav=amax1(avy(i,j),ycmav)
            ycmin=amin1(avy(i,j),ycmin)
        enddo
    enddo
    !
    if(abs(xcmin).gt.xcmax)xcmax=abs(xcmin)
    if(abs(ycmin).gt.ycmax)ycmax=abs(ycmin)
    cmax=amax1(xcmax,ycmax)
    if(cmax.gt.1.e-4)then
        title='maxspeed '
        call wtstr(.3,.96,title,1,0,0)
        write(title,'(f7.4)')cmax
        call wtstr(.4,.96,title,1,0,0)
        !      call ezvec(avx,avy,mx,my)
        !      call velvct(avx,mx,avy,mx,mx,my,0.1*cmax,0.,0,100,0,spv)
        call vvinit(avx,mx,avy,mx,dummy,0,mx,my,dummy,0)
        call vvsetr('vrl',0.20)
        call vvectr(avx,avy,dummy,dummy,dummy,dummy)
    endif
    !
    !     make constant latitude circles
    !
    dd=degrees/float(mx2-1)
    fmax=0.
    do j=1,my
        do i=1,mx
            dlat=sqrt(1.*(i-mx2)**2+1.*(j-my2)**2)
            work(i,j)=90.-dlat*dd
            fmax=amax1(fmax,work(i,j))
        enddo
    enddo
    !
    !     ilab=1
    call conrec(work,mx,mx,my,60.,90.,10.,0,-1,-1012)
    !
    !      output data to data file if necessary
    !
    if(write_dat)then
        write(recdt_f,*) time,avx,avy
    endif
    return
end
!
!
!	****************************************
!
!
subroutine plot1d(xary,yary,npt,time,label)
    !
    dimension xary(npt),yary(npt),xmax(2),ymax(2)
    character*12 label
    character*20 title
    !
    fmax=yary(1)
    fmin=fmax
    do i=2,npt
        fmax=amax1(yary(i),fmax)
        fmin=amin1(yary(i),fmin)
    enddo
    !
    if((fmax.gt.0.0).and.(fmin.gt.0.0))fmin=0.0
    if((fmax.lt.0.0).and.(fmin.lt.0.0))fmax=0.0
    !
    xmax(2)=xary(npt)
    xmax(1)=xary(1)
    ymax(2)=fmax
    ymax(1)=fmin
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    write(title,'(f7.3)')time
    title='t= '//title
    call wtstr(.75,.97,title,1,0,0)
    call ezxy(xmax,ymax,2,label)
    call curve(xary,yary,npt)
    call agseti('set.',1)
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    !
    return
end
!
!
!	****************************************
!
!
subroutine rungea(bx,by,bz,press,rho,nx,ny,nz, &
    box,xi,yi,zi,xmin,xmax,ymin,ymax, &
    zmin,zmax,add_dip,ergies,rhod,maxpts,dir, &
    delx,dely,delz,r_inner,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine uses a fourth order runge-kutta method
    !         to trace stream functions of current and magnetic field
    !         xi,yi,zi, and find maximum temperature along them
    !
    integer box
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz), &
    press(nx,ny,nz),rho(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    real xray(1000),yray(1000),zray(1000)
    !
    logical add_dip,roc
    !
    maxpts=1000
    !
    tstep=0.25*dir
    !
    xray(1)=xi
    yray(1)=yi
    zray(1)=zi
    !
    ts6=tstep/6.
    do n=2,maxpts
        !
        !           step 1
        !
        ax=xray(n-1)
        ay=yray(n-1)
        az=zray(n-1)
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx1,dy1,dz1,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 2
        !
        ax=xray(n-1)+tstep*0.5*dx1
        ay=yray(n-1)+tstep*0.5*dy1
        az=zray(n-1)+tstep*0.5*dz1
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx2,dy2,dz2,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 3
        !
        ax=xray(n-1)+tstep*0.5*dx2
        ay=yray(n-1)+tstep*0.5*dy2
        az=zray(n-1)+tstep*0.5*dz2
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx3,dy3,dz3,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 4
        !
        ax=xray(n-1)+tstep*dx3
        ay=yray(n-1)+tstep*dy3
        az=zray(n-1)+tstep*dz3
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx4,dy4,dz4,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       update new position and wavevector
        !
        xray(n)=xray(n-1)+ts6*(dx1+2.*dx2+2.*dx3+dx4)
        yray(n)=yray(n-1)+ts6*(dy1+2.*dy2+2.*dy3+dy4)
        zray(n)=zray(n-1)+ts6*(dz1+2.*dz2+2.*dz3+dz4)
        !
        !    find energy at the point
        !
        ax=xray(n)
        ay=yray(n)
        az=zray(n)
        !
        !          interpolate data to grid point
        !
        ak=1.+(az-grd_zmin(box))/delz
        k1=ak
        k2=k1+1
        dz=ak-k1
        !
        aj=1.+(ay-grd_ymin(box))/dely
        j1=aj
        j2=j1+1
        dy=aj-j1
        !
        ai=1.+(ax-grd_xmin(box))/delx
        i1=ai
        i2=i1+1
        dx=ai-i1
        !
        aerg=press(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
        +press(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
        +press(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
        +press(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
        +press(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
        +press(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
        +press(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
        +press(i2,j2,k2)*(dx)*(dy)*(dz)
        arho=rho(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
        +rho(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
        +rho(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
        +rho(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
        +rho(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
        +rho(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
        +rho(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
        +rho(i2,j2,k2)*(dx)*(dy)*(dz)
        !
        aergd=aerg/arho
        arho=amax1(arho,0.0001)
        !
        ergies=amax1(aergd,ergies)
        rhod=amin1(arho,rhod)
        !
        !       test to see if ray is within selected region
        !
        ar=sqrt((xray(n)-xdip)**2+(yray(n)-ydip)**2 &
        +(zray(n)-zdip)**2)
        radius=sqrt((xray(n)-xray(1))**2+(yray(n)-yray(1))**2 &
        +(zray(n)-zray(1))**2)
        if((xray(n).lt.xmin).or.(xray(n).gt.xmax).or. &
        (yray(n).lt.ymin).or.(yray(n).gt.ymax).or. &
        (zray(n).lt.zmin).or.(zray(n).gt.zmax).or. &
        (ar.lt.r_inner+1.).or.(radius.gt.4.*r_inner).or.(roc)) return
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
subroutine rungeb(efldx,efldy,efldz,bx,by,bz,nx,ny,nz, &
    n_grids,box,rx,xi,yi,zi,xmin,xmax,ymin,ymax, &
    zmin,zmax,add_dip,potential,maxpts,dir,r_inner, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine uses a fourth order runge-kutta method
    !         to trace stream functions of current and magnetic field
    !         xi,yi,zi, and integrates the field aligned potential drop
    !
    integer box
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz)
    real xray(1000),yray(1000),zray(1000)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    logical add_dip,roc
    !
    rx2=rx*3.
    potential=0.
    maxpts=1000
    delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    tstep=0.25*dir
    !
    xray(1)=xi
    yray(1)=yi
    zray(1)=zi
    !
    ts6=tstep/6.
    do n=2,maxpts
        !
        !           step 1
        !
        ax=xray(n-1)
        ay=yray(n-1)
        az=zray(n-1)
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx1,dy1,dz1,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 2
        !
        ax=xray(n-1)+tstep*0.5*dx1
        ay=yray(n-1)+tstep*0.5*dy1
        az=zray(n-1)+tstep*0.5*dz1
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx2,dy2,dz2,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 3
        !
        ax=xray(n-1)+tstep*0.5*dx2
        ay=yray(n-1)+tstep*0.5*dy2
        az=zray(n-1)+tstep*0.5*dz2
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx3,dy3,dz3,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 4
        !
        ax=xray(n-1)+tstep*dx3
        ay=yray(n-1)+tstep*dy3
        az=zray(n-1)+tstep*dz3
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx4,dy4,dz4,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       update new position and wavevector
        !
        xray(n)=xray(n-1)+ts6*(dx1+2.*dx2+2.*dx3+dx4)
        yray(n)=yray(n-1)+ts6*(dy1+2.*dy2+2.*dy3+dy4)
        zray(n)=zray(n-1)+ts6*(dz1+2.*dz2+2.*dz3+dz4)
        !
        !      define integral limits
        !
        xi=xray(n-1)
        yi=yray(n-1)
        zi=zray(n-1)
        !
        xf=xray(n)
        yf=yray(n)
        zf=zray(n)
        !
        !          interpolate data to grid point
        !
        ax=(xf+xi)/2.
        ay=(yf+yi)/2.
        az=(zf+zi)/2.
        !
        ak=1.+(az-grd_zmin(box))/delz
        k1=ak
        k2=k1+1
        dz=ak-k1
        ddz=1.-dz
        !
        aj=1.+(ay-grd_ymin(box))/dely
        j1=aj
        j2=j1+1
        dy=aj-j1
        ddy=1.-dy
        !
        ai=1.+(ax-grd_xmin(box))/delx
        i1=ai
        i2=i1+1
        dx=ai-i1
        ddx=1.-dx
        !
        aex=  efldx(i1,j1,k1)*ddx*ddy*ddz &
        +efldx(i1,j1,k2)*ddx*ddy*dz &
        +efldx(i1,j2,k1)*ddx*dy*ddz &
        +efldx(i1,j2,k2)*ddx*dy*dz &
        +efldx(i2,j1,k1)*dx*ddy*ddz &
        +efldx(i2,j1,k2)*dx*ddy*dz &
        +efldx(i2,j2,k1)*dx*dy*ddz &
        +efldx(i2,j2,k2)*dx*dy*dz
        aey=  efldy(i1,j1,k1)*ddx*ddy*ddz &
        +efldy(i1,j1,k2)*ddx*ddy*dz &
        +efldy(i1,j2,k1)*ddx*dy*ddz &
        +efldy(i1,j2,k2)*ddx*dy*dz &
        +efldy(i2,j1,k1)*dx*ddy*ddz &
        +efldy(i2,j1,k2)*dx*ddy*dz &
        +efldy(i2,j2,k1)*dx*dy*ddz &
        +efldy(i2,j2,k2)*dx*dy*dz
        aez=  efldz(i1,j1,k1)*ddx*ddy*ddz &
        +efldz(i1,j1,k2)*ddx*ddy*dz &
        +efldz(i1,j2,k1)*ddx*dy*ddz &
        +efldz(i1,j2,k2)*ddx*dy*dz &
        +efldz(i2,j1,k1)*dx*ddy*ddz &
        +efldz(i2,j1,k2)*dx*ddy*dz &
        +efldz(i2,j2,k1)*dx*dy*ddz &
        +efldz(i2,j2,k2)*dx*dy*dz
        !
        drx=xf-xi
        dry=yf-yi
        drz=zf-zi
        !
        potential=potential+drx*aex+dry*aey+drz*aez
        !
        !       test to see if ray is within selected region
        !
        ar=sqrt(xray(n)**2+yray(n)**2 +zray(n)**2)
    	!
        radius=sqrt((xray(n)-xray(1))**2+(yray(n)-yray(1))**2 &
        +(zray(n)-zray(1))**2)
        if((xray(n).lt.xmin).or.(xray(n).gt.xmax).or. &
        (yray(n).lt.ymin).or.(yray(n).gt.ymax).or. &
        (zray(n).lt.zmin).or.(zray(n).gt.zmax).or. &
        (ar.lt.r_inner+rx2).or.(radius.gt.3.*r_inner) &
        .or.(roc)) return
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
subroutine rungem(bx,by,bz,nx,ny,nz,box,rx, &
    xi,yi,zi,xmin,xmax,ymin,ymax,zmin,zmax,r_inner, &
    add_dip,xray,yray,zray,maxpts,npts,dir,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine uses a fourth order runge-kutta method
    !         to trace stream functions of current and magnetic field
    !         xi,yi,zi, until it hits boundary or maximum number
    !         of points is reached
    !
    integer box
    dimension bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)
    dimension xray(maxpts),yray(maxpts),zray(maxpts)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    !
    logical add_dip,roc
    !
    !      maxpts=1000
    !
    tstep=0.1*dir
    !
    xray(1)=xi
    yray(1)=yi
    zray(1)=zi
    !
    ts6=tstep/6.
    do n=2,maxpts
        !
        !           step 1
        !
        ax=xray(n-1)
        ay=yray(n-1)
        az=zray(n-1)
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx1,dy1,dz1,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 2
        !
        ax=xray(n-1)+tstep*0.5*dx1
        ay=yray(n-1)+tstep*0.5*dy1
        az=zray(n-1)+tstep*0.5*dz1
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx2,dy2,dz2,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 3
        !
        ax=xray(n-1)+tstep*0.5*dx2
        ay=yray(n-1)+tstep*0.5*dy2
        az=zray(n-1)+tstep*0.5*dz2
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx3,dy3,dz3,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !        step 4
        !
        ax=xray(n-1)+tstep*dx3
        ay=yray(n-1)+tstep*dy3
        az=zray(n-1)+tstep*dz3
        call intpol(bx,by,bz,nx,ny,nz,box,add_dip, &
        ax,ay,az,dx4,dy4,dz4,1.,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       update new position and wavevector
        !
        npts=n
        xray(n)=xray(n-1)+ts6*(dx1+2.*dx2+2.*dx3+dx4)
        yray(n)=yray(n-1)+ts6*(dy1+2.*dy2+2.*dy3+dy4)
        zray(n)=zray(n-1)+ts6*(dz1+2.*dz2+2.*dz3+dz4)
        !
        !       test to see if ray is within selected region
        !
        ar=sqrt(xray(n)**2+yray(n) **2+zray(n)**2)
        if((xray(n).lt.xmin).or.(xray(n).gt.xmax).or. &
            (yray(n).lt.ymin).or.(yray(n).gt.ymax).or. &
            (zray(n).lt.zmin).or.(zray(n).gt.zmax).or. &
            (ar.lt.r_inner+1.5*rx).or.(roc)) then
            !
            npts=npts-1
            return
        endif
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine step(bx,by,bz,n1,n2,n3,box,add_dip,x,y,z,ds,errin,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     re-calculates coords x,y,z for one step along field line.
    !     ds is step size,
    !     errin is permissible error value.
    !
    integer box
    dimension bx(n1,n2,n3),by(n1,n2,n3),bz(n1,n2,n3)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    logical roc,add_dip
    !
    do 
        ds3=-ds/3.
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x,y,z,r11,r12,r13,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        if(roc) return
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+r11,y+r12,z+r13,r21,r22,r23,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+.5*(r11+r21),y+.5*(r12+r22),z+.5* &
        (r13+r23),r31,r32,r33,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+.375*(r11+3.*r31),y+.375*(r12+3.*r32 &
        ),z+.375*(r13+3.*r33),r41,r42,r43,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+1.5*(r11-3.*r31+4.*r41),y+1.5*(r12- &
        3.*r32+4.*r42),z+1.5*(r13-3.*r33+4.*r43), &
        r51,r52,r53,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        errcur=abs(r11-4.5*r31+4.*r41-.5*r51)+abs(r12-4.5*r32+4.*r42-.5* &
        r52)+abs(r13-4.5*r33+4.*r43-.5*r53)
        if (errcur.lt.errin) exit
        ds=ds*.5
        if (ds.lt.0.2) exit
    enddo
    x=x+.5*(r11+4.*r41+r51)
    y=y+.5*(r12+4.*r42+r52)
    z=z+.5*(r13+4.*r43+r53)
    if(errcur.lt.errin*.04.and.abs(ds).lt.1.33) ds=ds*1.5
    return
end
!
!
!	****************************************
!
!
subroutine trace(bx,by,bz,n1,n2,n3,box,add_dip, &
    xi,yi,zi,dir,np,xf,yf,zf,xx,yy,zz, &
    xmin,xmax,ymin,ymax,zmin,zmax,l,r_inner,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !    3 field line tracing subroutines originally by tsyganenko.
    !    modified to trace field lines of bx,by,bz (nx,ny,nz)
    !
    !   traces field line from arbitrary point of space up to the earth
    !   surface or up to model limiting boundary.
    !-------------- input parameters:
    !   xi,yi,zi - gsm coords of initial point (in earth radii),
    !   dir - sign of tracing direction: if dir=1. then antiparallel to
    !     b vector (e.g. from northern to southern conjugate point),
    !     and if dir=-1. then parallel to b.
    !   np - upper estimate of number of steps along the field line
    !     (of the order of several hundreds).
    !-------------- output parameters:
    !   xf,yf,zf - gsm coords of final point
    !   xx,yy,zz - arrays (length np) containing coords of field line points
    !   l - actual number of field line points. if l exceeds np, tracing
    !     terminates, and a warning is displayed
    !
    !                   author: nikolai a. tsyganenko
    !                           institute of physics
    !                           st.-petersburg state university
    !                           stary petergof 198904
    !                           st.-petersburg
    !                           russia
    !
    integer box
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    dimension bx(n1,n2,n3),by(n1,n2,n3),bz(n1,n2,n3), &
    xx(np),yy(np),zz(np)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    logical roc,add_dip
    !
    err=0.05
    !     ds=0.5*dir
    x=xi
    y=yi
    z=zi
    ds=dir
    ds3=ds/3.
    l=0
    !
    call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
    x,y,z,r1,r2,r3,ds3,roc,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    if(roc) return
    !
    do l=1,np
        xx(l)=x
        yy(l)=y
        zz(l)=z
        !
        call step(bx,by,bz,n1,n2,n3,box,add_dip,x,y,z,ds,err,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        ar=sqrt((x-xdip)**2+(y-ydip)**2+(z-zdip)**2)
        if((x.lt.xmin).or.(x.gt.xmax).or. &
            (y.lt.ymin).or.(y.gt.ymax).or. &
            (z.lt.zmin).or.(z.gt.zmax).or. &
            (ar.lt.r_inner+1.)) then
            xf=x
            yf=y
            zf=z
            xx(l)=x
            yy(l)=y
            zz(l)=z
            err=0.0005
            return
        endif
    	!
    enddo
    !
    !     not enough points
    !
    l=l-1
    return
end
!
!
!	****************************************
!
!
subroutine vct3d(vx,vy,vz,nx,ny,nz,box,xmin,xmax,ymin, &
    ymax,zmin,zmax,time,label,iv,strtch,r_inner, &
    t,tx,ty,tz,work,mx,my,mz,muvwp2,mz2,n_grids, &
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
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    real eye(3)
    integer vpl,vpr,vpb,vpt
    integer box
    character*4 wd1,wd2,wd3
    character*12 label
    character*15 title
    !
    !      dimension for plotted array
	!
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
    axmax=amin1(xmax,grd_xmax(box)-.00001)
    aymax=amin1(ymax,grd_ymax(box)-.00001)
    azmax=amin1(zmax,grd_zmax(box)-.00001)
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
    ddx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    ddy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    ddz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
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
    tisom=r_inner
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
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
    call sflush
    !
    !     determine vector field
    !
    vm=0
    do k=1,mz
        az=zmin+(k-1.)*delz
        ak=1.+(az-grd_zmin(box))/ddz
        k1=ak
        k2=k1+1
        dz=ak-k1
        !
        do j=1,my
            !
            ay=ymin+(j-1.)*dely
            aj=1.+(ay-grd_ymin(box))/ddy
            j1=aj
            j2=j1+1
            dy=aj-j1
            !
            do i=1,mx
                ax=xmin+(i-1.)*delx
                ai=1.+(ax-grd_xmin(box))/ddx
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
                 if(radius.lt.r_inner+1.)then
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
    tisom=r_inner
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
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
    call sflush
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
    tisom=r_inner
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
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
    call sflush
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
    !     top view -- dusk side
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
    tisom=r_inner
    call isosrf(t,mx,mx,my,my,mz,eye,muvwp2,work,tisom,-3, &
    vpl,vpr,vpb,vpt)
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
    call sflush
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
!
!
!	************************************
!
!
subroutine filldatav(stuff,vx,vy,vz,nx,ny,nz,n_grids,mm,box, &
    t,tt,tx,ty,tz,mx,my,mz,xmin,ymin,zmin, &
    delx,dely,delz,vm,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	!
	!	Fill in diagnostic data cubes
	!    
    integer box
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension stuff(nx,ny,nz,n_grids),vx(nx,ny,nz), &
    vy(nx,ny,nz),vz(nx,ny,nz)
    dimension t(mx,my,mz),tt(mx,my,mz), &
    tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz)
    !
    !      load t stuff
    !
    ddx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    ddy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    ddz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    vm=0.00001
    do k=1,mz
        az=zmin+delz*(k-1)
        ak=1.+(az-grd_zmin(box))/ddz
        k1=ak
        k2=k1+1
        dz=ak-k1
        !
        do j=1,my
            ay=ymin+dely*(j-1)
            aj=1.+(ay-grd_ymin(box))/ddy
            j1=aj
            j2=j1+1
            dy=aj-j1
            !
            do i=1,mx
                ax=xmin+delx*(i-1)
                ai=1.+(ax-grd_xmin(box))/ddx
                i1=ai
                i2=i1+1
                dx=ai-i1
                !
                t(i,j,k)=stuff(i1,j1,k1,mm)*(1.-dx)*(1.-dy)*(1.-dz) &
                +stuff(i1,j1,k2,mm)*(1.-dx)*(1.-dy)*(dz) &
                +stuff(i1,j2,k1,mm)*(1.-dx)*(dy)*(1.-dz) &
                +stuff(i1,j2,k2,mm)*(1.-dx)*(dy)*(dz) &
                +stuff(i2,j1,k1,mm)*(dx)*(1.-dy)*(1.-dz) &
                +stuff(i2,j1,k2,mm)*(dx)*(1.-dy)*(dz) &
                +stuff(i2,j2,k1,mm)*(dx)*(dy)*(1.-dz) &
                +stuff(i2,j2,k2,mm)*(dx)*(dy)*(dz)
                !
                radius=sqrt(ax**2+ay**2+az**2)
                tt(i,j,k)=radius
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
                !
            enddo
        enddo
    enddo
    vm=amax1(vm,0.00004)
    !
    return
end subroutine filldatav
!
!
!	*****************************************
!
!
subroutine intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
    ax,ay,az,r1,r2,r3,ds3,roc,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     interpolate the magnetic to from nearest grid points
    !     to the actual point ax,ay,az and then calculates the
    !     direction for the increments in r1,r2,r3
    !
    integer box
    dimension bx(n1,n2,n3),by(n1,n2,n3),bz(n1,n2,n3)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    logical roc,add_dip
    !
    !       altered to make plots below equator assuming symmetry
    !
    ddx=(grd_xmax(box)-grd_xmin(box))/(n1-1.)
    ddy=(grd_ymax(box)-grd_ymin(box))/(n2-1.)
    ddz=(grd_zmax(box)-grd_zmin(box))/(n3-1.)
    !
    ak=1.+(az-grd_zmin(box))/ddz
    k1=ak
    k2=k1+1
    dz=ak-k1
    !
    aj=1.+(ay-grd_ymin(box))/ddy
    j1=aj
    j2=j1+1
    dy=aj-j1
    !
    ai=1.+(ax-grd_xmin(box))/ddx
    i1=ai
    i2=i1+1
    dx=ai-i1
    !
    abx=0.
    aby=0.
    abz=0.
    !
    if(add_dip)call dipole(abx,aby,abz,ax,ay,az)
    !
    abx=abx+bx(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
    +bx(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
    +bx(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
    +bx(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
    +bx(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
    +bx(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
    +bx(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
    +bx(i2,j2,k2)*(dx)*(dy)*(dz)
    !
    aby=aby+by(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
    +by(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
    +by(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
    +by(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
    +by(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
    +by(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
    +by(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
    +by(i2,j2,k2)*(dx)*(dy)*(dz)
    !
    abz=abz+bz(i1,j1,k1)*(1. - dx)*(1. - dy)*(1. - dz) &
    +bz(i1,j1,k2)*(1. - dx)*(1. - dy)*(dz) &
    +bz(i1,j2,k1)*(1. - dx)*(dy)*(1. - dz) &
    +bz(i1,j2,k2)*(1. - dx)*(dy)*(dz) &
    +bz(i2,j1,k1)*(dx)*(1. - dy)*(1. - dz) &
    +bz(i2,j1,k2)*(dx)*(1. - dy)*(dz) &
    +bz(i2,j2,k1)*(dx)*(dy)*(1. - dz) &
    +bz(i2,j2,k2)*(dx)*(dy)*(dz)
    !
    sqb=sqrt(abx**2+aby**2+abz**2)
    if(sqb.lt.1.e-4) then
        roc=.true.
        return
    endif
    !
    roc=.false.
    b=ds3/sqb
    r1=abx*b
    r2=aby*b
    r3=abz*b
    !
    return
end subroutine intpol
