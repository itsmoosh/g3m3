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
        call totfld(bx,bx0,bsx,nx,ny,nz,n_grids,box)
        call totfld(by,by0,bsy,nx,ny,nz,n_grids,box)
        call totfld(bz,bz0,bsz,nx,ny,nz,n_grids,box)
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
        if(box.le.3)then
            preslim=20.0/float(box)
            po=preslim*0.33
        else
            preslim=20.0/(box-1.)
            po=preslim*0.33
        endif
        !
        call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,n_grids,box)
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
        call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,n_grids,box)
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
        call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,n_grids,box)
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
        call totfld(bx,bx0,bsx,nx,ny,nz,n_grids,box)
        call totfld(by,by0,bsy,nx,ny,nz,n_grids,box)
        call totfld(bz,bz0,bsz,nx,ny,nz,n_grids,box)
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
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
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
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
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
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
