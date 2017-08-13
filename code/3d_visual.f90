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
    nx,ny,nz,ngrd,xspac, &
    cross,along,flat,xcraft,ncraft,re_equiv, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,ut,b_equiv,ti_te,rho_equiv)
    !
    real grd_xmin(ngrd),grd_xmax(ngrd),grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd),xspac(ngrd)
    real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    real bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
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
    do m=ngrd,1,-1
        rx=xspac(m)
        ry=xspac(m)
        rz=xspac(m)
        !
        ymin=grd_ymin(m)
        ymax=grd_ymax(m)
        zmin=grd_zmin(m)
        zmax=grd_zmax(m)
        xmin=grd_xmin(m)
        xmax=grd_xmax(m)
        !
        xcut=((xmin+xmax))/2.
        !
        add_dip=.false.
        !
        call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
        call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
        call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
        if(m.le.3)then
            preslim=20.0/float(m)
            po=preslim*0.33
        else
            preslim=20.0/(m-1.)
            po=preslim*0.33
        endif
        !
        call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=qpresx(i,j,k,m)
                    presy=qpresy(i,j,k,m)
                    presz=qpresz(i,j,k,m)
                    presxy=qpresxy(i,j,k,m)
                    presxz=qpresxz(i,j,k,m)
                    presyz=qpresyz(i,j,k,m)
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
        write(wd1,'(i1)')m
        label='qpres '//wd1
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    		!
        label='qpara '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='qcross '//wd1
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='qperp '//wd1
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
    	!
        call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=hpresx(i,j,k,m)
                    presy=hpresy(i,j,k,m)
                    presz=hpresz(i,j,k,m)
                    presxy=hpresxy(i,j,k,m)
                    presxz=hpresxz(i,j,k,m)
                    presyz=hpresyz(i,j,k,m)
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
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        label='h_para '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)

        label='h_cross '//wd1
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)

        label='h_perp '//wd1
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
    	!
        call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=opresx(i,j,k,m)
                    presy=opresy(i,j,k,m)
                    presz=opresz(i,j,k,m)
                    presxy=opresxy(i,j,k,m)
                    presxz=opresxz(i,j,k,m)
                    presyz=opresyz(i,j,k,m)
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
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='o_para '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='o_cross '//wd1
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
			!
        label='o_perp '//wd1
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=sqrt(epres(i,j,k,m))
                enddo
            enddo
        enddo
        !
        label='epres '//wd1
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=abs(qpresx(i,j,k,m)+hpresx(i,j,k,m) &
                        +opresx(i,j,k,m))/3.
                    efldy(i,j,k)=abs(qpresy(i,j,k,m)+hpresy(i,j,k,m) &
                        +opresy(i,j,k,m))/3.
                    efldz(i,j,k)=abs(qpresz(i,j,k,m)+hpresz(i,j,k,m) &
                        +opresz(i,j,k,m))/3.
                    tvx(i,j,k)=sqrt(efldx(i,j,k)+efldy(i,j,k)+efldz(i,j,k) &
                        +epres(i,j,k,m))
                enddo
            enddo
        enddo
		!
        label='tpres '//wd1
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,preslim*3., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !      trace velocity streams
        !
        call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m, &
            rmassq,rmassh,rmasso)
        !
        label='pres-vel '//wd1
        call conflow(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,11,1,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !     find total magnetic field
        !
        call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
        call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
        call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
    	!
        label='box '//wd1
        call concross(tvx,curx,cury,curz,bsx,bsy,bsz, &
            nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,rearth, &
            xmin,xmax,ymin,ymax,zmin,zmax, &
            ut,label,3,11,2.0,add_dip,1,-2,start, &
            tx,ty,tz,tt,tg1,tg2,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    		!
        call contop(tvx,curx,cury,curz,bsx,bsy,bsz, &
            nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,rearth, &
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
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
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
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,12,1,2.0,2.*blim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !      calculate alfven mach number
        !
        call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m, &
            rmassq,rmassh,rmasso)
        !
        dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
        dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
        dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    aden=qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m)+1.e-8
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
                    ay=grd_ymin(m)+dy*(j-1)-ydip
                    ax=grd_xmin(m)+dx*(i-1)-xdip
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
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        label='cs_mach '//wd1 
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        label='rot_mach '//wd1 
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
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
        if(m.gt.1) then
            rho_lim=10.
        else
            rho_lim=5.0
        endif
        rfrac=0.5
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    qden=qrho(i,j,k,m)/rmassq
                    if(qden.gt.0.001)then
                        apres=(qpresx(i,j,k,m)+qpresy(i,j,k,m) &
                            +qpresz(i,j,k,m))/3.
                        bsx(i,j,k)=amin1(tempx,sqrt(apres/qden))
                    else
                        bsx(i,j,k)=0.
                    endif
                    hden=hrho(i,j,k,m)/rmassh
                    if(hden.gt.0.0005)then
                        apres=(hpresx(i,j,k,m)+hpresy(i,j,k,m) &
                            +hpresz(i,j,k,m))/3.
                        bsy(i,j,k)=amin1(temph,sqrt(apres/hden))
                    else
                        bsy(i,j,k)=0.
                    endif
                    !
                    oden=orho(i,j,k,m)/rmasso
                    if(oden.gt.0.00002)then
                        apres=(opresx(i,j,k,m)+opresy(i,j,k,m) &
                            +opresz(i,j,k,m))/3.
                        bsz(i,j,k)=amin1(tempo,sqrt(apres/oden))
                    else
                        bsz(i,j,k)=0.
                    endif
                    !
                    efldx(i,j,k)=sqrt(epres(i,j,k,m)/(qden+hden+oden))
                    !
                enddo
            enddo
        enddo
        !
        label='q temp '//wd1 
        call conhot(bsx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,tempx, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='h temp '//wd1 
        call conhot(bsy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,temph, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='o temp '//wd1 
        call conhot(bsz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,tempo, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='e temp '//wd1 
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
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
        !       call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)
                    efldx(i,j,k)=alog10(efldx(i,j,k)*rho_equiv)+6. ! per m**3
                enddo
            enddo
        enddo
        !
        if(m.le.3)then
            plot_min=5.
            plot_max=8.5
        else
            plot_min=4.-0.5*(m-3)
            plot_max=7.5-0.5*(m-3)
        endif
        !
        label='q den '//wd1 
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : ionospheric h
        !
        !      call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldy(i,j,k)=hrho(i,j,k,m)/rmassh
                    efldy(i,j,k)=alog10(efldy(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        label='h den '//wd1 
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m,&
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows: ionospheric o
        !
        !       call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldz(i,j,k)=orho(i,j,k,m)/rmasso
                    efldz(i,j,k)=alog10(efldz(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        label='o den '//wd1 
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows: total
        !
        !      call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
        !    +       opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m,
        !    +      rmassq,rmassh,rmasso)
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,m)/rmasso+ &
                        hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
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
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,m)/rmasso+ &
                    hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
                    efldy(i,j,k)=(hrho(i,j,k,m)/rmassh)/(tden+1.e-6)
                    efldz(i,j,k)=(orho(i,j,k,m)/rmasso)/(tden+1.e-6)
                    efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)/(tden+1.e-6)
                enddo
            enddo
        enddo
        !
        label='rqdens x2 '//wd1 
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,0.5, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='rhdens '//wd1 
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,1.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        label='rodens x2 '//wd1 
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,label,3,18,1,2.0,0.5, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    enddo
    !
    return
end
