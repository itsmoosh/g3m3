!
!	This file contains two subroutines:
!	set_bndry_moon
!	set_bndry_moon_ram
!
subroutine set_bndry_moon(rmassq,rmassh,rmasso, &
    box,nx,ny,nz,n_grids, &
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
    integer box
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !      write(grdpt_f,*)'bndry_moon values', &
    !           qden_moon,hden_moon,oden_moon,tempi,gamma, &
    !           ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_m
    !      write(grdpt_f,*)mbndry,msrf,mmid,mzero
    !
    !      set scale lengths
    !
    dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    !      reset indices
    !
    numsrf(box)=0
    nummid(box)=0
    numzero(box)=0
    ntot=0
    !
    !       set conditions at moon's new position
    !
    do k=1,nz
        az=grd_zmin(box)+dz*(k-1)
        zm=az-zmoon
        do j=1,ny
            ay=grd_ymin(box)+dy*(j-1)
            ym=ay-ymoon
            do i=1,nx
                ax=grd_xmin(box)+dx*(i-1)
                xm=ax-xmoon
                ar_moon=sqrt(xm**2+ym**2+zm**2)
                !
                ra_moon=((ar_moon+0.3*rmoon)/(1.3*rmoon))**(-alpha_m)
                ra_moon=amin1(1.,ra_moon)
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
    			!
                if(ar_moon.le.rmoon+0.6*dx)then
                    ntot=ntot+1
                endif
                !
                if(ar_moon.le.rmoon-1.5*dx)then
                    numzero(box)=numzero(box)+1
                    if(numzero(box).gt.mzero)then
                        write(6,*)'numzero too large',box,numzero(box),mzero
                        stop
                    endif
                    ijzero(box,1,numzero(box))=i
                    ijzero(box,2,numzero(box))=j
                    ijzero(box,3,numzero(box))=k
                    !
                    parm_zero(box,1,numzero(box))=qden*rmassq
                    parm_zero(box,2,numzero(box))=hden*rmassh
                    parm_zero(box,3,numzero(box))=oden*rmasso
                    parm_zero(box,4,numzero(box))=qpress
                    parm_zero(box,5,numzero(box))=hpress
                    parm_zero(box,6,numzero(box))=opress
                    parm_zero(box,7,numzero(box))=epress
                    !
                else  if(ar_moon.le.rmoon-0.5*dx) then
                    nummid(box)=nummid(box)+1
                    if(nummid(box).gt.mmid)then
                        write(*,*) 'nummid too large. box, nummid(box), mmid:'
						write(*,*) box, nummid(box), mmid
                        stop
                    endif
                    ijmid(box,1,nummid(box))=i
                    ijmid(box,2,nummid(box))=j
                    ijmid(box,3,nummid(box))=k
                    !
                    parm_mid(box,1,nummid(box))=qden*rmassq
                    parm_mid(box,2,nummid(box))=hden*rmassh
                    parm_mid(box,3,nummid(box))=oden*rmasso
                    parm_mid(box,4,nummid(box))=qpress
                    parm_mid(box,5,nummid(box))=hpress
                    parm_mid(box,6,nummid(box))=opress
                    parm_mid(box,7,nummid(box))=epress
                    !
                else if(ar_moon.le.rmoon+0.6*dx) then
                    numsrf(box)=numsrf(box)+1
                    if(numsrf(box).gt.msrf)then
                        write(6,*)'numsrf too large',box,numsrf(box),msrf
                        stop
                        !                endif
                        ijsrf(box,1,numsrf(box))=i
                        ijsrf(box,2,numsrf(box))=j
                        ijsrf(box,3,numsrf(box))=k
                        !
                        parm_srf(box,1,numsrf(box))=qden*rmassq
                        parm_srf(box,2,numsrf(box))=hden*rmassh
                        parm_srf(box,3,numsrf(box))=oden*rmasso
                        parm_srf(box,4,numsrf(box))=qpress
                        parm_srf(box,5,numsrf(box))=hpress
                        parm_srf(box,6,numsrf(box))=opress
                        parm_srf(box,7,numsrf(box))=epress
                        !
                    endif
                endif
                !
            enddo
        enddo
    enddo
    !
    !         write(grdpt_f,*) 'Total pts: ', ntot
    !         write(grdpt_f,*) 'moon bndry_m pts: ', box, numsrf(box), nummid(box), numzero(box)
    !
    return
end
!
!
!	****************************************
!
!
subroutine set_bndry_moon_ram(rmassq,rmassh,rmasso, &
    box,nx,ny,nz,n_grids, &
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
    integer box
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !      write(grdpt_f,*)'bndry_moon values', &
    !           qden_moon,hden_moon,oden_moon,tempi,gamma, &
    !           ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_m
    !      write(grdpt_f,*)mbndry,msrf,mmid,mzero
    !
    !      set scale lengths
    !
    dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    !      reset indices
    !
    numsrf(box)=0
    nummid(box)=0
    numzero(box)=0
    ntot=0
    !
    !      speed for ram induced ionosphere
    !
    vtot=sqrt(vx_moon**2+vy_moon**2)+1.e-11
    !
    !       set conditions at moon's new position
    !
    do k=1,nz
        az=grd_zmin(box)+dz*(k-1)
        zm=az-zmoon
        do j=1,ny
            ay=grd_ymin(box)+dy*(j-1)
            ym=ay-ymoon
            do i=1,nx
                ax=grd_xmin(box)+dx*(i-1)
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
    			!
                if(ar_moon.le.rmoon+0.6*dx)then
                    ntot=ntot+1
                endif
                !
                if(ar_moon.le.rmoon-1.5*dx)then
                    numzero(box)=numzero(box)+1
                    if(numzero(box).gt.mzero)then
                        write(6,*)'numzero too large',box,numzero(box),mzero
                        stop
                    endif
                    ijzero(box,1,numzero(box))=i
                    ijzero(box,2,numzero(box))=j
                    ijzero(box,3,numzero(box))=k
                    !
                    parm_zero(box,1,numzero(box))=qden*rmassq
                    parm_zero(box,2,numzero(box))=hden*rmassh
                    parm_zero(box,3,numzero(box))=oden*rmasso
                    parm_zero(box,4,numzero(box))=qpress
                    parm_zero(box,5,numzero(box))=hpress
                    parm_zero(box,6,numzero(box))=opress
                    parm_zero(box,7,numzero(box))=epress
                    !
                else if(ar_moon.le.rmoon-0.5*dx) then
                    nummid(box)=nummid(box)+1
                    if(nummid(box).gt.mmid)then
                        write(6,*)'nummid too large',box,nummid(box),mmid
                        stop
                    endif
                    ijmid(box,1,nummid(box))=i
                    ijmid(box,2,nummid(box))=j
                    ijmid(box,3,nummid(box))=k
                    !
                    parm_mid(box,1,nummid(box))=qden*rmassq
                    parm_mid(box,2,nummid(box))=hden*rmassh
                    parm_mid(box,3,nummid(box))=oden*rmasso
                    parm_mid(box,4,nummid(box))=qpress
                    parm_mid(box,5,nummid(box))=hpress
                    parm_mid(box,6,nummid(box))=opress
                    parm_mid(box,7,nummid(box))=epress
                    !
                else if(ar_moon.le.rmoon+0.6*dx) then
                    numsrf(box)=numsrf(box)+1
                    if(numsrf(box).gt.msrf)then
                        write(6,*)'numsrf too large',box,numsrf(box),msrf
                        stop
                        !                endif
                        ijsrf(box,1,numsrf(box))=i
                        ijsrf(box,2,numsrf(box))=j
                        ijsrf(box,3,numsrf(box))=k
                        !
                        parm_srf(box,1,numsrf(box))=qden*rmassq
                        parm_srf(box,2,numsrf(box))=hden*rmassh
                        parm_srf(box,3,numsrf(box))=oden*rmasso
                        parm_srf(box,4,numsrf(box))=qpress
                        parm_srf(box,5,numsrf(box))=hpress
                        parm_srf(box,6,numsrf(box))=opress
                        parm_srf(box,7,numsrf(box))=epress
                        !
                    endif
                endif
                !
            enddo
        enddo
    enddo
    !
    !         write(grdpt_f,*) 'Total pts: ', ntot
    !         write(grdpt_f,*) 'moon bndry_m pts: ', box, numsrf(box), nummid(box), numzero(box)
    !
    return
end
