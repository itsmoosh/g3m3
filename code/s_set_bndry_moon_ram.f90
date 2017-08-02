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
    			!
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
