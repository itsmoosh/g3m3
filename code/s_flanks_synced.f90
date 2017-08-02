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
