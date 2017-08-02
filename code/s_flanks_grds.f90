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
