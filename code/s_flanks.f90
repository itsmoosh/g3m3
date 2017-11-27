!
!	This file contains three subroutines:
!	flanks
!	flanks_grds
!	flanks_synced
!
subroutine flanks(wrkrho,rho,oldrho,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
	integer, parameter :: dp = kind(1.d0)
    integer box
	real(dp) t_old(n_grids), t_new(n_grids), t
    dimension rho(nx,ny,nz,n_grids),oldrho(nx,ny,nz,n_grids), &
    wrkrho(nx,ny,nz,n_grids),work(nx,ny,nz)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    !     box is the suject array index
    !
    !      sets the outer boundaries of fine grid at position ns
    !      using data from coarse grid at position nb
    !
    mb=box+1
    ms=box
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
    !     do bottom panels
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
    !     do front panels
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
    !     do back panels
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
            enddo
        enddo
    enddo
    !
    !     do front-left panels
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
!
!	****************************************
!
!
subroutine flanks_grds( &
    rho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd, &
    rho,oldrho,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    integer box
    dimension rho(nx,ny,nz,n_grids),oldrho(nx,ny,nz,n_grids), &
    rho_n(nx_n,ny_n,nz_n,ngrd_n), &
    work(nx,ny,nz)
    dimension t_new(n_grids),t_old(n_grids)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    !     box is the suject array index
    !
    !      sets the outer boundaries of fine grid at position ns
    !      using data from coarse grid at position nb
    !
    mb=main_ngrd
    !
    t1=(t_new(mb)-t)/(t_new(mb)-t_old(mb))
    t2=1.-t1
    !
    !     write(grdpt_f,*)'flanks_grds',mb,t1,t2
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
subroutine flanks_synced(rho,nx,ny,nz,n_grids,box, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    integer box
    dimension rho(nx,ny,nz,n_grids)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    !      sets the outer boundaries of fine grid at position ns
    !      using data from coarse grid at position nb
    !
    mb=box+1
    nb=mb
    !
    ms=box
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
    !     do front and back panels
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
