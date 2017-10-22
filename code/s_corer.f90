!
!	This file contains two subroutines:
!	corer
!	corer_grds
!
subroutine corer(rho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      box is the coring subject grid index
    !
    integer box
    dimension rho(nx,ny,nz,n_grids)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    !    fills in the core of the coarser (big) grid in position nb
    !       using data from the finer (small) grid in position ns
    !
    !     set limits for coarse grid
    !
    mb=box+1
    ms=box
    sx=(grd_xmax(ms)-grd_xmin(ms))/(nx-1.)
    sy=(grd_ymax(ms)-grd_ymin(ms))/(ny-1.)
    sz=(grd_zmax(ms)-grd_zmin(ms))/(nz-1.)
    dx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
    dy=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
    dz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
    !
    nskip=dx/sx
    !
    ai=1.+(grd_xmin(ms)-grd_xmin(mb))/dx
    nbx1=ai+.5
    aj=1.+(grd_ymin(ms)-grd_ymin(mb))/dy
    nby1=aj+.5
    ak=1.+(grd_zmin(ms)-grd_zmin(mb))/dz
    nbz1=ak+.5
    !
    !$omp  parallel do
    do  ks=1+nskip,nz-nskip,nskip
        kb=nbz1+(ks-1)/nskip
        do  js=1+nskip,ny-nskip,nskip
            jb=nby1+(js-1)/nskip
            do  is=1+nskip,nx-nskip,nskip
                ib=nbx1+(is-1)/nskip
                !
                !       old version with no averaging
                !        rho(ib,jb,kb,mb)=rho(is,js,ks,ms)
                !
                !       smooth to make compatible with coarse grid
                !
                rho(ib,jb,kb,mb)=(6.*rho(is,js,ks,ms)+ &
                rho(is+1,js,ks,ms)+rho(is-1,js,ks,ms)+ &
                rho(is,js+1,ks,ms)+rho(is,js-1,ks,ms)+ &
                rho(is,js,ks+1,ms)+rho(is,js,ks-1,ms))/12.
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
subroutine corer_grds(rho,nx,ny,nz,n_grids,main_ngrd, &
    rho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    integer box
    dimension rho(nx,ny,nz,n_grids),rho_n(nx_n,ny_n,nz_n,ngrd_n)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    !    fills in the core of the main grid in position mb
    !       using data from the finer (small) grid in position ms
    !
    mb=main_ngrd
    ms=ngrd_n
    !
    !     set limits for coarse grid
    !
    dx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
    dy=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
    dz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
    sx=(grd_xmax_n(ms)-grd_xmin_n(ms))/(nx_n-1.)
    sy=(grd_ymax_n(ms)-grd_ymin_n(ms))/(ny_n-1.)
    sz=(grd_zmax_n(ms)-grd_zmin_n(ms))/(nz_n-1.)
    nskip=dx/sx
    !
    ai=1.+(grd_xmin_n(ms)-grd_xmin(mb))/dx
    nbx1=ai
    aj=1.+(grd_ymin_n(ms)-grd_ymin(mb))/dy
    nby1=aj
    ak=1.+(grd_zmin_n(ms)-grd_zmin(mb))/dz
    nbz1=ak
    !
    !$omp  parallel do
    do ks=1+nskip,nz_n-nskip,nskip
        kb=nbz1+(ks-1)/nskip
        do js=1+nskip,ny_n-nskip,nskip
            jb=nby1+(js-1)/nskip
            do is=1+nskip,nx_n-nskip,nskip
                ib=nbx1+(is-1)/nskip
                !
                !         rho(ib,jb,kb,mb)=rho_n(is,js,ks,ms)
                !
                rho(ib,jb,kb,mb)=(6.*rho_n(is,js,ks,ms)+ &
                rho_n(is+1,js,ks,ms)+rho_n(is-1,js,ks,ms)+ &
                rho_n(is,js+1,ks,ms)+rho_n(is,js-1,ks,ms)+ &
                rho_n(is,js,ks+1,ms)+rho_n(is,js,ks-1,ms))/12.
            enddo
        enddo
    enddo
    !
    return
end
