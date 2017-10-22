subroutine divb_cor(bx,by,bz,px,py,pz,rho,ppres, &
    nx,ny,nz,n_grids,box,srho,r_inner, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this subroutine corrects for possible non-divergence of b
    !         it integrates from the wind boundary
    !
    real bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids), &
    bz(nx,ny,nz,n_grids),rho(nx,ny,nz,n_grids),ppres(nx,ny,nz,n_grids), &
    px(nx,ny,nz,n_grids),py(nx,ny,nz,n_grids),pz(nx,ny,nz,n_grids)
    !
    real grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
	!
	real r_max, r_min, rad_lim
	real dx, dy, dz, dxt, dyt, dzt
	integer i, j, k, box
	real ai, aj, ak, ar, arho
    !
    !     set limit of x integration depending on when you
    !         hit magnetopause
    !
    r_max=4.5*srho
    r_min=0.15*srho
    rad_lim=3.*r_inner
    !
    dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    dxt=2.*dx
    dyt=2.*dy
    dzt=2.*dz
    !
    do k=2,nz-1
        ak=grd_zmin(box)+dz*(k-1)
        do j=2,ny-1
            aj=grd_ymin(box)+dy*(j-1)
            i=2
    	    !
            ai=grd_xmin(box)+dx*(i-1)
            ar=sqrt(ai**2+ak**2+aj**2)
            arho=rho(i,j,k,box)
            if((arho.gt.r_min).or.(ar.gt.rad_lim))then
                bx(i,j,k,box)=bx(i-1,j,k,box) &
                -(dx/dyt)*((by(i-1,j+1,k,box)+by(i,j+1,k,box))/2. &
                -(by(i-1,j-1,k,box)+by(i,j-1,k,box))/2. ) &
                -(dx/dzt)*((bz(i-1,j,k+1,box)+bz(i,j,k+1,box))/2. &
                -(bz(i-1,j,k-1,box)+bz(i,j,k-1,box))/2. )
            endif
            do i=3,nx-1
                ai=grd_xmin(box)+dx*(i-1)
                ar=sqrt(ai**2+ak**2+aj**2)
                arho=rho(i,j,k,box)
                if((arho.lt.r_min).or.(ar.lt.rad_lim)) exit
                bx(i,j,k,box)=bx(i-2,j,k,box) &
                -(dxt/dyt)*(by(i-1,j+1,k,box)-by(i-1,j-1,k,box)) &
                -(dxt/dzt)*(bz(i-1,j,k+1,box)-bz(i-1,j,k-1,box))
            enddo
        enddo
    enddo
    !
    return
end
