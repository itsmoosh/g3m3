subroutine divb_cor(bx,by,bz,px,py,pz,rho,ppres, &
    nx,ny,nz,ngrd,m,srho,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this subroutine corrects for possible non-divergence of b
    !         it integrates from the wind boundary
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd), &
    bz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd),ppres(nx,ny,nz,ngrd), &
    px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    !     set limit of x integration depending on when you
    !         hit magnetopause
    !
    r_max=4.5*srho
    r_min=0.15*srho
    rad_lim=3.*rearth
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    dxt=2.*dx
    dyt=2.*dy
    dzt=2.*dz
    !
    do k=2,nz-1
        ak=grd_zmin(m)+dz*(k-1)
        do j=2,ny-1
            aj=grd_ymin(m)+dy*(j-1)
            i=2
    	    !
            ai=grd_xmin(m)+dx*(i-1)
            ar=sqrt(ai**2+ak**2+aj**2)
            arho=rho(i,j,k,m)
            if((arho.gt.r_min).or.(ar.gt.rad_lim))then
                bx(i,j,k,m)=bx(i-1,j,k,m) &
                -(dx/dyt)*((by(i-1,j+1,k,m)+by(i,j+1,k,m))/2. &
                -(by(i-1,j-1,k,m)+by(i,j-1,k,m))/2. ) &
                -(dx/dzt)*((bz(i-1,j,k+1,m)+bz(i,j,k+1,m))/2. &
                -(bz(i-1,j,k-1,m)+bz(i,j,k-1,m))/2. )
            endif
            do i=3,nx-1
                ai=grd_xmin(m)+dx*(i-1)
                ar=sqrt(ai**2+ak**2+aj**2)
                arho=rho(i,j,k,m)
                if((arho.lt.r_min).or.(ar.lt.rad_lim)) exit
                bx(i,j,k,m)=bx(i-2,j,k,m) &
                -(dxt/dyt)*(by(i-1,j+1,k,m)-by(i-1,j-1,k,m)) &
                -(dxt/dzt)*(bz(i-1,j,k+1,m)-bz(i-1,j,k-1,m))
            enddo
        enddo
    enddo
    !
    return
end
