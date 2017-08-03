subroutine refinement(rho,nx,ny,nz,ngrd,main_grid, &
    rho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    dimension rho(nx,ny,nz,ngrd),rho_n(nx_n,ny_n,nz_n,ngrd_n)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    !    fills in the core of the coarser (big) grid in position nb
    !       using data from the finer (small) grid in position ns
    !
    !     set limits for coarse grid
    !
    rho_n=0.
    m=main_grid
    m_n=ngrd_n
    rho_max=0.
    rho_n_max=0.
    !
    rx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    ry=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    rz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    sx=(grd_xmax_n(m_n)-grd_xmin_n(m_n))/(nx_n-1.)
    sy=(grd_ymax_n(m_n)-grd_ymin_n(m_n))/(ny_n-1.)
    sz=(grd_zmax_n(m_n)-grd_zmin_n(m_n))/(nz_n-1.)
    !
    !    interpolate between main and sub grids
    !
    !$omp  parallel do
    do k_n=1,nz_n
        az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
        !
        !     find position on main grid
        !
        ak=1.+(az_n-grd_zmin(m))/rz
        k=ak
        kk=k+1
        dz=ak-k
        ddz=1.-dz
        !
        do j_n=1,ny_n
            ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
            aj=1.+(ay_n-grd_ymin(m))/ry
            j=aj
            jj=j+1
            dy=aj-j
            ddy=1.-dy
            !
            do i_n=1,nx_n
                ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                ai=1.+(ax_n-grd_xmin(m))/rx
                i=ai
                ii=i+1
                dx=ai-i
                ddx=1.-dx
                !
                rho_n(i_n,j_n,k_n,m_n)= &
                rho(i,j,k,m)*ddx*ddy*ddz+ &
                rho(i,j,kk,m)*ddx*ddy*dz+ &
                rho(i,jj,k,m)*ddx*dy*ddz+ &
                rho(i,jj,kk,m)*ddx*dy*dz+ &
                rho(ii,j,k,m)*dx*ddy*ddz+ &
                rho(ii,j,kk,m)*dx*ddy*dz+ &
                rho(ii,jj,k,m)*dx*dy*ddz+ &
                rho(ii,jj,kk,m)*dx*dy*dz
                !
                !       rho_n(i_n,j_n,k_n,m_n)= rho(i,j,k,m)
                !
                !       write(7,*)i,j,k,m
                !          rho_max=amax1(rho_max,rho(i,j,k,m))
                !          rho_n_max=amax1(rho_n_max,rho_n(i_n,j_n,k_n,m_n))
            enddo
        enddo
    enddo
    !      write(6,*)'refine',rho_max,rho_n_max
    !
    !      interpolate onto subgrids
    !
    do m_n=ngrd_n-1,1,-1
        m=m_n+1
        rx=(grd_xmax_n(m)-grd_xmin_n(m))/(nx_n-1.)
        ry=(grd_ymax_n(m)-grd_ymin_n(m))/(ny_n-1.)
        rz=(grd_zmax_n(m)-grd_zmin_n(m))/(nz_n-1.)
        sx=(grd_xmax_n(m_n)-grd_xmin_n(m_n))/(nx_n-1.)
        sy=(grd_ymax_n(m_n)-grd_ymin_n(m_n))/(ny_n-1.)
        sz=(grd_zmax_n(m_n)-grd_zmin_n(m_n))/(nz_n-1.)
        !
        !    interpolate between main and sub grids
        !
        do k_n=1,nz_n
            az_n=grd_zmin_n(m_n)+(k_n-1.)*sz
            !
            !     find position on bigger grid
            !
            ak=1.+(az_n-grd_zmin_n(m))/rz
            k=ak
            kk=k+1
            dz=ak-k
            ddz=1.-dz
            !
            do j_n=1,ny_n
                ay_n=grd_ymin_n(m_n)+(j_n-1.)*sy
                aj=1.+(ay_n-grd_ymin_n(m))/ry
                j=aj
                jj=j+1
                dy=aj-j
                ddy=1.-dy
                !
                do i_n=1,nx_n
                    ax_n=grd_xmin_n(m_n)+(i_n-1.)*sx
                    ai=1.+(ax_n-grd_xmin_n(m))/rx
                    i=ai
                    ii=i+1
                    dx=ai-i
                    ddx=1.-dx
                    !
                    rho_n(i_n,j_n,k_n,m_n)= &
                    rho_n(i,j,k,m)*ddx*ddy*ddz+ &
                    rho_n(i,j,kk,m)*ddx*ddy*dz+ &
                    rho_n(i,jj,k,m)*ddx*dy*ddz+ &
                    rho_n(i,jj,kk,m)*ddx*dy*dz+ &
                    rho_n(ii,j,k,m)*dx*ddy*ddz+ &
                    rho_n(ii,j,kk,m)*dx*ddy*dz+ &
                    rho_n(ii,jj,k,m)*dx*dy*ddz+ &
                    rho_n(ii,jj,kk,m)*dx*dy*dz
                    !
                enddo
            enddo
        enddo
    enddo
    !
    return
end
