subroutine refinement(rho,nx,ny,nz,n_grids,main_grid, &
    rho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    integer box, box_n
    dimension rho(nx,ny,nz,n_grids),rho_n(nx_n,ny_n,nz_n,ngrd_n)
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
    !    fills in the core of the coarser (big) grid in position nb
    !       using data from the finer (small) grid in position ns
    !
    !     set limits for coarse grid
    !
    rho_n=0.
    box=main_grid
    box_n=ngrd_n
    rho_max=0.
    rho_n_max=0.
    !
    rx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    ry=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    rz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    sx=(grd_xmax_n(box_n)-grd_xmin_n(box_n))/(nx_n-1.)
    sy=(grd_ymax_n(box_n)-grd_ymin_n(box_n))/(ny_n-1.)
    sz=(grd_zmax_n(box_n)-grd_zmin_n(box_n))/(nz_n-1.)
    !
    !    interpolate between main and sub grids
    !
    !$omp  parallel do
    do k_n=1,nz_n
        az_n=grd_zmin_n(box_n)+(k_n-1.)*sz
        !
        !     find position on main grid
        !
        ak=1.+(az_n-grd_zmin(box))/rz
        k=ak
        kk=k+1
        dz=ak-k
        ddz=1.-dz
        !
        do j_n=1,ny_n
            ay_n=grd_ymin_n(box_n)+(j_n-1.)*sy
            aj=1.+(ay_n-grd_ymin(box))/ry
            j=aj
            jj=j+1
            dy=aj-j
            ddy=1.-dy
            !
            do i_n=1,nx_n
                ax_n=grd_xmin_n(box_n)+(i_n-1.)*sx
                ai=1.+(ax_n-grd_xmin(box))/rx
                i=ai
                ii=i+1
                dx=ai-i
                ddx=1.-dx
                !
                rho_n(i_n,j_n,k_n,box_n)= &
                rho(i,j,k,box)*ddx*ddy*ddz+ &
                rho(i,j,kk,box)*ddx*ddy*dz+ &
                rho(i,jj,k,box)*ddx*dy*ddz+ &
                rho(i,jj,kk,box)*ddx*dy*dz+ &
                rho(ii,j,k,box)*dx*ddy*ddz+ &
                rho(ii,j,kk,box)*dx*ddy*dz+ &
                rho(ii,jj,k,box)*dx*dy*ddz+ &
                rho(ii,jj,kk,box)*dx*dy*dz
                !
                !       rho_n(i_n,j_n,k_n,box_n)= rho(i,j,k,box)
                !
                !       write(grdpt_f,*)i,j,k,box
                !          rho_max=amax1(rho_max,rho(i,j,k,box))
                !          rho_n_max=amax1(rho_n_max,rho_n(i_n,j_n,k_n,box_n))
            enddo
        enddo
    enddo
    !      write(fluxs_f,*)'refine',rho_max,rho_n_max
    !
    !      interpolate onto subgrids
    !
    do box_n=ngrd_n-1,1,-1
        box=box_n+1
        rx=(grd_xmax_n(box)-grd_xmin_n(box))/(nx_n-1.)
        ry=(grd_ymax_n(box)-grd_ymin_n(box))/(ny_n-1.)
        rz=(grd_zmax_n(box)-grd_zmin_n(box))/(nz_n-1.)
        sx=(grd_xmax_n(box_n)-grd_xmin_n(box_n))/(nx_n-1.)
        sy=(grd_ymax_n(box_n)-grd_ymin_n(box_n))/(ny_n-1.)
        sz=(grd_zmax_n(box_n)-grd_zmin_n(box_n))/(nz_n-1.)
        !
        !    interpolate between main and sub grids
        !
        do k_n=1,nz_n
            az_n=grd_zmin_n(box_n)+(k_n-1.)*sz
            !
            !     find position on bigger grid
            !
            ak=1.+(az_n-grd_zmin_n(box))/rz
            k=ak
            kk=k+1
            dz=ak-k
            ddz=1.-dz
            !
            do j_n=1,ny_n
                ay_n=grd_ymin_n(box_n)+(j_n-1.)*sy
                aj=1.+(ay_n-grd_ymin_n(box))/ry
                j=aj
                jj=j+1
                dy=aj-j
                ddy=1.-dy
                !
                do i_n=1,nx_n
                    ax_n=grd_xmin_n(box_n)+(i_n-1.)*sx
                    ai=1.+(ax_n-grd_xmin_n(box))/rx
                    i=ai
                    ii=i+1
                    dx=ai-i
                    ddx=1.-dx
                    !
                    rho_n(i_n,j_n,k_n,box_n)= &
                    rho_n(i,j,k,box)*ddx*ddy*ddz+ &
                    rho_n(i,j,kk,box)*ddx*ddy*dz+ &
                    rho_n(i,jj,k,box)*ddx*dy*ddz+ &
                    rho_n(i,jj,kk,box)*ddx*dy*dz+ &
                    rho_n(ii,j,k,box)*dx*ddy*ddz+ &
                    rho_n(ii,j,kk,box)*dx*ddy*dz+ &
                    rho_n(ii,jj,k,box)*dx*dy*ddz+ &
                    rho_n(ii,jj,kk,box)*dx*dy*dz
                    !
                enddo
            enddo
        enddo
    enddo
    !
    return
end
