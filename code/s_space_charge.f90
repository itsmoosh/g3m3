subroutine space_charge(charge,efldx,efldy,efldz, &
    bx0,by0,bz0,px,py,pz,rho,nx,ny,nz,n_grids,box,r_inner, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      calculates the space-charge fields in the plasma
    !
    integer box
    dimension bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids), &
    bz0(nx,ny,nz,n_grids),charge(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    px(nx,ny,nz,n_grids),py(nx,ny,nz,n_grids),pz(nx,ny,nz,n_grids), &
    rho(nx,ny,nz,n_grids)
    !
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    rx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    ry=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    rz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    !
    d_min=0.001
    !
    !     take the divergence of e and correct
    !          for curl b_0 errors
    !
    do k=2,nz-1
        km=k-1
        kp=k+1
        do j=2,ny-1
            jm=j-1
            jp=j+1
            do i=2,nx-1
                im=i-1
                ip=i+1
                !
                !     div. e
                !
                charge(i,j,k)= &
                (efldx(ip,j,k)-efldx(im,j,k))/dxt &
                +(efldy(i,jp,k)-efldy(i,jm,k))/dyt &
                +(efldz(i,j,kp)-efldz(i,j,km))/dzt
                !
                !     correct for curl b_0 errors
                !
                arho=rho(i,j,k,box)
                arho=amax1(arho,d_min)
                !
                apx=px(i,j,k,box)/arho
                apy=py(i,j,k,box)/arho
                apz=pz(i,j,k,box)/arho
                !
                curlbx=( (bz0(i,jp,k,box)-bz0(i,jm,k,box))/dyt &
                -(by0(i,j,kp,box)-by0(i,j,km,box))/dzt )/arho
                curlby=( (bx0(i,j,kp,box)-bx0(i,j,km,box))/dzt &
                -(bz0(ip,j,k,box)-bz0(im,j,k,box))/dxt )/arho
                curlbz=( (by0(ip,j,k,box)-by0(im,j,k,box))/dxt &
                -(bx0(i,jp,k,box)-bx0(i,jm,k,box))/dyt )/arho
                !
                charge(i,j,k)=charge(i,j,k) &
                -(apx*curlbx+apy*curlby+apz*curlbz)
                !
            enddo
        enddo
    enddo
    !
    !      set flank boundary conditions
    !
    nx1=nx-1
    ny1=ny-1
    nz1=nz-1
    !
    do k=1,nz
        do j=1,ny
            charge(1,j,k)=charge(2,j,k)
            charge(nx,j,k)=charge(nx1,j,k)
        enddo
    enddo
    !
    do j=1,ny
        do i=1,nx
            charge(i,j,1)=charge(i,j,2)
            charge(i,j,nz)=charge(i,j,nz1)
        enddo
    enddo
    !
    do k=1,nz
        do i=1,nx
            charge(i,1,k)=charge(i,2,k)
            charge(i,ny,k)=charge(i,ny1,k)
        enddo
    enddo
    !
    return
end
