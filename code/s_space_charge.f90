subroutine space_charge(charge,efldx,efldy,efldz, &
    bx0,by0,bz0,px,py,pz,rho,nx,ny,nz,ngrd,m,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      calculates the space-charge fields in the plasma
    !
    dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd), &
    bz0(nx,ny,nz,ngrd),charge(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd), &
    rho(nx,ny,nz,ngrd)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    rx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    ry=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    rz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
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
                arho=rho(i,j,k,m)
                arho=amax1(arho,d_min)
                !
                apx=px(i,j,k,m)/arho
                apy=py(i,j,k,m)/arho
                apz=pz(i,j,k,m)/arho
                !
                curlbx=( (bz0(i,jp,k,m)-bz0(i,jm,k,m))/dyt &
                -(by0(i,j,kp,m)-by0(i,j,km,m))/dzt )/arho
                curlby=( (bx0(i,j,kp,m)-bx0(i,j,km,m))/dzt &
                -(bz0(ip,j,k,m)-bz0(im,j,k,m))/dxt )/arho
                curlbz=( (by0(ip,j,k,m)-by0(im,j,k,m))/dxt &
                -(bx0(i,jp,k,m)-bx0(i,jm,k,m))/dyt )/arho
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
