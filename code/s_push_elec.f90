subroutine push_elec(epres,oldepres,wrkepres,evx,evy,evz, &
    gamma,gamma1,nx,ny,nz,ngrd,m,delt,rx,ry,rz)
    !
    !      evolves the electron pressure equation
    !
    !
    dimension epres(nx,ny,nz,ngrd),oldepres(nx,ny,nz,ngrd), &
    wrkepres(nx,ny,nz,ngrd)
    !
    dimension evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
    !
    !
    !       set physical grid spacing
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    !
    !$omp  parallel do
    do k=2,nz-1
        km=k-1
        kk=k+1
    
        do j=2,ny-1
            jm=j-1
            jj=j+1
    
            do i=2,nx-1
                ii=i+1
                im=i-1
                egradp_x=(wrkepres(ii,j,k,m)-wrkepres(im,j,k,m))/dxt
                egradp_y=(wrkepres(i,jj,k,m)-wrkepres(i,jm,k,m))/dyt
                egradp_z=(wrkepres(i,j,kk,m)-wrkepres(i,j,km,m))/dzt
                !
                !
                !       pressure equations:
                !
                epres(i,j,k,m)=oldepres(i,j,k,m)-delt*gamma* &
                ( ( (wrkepres(ii,j,k,m)*evx(ii,j,k) &
                -wrkepres(im,j,k,m)*evx(im,j,k))/dxt ) + &
                ( (wrkepres(i,jj,k,m)*evy(i,jj,k) &
                -wrkepres(i,jm,k,m)*evy(i,jm,k))/dyt ) + &
                ( (wrkepres(i,j,kk,m)*evz(i,j,kk) &
                -wrkepres(i,j,km,m)*evz(i,j,km))/dzt ) )
                epres(i,j,k,m)=epres(i,j,k,m) &
                + delt*gamma1*( &
                evx(i,j,k)*egradp_x+evy(i,j,k)*egradp_y &
                +evz(i,j,k)*egradp_z )
                !
            enddo
        enddo
    enddo
    !
    return
end
