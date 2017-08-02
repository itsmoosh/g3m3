subroutine push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
    efldx,efldy,efldz,nx,ny,nz,ngrd,m,delt, &
    rx,ry,rz)
    !
    !      standard runge-kutta time step
    !
    dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd),oldbz(nx,ny,nz,ngrd), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz)
    !
    !       set physical grid spacing
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    !
    ! parallelizes loop rw, aug. 17, 2004
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
                !
                !       induction equation
                !
                bx(i,j,k,m)=oldbx(i,j,k,m)-delt*( &
                (  (efldz(i,jj,k)-efldz(i,jm,k))/dyt ) &
                -( (efldy(i,j,kk)-efldy(i,j,km))/dzt )  )
                !
                by(i,j,k,m)=oldby(i,j,k,m)+delt*( &
                (  (efldz(ii,j,k)-efldz(im,j,k))/dxt ) &
                -( (efldx(i,j,kk)-efldx(i,j,km))/dzt )  )
                !
                bz(i,j,k,m)=oldbz(i,j,k,m)-delt*( &
                (  (efldy(ii,j,k)-efldy(im,j,k))/dxt) &
                -( (efldx(i,jj,k)-efldx(i,jm,k))/dyt )     )
                !
            enddo             ! k loop
        enddo             ! j loop
        !
    enddo             ! i loop
    !
    return
end
