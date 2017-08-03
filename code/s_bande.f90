subroutine bande(efldx,efldy,efldz,bsx,bsy,bsz, &
    curx,cury,curz,evx,evy,evz,btot, &
    epres,qrho,hrho,orho,rst,resist,reynolds, &
    nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
    ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
    rx,ry,rz)
    !
    !      calculates the surface magnetic field bs
    !      and the body electric field eb
    !
    dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
    efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz), &
    qrho(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd),btot(nx,ny,nz)
    dimension rst(nx,ny,nz,mbndry)
    integer ijmid(mbndry,3,mmid),ijzero(mbndry,3,mzero)
    integer nummid(mbndry),numzero(mbndry)
    !
    !      ohm's law: ve = vi-j/ne
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                !
                avx=evx(i,j,k)
                avy=evy(i,j,k)
                avz=evz(i,j,k)
                !
                abx=bsx(i,j,k)
                aby=bsy(i,j,k)
                abz=bsz(i,j,k)
                !
                efldx(i,j,k)=-(avy*abz-avz*aby)
                efldy(i,j,k)=-(avz*abx-avx*abz)
                efldz(i,j,k)=-(avx*aby-avy*abx)
                !
            enddo
        enddo
    enddo
    !
    !     add in grad p term
    !
    dxt=2.*rx
    dyt=2.*ry
    dzt=2.*rz
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=2,nz-1
        kp=k+1
        km=k-1
        do j=2,ny-1
            jp=j+1
            jm=j-1
            do i=2,nx-1
                ip=i+1
                im=i-1
                !
                arho=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
                orho(i,j,k,m)/rmasso)*reynolds
                !
                efldx(i,j,k)=efldx(i,j,k)- &
                ((epres(ip,j,k,m)-epres(im,j,k,m))/dxt)/arho
                efldy(i,j,k)=efldy(i,j,k)- &
                ((epres(i,jp,k,m)-epres(i,jm,k,m))/dyt)/arho
                efldz(i,j,k)=efldz(i,j,k)- &
                ((epres(i,j,kp,m)-epres(i,j,km,m))/dzt)/arho
                !
            enddo
        enddo
    enddo
    !
    !      add in ionospheric resistance
    !
    if((m.le.mbndry).and.(resist.lt.5000.))then
        ! parallelizes loop. RW, oct. 23, 2002
        !$omp  parallel do
        do n=1,nummid(m)
            i=ijmid(m,1,n)
            j=ijmid(m,2,n)
            k=ijmid(m,3,n)
            !
            eden=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
            orho(i,j,k,m)/rmasso)
            !
            efldx(i,j,k)=efldx(i,j,k)+rst(i,j,k,m)*curx(i,j,k)/eden
            efldy(i,j,k)=efldy(i,j,k)+rst(i,j,k,m)*cury(i,j,k)/eden
            efldz(i,j,k)=efldz(i,j,k)+rst(i,j,k,m)*curz(i,j,k)/eden
        enddo
        do n=1,numzero(m)
            i=ijzero(m,1,n)
            j=ijzero(m,2,n)
            k=ijzero(m,3,n)
            !
            eden=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
            orho(i,j,k,m)/rmasso)
            !
            efldx(i,j,k)=efldx(i,j,k)+rst(i,j,k,m)*curx(i,j,k)/eden
            efldy(i,j,k)=efldy(i,j,k)+rst(i,j,k,m)*cury(i,j,k)/eden
            efldz(i,j,k)=efldz(i,j,k)+rst(i,j,k,m)*curz(i,j,k)/eden
        enddo
    endif
    !
    !      boundary condtions for terrestrial surface currents
    !
    !     if(m.le.mbndry)then
    !      do n=1,nummid(m)
    !        i=ijmid(m,1,n)
    !        j=ijmid(m,2,n)
    !        k=ijmid(m,3,n)
    !        efldx(i,j,k)=0.
    !        efldy(i,j,k)=0.
    !        efldz(i,j,k)=0.
    !       enddo
    !      endif
    !
    !      set flank bounday conditions
    !
    nx1=nx-1
    ny1=ny-1
    nz1=nz-1
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            efldx(1,j,k)=efldx(2,j,k)
            efldx(nx,j,k)=efldx(nx1,j,k)
            efldy(1,j,k)=efldy(2,j,k)
            efldy(nx,j,k)=efldy(nx1,j,k)
            efldz(1,j,k)=efldz(2,j,k)
            efldz(nx,j,k)=efldz(nx1,j,k)
        enddo
    enddo
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do j=1,ny
        do i=1,nx
            efldx(i,j,1)=efldx(i,j,2)
            efldx(i,j,nz)=efldx(i,j,nz1)
            efldy(i,j,1)=efldy(i,j,2)
            efldy(i,j,nz)=efldy(i,j,nz1)
            efldz(i,j,1)=efldz(i,j,2)
            efldz(i,j,nz)=efldz(i,j,nz1)
        enddo
    enddo
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do k=1,nz
        do i=1,nx
            efldx(i,1,k)=efldx(i,2,k)
            efldx(i,ny,k)=efldx(i,ny1,k)
            efldy(i,1,k)=efldy(i,2,k)
            efldy(i,ny,k)=efldy(i,ny1,k)
            efldz(i,1,k)=efldz(i,2,k)
            efldz(i,ny,k)=efldz(i,ny1,k)
        enddo
    enddo
    !
    return
end
