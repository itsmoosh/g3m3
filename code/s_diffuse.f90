subroutine diffuse(qrho,qpres,qpx,qpy,qpz, &
    hrho,hpres,hpx,hpy,hpz, &
    orho,opres,opx,opy,opz,epres,bx,by,bz, &
    wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,wrkorho, &
    wrkopres,wrkopx,wrkopy,wrkopz, &
    wrkepres,wrkbx,wrkby,wrkbz,nx,ny,nz,ngrd,m, &
    difrho,difpxyz,diferg,delt,xspac,yspac,zspac)
    !
    !     apply straight artifical diffusion:
    !
    dimension qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpres(nx,ny,nz,ngrd), &
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd),wrkhrho(nx,ny,nz,ngrd), &
    wrkhpres(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd), &
    wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd), &
    wrkorho(nx,ny,nz,ngrd),wrkopres(nx,ny,nz,ngrd), &
    wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
    wrkopz(nx,ny,nz,ngrd),wrkbx(nx,ny,nz,ngrd), &
    wrkby(nx,ny,nz,ngrd),wrkbz(nx,ny,nz,ngrd), &
    wrkepres(nx,ny,nz,ngrd)
    !
    difb=diferg*delt/xspac
    dife=difrho*delt/xspac
    difp=difpxyz*delt/xspac
    difr=difrho*delt/xspac
    difo=0.25*difr
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
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
                !       species 1
                !
                qrho(i,j,k,m)=wrkqrho(i,j,k,m) +difr*( &
                +wrkqrho(ip,j,k,m)+wrkqrho(im,j,k,m) &
                +wrkqrho(i,jp,k,m)+wrkqrho(i,jm,k,m) &
                +wrkqrho(i,j,kp,m)+wrkqrho(i,j,km,m) &
                -6.*wrkqrho(i,j,k,m) )
    			!
                qpx(i,j,k,m)=wrkqpx(i,j,k,m) +difp*( &
                +wrkqpx(ip,j,k,m)+wrkqpx(im,j,k,m) &
                +wrkqpx(i,jp,k,m)+wrkqpx(i,jm,k,m) &
                +wrkqpx(i,j,kp,m)+wrkqpx(i,j,km,m) &
                -6.*wrkqpx(i,j,k,m) )
                !
                qpy(i,j,k,m)=wrkqpy(i,j,k,m) +difp*( &
                +wrkqpy(ip,j,k,m)+wrkqpy(im,j,k,m) &
                +wrkqpy(i,jp,k,m)+wrkqpy(i,jm,k,m) &
                +wrkqpy(i,j,kp,m)+wrkqpy(i,j,km,m) &
                -6.*wrkqpy(i,j,k,m) )
                !
                qpz(i,j,k,m)=wrkqpz(i,j,k,m) +difp*( &
                +wrkqpz(ip,j,k,m)+wrkqpz(im,j,k,m) &
                +wrkqpz(i,jp,k,m)+wrkqpz(i,jm,k,m) &
                +wrkqpz(i,j,kp,m)+wrkqpz(i,j,km,m) &
                -6.*wrkqpz(i,j,k,m) )
                !
                qpres(i,j,k,m)=wrkqpres(i,j,k,m) +dife*( &
                + wrkqpres(ip,j,k,m)+wrkqpres(im,j,k,m) &
                + wrkqpres(i,jp,k,m)+wrkqpres(i,jm,k,m) &
                + wrkqpres(i,j,kp,m)+wrkqpres(i,j,km,m) &
                -6.*wrkqpres(i,j,k,m) )
                !
                !       species 2
                !
                hrho(i,j,k,m)=wrkhrho(i,j,k,m) +difr*( &
                + wrkhrho(ip,j,k,m)+wrkhrho(im,j,k,m) &
                + wrkhrho(i,jp,k,m)+wrkhrho(i,jm,k,m) &
                + wrkhrho(i,j,kp,m)+wrkhrho(i,j,km,m) &
                -6.*wrkhrho(i,j,k,m) )
    			!
                hpx(i,j,k,m)=wrkhpx(i,j,k,m) +difp*( &
                +wrkhpx(ip,j,k,m)+wrkhpx(im,j,k,m) &
                +wrkhpx(i,jp,k,m)+wrkhpx(i,jm,k,m) &
                +wrkhpx(i,j,kp,m)+wrkhpx(i,j,km,m) &
                -6.*wrkhpx(i,j,k,m) )
                !
                hpy(i,j,k,m)=wrkhpy(i,j,k,m) +difp*( &
                +wrkhpy(ip,j,k,m)+wrkhpy(im,j,k,m) &
                +wrkhpy(i,jp,k,m)+wrkhpy(i,jm,k,m) &
                +wrkhpy(i,j,kp,m)+wrkhpy(i,j,km,m) &
                -6.*wrkhpy(i,j,k,m) )
                !
                hpz(i,j,k,m)=wrkhpz(i,j,k,m) +difp*( &
                +wrkhpz(ip,j,k,m)+wrkhpz(im,j,k,m) &
                +wrkhpz(i,jp,k,m)+wrkhpz(i,jm,k,m) &
                +wrkhpz(i,j,kp,m)+wrkhpz(i,j,km,m) &
                -6.*wrkhpz(i,j,k,m) )
                !
                hpres(i,j,k,m)=wrkhpres(i,j,k,m) +dife*( &
                + wrkhpres(ip,j,k,m)+wrkhpres(im,j,k,m) &
                + wrkhpres(i,jp,k,m)+wrkhpres(i,jm,k,m) &
                + wrkhpres(i,j,kp,m)+wrkhpres(i,j,km,m) &
                -6.*wrkhpres(i,j,k,m) )
                !
                !       species 3
                !
                orho(i,j,k,m)=wrkorho(i,j,k,m) +difo*( &
                + wrkorho(ip,j,k,m)+wrkorho(im,j,k,m) &
                + wrkorho(i,jp,k,m)+wrkorho(i,jm,k,m) &
                + wrkorho(i,j,kp,m)+wrkorho(i,j,km,m) &
                -6.*wrkorho(i,j,k,m) )
    			!
                opx(i,j,k,m)=wrkopx(i,j,k,m) +difo*( &
                +wrkopx(ip,j,k,m)+wrkopx(im,j,k,m) &
                +wrkopx(i,jp,k,m)+wrkopx(i,jm,k,m) &
                +wrkopx(i,j,kp,m)+wrkopx(i,j,km,m) &
                -6.*wrkopx(i,j,k,m) )
                !
                opy(i,j,k,m)=wrkopy(i,j,k,m) +difo*( &
                +wrkopy(ip,j,k,m)+wrkopy(im,j,k,m) &
                +wrkopy(i,jp,k,m)+wrkopy(i,jm,k,m) &
                +wrkopy(i,j,kp,m)+wrkopy(i,j,km,m) &
                -6.*wrkopy(i,j,k,m) )
                !
                opz(i,j,k,m)=wrkopz(i,j,k,m) +difo*( &
                +wrkopz(ip,j,k,m)+wrkopz(im,j,k,m) &
                +wrkopz(i,jp,k,m)+wrkopz(i,jm,k,m) &
                +wrkopz(i,j,kp,m)+wrkopz(i,j,km,m) &
                -6.*wrkopz(i,j,k,m) )
                !
                opres(i,j,k,m)=wrkopres(i,j,k,m) +dife*( &
                + wrkopres(ip,j,k,m)+wrkopres(im,j,k,m) &
                + wrkopres(i,jp,k,m)+wrkopres(i,jm,k,m) &
                + wrkopres(i,j,kp,m)+wrkopres(i,j,km,m) &
                -6.*wrkopres(i,j,k,m) )
                !
                epres(i,j,k,m)=wrkepres(i,j,k,m) +dife*( &
                + wrkepres(ip,j,k,m)+wrkepres(im,j,k,m) &
                + wrkepres(i,jp,k,m)+wrkepres(i,jm,k,m) &
                + wrkepres(i,j,kp,m)+wrkepres(i,j,km,m) &
                -6.*wrkepres(i,j,k,m) )
                !
                bx(i,j,k,m)=wrkbx(i,j,k,m)+difb*( &
                + wrkbx(ip,j,k,m)+wrkbx(im,j,k,m) &
                + wrkbx(i,jp,k,m)+wrkbx(i,jm,k,m) &
                + wrkbx(i,j,kp,m)+wrkbx(i,j,km,m) &
                -6.*wrkbx(i,j,k,m) )
                !
                by(i,j,k,m)=wrkby(i,j,k,m)+difb*( &
                + wrkby(ip,j,k,m)+wrkby(im,j,k,m) &
                + wrkby(i,jp,k,m)+wrkby(i,jm,k,m) &
                + wrkby(i,j,kp,m)+wrkby(i,j,km,m) &
                -6.*wrkby(i,j,k,m) )
                !
                bz(i,j,k,m)=wrkbz(i,j,k,m)+difb*( &
                + wrkbz(ip,j,k,m)+wrkbz(im,j,k,m) &
                + wrkbz(i,jp,k,m)+wrkbz(i,jm,k,m) &
                + wrkbz(i,j,kp,m)+wrkbz(i,j,km,m) &
                -6.*wrkbz(i,j,k,m) )
            enddo
        enddo
    enddo
    !
    !     reset time index so that desired quanties lie at nt1 rather than nt2
    !
    !     ntnew=nt2
    !     nt2=nt1
    !     nt1=ntnew
    !
    return
end
