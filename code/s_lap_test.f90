subroutine lap_test(qrho,qpres,qpx,qpy,qpz, &
    wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
    hrho,hpres,hpx,hpy,hpz, &
    wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz, &
    orho,opres,opx,opy,opz, &
    wrkorho,wrkopres,wrkopx,wrkopy,wrkopz, &
    nx,ny,nz,ngrd,m,chirho,chipxyz,chierg,delt, &
    rmassq,rmassh,rmasso)
    !
    dimension qrho(nx,ny,nz,ngrd),qpres(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpres(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opres(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd)
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpres(nx,ny,nz,ngrd), &
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd),wrkhrho(nx,ny,nz,ngrd), &
    wrkhpres(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd), &
    wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd), &
    wrkorho(nx,ny,nz,ngrd),wrkopres(nx,ny,nz,ngrd), &
    wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
    wrkopz(nx,ny,nz,ngrd)
    !
    write(6,*)'test lap entered'
    return
end
