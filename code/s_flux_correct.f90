subroutine flux_correct(qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
    wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
    wrkqpresxy,wrkqpresxz,wrkqpresyz, &
    wrkqpx,wrkqpy,wrkqpz, &
    oldqrho,oldqpresx,oldqpresy,oldqpresz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    oldqpx,oldqpy,oldqpz, &
    !
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
    wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
    wrkhpresxy,wrkhpresxz,wrkhpresyz, &
    wrkhpx,wrkhpy,wrkhpz, &
    oldhrho,oldhpresx,oldhpresy,oldhpresz, &
    oldhpresxy,oldhpresxz,oldhpresyz, &
    oldhpx,oldhpy,oldhpz, &
    !
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,opx,opy,opz, &
    wrkorho,wrkopresx,wrkopresy,wrkopresz, &
    wrkopresxy,wrkopresxz,wrkopresyz, &
    wrkopx,wrkopy,wrkopz, &
    oldorho,oldopresx,oldopresy,oldopresz, &
    oldopresxy,oldopresxz,oldopresyz, &
    oldopx,oldopy,oldopz, &
    !
    epres,wrkepres,oldepres, &
    bx,by,bz,wrkbx,wrkby,wrkbz, &
    oldbx,oldby,oldbz,vvx,vvy,vvz, &
    nx,ny,nz,ngrd,m,difrho,diferg,xspac, &
    isotropic)
    !
    !     apply flux correction smoothing technique
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd), &
    bz(nx,ny,nz,ngrd),epres(nx,ny,nz,ngrd), &
    vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz), &
    xspac(ngrd)
    !
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpx(nx,ny,nz,ngrd), &
    wrkqpy(nx,ny,nz,ngrd),wrkqpz(nx,ny,nz,ngrd), &
    wrkqpresx(nx,ny,nz,ngrd),wrkqpresy(nx,ny,nz,ngrd), &
    wrkqpresz(nx,ny,nz,ngrd),wrkqpresxy(nx,ny,nz,ngrd), &
    wrkqpresxz(nx,ny,nz,ngrd),wrkqpresyz(nx,ny,nz,ngrd), &
    !
    wrkhrho(nx,ny,nz,ngrd),wrkhpx(nx,ny,nz,ngrd), &
    wrkhpy(nx,ny,nz,ngrd),wrkhpz(nx,ny,nz,ngrd), &
    wrkhpresx(nx,ny,nz,ngrd),wrkhpresy(nx,ny,nz,ngrd), &
    wrkhpresz(nx,ny,nz,ngrd),wrkhpresxy(nx,ny,nz,ngrd), &
    wrkhpresxz(nx,ny,nz,ngrd),wrkhpresyz(nx,ny,nz,ngrd), &
    !
    wrkorho(nx,ny,nz,ngrd),wrkopx(nx,ny,nz,ngrd), &
    wrkopy(nx,ny,nz,ngrd),wrkopz(nx,ny,nz,ngrd), &
    wrkopresx(nx,ny,nz,ngrd),wrkopresy(nx,ny,nz,ngrd), &
    wrkopresz(nx,ny,nz,ngrd),wrkopresxy(nx,ny,nz,ngrd), &
    wrkopresxz(nx,ny,nz,ngrd),wrkopresyz(nx,ny,nz,ngrd), &
    !
    wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
    wrkbz(nx,ny,nz,ngrd),wrkepres(nx,ny,nz,ngrd)
    !
    !
    dimension oldqrho(nx,ny,nz,ngrd),oldqpx(nx,ny,nz,ngrd), &
    oldqpy(nx,ny,nz,ngrd),oldqpz(nx,ny,nz,ngrd), &
    oldqpresx(nx,ny,nz,ngrd),oldqpresy(nx,ny,nz,ngrd), &
    oldqpresz(nx,ny,nz,ngrd),oldqpresxy(nx,ny,nz,ngrd), &
    oldqpresxz(nx,ny,nz,ngrd),oldqpresyz(nx,ny,nz,ngrd), &
    !
    oldhrho(nx,ny,nz,ngrd),oldhpx(nx,ny,nz,ngrd), &
    oldhpy(nx,ny,nz,ngrd),oldhpz(nx,ny,nz,ngrd), &
    oldhpresx(nx,ny,nz,ngrd),oldhpresy(nx,ny,nz,ngrd), &
    oldhpresz(nx,ny,nz,ngrd),oldhpresxy(nx,ny,nz,ngrd), &
    oldhpresxz(nx,ny,nz,ngrd),oldhpresyz(nx,ny,nz,ngrd), &
    !
    oldorho(nx,ny,nz,ngrd),oldopx(nx,ny,nz,ngrd), &
    oldopy(nx,ny,nz,ngrd),oldopz(nx,ny,nz,ngrd), &
    oldopresx(nx,ny,nz,ngrd),oldopresy(nx,ny,nz,ngrd), &
    oldopresz(nx,ny,nz,ngrd),oldopresxy(nx,ny,nz,ngrd), &
    oldopresxz(nx,ny,nz,ngrd),oldopresyz(nx,ny,nz,ngrd), &
    !
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
    oldbz(nx,ny,nz,ngrd),oldepres(nx,ny,nz,ngrd)
    !
    !
    !      write(6,*)'in flux correct with',nx,ny,nz,ngrd,m,chirho,
    !    +              diferg,xspac
    !
    rx=xspac(m)
    ry=xspac(m)
    rz=xspac(m)
    chifcs=difrho
    call fcsmooth(qrho,oldqrho,wrkqrho,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    call fcsmooth(qpresx,oldqpresx,wrkqpresx,nx,ny,nz,ngrd,m, &
        chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresy,oldqpresy,wrkqpresy,nx,ny,nz,ngrd,m, &
        chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresz,oldqpresz,wrkqpresz,nx,ny,nz,ngrd,m, &
        chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresxy,oldqpresxy,wrkqpresxy,nx,ny,nz,ngrd, &
        m,chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresxz,oldqpresxz,wrkqpresxz,nx,ny,nz,ngrd, &
        m,chifcs,vvx,vvy,vvz)
    call fcsmooth(qpresyz,oldqpresyz,wrkqpresyz,nx,ny,nz,ngrd, &
        m,chifcs,vvx,vvy,vvz)
    !
    call fcsmooth(qpx,oldqpx,wrkqpx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(qpy,oldqpy,wrkqpy,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(qpz,oldqpz,wrkqpz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    call fcsmooth(hrho,oldhrho,wrkhrho,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    call fcsmooth(hpresx,oldhpresx,wrkhpresx,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresy,oldhpresy,wrkhpresy,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresz,oldhpresz,wrkhpresz,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresxy,oldhpresxy,wrkhpresxy,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresxz,oldhpresxz,wrkhpresxz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(hpresyz,oldhpresyz,wrkhpresyz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    !
    call fcsmooth(hpx,oldhpx,wrkhpx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(hpy,oldhpy,wrkhpy,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(hpz,oldhpz,wrkhpz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    call fcsmooth(orho,oldorho,wrkorho,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    call fcsmooth(opresx,oldopresx,wrkopresx,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(opresy,oldopresy,wrkopresy,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(opresz,oldopresz,wrkopresz,nx,ny,nz,ngrd,m, &
    chifcs,vvx,vvy,vvz)
    call fcsmooth(opresxy,oldopresxy,wrkopresxy,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(opresxz,oldopresxz,wrkopresxz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    call fcsmooth(opresyz,oldopresyz,wrkopresyz,nx,ny,nz,ngrd, &
    m,chifcs,vvx,vvy,vvz)
    !
    call fcsmooth(opx,oldopx,wrkopx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(opy,oldopy,wrkopy,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(opz,oldopz,wrkopz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    !
    call fcsmooth(epres,oldepres,wrkepres,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    chifcs=diferg
    call fcsmooth(bx,oldbx,wrkbx,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(by,oldby,wrkby,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    call fcsmooth(bz,oldbz,wrkbz,nx,ny,nz,ngrd,m,chifcs, &
    vvx,vvy,vvz)
    !
    return
end
