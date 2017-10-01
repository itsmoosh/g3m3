!
!	This file contains 7 subroutines:
!	bndry_corer
!	bndry_flanks
!	bndry_grd_core
!	bndry_grds
!	bndry_inner
!	bndry_moon
!	bndry_outer
!
subroutine bndry_corer( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz,opx,opy,opz, &
    epres,bx,by,bz, &
    qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz, &
    opresxy,opresxz,opresyz, &
    nx,ny,nz,ngrd,m, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    dimension &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    call corer(qrho,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresxy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresxz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresyz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call corer(hrho,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresxy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresxz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresyz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call corer(orho,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresxy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresxz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresyz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opy,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(epres,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(bx,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(by,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(bz,nx,ny,nz,ngrd,m,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    return
end
!
!
!	**************************************************
!
!
subroutine bndry_flanks( &
    wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
    wrkqpresx,wrkqpresy,wrkqpresz, &
    wrkqpresxy,wrkqpresxz,wrkqpresyz, &
    wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
    wrkhpresx,wrkhpresy,wrkhpresz, &
    wrkhpresxy,wrkhpresxz,wrkhpresyz, &
    wrkorho,wrkopx,wrkopy,wrkopz, &
    wrkopresx,wrkopresy,wrkopresz, &
    wrkopresxy,wrkopresxz,wrkopresyz, &
    wrkepres, wrkbx,wrkby,wrkbz, &
    qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz, &
    hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz, &
    orho,opx,opy,opz,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz, &
    epres, bx,by,bz, &
    oldqrho,oldqpx,oldqpy,oldqpz, &
    oldqpresx,oldqpresy,oldqpresz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    oldhrho,oldhpx,oldhpy,oldhpz, &
    oldhpresx,oldhpresy,oldhpresz, &
    oldhpresxy,oldhpresxz,oldhpresyz, &
    oldorho,oldopx,oldopy,oldopz, &
    oldopresx,oldopresy,oldopresz, &
    oldopresxy,oldopresxz,oldopresyz, &
    oldepres, oldbx,oldby,oldbz,work, &
    nx,ny,nz,ngrd,m,t_old,t_new,t, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !     this subroutine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    dimension wrkqrho(nx,ny,nz,ngrd),wrkqpresx(nx,ny,nz,ngrd), &
    wrkqpresy(nx,ny,nz,ngrd),wrkqpresz(nx,ny,nz,ngrd), &
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd), &
    wrkhrho(nx,ny,nz,ngrd),wrkhpresx(nx,ny,nz,ngrd), &
    wrkhpresy(nx,ny,nz,ngrd),wrkhpresz(nx,ny,nz,ngrd), &
    wrkhpx(nx,ny,nz,ngrd),wrkhpy(nx,ny,nz,ngrd), &
    wrkhpz(nx,ny,nz,ngrd), &
    wrkorho(nx,ny,nz,ngrd),wrkopresx(nx,ny,nz,ngrd), &
    wrkopresy(nx,ny,nz,ngrd),wrkopresz(nx,ny,nz,ngrd), &
    wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
    wrkopz(nx,ny,nz,ngrd), &
    wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
    wrkbz(nx,ny,nz,ngrd), &
    wrkepres(nx,ny,nz,ngrd)
    dimension oldqrho(nx,ny,nz,ngrd),oldqpresx(nx,ny,nz,ngrd), &
    oldqpresy(nx,ny,nz,ngrd),oldqpresz(nx,ny,nz,ngrd), &
    oldqpx(nx,ny,nz,ngrd),oldqpy(nx,ny,nz,ngrd), &
    oldqpz(nx,ny,nz,ngrd), &
    oldhrho(nx,ny,nz,ngrd),oldhpresx(nx,ny,nz,ngrd), &
    oldhpresy(nx,ny,nz,ngrd),oldhpresz(nx,ny,nz,ngrd), &
    oldhpx(nx,ny,nz,ngrd),oldhpy(nx,ny,nz,ngrd), &
    oldhpz(nx,ny,nz,ngrd), &
    oldorho(nx,ny,nz,ngrd),oldopresx(nx,ny,nz,ngrd), &
    oldopresy(nx,ny,nz,ngrd),oldopresz(nx,ny,nz,ngrd), &
    oldopx(nx,ny,nz,ngrd),oldopy(nx,ny,nz,ngrd), &
    oldopz(nx,ny,nz,ngrd), &
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
    oldbz(nx,ny,nz,ngrd), &
    oldepres(nx,ny,nz,ngrd)
    !
    dimension qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd),hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd)
    !
    dimension oldqpresxy(nx,ny,nz,ngrd),oldqpresxz(nx,ny,nz,ngrd), &
    oldqpresyz(nx,ny,nz,ngrd),oldhpresxy(nx,ny,nz,ngrd), &
    oldhpresxz(nx,ny,nz,ngrd),oldhpresyz(nx,ny,nz,ngrd), &
    oldopresxy(nx,ny,nz,ngrd),oldopresxz(nx,ny,nz,ngrd), &
    oldopresyz(nx,ny,nz,ngrd)
    dimension wrkqpresxy(nx,ny,nz,ngrd),wrkqpresxz(nx,ny,nz,ngrd), &
    wrkqpresyz(nx,ny,nz,ngrd),wrkhpresxy(nx,ny,nz,ngrd), &
    wrkhpresxz(nx,ny,nz,ngrd),wrkhpresyz(nx,ny,nz,ngrd), &
    wrkopresxy(nx,ny,nz,ngrd),wrkopresxz(nx,ny,nz,ngrd), &
    wrkopresyz(nx,ny,nz,ngrd)
    !
    dimension work(nx,ny,nz)
    !
    dimension t_old(ngrd),t_new(ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    call flanks(wrkqrho,qrho,oldqrho,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresx,qpresx,oldqpresx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresy,qpresy,oldqpresy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresz,qpresz,oldqpresz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresxy,qpresxy,oldqpresxy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresxz,qpresxz,oldqpresxz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresyz,qpresyz,oldqpresyz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpx,qpx,oldqpx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpy,qpy,oldqpy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpz,qpz,oldqpz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkhrho,hrho,oldhrho,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresx,hpresx,oldhpresx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresy,hpresy,oldhpresy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresz,hpresz,oldhpresz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresxy,hpresxy,oldhpresxy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresxz,hpresxz,oldhpresxz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresyz,hpresyz,oldhpresyz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpx,hpx,oldhpx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpy,hpy,oldhpy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpz,hpz,oldhpz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkorho,orho,oldorho,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresx,opresx,oldopresx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresy,opresy,oldopresy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresz,opresz,oldopresz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresxy,opresxy,oldopresxy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresxz,opresxz,oldopresxz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresyz,opresyz,oldopresyz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopx,opx,oldopx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopy,opy,oldopy,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopz,opz,oldopz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkepres,epres,oldepres,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkbx,bx,oldbx,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkby,by,oldby,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkbz,bz,oldbz,nx,ny,nz,ngrd,m, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    return
end
!
!
!	**************************************************
!
!
subroutine bndry_grd_core( &
    qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,opx,opy,opz, &
    epres,bx,by,bz, &
    nx,ny,nz,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    qrho_n,qpresx_n,qpresy_n,qpresz_n, &
    qpresxy_n,qpresxz_n,qpresyz_n, &
    qpx_n,qpy_n,qpz_n, &
    hrho_n,hpresx_n,hpresy_n,hpresz_n, &
    hpresxy_n,hpresxz_n,hpresyz_n, &
    hpx_n,hpy_n,hpz_n, &
    orho_n,opresx_n,opresy_n,opresz_n, &
    opresxy_n,opresxz_n,opresyz_n, &
    opx_n,opy_n,opz_n, &
    epres_n,bx_n,by_n,bz_n, &
    nx_n,ny_n,nz_n,ngrd_n,main_ngrd, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    !     this subroutine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    dimension &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    ! 
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd)
    dimension qrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpy_n(nx_n,ny_n,nz_n,ngrd_n),qpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    hrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpx_n(nx_n,ny_n,nz_n,ngrd_n),hpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    orho_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opx_n(nx_n,ny_n,nz_n,ngrd_n),opy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    bx_n(nx_n,ny_n,nz_n,ngrd_n),by_n(nx_n,ny_n,nz_n,ngrd_n), &
    bz_n(nx_n,ny_n,nz_n,ngrd_n), &
    epres_n(nx_n,ny_n,nz_n,ngrd_n)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    call corer_grds(qrho,nx,ny,nz,ngrd,main_ngrd, &
    qrho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresx,nx,ny,nz,ngrd,main_ngrd, &
    qpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresy,nx,ny,nz,ngrd,main_ngrd, &
    qpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresz,nx,ny,nz,ngrd,main_ngrd, &
    qpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresxy,nx,ny,nz,ngrd,main_ngrd, &
    qpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresxz,nx,ny,nz,ngrd,main_ngrd, &
    qpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresyz,nx,ny,nz,ngrd,main_ngrd, &
    qpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpx,nx,ny,nz,ngrd,main_ngrd, &
    qpx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpy,nx,ny,nz,ngrd,main_ngrd, &
    qpy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpz,nx,ny,nz,ngrd,main_ngrd, &
    qpz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(hrho,nx,ny,nz,ngrd,main_ngrd, &
    hrho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresx,nx,ny,nz,ngrd,main_ngrd, &
    hpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresy,nx,ny,nz,ngrd,main_ngrd, &
    hpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresz,nx,ny,nz,ngrd,main_ngrd, &
    hpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresxy,nx,ny,nz,ngrd,main_ngrd, &
    hpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresxz,nx,ny,nz,ngrd,main_ngrd, &
    hpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresyz,nx,ny,nz,ngrd,main_ngrd, &
    hpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpx,nx,ny,nz,ngrd,main_ngrd, &
    hpx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpy,nx,ny,nz,ngrd,main_ngrd, &
    hpy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpz,nx,ny,nz,ngrd,main_ngrd, &
    hpz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(orho,nx,ny,nz,ngrd,main_ngrd, &
    orho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresx,nx,ny,nz,ngrd,main_ngrd, &
    opresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresy,nx,ny,nz,ngrd,main_ngrd, &
    opresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresz,nx,ny,nz,ngrd,main_ngrd, &
    opresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresxy,nx,ny,nz,ngrd,main_ngrd, &
    opresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresxz,nx,ny,nz,ngrd,main_ngrd, &
    opresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresyz,nx,ny,nz,ngrd,main_ngrd, &
    opresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opx,nx,ny,nz,ngrd,main_ngrd, &
    opx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opy,nx,ny,nz,ngrd,main_ngrd, &
    opy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opz,nx,ny,nz,ngrd,main_ngrd, &
    opz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(epres,nx,ny,nz,ngrd,main_ngrd, &
    epres_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(bx,nx,ny,nz,ngrd,main_ngrd, &
    bx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(by,nx,ny,nz,ngrd,main_ngrd, &
    by_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(bz,nx,ny,nz,ngrd,main_ngrd, &
    bz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    return
end
!
!
!	**************************************************
!
!
subroutine bndry_grds( &
    qrho_n,qpresx_n,qpresy_n,qpresz_n, &
    qpresxy_n,qpresxz_n,qpresyz_n,qpx_n,qpy_n,qpz_n, &
    hrho_n,hpresx_n,hpresy_n,hpresz_n, &
    hpresxy_n,hpresxz_n,hpresyz_n,hpx_n,hpy_n,hpz_n, &
    orho_n,opresx_n,opresy_n,opresz_n, &
    opresxy_n,opresxz_n,opresyz_n,opx_n,opy_n,opz_n, &
    epres_n,bx_n,by_n,bz_n, &
    nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n,t, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n, &
    qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,opx,opy,opz, &
    epres,bx,by,bz, &
    oldqrho,oldqpresx,oldqpresy,oldqpresz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    oldqpx,oldqpy,oldqpz, &
    oldhrho,oldhpresx,oldhpresy,oldhpresz, &
    oldhpresxy,oldhpresxz,oldhpresyz, &
    oldhpx,oldhpy,oldhpz, &
    oldorho,oldopresx,oldopresy,oldopresz, &
    oldopresxy,oldopresxz,oldopresyz, &
    oldopx,oldopy,oldopz, &
    oldepres,oldbx,oldby,oldbz,work, &
    nx,ny,nz,ngrd,t_old,t_new, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !     this subroutine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    dimension qrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    hrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    orho_n(nx_n,ny_n,nz_n,ngrd_n), &
    opx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    bx_n(nx_n,ny_n,nz_n,ngrd_n),by_n(nx_n,ny_n,nz_n,ngrd_n), &
    bz_n(nx_n,ny_n,nz_n,ngrd_n),epres_n(nx_n,ny_n,nz_n,ngrd_n)
    !
    dimension qrho(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd),qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd),hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    !
    orho(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresy(nx,ny,nz,ngrd), &
    opresz(nx,ny,nz,ngrd),opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    dimension oldqrho(nx,ny,nz,ngrd), &
    oldqpresx(nx,ny,nz,ngrd),oldqpresy(nx,ny,nz,ngrd), &
    oldqpresz(nx,ny,nz,ngrd),oldqpresxy(nx,ny,nz,ngrd), &
    oldqpresxz(nx,ny,nz,ngrd),oldqpresyz(nx,ny,nz,ngrd), &
    oldqpx(nx,ny,nz,ngrd),oldqpy(nx,ny,nz,ngrd), &
    oldqpz(nx,ny,nz,ngrd), &
    !
    oldhrho(nx,ny,nz,ngrd), &
    oldhpresx(nx,ny,nz,ngrd),oldhpresy(nx,ny,nz,ngrd), &
    oldhpresz(nx,ny,nz,ngrd),oldhpresxy(nx,ny,nz,ngrd), &
    oldhpresxz(nx,ny,nz,ngrd),oldhpresyz(nx,ny,nz,ngrd), &
    oldhpx(nx,ny,nz,ngrd),oldhpy(nx,ny,nz,ngrd), &
    oldhpz(nx,ny,nz,ngrd), &
    !
    oldorho(nx,ny,nz,ngrd), &
    oldopresx(nx,ny,nz,ngrd),oldopresy(nx,ny,nz,ngrd), &
    oldopresz(nx,ny,nz,ngrd),oldopresxy(nx,ny,nz,ngrd), &
    oldopresxz(nx,ny,nz,ngrd),oldopresyz(nx,ny,nz,ngrd), &
    oldopx(nx,ny,nz,ngrd),oldopy(nx,ny,nz,ngrd), &
    oldopz(nx,ny,nz,ngrd), &
    oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
    oldbz(nx,ny,nz,ngrd),oldepres(nx,ny,nz,ngrd)
    !
    dimension work(nx,ny,nz)
    dimension t_old(ngrd),t_new(ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    call flanks_grds(qrho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qrho,oldqrho,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresx,oldqpresx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresy,oldqpresy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresz,oldqpresz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresxy,oldqpresxy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresxz,oldqpresxz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresyz,oldqpresyz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpx,oldqpx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpy,oldqpy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpz,oldqpz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(hrho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hrho,oldhrho,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresx,oldhpresx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresy,oldhpresy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresz,oldhpresz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresxy,oldhpresxy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresxz,oldhpresxz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresyz,oldhpresyz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpx,oldhpx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpy,oldhpy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpz,oldhpz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(orho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    orho,oldorho,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresx,oldopresx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresy,oldopresy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresz,oldopresz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresxy,oldopresxy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresxz,oldopresxz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresyz,oldopresyz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opx,oldopx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opy,oldopy,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opz,oldopz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(epres_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    epres,oldepres,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(bx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    bx,oldbx,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(by_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    by,oldby,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(bz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    bz,oldbz,work,nx,ny,nz,ngrd, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    return
end
!
!
!	**************************************************
!
!
subroutine bndry_inner( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
    orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
    epres,bx,by,bz, &
    qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz, &
    opresxy,opresxz,opresyz, &
    nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
    ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
    mbndry,msrf,mmid,mzero, &
    erho,epress,alpha_e,ti_te,o_conc, &
    re_equiv,rearth,sbx_wind,spress, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !     this subroutine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    dimension &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    ! 
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real  parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    !     set surface conditions around earth
    !      withinterior temperature as constant
    !
    !     write(6,*)'inside bndry inner with',mbndry
    !     write(6,*)  numsrf,nummid,numzero
    !
    aheight=rearth+2.
    !
    m=1
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    r_rot2=r_rot**2
    !
    !     set ionospheric boundary points - plasma quantities defined
    !
    do n=1,numsrf(m)
        !
        !        get coords of point
        !
        i=ijsrf(m,1,n)
        j=ijsrf(m,2,n)
        k=ijsrf(m,3,n)
        !
        !         set rotational speeds
        !
        ay=grd_ymin(m)+dy*(j-1)
        ax=grd_xmin(m)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        !        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        !        qden=parm_srf(1,n)
        qden=amax1(parm_srf(m,1,n),qrho(i,j,k,m))
        hden=amax1(parm_srf(m,2,n),hrho(i,j,k,m))
        oden=amax1(parm_srf(m,3,n),orho(i,j,k,m))
        !
        qrho(i,j,k,m)=qden
        hrho(i,j,k,m)=hden
        orho(i,j,k,m)=oden
        !
        qpresx(i,j,k,m)=parm_srf(m,4,n)
        hpresx(i,j,k,m)=parm_srf(m,5,n)
        opresx(i,j,k,m)=parm_srf(m,6,n)
        !
        qpresy(i,j,k,m)=parm_srf(m,4,n)
        hpresy(i,j,k,m)=parm_srf(m,5,n)
        opresy(i,j,k,m)=parm_srf(m,6,n)
        !
        qpresz(i,j,k,m)=parm_srf(m,4,n)
        hpresz(i,j,k,m)=parm_srf(m,5,n)
        opresz(i,j,k,m)=parm_srf(m,6,n)
        !
        epres(i,j,k,m)=parm_srf(m,7,n)
        !
        hpx(i,j,k,m)=hden*rvx
        hpy(i,j,k,m)=hden*rvy
        hpz(i,j,k,m)=hden*rvz
        opx(i,j,k,m)=oden*rvx
        opy(i,j,k,m)=oden*rvy
        opz(i,j,k,m)=oden*rvz
        qpx(i,j,k,m)=qden*rvx
        qpy(i,j,k,m)=qden*rvy
        qpz(i,j,k,m)=qden*rvz
        !
        !        iostropic conditions: set off axis pressures to zero
        !
        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
    enddo
    !
    !     set atmospheric boundary points - specify magnetic field profile
    !
    do n=1,nummid(m)
        !
        !        get coords of point
        !
        i=ijmid(m,1,n)
        j=ijmid(m,2,n)
        k=ijmid(m,3,n)
        !
        qrho(i,j,k,m)=parm_mid(m,1,n)
        hrho(i,j,k,m)=parm_mid(m,2,n)
        orho(i,j,k,m)=parm_mid(m,3,n)
        !
        qpresx(i,j,k,m)=parm_mid(m,4,n)
        hpresx(i,j,k,m)=parm_mid(m,5,n)
        opresx(i,j,k,m)=parm_mid(m,6,n)
        !
        qpresy(i,j,k,m)=parm_mid(m,4,n)
        hpresy(i,j,k,m)=parm_mid(m,5,n)
        opresy(i,j,k,m)=parm_mid(m,6,n)
        !
        qpresz(i,j,k,m)=parm_mid(m,4,n)
        hpresz(i,j,k,m)=parm_mid(m,5,n)
        opresz(i,j,k,m)=parm_mid(m,6,n)
        !
        epres(i,j,k,m)=parm_mid(m,7,n)
        !
        ay=grd_ymin(m)+dy*(j-1)
        ax=grd_xmin(m)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        !        bx(i,j,k,m)=sbx_wind
        !        by(i,j,k,m)=0.
        !        bz(i,j,k,m)=0.
        !
    enddo
    !
    !      reset interior points
    !
    do n=1,numzero(m)
        i=ijzero(m,1,n)
        j=ijzero(m,2,n)
        k=ijzero(m,3,n)
        !
        qrho(i,j,k,m)=parm_zero(m,1,n)
        hrho(i,j,k,m)=parm_zero(m,2,n)
        orho(i,j,k,m)=parm_zero(m,3,n)
        !
        qpresx(i,j,k,m)=parm_zero(m,4,n)
        hpresx(i,j,k,m)=parm_zero(m,5,n)
        opresx(i,j,k,m)=parm_zero(m,6,n)
        !
        qpresy(i,j,k,m)=parm_zero(m,4,n)
        hpresy(i,j,k,m)=parm_zero(m,5,n)
        opresy(i,j,k,m)=parm_zero(m,6,n)
        !
        qpresz(i,j,k,m)=parm_zero(m,4,n)
        hpresz(i,j,k,m)=parm_zero(m,5,n)
        opresz(i,j,k,m)=parm_zero(m,6,n)
        !
        epres(i,j,k,m)=parm_zero(m,7,n)
        !
        ay=grd_ymin(m)+dy*(j-1)
        ax=grd_xmin(m)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        bx(i,j,k,m)=0.
        by(i,j,k,m)=0.
        bz(i,j,k,m)=0.
        !
    enddo
    !
    return
end
!
!
!	**************************************************
!
!
subroutine bndry_moon(qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz, &
    qpx,qpy,qpz,rmassq, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz, &
    hpx,hpy,hpz,rmassh, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz, &
    opx,opy,opz,rmasso,epres,bx,by,bz, &
    m,nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
    ijsrf,numsrf,ijmid,nummid,ijzero, &
    numzero,mbndry,msrf,mmid,mzero, &
    qden_moon,hden_moon,oden_moon,cs_moon,gamma,ti_te, &
    vx_moon,vy_moon,vz_moon, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    isotropic)
    !
    !     this subroutine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
    tempi=cs_moon**2/gamma
    alpha_m=6.
    rvx=vx_moon
    rvy=vy_moon
    rvz=vz_moon
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !     set atmospheric boundary points -- specify magnetic field profile
    !
    do n=1,nummid(m)
        !
        !        get coords of point
        !
        i=ijmid(m,1,n)
        j=ijmid(m,2,n)
        k=ijmid(m,3,n)
        !
        qrho(i,j,k,m)=parm_mid(m,1,n)
        hrho(i,j,k,m)=parm_mid(m,2,n)
        orho(i,j,k,m)=parm_mid(m,3,n)
        !
        qpresx(i,j,k,m)=parm_mid(m,4,n)
        qpresy(i,j,k,m)=qpresx(i,j,k,m)
        qpresz(i,j,k,m)=qpresx(i,j,k,m)*0.9

        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        !
        hpresx(i,j,k,m)=parm_mid(m,5,n)
        hpresy(i,j,k,m)=hpresx(i,j,k,m)
        hpresz(i,j,k,m)=hpresx(i,j,k,m)*0.9
		!
        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        !
        opresx(i,j,k,m)=parm_mid(m,6,n)
        opresy(i,j,k,m)=opresx(i,j,k,m)
        opresz(i,j,k,m)=opresx(i,j,k,m)*0.9
		!
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        epres(i,j,k,m)=parm_mid(m,7,n)
        !
        !        qrho(i,j,k,m)=amax1(0.25*qrho(i,j,k,m),parm_mid(m,1,n))
        !        hrho(i,j,k,m)=amax1(0.25*hrho(i,j,k,m),parm_mid(m,2,n))
        !        orho(i,j,k,m)=amax1(0.25*orho(i,j,k,m),parm_mid(m,3,n))
        !
        !        qpres(i,j,k,m)=amax1(0.25*qpres(i,j,k,m),parm_mid(m,4,n))
        !        hpres(i,j,k,m)=amax1(0.25*hpres(i,j,k,m),parm_mid(m,5,n))
        !        opres(i,j,k,m)=amax1(0.25*opres(i,j,k,m),parm_mid(m,6,n))
        !        epres(i,j,k,m)=amax1(0.25*epres(i,j,k,m),parm_mid(m,7,n))
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
        !        bx(i,j,k,m)=sbx_wind
        !        by(i,j,k,m)=0.
        !        bz(i,j,k,m)=0.
        !
    enddo
    !
    !      reset interior points
    !
    do n=1,numzero(m)
        i=ijzero(m,1,n)
        j=ijzero(m,2,n)
        k=ijzero(m,3,n)
        !
        qrho(i,j,k,m)=parm_zero(m,1,n)
        hrho(i,j,k,m)=parm_zero(m,2,n)
        orho(i,j,k,m)=parm_zero(m,3,n)
        !
        qpresx(i,j,k,m)=parm_zero(m,4,n)
        qpresy(i,j,k,m)=qpresx(i,j,k,m)
        qpresz(i,j,k,m)=qpresx(i,j,k,m)*0.9
		!
        qpresxy(i,j,k,m)=0.
        qpresxz(i,j,k,m)=0.
        qpresyz(i,j,k,m)=0.
        !
        hpresx(i,j,k,m)=parm_zero(m,5,n)
        hpresy(i,j,k,m)=hpresx(i,j,k,m)
        hpresz(i,j,k,m)=hpresx(i,j,k,m)*0.9

        hpresxy(i,j,k,m)=0.
        hpresxz(i,j,k,m)=0.
        hpresyz(i,j,k,m)=0.
        !
        opresx(i,j,k,m)=parm_zero(m,6,n)
        opresy(i,j,k,m)=opresx(i,j,k,m)
        opresz(i,j,k,m)=opresx(i,j,k,m)*0.9
		!
        opresxy(i,j,k,m)=0.
        opresxz(i,j,k,m)=0.
        opresyz(i,j,k,m)=0.
        !
        epres(i,j,k,m)=parm_zero(m,7,n)
        !
        !        qrho(i,j,k,m)=amax1(0.10*qrho(i,j,k,m),parm_zero(m,1,n))
        !        hrho(i,j,k,m)=amax1(0.10*hrho(i,j,k,m),parm_zero(m,2,n))
        !        orho(i,j,k,m)=amax1(0.10*orho(i,j,k,m),parm_zero(m,3,n))
        !
        !        qpres(i,j,k,m)=amax1(0.05*qpres(i,j,k,m),parm_zero(m,4,n))
        !        hpres(i,j,k,m)=amax1(0.05*hpres(i,j,k,m),parm_zero(m,5,n))
        !        opres(i,j,k,m)=amax1(0.05*opres(i,j,k,m),parm_zero(m,6,n))
        !        epres(i,j,k,m)=amax1(0.05*epres(i,j,k,m),parm_zero(m,7,n))
        !
        opx(i,j,k,m)=orho(i,j,k,m)*rvx
        opy(i,j,k,m)=orho(i,j,k,m)*rvy
        opz(i,j,k,m)=orho(i,j,k,m)*rvz
        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
        hpz(i,j,k,m)=hrho(i,j,k,m)*rvz
        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
        qpz(i,j,k,m)=qrho(i,j,k,m)*rvz
        !
    enddo
    !
    return
end
!
!
!	**************************************************
!
!
subroutine bndry_outer( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
    orho,opresx,opresy,opresz,opx,opy,opz, &
    epres, &
    qpresxy,qpresxz,qpresyz, &
    hpresxy,hpresxz,hpresyz, &
    opresxy,opresxz,opresyz, &
    rmassq,rmassh,rmasso, &
    bx,by,bz,bx0,by0,bz0,nx,ny,nz,ngrd, &
    srho,rho_frac,o_conc,spress,spx,spy,spz, &
    sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)
    !
    !     this subroutine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    dimension &
    qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
    qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
    qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
    qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
    !
    hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
    hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
    hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
    hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
    ! 
    orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
    opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
    opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd), &
    !
    bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
    !
    !     do x boundary conditions: set oxygen in solar wind
    !       at 2.5%
    !
    m=ngrd
    ofrac=rho_frac*o_conc
    frac_h=rmassh/rmassq
    frac_o=rmasso/rmassq
    !
    i=1
    ii=2
    do k=1,nz
        do j=1,ny
            !
            !       wind boundary conditions at i=1
            !
            qrho(i,j,k,m)=srho
            qpx(i,j,k,m)=spx
            qpy(i,j,k,m)=spy
            qpz(i,j,k,m)=spz
            qpresx(i,j,k,m)=0.5*spress
            qpresy(i,j,k,m)=qpresx(i,j,k,m)
            qpresz(i,j,k,m)=qpresx(i,j,k,m)
            !
            hrho(i,j,k,m)=srho*rho_frac*frac_h
            hpx(i,j,k,m)=spx*rho_frac*frac_h
            hpy(i,j,k,m)=spy*rho_frac*frac_h
            hpz(i,j,k,m)=spz*rho_frac*frac_h
            hpresx(i,j,k,m)=0.5*spress*rho_frac
            hpresy(i,j,k,m)=hpresx(i,j,k,m)
            hpresz(i,j,k,m)=hpresx(i,j,k,m)
            !
            orho(i,j,k,m)=srho*ofrac*frac_o
            opx(i,j,k,m)=spx*ofrac*frac_o
            opy(i,j,k,m)=spy*ofrac*frac_o
            opz(i,j,k,m)=spz*ofrac*frac_o
            opresx(i,j,k,m)=0.5*spress*ofrac
            opresy(i,j,k,m)=opresx(i,j,k,m)
            opresz(i,j,k,m)=opresx(i,j,k,m)
            !
            epres(i,j,k,m)=qpresx(i,j,k,m)/ti_te
            !
            !        iostropic conditions: set off axis pressures to zero
            !
            qpresxy(i,j,k,m)=0.
            qpresxz(i,j,k,m)=0.
            qpresyz(i,j,k,m)=0.
            hpresxy(i,j,k,m)=0.
            hpresxz(i,j,k,m)=0.
            hpresyz(i,j,k,m)=0.
            opresxy(i,j,k,m)=0.
            opresxz(i,j,k,m)=0.
            opresyz(i,j,k,m)=0.
            !
            !       set bx of wind to zero, by to the wind values
            !          and subtract any geomagnetic field
            !
            bx(i,j,k,m)=sbx_wind-bx0(i,j,k,m)
            by(i,j,k,m)=sby_wind-by0(i,j,k,m)
            bz(i,j,k,m)=sbz_wind-bz0(i,j,k,m)
            !
            !       bx(i,j,k,m)=bx(2,j,k,m)
            !       by(i,j,k,m)=by(2,j,k,m)
            !       bz(i,j,k,m)=bz(2,j,k,m)
            !
        enddo
    enddo
    !
    !     do x boundary conditions - back wall
    !
    nx1=nx-1
    do k=1,nz
        do j=1,ny
            !
            qrho(nx,j,k,m)=qrho(nx1,j,k,m)
            qpx(nx,j,k,m)=abs(qpx(nx1,j,k,m))
            qpy(nx,j,k,m)=qpy(nx1,j,k,m)
            qpz(nx,j,k,m)=qpz(nx1,j,k,m)
            qpresx(nx,j,k,m)=qpresx(nx1,j,k,m)
            qpresy(nx,j,k,m)=qpresy(nx1,j,k,m)
            qpresz(nx,j,k,m)=qpresz(nx1,j,k,m)
            !
            hrho(nx,j,k,m)=hrho(nx1,j,k,m)
            hpx(nx,j,k,m)=abs(hpx(nx1,j,k,m))
            hpy(nx,j,k,m)=hpy(nx1,j,k,m)
            hpz(nx,j,k,m)=hpz(nx1,j,k,m)
            hpresx(nx,j,k,m)=hpresx(nx1,j,k,m)
            hpresy(nx,j,k,m)=hpresy(nx1,j,k,m)
            hpresz(nx,j,k,m)=hpresz(nx1,j,k,m)
            !
            orho(nx,j,k,m)=orho(nx1,j,k,m)
            opx(nx,j,k,m)=abs(opx(nx1,j,k,m))
            opy(nx,j,k,m)=opy(nx1,j,k,m)
            opz(nx,j,k,m)=opz(nx1,j,k,m)
            opresx(nx,j,k,m)=opresx(nx1,j,k,m)
	   		!
            opresy(nx,j,k,m)=opresy(nx1,j,k,m)
            opresz(nx,j,k,m)=opresz(nx1,j,k,m)
            !
            epres(nx,j,k,m)=epres(nx1,j,k,m)
            !
            !       continuous boundary condition
            !
            qpresxy(nx,j,k,m)=qpresxy(nx1,j,k,m)
            qpresxz(nx,j,k,m)=qpresxz(nx1,j,k,m)
            qpresyz(nx,j,k,m)=qpresyz(nx1,j,k,m)
            hpresxy(nx,j,k,m)=hpresxy(nx1,j,k,m)
            hpresxz(nx,j,k,m)=hpresxz(nx1,j,k,m)
            hpresyz(nx,j,k,m)=hpresyz(nx1,j,k,m)
            opresxy(nx,j,k,m)=opresxy(nx1,j,k,m)
            opresxz(nx,j,k,m)=opresxz(nx1,j,k,m)
            opresyz(nx,j,k,m)=opresyz(nx1,j,k,m)
            !
            bx(nx,j,k,m)=bx(nx1,j,k,m)
            by(nx,j,k,m)=by(nx1,j,k,m)
            bz(nx,j,k,m)=bz(nx1,j,k,m)
        enddo
    enddo
    !
    !     do y boundary conditions
    !
    ny1=ny-1
    do k=1,nz
        do i=2,nx
            i1=i
            !
            qrho(i,1,k,m)=qrho(i1,2,k,m)
            qpresx(i,1,k,m)=qpresx(i1,2,k,m)
            qpresy(i,1,k,m)=qpresy(i1,2,k,m)
            qpresz(i,1,k,m)=qpresz(i1,2,k,m)
		    !
            qpx(i,1,k,m)=qpx(i1,2,k,m)
            qpy(i,1,k,m)=-abs(qpy(i1,2,k,m))
            qpz(i,1,k,m)=qpz(i1,2,k,m)
            !
            hrho(i,1,k,m)=hrho(i1,2,k,m)
            hpresx(i,1,k,m)=hpresx(i1,2,k,m)
            hpresy(i,1,k,m)=hpresy(i1,2,k,m)
            hpresz(i,1,k,m)=hpresz(i1,2,k,m)
		    !
            hpx(i,1,k,m)=hpx(i1,2,k,m)
            hpy(i,1,k,m)=-abs(hpy(i1,2,k,m))
            hpz(i,1,k,m)=hpz(i1,2,k,m)
            !
            orho(i,1,k,m)=orho(i1,2,k,m)
            opresx(i,1,k,m)=opresx(i1,2,k,m)
            opresy(i,1,k,m)=opresy(i1,2,k,m)
            opresz(i,1,k,m)=opresz(i1,2,k,m)
		    !
            opx(i,1,k,m)=opx(i1,2,k,m)
            opy(i,1,k,m)=-abs(opy(i1,2,k,m))
            opz(i,1,k,m)=opz(i1,2,k,m)
            !
            epres(i,1,k,m)=epres(i1,2,k,m)
            !
            !       continuous boundary condition
            !
            qpresxy(i,1,k,m)=qpresxy(i1,2,k,m)
            qpresxz(i,1,k,m)=qpresxz(i1,2,k,m)
            qpresyz(i,1,k,m)=qpresyz(i1,2,k,m)
            hpresxy(i,1,k,m)=hpresxy(i1,2,k,m)
            hpresxz(i,1,k,m)=hpresxz(i1,2,k,m)
            hpresyz(i,1,k,m)=hpresyz(i1,2,k,m)
            opresxy(i,1,k,m)=opresxy(i1,2,k,m)
            opresxz(i,1,k,m)=opresxz(i1,2,k,m)
            opresyz(i,1,k,m)=opresyz(i1,2,k,m)
            !
            by(i,1,k,m)=by(i1,2,k,m)
            !    +             +(by0(i1,2,k,m)-by0(i,1,k,m))
            bx(i,1,k,m)=bx(i1,2,k,m)
            !    +              +(bx0(i1,2,k,m)-bx0(i,1,k,m))
            bz(i,1,k,m)=bz(i1,2,k,m)
            !    +              +(bz0(i1,2,k,m)-bz0(i,1,k,m))
            !
            qrho(i,ny,k,m)=qrho(i1,ny1,k,m)
            qpresx(i,ny,k,m)=qpresx(i1,ny1,k,m)
            qpresy(i,ny,k,m)=qpresy(i1,ny1,k,m)
            qpresz(i,ny,k,m)=qpresz(i1,ny1,k,m)
			!
            qpx(i,ny,k,m)=qpx(i1,ny1,k,m)
            qpy(i,ny,k,m)=abs(qpy(i1,ny1,k,m))
            qpz(i,ny,k,m)=qpz(i1,ny1,k,m)
            !
            hrho(i,ny,k,m)=hrho(i1,ny1,k,m)
            hpresx(i,ny,k,m)=hpresx(i1,ny1,k,m)
            hpresy(i,ny,k,m)=hpresy(i1,ny1,k,m)
            hpresz(i,ny,k,m)=hpresz(i1,ny1,k,m)
            hpx(i,ny,k,m)=hpx(i1,ny1,k,m)
            hpy(i,ny,k,m)=abs(hpy(i1,ny1,k,m))
            hpz(i,ny,k,m)=hpz(i1,ny1,k,m)
            !
            orho(i,ny,k,m)=orho(i1,ny1,k,m)
            opresx(i,ny,k,m)=opresx(i1,ny1,k,m)
            opresy(i,ny,k,m)=opresy(i1,ny1,k,m)
            opresz(i,ny,k,m)=opresz(i1,ny1,k,m)
			!
            opx(i,ny,k,m)=opx(i1,ny1,k,m)
            opy(i,ny,k,m)=abs(opy(i1,ny1,k,m))
            opz(i,ny,k,m)=opz(i1,ny1,k,m)
            !
            epres(i,ny,k,m)=epres(i1,ny1,k,m)
            !
            !       continuous boundary condition
            !
            qpresxy(i,ny,k,m)=qpresxy(i1,ny1,k,m)
            qpresxz(i,ny,k,m)=qpresxz(i1,ny1,k,m)
            qpresyz(i,ny,k,m)=qpresyz(i1,ny1,k,m)
            hpresxy(i,ny,k,m)=hpresxy(i1,ny1,k,m)
            hpresxz(i,ny,k,m)=hpresxz(i1,ny1,k,m)
            hpresyz(i,ny,k,m)=hpresyz(i1,ny1,k,m)
            opresxy(i,ny,k,m)=opresxy(i1,ny1,k,m)
            opresxz(i,ny,k,m)=opresxz(i1,ny1,k,m)
            opresyz(i,ny,k,m)=opresyz(i1,ny1,k,m)
			!
            by(i,ny,k,m)=by(i1,ny1,k,m)
            !    +              +(by0(i1,ny1,k,m)-by0(i,ny,k,m))
            bx(i,ny,k,m)=bx(i1,ny1,k,m)
            !    +               +(bx0(i1,ny1,k,m)-bx0(i,ny,k,m))
            bz(i,ny,k,m)=bz(i1,ny1,k,m)
            !    +              +(bz0(i1,ny1,k,m)-bz0(i,ny,k,m))
        enddo
    enddo
    !
    !     do z boundary conditions
    !
    nz1=nz-1
    do j=1,ny
        do i=2,nx
            i1=i
            !
            qrho(i,j,nz,m)=qrho(i1,j,nz1,m)
            qpresx(i,j,nz,m)=qpresx(i1,j,nz1,m)
            qpresy(i,j,nz,m)=qpresy(i1,j,nz1,m)
            qpresz(i,j,nz,m)=qpresz(i1,j,nz1,m)
			!
            qpx(i,j,nz,m)=qpx(i1,j,nz1,m)
            qpy(i,j,nz,m)=qpy(i1,j,nz1,m)
            qpz(i,j,nz,m)=abs(qpz(i1,j,nz1,m))
            !
            hrho(i,j,nz,m)=hrho(i1,j,nz1,m)
            hpresx(i,j,nz,m)=hpresx(i1,j,nz1,m)
            hpresy(i,j,nz,m)=hpresy(i1,j,nz1,m)
            hpresz(i,j,nz,m)=hpresz(i1,j,nz1,m)
			!
            hpx(i,j,nz,m)=hpx(i1,j,nz1,m)
            hpy(i,j,nz,m)=hpy(i1,j,nz1,m)
            hpz(i,j,nz,m)=abs(hpz(i1,j,nz1,m))
            !
            orho(i,j,nz,m)=orho(i1,j,nz1,m)
            opresx(i,j,nz,m)=opresx(i1,j,nz1,m)
            opresy(i,j,nz,m)=opresy(i1,j,nz1,m)
            opresz(i,j,nz,m)=opresz(i1,j,nz1,m)
			!
            opx(i,j,nz,m)=opx(i1,j,nz1,m)
            opy(i,j,nz,m)=opy(i1,j,nz1,m)
            opz(i,j,nz,m)=abs(opz(i1,j,nz1,m))
            !
            epres(i,j,nz,m)=epres(i1,j,nz1,m)
            !
            !       off diagnoal elements
            !
            qpresxy(i,j,nz,m)=qpresxy(i1,j,nz1,m)
            qpresxz(i,j,nz,m)=qpresxz(i1,j,nz1,m)
            qpresyz(i,j,nz,m)=qpresyz(i1,j,nz1,m)
            hpresxy(i,j,nz,m)=hpresxy(i1,j,nz1,m)
            hpresxz(i,j,nz,m)=hpresxz(i1,j,nz1,m)
            hpresyz(i,j,nz,m)=hpresyz(i1,j,nz1,m)
            opresxy(i,j,nz,m)=opresxy(i1,j,nz1,m)
            opresxz(i,j,nz,m)=opresxz(i1,j,nz1,m)
            opresyz(i,j,nz,m)=opresyz(i1,j,nz1,m)
            !
            bz(i,j,nz,m)=bz(i1,j,nz1,m)
            !    +            +bz0(i1,j,nz1,m)-bz0(i,j,nz,m)
            bx(i,j,nz,m)=bx(i1,j,nz1,m)
            !    +            +bx0(i1,j,nz1,m)-bx0(i,j,nz,m)
            by(i,j,nz,m)=by(i1,j,nz1,m)
            !    +            +by0(i1,j,nz1,m)-by0(i,j,nz,m)
        enddo
    enddo
    !
    !       symmetry boundary conditions at k=1 and k=2
    !
    do j=1,ny
        do i=2,nx
            !
            qrho(i,j,1,m)=qrho(i,j,2,m)
            qpx(i,j,1,m)=qpx(i,j,2,m)
            qpy(i,j,1,m)=qpy(i,j,2,m)
            qpz(i,j,1,m)=-abs(qpz(i,j,2,m))
            qpresx(i,j,1,m)=qpresx(i,j,2,m)
            qpresy(i,j,1,m)=qpresy(i,j,2,m)
            qpresz(i,j,1,m)=qpresz(i,j,2,m)
            !
            hrho(i,j,1,m)=hrho(i,j,2,m)
            hpx(i,j,1,m)=hpx(i,j,2,m)
            hpy(i,j,1,m)=hpy(i,j,2,m)
            hpz(i,j,1,m)=-abs(hpz(i,j,2,m))
            hpresx(i,j,1,m)=hpresx(i,j,2,m)
            hpresy(i,j,1,m)=hpresy(i,j,2,m)
            hpresz(i,j,1,m)=hpresz(i,j,2,m)
            !
            orho(i,j,1,m)=orho(i,j,2,m)
            opx(i,j,1,m)=opx(i,j,2,m)
            opy(i,j,1,m)=opy(i,j,2,m)
            opz(i,j,1,m)=-abs(opz(i,j,2,m))
            opresx(i,j,1,m)=opresx(i,j,2,m)
            opresy(i,j,1,m)=opresy(i,j,2,m)
            opresz(i,j,1,m)=opresz(i,j,2,m)
            !
            epres(i,j,1,m)=epres(i,j,2,m)
            !
            !       off diagonal elements
            !
            qpresxy(i,j,1,m)=qpresxy(i,j,2,m)
            qpresxz(i,j,1,m)=qpresxz(i,j,2,m)
            qpresyz(i,j,1,m)=qpresyz(i,j,2,m)
            hpresxy(i,j,1,m)=hpresxy(i,j,2,m)
            hpresxz(i,j,1,m)=hpresxz(i,j,2,m)
            hpresyz(i,j,1,m)=hpresyz(i,j,2,m)
            opresxy(i,j,1,m)=opresxy(i,j,2,m)
            opresxz(i,j,1,m)=opresxz(i,j,2,m)
            opresyz(i,j,1,m)=opresyz(i,j,2,m)
            !
            !       set bx of wind to zero, by to the wind values
            !          and subtract any geomagnetic field
            !       force bz to be positive to attain equilibrium
            !
            bx(i,j,1,m)=bx(i,j,2,m)
            by(i,j,1,m)=by(i,j,2,m)
            bz(i,j,1,m)=bz(i,j,2,m)
        enddo
    enddo
    !
    !     do corner lines at left-right back wall
    !
    do k=2,nz-1
        qrho(nx,1,k,m)=(qrho(nx1,1,k,m)+qrho(nx,2,k,m))/2.
        qrho(nx,ny,k,m)=(qrho(nx1,ny,k,m)+qrho(nx,ny1,k,m))/2.
        hrho(nx,1,k,m)=(hrho(nx1,1,k,m)+hrho(nx,2,k,m))/2.
        hrho(nx,ny,k,m)=(hrho(nx1,ny,k,m)+hrho(nx,ny1,k,m))/2.
        orho(nx,1,k,m)=(orho(nx1,1,k,m)+orho(nx,2,k,m))/2.
        orho(nx,ny,k,m)=(orho(nx1,ny,k,m)+orho(nx,ny1,k,m))/2.
        !
        qpresx(nx,1,k,m)=(qpresx(nx1,1,k,m) &
        +qpresx(nx,2,k,m))/2.
        qpresx(nx,ny,k,m)=(qpresx(nx1,ny,k,m) &
        +qpresx(nx,ny1,k,m))/2.
        qpresy(nx,1,k,m)=(qpresy(nx1,1,k,m) &
            +qpresy(nx,2,k,m))/2.
        qpresy(nx,ny,k,m)=(qpresy(nx1,ny,k,m) &
            +qpresy(nx,ny1,k,m))/2.
        qpresz(nx,1,k,m)=(qpresz(nx1,1,k,m) &
            +qpresz(nx,2,k,m))/2.
        qpresz(nx,ny,k,m)=(qpresz(nx1,ny,k,m) &
            +qpresz(nx,ny1,k,m))/2.
        !
        hpresx(nx,1,k,m)=(hpresx(nx1,1,k,m) &
        +hpresx(nx,2,k,m))/2.
        hpresx(nx,ny,k,m)=(hpresx(nx1,ny,k,m) &
        +hpresx(nx,ny1,k,m))/2.
        hpresy(nx,1,k,m)=(hpresy(nx1,1,k,m) &
            +hpresy(nx,2,k,m))/2.
        hpresy(nx,ny,k,m)=(hpresy(nx1,ny,k,m) &
            +hpresy(nx,ny1,k,m))/2.
        hpresz(nx,1,k,m)=(hpresz(nx1,1,k,m) &
            +hpresz(nx,2,k,m))/2.
        hpresz(nx,ny,k,m)=(hpresz(nx1,ny,k,m) &
            +hpresz(nx,ny1,k,m))/2.
        !
        opresx(nx,1,k,m)=(opresx(nx1,1,k,m) &
        +opresx(nx,2,k,m))/2.
        opresx(nx,ny,k,m)=(opresx(nx1,ny,k,m) &
        +opresx(nx,ny1,k,m))/2.
        opresy(nx,1,k,m)=(opresy(nx1,1,k,m) &
            +opresy(nx,2,k,m))/2.
        opresy(nx,ny,k,m)=(opresy(nx1,ny,k,m) &
            +opresy(nx,ny1,k,m))/2.
        opresz(nx,1,k,m)=(opresz(nx1,1,k,m) &
            +opresz(nx,2,k,m))/2.
        opresz(nx,ny,k,m)=(opresz(nx1,ny,k,m) &
            +opresz(nx,ny1,k,m))/2.
        !
        epres(nx,1,k,m)=(epres(nx1,1,k,m) &
        +epres(nx,2,k,m))/2.
        epres(nx,ny,k,m)=(epres(nx1,ny,k,m) &
        +epres(nx,ny1,k,m))/2.
        !
        !       off-diagonal elements
        !
        qpresxy(nx,1,k,m)=(qpresxy(nx1,1,k,m) &
            +qpresxy(nx,2,k,m))/2.
        qpresxy(nx,ny,k,m)=(qpresxy(nx1,ny,k,m) &
            +qpresxy(nx,ny1,k,m))/2.
        qpresxz(nx,1,k,m)=(qpresxz(nx1,1,k,m) &
            +qpresxz(nx,2,k,m))/2.
        qpresxz(nx,ny,k,m)=(qpresxz(nx1,ny,k,m) &
            +qpresxz(nx,ny1,k,m))/2.
        qpresyz(nx,1,k,m)=(qpresyz(nx1,1,k,m) &
            +qpresyz(nx,2,k,m))/2.
        qpresyz(nx,ny,k,m)=(qpresyz(nx1,ny,k,m) &
            +qpresyz(nx,ny1,k,m))/2.
        hpresxy(nx,1,k,m)=(hpresxy(nx1,1,k,m) &
            +hpresxy(nx,2,k,m))/2.
        hpresxy(nx,ny,k,m)=(hpresxy(nx1,ny,k,m) &
            +hpresxy(nx,ny1,k,m))/2.
        hpresxz(nx,1,k,m)=(hpresxz(nx1,1,k,m) &
            +hpresxz(nx,2,k,m))/2.
        hpresxz(nx,ny,k,m)=(hpresxz(nx1,ny,k,m) &
            +hpresxz(nx,ny1,k,m))/2.
        hpresyz(nx,1,k,m)=(hpresyz(nx1,1,k,m) &
            +hpresyz(nx,2,k,m))/2.
        hpresyz(nx,ny,k,m)=(hpresyz(nx1,ny,k,m) &
            +hpresyz(nx,ny1,k,m))/2.
        opresxy(nx,1,k,m)=(opresxy(nx1,1,k,m) &
            +opresxy(nx,2,k,m))/2.
        opresxy(nx,ny,k,m)=(opresxy(nx1,ny,k,m) &
            +opresxy(nx,ny1,k,m))/2.
        opresxz(nx,1,k,m)=(opresxz(nx1,1,k,m) &
            +opresxz(nx,2,k,m))/2.
        opresxz(nx,ny,k,m)=(opresxz(nx1,ny,k,m) &
            +opresxz(nx,ny1,k,m))/2.
        opresyz(nx,1,k,m)=(opresyz(nx1,1,k,m) &
            +opresyz(nx,2,k,m))/2.
        opresyz(nx,ny,k,m)=(opresyz(nx1,ny,k,m) &
            +opresyz(nx,ny1,k,m))/2.
        !
        qpx(nx,1,k,m)=(qpx(nx1,1,k,m)+qpx(nx,2,k,m))/2.
        qpx(nx,ny,k,m)=(qpx(nx1,ny,k,m)+qpx(nx,ny1,k,m))/2.
        !
        qpy(nx,1,k,m)=(qpy(nx1,1,k,m)+qpy(nx,2,k,m))/2.
        qpy(nx,ny,k,m)=(qpy(nx1,ny,k,m)+qpy(nx,ny1,k,m))/2.
        !
        qpz(nx,1,k,m)=(qpz(nx1,1,k,m)+qpz(nx,2,k,m))/2.
        qpz(nx,ny,k,m)=(qpz(nx1,ny,k,m)+qpz(nx,ny1,k,m))/2.
        !
        hpx(nx,1,k,m)=(hpx(nx1,1,k,m)+hpx(nx,2,k,m))/2.
        hpx(nx,ny,k,m)=(hpx(nx1,ny,k,m)+hpx(nx,ny1,k,m))/2.
        !
        hpy(nx,1,k,m)=(hpy(nx1,1,k,m)+hpy(nx,2,k,m))/2.
        hpy(nx,ny,k,m)=(hpy(nx1,ny,k,m)+hpy(nx,ny1,k,m))/2.
        !
        hpz(nx,1,k,m)=(hpz(nx1,1,k,m)+hpz(nx,2,k,m))/2.
        hpz(nx,ny,k,m)=(hpz(nx1,ny,k,m)+hpz(nx,ny1,k,m))/2.
        !
        opx(nx,1,k,m)=(opx(nx1,1,k,m)+opx(nx,2,k,m))/2.
        opx(nx,ny,k,m)=(opx(nx1,ny,k,m)+opx(nx,ny1,k,m))/2.
        !
        opy(nx,1,k,m)=(opy(nx1,1,k,m)+opy(nx,2,k,m))/2.
        opy(nx,ny,k,m)=(opy(nx1,ny,k,m)+opy(nx,ny1,k,m))/2.
        !
        opz(nx,1,k,m)=(opz(nx1,1,k,m)+opz(nx,2,k,m))/2.
        opz(nx,ny,k,m)=(opz(nx1,ny,k,m)+opz(nx,ny1,k,m))/2.
        !
        bx(nx,1,k,m)=bx(nx1,1,k,m)
        bx(nx,ny,k,m)=bx(nx1,ny,k,m)
        !
        by(nx,1,k,m)=(by(nx1,1,k,m) &
        +by(nx,2,k,m))/2.
        by(nx,ny,k,m)=(by(nx1,ny,k,m) &
        +by(nx,ny1,k,m))/2.
        !
        bz(nx,1,k,m)=(bz(nx1,1,k,m) &
        +bz(nx,2,k,m))/2.
        bz(nx,ny,k,m)=(bz(nx1,ny,k,m) &
        +bz(nx,ny1,k,m))/2.
        !
    enddo
    !
    !     do corner lines at top and bottom back wall
    !
    do j=2,ny-1
        qrho(nx,j,1,m)=qrho(nx,j,2,m)
        qrho(nx,j,nz,m)=qrho(nx1,j,nz,m)
        hrho(nx,j,1,m)=hrho(nx,j,2,m)
        hrho(nx,j,nz,m)=hrho(nx1,j,nz,m)
        orho(nx,j,1,m)=orho(nx,j,2,m)
        orho(nx,j,nz,m)=orho(nx1,j,nz,m)
        !
        qpresx(nx,j,1,m)=qpresx(nx,j,2,m)
        qpresx(nx,j,nz,m)=qpresx(nx1,j,nz,m)
        qpresy(nx,j,1,m)=qpresy(nx,j,2,m)
        qpresy(nx,j,nz,m)=qpresy(nx1,j,nz,m)
        qpresz(nx,j,1,m)=qpresz(nx,j,2,m)
        qpresz(nx,j,nz,m)=qpresz(nx1,j,nz,m)
        !
        hpresx(nx,j,1,m)=hpresx(nx,j,2,m)
        hpresx(nx,j,nz,m)=hpresx(nx1,j,nz,m)
        hpresy(nx,j,1,m)=hpresy(nx,j,2,m)
        hpresy(nx,j,nz,m)=hpresy(nx1,j,nz,m)
        hpresz(nx,j,1,m)=hpresz(nx,j,2,m)
        hpresz(nx,j,nz,m)=hpresz(nx1,j,nz,m)
        !
        opresx(nx,j,1,m)=opresx(nx,j,2,m)
        opresx(nx,j,nz,m)=opresx(nx1,j,nz,m)
        opresy(nx,j,1,m)=opresy(nx,j,2,m)
        opresy(nx,j,nz,m)=opresy(nx1,j,nz,m)
        opresz(nx,j,1,m)=opresz(nx,j,2,m)
        opresz(nx,j,nz,m)=opresz(nx1,j,nz,m)
        !
        epres(nx,j,1,m)=epres(nx,j,2,m)
        !
        !       off diagnonal elements
        !
        qpresxy(nx,j,1,m)=qpresxy(nx,j,2,m)
        qpresxy(nx,j,nz,m)=qpresxy(nx1,j,nz,m)
        qpresxz(nx,j,1,m)=qpresxz(nx,j,2,m)
        qpresxz(nx,j,nz,m)=qpresxz(nx1,j,nz,m)
        qpresyz(nx,j,1,m)=qpresyz(nx,j,2,m)
        qpresyz(nx,j,nz,m)=qpresyz(nx1,j,nz,m)
        hpresxy(nx,j,1,m)=hpresxy(nx,j,2,m)
        hpresxy(nx,j,nz,m)=hpresxy(nx1,j,nz,m)
        hpresxz(nx,j,1,m)=hpresxz(nx,j,2,m)
        hpresxz(nx,j,nz,m)=hpresxz(nx1,j,nz,m)
        hpresyz(nx,j,1,m)=hpresyz(nx,j,2,m)
        hpresyz(nx,j,nz,m)=hpresyz(nx1,j,nz,m)
        opresxy(nx,j,1,m)=opresxy(nx,j,2,m)
        opresxy(nx,j,nz,m)=opresxy(nx1,j,nz,m)
        opresxz(nx,j,1,m)=opresxz(nx,j,2,m)
        opresxz(nx,j,nz,m)=opresxz(nx1,j,nz,m)
        opresyz(nx,j,1,m)=opresyz(nx,j,2,m)
        opresyz(nx,j,nz,m)=opresyz(nx1,j,nz,m)
        !
        qpx(nx,j,1,m)=qpx(nx,j,2,m)
        qpx(nx,j,nz,m)=qpx(nx1,j,nz,m)
        !
        qpy(nx,j,1,m)=qpy(nx,j,2,m)
        qpy(nx,j,nz,m)=qpy(nx1,j,nz,m)
        !
        qpz(nx,j,1,m)=qpz(nx,j,2,m)
        qpz(nx,j,nz,m)=qpz(nx1,j,nz,m)
        !
        hpx(nx,j,1,m)=hpx(nx,j,2,m)
        hpx(nx,j,nz,m)=hpx(nx1,j,nz,m)
        !
        hpy(nx,j,1,m)=hpy(nx,j,2,m)
        hpy(nx,j,nz,m)=hpy(nx1,j,nz,m)
        !
        hpz(nx,j,1,m)=hpz(nx,j,2,m)
        hpz(nx,j,nz,m)=hpz(nx1,j,nz,m)
        !
        opx(nx,j,1,m)=opx(nx,j,2,m)
        opx(nx,j,nz,m)=opx(nx1,j,nz,m)
        !
        opy(nx,j,1,m)=opy(nx,j,2,m)
        opy(nx,j,nz,m)=opy(nx1,j,nz,m)
        !
        opz(nx,j,1,m)=opz(nx,j,2,m)
        opz(nx,j,nz,m)=opz(nx1,j,nz,m)
        !
        bx(nx,j,1,m)=bx(nx,j,2,m)
        bx(nx,j,nz,m)=bx(nx1,j,nz,m)
        !
        by(nx,j,1,m)=by(nx,j,2,m)
        by(nx,j,nz,m)=by(nx1,j,nz,m)
        !
        bz(nx,j,1,m)=bz(nx,j,2,m)
        bz(nx,j,nz,m)=bz(nx1,j,nz,m)
        !
    enddo
    !
    do i=2,nx
        qrho(i,1,1,m)=qrho(i,1,2,m)
        qrho(i,ny,1,m)=qrho(i,ny,2,m)
        qrho(i,1,nz,m)=qrho(i,1,nz1,m)
        qrho(i,ny,nz,m)=qrho(i,ny,nz1,m)
        !
        hrho(i,1,1,m)=hrho(i,1,2,m)
        hrho(i,ny,1,m)=hrho(i,ny,2,m)
        hrho(i,1,nz,m)=hrho(i,1,nz1,m)
        hrho(i,ny,nz,m)=hrho(i,ny,nz1,m)
        !
        orho(i,1,1,m)=orho(i,1,2,m)
        orho(i,ny,1,m)=orho(i,ny,2,m)
        orho(i,1,nz,m)=orho(i,1,nz1,m)
        orho(i,ny,nz,m)=orho(i,ny,nz1,m)
        !
        qpresx(i,1,1,m)=qpresx(i,1,2,m)
        qpresx(i,ny,1,m)=qpresx(i,ny,2,m)
        qpresx(i,1,nz,m)=qpresx(i,1,nz1,m)
        qpresx(i,ny,nz,m)=qpresx(i,ny,nz1,m)
        qpresy(i,1,1,m)=qpresy(i,1,2,m)
        qpresy(i,ny,1,m)=qpresy(i,ny,2,m)
        qpresy(i,1,nz,m)=qpresy(i,1,nz1,m)
        qpresy(i,ny,nz,m)=qpresy(i,ny,nz1,m)
        qpresz(i,1,1,m)=qpresz(i,1,2,m)
        qpresz(i,ny,1,m)=qpresz(i,ny,2,m)
        qpresz(i,1,nz,m)=qpresz(i,1,nz1,m)
        qpresz(i,ny,nz,m)=qpresz(i,ny,nz1,m)
        !
        hpresx(i,1,1,m)=hpresx(i,1,2,m)
        hpresx(i,ny,1,m)=hpresx(i,ny,2,m)
        hpresx(i,1,nz,m)=hpresx(i,1,nz1,m)
        hpresx(i,ny,nz,m)=hpresx(i,ny,nz1,m)
        hpresy(i,1,1,m)=hpresy(i,1,2,m)
        hpresy(i,ny,1,m)=hpresy(i,ny,2,m)
        hpresy(i,1,nz,m)=hpresy(i,1,nz1,m)
        hpresy(i,ny,nz,m)=hpresy(i,ny,nz1,m)
        hpresz(i,1,1,m)=hpresz(i,1,2,m)
        hpresz(i,ny,1,m)=hpresz(i,ny,2,m)
        hpresz(i,1,nz,m)=hpresz(i,1,nz1,m)
        hpresz(i,ny,nz,m)=hpresz(i,ny,nz1,m)
        !
        opresx(i,1,1,m)=opresx(i,1,2,m)
        opresx(i,ny,1,m)=opresx(i,ny,2,m)
        opresx(i,1,nz,m)=opresx(i,1,nz1,m)
        opresx(i,ny,nz,m)=opresx(i,ny,nz1,m)
        opresy(i,1,1,m)=opresy(i,1,2,m)
        opresy(i,ny,1,m)=opresy(i,ny,2,m)
        opresy(i,1,nz,m)=opresy(i,1,nz1,m)
        opresy(i,ny,nz,m)=opresy(i,ny,nz1,m)
        opresz(i,1,1,m)=opresz(i,1,2,m)
        opresz(i,ny,1,m)=opresz(i,ny,2,m)
        opresz(i,1,nz,m)=opresz(i,1,nz1,m)
        opresz(i,ny,nz,m)=opresz(i,ny,nz1,m)
        !
        epres(i,1,1,m)=epres(i,1,2,m)
        epres(i,ny,1,m)=epres(i,ny,2,m)
        epres(i,1,nz,m)=epres(i,1,nz1,m)
        epres(i,ny,nz,m)=epres(i,ny,nz1,m)
        !
        !       off-diagonal elements
        !
        qpresxy(i,1,1,m)=qpresxy(i,1,2,m)
        qpresxy(i,ny,1,m)=qpresxy(i,ny,2,m)
        qpresxy(i,1,nz,m)=qpresxy(i,1,nz1,m)
        qpresxy(i,ny,nz,m)=qpresxy(i,ny,nz1,m)
        qpresxz(i,1,1,m)=qpresxz(i,1,2,m)
        qpresxz(i,ny,1,m)=qpresxz(i,ny,2,m)
        qpresxz(i,1,nz,m)=qpresxz(i,1,nz1,m)
        qpresxz(i,ny,nz,m)=qpresxz(i,ny,nz1,m)
        qpresyz(i,1,1,m)=qpresyz(i,1,2,m)
        qpresyz(i,ny,1,m)=qpresyz(i,ny,2,m)
        qpresyz(i,1,nz,m)=qpresyz(i,1,nz1,m)
        qpresyz(i,ny,nz,m)=qpresyz(i,ny,nz1,m)
        !
        hpresxy(i,1,1,m)=hpresxy(i,1,2,m)
        hpresxy(i,ny,1,m)=hpresxy(i,ny,2,m)
        hpresxy(i,1,nz,m)=hpresxy(i,1,nz1,m)
        hpresxy(i,ny,nz,m)=hpresxy(i,ny,nz1,m)
        hpresxz(i,1,1,m)=hpresxz(i,1,2,m)
        hpresxz(i,ny,1,m)=hpresxz(i,ny,2,m)
        hpresxz(i,1,nz,m)=hpresxz(i,1,nz1,m)
        hpresxz(i,ny,nz,m)=hpresxz(i,ny,nz1,m)
        hpresyz(i,1,1,m)=hpresyz(i,1,2,m)
        hpresyz(i,ny,1,m)=hpresyz(i,ny,2,m)
        hpresyz(i,1,nz,m)=hpresyz(i,1,nz1,m)
        hpresyz(i,ny,nz,m)=hpresyz(i,ny,nz1,m)
        !
        opresxy(i,1,1,m)=opresxy(i,1,2,m)
        opresxy(i,ny,1,m)=opresxy(i,ny,2,m)
        opresxy(i,1,nz,m)=opresxy(i,1,nz1,m)
        opresxy(i,ny,nz,m)=opresxy(i,ny,nz1,m)
        opresxz(i,1,1,m)=opresxz(i,1,2,m)
        opresxz(i,ny,1,m)=opresxz(i,ny,2,m)
        opresxz(i,1,nz,m)=opresxz(i,1,nz1,m)
        opresxz(i,ny,nz,m)=opresxz(i,ny,nz1,m)
        opresyz(i,1,1,m)=opresyz(i,1,2,m)
        opresyz(i,ny,1,m)=opresyz(i,ny,2,m)
        opresyz(i,1,nz,m)=opresyz(i,1,nz1,m)
        opresyz(i,ny,nz,m)=opresyz(i,ny,nz1,m)
        !
        qpx(i,1,1,m)=qpx(i,1,2,m)
        qpx(i,ny,1,m)=qpx(i,ny,2,m)
        qpx(i,1,nz,m)=qpx(i,1,nz1,m)
        qpx(i,ny,nz,m)=qpx(i,ny,nz1,m)
        !
        qpy(i,1,1,m)=qpy(i,1,2,m)
        qpy(i,ny,1,m)=qpy(i,ny,2,m)
        qpy(i,1,nz,m)=qpy(i,1,nz1,m)
        qpy(i,ny,nz,m)=qpy(i,ny,nz1,m)
        !
        qpz(i,1,1,m)=qpz(i,1,2,m)
        qpz(i,ny,1,m)=qpz(i,ny,2,m)
        qpz(i,1,nz,m)=qpz(i,1,nz1,m)
        qpz(i,ny,nz,m)=qpz(i,ny,nz1,m)
        !
        hpx(i,1,1,m)=hpx(i,1,2,m)
        hpx(i,ny,1,m)=hpx(i,ny,2,m)
        hpx(i,1,nz,m)=hpx(i,1,nz1,m)
        hpx(i,ny,nz,m)=hpx(i,ny,nz1,m)
        !
        hpy(i,1,1,m)=hpy(i,1,2,m)
        hpy(i,ny,1,m)=hpy(i,ny,2,m)
        hpy(i,1,nz,m)=hpy(i,1,nz1,m)
        hpy(i,ny,nz,m)=hpy(i,ny,nz1,m)
        !
        hpz(i,1,1,m)=hpz(i,1,2,m)
        hpz(i,ny,1,m)=hpz(i,ny,2,m)
        hpz(i,1,nz,m)=hpz(i,1,nz1,m)
        hpz(i,ny,nz,m)=hpz(i,ny,nz1,m)
        !
        opx(i,1,1,m)=opx(i,1,2,m)
        opx(i,ny,1,m)=opx(i,ny,2,m)
        opx(i,1,nz,m)=opx(i,1,nz1,m)
        opx(i,ny,nz,m)=opx(i,ny,nz1,m)
        !
        opy(i,1,1,m)=opy(i,1,2,m)
        opy(i,ny,1,m)=opy(i,ny,2,m)
        opy(i,1,nz,m)=opy(i,1,nz1,m)
        opy(i,ny,nz,m)=opy(i,ny,nz1,m)
        !
        opz(i,1,1,m)=opz(i,1,2,m)
        opz(i,ny,1,m)=opz(i,ny,2,m)
        opz(i,1,nz,m)=opz(i,1,nz1,m)
        opz(i,ny,nz,m)=opz(i,ny,nz1,m)
        !
        bx(i,1,1,m)=bx(i,1,2,m)
        bx(i,ny,1,m)=bx(i,ny,2,m)
        bx(i,1,nz,m)=bx(i,1,nz1,m)
        bx(i,ny,nz,m)=bx(i,ny,nz1,m)
        !
        by(i,1,1,m)=by(i,1,2,m)
        by(i,ny,1,m)=by(i,ny,2,m)
        by(i,1,nz,m)=by(i,1,nz1,m)
        by(i,ny,nz,m)=by(i,ny,nz1,m)
        !
        bz(i,1,1,m)=bz(i,1,2,m)
        bz(i,ny,1,m)=bz(i,ny,2,m)
        bz(i,1,nz,m)=bz(i,1,nz1,m)
        bz(i,ny,nz,m)=bz(i,ny,nz1,m)
        !
    enddo
    !
    return
end
