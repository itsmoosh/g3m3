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
    nx,ny,nz,n_grids,box, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     this routine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    integer box
    dimension &
    bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    qrho(nx,ny,nz,n_grids),qpx(nx,ny,nz,n_grids), &
    qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    qpresx(nx,ny,nz,n_grids),qpresy(nx,ny,nz,n_grids), &
    qpresz(nx,ny,nz,n_grids), qpresxy(nx,ny,nz,n_grids), &
    qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
    !
    hrho(nx,ny,nz,n_grids),hpx(nx,ny,nz,n_grids), &
    hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    hpresx(nx,ny,nz,n_grids),hpresy(nx,ny,nz,n_grids), &
    hpresz(nx,ny,nz,n_grids), hpresxy(nx,ny,nz,n_grids), &
    hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
    !
    orho(nx,ny,nz,n_grids),opx(nx,ny,nz,n_grids), &
    opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    opresx(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
    opresy(nx,ny,nz,n_grids), opresxy(nx,ny,nz,n_grids), &
    opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
    !
    epres(nx,ny,nz,n_grids)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    call corer(qrho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresxy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresxz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpresyz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(qpz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call corer(hrho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresxy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresxz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpresyz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(hpz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call corer(orho,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresxy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresxz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opresyz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opy,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(opz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(epres,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(bx,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(by,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call corer(bz,nx,ny,nz,n_grids,box,grd_xmin,grd_xmax, &
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
    nx,ny,nz,n_grids,box,t_old,t_new,t, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !     this subroutine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    integer box
    dimension qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
    qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
    qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
    hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
    hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
    opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
    opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    epres(nx,ny,nz,n_grids)
    dimension wrkqrho(nx,ny,nz,n_grids),wrkqpresx(nx,ny,nz,n_grids), &
    wrkqpresy(nx,ny,nz,n_grids),wrkqpresz(nx,ny,nz,n_grids), &
    wrkqpx(nx,ny,nz,n_grids),wrkqpy(nx,ny,nz,n_grids), &
    wrkqpz(nx,ny,nz,n_grids), &
    wrkhrho(nx,ny,nz,n_grids),wrkhpresx(nx,ny,nz,n_grids), &
    wrkhpresy(nx,ny,nz,n_grids),wrkhpresz(nx,ny,nz,n_grids), &
    wrkhpx(nx,ny,nz,n_grids),wrkhpy(nx,ny,nz,n_grids), &
    wrkhpz(nx,ny,nz,n_grids), &
    wrkorho(nx,ny,nz,n_grids),wrkopresx(nx,ny,nz,n_grids), &
    wrkopresy(nx,ny,nz,n_grids),wrkopresz(nx,ny,nz,n_grids), &
    wrkopx(nx,ny,nz,n_grids),wrkopy(nx,ny,nz,n_grids), &
    wrkopz(nx,ny,nz,n_grids), &
    wrkbx(nx,ny,nz,n_grids),wrkby(nx,ny,nz,n_grids), &
    wrkbz(nx,ny,nz,n_grids), &
    wrkepres(nx,ny,nz,n_grids)
    dimension oldqrho(nx,ny,nz,n_grids),oldqpresx(nx,ny,nz,n_grids), &
    oldqpresy(nx,ny,nz,n_grids),oldqpresz(nx,ny,nz,n_grids), &
    oldqpx(nx,ny,nz,n_grids),oldqpy(nx,ny,nz,n_grids), &
    oldqpz(nx,ny,nz,n_grids), &
    oldhrho(nx,ny,nz,n_grids),oldhpresx(nx,ny,nz,n_grids), &
    oldhpresy(nx,ny,nz,n_grids),oldhpresz(nx,ny,nz,n_grids), &
    oldhpx(nx,ny,nz,n_grids),oldhpy(nx,ny,nz,n_grids), &
    oldhpz(nx,ny,nz,n_grids), &
    oldorho(nx,ny,nz,n_grids),oldopresx(nx,ny,nz,n_grids), &
    oldopresy(nx,ny,nz,n_grids),oldopresz(nx,ny,nz,n_grids), &
    oldopx(nx,ny,nz,n_grids),oldopy(nx,ny,nz,n_grids), &
    oldopz(nx,ny,nz,n_grids), &
    oldbx(nx,ny,nz,n_grids),oldby(nx,ny,nz,n_grids), &
    oldbz(nx,ny,nz,n_grids), &
    oldepres(nx,ny,nz,n_grids)
    !
    dimension qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
    qpresyz(nx,ny,nz,n_grids),hpresxy(nx,ny,nz,n_grids), &
    hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
    opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
    opresyz(nx,ny,nz,n_grids)
    !
    dimension oldqpresxy(nx,ny,nz,n_grids),oldqpresxz(nx,ny,nz,n_grids), &
    oldqpresyz(nx,ny,nz,n_grids),oldhpresxy(nx,ny,nz,n_grids), &
    oldhpresxz(nx,ny,nz,n_grids),oldhpresyz(nx,ny,nz,n_grids), &
    oldopresxy(nx,ny,nz,n_grids),oldopresxz(nx,ny,nz,n_grids), &
    oldopresyz(nx,ny,nz,n_grids)
    dimension wrkqpresxy(nx,ny,nz,n_grids),wrkqpresxz(nx,ny,nz,n_grids), &
    wrkqpresyz(nx,ny,nz,n_grids),wrkhpresxy(nx,ny,nz,n_grids), &
    wrkhpresxz(nx,ny,nz,n_grids),wrkhpresyz(nx,ny,nz,n_grids), &
    wrkopresxy(nx,ny,nz,n_grids),wrkopresxz(nx,ny,nz,n_grids), &
    wrkopresyz(nx,ny,nz,n_grids)
    !
    dimension work(nx,ny,nz)
    !
    dimension t_old(n_grids),t_new(n_grids)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    call flanks(wrkqrho,qrho,oldqrho,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresx,qpresx,oldqpresx,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresy,qpresy,oldqpresy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresz,qpresz,oldqpresz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresxy,qpresxy,oldqpresxy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresxz,qpresxz,oldqpresxz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpresyz,qpresyz,oldqpresyz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpx,qpx,oldqpx,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpy,qpy,oldqpy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkqpz,qpz,oldqpz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkhrho,hrho,oldhrho,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresx,hpresx,oldhpresx,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresy,hpresy,oldhpresy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresz,hpresz,oldhpresz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresxy,hpresxy,oldhpresxy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresxz,hpresxz,oldhpresxz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpresyz,hpresyz,oldhpresyz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpx,hpx,oldhpx,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpy,hpy,oldhpy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkhpz,hpz,oldhpz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkorho,orho,oldorho,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresx,opresx,oldopresx,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresy,opresy,oldopresy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresz,opresz,oldopresz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresxy,opresxy,oldopresxy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresxz,opresxz,oldopresxz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopresyz,opresyz,oldopresyz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopx,opx,oldopx,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopy,opy,oldopy,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkopz,opz,oldopz,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    call flanks(wrkepres,epres,oldepres,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkbx,bx,oldbx,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkby,by,oldby,nx,ny,nz,n_grids,box, &
    t_new,t_old,t,work,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    call flanks(wrkbz,bz,oldbz,nx,ny,nz,n_grids,box, &
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
    nx,ny,nz,n_grids, &
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
    integer box
    dimension &
    bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    qrho(nx,ny,nz,n_grids),qpx(nx,ny,nz,n_grids), &
    qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    qpresx(nx,ny,nz,n_grids),qpresy(nx,ny,nz,n_grids), &
    qpresz(nx,ny,nz,n_grids), qpresxy(nx,ny,nz,n_grids), &
    qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
    !
    hrho(nx,ny,nz,n_grids),hpx(nx,ny,nz,n_grids), &
    hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    hpresx(nx,ny,nz,n_grids),hpresy(nx,ny,nz,n_grids), &
    hpresz(nx,ny,nz,n_grids), hpresxy(nx,ny,nz,n_grids), &
    hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
    ! 
    orho(nx,ny,nz,n_grids),opx(nx,ny,nz,n_grids), &
    opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    opresx(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
    opresy(nx,ny,nz,n_grids), opresxy(nx,ny,nz,n_grids), &
    opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
    !
    epres(nx,ny,nz,n_grids)
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
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    call corer_grds(qrho,nx,ny,nz,n_grids,main_ngrd, &
    qrho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresx,nx,ny,nz,n_grids,main_ngrd, &
    qpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresy,nx,ny,nz,n_grids,main_ngrd, &
    qpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresz,nx,ny,nz,n_grids,main_ngrd, &
    qpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresxy,nx,ny,nz,n_grids,main_ngrd, &
    qpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresxz,nx,ny,nz,n_grids,main_ngrd, &
    qpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpresyz,nx,ny,nz,n_grids,main_ngrd, &
    qpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpx,nx,ny,nz,n_grids,main_ngrd, &
    qpx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpy,nx,ny,nz,n_grids,main_ngrd, &
    qpy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(qpz,nx,ny,nz,n_grids,main_ngrd, &
    qpz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(hrho,nx,ny,nz,n_grids,main_ngrd, &
    hrho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresx,nx,ny,nz,n_grids,main_ngrd, &
    hpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresy,nx,ny,nz,n_grids,main_ngrd, &
    hpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresz,nx,ny,nz,n_grids,main_ngrd, &
    hpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresxy,nx,ny,nz,n_grids,main_ngrd, &
    hpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresxz,nx,ny,nz,n_grids,main_ngrd, &
    hpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpresyz,nx,ny,nz,n_grids,main_ngrd, &
    hpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpx,nx,ny,nz,n_grids,main_ngrd, &
    hpx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpy,nx,ny,nz,n_grids,main_ngrd, &
    hpy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(hpz,nx,ny,nz,n_grids,main_ngrd, &
    hpz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(orho,nx,ny,nz,n_grids,main_ngrd, &
    orho_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresx,nx,ny,nz,n_grids,main_ngrd, &
    opresx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresy,nx,ny,nz,n_grids,main_ngrd, &
    opresy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresz,nx,ny,nz,n_grids,main_ngrd, &
    opresz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresxy,nx,ny,nz,n_grids,main_ngrd, &
    opresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresxz,nx,ny,nz,n_grids,main_ngrd, &
    opresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opresyz,nx,ny,nz,n_grids,main_ngrd, &
    opresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opx,nx,ny,nz,n_grids,main_ngrd, &
    opx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opy,nx,ny,nz,n_grids,main_ngrd, &
    opy_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(opz,nx,ny,nz,n_grids,main_ngrd, &
    opz_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    !
    call corer_grds(epres,nx,ny,nz,n_grids,main_ngrd, &
    epres_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(bx,nx,ny,nz,n_grids,main_ngrd, &
    bx_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(by,nx,ny,nz,n_grids,main_ngrd, &
    by_n,nx_n,ny_n,nz_n,ngrd_n, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
    grd_zmin_n,grd_zmax_n)
    call corer_grds(bz,nx,ny,nz,n_grids,main_ngrd, &
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
    nx,ny,nz,n_grids,t_old,t_new, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !     this subroutine applies boundary conditions to all grid types
    !     of the system and at any irregular boundaries
    !
    integer box
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
    dimension qrho(nx,ny,nz,n_grids), &
    qpresx(nx,ny,nz,n_grids),qpresy(nx,ny,nz,n_grids), &
    qpresz(nx,ny,nz,n_grids),qpresxy(nx,ny,nz,n_grids), &
    qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
    qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    !
    hrho(nx,ny,nz,n_grids), &
    hpresx(nx,ny,nz,n_grids),hpresy(nx,ny,nz,n_grids), &
    hpresz(nx,ny,nz,n_grids),hpresxy(nx,ny,nz,n_grids), &
    hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
    hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    !
    orho(nx,ny,nz,n_grids), &
    opresx(nx,ny,nz,n_grids),opresy(nx,ny,nz,n_grids), &
    opresz(nx,ny,nz,n_grids),opresxy(nx,ny,nz,n_grids), &
    opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
    opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    !
    bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    epres(nx,ny,nz,n_grids)
    !
    dimension oldqrho(nx,ny,nz,n_grids), &
    oldqpresx(nx,ny,nz,n_grids),oldqpresy(nx,ny,nz,n_grids), &
    oldqpresz(nx,ny,nz,n_grids),oldqpresxy(nx,ny,nz,n_grids), &
    oldqpresxz(nx,ny,nz,n_grids),oldqpresyz(nx,ny,nz,n_grids), &
    oldqpx(nx,ny,nz,n_grids),oldqpy(nx,ny,nz,n_grids), &
    oldqpz(nx,ny,nz,n_grids), &
    !
    oldhrho(nx,ny,nz,n_grids), &
    oldhpresx(nx,ny,nz,n_grids),oldhpresy(nx,ny,nz,n_grids), &
    oldhpresz(nx,ny,nz,n_grids),oldhpresxy(nx,ny,nz,n_grids), &
    oldhpresxz(nx,ny,nz,n_grids),oldhpresyz(nx,ny,nz,n_grids), &
    oldhpx(nx,ny,nz,n_grids),oldhpy(nx,ny,nz,n_grids), &
    oldhpz(nx,ny,nz,n_grids), &
    !
    oldorho(nx,ny,nz,n_grids), &
    oldopresx(nx,ny,nz,n_grids),oldopresy(nx,ny,nz,n_grids), &
    oldopresz(nx,ny,nz,n_grids),oldopresxy(nx,ny,nz,n_grids), &
    oldopresxz(nx,ny,nz,n_grids),oldopresyz(nx,ny,nz,n_grids), &
    oldopx(nx,ny,nz,n_grids),oldopy(nx,ny,nz,n_grids), &
    oldopz(nx,ny,nz,n_grids), &
    oldbx(nx,ny,nz,n_grids),oldby(nx,ny,nz,n_grids), &
    oldbz(nx,ny,nz,n_grids),oldepres(nx,ny,nz,n_grids)
    !
    dimension work(nx,ny,nz)
    dimension t_old(n_grids),t_new(n_grids)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n)
    !
    call flanks_grds(qrho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qrho,oldqrho,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresx,oldqpresx,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresy,oldqpresy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresz,oldqpresz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresxy,oldqpresxy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresxz,oldqpresxz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpresyz,oldqpresyz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpx,oldqpx,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpy,oldqpy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(qpz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    qpz,oldqpz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(hrho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hrho,oldhrho,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresx,oldhpresx,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresy,oldhpresy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresz,oldhpresz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresxy,oldhpresxy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresxz,oldhpresxz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpresyz,oldhpresyz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpx,oldhpx,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpy,oldhpy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(hpz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    hpz,oldhpz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(orho_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    orho,oldorho,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresx,oldopresx,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresy,oldopresy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresz,oldopresz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresxy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresxy,oldopresxy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresxz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresxz,oldopresxz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opresyz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opresyz,oldopresyz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opx,oldopx,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opy_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opy,oldopy,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(opz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    opz,oldopz,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    !
    call flanks_grds(epres_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    epres,oldepres,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(bx_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    bx,oldbx,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(by_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    by,oldby,work,nx,ny,nz,n_grids, &
    t_new,t_old,t,grd_xmin,grd_xmax, &
    grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
    grd_xmin_n,grd_xmax_n,grd_ymin_n, &
    grd_ymax_n,grd_zmin_n,grd_zmax_n)
    call flanks_grds(bz_n,nx_n,ny_n,nz_n,ngrd_n,main_ngrd_n, &
    bz,oldbz,work,nx,ny,nz,n_grids, &
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
    nx,ny,nz,n_grids,parm_srf,parm_mid,parm_zero, &
    ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
    mbndry,msrf,mmid,mzero, &
    erho,epress,alpha_e,ti_te,o_conc, &
    re_equiv,r_inner,sbx_wind,spress, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !     this subroutine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    integer box
    dimension &
    qrho(nx,ny,nz,n_grids),qpx(nx,ny,nz,n_grids), &
    qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    qpresx(nx,ny,nz,n_grids),qpresy(nx,ny,nz,n_grids), &
    qpresz(nx,ny,nz,n_grids), qpresxy(nx,ny,nz,n_grids), &
    qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
    !
    hrho(nx,ny,nz,n_grids),hpx(nx,ny,nz,n_grids), &
    hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    hpresx(nx,ny,nz,n_grids),hpresy(nx,ny,nz,n_grids), &
    hpresz(nx,ny,nz,n_grids), hpresxy(nx,ny,nz,n_grids), &
    hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
    ! 
    orho(nx,ny,nz,n_grids),opx(nx,ny,nz,n_grids), &
    opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    opresx(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
    opresy(nx,ny,nz,n_grids), opresxy(nx,ny,nz,n_grids), &
    opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
    !
    epres(nx,ny,nz,n_grids), &
    !
    bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids)
    !
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real  parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
	integer fluxs_f
	integer speed_f
	integer concs_f
	integer grdpt_f
	integer recdt_f
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	!
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    !     set surface conditions around earth
    !      withinterior temperature as constant
    !
    !     write(grdpt_f,*)'inside bndry inner with',mbndry
    !     write(grdpt_f,*)  numsrf,nummid,numzero
    !
    aheight=r_inner+2.
    !
    box=1
    dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    r_rot2=r_rot**2
    !
    !     set ionospheric boundary points - plasma quantities defined
    !
    do n=1,numsrf(box)
        !
        !        get coords of point
        !
        i=ijsrf(box,1,n)
        j=ijsrf(box,2,n)
        k=ijsrf(box,3,n)
        !
        !         set rotational speeds
        !
        ay=grd_ymin(box)+dy*(j-1)
        ax=grd_xmin(box)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        !        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        !        qden=parm_srf(1,n)
        qden=amax1(parm_srf(box,1,n),qrho(i,j,k,box))
        hden=amax1(parm_srf(box,2,n),hrho(i,j,k,box))
        oden=amax1(parm_srf(box,3,n),orho(i,j,k,box))
        !
        qrho(i,j,k,box)=qden
        hrho(i,j,k,box)=hden
        orho(i,j,k,box)=oden
        !
        qpresx(i,j,k,box)=parm_srf(box,4,n)
        hpresx(i,j,k,box)=parm_srf(box,5,n)
        opresx(i,j,k,box)=parm_srf(box,6,n)
        !
        qpresy(i,j,k,box)=parm_srf(box,4,n)
        hpresy(i,j,k,box)=parm_srf(box,5,n)
        opresy(i,j,k,box)=parm_srf(box,6,n)
        !
        qpresz(i,j,k,box)=parm_srf(box,4,n)
        hpresz(i,j,k,box)=parm_srf(box,5,n)
        opresz(i,j,k,box)=parm_srf(box,6,n)
        !
        epres(i,j,k,box)=parm_srf(box,7,n)
        !
        hpx(i,j,k,box)=hden*rvx
        hpy(i,j,k,box)=hden*rvy
        hpz(i,j,k,box)=hden*rvz
        opx(i,j,k,box)=oden*rvx
        opy(i,j,k,box)=oden*rvy
        opz(i,j,k,box)=oden*rvz
        qpx(i,j,k,box)=qden*rvx
        qpy(i,j,k,box)=qden*rvy
        qpz(i,j,k,box)=qden*rvz
        !
        !        iostropic conditions: set off axis pressures to zero
        !
        qpresxy(i,j,k,box)=0.
        qpresxz(i,j,k,box)=0.
        qpresyz(i,j,k,box)=0.
        hpresxy(i,j,k,box)=0.
        hpresxz(i,j,k,box)=0.
        hpresyz(i,j,k,box)=0.
        opresxy(i,j,k,box)=0.
        opresxz(i,j,k,box)=0.
        opresyz(i,j,k,box)=0.
        !
    enddo
    !
    !     set atmospheric boundary points - specify magnetic field profile
    !
    do n=1,nummid(box)
        !
        !        get coords of point
        !
        i=ijmid(box,1,n)
        j=ijmid(box,2,n)
        k=ijmid(box,3,n)
        !
        qrho(i,j,k,box)=parm_mid(box,1,n)
        hrho(i,j,k,box)=parm_mid(box,2,n)
        orho(i,j,k,box)=parm_mid(box,3,n)
        !
        qpresx(i,j,k,box)=parm_mid(box,4,n)
        hpresx(i,j,k,box)=parm_mid(box,5,n)
        opresx(i,j,k,box)=parm_mid(box,6,n)
        !
        qpresy(i,j,k,box)=parm_mid(box,4,n)
        hpresy(i,j,k,box)=parm_mid(box,5,n)
        opresy(i,j,k,box)=parm_mid(box,6,n)
        !
        qpresz(i,j,k,box)=parm_mid(box,4,n)
        hpresz(i,j,k,box)=parm_mid(box,5,n)
        opresz(i,j,k,box)=parm_mid(box,6,n)
        !
        epres(i,j,k,box)=parm_mid(box,7,n)
        !
        ay=grd_ymin(box)+dy*(j-1)
        ax=grd_xmin(box)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        opx(i,j,k,box)=orho(i,j,k,box)*rvx
        opy(i,j,k,box)=orho(i,j,k,box)*rvy
        opz(i,j,k,box)=orho(i,j,k,box)*rvz
        hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
        hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
        hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
        qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
        qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
        qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
        !
        qpresxy(i,j,k,box)=0.
        qpresxz(i,j,k,box)=0.
        qpresyz(i,j,k,box)=0.
        hpresxy(i,j,k,box)=0.
        hpresxz(i,j,k,box)=0.
        hpresyz(i,j,k,box)=0.
        opresxy(i,j,k,box)=0.
        opresxz(i,j,k,box)=0.
        opresyz(i,j,k,box)=0.
        !
        !        bx(i,j,k,box)=sbx_wind
        !        by(i,j,k,box)=0.
        !        bz(i,j,k,box)=0.
        !
    enddo
    !
    !      reset interior points
    !
    do n=1,numzero(box)
        i=ijzero(box,1,n)
        j=ijzero(box,2,n)
        k=ijzero(box,3,n)
        !
        qrho(i,j,k,box)=parm_zero(box,1,n)
        hrho(i,j,k,box)=parm_zero(box,2,n)
        orho(i,j,k,box)=parm_zero(box,3,n)
        !
        qpresx(i,j,k,box)=parm_zero(box,4,n)
        hpresx(i,j,k,box)=parm_zero(box,5,n)
        opresx(i,j,k,box)=parm_zero(box,6,n)
        !
        qpresy(i,j,k,box)=parm_zero(box,4,n)
        hpresy(i,j,k,box)=parm_zero(box,5,n)
        opresy(i,j,k,box)=parm_zero(box,6,n)
        !
        qpresz(i,j,k,box)=parm_zero(box,4,n)
        hpresz(i,j,k,box)=parm_zero(box,5,n)
        opresz(i,j,k,box)=parm_zero(box,6,n)
        !
        epres(i,j,k,box)=parm_zero(box,7,n)
        !
        ay=grd_ymin(box)+dy*(j-1)
        ax=grd_xmin(box)+dx*(i-1)
        rx=ax*re_equiv
        ry=ay*re_equiv
        vfrac=(r_rot2)/(r_rot2+ rx**2+ry**2)
        vfrac=1.
        rvy=rx*v_rot*vfrac
        rvx=-ry*v_rot*vfrac
        rvz=0.
        !
        opx(i,j,k,box)=orho(i,j,k,box)*rvx
        opy(i,j,k,box)=orho(i,j,k,box)*rvy
        opz(i,j,k,box)=orho(i,j,k,box)*rvz
        hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
        hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
        hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
        qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
        qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
        qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
        !
        qpresxy(i,j,k,box)=0.
        qpresxz(i,j,k,box)=0.
        qpresyz(i,j,k,box)=0.
        hpresxy(i,j,k,box)=0.
        hpresxz(i,j,k,box)=0.
        hpresyz(i,j,k,box)=0.
        opresxy(i,j,k,box)=0.
        opresxz(i,j,k,box)=0.
        opresyz(i,j,k,box)=0.
        !
        bx(i,j,k,box)=0.
        by(i,j,k,box)=0.
        bz(i,j,k,box)=0.
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
    box,nx,ny,nz,n_grids,parm_srf,parm_mid,parm_zero, &
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
    integer box
    dimension qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
    qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
    qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
    qpresyz(nx,ny,nz,n_grids), &
    qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
    hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
    hpresxy(nx,ny,nz,n_grids),hpresxz(nx,ny,nz,n_grids), &
    hpresyz(nx,ny,nz,n_grids), &
    hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
    opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
    opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
    opresyz(nx,ny,nz,n_grids), &
    opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    epres(nx,ny,nz,n_grids)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
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
    dx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
    dy=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
    dz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
    !
    !     set atmospheric boundary points -- specify magnetic field profile
    !
    do n=1,nummid(box)
        !
        !        get coords of point
        !
        i=ijmid(box,1,n)
        j=ijmid(box,2,n)
        k=ijmid(box,3,n)
        !
        qrho(i,j,k,box)=parm_mid(box,1,n)
        hrho(i,j,k,box)=parm_mid(box,2,n)
        orho(i,j,k,box)=parm_mid(box,3,n)
        !
        qpresx(i,j,k,box)=parm_mid(box,4,n)
        qpresy(i,j,k,box)=qpresx(i,j,k,box)
        qpresz(i,j,k,box)=qpresx(i,j,k,box)*0.9

        qpresxy(i,j,k,box)=0.
        qpresxz(i,j,k,box)=0.
        qpresyz(i,j,k,box)=0.
        !
        hpresx(i,j,k,box)=parm_mid(box,5,n)
        hpresy(i,j,k,box)=hpresx(i,j,k,box)
        hpresz(i,j,k,box)=hpresx(i,j,k,box)*0.9
		!
        hpresxy(i,j,k,box)=0.
        hpresxz(i,j,k,box)=0.
        hpresyz(i,j,k,box)=0.
        !
        opresx(i,j,k,box)=parm_mid(box,6,n)
        opresy(i,j,k,box)=opresx(i,j,k,box)
        opresz(i,j,k,box)=opresx(i,j,k,box)*0.9
		!
        opresxy(i,j,k,box)=0.
        opresxz(i,j,k,box)=0.
        opresyz(i,j,k,box)=0.
        !
        epres(i,j,k,box)=parm_mid(box,7,n)
        !
        !        qrho(i,j,k,box)=amax1(0.25*qrho(i,j,k,box),parm_mid(box,1,n))
        !        hrho(i,j,k,box)=amax1(0.25*hrho(i,j,k,box),parm_mid(box,2,n))
        !        orho(i,j,k,box)=amax1(0.25*orho(i,j,k,box),parm_mid(box,3,n))
        !
        !        qpres(i,j,k,box)=amax1(0.25*qpres(i,j,k,box),parm_mid(box,4,n))
        !        hpres(i,j,k,box)=amax1(0.25*hpres(i,j,k,box),parm_mid(box,5,n))
        !        opres(i,j,k,box)=amax1(0.25*opres(i,j,k,box),parm_mid(box,6,n))
        !        epres(i,j,k,box)=amax1(0.25*epres(i,j,k,box),parm_mid(box,7,n))
        !
        opx(i,j,k,box)=orho(i,j,k,box)*rvx
        opy(i,j,k,box)=orho(i,j,k,box)*rvy
        opz(i,j,k,box)=orho(i,j,k,box)*rvz
        hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
        hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
        hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
        qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
        qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
        qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
        !
        !        bx(i,j,k,box)=sbx_wind
        !        by(i,j,k,box)=0.
        !        bz(i,j,k,box)=0.
        !
    enddo
    !
    !      reset interior points
    !
    do n=1,numzero(box)
        i=ijzero(box,1,n)
        j=ijzero(box,2,n)
        k=ijzero(box,3,n)
        !
        qrho(i,j,k,box)=parm_zero(box,1,n)
        hrho(i,j,k,box)=parm_zero(box,2,n)
        orho(i,j,k,box)=parm_zero(box,3,n)
        !
        qpresx(i,j,k,box)=parm_zero(box,4,n)
        qpresy(i,j,k,box)=qpresx(i,j,k,box)
        qpresz(i,j,k,box)=qpresx(i,j,k,box)*0.9
		!
        qpresxy(i,j,k,box)=0.
        qpresxz(i,j,k,box)=0.
        qpresyz(i,j,k,box)=0.
        !
        hpresx(i,j,k,box)=parm_zero(box,5,n)
        hpresy(i,j,k,box)=hpresx(i,j,k,box)
        hpresz(i,j,k,box)=hpresx(i,j,k,box)*0.9

        hpresxy(i,j,k,box)=0.
        hpresxz(i,j,k,box)=0.
        hpresyz(i,j,k,box)=0.
        !
        opresx(i,j,k,box)=parm_zero(box,6,n)
        opresy(i,j,k,box)=opresx(i,j,k,box)
        opresz(i,j,k,box)=opresx(i,j,k,box)*0.9
		!
        opresxy(i,j,k,box)=0.
        opresxz(i,j,k,box)=0.
        opresyz(i,j,k,box)=0.
        !
        epres(i,j,k,box)=parm_zero(box,7,n)
        !
        !        qrho(i,j,k,box)=amax1(0.10*qrho(i,j,k,box),parm_zero(box,1,n))
        !        hrho(i,j,k,box)=amax1(0.10*hrho(i,j,k,box),parm_zero(box,2,n))
        !        orho(i,j,k,box)=amax1(0.10*orho(i,j,k,box),parm_zero(box,3,n))
        !
        !        qpres(i,j,k,box)=amax1(0.05*qpres(i,j,k,box),parm_zero(box,4,n))
        !        hpres(i,j,k,box)=amax1(0.05*hpres(i,j,k,box),parm_zero(box,5,n))
        !        opres(i,j,k,box)=amax1(0.05*opres(i,j,k,box),parm_zero(box,6,n))
        !        epres(i,j,k,box)=amax1(0.05*epres(i,j,k,box),parm_zero(box,7,n))
        !
        opx(i,j,k,box)=orho(i,j,k,box)*rvx
        opy(i,j,k,box)=orho(i,j,k,box)*rvy
        opz(i,j,k,box)=orho(i,j,k,box)*rvz
        hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
        hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
        hpz(i,j,k,box)=hrho(i,j,k,box)*rvz
        qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
        qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
        qpz(i,j,k,box)=qrho(i,j,k,box)*rvz
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
    bx,by,bz,bx0,by0,bz0,nx,ny,nz,n_grids, &
    srho,rho_frac,o_conc,spress,spx,spy,spz, &
    sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)
    !
    !     this subroutine applies boundary conditions around the edges
    !     of the system and at any irregular boundaries
    !
    integer box
    dimension &
    qrho(nx,ny,nz,n_grids),qpx(nx,ny,nz,n_grids), &
    qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
    qpresx(nx,ny,nz,n_grids),qpresy(nx,ny,nz,n_grids), &
    qpresz(nx,ny,nz,n_grids), qpresxy(nx,ny,nz,n_grids), &
    qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
    !
    hrho(nx,ny,nz,n_grids),hpx(nx,ny,nz,n_grids), &
    hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
    hpresx(nx,ny,nz,n_grids),hpresy(nx,ny,nz,n_grids), &
    hpresz(nx,ny,nz,n_grids), hpresxy(nx,ny,nz,n_grids), &
    hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
    ! 
    orho(nx,ny,nz,n_grids),opx(nx,ny,nz,n_grids), &
    opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
    opresx(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
    opresy(nx,ny,nz,n_grids), opresxy(nx,ny,nz,n_grids), &
    opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
    !
    epres(nx,ny,nz,n_grids), &
    !
    bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
    bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids),bz0(nx,ny,nz,n_grids)
    !
    !     do x boundary conditions: set oxygen in solar wind
    !       at 2.5%
    !
    box=n_grids
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
            qrho(i,j,k,box)=srho
            qpx(i,j,k,box)=spx
            qpy(i,j,k,box)=spy
            qpz(i,j,k,box)=spz
            qpresx(i,j,k,box)=0.5*spress
            qpresy(i,j,k,box)=qpresx(i,j,k,box)
            qpresz(i,j,k,box)=qpresx(i,j,k,box)
            !
            hrho(i,j,k,box)=srho*rho_frac*frac_h
            hpx(i,j,k,box)=spx*rho_frac*frac_h
            hpy(i,j,k,box)=spy*rho_frac*frac_h
            hpz(i,j,k,box)=spz*rho_frac*frac_h
            hpresx(i,j,k,box)=0.5*spress*rho_frac
            hpresy(i,j,k,box)=hpresx(i,j,k,box)
            hpresz(i,j,k,box)=hpresx(i,j,k,box)
            !
            orho(i,j,k,box)=srho*ofrac*frac_o
            opx(i,j,k,box)=spx*ofrac*frac_o
            opy(i,j,k,box)=spy*ofrac*frac_o
            opz(i,j,k,box)=spz*ofrac*frac_o
            opresx(i,j,k,box)=0.5*spress*ofrac
            opresy(i,j,k,box)=opresx(i,j,k,box)
            opresz(i,j,k,box)=opresx(i,j,k,box)
            !
            epres(i,j,k,box)=qpresx(i,j,k,box)/ti_te
            !
            !        iostropic conditions: set off axis pressures to zero
            !
            qpresxy(i,j,k,box)=0.
            qpresxz(i,j,k,box)=0.
            qpresyz(i,j,k,box)=0.
            hpresxy(i,j,k,box)=0.
            hpresxz(i,j,k,box)=0.
            hpresyz(i,j,k,box)=0.
            opresxy(i,j,k,box)=0.
            opresxz(i,j,k,box)=0.
            opresyz(i,j,k,box)=0.
            !
            !       set bx of wind to zero, by to the wind values
            !          and subtract any geomagnetic field
            !
            bx(i,j,k,box)=sbx_wind-bx0(i,j,k,box)
            by(i,j,k,box)=sby_wind-by0(i,j,k,box)
            bz(i,j,k,box)=sbz_wind-bz0(i,j,k,box)
            !
            !       bx(i,j,k,box)=bx(2,j,k,box)
            !       by(i,j,k,box)=by(2,j,k,box)
            !       bz(i,j,k,box)=bz(2,j,k,box)
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
            qrho(nx,j,k,box)=qrho(nx1,j,k,box)
            qpx(nx,j,k,box)=abs(qpx(nx1,j,k,box))
            qpy(nx,j,k,box)=qpy(nx1,j,k,box)
            qpz(nx,j,k,box)=qpz(nx1,j,k,box)
            qpresx(nx,j,k,box)=qpresx(nx1,j,k,box)
            qpresy(nx,j,k,box)=qpresy(nx1,j,k,box)
            qpresz(nx,j,k,box)=qpresz(nx1,j,k,box)
            !
            hrho(nx,j,k,box)=hrho(nx1,j,k,box)
            hpx(nx,j,k,box)=abs(hpx(nx1,j,k,box))
            hpy(nx,j,k,box)=hpy(nx1,j,k,box)
            hpz(nx,j,k,box)=hpz(nx1,j,k,box)
            hpresx(nx,j,k,box)=hpresx(nx1,j,k,box)
            hpresy(nx,j,k,box)=hpresy(nx1,j,k,box)
            hpresz(nx,j,k,box)=hpresz(nx1,j,k,box)
            !
            orho(nx,j,k,box)=orho(nx1,j,k,box)
            opx(nx,j,k,box)=abs(opx(nx1,j,k,box))
            opy(nx,j,k,box)=opy(nx1,j,k,box)
            opz(nx,j,k,box)=opz(nx1,j,k,box)
            opresx(nx,j,k,box)=opresx(nx1,j,k,box)
	   		!
            opresy(nx,j,k,box)=opresy(nx1,j,k,box)
            opresz(nx,j,k,box)=opresz(nx1,j,k,box)
            !
            epres(nx,j,k,box)=epres(nx1,j,k,box)
            !
            !       continuous boundary condition
            !
            qpresxy(nx,j,k,box)=qpresxy(nx1,j,k,box)
            qpresxz(nx,j,k,box)=qpresxz(nx1,j,k,box)
            qpresyz(nx,j,k,box)=qpresyz(nx1,j,k,box)
            hpresxy(nx,j,k,box)=hpresxy(nx1,j,k,box)
            hpresxz(nx,j,k,box)=hpresxz(nx1,j,k,box)
            hpresyz(nx,j,k,box)=hpresyz(nx1,j,k,box)
            opresxy(nx,j,k,box)=opresxy(nx1,j,k,box)
            opresxz(nx,j,k,box)=opresxz(nx1,j,k,box)
            opresyz(nx,j,k,box)=opresyz(nx1,j,k,box)
            !
            bx(nx,j,k,box)=bx(nx1,j,k,box)
            by(nx,j,k,box)=by(nx1,j,k,box)
            bz(nx,j,k,box)=bz(nx1,j,k,box)
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
            qrho(i,1,k,box)=qrho(i1,2,k,box)
            qpresx(i,1,k,box)=qpresx(i1,2,k,box)
            qpresy(i,1,k,box)=qpresy(i1,2,k,box)
            qpresz(i,1,k,box)=qpresz(i1,2,k,box)
		    !
            qpx(i,1,k,box)=qpx(i1,2,k,box)
            qpy(i,1,k,box)=-abs(qpy(i1,2,k,box))
            qpz(i,1,k,box)=qpz(i1,2,k,box)
            !
            hrho(i,1,k,box)=hrho(i1,2,k,box)
            hpresx(i,1,k,box)=hpresx(i1,2,k,box)
            hpresy(i,1,k,box)=hpresy(i1,2,k,box)
            hpresz(i,1,k,box)=hpresz(i1,2,k,box)
		    !
            hpx(i,1,k,box)=hpx(i1,2,k,box)
            hpy(i,1,k,box)=-abs(hpy(i1,2,k,box))
            hpz(i,1,k,box)=hpz(i1,2,k,box)
            !
            orho(i,1,k,box)=orho(i1,2,k,box)
            opresx(i,1,k,box)=opresx(i1,2,k,box)
            opresy(i,1,k,box)=opresy(i1,2,k,box)
            opresz(i,1,k,box)=opresz(i1,2,k,box)
		    !
            opx(i,1,k,box)=opx(i1,2,k,box)
            opy(i,1,k,box)=-abs(opy(i1,2,k,box))
            opz(i,1,k,box)=opz(i1,2,k,box)
            !
            epres(i,1,k,box)=epres(i1,2,k,box)
            !
            !       continuous boundary condition
            !
            qpresxy(i,1,k,box)=qpresxy(i1,2,k,box)
            qpresxz(i,1,k,box)=qpresxz(i1,2,k,box)
            qpresyz(i,1,k,box)=qpresyz(i1,2,k,box)
            hpresxy(i,1,k,box)=hpresxy(i1,2,k,box)
            hpresxz(i,1,k,box)=hpresxz(i1,2,k,box)
            hpresyz(i,1,k,box)=hpresyz(i1,2,k,box)
            opresxy(i,1,k,box)=opresxy(i1,2,k,box)
            opresxz(i,1,k,box)=opresxz(i1,2,k,box)
            opresyz(i,1,k,box)=opresyz(i1,2,k,box)
            !
            by(i,1,k,box)=by(i1,2,k,box)
            !    +             +(by0(i1,2,k,box)-by0(i,1,k,box))
            bx(i,1,k,box)=bx(i1,2,k,box)
            !    +              +(bx0(i1,2,k,box)-bx0(i,1,k,box))
            bz(i,1,k,box)=bz(i1,2,k,box)
            !    +              +(bz0(i1,2,k,box)-bz0(i,1,k,box))
            !
            qrho(i,ny,k,box)=qrho(i1,ny1,k,box)
            qpresx(i,ny,k,box)=qpresx(i1,ny1,k,box)
            qpresy(i,ny,k,box)=qpresy(i1,ny1,k,box)
            qpresz(i,ny,k,box)=qpresz(i1,ny1,k,box)
			!
            qpx(i,ny,k,box)=qpx(i1,ny1,k,box)
            qpy(i,ny,k,box)=abs(qpy(i1,ny1,k,box))
            qpz(i,ny,k,box)=qpz(i1,ny1,k,box)
            !
            hrho(i,ny,k,box)=hrho(i1,ny1,k,box)
            hpresx(i,ny,k,box)=hpresx(i1,ny1,k,box)
            hpresy(i,ny,k,box)=hpresy(i1,ny1,k,box)
            hpresz(i,ny,k,box)=hpresz(i1,ny1,k,box)
            hpx(i,ny,k,box)=hpx(i1,ny1,k,box)
            hpy(i,ny,k,box)=abs(hpy(i1,ny1,k,box))
            hpz(i,ny,k,box)=hpz(i1,ny1,k,box)
            !
            orho(i,ny,k,box)=orho(i1,ny1,k,box)
            opresx(i,ny,k,box)=opresx(i1,ny1,k,box)
            opresy(i,ny,k,box)=opresy(i1,ny1,k,box)
            opresz(i,ny,k,box)=opresz(i1,ny1,k,box)
			!
            opx(i,ny,k,box)=opx(i1,ny1,k,box)
            opy(i,ny,k,box)=abs(opy(i1,ny1,k,box))
            opz(i,ny,k,box)=opz(i1,ny1,k,box)
            !
            epres(i,ny,k,box)=epres(i1,ny1,k,box)
            !
            !       continuous boundary condition
            !
            qpresxy(i,ny,k,box)=qpresxy(i1,ny1,k,box)
            qpresxz(i,ny,k,box)=qpresxz(i1,ny1,k,box)
            qpresyz(i,ny,k,box)=qpresyz(i1,ny1,k,box)
            hpresxy(i,ny,k,box)=hpresxy(i1,ny1,k,box)
            hpresxz(i,ny,k,box)=hpresxz(i1,ny1,k,box)
            hpresyz(i,ny,k,box)=hpresyz(i1,ny1,k,box)
            opresxy(i,ny,k,box)=opresxy(i1,ny1,k,box)
            opresxz(i,ny,k,box)=opresxz(i1,ny1,k,box)
            opresyz(i,ny,k,box)=opresyz(i1,ny1,k,box)
			!
            by(i,ny,k,box)=by(i1,ny1,k,box)
            !    +              +(by0(i1,ny1,k,box)-by0(i,ny,k,box))
            bx(i,ny,k,box)=bx(i1,ny1,k,box)
            !    +               +(bx0(i1,ny1,k,box)-bx0(i,ny,k,box))
            bz(i,ny,k,box)=bz(i1,ny1,k,box)
            !    +              +(bz0(i1,ny1,k,box)-bz0(i,ny,k,box))
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
            qrho(i,j,nz,box)=qrho(i1,j,nz1,box)
            qpresx(i,j,nz,box)=qpresx(i1,j,nz1,box)
            qpresy(i,j,nz,box)=qpresy(i1,j,nz1,box)
            qpresz(i,j,nz,box)=qpresz(i1,j,nz1,box)
			!
            qpx(i,j,nz,box)=qpx(i1,j,nz1,box)
            qpy(i,j,nz,box)=qpy(i1,j,nz1,box)
            qpz(i,j,nz,box)=abs(qpz(i1,j,nz1,box))
            !
            hrho(i,j,nz,box)=hrho(i1,j,nz1,box)
            hpresx(i,j,nz,box)=hpresx(i1,j,nz1,box)
            hpresy(i,j,nz,box)=hpresy(i1,j,nz1,box)
            hpresz(i,j,nz,box)=hpresz(i1,j,nz1,box)
			!
            hpx(i,j,nz,box)=hpx(i1,j,nz1,box)
            hpy(i,j,nz,box)=hpy(i1,j,nz1,box)
            hpz(i,j,nz,box)=abs(hpz(i1,j,nz1,box))
            !
            orho(i,j,nz,box)=orho(i1,j,nz1,box)
            opresx(i,j,nz,box)=opresx(i1,j,nz1,box)
            opresy(i,j,nz,box)=opresy(i1,j,nz1,box)
            opresz(i,j,nz,box)=opresz(i1,j,nz1,box)
			!
            opx(i,j,nz,box)=opx(i1,j,nz1,box)
            opy(i,j,nz,box)=opy(i1,j,nz1,box)
            opz(i,j,nz,box)=abs(opz(i1,j,nz1,box))
            !
            epres(i,j,nz,box)=epres(i1,j,nz1,box)
            !
            !       off diagnoal elements
            !
            qpresxy(i,j,nz,box)=qpresxy(i1,j,nz1,box)
            qpresxz(i,j,nz,box)=qpresxz(i1,j,nz1,box)
            qpresyz(i,j,nz,box)=qpresyz(i1,j,nz1,box)
            hpresxy(i,j,nz,box)=hpresxy(i1,j,nz1,box)
            hpresxz(i,j,nz,box)=hpresxz(i1,j,nz1,box)
            hpresyz(i,j,nz,box)=hpresyz(i1,j,nz1,box)
            opresxy(i,j,nz,box)=opresxy(i1,j,nz1,box)
            opresxz(i,j,nz,box)=opresxz(i1,j,nz1,box)
            opresyz(i,j,nz,box)=opresyz(i1,j,nz1,box)
            !
            bz(i,j,nz,box)=bz(i1,j,nz1,box)
            !    +            +bz0(i1,j,nz1,box)-bz0(i,j,nz,box)
            bx(i,j,nz,box)=bx(i1,j,nz1,box)
            !    +            +bx0(i1,j,nz1,box)-bx0(i,j,nz,box)
            by(i,j,nz,box)=by(i1,j,nz1,box)
            !    +            +by0(i1,j,nz1,box)-by0(i,j,nz,box)
        enddo
    enddo
    !
    !       symmetry boundary conditions at k=1 and k=2
    !
    do j=1,ny
        do i=2,nx
            !
            qrho(i,j,1,box)=qrho(i,j,2,box)
            qpx(i,j,1,box)=qpx(i,j,2,box)
            qpy(i,j,1,box)=qpy(i,j,2,box)
            qpz(i,j,1,box)=-abs(qpz(i,j,2,box))
            qpresx(i,j,1,box)=qpresx(i,j,2,box)
            qpresy(i,j,1,box)=qpresy(i,j,2,box)
            qpresz(i,j,1,box)=qpresz(i,j,2,box)
            !
            hrho(i,j,1,box)=hrho(i,j,2,box)
            hpx(i,j,1,box)=hpx(i,j,2,box)
            hpy(i,j,1,box)=hpy(i,j,2,box)
            hpz(i,j,1,box)=-abs(hpz(i,j,2,box))
            hpresx(i,j,1,box)=hpresx(i,j,2,box)
            hpresy(i,j,1,box)=hpresy(i,j,2,box)
            hpresz(i,j,1,box)=hpresz(i,j,2,box)
            !
            orho(i,j,1,box)=orho(i,j,2,box)
            opx(i,j,1,box)=opx(i,j,2,box)
            opy(i,j,1,box)=opy(i,j,2,box)
            opz(i,j,1,box)=-abs(opz(i,j,2,box))
            opresx(i,j,1,box)=opresx(i,j,2,box)
            opresy(i,j,1,box)=opresy(i,j,2,box)
            opresz(i,j,1,box)=opresz(i,j,2,box)
            !
            epres(i,j,1,box)=epres(i,j,2,box)
            !
            !       off diagonal elements
            !
            qpresxy(i,j,1,box)=qpresxy(i,j,2,box)
            qpresxz(i,j,1,box)=qpresxz(i,j,2,box)
            qpresyz(i,j,1,box)=qpresyz(i,j,2,box)
            hpresxy(i,j,1,box)=hpresxy(i,j,2,box)
            hpresxz(i,j,1,box)=hpresxz(i,j,2,box)
            hpresyz(i,j,1,box)=hpresyz(i,j,2,box)
            opresxy(i,j,1,box)=opresxy(i,j,2,box)
            opresxz(i,j,1,box)=opresxz(i,j,2,box)
            opresyz(i,j,1,box)=opresyz(i,j,2,box)
            !
            !       set bx of wind to zero, by to the wind values
            !          and subtract any geomagnetic field
            !       force bz to be positive to attain equilibrium
            !
            bx(i,j,1,box)=bx(i,j,2,box)
            by(i,j,1,box)=by(i,j,2,box)
            bz(i,j,1,box)=bz(i,j,2,box)
        enddo
    enddo
    !
    !     do corner lines at left-right back wall
    !
    do k=2,nz-1
        qrho(nx,1,k,box)=(qrho(nx1,1,k,box)+qrho(nx,2,k,box))/2.
        qrho(nx,ny,k,box)=(qrho(nx1,ny,k,box)+qrho(nx,ny1,k,box))/2.
        hrho(nx,1,k,box)=(hrho(nx1,1,k,box)+hrho(nx,2,k,box))/2.
        hrho(nx,ny,k,box)=(hrho(nx1,ny,k,box)+hrho(nx,ny1,k,box))/2.
        orho(nx,1,k,box)=(orho(nx1,1,k,box)+orho(nx,2,k,box))/2.
        orho(nx,ny,k,box)=(orho(nx1,ny,k,box)+orho(nx,ny1,k,box))/2.
        !
        qpresx(nx,1,k,box)=(qpresx(nx1,1,k,box) &
        +qpresx(nx,2,k,box))/2.
        qpresx(nx,ny,k,box)=(qpresx(nx1,ny,k,box) &
        +qpresx(nx,ny1,k,box))/2.
        qpresy(nx,1,k,box)=(qpresy(nx1,1,k,box) &
            +qpresy(nx,2,k,box))/2.
        qpresy(nx,ny,k,box)=(qpresy(nx1,ny,k,box) &
            +qpresy(nx,ny1,k,box))/2.
        qpresz(nx,1,k,box)=(qpresz(nx1,1,k,box) &
            +qpresz(nx,2,k,box))/2.
        qpresz(nx,ny,k,box)=(qpresz(nx1,ny,k,box) &
            +qpresz(nx,ny1,k,box))/2.
        !
        hpresx(nx,1,k,box)=(hpresx(nx1,1,k,box) &
        +hpresx(nx,2,k,box))/2.
        hpresx(nx,ny,k,box)=(hpresx(nx1,ny,k,box) &
        +hpresx(nx,ny1,k,box))/2.
        hpresy(nx,1,k,box)=(hpresy(nx1,1,k,box) &
            +hpresy(nx,2,k,box))/2.
        hpresy(nx,ny,k,box)=(hpresy(nx1,ny,k,box) &
            +hpresy(nx,ny1,k,box))/2.
        hpresz(nx,1,k,box)=(hpresz(nx1,1,k,box) &
            +hpresz(nx,2,k,box))/2.
        hpresz(nx,ny,k,box)=(hpresz(nx1,ny,k,box) &
            +hpresz(nx,ny1,k,box))/2.
        !
        opresx(nx,1,k,box)=(opresx(nx1,1,k,box) &
        +opresx(nx,2,k,box))/2.
        opresx(nx,ny,k,box)=(opresx(nx1,ny,k,box) &
        +opresx(nx,ny1,k,box))/2.
        opresy(nx,1,k,box)=(opresy(nx1,1,k,box) &
            +opresy(nx,2,k,box))/2.
        opresy(nx,ny,k,box)=(opresy(nx1,ny,k,box) &
            +opresy(nx,ny1,k,box))/2.
        opresz(nx,1,k,box)=(opresz(nx1,1,k,box) &
            +opresz(nx,2,k,box))/2.
        opresz(nx,ny,k,box)=(opresz(nx1,ny,k,box) &
            +opresz(nx,ny1,k,box))/2.
        !
        epres(nx,1,k,box)=(epres(nx1,1,k,box) &
        +epres(nx,2,k,box))/2.
        epres(nx,ny,k,box)=(epres(nx1,ny,k,box) &
        +epres(nx,ny1,k,box))/2.
        !
        !       off-diagonal elements
        !
        qpresxy(nx,1,k,box)=(qpresxy(nx1,1,k,box) &
            +qpresxy(nx,2,k,box))/2.
        qpresxy(nx,ny,k,box)=(qpresxy(nx1,ny,k,box) &
            +qpresxy(nx,ny1,k,box))/2.
        qpresxz(nx,1,k,box)=(qpresxz(nx1,1,k,box) &
            +qpresxz(nx,2,k,box))/2.
        qpresxz(nx,ny,k,box)=(qpresxz(nx1,ny,k,box) &
            +qpresxz(nx,ny1,k,box))/2.
        qpresyz(nx,1,k,box)=(qpresyz(nx1,1,k,box) &
            +qpresyz(nx,2,k,box))/2.
        qpresyz(nx,ny,k,box)=(qpresyz(nx1,ny,k,box) &
            +qpresyz(nx,ny1,k,box))/2.
        hpresxy(nx,1,k,box)=(hpresxy(nx1,1,k,box) &
            +hpresxy(nx,2,k,box))/2.
        hpresxy(nx,ny,k,box)=(hpresxy(nx1,ny,k,box) &
            +hpresxy(nx,ny1,k,box))/2.
        hpresxz(nx,1,k,box)=(hpresxz(nx1,1,k,box) &
            +hpresxz(nx,2,k,box))/2.
        hpresxz(nx,ny,k,box)=(hpresxz(nx1,ny,k,box) &
            +hpresxz(nx,ny1,k,box))/2.
        hpresyz(nx,1,k,box)=(hpresyz(nx1,1,k,box) &
            +hpresyz(nx,2,k,box))/2.
        hpresyz(nx,ny,k,box)=(hpresyz(nx1,ny,k,box) &
            +hpresyz(nx,ny1,k,box))/2.
        opresxy(nx,1,k,box)=(opresxy(nx1,1,k,box) &
            +opresxy(nx,2,k,box))/2.
        opresxy(nx,ny,k,box)=(opresxy(nx1,ny,k,box) &
            +opresxy(nx,ny1,k,box))/2.
        opresxz(nx,1,k,box)=(opresxz(nx1,1,k,box) &
            +opresxz(nx,2,k,box))/2.
        opresxz(nx,ny,k,box)=(opresxz(nx1,ny,k,box) &
            +opresxz(nx,ny1,k,box))/2.
        opresyz(nx,1,k,box)=(opresyz(nx1,1,k,box) &
            +opresyz(nx,2,k,box))/2.
        opresyz(nx,ny,k,box)=(opresyz(nx1,ny,k,box) &
            +opresyz(nx,ny1,k,box))/2.
        !
        qpx(nx,1,k,box)=(qpx(nx1,1,k,box)+qpx(nx,2,k,box))/2.
        qpx(nx,ny,k,box)=(qpx(nx1,ny,k,box)+qpx(nx,ny1,k,box))/2.
        !
        qpy(nx,1,k,box)=(qpy(nx1,1,k,box)+qpy(nx,2,k,box))/2.
        qpy(nx,ny,k,box)=(qpy(nx1,ny,k,box)+qpy(nx,ny1,k,box))/2.
        !
        qpz(nx,1,k,box)=(qpz(nx1,1,k,box)+qpz(nx,2,k,box))/2.
        qpz(nx,ny,k,box)=(qpz(nx1,ny,k,box)+qpz(nx,ny1,k,box))/2.
        !
        hpx(nx,1,k,box)=(hpx(nx1,1,k,box)+hpx(nx,2,k,box))/2.
        hpx(nx,ny,k,box)=(hpx(nx1,ny,k,box)+hpx(nx,ny1,k,box))/2.
        !
        hpy(nx,1,k,box)=(hpy(nx1,1,k,box)+hpy(nx,2,k,box))/2.
        hpy(nx,ny,k,box)=(hpy(nx1,ny,k,box)+hpy(nx,ny1,k,box))/2.
        !
        hpz(nx,1,k,box)=(hpz(nx1,1,k,box)+hpz(nx,2,k,box))/2.
        hpz(nx,ny,k,box)=(hpz(nx1,ny,k,box)+hpz(nx,ny1,k,box))/2.
        !
        opx(nx,1,k,box)=(opx(nx1,1,k,box)+opx(nx,2,k,box))/2.
        opx(nx,ny,k,box)=(opx(nx1,ny,k,box)+opx(nx,ny1,k,box))/2.
        !
        opy(nx,1,k,box)=(opy(nx1,1,k,box)+opy(nx,2,k,box))/2.
        opy(nx,ny,k,box)=(opy(nx1,ny,k,box)+opy(nx,ny1,k,box))/2.
        !
        opz(nx,1,k,box)=(opz(nx1,1,k,box)+opz(nx,2,k,box))/2.
        opz(nx,ny,k,box)=(opz(nx1,ny,k,box)+opz(nx,ny1,k,box))/2.
        !
        bx(nx,1,k,box)=bx(nx1,1,k,box)
        bx(nx,ny,k,box)=bx(nx1,ny,k,box)
        !
        by(nx,1,k,box)=(by(nx1,1,k,box) &
        +by(nx,2,k,box))/2.
        by(nx,ny,k,box)=(by(nx1,ny,k,box) &
        +by(nx,ny1,k,box))/2.
        !
        bz(nx,1,k,box)=(bz(nx1,1,k,box) &
        +bz(nx,2,k,box))/2.
        bz(nx,ny,k,box)=(bz(nx1,ny,k,box) &
        +bz(nx,ny1,k,box))/2.
        !
    enddo
    !
    !     do corner lines at top and bottom back wall
    !
    do j=2,ny-1
        qrho(nx,j,1,box)=qrho(nx,j,2,box)
        qrho(nx,j,nz,box)=qrho(nx1,j,nz,box)
        hrho(nx,j,1,box)=hrho(nx,j,2,box)
        hrho(nx,j,nz,box)=hrho(nx1,j,nz,box)
        orho(nx,j,1,box)=orho(nx,j,2,box)
        orho(nx,j,nz,box)=orho(nx1,j,nz,box)
        !
        qpresx(nx,j,1,box)=qpresx(nx,j,2,box)
        qpresx(nx,j,nz,box)=qpresx(nx1,j,nz,box)
        qpresy(nx,j,1,box)=qpresy(nx,j,2,box)
        qpresy(nx,j,nz,box)=qpresy(nx1,j,nz,box)
        qpresz(nx,j,1,box)=qpresz(nx,j,2,box)
        qpresz(nx,j,nz,box)=qpresz(nx1,j,nz,box)
        !
        hpresx(nx,j,1,box)=hpresx(nx,j,2,box)
        hpresx(nx,j,nz,box)=hpresx(nx1,j,nz,box)
        hpresy(nx,j,1,box)=hpresy(nx,j,2,box)
        hpresy(nx,j,nz,box)=hpresy(nx1,j,nz,box)
        hpresz(nx,j,1,box)=hpresz(nx,j,2,box)
        hpresz(nx,j,nz,box)=hpresz(nx1,j,nz,box)
        !
        opresx(nx,j,1,box)=opresx(nx,j,2,box)
        opresx(nx,j,nz,box)=opresx(nx1,j,nz,box)
        opresy(nx,j,1,box)=opresy(nx,j,2,box)
        opresy(nx,j,nz,box)=opresy(nx1,j,nz,box)
        opresz(nx,j,1,box)=opresz(nx,j,2,box)
        opresz(nx,j,nz,box)=opresz(nx1,j,nz,box)
        !
        epres(nx,j,1,box)=epres(nx,j,2,box)
        !
        !       off diagnonal elements
        !
        qpresxy(nx,j,1,box)=qpresxy(nx,j,2,box)
        qpresxy(nx,j,nz,box)=qpresxy(nx1,j,nz,box)
        qpresxz(nx,j,1,box)=qpresxz(nx,j,2,box)
        qpresxz(nx,j,nz,box)=qpresxz(nx1,j,nz,box)
        qpresyz(nx,j,1,box)=qpresyz(nx,j,2,box)
        qpresyz(nx,j,nz,box)=qpresyz(nx1,j,nz,box)
        hpresxy(nx,j,1,box)=hpresxy(nx,j,2,box)
        hpresxy(nx,j,nz,box)=hpresxy(nx1,j,nz,box)
        hpresxz(nx,j,1,box)=hpresxz(nx,j,2,box)
        hpresxz(nx,j,nz,box)=hpresxz(nx1,j,nz,box)
        hpresyz(nx,j,1,box)=hpresyz(nx,j,2,box)
        hpresyz(nx,j,nz,box)=hpresyz(nx1,j,nz,box)
        opresxy(nx,j,1,box)=opresxy(nx,j,2,box)
        opresxy(nx,j,nz,box)=opresxy(nx1,j,nz,box)
        opresxz(nx,j,1,box)=opresxz(nx,j,2,box)
        opresxz(nx,j,nz,box)=opresxz(nx1,j,nz,box)
        opresyz(nx,j,1,box)=opresyz(nx,j,2,box)
        opresyz(nx,j,nz,box)=opresyz(nx1,j,nz,box)
        !
        qpx(nx,j,1,box)=qpx(nx,j,2,box)
        qpx(nx,j,nz,box)=qpx(nx1,j,nz,box)
        !
        qpy(nx,j,1,box)=qpy(nx,j,2,box)
        qpy(nx,j,nz,box)=qpy(nx1,j,nz,box)
        !
        qpz(nx,j,1,box)=qpz(nx,j,2,box)
        qpz(nx,j,nz,box)=qpz(nx1,j,nz,box)
        !
        hpx(nx,j,1,box)=hpx(nx,j,2,box)
        hpx(nx,j,nz,box)=hpx(nx1,j,nz,box)
        !
        hpy(nx,j,1,box)=hpy(nx,j,2,box)
        hpy(nx,j,nz,box)=hpy(nx1,j,nz,box)
        !
        hpz(nx,j,1,box)=hpz(nx,j,2,box)
        hpz(nx,j,nz,box)=hpz(nx1,j,nz,box)
        !
        opx(nx,j,1,box)=opx(nx,j,2,box)
        opx(nx,j,nz,box)=opx(nx1,j,nz,box)
        !
        opy(nx,j,1,box)=opy(nx,j,2,box)
        opy(nx,j,nz,box)=opy(nx1,j,nz,box)
        !
        opz(nx,j,1,box)=opz(nx,j,2,box)
        opz(nx,j,nz,box)=opz(nx1,j,nz,box)
        !
        bx(nx,j,1,box)=bx(nx,j,2,box)
        bx(nx,j,nz,box)=bx(nx1,j,nz,box)
        !
        by(nx,j,1,box)=by(nx,j,2,box)
        by(nx,j,nz,box)=by(nx1,j,nz,box)
        !
        bz(nx,j,1,box)=bz(nx,j,2,box)
        bz(nx,j,nz,box)=bz(nx1,j,nz,box)
        !
    enddo
    !
    do i=2,nx
        qrho(i,1,1,box)=qrho(i,1,2,box)
        qrho(i,ny,1,box)=qrho(i,ny,2,box)
        qrho(i,1,nz,box)=qrho(i,1,nz1,box)
        qrho(i,ny,nz,box)=qrho(i,ny,nz1,box)
        !
        hrho(i,1,1,box)=hrho(i,1,2,box)
        hrho(i,ny,1,box)=hrho(i,ny,2,box)
        hrho(i,1,nz,box)=hrho(i,1,nz1,box)
        hrho(i,ny,nz,box)=hrho(i,ny,nz1,box)
        !
        orho(i,1,1,box)=orho(i,1,2,box)
        orho(i,ny,1,box)=orho(i,ny,2,box)
        orho(i,1,nz,box)=orho(i,1,nz1,box)
        orho(i,ny,nz,box)=orho(i,ny,nz1,box)
        !
        qpresx(i,1,1,box)=qpresx(i,1,2,box)
        qpresx(i,ny,1,box)=qpresx(i,ny,2,box)
        qpresx(i,1,nz,box)=qpresx(i,1,nz1,box)
        qpresx(i,ny,nz,box)=qpresx(i,ny,nz1,box)
        qpresy(i,1,1,box)=qpresy(i,1,2,box)
        qpresy(i,ny,1,box)=qpresy(i,ny,2,box)
        qpresy(i,1,nz,box)=qpresy(i,1,nz1,box)
        qpresy(i,ny,nz,box)=qpresy(i,ny,nz1,box)
        qpresz(i,1,1,box)=qpresz(i,1,2,box)
        qpresz(i,ny,1,box)=qpresz(i,ny,2,box)
        qpresz(i,1,nz,box)=qpresz(i,1,nz1,box)
        qpresz(i,ny,nz,box)=qpresz(i,ny,nz1,box)
        !
        hpresx(i,1,1,box)=hpresx(i,1,2,box)
        hpresx(i,ny,1,box)=hpresx(i,ny,2,box)
        hpresx(i,1,nz,box)=hpresx(i,1,nz1,box)
        hpresx(i,ny,nz,box)=hpresx(i,ny,nz1,box)
        hpresy(i,1,1,box)=hpresy(i,1,2,box)
        hpresy(i,ny,1,box)=hpresy(i,ny,2,box)
        hpresy(i,1,nz,box)=hpresy(i,1,nz1,box)
        hpresy(i,ny,nz,box)=hpresy(i,ny,nz1,box)
        hpresz(i,1,1,box)=hpresz(i,1,2,box)
        hpresz(i,ny,1,box)=hpresz(i,ny,2,box)
        hpresz(i,1,nz,box)=hpresz(i,1,nz1,box)
        hpresz(i,ny,nz,box)=hpresz(i,ny,nz1,box)
        !
        opresx(i,1,1,box)=opresx(i,1,2,box)
        opresx(i,ny,1,box)=opresx(i,ny,2,box)
        opresx(i,1,nz,box)=opresx(i,1,nz1,box)
        opresx(i,ny,nz,box)=opresx(i,ny,nz1,box)
        opresy(i,1,1,box)=opresy(i,1,2,box)
        opresy(i,ny,1,box)=opresy(i,ny,2,box)
        opresy(i,1,nz,box)=opresy(i,1,nz1,box)
        opresy(i,ny,nz,box)=opresy(i,ny,nz1,box)
        opresz(i,1,1,box)=opresz(i,1,2,box)
        opresz(i,ny,1,box)=opresz(i,ny,2,box)
        opresz(i,1,nz,box)=opresz(i,1,nz1,box)
        opresz(i,ny,nz,box)=opresz(i,ny,nz1,box)
        !
        epres(i,1,1,box)=epres(i,1,2,box)
        epres(i,ny,1,box)=epres(i,ny,2,box)
        epres(i,1,nz,box)=epres(i,1,nz1,box)
        epres(i,ny,nz,box)=epres(i,ny,nz1,box)
        !
        !       off-diagonal elements
        !
        qpresxy(i,1,1,box)=qpresxy(i,1,2,box)
        qpresxy(i,ny,1,box)=qpresxy(i,ny,2,box)
        qpresxy(i,1,nz,box)=qpresxy(i,1,nz1,box)
        qpresxy(i,ny,nz,box)=qpresxy(i,ny,nz1,box)
        qpresxz(i,1,1,box)=qpresxz(i,1,2,box)
        qpresxz(i,ny,1,box)=qpresxz(i,ny,2,box)
        qpresxz(i,1,nz,box)=qpresxz(i,1,nz1,box)
        qpresxz(i,ny,nz,box)=qpresxz(i,ny,nz1,box)
        qpresyz(i,1,1,box)=qpresyz(i,1,2,box)
        qpresyz(i,ny,1,box)=qpresyz(i,ny,2,box)
        qpresyz(i,1,nz,box)=qpresyz(i,1,nz1,box)
        qpresyz(i,ny,nz,box)=qpresyz(i,ny,nz1,box)
        !
        hpresxy(i,1,1,box)=hpresxy(i,1,2,box)
        hpresxy(i,ny,1,box)=hpresxy(i,ny,2,box)
        hpresxy(i,1,nz,box)=hpresxy(i,1,nz1,box)
        hpresxy(i,ny,nz,box)=hpresxy(i,ny,nz1,box)
        hpresxz(i,1,1,box)=hpresxz(i,1,2,box)
        hpresxz(i,ny,1,box)=hpresxz(i,ny,2,box)
        hpresxz(i,1,nz,box)=hpresxz(i,1,nz1,box)
        hpresxz(i,ny,nz,box)=hpresxz(i,ny,nz1,box)
        hpresyz(i,1,1,box)=hpresyz(i,1,2,box)
        hpresyz(i,ny,1,box)=hpresyz(i,ny,2,box)
        hpresyz(i,1,nz,box)=hpresyz(i,1,nz1,box)
        hpresyz(i,ny,nz,box)=hpresyz(i,ny,nz1,box)
        !
        opresxy(i,1,1,box)=opresxy(i,1,2,box)
        opresxy(i,ny,1,box)=opresxy(i,ny,2,box)
        opresxy(i,1,nz,box)=opresxy(i,1,nz1,box)
        opresxy(i,ny,nz,box)=opresxy(i,ny,nz1,box)
        opresxz(i,1,1,box)=opresxz(i,1,2,box)
        opresxz(i,ny,1,box)=opresxz(i,ny,2,box)
        opresxz(i,1,nz,box)=opresxz(i,1,nz1,box)
        opresxz(i,ny,nz,box)=opresxz(i,ny,nz1,box)
        opresyz(i,1,1,box)=opresyz(i,1,2,box)
        opresyz(i,ny,1,box)=opresyz(i,ny,2,box)
        opresyz(i,1,nz,box)=opresyz(i,1,nz1,box)
        opresyz(i,ny,nz,box)=opresyz(i,ny,nz1,box)
        !
        qpx(i,1,1,box)=qpx(i,1,2,box)
        qpx(i,ny,1,box)=qpx(i,ny,2,box)
        qpx(i,1,nz,box)=qpx(i,1,nz1,box)
        qpx(i,ny,nz,box)=qpx(i,ny,nz1,box)
        !
        qpy(i,1,1,box)=qpy(i,1,2,box)
        qpy(i,ny,1,box)=qpy(i,ny,2,box)
        qpy(i,1,nz,box)=qpy(i,1,nz1,box)
        qpy(i,ny,nz,box)=qpy(i,ny,nz1,box)
        !
        qpz(i,1,1,box)=qpz(i,1,2,box)
        qpz(i,ny,1,box)=qpz(i,ny,2,box)
        qpz(i,1,nz,box)=qpz(i,1,nz1,box)
        qpz(i,ny,nz,box)=qpz(i,ny,nz1,box)
        !
        hpx(i,1,1,box)=hpx(i,1,2,box)
        hpx(i,ny,1,box)=hpx(i,ny,2,box)
        hpx(i,1,nz,box)=hpx(i,1,nz1,box)
        hpx(i,ny,nz,box)=hpx(i,ny,nz1,box)
        !
        hpy(i,1,1,box)=hpy(i,1,2,box)
        hpy(i,ny,1,box)=hpy(i,ny,2,box)
        hpy(i,1,nz,box)=hpy(i,1,nz1,box)
        hpy(i,ny,nz,box)=hpy(i,ny,nz1,box)
        !
        hpz(i,1,1,box)=hpz(i,1,2,box)
        hpz(i,ny,1,box)=hpz(i,ny,2,box)
        hpz(i,1,nz,box)=hpz(i,1,nz1,box)
        hpz(i,ny,nz,box)=hpz(i,ny,nz1,box)
        !
        opx(i,1,1,box)=opx(i,1,2,box)
        opx(i,ny,1,box)=opx(i,ny,2,box)
        opx(i,1,nz,box)=opx(i,1,nz1,box)
        opx(i,ny,nz,box)=opx(i,ny,nz1,box)
        !
        opy(i,1,1,box)=opy(i,1,2,box)
        opy(i,ny,1,box)=opy(i,ny,2,box)
        opy(i,1,nz,box)=opy(i,1,nz1,box)
        opy(i,ny,nz,box)=opy(i,ny,nz1,box)
        !
        opz(i,1,1,box)=opz(i,1,2,box)
        opz(i,ny,1,box)=opz(i,ny,2,box)
        opz(i,1,nz,box)=opz(i,1,nz1,box)
        opz(i,ny,nz,box)=opz(i,ny,nz1,box)
        !
        bx(i,1,1,box)=bx(i,1,2,box)
        bx(i,ny,1,box)=bx(i,ny,2,box)
        bx(i,1,nz,box)=bx(i,1,nz1,box)
        bx(i,ny,nz,box)=bx(i,ny,nz1,box)
        !
        by(i,1,1,box)=by(i,1,2,box)
        by(i,ny,1,box)=by(i,ny,2,box)
        by(i,1,nz,box)=by(i,1,nz1,box)
        by(i,ny,nz,box)=by(i,ny,nz1,box)
        !
        bz(i,1,1,box)=bz(i,1,2,box)
        bz(i,ny,1,box)=bz(i,ny,2,box)
        bz(i,1,nz,box)=bz(i,1,nz1,box)
        bz(i,ny,nz,box)=bz(i,ny,nz1,box)
        !
    enddo
    !
    return
end
