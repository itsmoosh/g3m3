!
!	This subroutine applies boundary conditions to all grid types
!	of the system and at any irregular boundaries
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

	integer, parameter :: dp = kind(1.d0)
	integer box
	real(dp) t_old(n_grids), t_new(n_grids), t_grid
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
	
	dimension qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
	qpresyz(nx,ny,nz,n_grids),hpresxy(nx,ny,nz,n_grids), &
	hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
	opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
	opresyz(nx,ny,nz,n_grids)
	
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
	
	dimension work(nx,ny,nz)
	
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)
	
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
	
	return
end subroutine bndry_flanks
