!
!	Apply flux correction smoothing technique
!
subroutine flux_correct(qrho,qpresx,qpresy,qpresz, &
	qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
	wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
	wrkqpresxy,wrkqpresxz,wrkqpresyz, &
	wrkqpx,wrkqpy,wrkqpz, &
	oldqrho,oldqpresx,oldqpresy,oldqpresz, &
	oldqpresxy,oldqpresxz,oldqpresyz, &
	oldqpx,oldqpy,oldqpz, &

	hrho,hpresx,hpresy,hpresz, &
	hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
	wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
	wrkhpresxy,wrkhpresxz,wrkhpresyz, &
	wrkhpx,wrkhpy,wrkhpz, &
	oldhrho,oldhpresx,oldhpresy,oldhpresz, &
	oldhpresxy,oldhpresxz,oldhpresyz, &
	oldhpx,oldhpy,oldhpz, &

	orho,opresx,opresy,opresz, &
	opresxy,opresxz,opresyz,opx,opy,opz, &
	wrkorho,wrkopresx,wrkopresy,wrkopresz, &
	wrkopresxy,wrkopresxz,wrkopresyz, &
	wrkopx,wrkopy,wrkopz, &
	oldorho,oldopresx,oldopresy,oldopresz, &
	oldopresxy,oldopresxz,oldopresyz, &
	oldopx,oldopy,oldopz, &

	epres,wrkepres,oldepres, &
	bx,by,bz,wrkbx,wrkby,wrkbz, &
	oldbx,oldby,oldbz,vvx,vvy,vvz, &
	nx,ny,nz,n_grids,box,difrho,diferg,xspac, &
	isotropic)

	integer box
	dimension qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
	qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
	qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
	qpresyz(nx,ny,nz,n_grids),qpx(nx,ny,nz,n_grids), &
	qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &

	hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
	hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
	hpresxy(nx,ny,nz,n_grids),hpresxz(nx,ny,nz,n_grids), &
	hpresyz(nx,ny,nz,n_grids),hpx(nx,ny,nz,n_grids), &
	hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &

	orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
	opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
	opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
	opresyz(nx,ny,nz,n_grids),opx(nx,ny,nz,n_grids), &
	opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &

	bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids), &
	bz(nx,ny,nz,n_grids),epres(nx,ny,nz,n_grids), &
	vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz), &
	xspac(n_grids)

	dimension wrkqrho(nx,ny,nz,n_grids),wrkqpx(nx,ny,nz,n_grids), &
	wrkqpy(nx,ny,nz,n_grids),wrkqpz(nx,ny,nz,n_grids), &
	wrkqpresx(nx,ny,nz,n_grids),wrkqpresy(nx,ny,nz,n_grids), &
	wrkqpresz(nx,ny,nz,n_grids),wrkqpresxy(nx,ny,nz,n_grids), &
	wrkqpresxz(nx,ny,nz,n_grids),wrkqpresyz(nx,ny,nz,n_grids), &

	wrkhrho(nx,ny,nz,n_grids),wrkhpx(nx,ny,nz,n_grids), &
	wrkhpy(nx,ny,nz,n_grids),wrkhpz(nx,ny,nz,n_grids), &
	wrkhpresx(nx,ny,nz,n_grids),wrkhpresy(nx,ny,nz,n_grids), &
	wrkhpresz(nx,ny,nz,n_grids),wrkhpresxy(nx,ny,nz,n_grids), &
	wrkhpresxz(nx,ny,nz,n_grids),wrkhpresyz(nx,ny,nz,n_grids), &

	wrkorho(nx,ny,nz,n_grids),wrkopx(nx,ny,nz,n_grids), &
	wrkopy(nx,ny,nz,n_grids),wrkopz(nx,ny,nz,n_grids), &
	wrkopresx(nx,ny,nz,n_grids),wrkopresy(nx,ny,nz,n_grids), &
	wrkopresz(nx,ny,nz,n_grids),wrkopresxy(nx,ny,nz,n_grids), &
	wrkopresxz(nx,ny,nz,n_grids),wrkopresyz(nx,ny,nz,n_grids), &

	wrkbx(nx,ny,nz,n_grids),wrkby(nx,ny,nz,n_grids), &
	wrkbz(nx,ny,nz,n_grids),wrkepres(nx,ny,nz,n_grids)

	dimension oldqrho(nx,ny,nz,n_grids),oldqpx(nx,ny,nz,n_grids), &
	oldqpy(nx,ny,nz,n_grids),oldqpz(nx,ny,nz,n_grids), &
	oldqpresx(nx,ny,nz,n_grids),oldqpresy(nx,ny,nz,n_grids), &
	oldqpresz(nx,ny,nz,n_grids),oldqpresxy(nx,ny,nz,n_grids), &
	oldqpresxz(nx,ny,nz,n_grids),oldqpresyz(nx,ny,nz,n_grids), &

	oldhrho(nx,ny,nz,n_grids),oldhpx(nx,ny,nz,n_grids), &
	oldhpy(nx,ny,nz,n_grids),oldhpz(nx,ny,nz,n_grids), &
	oldhpresx(nx,ny,nz,n_grids),oldhpresy(nx,ny,nz,n_grids), &
	oldhpresz(nx,ny,nz,n_grids),oldhpresxy(nx,ny,nz,n_grids), &
	oldhpresxz(nx,ny,nz,n_grids),oldhpresyz(nx,ny,nz,n_grids), &

	oldorho(nx,ny,nz,n_grids),oldopx(nx,ny,nz,n_grids), &
	oldopy(nx,ny,nz,n_grids),oldopz(nx,ny,nz,n_grids), &
	oldopresx(nx,ny,nz,n_grids),oldopresy(nx,ny,nz,n_grids), &
	oldopresz(nx,ny,nz,n_grids),oldopresxy(nx,ny,nz,n_grids), &
	oldopresxz(nx,ny,nz,n_grids),oldopresyz(nx,ny,nz,n_grids), &

	oldbx(nx,ny,nz,n_grids),oldby(nx,ny,nz,n_grids), &
	oldbz(nx,ny,nz,n_grids),oldepres(nx,ny,nz,n_grids)

	!	write(*,*) 'In flux correct with: nx, ny, nz, n_grids, box, &
	!		chirho, diferg, xspac'
	!	write(*,*) nx, ny, nz, n_grids, box, chirho, diferg, xspac

	rx=xspac(box)
	ry=xspac(box)
	rz=xspac(box)
	chifcs=difrho
	call fcsmooth(qrho,oldqrho,wrkqrho,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)

	call fcsmooth(qpresx,oldqpresx,wrkqpresx,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(qpresy,oldqpresy,wrkqpresy,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(qpresz,oldqpresz,wrkqpresz,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(qpresxy,oldqpresxy,wrkqpresxy,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	call fcsmooth(qpresxz,oldqpresxz,wrkqpresxz,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	call fcsmooth(qpresyz,oldqpresyz,wrkqpresyz,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	
	call fcsmooth(qpx,oldqpx,wrkqpx,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(qpy,oldqpy,wrkqpy,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(qpz,oldqpz,wrkqpz,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	
	call fcsmooth(hrho,oldhrho,wrkhrho,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	
	call fcsmooth(hpresx,oldhpresx,wrkhpresx,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(hpresy,oldhpresy,wrkhpresy,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(hpresz,oldhpresz,wrkhpresz,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(hpresxy,oldhpresxy,wrkhpresxy,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	call fcsmooth(hpresxz,oldhpresxz,wrkhpresxz,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	call fcsmooth(hpresyz,oldhpresyz,wrkhpresyz,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	
	call fcsmooth(hpx,oldhpx,wrkhpx,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(hpy,oldhpy,wrkhpy,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(hpz,oldhpz,wrkhpz,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	
	call fcsmooth(orho,oldorho,wrkorho,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	
	call fcsmooth(opresx,oldopresx,wrkopresx,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(opresy,oldopresy,wrkopresy,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(opresz,oldopresz,wrkopresz,nx,ny,nz,n_grids,box, &
		chifcs,vvx,vvy,vvz)
	call fcsmooth(opresxy,oldopresxy,wrkopresxy,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	call fcsmooth(opresxz,oldopresxz,wrkopresxz,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	call fcsmooth(opresyz,oldopresyz,wrkopresyz,nx,ny,nz,n_grids, &
		box,chifcs,vvx,vvy,vvz)
	
	call fcsmooth(opx,oldopx,wrkopx,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(opy,oldopy,wrkopy,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(opz,oldopz,wrkopz,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	
	call fcsmooth(epres,oldepres,wrkepres,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	
	chifcs=diferg
	call fcsmooth(bx,oldbx,wrkbx,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(by,oldby,wrkby,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	call fcsmooth(bz,oldbz,wrkbz,nx,ny,nz,n_grids,box,chifcs, &
		vvx,vvy,vvz)
	
	return
end subroutine lap_smooth
