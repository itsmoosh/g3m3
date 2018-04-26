!
!   This subroutine prints simulation data for plotting. It is derived
!	from a cutter file, and prints ASCII text files intended for matplotlib.
!	
subroutine write_graphing_data( &
	nx,ny,nz,n_grids, limit, &
	mbndry, num_zqt, msrf, mmid, mzero, &
	xspac, grid_minvals, grid_maxvals, ut, &
	t, &
	qrho,qpx,qpy,qpz, &
	qpresx,qpresy,qpresz, &
	qpresxy,qpresxz,qpresyz, &
	hrho,hpx,hpy,hpz, &
	hpresx,hpresy,hpresz, &
	hpresxy,hpresxz,hpresyz, &
	orho,opx,opy,opz, &
	opresx,opresy,opresz, &
	opresxy,opresxz,opresyz, &
	bx,by,bz,epres,bx0,by0,bz0, &
	parm_srf,parm_mid,parm_zero, &
	ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
	bsx,bsy,bsz, &
	rmassq,rmassh,rmasso, &
	reynolds, resistive, resist, &
	curx,cury,curz, &
	ncraft, xcraft, re_equiv, b_equiv, v_equiv, t_equiv, &
	ti_te, rho_equiv, planet_rad, planet_per, moon_rad, &
	run_name, dummy_f, nplots)

	use astrometry
	implicit none

	!	double precision
	integer, parameter :: dp = kind(1.d0)

	!	graphing parameters
	integer, parameter :: n_cuts=3	! Number of cuts to make
	integer, parameter :: nlines=479 ! Synthetic trajectory length
	integer, parameter :: skip=1 ! Sample skip for flythrough data
	integer, parameter :: vbins=100 ! Number of energy/velocity bins on Y-axis
	integer, parameter :: phi_res=10 ! Number of angle bins on Y-axis

	!	file naming
	character*32, parameter :: python_dir = "figures/"
	character*32, parameter :: python_data = "figures/data/"
	character*32, parameter :: python_plotter = "plot_data.py"
	character*2, parameter :: cut_label(n_cuts) = ['xy', 'xz', 'yz']
	character*5, parameter :: fname1 = 'plas'
	character*5, parameter :: fname2 = 'flow'
	character*5, parameter :: fname3 = 'pres'
	character*5, parameter :: fname4 = 'bande'
	character*5, parameter :: fname5 = 'model'

	!	file IDs
	integer, parameter :: plas_f(n_cuts) = [111, 121, 131]
	integer, parameter :: flow_f(n_cuts) = [112, 122, 132]
	integer, parameter :: pres_f(n_cuts) = [113, 123, 133]
	integer, parameter :: bande_f(n_cuts) = [114, 124, 134]
	integer, parameter :: model_f(n_cuts) = [115, 125, 135]

!       header information
	character*300 plas_header, flow_header, bande_header, model_header, &
		efld_header, pres_header
	character*53, parameter :: header_intro = "This file contains data for the following quantities:"
	character*10 ut_string

	!**********************
	!	Input variables
	!**********************
	!	grid quantities
	integer, intent(in) :: nx, ny, nz, n_grids
	real, intent(in) :: limit, xspac(n_grids), &
		grid_minvals(3,n_grids), grid_maxvals(3,n_grids)
	integer, intent(in) :: mbndry, num_zqt, msrf, mmid, mzero
	!	physics sim quantities
	real(dp), intent(in) :: ut, t
	real, intent(in) ::	qrho(nx,ny,nz,n_grids), qpx(nx,ny,nz,n_grids), &
		qpy(nx,ny,nz,n_grids), qpz(nx,ny,nz,n_grids), &
		qpresx(nx,ny,nz,n_grids), qpresy(nx,ny,nz,n_grids), &
		qpresz(nx,ny,nz,n_grids), qpresxy(nx,ny,nz,n_grids), &
		qpresxz(nx,ny,nz,n_grids), qpresyz(nx,ny,nz,n_grids)
	real, intent(in) :: hrho(nx,ny,nz,n_grids), hpx(nx,ny,nz,n_grids), &
		hpy(nx,ny,nz,n_grids), hpz(nx,ny,nz,n_grids), &
		hpresx(nx,ny,nz,n_grids), hpresy(nx,ny,nz,n_grids), &
		hpresz(nx,ny,nz,n_grids), hpresxy(nx,ny,nz,n_grids), &
		hpresxz(nx,ny,nz,n_grids), hpresyz(nx,ny,nz,n_grids)
	real, intent(in) :: orho(nx,ny,nz,n_grids), opx(nx,ny,nz,n_grids), &
		opy(nx,ny,nz,n_grids), opz(nx,ny,nz,n_grids), &
		opresx(nx,ny,nz,n_grids), opresy(nx,ny,nz,n_grids), &
		opresz(nx,ny,nz,n_grids), opresxy(nx,ny,nz,n_grids), &
		opresxz(nx,ny,nz,n_grids), opresyz(nx,ny,nz,n_grids)
	real, intent(in) :: bx(nx,ny,nz,n_grids), by(nx,ny,nz,n_grids), &
		bz(nx,ny,nz,n_grids), epres(nx,ny,nz,n_grids), &
		bx0(nx,ny,nz,n_grids), by0(nx,ny,nz,n_grids), bz0(nx,ny,nz,n_grids)
	!	zeroed quantities
	real, intent(in) :: parm_srf(mbndry,num_zqt,msrf), &
		parm_mid(mbndry,num_zqt,mmid), parm_zero(mbndry,num_zqt,mzero)
	integer, intent(in) :: ijzero(mbndry,3,mzero), numzero(mbndry), &
		ijmid(mbndry,3,mmid), nummid(mbndry), &
		ijsrf(mbndry,3,msrf), numsrf(mbndry)
	real, intent(in) :: bsx(nx,ny,nz), bsy(nx,ny,nz), bsz(nx,ny,nz)
	real, intent(in) :: rmassq, rmassh, rmasso, &
		reynolds, resistive(nx,ny,nz,mbndry), resist
	real, intent(in) :: curx(nx,ny,nz), cury(nx,ny,nz), curz(nx,ny,nz)
	integer, intent(in) :: ncraft
	real, intent(in) :: xcraft(4,ncraft), re_equiv, b_equiv, v_equiv, &
		t_equiv, ti_te, rho_equiv, planet_rad, planet_per, moon_rad
	character*8, intent(in) :: run_name
	integer, intent(in) :: dummy_f

	integer, intent(out) :: nplots

!     --------------------------------------------------------------

	!	grid spacing working values
	real rx,ry,rz

	!	string dump working variables
	character*40 fname_ending(n_cuts)
	character*40 fname_starting

!      grid limits now set by grd_min grd_max arrays
!      ncore denotes couser grid to be hollowed out by fine grid
!      nbndry denotes finer grid to which coarser grid sets flanks
!      xspac is the relative grid spacing relative to inner grid system
!      main_n_grids gives the box number from which the fine gridding is
!                   interpolated

	integer ncore(n_grids), nbndry(n_grids)
	integer box

	!      xcraft is the actual position of the spacecraft in RE
	!          4th dimension of the actual time
	!      zcraft is the future position of the spacecraft in RE
	!      rcraft is the position of the spacecraft for which
	!           IMF is reference. NO alteration from boundary conditions applied
	!

	!	physics plasma quantities
	real qpara(nx,ny,nz,n_grids), qperp(nx,ny,nz,n_grids), &
		qcross(nx,ny,nz,n_grids)
	real hpara(nx,ny,nz,n_grids), hperp(nx,ny,nz,n_grids), &
		hcross(nx,ny,nz,n_grids)
	real opara(nx,ny,nz,n_grids), operp(nx,ny,nz,n_grids), &
		ocross(nx,ny,nz,n_grids)
	real qvx, qvy, qvz, ovx, ovy, ovz, hvx, hvy, hvz
	real tvx(nx,ny,nz), tvy(nx,ny,nz), tvz(nx,ny,nz)
	real evelx(nx,ny,nz,n_grids), evely(nx,ny,nz,n_grids), &
		evelz(nx,ny,nz,n_grids)
	real evx(nx,ny,nz), evy(nx,ny,nz), evz(nx,ny,nz)
	real vvx(nx,ny,nz), vvy(nx,ny,nz), vvz(nx,ny,nz)
	real ex(nx,ny,nz), ey(nx,ny,nz), ez(nx,ny,nz)

	!	Unperturbed quantities
	real bxs(nx,ny,nz,n_grids),bys(nx,ny,nz,n_grids),bzs(nx,ny,nz,n_grids)
	real br0(nx,ny,nz,n_grids),bt0(nx,ny,nz,n_grids),bp0(nx,ny,nz,n_grids)
	real bxt(nx,ny,nz),byt(nx,ny,nz),bzt(nx,ny,nz),btot(nx,ny,nz)

	real curx_all, cury_all, curz_all

	!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*

	integer a,b,i,j,k,m,n

	real aqpres, ahpres, aopres, aepres, qden, hden, oden, tden, &
		eden, qtemp, htemp, otemp, etemp, qpress, hpress, &
		opress, epress, qdens, hdens, odens, edens, lqpara, lqperp, lqcross, &
		lhpara, lhperp, lhcross, lopara, loperp, locross, &
		qtemps, htemps, otemps, etemps, abx, aby, abz, bdx, bdy, bdz, &
		ri, rj, rk, rad, bmag, bdipmag, time, ub0x, ub0y, ub0z, &
		ri_moon, rj_moon, rk_moon, aoden, aqden, ahden, aeden, &
		aotemp, ahtemp, aqtemp, aetemp, aefldx, aefldy, aefldz, &
		net_flow, v_alf, alfven_mach, unetx, unety, unetz, emag

	character*1 boxchar
	character*3 nplots_char

	character*64 wd1(n_cuts), wd2(n_cuts), wd3(n_cuts), wd4(n_cuts), &
		wd5(n_cuts), wd6(n_cuts), wd7(n_cuts), wd8(n_cuts), &
		wd9(n_cuts), wd0(n_cuts)
	character*64 dummy_fname, data_grep
	logical dummy_exists
	logical cuts(n_cuts)
	integer cuts_vals(n_cuts)

	!     rho,pres,erg,px,py are the density, pressure, energy and momentum
	!              in the x and y directions, respectively, of the fluid
	!       the indices 1, 2 is required to store old and new values
	!
	!     ijsrf give position of ionosphere in grid units - plasma
	!           parameters stored in parm_srf
	!     ijmid gives intermediate boundary - say representing the
	!            atmosphere - plasma parameters stored in parm_mid
	!     ijzero gives position of all grid units interior to surface
	!
	!     frho,ferg and fpx,fpy are the estimates of the fluid quantities
	!           at n+1/2 as defined by the Lax-Wendroff scheme
	!       the index of 3 is needed to store adjacent x values for these fns
	!
	!     d_min is the minimum allowable density
	!     stepsz is the size of the time step in terms of delta_x,y/(|umax|+c_s)
	!
	!     system dimensions are nx,ny,nz
	!     variable grid spacing enabled with rxyz >1
	!           rx,ry,rz should not be set identically to zero

	real temp_equiv	!1.0432e4 with v_equiv = 1000.
	real pres_equiv	!3.34e-10 with rho_equiv = 0.2, v_equiv = 1000.
	real cur_equiv	!1.3725e-9 with re_equiv = 0.2, planet_rad = 60268., b_equiv = 20.79 (rho_equiv = 0.2, v_equiv = 1000.)
	real e_equiv	!20.79e6 with b_equiv = 20.79 and v_equiv = 1000.


!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
!
!	Executions
!
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~


	! These multipliers are in SI units but are calculated from sim units.
	! Numerical factors appear for the following conversions:
	! rho_equiv is in cm^-3 of rmassq species
	! v_equiv is in km/s

	temp_equiv = m_prot * (v_equiv*1.e3)**2 / q_elec	!	In eV
	pres_equiv = rho_equiv*1.e6*m_prot * (v_equiv*1.e3)**2	! In Pa
	cur_equiv = b_equiv / mu0 / (planet_rad*1.e3 * re_equiv)	! In nA/m^2
	e_equiv = b_equiv * v_equiv*1.e3
	write(ut_string,'(F9.4)') ut

!
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
!	large grid system
!	(Not currently implemented)
!		ncore -> grid gets all the information from this grid no.
!		nbdry -> which grid data will be applied to this grid no. for bndry


	!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
!	do box=1,n_grids

!	enddo !loop over box

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	!	Check number of data files printed for this run_name,
	!		and increment it. This number is the only output from this subroutine.

	fname_starting = trim(python_data)//trim(run_name)//'_'

	dummy_fname = trim(python_data)//"dummy.txt"
	inquire(file=trim(dummy_fname), exist=dummy_exists)
	if(dummy_exists) call system ('rm '//trim(dummy_fname))
	
	data_grep = trim(fname_starting)//trim(fname1)//'_1'//trim(cut_label(1))//'_t???.dat'
	call system("ls 2>/dev/null -1q "//trim(data_grep)//" | wc -l > "//trim(dummy_fname))

	open(dummy_f,file=trim(dummy_fname),status='unknown',form='formatted')
		read(dummy_f,*) nplots
	close(dummy_f)
	call system ('rm '//trim(dummy_fname))

	nplots = nplots + 1
	if(nplots .gt. 999) write(*,*) "ERROR: 1000 plots created for this run_name. Max is 999."
	
	!	****************************************************
	!	Make a data cube for post-processing with matplotlib
	!	****************************************************

	write(*,'(A)',advance='no') "Writing data to files..."

	!	********************
	!		Low Res Data
	!	********************

	do box=1,n_grids
		write(boxchar,'(I1)') box
		write(nplots_char,'(I0.3)') nplots

		do m=1,n_cuts
			fname_ending(m) = '_'//trim(adjustl(boxchar))//trim(cut_label(m))//'_t'//trim(adjustl(nplots_char))//'.dat'

			!	************
			!	Create files
			!	************

			wd1(m) = trim(fname_starting)//trim(fname1)//trim(fname_ending(m))
			wd2(m) = trim(fname_starting)//trim(fname2)//trim(fname_ending(m))
			wd3(m) = trim(fname_starting)//trim(fname3)//trim(fname_ending(m))
			wd4(m) = trim(fname_starting)//trim(fname4)//trim(fname_ending(m))
			wd5(m) = trim(fname_starting)//trim(fname5)//trim(fname_ending(m))

			open(plas_f(m),file=wd1(m),status='replace',form='formatted')
			open(flow_f(m),file=wd2(m),status='replace',form='formatted')
			open(pres_f(m),file=wd3(m),status='replace',form='formatted')
			open(bande_f(m),file=wd4(m),status='replace',form='formatted')
			open(model_f(m),file=wd5(m),status='replace',form='formatted')

			!	****************
			!	Generate headers
			!	****************

			write(plas_header,'(A)') header_intro//new_line('A')// &
				"  Plasma densities and temperatures for each species, in log(1/cm^3) and eV."//new_line('A')// &
				'  run name:      '// run_name//new_line('A')// &
				'  box, cut:      '// boxchar//', '// cut_label(m)//'-plane'//new_line('A')// &
				'  ut, #in seq: '// ut_string// ', '//nplots_char
			write(flow_header,'(A)') header_intro//new_line('A')// &
				"  Flow velocity components for each species and net flow, in km/s, and Alfven Mach number."//new_line('A')// &
				'  run name:      '// run_name//new_line('A')// &
				'  box, cut:      '// boxchar//', '// cut_label(m)//'-plane'//new_line('A')// &
				'  ut, #in seq: '// ut_string// ', '//nplots_char
			write(pres_header,'(A)') header_intro//new_line('A')// &
				"  Pressure tensor terms for each species: parallel, perpendicular, and cross, and isotropic electron pressure, all in log(nPa)."//new_line('A')// &
				'  run name:      '// run_name//new_line('A')// &
				'  box, cut:      '// boxchar//', '// cut_label(m)//'-plane'//new_line('A')// &
				'  ut, #in seq: '// ut_string// ', '//nplots_char
			write(bande_header,'(A)') header_intro//new_line('A')// &
				"  Net magnetic and electric field components and magnitudes, and current density components, in nT, nV/m, and nA/m^2,."//new_line('A')// &
				'  run name:      '// run_name//new_line('A')// &
				'  box, cut:      '// boxchar//', '// cut_label(m)//'-plane'//new_line('A')// &
				'  ut, #in seq: '// ut_string// ', '//nplots_char
			write(model_header,'(A)') header_intro//new_line('A')// &
				"  Planetocentric distance in planetary radii, magnetic field comps for perturbed and background, and BG magnitude, all in nT."//new_line('A')// &
				'  run name:      '// run_name//new_line('A')// &
				'  box, cut:      '// boxchar//', '// cut_label(m)//'-plane'//new_line('A')// &
				'  ut, #in seq: '// ut_string// ', '//nplots_char

			write(plas_f(m),'(A)') trim(plas_header)
			write(flow_f(m),'(A)') trim(flow_header)
			write(pres_f(m),'(A)') trim(pres_header)
			write(bande_f(m),'(A)') trim(bande_header)
			write(model_f(m),'(A)') trim(model_header)

			write(plas_f(m),'(8(A14))') &
				'log(qden/cc)','qtemp(eV)', 'log(hden/cc)','htemp(eV)', &	! Densities are in units of # per cc
				'log(oden/cc)','otemp(eV)', 'log(eden/cc)','etemp(eV)'
			write(flow_f(m),'(13(A14))') &
				'qvx(km/s)','qvy(km/s)','qvz(km/s)', &
				'hvx(km/s)','hvy(km/s)','hvz(km/s)', &
				'ovx(km/s)','ovy(km/s)','ovz(km/s)', &
				'u_netx(km/s)','u_nety(km/s)','u_netz(km/s)', &
				'|v|/V_Alf'
			write(pres_f(m),'(10(A14))') &
				'log_qpar(nPa)','log_qprp(nPa)','log_qcrs(nPa)', &
				'log_hpar(nPa)','log_hprp(nPa)','log_hcrs(nPa)', &
				'log_opar(nPa)','log_oprp(nPa)','log_ocrs(nPa)', &
				'log_eprs(nPa)'
			write(bande_f(m),'(11(A14))') &
				'Bx(nT)','By(nT)','Bz(nT)','|Btot|(nT)', &
				'Ex(nV/m)','Ey(nV/m)','Ez(nV/m)','|Etot|(nV/m)', &
				'Jx(nA/m^2)','Jy(nA/m^2)','Jz(nA/m^2)'
			write(model_f(m),'(4(A9),7(A14))') 'x(R_E)','y(R_E)','z(R_E)','rad(R_E)', &
				'Bdx(nT)','Bdy(nT)','Bdz(nT)', &
				'Bdipx(nT)','Bdipy(nT)','Bdipz(nT)','|Bdip|(nT)'
		enddo

		!	Whole-box calcs
		rx=xspac(box)
		ry=xspac(box)
		rz=xspac(box)

		call calcur(bx,by,bz,nx,ny,nz,n_grids,box, &
			curx,cury,curz,rx,ry,rz)

		bxt(:,:,:) = bx0(:,:,:,box) + bx(:,:,:,box)
		byt(:,:,:) = by0(:,:,:,box) + by(:,:,:,box)
		bzt(:,:,:) = bz0(:,:,:,box) + bz(:,:,:,box)
		!
		!	Find magnitude of B
		!
		call tot_b(btot,bxt,byt,bzt,nx,ny,nz)

		call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
			opx,opy,opz,orho,curx,cury,curz, &
			evx,evy,evz,tvx,tvy,tvz, &
			nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso,reynolds)

		call bande(ex,ey,ez,bxt,byt,bzt, &
			curx,cury,curz,evx,evy,evz,btot, &
			epres,qrho,hrho,orho,resistive,resist,reynolds, &
			nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso, &
			ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
			rx,ry,rz)

		!	Find anisotropy values for each species:
		!
		!	*********
		!	Species q
		!	*********

		call fnd_vel(qpx,qpy,qpz,qrho,vvx,vvy,vvz,nx,ny,nz,n_grids,box)
		call fnd_pres(qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz,&
			qpara,qperp,qcross, &
			vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,n_grids,box)
		qvx = vvx(i,j,k)*v_equiv
		qvy = vvy(i,j,k)*v_equiv
		qvz = vvz(i,j,k)*v_equiv

		!	*********
		!	Species h
		!	*********

		call fnd_vel(hpx,hpy,hpz,hrho,vvx,vvy,vvz,nx,ny,nz,n_grids,box)
		call fnd_pres(hpresx,hpresy,hpresz,hpresxy,hpresxz,hpresyz,&
			hpara,hperp,hcross, &
			vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,n_grids,box)
		hvx = vvx(i,j,k)*v_equiv
		hvy = vvy(i,j,k)*v_equiv
		hvz = vvz(i,j,k)*v_equiv

		!	*********
		!	Species o
		!	*********

		call fnd_vel(opx,opy,opz,orho,vvx,vvy,vvz,nx,ny,nz,n_grids,box)
		call fnd_pres(opresx,opresy,opresz,opresxy,opresxz,opresyz,&
			opara,operp,ocross, &
			vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,n_grids,box)
		ovx = vvx(i,j,k)*v_equiv
		ovy = vvy(i,j,k)*v_equiv
		ovz = vvz(i,j,k)*v_equiv


		!	******************
		!	Choose slice here:
		!	******************

		!	For limit=60., choosing k=31 gives the z=0 plane.
		!	For limit=60., choosing i=61 or j=61 gives the x=0
		!		or y=0 plane, respectively, assuming the box
		!		is symmetric about the origin.

		cuts_vals(1) = int(limit/2.+1)
		cuts_vals(2) = int(limit+1)
		cuts_vals(3) = int(limit-14)	!	Box is not centered on the origin

		!	*****************
		!		Cut 1: xy
		!	*****************

		k=cuts_vals(1)
		do j=1, ny
			do i=1, nx
				unetx = tvx(i,j,k) * v_equiv
				unety = tvy(i,j,k) * v_equiv
				unetz = tvz(i,j,k) * v_equiv
				net_flow = sqrt(unetx**2 + unety**2 + unetz**2)

				curx_all = curx(i,j,k)*cur_equiv
				cury_all = cury(i,j,k)*cur_equiv
				curz_all = curz(i,j,k)*cur_equiv

				aqpres = amax1(  0.0000001, 0.333*( qpresx(i,j,k,box) +  &
					qpresy(i,j,k,box) + qpresz(i,j,k,box) )  )
				ahpres = amax1(  0.0000001, 0.333*( hpresx(i,j,k,box) + &
					hpresy(i,j,k,box) + hpresz(i,j,k,box) )  )
				aopres = amax1(  0.0000001, 0.333*( opresx(i,j,k,box) + &
					opresy(i,j,k,box) + opresz(i,j,k,box) )  )
				aepres = amax1( 0.0000001, epres(i,j,k,box) )

				lqpara = alog10( qpara(i,j,k,box) * pres_equiv )
				lqperp = alog10( qperp(i,j,k,box) * pres_equiv )
				lqcross = alog10( qcross(i,j,k,box) * pres_equiv )
				lhpara = alog10( hpara(i,j,k,box) * pres_equiv )
				lhperp = alog10( hperp(i,j,k,box) * pres_equiv )
				lhcross = alog10( hcross(i,j,k,box) * pres_equiv )
				lopara = alog10( opara(i,j,k,box) * pres_equiv )
				loperp = alog10( operp(i,j,k,box) * pres_equiv )
				locross = alog10( ocross(i,j,k,box) * pres_equiv )

				qden = amax1(qrho(i,j,k,box),0.000001)
				hden = amax1(0.0001,hrho(i,j,k,box))
				oden = amax1(0.0001,orho(i,j,k,box))

				tden = ( qrho(i,j,k,box) + hrho(i,j,k,box) + orho(i,j,k,box) ) * rho_equiv

				qden = qden / rmassq * rho_equiv
				hden = hden / rmassh * rho_equiv
				oden = oden / rmasso * rho_equiv
				eden = oden + hden + qden + 0.0001

				qtemp = amax1(0.00001,aqpres/qden) * temp_equiv
				htemp = ahpres / hden * temp_equiv
				otemp = aopres / oden * temp_equiv
				etemp = aepres / eden * temp_equiv

				qdens = alog10(qden)
				hdens = alog10(hden)
				odens = alog10(oden)
				edens = alog10(eden)

				aqpres = alog10( aqpres * pres_equiv )
				ahpres = alog10( ahpres * pres_equiv )
				aopres = alog10( aopres * pres_equiv )
				aepres = alog10( aepres * pres_equiv )

				abx = b_equiv * ( bx(i,j,k,box) + bx0(i,j,k,box) )
				aby = b_equiv * ( by(i,j,k,box) + by0(i,j,k,box) )
				abz = b_equiv * ( bz(i,j,k,box) + bz0(i,j,k,box) )

				aefldx = ex(i,j,k) * e_equiv
				aefldy = ey(i,j,k) * e_equiv
				aefldz = ez(i,j,k) * e_equiv
				emag = sqrt( aefldx**2 + aefldy**2 + aefldz**2 )

				ri = grid_minvals(1,box) + (xspac(box)*real(i-1))
				ri = ri*re_equiv
				rj = grid_minvals(2,box) + (xspac(box)*real(j-1))
				rj = rj*re_equiv
				rk = grid_minvals(3,box) + (xspac(box)*real(k-1))
				rk = rk*re_equiv

				rad = sqrt( ri**2 + rj**2 + rk**2 )

				bdipmag = sqrt( (bx0(i,j,k,box)**2) &
					+ (by0(i,j,k,box)**2) &
					+ (bz0(i,j,k,box)**2) )
				bdipmag = bdipmag*b_equiv
				bmag = sqrt( abx**2 + aby**2 + abz**2 )

				v_alf = bmag*1.e-9 / sqrt(mu0 * tden*1.e6 * m_prot)
				alfven_mach = net_flow*1.e3 / v_alf

				! Perturbed and unperturbed quantities
				bdx=b_equiv*(bx(i,j,k,box))
				ub0x=b_equiv*(bx0(i,j,k,box))
				bdy=b_equiv*(by(i,j,k,box))
				ub0y=b_equiv*(by0(i,j,k,box))
				bdz=b_equiv*(bz(i,j,k,box))
				ub0z=b_equiv*(bz0(i,j,k,box))

				!	~~~~~~~~~~~
				!	Cut 1 print
				!	~~~~~~~~~~~

				write(plas_f(1),'(8(es14.6))') &
					qdens,qtemp, hdens,htemp, &
					odens,otemp, edens,etemp

				write(flow_f(1),'(13(es14.6))') &
					qvx,qvy,qvz, hvx,hvy,hvz, ovx,ovy,ovz, &
					unetx,unety,unetz, alfven_mach

				write(pres_f(1),'(10(es14.6))') &
					lqpara,lqperp,lqcross, &
					lhpara,lhperp,lhcross, &
					lopara,loperp,locross, &
					aepres

				write(bande_f(1),'(11(es14.6))') &
					abx,aby,abz,bmag, aefldx,aefldy,aefldz,emag, &
					curx_all,cury_all,curz_all
					
				write(model_f(1),'(4(f9.2),7(es14.6))') &
					ri,rj,rk,rad, &
					bdx,bdy,bdz, ub0x,ub0y,ub0z,bdipmag
			enddo
		enddo

		!	*****************
		!		Cut 2: xz
		!	*****************

		j=cuts_vals(2)
		do k=1, nz
			do i=1, nx
				unetx = tvx(i,j,k) * v_equiv
				unety = tvy(i,j,k) * v_equiv
				unetz = tvz(i,j,k) * v_equiv
				net_flow = sqrt(unetx**2 + unety**2 + unetz**2)

				curx_all = curx(i,j,k)*cur_equiv
				cury_all = cury(i,j,k)*cur_equiv
				curz_all = curz(i,j,k)*cur_equiv

				aqpres = amax1(  0.0000001, 0.333*( qpresx(i,j,k,box) +  &
					qpresy(i,j,k,box) + qpresz(i,j,k,box) )  )
				ahpres = amax1(  0.0000001, 0.333*( hpresx(i,j,k,box) + &
					hpresy(i,j,k,box) + hpresz(i,j,k,box) )  )
				aopres = amax1(  0.0000001, 0.333*( opresx(i,j,k,box) + &
					opresy(i,j,k,box) + opresz(i,j,k,box) )  )
				aepres = amax1( 0.0000001, epres(i,j,k,box) )

				lqpara = alog10( qpara(i,j,k,box) * pres_equiv )
				lqperp = alog10( qperp(i,j,k,box) * pres_equiv )
				lqcross = alog10( qcross(i,j,k,box) * pres_equiv )
				lhpara = alog10( hpara(i,j,k,box) * pres_equiv )
				lhperp = alog10( hperp(i,j,k,box) * pres_equiv )
				lhcross = alog10( hcross(i,j,k,box) * pres_equiv )
				lopara = alog10( opara(i,j,k,box) * pres_equiv )
				loperp = alog10( operp(i,j,k,box) * pres_equiv )
				locross = alog10( ocross(i,j,k,box) * pres_equiv )

				qden = amax1(qrho(i,j,k,box),0.000001)
				hden = amax1(0.0001,hrho(i,j,k,box))
				oden = amax1(0.0001,orho(i,j,k,box))

				tden = ( qrho(i,j,k,box) + hrho(i,j,k,box) + orho(i,j,k,box) ) * rho_equiv

				qden = qden / rmassq * rho_equiv
				hden = hden / rmassh * rho_equiv
				oden = oden / rmasso * rho_equiv
				eden = oden + hden + qden + 0.0001

				qtemp = amax1(0.00001,aqpres/qden) * temp_equiv
				htemp = ahpres / hden * temp_equiv
				otemp = aopres / oden * temp_equiv
				etemp = aepres / eden * temp_equiv

				qdens = alog10(qden)
				hdens = alog10(hden)
				odens = alog10(oden)
				edens = alog10(eden)

				aqpres = alog10( aqpres * pres_equiv )
				ahpres = alog10( ahpres * pres_equiv )
				aopres = alog10( aopres * pres_equiv )
				aepres = alog10( aepres * pres_equiv )

				abx = b_equiv * ( bx(i,j,k,box) + bx0(i,j,k,box) )
				aby = b_equiv * ( by(i,j,k,box) + by0(i,j,k,box) )
				abz = b_equiv * ( bz(i,j,k,box) + bz0(i,j,k,box) )

				aefldx = ex(i,j,k) * e_equiv
				aefldy = ey(i,j,k) * e_equiv
				aefldz = ez(i,j,k) * e_equiv
				emag = sqrt( aefldx**2 + aefldy**2 + aefldz**2 )

				ri = grid_minvals(1,box) + (xspac(box)*real(i-1))
				ri = ri*re_equiv
				rj = grid_minvals(2,box) + (xspac(box)*real(j-1))
				rj = rj*re_equiv
				rk = grid_minvals(3,box) + (xspac(box)*real(k-1))
				rk = rk*re_equiv

				rad = sqrt( ri**2 + rj**2 + rk**2 )

				bdipmag = sqrt( (bx0(i,j,k,box)**2) &
					+ (by0(i,j,k,box)**2) &
					+ (bz0(i,j,k,box)**2) )
				bdipmag = bdipmag*b_equiv
				bmag = sqrt( abx**2 + aby**2 + abz**2 )

				v_alf = bmag*1.e-9 / sqrt(mu0 * tden*1.e6 * m_prot)
				alfven_mach = net_flow*1.e3 / v_alf

				! Perturbed and unperturbed quantities
				bdx=b_equiv*(bx(i,j,k,box))
				ub0x=b_equiv*(bx0(i,j,k,box))
				bdy=b_equiv*(by(i,j,k,box))
				ub0y=b_equiv*(by0(i,j,k,box))
				bdz=b_equiv*(bz(i,j,k,box))
				ub0z=b_equiv*(bz0(i,j,k,box))

				!	~~~~~~~~~~~
				!	Cut 2 print
				!	~~~~~~~~~~~

				write(plas_f(2),'(8(es14.6))') &
					qdens,qtemp, hdens,htemp, &
					odens,otemp, edens,etemp

				write(flow_f(2),'(13(es14.6))') &
					qvx,qvy,qvz, hvx,hvy,hvz, ovx,ovy,ovz, &
					unetx,unety,unetz, alfven_mach

				write(pres_f(2),'(10(es14.6))') &
					lqpara,lqperp,lqcross, &
					lhpara,lhperp,lhcross, &
					lopara,loperp,locross, &
					aepres

				write(bande_f(2),'(11(es14.6))') &
					abx,aby,abz,bmag, aefldx,aefldy,aefldz,emag, &
					curx_all,cury_all,curz_all
					
				write(model_f(2),'(4(f9.2),7(es14.6))') &
					ri,rj,rk,rad, &
					bdx,bdy,bdz, ub0x,ub0y,ub0z,bdipmag
			enddo
		enddo

		!	*****************
		!		Cut 3: yz
		!	*****************

		i=cuts_vals(3)
		do k=1, nz
			do j=1, ny
				unetx = tvx(i,j,k) * v_equiv
				unety = tvy(i,j,k) * v_equiv
				unetz = tvz(i,j,k) * v_equiv
				net_flow = sqrt(unetx**2 + unety**2 + unetz**2)

				curx_all = curx(i,j,k)*cur_equiv
				cury_all = cury(i,j,k)*cur_equiv
				curz_all = curz(i,j,k)*cur_equiv

				aqpres = amax1(  0.0000001, 0.333*( qpresx(i,j,k,box) +  &
					qpresy(i,j,k,box) + qpresz(i,j,k,box) )  )
				ahpres = amax1(  0.0000001, 0.333*( hpresx(i,j,k,box) + &
					hpresy(i,j,k,box) + hpresz(i,j,k,box) )  )
				aopres = amax1(  0.0000001, 0.333*( opresx(i,j,k,box) + &
					opresy(i,j,k,box) + opresz(i,j,k,box) )  )
				aepres = amax1( 0.0000001, epres(i,j,k,box) )

				lqpara = alog10( qpara(i,j,k,box) * pres_equiv )
				lqperp = alog10( qperp(i,j,k,box) * pres_equiv )
				lqcross = alog10( qcross(i,j,k,box) * pres_equiv )
				lhpara = alog10( hpara(i,j,k,box) * pres_equiv )
				lhperp = alog10( hperp(i,j,k,box) * pres_equiv )
				lhcross = alog10( hcross(i,j,k,box) * pres_equiv )
				lopara = alog10( opara(i,j,k,box) * pres_equiv )
				loperp = alog10( operp(i,j,k,box) * pres_equiv )
				locross = alog10( ocross(i,j,k,box) * pres_equiv )

				qden = amax1(qrho(i,j,k,box),0.000001)
				hden = amax1(0.0001,hrho(i,j,k,box))
				oden = amax1(0.0001,orho(i,j,k,box))

				tden = ( qrho(i,j,k,box) + hrho(i,j,k,box) + orho(i,j,k,box) ) * rho_equiv

				qden = qden / rmassq * rho_equiv
				hden = hden / rmassh * rho_equiv
				oden = oden / rmasso * rho_equiv
				eden = oden + hden + qden + 0.0001

				qtemp = amax1(0.00001,aqpres/qden) * temp_equiv
				htemp = ahpres / hden * temp_equiv
				otemp = aopres / oden * temp_equiv
				etemp = aepres / eden * temp_equiv

				qdens = alog10(qden)
				hdens = alog10(hden)
				odens = alog10(oden)
				edens = alog10(eden)

				aqpres = alog10( aqpres * pres_equiv )
				ahpres = alog10( ahpres * pres_equiv )
				aopres = alog10( aopres * pres_equiv )
				aepres = alog10( aepres * pres_equiv )

				abx = b_equiv * ( bx(i,j,k,box) + bx0(i,j,k,box) )
				aby = b_equiv * ( by(i,j,k,box) + by0(i,j,k,box) )
				abz = b_equiv * ( bz(i,j,k,box) + bz0(i,j,k,box) )

				aefldx = ex(i,j,k) * e_equiv
				aefldy = ey(i,j,k) * e_equiv
				aefldz = ez(i,j,k) * e_equiv
				emag = sqrt( aefldx**2 + aefldy**2 + aefldz**2 )

				ri = grid_minvals(1,box) + (xspac(box)*real(i-1))
				ri = ri*re_equiv
				rj = grid_minvals(2,box) + (xspac(box)*real(j-1))
				rj = rj*re_equiv
				rk = grid_minvals(3,box) + (xspac(box)*real(k-1))
				rk = rk*re_equiv

				rad = sqrt( ri**2 + rj**2 + rk**2 )

				bdipmag = sqrt( (bx0(i,j,k,box)**2) &
					+ (by0(i,j,k,box)**2) &
					+ (bz0(i,j,k,box)**2) )
				bdipmag = bdipmag*b_equiv
				bmag = sqrt( abx**2 + aby**2 + abz**2 )

				v_alf = bmag*1.e-9 / sqrt(mu0 * tden*1.e6 * m_prot)
				alfven_mach = net_flow*1.e3 / v_alf

				! Perturbed and unperturbed quantities
				bdx=b_equiv*(bx(i,j,k,box))
				ub0x=b_equiv*(bx0(i,j,k,box))
				bdy=b_equiv*(by(i,j,k,box))
				ub0y=b_equiv*(by0(i,j,k,box))
				bdz=b_equiv*(bz(i,j,k,box))
				ub0z=b_equiv*(bz0(i,j,k,box))

				!	~~~~~~~~~~~
				!	Cut 3 print
				!	~~~~~~~~~~~

				write(plas_f(3),'(8(es14.6))') &
					qdens,qtemp, hdens,htemp, &
					odens,otemp, edens,etemp

				write(flow_f(3),'(13(es14.6))') &
					qvx,qvy,qvz, hvx,hvy,hvz, ovx,ovy,ovz, &
					unetx,unety,unetz, alfven_mach

				write(pres_f(3),'(10(es14.6))') &
					lqpara,lqperp,lqcross, &
					lhpara,lhperp,lhcross, &
					lopara,loperp,locross, &
					aepres

				write(bande_f(3),'(11(es14.6))') &
					abx,aby,abz,bmag, aefldx,aefldy,aefldz,emag, &
					curx_all,cury_all,curz_all
					
				write(model_f(3),'(4(f9.2),7(es14.6))') &
					ri,rj,rk,rad, &
					bdx,bdy,bdz, ub0x,ub0y,ub0z,bdipmag
			enddo
		enddo

		do m=1, n_cuts
			close(plas_f(m))
			close(flow_f(m))
			close(pres_f(m))
			close(bande_f(m))
			close(model_f(m))
		enddo

	enddo ! loop over box

	write(*,*) "Done."

!	call system("python3 "//trim(python_dir)//trim(python_plotter))

	if(mod(nplots,10) .eq. 0) then
		! MJS: Make gifs out of time sequences here (and overwrite)
		write(*,*) "gifs updated."
	endif

	return
end subroutine write_graphing_data
