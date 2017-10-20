!
!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!		@											@
!		@		MULTIFLUID SIM MAIN SEQUENCE		@
!		@											@
!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!   This is a 3-d modified three fluid simulation using
!         electrons  : arrays starting with e
!         solar wind : arrays starting with q (protons)
!         ionospheric: arrays starting with o oxygen,
!                                           h hydrogen
!
!   ?[ and change set33(0.1,1.,0.1,1.0  to
!               set3(0.,1.,0.,1.
!   ]?
!	Warnings: make sure space is compatible with graphics
!               routines
!             the arrays ijzero,ijmid  and ijsrf have to
!               be modified in bsubs.f if modified in main
!             add in mirror dipole so bx is zero at wind boundary
!             graphics - contour and flows - have to be manually set
!               for right aspect ratios
!             plasma and magnetic field data must be aligned in time
!
!	Grid within grid size nt = division*n_grids
!	ncts is the size of the data array for the IMF data file
!
!	File indices:
!		1-29	Top-level I/O, fluid files, primary output data, etc
!		30-59	Spacecraft position/trajectory input files
!		60-89	Spacecraft data recording
!
program multifluid
	!
	!	***********************
	!			Modules
	!	***********************
	!
	!	Contains planetary constants
	!	And some physical constants, including pi
	use astrometry
	!
	!	*******************
	!	Critical parameters
	!	*******************
		integer,parameter :: nx=121,ny=121,nz=61,n_grids=5,division=2, &
		mbndry=1, ncts=281, num_pts(3)=[nx,ny,nz]
		integer, parameter :: dp = kind(1.d0)  !	Double precision
		!
		real,parameter	:: d_min=0.01	! Minimum density
		real,parameter	:: o_conc_min=0.01
		real,parameter	:: o_conc_max=1.0
		real,parameter	:: cur_min=0.75
		real,parameter	:: cur_max=20.0
		!
	    !	Scale lengths for plasma sheet population
		real,parameter	:: sheet_den=0.25
		real,parameter	:: alpha_s=4.
		!
		integer,parameter :: num_quants = 7	! Number of quantities to zero inside inner boundary
		integer mzero	!	Number of grid points to zero inside innermost boundary
		integer mmid	!	Number of grid points in 'mid' region between inner boundary and zero region
		integer msrf	!	Number of grid points on inner boundary surface
		real,parameter	:: zero_bndry_offset = -1.5
		real,parameter	:: mid_bndry_offset = -0.5
		real,parameter	:: srf_bndry_offset = 0.6
		real zero_bndry
		real mid_bndry
		real srf_bndry
		integer,allocatable	:: ijsrf(:,:,:)
		integer,allocatable	:: ijmid(:,:,:)
		integer,allocatable	:: ijzero(:,:,:)
		integer	numsrf(mbndry)
		integer	nummid(mbndry)
		integer	numzero(mbndry)
		real,allocatable	:: parm_srf(:,:,:)
		real,allocatable	:: parm_mid(:,:,:)
		real,allocatable	:: parm_zero(:,:,:)
		!
	!
	!	********************
	!	Graphics parameters:
	!	********************
		integer,parameter :: mx=61
		integer,parameter :: my=61
		integer,parameter :: mz=31
		integer,parameter :: muvwp2=63
		integer,parameter :: mz2=16
		!	Note: muvwp2=amax(mx,my,mz)+2,mz2=(mz-1)/2+1
		!
	!
	!	*********************
	!	Coordinate variables:
	!	*********************
		real dx,dy,dz		!	Used for differences between points in various checks (stand-in for xspac)
		real ax,ay,az		!	Used in loops over grid points--x,y,z locations in re_equiv relative to origin
		real xr,yr,zr		!	Rotated coordinates
		real xp,yp,zp		!	Dipole coordinates
		real ar, ar2		!	Radial distance combining ax, ay, az
		real ra				!	
		real avx,avy,avz	!	Active loop values for velocity components
		real rx, ry, rz, rd	!	Distances relative to the origin
		real rvx,rvy,rvz	!	Velocity variables derived from rotational motion
		!
		real corotate		!	Corotation speed for current distance
		real ar_tmid		!	Radial distance from origin to torus center in sim units
		real dr_tmid		!	Distance from center of torus to current point
		real ar_iono
		real r_equat
		real rerg_sphere
		real rden_sphere
		real zheight
		real rho_iono
		real vt
		real abtot
		!
		integer nrot_p, nrot_m	!	Number of full rotations for planet and moon
		real rot_hrs, rot_hrs_m	!	Number of hours into current local day/lunar orbit
		real rot_angle, rot_angle_m
		real sin_rot, cos_rot, sin_rot_m, cos_rot_m
		!
    !
	!
	!	*********************
	!	Input file parameters
	!	*********************
		! group 'option':
		real tmax, stepsz, tsave
		integer ntgraph, ntinj
		logical start, isotropic
		! group 'planet'
		character*10 bodyname, moonname
		real xdip, ydip, zdip, r_inner, torus_rad, &
		tilt1, tilt2, &
		rmassq, rmassh, rmasso
		logical tilting
		! group 'speeds'
		real cs_inner, alf_inner1, alf_inner2, &
		alpha_e, denh_inner, denq_torus, denh_torus, deno_torus, &
		gravity, ti_te, gamma, reduct, t_torus, aniso_factor, &
		ani_q, ani_h, ani_o
		logical ringo, update, reload, divb_lores, divb_hires
		! group 'windy'
		real re_wind, cs_wind, &
		vx_wind1, vx_wind2, vy_wind1, vy_wind2, vz_wind1, vz_wind2, &
		alfx_wind1, alfx_wind2, alfy_wind1, alfy_wind2, alfz_wind1, alfz_wind2, &
		den_wind1, den_wind2, &
		reynolds, resist, o_conc, rho_frac, bfrac, vfrac
		! group 'lunar'
		real orbit_moon, theta_moon, cs_moon, &
		qden_moon, hden_moon, oden_moon, &
		alf_moon, ti_te_moon, &
		xdip_moon, ydip_moon, zdip_moon, offset
		! group 'physical' 
		real re_equiv, b_equiv, v_equiv, rho_equiv
		!real(dp) utstart	!	Some day.
		real utstart
		logical spacecraft, craft_input, warp
		! group 'smooth'
		real chirho, chipxyz, chierg, difrho, difpxyz, diferg
		!
	!	**********************
	!	Output file parameters
	!	**********************
		! group 'crafthead'
		character*8	:: cname = 'default'
		integer		:: num_vals = 0	!	FIX NEEDED: Write the ntimes value for each craft to this
		real		:: rot_closest = 0.0	!	FIX NEEDED: Use utc_to_jd subroutine to find these numbers so they can print to file headers
		character*40 :: git_hash = 'deadbeef'
		!
		real,parameter :: wind_adjust=4./3., limit=60.
		real xspac(n_grids), grid_minvals(3,n_grids), grid_maxvals(3,n_grids)
!		real(dp) ut	!	Some day.
		real ut
		!
		!	Notes about grid variables:
		!		
		!	Grid limits are set by grid_minvals grid_maxvals arrays. 'limit' is the
		!		base min/max value for the smallest xy box. z is half as much.
		!	xspac is the relative grid spacing; the distance between grid
		!		points in units of the smallest grid step size.
		!	wind_adjust is the factor to stretch/compress the biggest box
		!		along the wind direction to capture the magnetotail.
		!	UT is in hours; we use higher precision so that we can
		!		track spacecraft times seconds apart after many hours of simulation time.
		!	UT = 0 is chosen to keep values relatively low for Europa sims, 1996-12-02T14:24:00.000
		!
		!
	!	*******************
	!	Spacecraft file I/O
	!	*******************
		!	File indices:
		!		1-29	Top-level I/O, fluid files, primary output data, etc
		!		30-59	Spacecraft position/trajectory input files
		!		60-89	Spacecraft data recording
		!		100		Git hash file
		!		101		Dummy craft file
		integer,parameter :: scin=30, scout=60, git_f=100, dummy_f=101
		character,parameter :: tab=char(9)
		character*120,parameter :: dat_header='time'//tab//'xpos'//tab// &
		'ypos'//tab//'zpos'//tab//'Bxval'//tab//'Byval'//tab//'Bzval'//tab//'temp'
		character*80,parameter :: flux_header='main grid'//tab//'UT'//tab// &
		'box'//tab//'qflux'//tab//'hflux'//tab//'oflux'
		character*12,parameter :: git_hash_file='git_hash.txt'
		!
	!
    !	************************************
    !	Physics plasma quantities: Main grid
    !	************************************
		real bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
		qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
		qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
		qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
		qpresxy(nx,ny,nz,n_grids),qpresxz(nx,ny,nz,n_grids), &
		qpresyz(nx,ny,nz,n_grids), &
		!
		hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
		hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
		hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
		hpresxy(nx,ny,nz,n_grids),hpresxz(nx,ny,nz,n_grids), &
		hpresyz(nx,ny,nz,n_grids), &
		!
		opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
		orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
		opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
		opresxy(nx,ny,nz,n_grids),opresxz(nx,ny,nz,n_grids), &
		opresyz(nx,ny,nz,n_grids), &
		!
		epres(nx,ny,nz,n_grids)
		!
	!
	!	****************************************************
    !	Work arrays for Runge-Kutta and smoothing: Main grid
	!	****************************************************
		real, allocatable, dimension(:,:,:,:) :: &
		oldbx,oldby,oldbz, &
		oldqrho,oldqpx,oldqpy,oldqpz, &
		oldqpresx,oldqpresy,oldqpresz, &
		oldqpresxy,oldqpresxz,oldqpresyz, &
		oldhrho,oldhpx,oldhpy,oldhpz, &
		oldhpresx,oldhpresy,oldhpresz, &
		oldhpresxy,oldhpresxz,oldhpresyz, &
		oldorho,oldopx,oldopy,oldopz, &
		oldopresx,oldopresy,oldopresz, &
		oldopresxy,oldopresxz,oldopresyz, &
		oldepres
		!
		real, allocatable, dimension(:,:,:,:) :: &
		wrkbx,wrkby,wrkbz, &
		wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
		wrkqpresx,wrkqpresy,wrkqpresz, &
		wrkqpresxy,wrkqpresxz,wrkqpresyz, &
		wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
		wrkhpresx,wrkhpresy,wrkhpresz, &
		wrkhpresxy,wrkhpresxz,wrkhpresyz, &
		wrkorho,wrkopx,wrkopy,wrkopz, &
		wrkopresx,wrkopresy,wrkopresz, &
		wrkopresxy,wrkopresxz,wrkopresyz, &
		wrkepres
		!
	!
	!	**********************
    !	Unperturbed quantities
	!	**********************
		real bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids),bz0(nx,ny,nz,n_grids)
		!
		real efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
		curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
		bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz), &
		resistive(nx,ny,nz,mbndry)
		!
	!
	!	*************************
    !	Variable time step arrays
	!	*************************
		!
!		real(dp) t_old(n_grids)	!	Some day we will do this. But not today.
!		real(dp) t_new(n_grids)
!		real(dp) t_step(n_grids)
!		real(dp) t_stepnew(n_grids)
!		real(dp) t_equiv
!		real(dp) t	!	We give t extra precision to make sure we don't easily hit an upper limit on runtime.
		real t_old(n_grids)
		real t_new(n_grids)
		real t_step(n_grids)
		real t_stepnew(n_grids)
		real t_equiv
		real t
		!
	!
	!	*************************
    !	Boundary condition arrays
	!	*************************
		real bxf(ny,nz), byf(ny,nz), bzf(ny,nz), &
		rhof(ny,nz), svxf(ny,nz), svyf(ny,nz), svzf(ny,nz), &
		bxp(ny,nz), byp(ny,nz), bzp(ny,nz), &
		rhop(ny,nz), svxp(ny,nz), svyp(ny,nz), svzp(ny,nz), &
		future(ny,nz), past(ny,nz), &
		bfld(ncts,4), rplas(ncts), svel(ncts,3)
		integer ncount(ny,nz)
		!
		real tx(mx,my,mz), ty(mx,my,mz), tz(mx,my,mz), tg1(mx,my,mz), &
		tg2(mx,my,mz2), tt(mx,my,mz), work(muvwp2,muvwp2), &
		cross(ny,nz), along(nx,nz), flat(nx,ny)
		!
	!
	!	***************************
	!	Labels for graphics outputs
	!	***************************
		character*8 wd1,wd2,wd3,wd4
		character*8 label
		character*15 title
		!
	!
	!	*************
	!	Miscellaneous
	!	*************
		integer n, m, ms, nn, box	!	Loop counters
		integer m_smallest, m_step, cbox	!	Placeholders
		integer mating_index, next_box	!	For checking grid mating
		real grid_diff	!	For checking grid mating
		!
		logical add_dip, save_dat, write_dat, &
		yes_step, yes_step_n, grid_reset
		!
	!
	!	************************
	!	Adjusted data parameters
	!	************************
		!
		real gamma1
		real epress, eerg, spress, serg
		real rho_wind1, rho_wind2
		real srho, delrho
		real svelx, svely, svelz
		real spx, spy, spz
		!
		real delvx_wind, delvy_wind, delvz_wind
		real delbx_wind, delby_wind, delbz_wind
		real sbx_wind, sby_wind, sbz_wind
		!
		real deltg, deltinj, delt
		!
		real	:: tgraph = 0.
		real	:: tinj = 0.
		real	:: tdiv = 10.
		real	:: del_tdiv = 5.
	!
	!		spacecraft decides whether to include any spacecraft at all
	!			(input or output)
	!		craft_input must be true to use input data from spacecraft
	!			(craft_input is not currently implemented)
	!			(spacecraft must be set to true to use craft_input)
    !
    !      ringo decides if you you want to plot 1 set of diagnostics
    !              with no time stepping
    !      update decides if time setting is to be reset
    !
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
    !           at n+1/2 as defined by the lax-wendroff scheme
    !       the index of 3 is needed to store adjacent x values for these fns
    !
    !     d_min is the minimum allowable density
    !     stepsz is the size of the time step in terms of delta_x,y/(|umax|+c_s)
    !
    !     system dimensions are nx,ny,nz
    !     variable grid spacing enabled with rxyz >1
    !           rx,ry,rz should not be set identically to zero
    !
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	@		SPACECRAFT INITIALIZATION		@
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!
	!       xcraft is the actual position of the spacecraft in re
    !       4th dimension is UT in hours
    !       zcraft is the future position of the spacecraft in re
    !       rcraft is the position of the spacecraft for which
    !           IMF is reference. No alteration from boundary conditions applied,
	!			and no data recorded. Not included in ncraft or ndef_craft.
	!		craft_info is directory for spacecraft position/times (relative to multifluid directory)
	!		craft_data is directory for output data (relative to multifluid directory)
	!
	!		craftpos contains all (x,y,z,t) for all aux craft 
	!		sxyz is linearly interpolated spacecraft xyz for a given UT.
	!		gridpts holds two vertices of the cube containing the 8 grid coordinates closest to the spacecraft for a given UT
	!		meas_qty(i) are scalar physical quantities to be measured by a spacecraft
	!		scdata(i) are the interpolated values corresponding to meas_qty(i) that a spacecraft actually records
	!		num_inst is how many 'instruments' the spacecraft has, i.e., how many scalar quantities are output to .dat
	!			When new measurements are to be added to spacecraft .dat files, num_inst must be incremented and
	!			'SPACECRAFT DATA RECORDING' section also needs a new line for the new physical quantity: scdata(i) = physical_qty.
	!			This must be done in two places, one for default craft and one for aux craft.
	!
	!		craftstat holds IOSTAT for reading in .craft files. Goes negative if EOF is reached before read() is done reading values
	!		deflt_rec_skp is the number of time stepping loops to skip between default craft measurements
	!		nskipped is the number of time stepping loops during which we skipped default craft recording
	!		ntimes(:,1) holds the number of measurements an aux craft has made
	!		ntimes(:,2) holds the number of (x,y,z,t) points we have for each aux craft
	!		craft_gridpt(1:3) holds the xyz grid indices of the nearest grid point for a given craft. Default craft must be placed exactly on a grid point.
	!			craft_gridpt(4) holds the box number for the smallest box our craft.
	!		fname is a path to file, relative to the multifluid directory
	!		recording is a flag, true by default, which is set to false when a spacecraft has recorded data for its entire trajectory
	!					
	!		cname, num_vals, and rot_closest are craft name, number of measurements, and rot_hrs at closest approach. These values are read from .craft files via a namelist and output as header information in spacecraft data files.
	!		Spacecraft data recording has a variable number of header lines. A dummy_craft file is created to print a header and count the number of lines to skip when seeking to the end of these files. The number of header lines is stored in nheadlines.
    !
	character*32, parameter	:: craft_info='spacecraft_info/'
	character*32, parameter	:: craft_data='data/spacecraft_data/'
	character*11, parameter :: dummy_craft='dummy.craft'
	character*120	junkline
	character*120	fname
	character*8		numstring
	!
	integer, parameter	:: deflt_rec_skp=20
	integer, parameter	:: num_inst=3	! (Bx, By, Bz)
	integer	:: craftstat=0
	integer	:: nskipped=0
	integer	ncraft
	integer	ndef_craft
	integer	naux_craft
	integer nheadlines
	!
	logical	:: dat_exists=.false.
	logical :: dummy_exists=.false.
	!
	real	meas_qty(2,2,2,num_inst)
	real	scdata(num_inst)
	real	gridpts(3,2)
	real	rcraft(3)
	real	sxyz(3)
	!
	character*8, allocatable	:: craftnames(:)
	integer, allocatable		:: craft_gridpt(:,:)
	integer, allocatable		:: ntimes(:,:)
	logical, allocatable		:: recording(:)
    real, allocatable			:: xcraft(:,:)
    real, allocatable			:: zcraft(:,:)
	real, allocatable			:: craftpos(:,:,:)
	!
	!	**********************
	!	Namelists for file I/O
	!	**********************
		namelist/option/tmax,ntgraph,stepsz,start,tsave,ntinj,isotropic
		namelist/planet/bodyname,moonname,xdip,ydip,zdip,r_inner,torus_rad, &
		tilt1,tilt2,tilting,rmassq,rmassh,rmasso
		namelist/speeds/cs_inner,alf_inner1,alf_inner2, &
		alpha_e,denh_inner,denq_torus,denh_torus,deno_torus, &
		gravity,ti_te,gamma,ringo,update,reload, &
		divb_lores,divb_hires,reduct,t_torus,aniso_factor, &
		ani_q,ani_h,ani_o
		namelist/windy/re_wind,cs_wind,vx_wind1,vx_wind2, &
		vy_wind1,vy_wind2,vz_wind1,vz_wind2, &
		alfx_wind1,alfx_wind2, &
		alfy_wind1,alfy_wind2, &
		alfz_wind1,alfz_wind2, &
		den_wind1,den_wind2, &
		reynolds,resist,o_conc,rho_frac,bfrac,vfrac
		namelist/lunar/orbit_moon,theta_moon,cs_moon, &
		qden_moon,hden_moon,oden_moon, &
		alf_moon,ti_te_moon, &
		xdip_moon,ydip_moon,zdip_moon,offset
		namelist/physical/re_equiv,b_equiv,v_equiv,rho_equiv, &
		spacecraft,craft_input,warp,utstart
		namelist/smooth/chirho,chipxyz,chierg, &
		difrho,difpxyz,diferg
		!
		namelist/crafthead/cname,num_vals,rot_closest,git_hash
	!
	!
	!	************************
	!	Planet & moon parameters
	!	************************
	!
	!	Notes:	
	!			Calculations continued after input file is read in
		!
		!	Adapted parameters for simulation use (boundaries etc., see 'Planet & moon calculations')
		!	torus_inj_dist is used for making a correction to plasma torus injection
			real lunar_rad, torus_inj_dist, grav	
			real r_orbit, v_orbit, tilt
			real xmoon, ymoon, zmoon, rmoon, b0_moon
			!
		!
	!	*************
	!	Common blocks
	!	*************
		real	vvx(nx,ny,nz), vvy(nx,ny,nz), vvz(nx,ny,nz), &
				tvx(nx,ny,nz), tvy(nx,ny,nz), tvz(nx,ny,nz), &
				evx(nx,ny,nz), evy(nx,ny,nz), evz(nx,ny,nz)
		real	v_rot, r_rot, sin_tilt, cos_tilt, b0
		real	planet_orbit_rad, planet_year, planet_rad, &
				planet_per, planet_mass, planet_obliq, planet_incl, &
				r_lim, torus_infall, planet_tilt, planet_init_long, &
				planet_xdip, planet_ydip, planet_zdip, torus_dist, &
				moon_orbit_rad, moon_per, moon_rad, moon_mass, &
				moon_incl, moon_init_rot
		!
		common /space/vvx,vvy,vvz, tvx,tvy,tvz, evx,evy,evz
		!
		common /rotation/v_rot,r_rot,rot_angle, xdip,ydip,zdip, &
		sin_tilt,cos_tilt,b0
		!
		common /planetary/planet_orbit_rad, planet_year, planet_rad, &
		planet_per, planet_mass, planet_obliq, planet_incl, &
		r_lim, torus_infall, planet_tilt, planet_init_long, &
		planet_xdip, planet_ydip, planet_zdip, torus_dist, &
		moon_orbit_rad, moon_per, moon_rad, moon_mass, moon_incl, &
		moon_init_rot
	!
	!
	!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!							EXECUTION SECTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	write(*,*) '-'
	call system('date +"START: %T-%Y-%m-%d"')
	write(*,*) '-'
	!
	!	*************
	!	Grab git hash
	!	*************
		!	
		call system ('git rev-parse --short HEAD > '//trim(git_hash_file))
		open(git_f,file=trim(git_hash_file),status='unknown',form='formatted')
			read(git_f,*) git_hash
			write(*,*) 'Git hash: ', git_hash
		close(git_f)
		!
	!
	write(*,*) '-'
	!
	!	******************
	!	Count header lines
	!	******************
		!
		inquire(file=trim(dummy_craft), exist=dummy_exists)
		if(dummy_exists) call system ('rm '//trim(dummy_craft))
		!
		open(dummy_f,file=trim(dummy_craft),status='unknown',form='formatted')
			write(dummy_f,crafthead)
			write(dummy_f,*) dat_header
			call system ('wc -l '//trim(dummy_craft)//' >> '//trim(dummy_craft))
			read(dummy_f,*) nheadlines, junkline
		close(dummy_f)
		call system ('rm '//trim(dummy_craft))
		!
	!
	!	********************************
    !	Open input and output data files
	!	********************************
    !
    open(1,file='input',status='old',form='formatted')
    open(2,file='input_out',status='unknown',form='formatted')
    open(3,file='fluxes.dat',status='unknown',form='formatted')
    open(6,file='speeds.dat',status='unknown',form='formatted')
    !open(7,file='grid.dat',status='unknown',form='formatted')
    !open(8,file='cur.dat',status='unknown',form='unformatted')
    !open(9,file='pot.dat',status='unknown',form='unformatted')
    open(10,file='conc.dat',status='unknown',form='formatted')
    !
	!	*****************
    !	Open ncargraphics
	!	*****************
    !
    call opngks
    call gsclip(0)
    !
	!	***************
    !	Set color table
	!	***************
    !
    call cpclrs
    !
    call gselnt(0)
    call gsplci(1)
    call wtstr(.4,.975,'3d mutant code',2,0,0)
    !
	!	************************************************************
    !	Read input parameters. We write to output at end of program.
	!	************************************************************
    !
    read(1,option)
    read(1,planet)
    read(1,speeds)
    read(1,windy)
    read(1,lunar)
    read(1,physical)
    read(1,smooth)
	!
	!	***********************************
    !	Write input data to stdout/log file
	!	***********************************
    !
    write(*,option)
    write(*,planet)
    write(*,speeds)
    write(*,windy)
    write(*,lunar)
    write(*,physical)
    write(*,smooth)
	!
	write(*,*) 'Grid xmin: ','Grid xmax: ','Grid ymin: ','Grid ymax: ', &
		'Grid zmin: ','Grid zmax: ','Grid spacing: '
    do i=1,n_grids
		xspac(i) = division**(i-1)
		grid_minvals(1,i) = -limit*xspac(i)
		grid_maxvals(1,i) = limit*xspac(i)
        grid_minvals(2,i) = -limit*xspac(i)
		grid_maxvals(2,i) = limit*xspac(i)        
		grid_minvals(3,i) = -limit/2*xspac(i)
		grid_maxvals(3,i) = limit/2*xspac(i)
		!
		if(i.eq.n_grids) then
			grid_diff = grid_maxvals(1,i) * wind_adjust - grid_maxvals(1,i)
			grid_maxvals(1,i) = grid_maxvals(1,i) + grid_diff
			grid_minvals(1,i) = grid_minvals(1,i) + grid_diff
		endif
		!
		write(*,*) grid_minvals(1,i), grid_maxvals(1,i), &
		grid_minvals(2,i), grid_maxvals(2,i), &
		grid_minvals(3,i), grid_maxvals(3,i), xspac(i)
		ix = 1 + ( grid_maxvals(1,i) - grid_minvals(1,i) ) / xspac(i)
		iy = 1 + ( grid_maxvals(2,i) - grid_minvals(2,i) ) / xspac(i)
        iz = 1 + ( grid_maxvals(3,i) - grid_minvals(3,i) ) / xspac(i)
        if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz)) then
            write(*,*)'ERROR: grid spacing does not match number of points.', &
				' Actual x,y,z: ',ix,iy,iz,'Expected x,y,z: ',nx,ny,nz
            stop
        endif
    enddo
    !
	!	**********************************************
    !	Check if subgrids mate properly to larger grid
	!	**********************************************
    !
    do box=1,n_grids-1
        next_box = box + 1
		!
        grid_diff = 1. + (grid_minvals(1,box) - grid_minvals(1,next_box)) / xspac(next_box)
        mating_index = int(grid_diff)
        dx = grid_diff - mating_index
        if(abs(dx).gt.0.001) then
            write(*,*)'Main grid: xmin dont match. ',box,next_box,mating_index,grid_diff
            stop
        endif
        !
        grid_diff = 1. + (grid_minvals(2,box) - grid_minvals(2,next_box)) / xspac(next_box)
        mating_index = int(grid_diff)
        dy = grid_diff - mating_index
        if(abs(dy).gt.0.001) then
            write(*,*)'Main grid: ymin dont match. ',box,next_box,mating_index,grid_diff
            stop
        endif
        !
        grid_diff = 1. + (grid_minvals(3,box) - grid_minvals(3,next_box)) / xspac(next_box)
        mating_index = int(grid_diff)
        dz = grid_diff - mating_index
        if(abs(dz).gt.0.001) then
            write(*,*)'Main grid: zmin dont match. ',box,next_box,mating_index,grid_diff
            stop
        endif
        !
    enddo
	!
	!	**************************
	!	Set inner boundary regions
	!	**************************
		!
		zero_bndry = r_inner + zero_bndry_offset
		mid_bndry = r_inner + mid_bndry_offset
		srf_bndry = r_inner + srf_bndry_offset
		!
	!	**************************
	!	Planet & moon calculations
	!	**************************
		!
		!	Use astrometry module to set planetary parameters. Import list:
		!	planet_orbit_rad,	planet_year,
		!	planet_rad,			planet_mass,
		!	planet_obliq,		planet_incl,
		!	r_lim,				torus_infall,
		!	planet_tilt,		planet_init_long,
		!	planet_xdip,		planet_ydip,
		!	planet_zdip,		torus_dist,
		!	moon_orbit_rad,		moon_per,
		!	moon_rad,			moon_mass,
		!	moon_incl,			moon_init_rot,
		call choose_system(bodyname,moonname)
		!
		r_rot = r_lim
		v_rot = 2.*pi*planet_rad/(planet_per*3600.)/v_equiv  !	Normalized rotation rate
		torus_inj_dist = torus_dist + torus_infall
		lunar_rad = 1.25 * moon_rad	!	Start exobase at 1.25 rt
		!
		rmoon = ( lunar_rad / planet_rad ) / re_equiv   !	In grid points
		r_orbit = moon_orbit_rad / planet_rad / re_equiv	!	In grid pts (Never used, 10/08/2017 MJS)
		v_orbit = (moon_orbit_rad*2.*pi)/(moon_per*3600.)/v_equiv    !	Sim units
		!
		tilt = planet_tilt
		sin_tilt = sin(tilt)
		cos_tilt = cos(tilt)
		dtilt = ( tilt2 - tilt1 ) / tmax
		if(tilting) then
			xdip = planet_xdip
			ydip = planet_ydip
			zdip = planet_zdip
		else
			xdip = 0.00001
			ydip = 0.
			zdip = 0.
		endif
		!
		!   Gravity in m/s**2 at planet's surface 
		!	Need t_equiv in normalization
		t_equiv = planet_rad * re_equiv / v_equiv
		grav = gravity * (planet_rad*re_equiv*1000.) / (1000.*v_equiv)**2
		!
	!	******************************
	!	Calculate inner boundary sizes
	!	******************************
		!
		!	Calculates approx. number of grid point volume units in each region
		msrf  = 1.2 * 1.33*pi*(srf_bndry**3 - mid_bndry**3)/xspac(mbndry)**3
		mmid  = 1.2 * 1.33*pi*(mid_bndry**3 - zero_bndry**3)/xspac(mbndry)**3
		mzero = 1.2 * 1.33*pi*(zero_bndry**3)/xspac(mbndry)**3
		!
	!
	!
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	@				ALLOCATE ARRAYS					@
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!
	!
		!
		!	Numbers of grid points to zero inside inner boundary
		!
		allocate(ijsrf(mbndry,3,msrf), ijmid(mbndry,3,mmid), &
		ijzero(mbndry,3,mzero), &
		!
		parm_srf(mbndry,num_quants,msrf), &
		parm_mid(mbndry,num_quants,mmid), &
		parm_zero(mbndry,num_quants,mzero) )
		!
		!	Work arrays for runge-kutta and smoothing
		!
		allocate(oldbx(nx,ny,nz,n_grids),oldby(nx,ny,nz,n_grids), &
		oldbz(nx,ny,nz,n_grids), &
		!
		oldqpx(nx,ny,nz,n_grids),oldqpy(nx,ny,nz,n_grids), &
		oldqpz(nx,ny,nz,n_grids),oldqrho(nx,ny,nz,n_grids), &
		oldqpresx(nx,ny,nz,n_grids),oldqpresy(nx,ny,nz,n_grids), &
		oldqpresz(nx,ny,nz,n_grids),oldqpresxy(nx,ny,nz,n_grids), &
		oldqpresxz(nx,ny,nz,n_grids),oldqpresyz(nx,ny,nz,n_grids), &
		!
		oldhpx(nx,ny,nz,n_grids),oldhpy(nx,ny,nz,n_grids), &
		oldhpz(nx,ny,nz,n_grids),oldhrho(nx,ny,nz,n_grids), &
		oldhpresx(nx,ny,nz,n_grids),oldhpresy(nx,ny,nz,n_grids), &
		oldhpresz(nx,ny,nz,n_grids),oldhpresxy(nx,ny,nz,n_grids), &
		oldhpresxz(nx,ny,nz,n_grids),oldhpresyz(nx,ny,nz,n_grids), &
		!
		oldopx(nx,ny,nz,n_grids),oldopy(nx,ny,nz,n_grids), &
		oldopz(nx,ny,nz,n_grids),oldorho(nx,ny,nz,n_grids), &
		oldopresx(nx,ny,nz,n_grids),oldopresy(nx,ny,nz,n_grids), &
		oldopresz(nx,ny,nz,n_grids),oldopresxy(nx,ny,nz,n_grids), &
		oldopresxz(nx,ny,nz,n_grids),oldopresyz(nx,ny,nz,n_grids), &
		!
		oldepres(nx,ny,nz,n_grids))
		!
		allocate(wrkbx(nx,ny,nz,n_grids),wrkby(nx,ny,nz,n_grids), &
		wrkbz(nx,ny,nz,n_grids), &
		!
		wrkqpx(nx,ny,nz,n_grids),wrkqpy(nx,ny,nz,n_grids), &
		wrkqpz(nx,ny,nz,n_grids),wrkqrho(nx,ny,nz,n_grids), &
		wrkqpresx(nx,ny,nz,n_grids),wrkqpresy(nx,ny,nz,n_grids), &
		wrkqpresz(nx,ny,nz,n_grids),wrkqpresxy(nx,ny,nz,n_grids), &
		wrkqpresxz(nx,ny,nz,n_grids),wrkqpresyz(nx,ny,nz,n_grids), &
		!
		wrkhpx(nx,ny,nz,n_grids),wrkhpy(nx,ny,nz,n_grids), &
		wrkhpz(nx,ny,nz,n_grids),wrkhrho(nx,ny,nz,n_grids), &
		wrkhpresx(nx,ny,nz,n_grids),wrkhpresy(nx,ny,nz,n_grids), &
		wrkhpresz(nx,ny,nz,n_grids),wrkhpresxy(nx,ny,nz,n_grids), &
		wrkhpresxz(nx,ny,nz,n_grids),wrkhpresyz(nx,ny,nz,n_grids), &
		!
		wrkopx(nx,ny,nz,n_grids),wrkopy(nx,ny,nz,n_grids), &
		wrkopz(nx,ny,nz,n_grids),wrkorho(nx,ny,nz,n_grids), &
		wrkopresx(nx,ny,nz,n_grids),wrkopresy(nx,ny,nz,n_grids), &
		wrkopresz(nx,ny,nz,n_grids),wrkopresxy(nx,ny,nz,n_grids), &
		wrkopresxz(nx,ny,nz,n_grids),wrkopresyz(nx,ny,nz,n_grids), &
		!
		wrkepres(nx,ny,nz,n_grids))
		!
	!
	!
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	@		PARAMETER DECLARATIONS COMPLETE			@
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!
	!
	!	*************************************
    !	Write important data to graphics file
	!	*************************************
		write(wd4,'(f6.3)')stepsz
		!
		title='stepsz = '//wd4
		call wtstr(.75,.85,title,1,0,0)
		!
		write(wd1,'(f6.3)')xdip
		write(wd2,'(f6.3)')ydip
		write(wd3,'(f6.3)')zdip
		!
		title='xdip = '//wd1
		call wtstr(.15,.82,title,1,0,0)
		title='ydip = '//wd2
		call wtstr(.35,.82,title,1,0,0)
		title='zdip = '//wd3
		call wtstr(.55,.82,title,1,0,0)
		!
		write(wd1,'(f5.1)')r_inner
		!
		title='r_inner = '//wd1
		call wtstr(.15,.79,title,1,0,0)
		!
		!	Calculate effective magnetic field strength
		!
		erho = denh_inner * rmassh
		b01 = alf_inner1 * sqrt(erho) * r_inner**3
		b02 = alf_inner2 * sqrt(erho) * r_inner**3
		alf_lim = 6.00 * alf_inner1
		b0 = b01
		write(*,*)'b0: ',b0
		bmax = b02
		delb0 = (b02 - b01) / tmax	! Never used as of 08/27/17. MJS
		!
		write(wd1,'(f6.3)')cs_inner
		write(wd2,'(f6.3)')alf_inner1
		write(wd3,'(f6.3)')alf_inner2
		write(wd4,'(f5.1)')denh_inner
		!
		title='cs_inner = '//wd1
		call wtstr(.15,.76,title,1,0,0)
		title='alf_inner1= '//wd2
		call wtstr(.35,.76,title,1,0,0)
		title='alf_inner2= '//wd3
		call wtstr(.55,.76,title,1,0,0)
		title='denh_inner = '//wd4
		call wtstr(.75,.76,title,1,0,0)
		!
		write(wd1,'(f6.3)')o_conc
		write(wd2,'(f6.3)')gravity
		write(wd3,'(f6.3)')rmasso
		!
		title='o_conc = '//wd1
		call wtstr(.15,.73,title,1,0,0)
		title='gravity= '//wd2
		call wtstr(.35,.73,title,1,0,0)
		title='rmasso= '//wd2
		call wtstr(.35,.73,title,1,0,0)
		!
		write(wd1,'(f6.3)')rho_wind1
		write(wd2,'(f6.3)')rho_wind2
		write(wd3,'(f6.3)')vx_wind1
		write(wd4,'(f6.3)')vx_wind2
		!
		title='rho_wind1= '//wd1
		call wtstr(.15,.7,title,1,0,0)
		title='rho_wind2= '//wd2
		call wtstr(.35,.7,title,1,0,0)
		title='vx_wind1 = '//wd3
		call wtstr(.55,.7,title,1,0,0)
		title='vx_wind2 = '//wd4
		call wtstr(.75,.7,title,1,0,0)
		!
		write(wd1,'(f6.3)')vy_wind1
		write(wd2,'(f6.3)')vy_wind2
		write(wd3,'(f6.3)')vz_wind1
		write(wd4,'(f6.3)')vz_wind2
		!
		title='vy_wind1= '//wd1
		call wtstr(.15,.67,title,1,0,0)
		title='vy_wind2= '//wd2
		call wtstr(.35,.67,title,1,0,0)
		title='vz_wind1 = '//wd3
		call wtstr(.55,.67,title,1,0,0)
		title='vz_wind2 = '//wd4
		call wtstr(.75,.67,title,1,0,0)
		!
		write(wd1,'(f6.3)')alfx_wind1
		write(wd2,'(f6.3)')alfx_wind2
		write(wd3,'(f6.3)')alfy_wind1
		write(wd4,'(f6.3)')alfy_wind2
		!
		title='alfx1 = '//wd1
		call wtstr(.15,.64,title,1,0,0)
		title='alfx2 = '//wd2
		call wtstr(.35,.64,title,1,0,0)
		title='alfy1 = '//wd3
		call wtstr(.55,.64,title,1,0,0)
		title='alfy2 = '//wd4
		call wtstr(.75,.64,title,1,0,0)
		!
		write(wd3,'(f6.3)')alfz_wind1
		write(wd4,'(f6.3)')alfz_wind2
		title='alfz1 = '//wd3
		call wtstr(.55,.61,title,1,0,0)
		title='alfz2 = '//wd4
		call wtstr(.75,.61,title,1,0,0)
		!
		write(wd1,'(f5.1)')re_wind
		write(wd2,'(f5.0)')reynolds
		write(wd3,'(f5.0)')resist
		!	re_wind sets radius from earth where initial wind placed
		!	reynolds coefficient for surface currents
		!	resist equivalent if you wish to run anomalous resistivity
		!	bfrac determines the percentage of the tangential magnetic
		!		field allowed at earth's surface
		!
		title='re_wind = '//wd1
		call wtstr(.15,.58,title,1,0,0)
		title='reynolds = '//wd2
		call wtstr(.35,.58,title,1,0,0)
		title='resist = '//wd3
		call wtstr(.55,.58,title,1,0,0)
		!
		write(wd1,'(f6.3)')bfrac
		write(wd2,'(f6.3)')vfrac
		title='bfrac = '//wd1
		call wtstr(.35,.55,title,1,0,0)
		title='vfrac = '//wd2
		call wtstr(.55,.55,title,1,0,0)
		!
		write(wd1,'(f6.3)')chirho
		write(wd2,'(f6.3)')chipxyz
		write(wd3,'(f6.3)')chierg
		!
		title='chirho = '//wd1
		call wtstr(.15,.52,title,1,0,0)
		title='chipxyz = '//wd2
		call wtstr(.35,.52,title,1,0,0)
		title='chierg = '//wd3
		call wtstr(.55,.52,title,1,0,0)
		!
		write(wd1,'(f6.3)')difrho
		write(wd2,'(f6.3)')difpxyz
		write(wd3,'(f6.3)')diferg
		!
		title='chirho = '//wd1
		call wtstr(.15,.49,title,1,0,0)
		title='chipxyz = '//wd2
		call wtstr(.35,.49,title,1,0,0)
		title='chierg = '//wd3
		call wtstr(.55,.49,title,1,0,0)
		!
		write(wd1,'(f4.1)')tilt1
		title='tilt1 = '//wd1
		call wtstr(.15,.40,title,1,0,0)
		write(wd1,'(f4.1)')tilt2
		title='tilt2 = '//wd1
		call wtstr(.30,.40,title,1,0,0)
		!
	!
	!
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!		@		RECONFIGURE INPUT DATA		@
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
	!
    !	The magnetic field, in dimensionless units, is actually in Alfven speeds
    !		relative to the normalizing velocity
    !	The temperature is in sound speeds, relative to the normalizing speed
    !		all squared
    !
    !	The magnetospheric plasma densities are assumed to vary as
    !		rho proportional to (r)**-alpha_e
    !		temperature proportional to (r)**alpha_e
    !		with total pressure constant which is needed for equilibrium
    !
    !	Now, find the equivalent pressure of magnetosphere for the given
    !		sound speed
    !
    gamma1 = gamma - 1.
    epress = cs_inner**2 * erho / gamma
    eerg = epress / gamma1
    !
    !	Do the same for the solar wind
    !
    rho_wind1 = den_wind1 * rmassq
    rho_wind2 = den_wind2 * rmassq
    srho = rho_wind1
    delrho = (rho_wind2 - rho_wind1) / tmax
    spress = (cs_wind**2 * srho / gamma) / gamma1
    svelx = vx_wind1
    svely = vy_wind1
    svelz = vz_wind1
    !
    spx = srho * svelx
    spy = srho * svely
    spz = srho * svelz
    serg = 0.5 * (svelx**2 + svely**2 + svelz**2) * srho + spress
    delvx_wind = (vx_wind2 - vx_wind1) / tmax
    delvy_wind = (vy_wind2 - vy_wind1) / tmax
    delvz_wind = (vz_wind2 - vz_wind1) / tmax
    !
    delbx_wind = ( alfx_wind2 * sqrt(rho_wind2) - alfx_wind1 * sqrt(rho_wind1) ) / tmax
    delby_wind = ( alfy_wind2 * sqrt(rho_wind2) - alfy_wind1 * sqrt(rho_wind1) ) / tmax
    delbz_wind = ( alfz_wind2 * sqrt(rho_wind2) - alfz_wind1 * sqrt(rho_wind1) ) / tmax
    sbx_wind = alfx_wind1 * sqrt(rho_wind1)
    sby_wind = alfy_wind1 * sqrt(rho_wind1)
    sbz_wind = alfz_wind1 * sqrt(rho_wind1)
    !
    deltg = tmax / float(ntgraph)
    deltinj = tmax / float(ntinj)
    delt = stepsz
	write_dat = .true.
    !
    !	*****************
    !	Check for restart
    !	*****************
    nchf=11
    if(.not.start) then
        !
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        !	@			RESTART FROM FLUID FILES			@
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        !
        if(reload) then
            nchf=15
            open(nchf,file='fluid00',status='unknown',form='unformatted')
        else
            open(11,file='fluid11',status='unknown',form='unformatted')
            open(12,file='fluid12',status='unknown',form='unformatted')
            read(nchf) t
            rewind nchf
            inchf=23-nchf
            read(inchf) t2
            rewind inchf
            !
            if(t.lt.t2) then
                close(nchf)
                nchf=inchf
            else
                close(inchf)
            endif
        endif
        !
		!	*****************
        !	Read restart data
		!	*****************
		    read(nchf)t
		    read(nchf)qrho
		    read(nchf)qpx
		    read(nchf)qpy
		    read(nchf)qpz
		    read(nchf)qpresx
		    read(nchf)qpresy
		    read(nchf)qpresz
		    read(nchf)qpresxy
		    read(nchf)qpresxz
		    read(nchf)qpresyz
		    read(nchf)hrho
		    read(nchf)hpx
		    read(nchf)hpy
		    read(nchf)hpz
		    read(nchf)hpresx
		    read(nchf)hpresy
		    read(nchf)hpresz
		    read(nchf)hpresxy
		    read(nchf)hpresxz
		    read(nchf)hpresyz
		    read(nchf)orho
		    read(nchf)opx
		    read(nchf)opy
		    read(nchf)opz
		    read(nchf)opresx
		    read(nchf)opresy
		    read(nchf)opresz
		    read(nchf)opresxy
		    read(nchf)opresxz
		    read(nchf)opresyz
		    read(nchf)bx
		    read(nchf)by
		    read(nchf)bz
		    read(nchf)epres
		    read(nchf)bx0
		    read(nchf)by0
		    read(nchf)bz0
		    read(nchf)parm_srf,parm_mid,parm_zero,ijzero,numzero,ijmid,nummid,ijsrf,numsrf
		close(nchf)
        !
		!###########################################
		ut=utstart	!	utstart read from input file
		!###########################################
		!
		!	***************
		!	Rotation timing
		!	***************
			nrot_p = int( ut/planet_per )
			rot_hrs = ut - planet_per*nrot_p
			rot_angle = 2.*pi*rot_hrs/planet_per
			nrot_m = int( ut/moon_per )
			rot_hrs_m = ut - moon_per*nrot_m
			rot_angle_m = 2.*pi*rot_hrs_m/moon_per + moon_init_rot
			!
			xmoon = moon_orbit_rad * cos(rot_angle_m)
			ymoon = moon_orbit_rad * sin(rot_angle_m)
			zmoon = 0.
			!
		!
        !	Check for div b errors
        !
        if(divb_lores) then
            range=1.33*torus_inj_dist/re_equiv
            write(*,*)'Range for divb taper: ',torus_inj_dist,range
            do box=n_grids,1,-1
                write(*,*)'divb on box: ',box
                call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                    curx,cury,curz,efldx,efldy,efldz, &
                    7,nx*ny*nz,nx,ny,nz,n_grids,box,xspac)
                call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                    curx,cury,curz,efldx,efldy,efldz, &
                    7,nx*ny*nz,nx,ny,nz,n_grids,box,xspac)
                write(*,*)'Completed divb on box: ',box
            enddo !end box loop
            !
            !	Apply boundary conditions
            !
            do box=n_grids-1,1,-1
                call flanks_synced(bx,nx,ny,nz,n_grids,box, &
                    grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
                call flanks_synced(by,nx,ny,nz,n_grids,box, &
                    grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
                call flanks_synced(bz,nx,ny,nz,n_grids,box, &
                    grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
            enddo
            do box=1,n_grids-1
                call bndry_corer( &
                    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                    orho,opresx,opresy,opresz,opx,opy,opz, &
                    epres,bx,by,bz, &
                    qpresxy,qpresxz,qpresyz, &
                    hpresxy,hpresxz,hpresyz, &
                    opresxy,opresxz,opresyz, &
                    nx,ny,nz,n_grids,box, &
                    grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
            enddo  !end bndry_corer
        endif  ! end divb_lores
		!
        !	if(update)t=0.
        write(*,*)'Entering lores visual'
        ut=utstart
        call visual(qrho,qpresx,qpresy,qpresz,qpresxy, &
            qpresxz,qpresyz,qpx,qpy,qpz,rmassq, &
            hrho,hpresx,hpresy,hpresz,hpresxy, &
            hpresxz,hpresyz,hpx,hpy,hpz,rmassh, &
            orho,opresx,opresy,opresz,opresxy, &
            opresxz,opresyz,opx,opy,opz,rmasso, &
            epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
            curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
            tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
            nx,ny,nz,n_grids,xspac, &
            cross,along,flat,xcraft,ncraft,re_equiv, &
            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
            grid_minvals(3,:), grid_maxvals(3,:), ut,b_equiv,ti_te,rho_equiv)
		!
        ts1 = t + tsave
        tstep = tmax
        tmax = t + tmax
        tgraph = t + deltg
        tinj = t + deltinj
        tdiv = t
        nchf = 11
        !
        !	Initialize plasma resistivity
        !
        call set_resist(resistive,nx,ny,nz,mbndry,resist, &
            ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
            msrf,mmid,mzero,1.)
        !
        write(*,*)'Entering lores visual'
        call visual(qrho,qpresx,qpresy,qpresz,qpresxy, &
            qpresxz,qpresyz,qpx,qpy,qpz,rmassq, &
            hrho,hpresx,hpresy,hpresz,hpresxy, &
            hpresxz,hpresyz,hpx,hpy,hpz,rmassh, &
            orho,opresx,opresy,opresz,opresxy, &
            opresxz,opresyz,opx,opy,opz,rmasso, &
            epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
            curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
            tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
            nx,ny,nz,n_grids,xspac, &
            cross,along,flat,xcraft,ncraft,re_equiv, &
            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
            grid_minvals(3,:), grid_maxvals(3,:), ut,b_equiv,ti_te,rho_equiv)
        !
        write(*,79) nchf
        79 format('  Restart from fluid_',i2)
        rewind nchf
        nchf=11
    else
        !
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!	@			INITIALIZATION			@
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!
        !###########################################
		ut=utstart	!	utstart read from input file
		!###########################################
		!
		ts1 = tsave
		!
		!	***************
		!	Rotation timing
		!	***************
			nrot_p = int( ut/planet_per )
			rot_hrs = ut - planet_per*nrot_p
			rot_angle = 2.*pi*rot_hrs/planet_per
			sin_rot = sin(rot_angle)
			cos_rot = cos(rot_angle)
			!
			nrot_m = int( ut/moon_per )
			rot_hrs_m = ut - moon_per*nrot_m
			rot_angle_m = moon_init_rot
			!
			xmoon = moon_orbit_rad * cos(rot_angle_m)
			ymoon = moon_orbit_rad * sin(rot_angle_m)
			zmoon = 0.
			!
		!
		!	**************
        !	Ambient plasma
		!	**************
			!	Initialize indices and surface points
		    !
		    numsrf=0
		    nummid=0
		    numzero=0
		    !
		    write(*,*)'Initial torus_inj_dist, rot_angle: ',torus_inj_dist, rot_angle
		    !
		!
        do box=1,n_grids
			dx = ( grid_maxvals(1,box) - grid_minvals(1,box) ) / ( nx - 1. )
			dy = ( grid_maxvals(2,box) - grid_minvals(2,box) ) / ( ny - 1. )
			dz = ( grid_maxvals(3,box) - grid_minvals(3,box) ) / ( nz - 1. )
            !
            !	Create dipole magnetic field and load magnetospheric plasma
            !
            do k=1,nz
                az = grid_minvals(3,box) + dz * (k-1) - zdip
                !
                do j=1,ny
                    ay = grid_minvals(2,box) + dy * (j-1) - ydip
                    !
                    do i=1,nx
                        ax = grid_minvals(1,box) + dx * (i-1) - xdip
                        !
                        !	Rotate coordinates for planet motion
                        !
                        xr=ax*cos_rot+ay*sin_rot
                        yr=-ax*sin_rot+ay*cos_rot
                        !
                        !	Tilt space to dipole space
                        !
                        xp=xr*cos_tilt-az*sin_tilt
                        yp=yr
                        zp=xr*sin_tilt+az*cos_tilt
                        ar=sqrt(xp**2+yp**2+zp**2)
                        !
                        !	Determine magnetic dipole field
                        !
                        call dipole(bx0(i,j,k,box),by0(i,j,k,box), &
                            bz0(i,j,k,box),ax,ay,az)
                        !
                        !	Zero interior magnetic field so alfven speed small
                        !
                        if (ar.lt.zero_bndry) then
                            bx0(i,j,k,box)=0.
                            by0(i,j,k,box)=0.
                            bz0(i,j,k,box)=0.
                        endif
                        !
                        !	Set up rotational properties
                        !
                        rx=ax*re_equiv
                        ry=ay*re_equiv
                        rd=sqrt(rx**2+ry**2)
                        if(rd.lt.r_rot) then
                            vfrac=1.
                        else
                            vfrac=((2.*r_rot**2)/(r_rot**2+ rd**2))**2
                        endif
                        rvy=rx*v_rot*vfrac
                        rvx=-ry*v_rot*vfrac
                        rvz=0.
                        corotate=sqrt(rvx**2+rvy**2)
                        !
                        !	For Jupiter top hat distribution
    					!
                        ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
                        !	isotropic
                        !	ar_iono=sqrt(xp**2+ay**2+zp**2)
                        !
                        ar_iono=amax1(0.0001,ar_iono)
                        ra=((ar_iono+0.5*r_inner)/(1.5*r_inner))**(-alpha_e)
                        zheight=amax1(1.,(zp**2+(1.5*r_inner)**2)/(3.0*r_inner)**2)
                        ra=ra/zheight**2
                        rho_iono=amax1(erho*ra,d_min)
                        !
                        r_equat=(ar**3+0.001)/(xp**2+ay**2+0.001)
                        r_equat=amax1(r_equat,r_inner)
                        rerg_sphere=(0.001+(r_inner/r_equat)**4)
                        !
                        if(r_equat.gt.5.*r_inner) then
                            rerg_sphere=(0.001+(r_inner/r_equat)**alpha_e)!constant pressure flux
                            rden_sphere=1.
                        else
                            rerg_sphere=0.1+0.9*(r_equat-r_inner)/(4.*r_inner)
                            rden_sphere=(rerg_sphere**2)          !reduce o+ in plasmasphere
                        endif
                        !
                        arho = amax1( rho_iono, d_min )
                        vt = amax1( sqrt(rvy**2 + rvx**2), cs_inner )
                        !
                        !   MAT - Corrected according to
                        !       doi:10.1029/JA094iA12p17287
                        !       and 
                        !       doi:10.1029/2008GL035433
                        !
                        qrho(i,j,k,box)=arho*rmassq/zheight
                        qpresx(i,j,k,box)=arho*vt
                        qpresy(i,j,k,box)=qpresx(i,j,k,box)
                        qpresz(i,j,k,box)=qpresx(i,j,k,box)
                        qpresxy(i,j,k,box)=0.
                        qpresxz(i,j,k,box)=0.
                        qpresyz(i,j,k,box)=0.
                        qpx(i,j,k,box)=qrho(i,j,k,box)*rvx
                        qpy(i,j,k,box)=qrho(i,j,k,box)*rvy
                        qpz(i,j,k,box)=0.
                        !
                        !	Add small amount of heavies everywhere
                        !
                        hrho(i,j,k,box)=0.05*qrho(i,j,k,box)*rmassh/rmassq
                        hpresx(i,j,k,box)=0.05*qpresx(i,j,k,box)
                        hpresy(i,j,k,box)=hpresx(i,j,k,box)
                        hpresz(i,j,k,box)=hpresx(i,j,k,box)
                        hpresxy(i,j,k,box)=0.
                        hpresxz(i,j,k,box)=0.
                        hpresyz(i,j,k,box)=0.
                        hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
                        hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
                        hpz(i,j,k,box)=0.
                        !
                        orho(i,j,k,box)=qrho(i,j,k,box)*rmasso/rmassq
                        opresx(i,j,k,box)=qpresx(i,j,k,box)
                        opresy(i,j,k,box)=opresx(i,j,k,box)
                        opresz(i,j,k,box)=opresx(i,j,k,box)
                        opresxy(i,j,k,box)=0.
                        opresxz(i,j,k,box)=0.
                        opresyz(i,j,k,box)=0.
                        opx(i,j,k,box)=orho(i,j,k,box)*rvx
                        opy(i,j,k,box)=orho(i,j,k,box)*rvy
                        opz(i,j,k,box)=0.
                        !
                        epres(i,j,k,box)=(qpresx(i,j,k,box)+hpresx(i,j,k,box)+ &
                            opresx(i,j,k,box))/ti_te
                        !
                        bx(i,j,k,box)=0.
                        by(i,j,k,box)=0.
                        bz(i,j,k,box)=0.
                        !
                        !	Test for boundary point of planets surface or
                        !		interior to planet
                        !
                        if((ar.le.srf_bndry).and.(box.le.mbndry)) then
                            !
                            if(ar.lt.zero_bndry) then
                                numzero(box)=numzero(box)+1
                                ijzero(box,1,numzero(box))=i
                                ijzero(box,2,numzero(box))=j
                                ijzero(box,3,numzero(box))=k
                                !
                                parm_zero(box,1,numzero(box))=qrho(i,j,k,box)
                                parm_zero(box,2,numzero(box))=hrho(i,j,k,box)
                                parm_zero(box,3,numzero(box))=orho(i,j,k,box)
                                parm_zero(box,4,numzero(box))=qpresx(i,j,k,box)
                                parm_zero(box,5,numzero(box))=hpresx(i,j,k,box)
                                parm_zero(box,6,numzero(box))=opresx(i,j,k,box)
                                parm_zero(box,7,numzero(box))=epres(i,j,k,box)
                                !
                            else if(ar.lt.mid_bndry) then
                                nummid(box)=nummid(box)+1
                                ijmid(box,1,nummid(box))=i
                                ijmid(box,2,nummid(box))=j
                                ijmid(box,3,nummid(box))=k
                                ar2=sqrt(xp**2+ay**2)
                                parm_mid(box,1,nummid(box))=qrho(i,j,k,box)
                                parm_mid(box,2,nummid(box))=hrho(i,j,k,box)
                                parm_mid(box,3,nummid(box))=orho(i,j,k,box)
                                parm_mid(box,4,nummid(box))=qpresx(i,j,k,box)
                                parm_mid(box,5,nummid(box))=hpresx(i,j,k,box)
                                parm_mid(box,6,nummid(box))=opresx(i,j,k,box)
                                parm_mid(box,7,nummid(box))=epres(i,j,k,box)
                                !
                            else
                                numsrf(box)=numsrf(box)+1
                                ijsrf(box,1,numsrf(box))=i
                                ijsrf(box,2,numsrf(box))=j
                                ijsrf(box,3,numsrf(box))=k
                                ar2=sqrt(xp**2+ay**2)
                                parm_srf(box,1,numsrf(box))=qrho(i,j,k,box)
                                parm_srf(box,2,numsrf(box))=hrho(i,j,k,box)
                                parm_srf(box,3,numsrf(box))=orho(i,j,k,box)
                                parm_srf(box,4,numsrf(box))=qpresx(i,j,k,box)
                                parm_srf(box,5,numsrf(box))=hpresx(i,j,k,box)
                                parm_srf(box,6,numsrf(box))=opresx(i,j,k,box)
                                parm_srf(box,7,numsrf(box))=epres(i,j,k,box)
                            endif
                        endif  ! end bndry condition test
                        !
                    enddo	! end x loop
                enddo	! end y loop
            enddo	! end z loop
        enddo	! end boxes loop
        !
        write(*,*) 'Interior and zero points:'
		write(*,*) 'box ','srf ','mid ','zero'
        do box=1,mbndry
            write(*,*) box,numsrf(box),nummid(box),numzero(box)
        enddo
        !
        !	Initialize solar wind plasma. Inserted beyond
        !		the planet at a radius of re_wind.
        !
        wind_bnd = r_rot / re_equiv
        ofrac = rho_frac
        do box=1,n_grids
            !
            dx = ( grid_maxvals(1,box) - grid_minvals(1,box) ) / (nx-1.)
            dy = ( grid_maxvals(2,box) - grid_minvals(2,box) ) / (ny-1.)
            dz = ( grid_maxvals(3,box) - grid_minvals(3,box) ) / (nz-1.)
            !
            do k=1,nz
                az = grid_minvals(3,box) + dz * (k-1) - zdip
                do j=1,ny
                    ay = grid_minvals(2,box) + dy * (j-1) - ydip
                    do i=1,nx
                        ax = grid_minvals(1,box) + dx * (i-1) - xdip
    					!
                        ar=sqrt(ax**2+ay**2+az**2)
                        if( (ar.ge.1.5*wind_bnd) .and. (ax.lt.0.) ) then
                            qrho(i,j,k,box)=srho+qrho(i,j,k,box)
                            qpx(i,j,k,box)=spx+qpx(i,j,k,box)
                            qpy(i,j,k,box)=spy+qpy(i,j,k,box)
                            qpz(i,j,k,box)=spz+qpz(i,j,k,box)
                            !
                            avx=qpx(i,j,k,box)/qrho(i,j,k,box)
                            avy=qpy(i,j,k,box)/qrho(i,j,k,box)
                            avz=qpz(i,j,k,box)/qrho(i,j,k,box)
                            !
                            apres=cs_wind**2*qrho(i,j,k,box)/gamma
                            qpresx(i,j,k,box)=0.5*apres+qpresx(i,j,k,box)
                            qpresy(i,j,k,box)=qpresx(i,j,k,box)
                            qpresz(i,j,k,box)=qpresx(i,j,k,box)
                            qpresxy(i,j,k,box)=0.
                            qpresxz(i,j,k,box)=0.
                            qpresyz(i,j,k,box)=0.
                            !
                            epres(i,j,k,box)=0.5*apres/ti_te+epres(i,j,k,box)
                            !
                            hrho(i,j,k,box)=srho*rho_frac*rmassh/rmassq+ &
                                hrho(i,j,k,box)
                            hpx(i,j,k,box)=avx*hrho(i,j,k,box)
                            hpy(i,j,k,box)=avy*hrho(i,j,k,box)
                            hpz(i,j,k,box)=avz*hrho(i,j,k,box)
                            hpresx(i,j,k,box)=0.5*spress*rho_frac &
                                +hpresx(i,j,k,box)
                            hpresy(i,j,k,box)=hpresx(i,j,k,box)
                            hpresz(i,j,k,box)=hpresx(i,j,k,box)
                            hpresxy(i,j,k,box)=0.
                            hpresxz(i,j,k,box)=0.
                            hpresyz(i,j,k,box)=0.
                            !
                            orho(i,j,k,box)=srho*ofrac*rmasso/rmassq+ &
                                orho(i,j,k,box)
                            opx(i,j,k,box)=avx*orho(i,j,k,box)
                            opy(i,j,k,box)=avy*orho(i,j,k,box)
                            opz(i,j,k,box)=avz*orho(i,j,k,box)
                            opresx(i,j,k,box)=0.5*spress*ofrac &
                                +opresx(i,j,k,box)
                            opresy(i,j,k,box)=opresx(i,j,k,box)
                            opresz(i,j,k,box)=opresx(i,j,k,box)
                            opresxy(i,j,k,box)=0.
                            opresxz(i,j,k,box)=0.
                            opresyz(i,j,k,box)=0.
                        endif
						!
                        if((ar.gt.wind_bnd).and.(ar.lt.1.5*wind_bnd) &
                            .and.(ax.lt.0.0)) then
                            dfrac=(ar-wind_bnd)/(0.5*wind_bnd)
                            qrho(i,j,k,box)=srho*dfrac+qrho(i,j,k,box)
                            qpx(i,j,k,box)=spx*dfrac+qpx(i,j,k,box)
                            qpy(i,j,k,box)=spy*dfrac+qpy(i,j,k,box)
                            qpz(i,j,k,box)=spz*dfrac+qpz(i,j,k,box)
                            avx=qpx(i,j,k,box)/qrho(i,j,k,box)
                            avy=qpy(i,j,k,box)/qrho(i,j,k,box)
                            avz=qpz(i,j,k,box)/qrho(i,j,k,box)
                            !
                            apres=cs_wind**2*qrho(i,j,k,box)/gamma
                            qpresx(i,j,k,box)=0.5*apres*dfrac+qpresx(i,j,k,box)
                            qpresy(i,j,k,box)=qpresx(i,j,k,box)
                            qpresz(i,j,k,box)=qpresx(i,j,k,box)
                            qpresxy(i,j,k,box)=0.
                            qpresxz(i,j,k,box)=0.
                            qpresyz(i,j,k,box)=0.
							!
                            epres(i,j,k,box)=0.5*apres/ti_te*dfrac+epres(i,j,k,box)
							!
                            hrho(i,j,k,box)=srho*rho_frac*rmassh/rmassq*dfrac+ &
                                hrho(i,j,k,box)
                            hpx(i,j,k,box)=avx*hrho(i,j,k,box)
                            hpy(i,j,k,box)=avy*hrho(i,j,k,box)
                            hpz(i,j,k,box)=avz*hrho(i,j,k,box)
                            hpresx(i,j,k,box)=0.5*spress*rho_frac*dfrac &
                                +hpresx(i,j,k,box)
                            hpresy(i,j,k,box)=hpresx(i,j,k,box)
                            hpresz(i,j,k,box)=hpresx(i,j,k,box)
                            hpresxy(i,j,k,box)=0.
                            hpresxz(i,j,k,box)=0.
                            hpresyz(i,j,k,box)=0.
							!
                            orho(i,j,k,box)=srho*ofrac*rmasso/rmassq*dfrac+ &
                                orho(i,j,k,box)
                            opx(i,j,k,box)=avx*orho(i,j,k,box)
                            opy(i,j,k,box)=avy*orho(i,j,k,box)
                            opz(i,j,k,box)=avz*orho(i,j,k,box)
                            opresx(i,j,k,box)=0.5*spress*ofrac*dfrac &
                                +opresx(i,j,k,box)
                            opresy(i,j,k,box)=opresx(i,j,k,box)
                            opresz(i,j,k,box)=opresx(i,j,k,box)
                            opresxy(i,j,k,box)=0.
                            opresxz(i,j,k,box)=0.
                            opresyz(i,j,k,box)=0.
                        endif
                        !
                    enddo
                enddo
            enddo
        enddo
        !
        !	Check speeds
        !
        do  box=n_grids,1,-1
            do  k=1,nz
                do  j=1,ny
                    do  i=1,nx
                        avx=qpx(i,j,k,box)/qrho(i,j,k,box)
                        avy=qpy(i,j,k,box)/qrho(i,j,k,box)
                        avz=qpz(i,j,k,box)/qrho(i,j,k,box)
                        spd=sqrt(avx**2+avy**2+avz**2)
                        if(spd.gt.1.) then
                            qpx(i,j,k,box)=0.
                            qpy(i,j,k,box)=0.
                            qpz(i,j,k,box)=0.
                        endif
                        !
                        avx=hpx(i,j,k,box)/hrho(i,j,k,box)
                        avy=hpy(i,j,k,box)/hrho(i,j,k,box)
                        avz=hpz(i,j,k,box)/hrho(i,j,k,box)
                        spd=sqrt(avx**2+avy**2+avz**2)
                        if(spd.gt.1.) then
                            hpx(i,j,k,box)=0.
                            hpy(i,j,k,box)=0.
                            hpz(i,j,k,box)=0.
                        endif
                        !
                        avx=opx(i,j,k,box)/orho(i,j,k,box)
                        avy=opy(i,j,k,box)/orho(i,j,k,box)
                        avz=opz(i,j,k,box)/orho(i,j,k,box)
                        spd=sqrt(avx**2+avy**2+avz**2)
                        if(spd.gt.1.) then
                            opx(i,j,k,box)=0.
                            opy(i,j,k,box)=0.
                            opz(i,j,k,box)=0.
                        endif
                        !
                    enddo
                enddo
            enddo
        enddo
		!
		!	Initialize plasma resistivity
		!
		call set_resist(resistive,nx,ny,nz,mbndry,resist, &
		    ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
		    msrf,mmid,mzero,1.)
		!
    endif ! if(.not.start) then
	!
	!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!						FINISHED INITIAL SETUP!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
    !
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) 'Did we start from scratch?', start
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) 'UT = ', ut
		!
	!
	!
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !		@		SPACECRAFT DATA INITIALIZATION		@
    !		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!
    !
	if(spacecraft) then
		!
		!	Count crafts named in craft_info directory.
		call countcraft(craft_info,scin,ncraft,ndef_craft,naux_craft)
		!
		allocate (xcraft(4,ncraft),zcraft(4,ncraft))
		allocate (craftnames(ncraft),ntimes(ncraft,2),recording(ncraft))
		allocate (craft_gridpt(4,ncraft))
		!
		!	Reference craft is a copy of 'wind' default craft
		rcraft(1)=xcraft(1,1)
		rcraft(2)=xcraft(2,1)
		rcraft(3)=xcraft(3,1)
		!
		!	Read auxiliary spacecraft names into craftnames array
		!	We do aux first so we can allocate craftpos array before reading default craft data.
		!	Step 1: ls .craft names into crafts.txt file. 
		!		This file will contain only .craft names, not default craft.
		!		We use wc -l *.craft so that we get numbers of lines as well as file names.
		call system ('wc -l '//trim(craft_info)//'*.craft > '//trim(craft_info)//'crafts.txt')
		!
		open(scin,file=trim(craft_info)//'crafts.txt',status='unknown',form='formatted')
			do n=ndef_craft+1,ncraft
				!	Step 2: Read .craft filenames into placeholder string
				read(scin,*) ntimes(n,2), junkline
				ntimes(n,2) = ntimes(n,2) - nheadlines	!	Excludes header line in craft file
				!	Step 3: Print everything before the .craft extension into craftnames array
				craftnames(n) = junkline(1:index(junkline,'.')-1)
			enddo
		close(scin)
		write(*,*) 'Spacecraft found: ', craftnames
		!
		allocate(craftpos(4,ncraft,maxval(ntimes(:,2))))	! Allocate array which stores all spacecraft position/time input values
		!
		!	Now we can read in default craft info
		!
		open(scin+1,file=trim(craft_info)//'defaults.pos',status='unknown',form='formatted')
			read(scin+1,*) junkline	!	Skip header line, only one line in default craft file
			!
			do n=1,ndef_craft
				recording(n)=.true.
				read(scin+1,*) craftnames(n), craftpos(1,n,1), &
					craftpos(2,n,1), craftpos(3,n,1)
				craftpos(4,n,1) = utstart	!	FIX NEEDED: utstart is not correct ut reference
				!
				!	Initial positions of default spacecraft in re but simulation directions
				craftpos(1:3,n,1) = craftpos(1:3,n,1) * re_equiv
				ntimes(n,2) = -1	!	We don't have a set number of recordings for default craft
				!
				!	Default craft must be placed on grid points.
				!	If not, snap to nearest grid point.
				call findgrid(craftpos(1:3,n,1),n_grids,grid_minvals,grid_maxvals,craft_gridpt(:,n))
				cbox = craft_gridpt(4,n)					
				do axis=1,3
					craftpos(axis,n,1) = grid_minvals(axis,cbox) + xspac(cbox)*re_equiv*craft_gridpt(axis,n)
				enddo
				xcraft(:,n) = craftpos(:,n,1)	!	Set spacecraft actual position to first value in craftpos array
			enddo
		close(scin+1)
		!
		do n=ndef_craft+1, ncraft
			recording=.true.
			craftstat=0
			open(scin+n,file=trim(craft_info)//trim(craftnames(n))//'.craft', &
				status='unknown',iostat=craftstat,form='formatted')
				do headnum=1, nheadlines
					read(scin+n,*) junkline		!	Skip header lines
				enddo
				!
				do m=1,ntimes(n,2)
					!	.craft files are formatted as t, x, y, z
					read(scin+n,*) craftpos(4,n,m), craftpos(1,n,m), &
						craftpos(2,n,m), craftpos(3,n,m)
				enddo
				!
				if (craftstat .ne. 0) then
					write(*,*) 'Error reading craft file: ', craftnames(n)
					ntimes(n,2) = 0
				elseif (ntimes(n,2).eq.1) then
					write(*,*) 'Craft file has only one xyzt coordinate:', craftnames(n)
					ntimes(n,2) = 0
				endif
				!
			close(scin+n)
		enddo
		!
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!	@		FINISHED READING IN SPACECRAFT INPUTS		@
    	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!
	    !xcraft(4,n) = zut	!	I don't know why this is being done. Removed because it will interfere with read-in locations otherwise. MJS 08/11/17
		!	This is what was being read in from spacecraft input files before my upgrade: MJS 08/12/17
		!	do m=1,ncts
        !	    read(29,*)bfld(m,4),bmag,bfld(m,1),bfld(m,2),bfld(m,3)
        !	    read(27,*)rut,rplas(m),svel(m,1),svel(m,2),svel(m,3)
		!	enddo
		!	It doesn't seem like any of it was being used.
		!
		call limcraft(xcraft,ncraft,re_equiv,n_grids, &
			grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
			grid_minvals(3,:), grid_maxvals(3,:) )
		!
        !zcraft(:,:)=xcraft(:,:)
		!
	else	!	If we are not using spacecraft, we can't use spacecraft-only functions.
		craft_input = .false.
		warp = .false.	
    endif	! endif(spacecraft)
	!
	!
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!		@		INITIALIZE SPACECRAFT OUTPUT FILES		@
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!
	!
    if(spacecraft) then
        !
        !	Open and go to end of .dat files for each craft
        !
        do n=1,ncraft
			fname=trim(craft_data)//trim(craftnames(n))//'.dat'
			inquire(file=fname, exist=dat_exists)
			craftstat=0
			!	Open .dat files for each craft. Files will be closed automatically when process completes or dies.
	    	if(.not.dat_exists) then
				open(scout+n,file=fname,iostat=craftstat,status='unknown',form='formatted')
				write(scout+n,crafthead)
				write(scout+n,*) dat_header
				zcraft(:,n) = craftpos(:,n,2)	!	If we haven't created a .dat file for this spacecraft yet, 
				xcraft(:,n) = craftpos(:,n,1)		!		its next recording should be at the first position in its xyz sequence
				ntimes(n,1) = 0
			else
				call system ('wc -l '//fname//' > '//trim(craft_data)//'dat_working.txt')
				open(scout,file=trim(craft_data)//'dat_working.txt',status='unknown',form='formatted')
					read(scout,*) ntimes(n,1), junkline
					ntimes(n,1) = ntimes(n,1) - nheadlines	!	Subtract header lines
				close(scout)
				!
				if(ntimes(n,2).le.ntimes(n,1) .and. n.gt.ndef_craft) then	!	If true, this is an aux. craft which has finished recording its data
					recording(n) = .false.
					write(*,*) 'Spacecraft ',craftnames(n),' was already finished recording before start.'
				else
					open(scout+n,file=fname,iostat=craftstat,status='unknown',form='formatted')
					do m=1, (ntimes(n,1) + nheadlines)	!	Step through the file until the end, so we are ready to write new measurements. We have to skip past header lines, too.
						read(scout+n,*) junkline
					enddo
					xcraft(:,n) = craftpos(:,n,ntimes(n,1))
					zcraft(:,n) = craftpos(:,n,ntimes(n,1)+1)
				endif
				if(craftstat.ne.0) then
					write(*,*) 'Problem reading .dat file for craft: ', craftnames(n)
					recording(n) = .false.
				endif
			endif
		enddo
        !
        !	Keep spacecraft within grid boundaries
        call limcraft(zcraft,ncraft,re_equiv,n_grids, &
            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
            grid_minvals(3,:), grid_maxvals(3,:) )
        !
        !	Calculate the distance between the reference spacecraft and
        !		solar wind boundary
        !
        distance=grid_minvals(1,n_grids)-rcraft(1)/re_equiv
        !write(*,*)'Wind displaced by: ',distance	! Reference craft used to be used to track leading edge of input solar wind. Not any more. MJS 08/18/17
		!
		!	Set current spacecraft location to init or last read in from .dat
        xcraft(:,:)=zcraft(:,:)
        !
        do  k=1,nz
            do  j=1,ny
                ncount(j,k)=0	!	?[]? What are these variables for?
                future(j,k)=ut-0.01
            enddo
        enddo
        !
		if(craft_input) then
			!
		    !	Read all the magnetic field data to minimize data sorting
			!	?[]? What data sorting is this even talking about?
		    !	?[]? Untraceable file units, unknown previous use, obsolete. MJS 08/10/17
		    !do m=1,ncts
		    !    read(29,*)bfld(m,4),bmag,bfld(m,1),bfld(m,2),bfld(m,3)
		    !    read(27,*)rut,rplas(m),svel(m,1),svel(m,2),svel(m,3)
		    !    !      read(27,*)rut,rplas(m)
		    !    !      read(28,*)vut,svel(m,1),svel(m,2),svel(m,3)
		    !    !      warning recalibration
		    !          !   keep bx in solar wind constant
		    !    !      bfld(m,1)=-sbx_wind*b_equiv
		    !enddo
		    !
		    !	Set timing
		    !		
		    nvx=0
		    vut=-999.
		    do while(ut.gt.vut)
		        svelx=zvelx
		        nvx=nvx+1
		        zvelx=-svel(nvx,1)/v_equiv
		        vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.	!	REQUIRES SPACECRAFT INPUT DATA FOR bfld!!!	MJS 08/12/17
		    enddo
		    !
		    write(*,*)'UT=',ut,' Wind time: ',bfld(nvx,4)
		    !
		    displace=0.
		    dx = ( grid_maxvals(1,n_grids) - grid_minvals(1,n_grids) ) / (nx-1.)
		    dy = ( grid_maxvals(2,n_grids) - grid_minvals(2,n_grids) ) / (ny-1.)
		    dz = ( grid_maxvals(3,n_grids) - grid_minvals(3,n_grids) ) / (nz-1.)
		    !
		    do k=1,nz
		        do j=1,ny
		            do while((ut.gt.future(j,k)) &
		                .and.(ncount(j,k)+1.le.ncts))
		                nc=ncount(j,k)+1
		                bxp(j,k)=bxf(j,k)
		                byp(j,k)=byf(j,k)
		                bzp(j,k)=bzf(j,k)
		                rhop(j,k)=rhof(j,k)
		                svxp(j,k)=svxf(j,k)
		                svyp(j,k)=svyf(j,k)
		                svzp(j,k)=svzf(j,k)
		                past(j,k)=future(j,k)
		                !
		                future(j,k)=bfld(nc,4)
		                bxf(j,k)=-bfld(nc,1)/b_equiv
		                byf(j,k)=-bfld(nc,2)/b_equiv
		                bzf(j,k)=bfld(nc,3)/b_equiv
		                rhof(j,k)=rplas(nc)/rho_equiv
		                svxf(j,k)=-svel(nc,1)/v_equiv
		                svyf(j,k)=0.
		                svzf(j,k)=0.
		                ncount(j,k)=nc
		                avx=svxf(j,k)
		                !
		                !	Calculate delay
		                !
		                if(warp) then
							stop
							!	'warp' is not implemented.
		                    b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
		                    b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
							!	WARNING: This is not in a loop over 'box'.
							ay = grid_minvals(2,box) + dy * (j-1) - rcraft(2) / re_equiv
		                    az = grid_minvals(3,box) + dz * (k-1) - rcraft(3) / re_equiv
		                    !
		                    !       Assuming Bz IMF is positive on average 
							!			and that we can ignore transients
		                    !       Also assuming By IMF is negative
		                    displace=-bxf(j,k)* &
		                        (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
		                endif
		                ar=distance+displace
		                future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
		            enddo
		        enddo
		    enddo
        	!
		    call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
		        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
		        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
		        orho,opresx,opresy,opresz,opx,opy,opz, &
		        rmassq,rmassh,rmasso,epres, &
		        qpresxy,qpresxz,qpresyz, &
		        hpresxy,hpresxz,hpresyz, &
		        opresxy,opresxz,opresyz, &
		        rhop,svxp,svyp,svzp,svelx,spress, &
		        ti_te,rho_frac,nx,ny,nz,n_grids)
		endif	!	end craft_input if
    endif   !	end spacecraft if
    !
    !	Check initial conditions
    !
    !      write(*,*) 'checking set speed'
    do box=n_grids,1,-1
        !
        !	Sync time steps
        !
        t_old(box)=0.
        t_new(box)=0.
        !
        !	Check density
        !
        call set_rho(qrho,qpresx,qpresy,qpresz, &
            qpresxy,qpresxz,qpresyz,rmassq, &
            hrho,hpresx,hpresy,hpresz, &
            hpresxy,hpresxz,hpresyz,rmassh, &
            orho,opresx,opresy,opresz, &
            opresxy,opresxz,opresyz,rmasso, &
            epres,nx,ny,nz,n_grids,box,o_conc)
        !
        !	Check speeds of individual grids
        !
        vlim=0.6*box
        !	Find maximum velocity to determine time step
        call set_speed_agrd( &
            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
            orho,opresx,opresy,opresz,opx,opy,opz, &
            epres,qpresxy,qpresxz,qpresyz, &
            hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
            bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
            vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
            rmassq,rmassh,rmasso,nx,ny,nz,n_grids,box, &
            pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
            vlim,alf_lim,o_conc,fastest,isotropic)
        write(*,195)box,csmax,alfmax,pxmax,pymax,pzmax
        195   format(1x,i2,5(1x,1pe12.5))
        !
        t_stepnew(box)=stepsz*xspac(box)/fastest
    enddo
    write(*,*)'Speeds checked.'
    !
    delt_old=delt
    delt=stepsz/fastest
    delt=amin1(delt,1.25*delt_old)
    delt=amax1(3.e-3,delt)
    write(*,163)t,delt,ut,bfld(nvx,4)
    163 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=', &
        1pe12.5,' wind time',1pe12.5)
    !
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    write(*,*) '    INITIALIZATION COMPLETED!    '
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
    !
	!
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !		@		********************************		@
    !		@			START MAIN TIME SEQUENCE			@
    !     	@		********************************		@
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
	!
    do while(t.lt.tmax)
        !
        write(*,*)'Start main loop, t=',t
        !
        do box=n_grids,1,-1
            t_step(box)=t_stepnew(box)
            !
            !	Sync time steps
            !
            t_old(box)=0.
            t_new(box)=0.
			!
            if(box.eq.n_grids) then
                smallest_step=t_step(box)
                m_smallest=box
                !write(*,*)'Syncing box, t_step(box): ',box,t_step(box)
            else
                !
                !	Check to see if box grid doesnt outstep grid box+1
                !write(*,*)'Unsync box, t_step(box): ',box,t_step(box)
                if(t_step(box).gt.t_step(box+1))t_step(box)=t_step(box+1)
                !
                if (smallest_step.gt.t_step(box)) then
                    smallest_step=t_step(box)
                    m_smallest=box
                endif
            endif
        enddo
        !
        !     set variable steps for each grid box
        !       round_step size off
        !
        !       i_step=smallest_step*10000
        !       smallest_step=i_step/10000.
        !
        !     write(*,*)'unsync steps',t_step,mallest_step,smallest_step
        !
        do box=1,n_grids
            if(box.le.m_smallest) then
                t_step(box)=smallest_step
                !write(*,*)'Sync step: box, t_step(box): ',box,t_step(box)
            else
                astep=(t_step(box)/t_step(box-1)+.50)  !round up as needed
                nsteps=astep                          !nearest integer
                !           write(*,*)astep,nsteps
                if(nsteps.gt.2)nsteps=2
                if(nsteps.lt.1)nsteps=1
                t_step(box)=t_step(box-1)*nsteps
                !write(*,*)'Sync step: box, t_step(box): ',box,t_step(box)
            endif
        enddo
        !
        m_step=t_step(n_grids)/t_step(m_smallest)
        !write(*,*)'lores steps: ',m_smallest,m_step,t_step
        !
        delt=t_step(n_grids)
        told=t
        utold=ut
        old_tilt=tilt
        !
        t=t+delt
        ut=t*t_equiv/3600.
        tilt=tilt+dtilt*delt
        !delay=t_equiv*distance/svelx/3600.	!	Does not appear to be used. MJS 08/18/17 
        !
        !write(*,201)t,delt,ut,bfld(nvx,4)
        !201 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=', &
        !    1pe12.5,' wind time',1pe12.5)
        !
        nrot=ut/planet_per
        rot_hrs=ut-nrot*planet_per
        rot_angle=2.*pi*rot_hrs/planet_per
		xmoon=cos(rot_angle)	!	?[]? Needs to be re-examined. MJS 09/26/17
		ymoon=sin(rot_angle)
		zmoon=0
        !
		!	Write to fluxes.dat file
        write(3,*) flux_header
        do box=n_grids,1,-1
            scale=rho_equiv*v_equiv*1.e5* &
            (xspac(box)*planet_rad*re_equiv*1.e5)**2
	        !	Calculate flux from torus
            call flux_counter(qpx,qpy,qpz,hpx,hpy,hpz, &
                opx,opy,opz,nx,ny,nz,n_grids,box, &
                rmassq,rmassh,rmasso, &
                scale,qflux_in,qflux_out,hflux_in,hflux_out, &
                oflux_in,oflux_out)
            write(3,*)ut,box,qflux_in,hflux_in,oflux_in,'grid_in'
            write(3,*)ut,box,qflux_out,hflux_out,oflux_out,'grid_out'
        enddo
        !
        !     ?[write out data if necessary - only using coarse gridding
        !        the moment]?
        !
		if(spacecraft) then
			!
		    !	Update position of moon diagnostic craft
		    if(trim(craftnames(2)) == 'moondgn') then
		    	xcraft(1,2)=xmoon*re_equiv*1.05
		    	xcraft(2,2)=ymoon*re_equiv*1.05
		    	xcraft(3,2)=zmoon*re_equiv*1.05
				xcraft(4,2)=ut
			endif
		    !
		    !do n=1,ncraft
		    !    call qvset(0.,bsx,nx*ny*nz)	! Sets spacecraft B components (bsx,bsy,bsz) to be zero initially
		    !    call qvset(0.,bsy,nx*ny*nz)
		    !    call qvset(0.,bsz,nx*ny*nz)
		    !    !
		    !    !
		    !    box=1
		    !    do while ( ((xcraft(1,n).gt.grd_xmax(box)*re_equiv).or. &
		    !        (xcraft(1,n).lt.grd_xmin(box)*re_equiv).or. &
		    !        (xcraft(2,n).gt.grd_ymax(box)*re_equiv).or. &
		    !        (xcraft(2,n).lt.grd_ymin(box)*re_equiv).or. &
		    !        (xcraft(3,n).gt.grd_zmax(box)*re_equiv).or. &
		    !        (xcraft(3,n).lt.grd_zmin(box)*re_equiv)) &
			!		.and. (box+1.le.n_grids) )
		    !        box=box+1	!	Finds the smallest box size the spacecraft location is inside of
		    !    enddo
		    !    rx=xspac(box)	!	Sets ?[rx]? to be spacing between grid points of smallest matching grid
		    !    ry=rx
		    !    rz=rz
		    !    !
		    !    call totfld(bx,bx0,bsx,nx,ny,nz,n_grids,box)
		    !    call totfld(by,by0,bsy,nx,ny,nz,n_grids,box)
		    !    call totfld(bz,bz0,bsz,nx,ny,nz,n_grids,box)
		    !    !
		    !    add_dip=.false.
		    !    !call crafdatv(bsx,bsy,bsz, &
		    !    !    qpx,qpy,qpz,qrho,qpresx,qpresy,qpresz,rmassq, &
		    !    !    hpx,hpy,hpz,hrho,hpresx,hpresy,hpresz,rmassh, &
		    !    !    opx,opy,opz,orho,opresx,opresy,opresz,rmasso, &
		    !    !    epres,nx,ny,nz,n_grids,box,xcraft,ncraft,n,ut, &
		    !    !    re_equiv,b_equiv,v_equiv,rho_equiv,gamma, &
		    !    !    grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
		    !    !    grid_minvals(3,:), grid_maxvals(3,:) )
		    !enddo
			!
! This entire block seems redundant and unnecessary. MJS 08/13/17 ?[]?
!	If data is read in from .craft files this step will break or do nothing.
			!
		    !     test to see whether spacecraft positions need to be updated
		    !
		    !    n=1	
		    !    if(ut.ge.zcraft(4,n)) then
		    !        read(scin+n,*)zcraft(4,n),zcraft(1,n), &
		    !            zcraft(2,n),zcraft(3,n)
		    !        !
		    !        !          Change direction to get from GSM to simulation coords
		    !        !
		    !        zcraft(1,n)=-zcraft(1,n)
		    !        zcraft(2,n)=-zcraft(2,n)
		    !        !
		    !        endif
		    !        !
		    !    endif
			!	rcraft(1)=zcraft(1,1)
		    !    rcraft(2)=zcraft(2,1)
		    !    rcraft(3)=zcraft(3,1)
		    !    !
		    !    !zcraft(4,3)=zcraft(4,2)	!	These 2 lines don't seem necessary. MJS 08/11/17
		    !    !zcraft(4,4)=zcraft(4,2)
		    !    !
		    !    !         set reference spacecraft position
		    !    !                   and spacecraft limits
		    !    !
		    !    call limcraft(zcraft,ncraft,re_equiv,n_grids, &
		    !        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
		    !        grid_minvals(3,:), grid_maxvals(3,:) )
		    !    !
		    !    !         set density and velocity
		    !    !
		    !    distance=grd_xmin(n_grids)-rcraft(1)/re_equiv
! End redundant block
			!
            !do while (ut.ge.vut)
            !    nvx=nvx+1
            !    zvelx=-svel(nvx,1)/v_equiv
            !    zvely=-svel(nvx,2)/v_equiv
            !    !         zvely=-svel(nvx,2)/v_equiv+0.03
            !    zvelz=svel(nvx,3)/v_equiv
            !    vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
            !end do
            !!
            !!	?[Fix up]? magnetic field
            !!
            !displace=0.
            !do k=1,nz
            !    do j=1,ny
            !        do while((ut.ge.future(j,k)) &
            !            .and.(ncount(j,k)+1.le.ncts))
            !            nc=ncount(j,k)+1
            !            future(j,k)=bfld(nc,4)
            !            bxf(j,k)=-bfld(nc,1)/b_equiv
            !            byf(j,k)=-bfld(nc,2)/b_equiv
            !            bzf(j,k)=bfld(nc,3)/b_equiv
            !            rhof(j,k)=rplas(nc)/rho_equiv
            !            svxf(j,k)=-svel(nc,1)/v_equiv
            !            svyf(j,k)=0.00
            !            svzf(j,k)=0.0
            !            !          svyf(j,k)=-svel(nc,2)/v_equiv
            !            !          svyf(j,k)=-svel(nc,2)/v_equiv+0.03
            !            !          svzf(j,k)=svel(nc,3)/v_equiv
            !            avx=svxf(j,k)
            !            ncount(j,k)=nc
            !            !
            !            !      Calculate delay
            !            !
                        if(warp) then
                            !           b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
                            !           b_perp=amax1(b_perp,0.1*abs(bxf(j,k)))
                            ay = (j-1.) * xspac(n_grids) + grid_minvals(2,n_grids) - rcraft(2) / re_equiv
                            az = (k-1.) * xspac(n_grids) + grid_minvals(3,n_grids) - rcraft(3) / re_equiv
                            !
                            !       Assuming Bz IMF is positive on average 
							!			and that we can ignore transients
		                    !       Also assuming By IMF is negative
                            b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
                            b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
                            displace=-bxf(j,k)* &
                                (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
							ar=distance+displace
	                        future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
                        endif
            !          enddo
            !      enddo
            !  enddo
            !	Update spacecraft position and magnetic field
			!dut=(ut-xcraft(4,n))/(zcraft(4,n)-utold)
            !xcraft(1,n)=xcraft(1,n)+dut*(zcraft(1,n)-xcraft(1,n))
            !xcraft(2,n)=xcraft(2,n)+dut*(zcraft(2,n)-xcraft(2,n))
            !xcraft(3,n)=xcraft(3,n)+dut*(zcraft(3,n)-xcraft(3,n))
            !xcraft(4,n)=ut
            !
			!
			!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            !		@		SPACECRAFT DATA RECORDING		@
			!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			!
			!	---------------------
			!		DEFAULT CRAFT
			!	---------------------
			!
			!	Only record values for default craft every deflt_rec_skp loops
			if(nskipped.ge.deflt_rec_skp .or. t.ge.tgraph) then
				write(*,*) 'Recording default spacecraft measurements.'
				do n=1,ndef_craft
					if(recording(n)) then
						!	WARNING: A line must be added below each time num_inst is incremented
						!		to add the new quantity to the array meas_qty
						scdata(1) = bx(craft_gridpt(1,n),craft_gridpt(2,n),craft_gridpt(3,n),craft_gridpt(4,n))
						scdata(2) = by(craft_gridpt(1,n),craft_gridpt(2,n),craft_gridpt(3,n),craft_gridpt(4,n))
						scdata(3) = bz(craft_gridpt(1,n),craft_gridpt(2,n),craft_gridpt(3,n),craft_gridpt(4,n))
						!	Default craft xyz are in simulation coordinates, including measurements
						write(scout+n) ut, xcraft(1:3,n), scdata
					endif
				enddo
				!
				ntimes(n,1) = ntimes(n,1) + 1
				nskipped = 0
			else
				nskipped = nskipped + 1
			endif
			!
			!	-----------------
			!		AUX CRAFT
			!	-----------------
			!
            do n=ndef_craft+1,ncraft
				if( recording(n) .and. (utold.le.xcraft(4,n) .and. ut.ge.xcraft(4,n)) ) then
					!	NEEDED: Edge case for when xcraft is the last value. Perhaps just make zcraft stay the same upon read when EOF is reached, when zcraft values are updated from craftpos
					call linterp(xcraft(:,n),zcraft(:,n),ut,sxyz)
					!	Convert from GSM to simulation coordinates, for grid ID:
					sxyz(1) = -sxyz(1)
					sxyz(2) = -sxyz(2)
					!	Plan: find indices of nearest grid point, then check if that point is above/below along each axis. Then go the other way to find next nearest point along that axis.
					!	Find closest xyz values below
					call findgrid(sxyz,n_grids,grid_minvals,grid_maxvals,num_pts,craft_gridpt)
					cbox = craft_gridpt(4,n)
					!	craft_gridpt contains x,y,z,box indices: the indices of the closest xyz LESS THAN the craft location, and the smallest box the craft fits within.
					!	x index is set to zero if there is a problem.
					if(craft_gridpt(1,n) .le. 0) then
						recording(n) = .false.
						close(scout+n)
						write(*,*) 'Gridding problem with craft ', craftnames(n), &
						'at UT = ', ut, ' Future points skipped.'
						cycle
					endif
					!	Indices in craft_gridpt are used to identify the physical parameters at those points
					!	Position values for the vertices of the grid point cube surrounding sxyz are stored in gridpts
					do axis=1,3
						gridpts(axis,1) = grid_minvals(axis,cbox) + (craft_gridpt(axis,n)-1.)*xspac(cbox)*re_equiv
						gridpts(axis,2) = gridpts(axis,1) + xspac(cbox)*re_equiv
					enddo
					!	WARNING: A line must be added below each time num_inst is incremented
					!		to add the new quantity to the array meas_qty. meas_qty values are in GSM
					!	xy directions are negated to get from sim coordinates to GSM
					meas_qty(:,:,:,1) = -bx( craft_gridpt(1,n):craft_gridpt(1,n)+1, craft_gridpt(2,n):craft_gridpt(2,n)+1, craft_gridpt(3,n):craft_gridpt(3,n)+1, craft_gridpt(4,n) )
					meas_qty(:,:,:,2) = -by( craft_gridpt(1,n):craft_gridpt(1,n)+1, craft_gridpt(2,n):craft_gridpt(2,n)+1, craft_gridpt(3,n):craft_gridpt(3,n)+1, craft_gridpt(4,n) )
					meas_qty(:,:,:,3) = bz( craft_gridpt(1,n):craft_gridpt(1,n)+1, craft_gridpt(2,n):craft_gridpt(2,n)+1, craft_gridpt(3,n):craft_gridpt(3,n)+1, craft_gridpt(4,n) )
					!	Separate trilin_interp calls needed for each scalar measurement and each component of vector quantities
					do nn=1,num_inst
						call trilin_interp( sxyz, gridpts, meas_qty(:,:,:,nn), scdata(nn) )
					enddo
					!	Convert back to GSM now that we have finished with calculations
					sxyz(1) = -sxyz(1)
					sxyz(2) = -sxyz(2)
					!	Write all spacecraft instruments to .dat file
					write(scout+n,*) ut, sxyz, scdata
					!
					ntimes(n,1) = ntimes(n,1) + 1
					!	Check if this spacecraft has finished recording
					if(ntimes(n,1).ge.ntimes(n,2)) then
						recording = .false.
					else
						xcraft(:,n) = zcraft(:,n)
						zcraft(:,n) = craftpos(:,n,ntimes(n,1))
						!	Change direction to get from GSM to simulation coords
						zcraft(1,n)=-zcraft(1,n)
						zcraft(2,n)=-zcraft(2,n)
					endif	!end if(update zcraft)
				endif	!end if(time to record)
            enddo
            !
			if(craft_input) then
		    	!	Perhaps we should call limcraft on the craftpos values when we read them in.
				!	Then, we won't need to call limcraft repeatedly during time stepping.	MJS 08/13/17
		        call limcraft(xcraft,ncraft,re_equiv,n_grids, &
		            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
		            grid_minvals(3,:), grid_maxvals(3,:) )
		        !
		        srho=0.
		        do k=1,nz
		            do j=1,ny
		                dut=(ut-utold)/(future(j,k)-utold)
		                bxp(j,k)=bxp(j,k)+dut*(bxf(j,k)-bxp(j,k))
		                byp(j,k)=byp(j,k)+dut*(byf(j,k)-byp(j,k))
		                bzp(j,k)=bzp(j,k)+dut*(bzf(j,k)-bzp(j,k))
		                rhop(j,k)=rhop(j,k)+dut*(rhof(j,k)-rhop(j,k))
		                svxp(j,k)=svxp(j,k)+dut*(svxf(j,k)-svxp(j,k))
		                svyp(j,k)=svyp(j,k)+dut*(svyf(j,k)-svyp(j,k))
		                svzp(j,k)=svzp(j,k)+dut*(svzf(j,k)-svzp(j,k))
		                srho=srho+rhop(j,k)
		            enddo
		        enddo
	            !
        	    dut=(ut-utold)/(vut-utold)
		        svelx=svelx+dut*(zvelx-svelx)
		        svely=0.
		        srho=srho/float(nz*ny)
		        !
		        svelx=svelx+delvx_wind*delt
		        svely=svely+delvy_wind*delt
		        srho=srho+delrho*delt
				!
				svelz=svelz+delvz_wind*delt
				!
				sbx_wind=sbx_wind+delbx_wind*delt
				sby_wind=sby_wind+delby_wind*delt
				sbz_wind=sbz_wind+delbz_wind*delt
				!
				spx=srho*svelx
				spy=srho*svely
				spz=srho*svelz
				!
				spress=(cs_wind**2*srho/gamma)/gamma1
				serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
				!
				delay=t_equiv*distance/svelx/3600.
			endif	!end if(craft_input)
		endif	!end if(spacecraft)
		!
        !     *******************************
        !     Main grid loop over delt/m_step
        !     *******************************
        !
        do ms=1,m_step
            do box=n_grids,1,-1
                !
                !	Test if grid sector needs to be moved in time
				!
                yes_step=.false.
                !
                if(box.eq.1) then
                    yes_step=.true.
                else
                    if ((t_old(box).eq.t_new(box)).and.(t_new(box).lt.t_step(n_grids)) &
                        .and.(abs(t_new(box)-t_new(1)).le.0.005*t_step(1)) ) then
                        yes_step=.true.
                    endif
                endif
                !
                !	Time step grid
                !
                if(yes_step) then
                    t_old(box)=t_new(box)
                    t_new(box)=t_old(box)+t_step(box)
                    !
                    write(6,*)'lores t step. box/t/old/new: ',box, t, t_old(box),t_new(box)
                    !
                    delt = t_step(box)
                    delt2=delt/2.
                    !
                    if(tilting) then
                        atilt=old_tilt+(t_old(box)+delt2)*dtilt
                        sin_tilt=sin(atilt*pi/180.)
                        cos_tilt=cos(atilt*pi/180.)
                        ut2=utold+(t_old(box)+delt2)*t_equiv/3600.
                        nrot2=ut2/planet_per
                        rot_hr2=ut2-nrot2*planet_per
                        rot_angle=2.*pi*rot_hr2/planet_per
                        !        write(*,*)'mak_dip half with', t_old(box),delt2
                        call mak_dip_grd(bx0,by0,bz0,nx,ny,nz, &
                            n_grids,mbndry,box,ijzero,numzero,mzero,r_inner, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    !
                    !    *******************************
                    !    Store initial plasma parameters
                    !    *******************************
					!
					call store_array(oldqrho,qrho,nx,ny,nz,n_grids,box)
                    call store_array(oldqpresx,qpresx,nx,ny,nz,n_grids,box)
                    call store_array(oldqpresy,qpresy,nx,ny,nz,n_grids,box)
                    call store_array(oldqpresz,qpresz,nx,ny,nz,n_grids,box)
                    call store_array(oldqpresxy,qpresxy,nx,ny,nz,n_grids,box)
                    call store_array(oldqpresxz,qpresxz,nx,ny,nz,n_grids,box)
                    call store_array(oldqpresyz,qpresyz,nx,ny,nz,n_grids,box)
                    call store_array(oldqpx,qpx,nx,ny,nz,n_grids,box)
                    call store_array(oldqpy,qpy,nx,ny,nz,n_grids,box)
                    call store_array(oldqpz,qpz,nx,ny,nz,n_grids,box)
                    !
                    call store_array(oldhrho,hrho,nx,ny,nz,n_grids,box)
                    call store_array(oldhpresx,hpresx,nx,ny,nz,n_grids,box)
                    call store_array(oldhpresy,hpresy,nx,ny,nz,n_grids,box)
                    call store_array(oldhpresz,hpresz,nx,ny,nz,n_grids,box)
                    call store_array(oldhpresxy,hpresxy,nx,ny,nz,n_grids,box)
                    call store_array(oldhpresxz,hpresxz,nx,ny,nz,n_grids,box)
                    call store_array(oldhpresyz,hpresyz,nx,ny,nz,n_grids,box)
                    call store_array(oldhpx,hpx,nx,ny,nz,n_grids,box)
                    call store_array(oldhpy,hpy,nx,ny,nz,n_grids,box)
                    call store_array(oldhpz,hpz,nx,ny,nz,n_grids,box)
                    !
                    call store_array(oldorho,orho,nx,ny,nz,n_grids,box)
                    call store_array(oldopresx,opresx,nx,ny,nz,n_grids,box)
                    call store_array(oldopresy,opresy,nx,ny,nz,n_grids,box)
                    call store_array(oldopresz,opresz,nx,ny,nz,n_grids,box)
                    call store_array(oldopresxy,opresxy,nx,ny,nz,n_grids,box)
                    call store_array(oldopresxz,opresxz,nx,ny,nz,n_grids,box)
                    call store_array(oldopresyz,opresyz,nx,ny,nz,n_grids,box)
                    call store_array(oldopx,opx,nx,ny,nz,n_grids,box)
                    call store_array(oldopy,opy,nx,ny,nz,n_grids,box)
                    call store_array(oldopz,opz,nx,ny,nz,n_grids,box)
                    !
                    call store_array(oldepres,epres,nx,ny,nz,n_grids,box)
                    call store_array(oldbx,bx,nx,ny,nz,n_grids,box)
                    call store_array(oldby,by,nx,ny,nz,n_grids,box)
                    call store_array(oldbz,bz,nx,ny,nz,n_grids,box)
                    !
                    !
                    !	**************************************
                    !	2-step Lax-Wendroff: Step 1
                    !		Estimate fluid quantites at n+1/2
                    !	**************************************
                    !
                    !	Calculate standard MHD current: j = curl B
                    !
                    rx=xspac(box)
                    ry=xspac(box)
                    rz=xspac(box)
                    !
                    !write(*,*)'Entering subroutine: calcur'
                    call calcur(bx,by,bz,nx,ny,nz,n_grids,box,curx,cury,curz, &
                        rx,ry,rz)
                    !
                    !     Find total magnetic field
                    !
                    !write(*,*)'Entering subroutine: totbfld'
                    call totfld(bx,bx0,bsx,nx,ny,nz,n_grids,box)
                    call totfld(by,by0,bsy,nx,ny,nz,n_grids,box)
                    call totfld(bz,bz0,bsz,nx,ny,nz,n_grids,box)
                    !
                    !     Find magnitude of B
                    !
                    call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
                    !
                    !write(*,*)'Entering subroutine: fnd_evel'
                    call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
                        opx,opy,opz,orho,curx,cury,curz,evx,evy,evz, &
                        tvx,tvy,tvz,nx,ny,nz,n_grids,box, &
                        rmassq,rmassh,rmasso,reynolds)
                    !
                    !write(*,*)'Entering subroutine: bande'
                    call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
                        curx,cury,curz,evx,evy,evz,btot, &
                        epres,qrho,hrho,orho,resistive,resist,reynolds, &
                        nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso, &
                        ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
                        rx,ry,rz)
                    !
					!write(*,*)'Entering subroutine: push_elec, first pass'
                    call push_elec(wrkepres,oldepres,epres,evx,evy,evz, &
                        gamma,gamma1,nx,ny,nz,n_grids,box,0.5*delt, &
                        rx,ry,rz)
                    !
                    !write(*,*)'Entering subroutine: push_ion (q), first pass'
                    call push_ion(wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                        oldqpx,oldqpy,oldqpz, &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        oldqpresxy,oldqpresxz,oldqpresyz, &
                        qpresxy,qpresxz,qpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,n_grids,box,0.5*delt,grav,re_equiv,reynolds, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                        grid_minvals(3,:), grid_maxvals(3,:), ani_q,isotropic)
					!
                    !write(*,*)'Entering subrout: push_ion (h), first pass'
                    call push_ion(wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                        oldhpx,oldhpy,oldhpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        oldhpresxy,oldhpresxz,oldhpresyz, &
                        hpresxy,hpresxz,hpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,n_grids,box,0.5*delt,grav,re_equiv,reynolds, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                        grid_minvals(3,:), grid_maxvals(3,:), ani_h,isotropic)
                    !
                    !write(*,*)'Entering subroutine: push_ion (o), first pass'
                    call push_ion(wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        oldorho,oldopresx,oldopresy,oldopresz, &
                        oldopx,oldopy,oldopz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        oldopresxy,oldopresxz,oldopresyz, &
                        opresxy,opresxz,opresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,n_grids,box,0.5*delt,grav,re_equiv,reynolds, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                        grid_minvals(3,:), grid_maxvals(3,:), ani_o,isotropic)
                    !
                    !write(*,*)'Entering subroutine: push_bfld'
                    call push_bfld(wrkbx,wrkby,wrkbz,oldbx,oldby,oldbz, &
                        efldx,efldy,efldz,nx,ny,nz,n_grids,box,0.5*delt, &
                        rx,ry,rz)
                    !
                    !    *****************************************
                    !		     Apply boundary conditions
                    !    *****************************************
					!
					!write(*,*)'Entering first main Lax-Wendroff loop.'
                    !
                    if(box.eq.n_grids) then
	                    !write(*,*)'Applying boundary conditions.'
                        call bndry_outer( &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz,wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rmassq,rmassh,rmasso,wrkbx,wrkby,wrkbz, &
                            bx0,by0,bz0,nx,ny,nz,n_grids, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)
                    else
                        t_grid=t_old(box)+delt2
                        if(t_grid.gt.t_new(box+1)) then
                            t_grid=t_new(box+1)
                        endif
                        if (t_grid.lt.t_old(box+1)) then
                            t_grid=t_old(box+1)
                        endif
						!
                        !write(*,*) 'Entering subroutine: bndry_flanks'
                        call bndry_flanks( &
                            wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
                            wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
                            wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkorho,wrkopx,wrkopy,wrkopz, &
                            wrkopresx,wrkopresy,wrkopresz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            wrkepres,wrkbx,wrkby,wrkbz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres,oldbx,oldby,oldbz,vvx, &
                            nx,ny,nz,n_grids,box,t_old,t_new,t_grid, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    if(box.le.mbndry) then
                        !write(*,*)'Entering subroutine: bndry_inner'
                        call bndry_inner( &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz,rmassq, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz,rmassh, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz,rmasso, &
                            wrkepres,wrkbx,wrkby,wrkbz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            nx,ny,nz,n_grids,parm_srf,parm_mid,parm_zero, &
                            ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
                            mbndry,msrf,mmid,mzero, &
                            erho,epress,alpha_e,ti_te,o_conc, &
                            re_equiv,r_inner,sbx_wind,spress, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    !
                    !      write(*,*)'lax 1 speeds'
                    !
                    !	Check fluid parameters
                    !
                    if(craft_input) then
                        call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp, &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz, &
                            rmassq,rmassh,rmasso,wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,n_grids)
                    endif
					!
                    call set_rho(wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz,rmassq, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz,rmassh, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopresxy,wrkopresxz,wrkopresyz,rmasso, &
                        wrkepres,nx,ny,nz,n_grids,box,o_conc)
                    !
                    vlim=0.6*box
                    call set_speed_agrd( &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        wrkepres,wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        wrkbx,wrkby,wrkbz,bx0,by0,bz0, &
                        bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,n_grids,box, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest,isotropic)
                    !
                    !
                    !	***********************************************************
                    !	Lax-Wendroff: Step 2
                    !		Use the predicted value to find corrected value for n+1
                    !	***********************************************************
                    !
                    if(tilting) then
                        atilt=old_tilt+(t_old(box)+delt)*dtilt
                        sin_tilt=sin(atilt*pi/180.)
                        cos_tilt=cos(atilt*pi/180.)
                        ut2=utold+(t_old(box)+delt)*t_equiv/3600.
                        nrot2=ut2/planet_per
                        rot_hr2=ut2-nrot2*planet_per
                        rot_angle=2.*pi*rot_hr2/planet_per
                        !         write(*,*)'mak_dip full with',t_old(box),delt
                        call mak_dip_grd(bx0,by0,bz0,nx,ny,nz,n_grids, &
                            mbndry,box,ijzero,numzero,mzero,r_inner, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
						!
                    endif
                    rx=xspac(box)
                    ry=xspac(box)
                    rz=xspac(box)
                    !
                    !     Calculate standard MHD current: j = curl B
                    !
                    call calcur(wrkbx,wrkby,wrkbz,nx,ny,nz,n_grids,box,curx,cury,curz, &
                        rx,ry,rz)
                    !
                    !     Find total magnetic field
                    !
                    call totfld(wrkbx,bx0,bsx,nx,ny,nz,n_grids,box)
                    call totfld(wrkby,by0,bsy,nx,ny,nz,n_grids,box)
                    call totfld(wrkbz,bz0,bsz,nx,ny,nz,n_grids,box)
                    !
                    !     Find magnitude of B
                    !
                    call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
                    !
                    !
                    !     Find the electric field from electron momentum eqn
                    !
                    call fnd_evel(wrkqpx,wrkqpy,wrkqpz,wrkqrho, &
                        wrkhpx,wrkhpy,wrkhpz,wrkhrho, &
                        wrkopx,wrkopy,wrkopz,wrkorho, &
                        curx,cury,curz,evx,evy,evz, &
                        tvx,tvy,tvz,nx,ny,nz,n_grids,box, &
                        rmassq,rmassh,rmasso,reynolds)
                    !
                    call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
                        curx,cury,curz,evx,evy,evz,btot, &
                        wrkepres,wrkqrho,wrkhrho,wrkorho, &
                        resistive,resist,reynolds, &
                        nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso, &
                        ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
                        rx,ry,rz)
                    !
                    !write(*,*)'Entering subroutine: push_*, second pass'
                    call push_elec(epres,oldepres,wrkepres,evx,evy,evz, &
                        gamma,gamma1,nx,ny,nz,n_grids,box,delt, &
                        rx,ry,rz)
                    !
                    call push_ion(qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                        oldqpx,oldqpy,oldqpz, &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        qpresxy,qpresxz,qpresyz, &
                        oldqpresxy,oldqpresxz,oldqpresyz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,n_grids,box,delt,grav,re_equiv,reynolds, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                        grid_minvals(3,:), grid_maxvals(3,:), ani_q,isotropic)
                    !
					call push_ion(hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                        oldhpx,oldhpy,oldhpz, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        hpresxy,hpresxz,hpresyz, &
                        oldhpresxy,oldhpresxz,oldhpresyz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,n_grids,box,delt,grav,re_equiv,reynolds, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                        grid_minvals(3,:), grid_maxvals(3,:), ani_h,isotropic)
                    !
                    call push_ion(orho,opresx,opresy,opresz,opx,opy,opz, &
                        oldorho,oldopresx,oldopresy,oldopresz, &
                        oldopx,oldopy,oldopz, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        opresxy,opresxz,opresyz, &
                        oldopresxy,oldopresxz,oldopresyz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,n_grids,box,delt,grav,re_equiv,reynolds, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                        grid_minvals(3,:), grid_maxvals(3,:), ani_o,isotropic)
                    !
                    call push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
                        efldx,efldy,efldz,nx,ny,nz,n_grids,box,delt, &
                        rx,ry,rz)
                    !
                    !    *****************************************
                    !		     Apply boundary conditions
                    !    *****************************************
                    !
                    !write(*,*)'Entering second main Lax-Wendroff loop.'
                    !
                    if(box.eq.n_grids) then
                        call bndry_outer( &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rmassq,rmassh,rmasso,bx,by,bz, &
                            bx0,by0,bz0,nx,ny,nz,n_grids, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)
                    else
                        t_grid=t_old(box)+delt
                        if(t_grid.gt.t_new(box+1)) then
                            !        write(*,*)'warning on lax2 max',box,t_grid,t_new(box+1),
                            !    +                         t_grid-t_new(box+1)
                            t_grid=t_new(box+1)
                        endif
                        if (t_grid.lt.t_old(box+1)) then
                            t_grid=t_old(box+1)
                        endif
						!
                        call bndry_flanks( &
                            qrho,qpx,qpy,qpz, &
                            qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz, &
                            hpresx,hpresy,hpresz,hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz, &
                            opresx,opresy,opresz,opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres,oldbx,oldby,oldbz,vvx, &
                            nx,ny,nz,n_grids,box,t_old,t_new,t_grid, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    !
                    if(box.le.mbndry) then
                        call bndry_inner( &
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
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    !
                    if(craft_input) then
                        call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            rmassq,rmassh,rmasso,epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,n_grids)
                    endif
                    !
                    !      write(*,*)'lax 2 speeds'
                    !
                    call set_rho(qrho,qpresx,qpresy,qpresz, &
                        qpresxy,qpresxz,qpresyz,rmassq, &
                        hrho,hpresx,hpresy,hpresz, &
                        hpresxy,hpresxz,hpresyz,rmassh, &
                        orho,opresx,opresy,opresz, &
                        opresxy,opresxz,opresyz,rmasso, &
                        epres,nx,ny,nz,n_grids,box,o_conc)
                    !
                    vlim=0.6*box
                    call set_speed_agrd( &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        epres,qpresxy,qpresxz,qpresyz, &
                        hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
                        bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,n_grids,box, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest,isotropic)
					!
                    !     ...........................................................
                    !     Try Lapidus smoothing - smoothed results will appear in nt2
                    !     ...........................................................
                    !
                    !     write(*,*)' doing smoothing okay'
                    rx=xspac(box)
                    ry=xspac(box)
                    rz=xspac(box)
                    !      write(*,*)'calling lapidus'
                    !
                    !     Lapidus smoothing
                    !
                    !      species q
                    !
                    !     write(*,*)'lap ion1'
                    call fnd_vel(qpx,qpy,qpz,qrho, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box)
                    call lap_plasma(qrho,qpx,qpy,qpz, &
                        qpresx,qpresy,qpresz, &
                        qpresxy,qpresxz,qpresyz, &
                        wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
                        wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box, &
                        chirho,chipxyz,chierg,delt,isotropic)
                    !
                    !     species h
                    !
                    !     write(*,*)'lap ion2'
                    call fnd_vel(hpx,hpy,hpz,hrho, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box)
                    call lap_plasma(hrho,hpx,hpy,hpz, &
                        hpresx,hpresy,hpresz, &
                        hpresxy,hpresxz,hpresyz, &
                        wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
                        wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box, &
                        chirho,chipxyz,chierg,delt,isotropic)
                    !
                    !     species o
                    !
                    !     write(*,*)'lap ion3'
                    call fnd_vel(opx,opy,opz,orho, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box)
                    call lap_plasma(orho,opx,opy,opz, &
                        opresx,opresy,opresz, &
                        opresxy,opresxz,opresyz, &
                        wrkorho,wrkopx,wrkopy,wrkopz, &
                        wrkopresx,wrkopresy,wrkopresz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box, &
                        chirho,chipxyz,chierg,delt,isotropic)
                    !
                    !     electrons
                    !
                    !     write(*,*)'lap electrons'
                    call fnd_vtot(qpx,qpy,qpz,qrho, &
                        hpx,hpy,hpz,hrho, &
                        opx,opy,opz,orho, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box, &
                        rmassq,rmassh,rmasso)
                    call lap_elec(epres,wrkepres, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box, &
                        chirho,chipxyz,chierg,delt)
                    !
                    !     magnetic fields
                    !
                    !     write(*,*)'lap bfld'
                    call lap_bfld(bx,by,bz, &
                    wrkbx,wrkby,wrkbz, &
                        curx,cury,curz, &
                        nx,ny,nz,n_grids,box, &
                        chirho,chipxyz,chierg,delt)
                    !
                    !     reset boundary conditions
                    !
                    if(box.eq.n_grids) then
                        call bndry_outer( &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz, &
                            wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rmassq,rmassh,rmasso,wrkbx,wrkby,wrkbz, &
                            bx0,by0,bz0,nx,ny,nz,n_grids, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)
                    else
                        call bndry_flanks( &
                            wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
                            wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
                            wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkorho,wrkopx,wrkopy,wrkopz, &
                            wrkopresx,wrkopresy,wrkopresz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            wrkepres,wrkbx,wrkby,wrkbz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres,oldbx,oldby,oldbz, vvx, &
                            nx,ny,nz,n_grids,box,t_old,t_new,t_grid, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    if(box.le.mbndry) then
                        !       write(*,*)'calling bndry inner lap'
                        call bndry_inner( &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz,rmassq, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz,rmassh, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz,rmasso, &
                            wrkepres,  wrkbx, wrkby, wrkbz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            nx,ny,nz,n_grids,parm_srf,parm_mid,parm_zero, &
                            ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
                            mbndry,msrf,mmid,mzero, &
                            erho,epress,alpha_e,ti_te,o_conc, &
                            sbx_wind,spress, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    !
                    if(craft_input) then
                        call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp, &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz, &
                            rmassq,rmassh,rmasso,wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,n_grids)
                    endif
                    call set_rho(wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz,rmassq, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz,rmassh, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopresxy,wrkopresxz,wrkopresyz,rmasso, &
                        wrkepres,nx,ny,nz,n_grids,box,o_conc)
                    !
                    !      write(*,*)'lapidus speeds'
                    !
                    vlim=0.6*box
                    call set_speed_agrd( &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        wrkepres,wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        wrkbx,wrkby,wrkbz,bx0,by0,bz0, &
                        bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,n_grids,box, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest,isotropic)
                    !
                    !write(*,*) 'Lapidus smoothing complete.'
                    !
                    !     ....................................
                    !     Add a little bit of flux correction:
                    !     ....................................
                    !
                    !write(*,*)'Entering subroutine: flux_correct'
                    call flux_correct(qrho,qpresx,qpresy,qpresz, &
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
                        nx,ny,nz,n_grids,box,difrho,diferg,xspac, &
                        isotropic)
                    !
                    !      Set boundary conditions
                    !
                    if(box.eq.n_grids) then
                        call bndry_outer( &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rmassq,rmassh,rmasso,bx,by,bz, &
                            bx0,by0,bz0,nx,ny,nz,n_grids, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te,isotropic)
                    else
                        call bndry_flanks( &
                            qrho,qpx,qpy,qpz, &
                            qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz, &
                            hpresx,hpresy,hpresz,hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz, &
                            opresx,opresy,opresz,opresxy,opresxz,opresyz, &
                            epres, &
                            bx,by,bz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres, &
                            bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres, &
                            oldbx,oldby,oldbz, vvx, &
                            nx,ny,nz,n_grids,box,t_old,t_new,t_grid, &
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                    endif
                    !
                    if(box.le.mbndry) then
                        !       write(*,*)'calling bndry inner'
                        call bndry_inner( &
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
                            grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                            grid_minvals(3,:), grid_maxvals(3,:) )
                        !
                    endif   ! end mbndry
           			!
                    if(craft_input) then
                        call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
                        	qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            rmassq,rmassh,rmasso,epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,n_grids)
                    endif
                    !
                    !	Write(*,*)'fcsmooth speeds'
                    !
                    call set_rho(qrho,qpresx,qpresy,qpresz, &
                        qpresxy,qpresxz,qpresyz,rmassq, &
                        hrho,hpresx,hpresy,hpresz, &
                        hpresxy,hpresxz,hpresyz,rmassh, &
                        orho,opresx,opresy,opresz, &
                        opresxy,opresxz,opresyz,rmasso, &
                        epres,nx,ny,nz,n_grids,box,o_conc)
                    !
                    vlim=0.6*box
                    call set_speed_agrd( &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        epres,qpresxy,qpresxz,qpresyz, &
                        hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
                        bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,n_grids,box, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest,isotropic)
					!
                    t_stepnew(box)=stepsz*xspac(box)/fastest
                    !     write(*,*)'needed step of',t_step(box),t_stepnew(box)
                    !
                    !     sync time steps if needed and apply core conditions
                    !
                    if(box.lt.n_grids) then
                        if(abs(t_new(box)-t_new(box+1)).le.0.005*t_step(box)) then
                            call bndry_corer( &
                                qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                                hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                                orho,opresx,opresy,opresz,opx,opy,opz, &
                                epres,bx,by,bz, &
                                qpresxy,qpresxz,qpresyz, &
                                hpresxy,hpresxz,hpresyz, &
                                opresxy,opresxz,opresyz, &
                                nx,ny,nz,n_grids,box, &
                                grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                                grid_minvals(3,:), grid_maxvals(3,:) )
                            !         write(*,*)'corer',box,t_new(box),box+1,t_new(box+1)
                            t_old(box+1)=t_new(box+1)
                        endif
                    endif
                    !
                    !
                    !
                endif	!	yes_step of main grid
            enddo	!	lores box increment
        enddo	!	box sweep of lores grid            
        !
        !write(*,*) 'Final sync on boundary conditions.'
        t_old(1)=t_new(1)
        !
        do box=1,n_grids-1
            call bndry_corer( &
                qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                orho,opresx,opresy,opresz,opx,opy,opz, &
                epres,bx,by,bz, &
                qpresxy,qpresxz,qpresyz, &
                hpresxy,hpresxz,hpresyz, &
                opresxy,opresxz,opresyz, &
                nx,ny,nz,n_grids,box, &
                grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                grid_minvals(3,:), grid_maxvals(3,:) )
            !
            call set_rho(qrho,qpresx,qpresy,qpresz, &
                qpresxy,qpresxz,qpresyz,rmassq, &
                hrho,hpresx,hpresy,hpresz, &
                hpresxy,hpresxz,hpresyz,rmassh, &
                orho,opresx,opresy,opresz, &
                opresxy,opresxz,opresyz,rmasso, &
                epres,nx,ny,nz,n_grids,box,o_conc)
                t_old(box)=t_old(box-1)
                t_new(box)=t_new(box-1)
            !
        enddo
		!
        !     write(*,*)'end loop 999'
        !
        !     write(*,999)t
        ! 999 format(' step 2 complete at t= ',1pe12.5)
        !
        if(t.ge.tgraph) then
            !
            !	Plot plasma properties
            !
            !	Calculate size of plotting stuff and ensure no distortions
            !		over desired scales
            !
            call visual(qrho,qpresx,qpresy,qpresz,qpresxy, &
                qpresxz,qpresyz,qpx,qpy,qpz,rmassq, &
                hrho,hpresx,hpresy,hpresz,hpresxy, &
                hpresxz,hpresyz,hpx,hpy,hpz,rmassh, &
                orho,opresx,opresy,opresz,opresxy, &
                opresxz,opresyz,opx,opy,opz,rmasso, &
                epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
                curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
                tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
                nx,ny,nz,n_grids,xspac, &
                cross,along,flat,xcraft,ncraft,re_equiv, &
                grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
                grid_minvals(3,:), grid_maxvals(3,:), ut,b_equiv,ti_te,rho_equiv)
            !
            tgraph=tgraph+deltg
			write(*,*)'Graphics plotted.'
        endif !if(t.ge.tgraph)
		!
        if(t.ge.tinj.and.ringo) then
            !
            !	Inject torus plasma
            !
            write(*,*) 'Injecting plasma torus, tinj:', tinj
			!
            box=1
            dx = ( grid_maxvals(1,box) - grid_minvals(1,box) ) / (nx-1.)
            dy = ( grid_maxvals(2,box) - grid_minvals(2,box) ) / (ny-1.)
            dz = ( grid_maxvals(3,box) - grid_minvals(3,box) ) / (nz-1.)
            !
            tot_o=0.
            tot_h=0.
            tot_q=0.
            !
            call totfld(bx,bx0,bsx,nx,ny,nz,n_grids,box)
            call totfld(by,by0,bsy,nx,ny,nz,n_grids,box)
            call totfld(bz,bz0,bsz,nx,ny,nz,n_grids,box)
            !
            do k=1,nz
                az = grid_minvals(3,box) + dz * (k-1) - zdip
                do j=1,ny
                    ay = grid_minvals(2,box) + dy * (j-1) - ydip
                    do i=1,nx
                        ax = grid_minvals(1,box) + dx * (i-1) - xdip
                        !
                        rx=ax*re_equiv
                        ry=ay*re_equiv
                        rd=sqrt(rx**2+ry**2)
						ar=sqrt( (ax+xdip)**2 + (ay+ydip)**2 + (az+zdip)**2 )	!	Radial distance from planet center
                        !
                        rvy=rx*v_rot
                        rvx=-ry*v_rot
                        rvz=0.
                        corotate=sqrt(rvx**2+rvy**2)
                        !
                        ar_tmid=sqrt(ax**2+ay**2)*re_equiv
                        dr_tmid=abs(ar_tmid - torus_inj_dist)
  						!
                        ! MAT - Corrected, according to 
                        !   doi:10.1029/2009JE003372
                        !   and
                        !   doi:10.1029/JA091iA08p08749
                        !
                        if( (abs(dr_tmid).lt.2.*torus_rad) .and. (ar.gt.srf_bndry) ) then
                            ! scale height in # planet_rad
                            rscale=exp(-((dr_tmid)/(0.5*torus_rad))**2) 
                            zscale=exp(-((az*re_equiv)/(0.5*torus_rad))**2)
                            abtot=(bsx(i,j,k)**2+bsy(i,j,k)**2)/bsz(i,j,k)**2
                            dscale=amax1(0.,1.-10.*abtot)
                            !
                            ! Inject h species
                            !
                            hden=denh_torus*rmassh*rscale*zscale*dscale
                            hrho(i,j,k,box)=hrho(i,j,k,box)+hden
                            !temp goes as rho*v**2
                            del_hp=hden*(corotate**2)*t_torus
                            hpresx(i,j,k,box)=hpresx(i,j,k,box)+del_hp
                            hpresy(i,j,k,box)=hpresy(i,j,k,box)+del_hp
                            hpresz(i,j,k,box)=hpresz(i,j,k,box)+del_hp*aniso_factor
                            hpresxy(i,j,k,box)=hpresxy(i,j,k,box)+del_hp
							!
                            hpx(i,j,k,box)=hpx(i,j,k,box)+reduct*hden*rvx
                            hpy(i,j,k,box)=hpy(i,j,k,box)+reduct*hden*rvy
                            !
                            ! Inject q species
                            !
                            qden=denq_torus*rmassq*rscale*zscale*dscale 
                            qrho(i,j,k,box)=qrho(i,j,k,box)+qden
                            !temp goes as v**2
                            del_qp=(qden/rmassq)*(corotate**2)*t_torus
                            qpresx(i,j,k,box)=qpresx(i,j,k,box)+del_qp
                            qpresy(i,j,k,box)=qpresy(i,j,k,box)+del_qp
                            qpresz(i,j,k,box)=qpresz(i,j,k,box)+del_qp*aniso_factor
                            qpresxy(i,j,k,box)=qpresxy(i,j,k,box)+del_qp
                            !
                            qpx(i,j,k,box)=qpx(i,j,k,box)+reduct*qden*rvx
                            qpy(i,j,k,box)=qpy(i,j,k,box)+reduct*qden*rvy
                            !
                            ! Inject o species 
                            !
                            oden=deno_torus*rmasso*rscale*zscale*dscale 
                            orho(i,j,k,box)=orho(i,j,k,box)+oden
                            !	temp goes as v**2
                            del_op=(oden/rmasso)*(corotate**2)*t_torus    
                            opresx(i,j,k,box)=opresx(i,j,k,box)+del_op
                            opresy(i,j,k,box)=opresy(i,j,k,box)+del_op
                            opresz(i,j,k,box)=opresz(i,j,k,box)+del_op*aniso_factor
                            opresxy(i,j,k,box)=opresxy(i,j,k,box)+del_op
							!
                            opx(i,j,k,box)=opx(i,j,k,box)+reduct*oden*rvx
                            opy(i,j,k,box)=opy(i,j,k,box)+reduct*oden*rvy
                            !	equal temps
                            epres(i,j,k,box)=epres(i,j,k,box)+(del_op+del_hp+del_qp)
                            !
                            tot_q=tot_q+qden*dx*dy*dz
                            tot_h=tot_h+hden*dx*dy*dz
                            tot_o=tot_o+oden*dx*dy*dz
                            !
                        endif
                    enddo
                enddo
            enddo
            !
            !	Scale factors to kg/s
            !
            volume=(re_equiv*planet_rad*1.e3)**3  !(cubic meters)
            atime=deltinj*t_equiv
            proton_mass=1.67e-27
            write(*,*)'ut,volume,t_equiv,atime',ut,volume,t_equiv,atime
            !
            tot_q=tot_q*volume/atime*rho_equiv*1.e6/rmassq
            tot_h=tot_h*volume/atime*rho_equiv*1.e6/rmassh
            tot_o=tot_o*volume/atime*rho_equiv*1.e6/rmasso
            write(*,*)'tot torus ions/s q,h,o',ut,tot_q,tot_h,tot_o
            write(10,*)'tot torus ions/s q,h,o',ut,tot_q,tot_h,tot_o
            !
            tot_q=tot_q*proton_mass*rmassq
            tot_h=tot_h*proton_mass*rmassh
            tot_o=tot_o*proton_mass*rmasso
            write(*,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
            write(10,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
            !
            tinj=tinj+deltinj
			!
        endif !if(t.ge.tinj.and.ringo)
		!
        if(t.ge.ts1) then
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
			write(*,*) 'Writing to fluid file. t =', t
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
            if(nchf.eq.11) &
                open(11,file='fluid11',status='unknown',form='unformatted')
            if(nchf.eq.12) &
                open(12,file='fluid12',status='unknown',form='unformatted')
            if(nchf.eq.13) &
                open(13,file='fluid13',status='unknown',form='unformatted')
            if(nchf.eq.14) &
                open(14,file='fluid14',status='unknown',form='unformatted')
            if(nchf.eq.15) &
                open(15,file='fluid15',status='unknown',form='unformatted')
            !
            !	Write restart data
            !
            write(nchf)t
            write(nchf)qrho
            write(nchf)qpx
            write(nchf)qpy
            write(nchf)qpz
            write(nchf)qpresx
            write(nchf)qpresy
            write(nchf)qpresz
            write(nchf)qpresxy
            write(nchf)qpresxz
            write(nchf)qpresyz
            write(nchf)hrho
            write(nchf)hpx
            write(nchf)hpy
            write(nchf)hpz
            write(nchf)hpresx
            write(nchf)hpresy
            write(nchf)hpresz
            write(nchf)hpresxy
            write(nchf)hpresxz
            write(nchf)hpresyz
            write(nchf)orho
            write(nchf)opx
            write(nchf)opy
            write(nchf)opz
            write(nchf)opresx
            write(nchf)opresy
            write(nchf)opresz
            write(nchf)opresxy
            write(nchf)opresxz
            write(nchf)opresyz
            write(nchf)bx
            write(nchf)by
            write(nchf)bz
            write(nchf)epres
            write(nchf)bx0
            write(nchf)by0
            write(nchf)bz0
            write(nchf)parm_srf,parm_mid,parm_zero, &
                ijzero,numzero,ijmid,nummid,ijsrf,numsrf
            close(nchf)
			!
            nchf=nchf+1
            if(nchf.gt.15)nchf=11
            ts1=ts1+tsave
            !
            if(divb_lores) then
                range=1.33*torus_inj_dist/re_equiv
                write(*,*)'Range for divb taper: ',torus_inj_dist,range
                do box=n_grids,1,-1
                    write(*,*)'divb on box: ',box
                    call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                        curx,cury,curz,efldx,efldy,efldz, &
                        7,nx*ny*nz,nx,ny,nz,n_grids,box,xspac)
                    call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                        curx,cury,curz,efldx,efldy,efldz, &
                        7,nx*ny*nz,nx,ny,nz,n_grids,box,xspac)
                    write(*,*)'Completed divb on box: ',box
                enddo	!end box loop
                !
                !	Apply boundary conditions
                !
                do box=n_grids-1,1,-1
                    call flanks_synced(bx,nx,ny,nz,n_grids,box, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:) )
                    call flanks_synced(by,nx,ny,nz,n_grids,box, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:) )
                    call flanks_synced(bz,nx,ny,nz,n_grids,box, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:) )
                enddo
                do box=1,n_grids-1
                    call bndry_corer( &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        epres,bx,by,bz, &
                        qpresxy,qpresxz,qpresyz, &
                        hpresxy,hpresxz,hpresyz, &
                        opresxy,opresxz,opresyz, &
                        nx,ny,nz,n_grids,box, &
                        grid_minvals(1,:), grid_maxvals(1,:), grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:) )
                enddo  !end bndry_corer
            endif  ! end divb_lores
            !
        endif	!end if(t.ge.ts1)
    enddo	!end do while(t.lt.tmax)
    call clsgks
	!
	utstart = ut
	write(2,option)
    write(2,planet)
    write(2,speeds)
    write(2,windy)
    write(2,lunar)
    write(2,physical)
    write(2,smooth)
	close(2)
	!
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) 'Run complete! t =', t
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) '-'
	call system('date +"FINISH: %T-%Y-%m-%d"')
	write(*,*) '-'
end
