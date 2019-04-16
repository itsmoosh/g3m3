!
!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!		@														@
!		@		GLOBAL 3D MAGNETOSPHERIC MULTIFLUID MODEL		@
!		@				(G3M^3) MAIN SEQUENCE					@
!		@														@
!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!
!
!	This is a 3D, three-fluid simulation using
!		electrons	: arrays starting with e
!		[solar] wind: arrays starting with q (protons)
!		ionospheric	: arrays starting with o oxygen,
!		h hydrogen
!
!	Warnings:
!		@	Add in mirror dipole so bx is zero at wind boundary
!		@	NCAR graphics - contour and flows - have to be manually set
!		for correct aspect ratios. Deprecated, as NCAR graphics is no
!		longer used.
!		@	Plasma and magnetic field data must be aligned in time
!
!	Grid within grid size nt = division*n_grids
!
!	File indices:
!	1-9		Top-level I/O
!	10-20	Fluid files
!	21-29	Primary output data
!	30-59	Spacecraft position/trajectory input files
!	60-89	Spacecraft data recording
!	100		git hash file
!	100-109	dummy files
!	110-150	Figure plotting
!
program multifluid
	
	!	*******
	!	Modules
	!	*******
	
	!	Contains planetary constants
	!	And some physical constants, including pi
	use astrometry
	
	implicit none

	!	*******************
	!	Critical parameters
	!	*******************

	!	Double precision
	integer, parameter :: dp = kind(1.d0)
	!	Small, arbitrary amount to prevent small denominators
	real, parameter :: smallbit = 1.e-5
	!	Tab character for output formatting
	character,parameter	:: tab=char(9)
	
	!	*******************
	!	Gridding parameters
	!	*******************

	real,parameter :: limit=60.	! Typically 60
	integer,parameter :: nx=2*limit+1	! Typically 121
	integer,parameter :: ny=2*limit+1
	integer,parameter :: nz=limit+1
	integer,parameter :: num_pts(3)=[nx,ny,nz]
	integer,parameter :: n_grids=5
	integer,parameter :: division=2	! Scaling between adjacent grids
	integer,parameter :: mbndry=1	!	Box for inner boundary
	integer,parameter :: ncts=281	!	Number of set_imf points
	!	wind_adjust * limit must yield a whole number
	real,parameter :: wind_adjust=5./4.
	real xspac(n_grids)
	real grid_minvals(3,n_grids), grid_maxvals(3,n_grids)
	real csmax, alfmax, pxmax, pymax, pzmax

	integer i,ix,iy,iz, j,k

	!	***********************
	!	Log file format strings
	!	***********************

	!	speeds_fmt: box, csmax, alfmax, pxmax, pymax, pzmax
	character*21, parameter :: speeds_fmt = '(1x,i2,5(1x,1pe12.5))'
	!	init_fmt: 'Init: t=' t, tab//'dt=', delt, tab//'ut=' ut
	character*37, parameter :: init_fmt = '(1X,A8,1pe12.5,A4,' &
		// '1pe12.5,A4,1pe12.5)'
	
	!	**************
	!	Grid variables
	!	**************
	
	!	Grid limits are set by grid_minvals grid_maxvals arrays.
	!	'limit' is the base min/max value for the smallest xy box.
	!	z is half as much.
	!	xspac is the relative grid spacing; the distance between grid
	!	points in units of the smallest grid step size.
	!	wind_adjust is the factor to stretch/compress the biggest box
	!	along the wind direction to capture the magnetotail.
	!	wind_adjust is the # of times larger to make the magnetotail
	!	side than the y min/max values.
	!	UT is in hours; we use higher precision so that we can
	!	track spacecraft times seconds apart after many hours of
	!	simulation time.
	!	UT = 0 is chosen to keep values relatively low for Europa sims,
	!	1996-12-02T14:24:00.000
	
	real,parameter	:: d_min=0.01	!	Minimum plasma density
	real,parameter	:: o_conc_min=0.01
	real,parameter	:: o_conc_max=1.0
	real,parameter	:: cur_min=0.75
	real,parameter	:: cur_max=20.0
	real,parameter	:: vlim_factor=0.6	!	Multiply by box to get max CFL speed
	!	Minimum time step, in sim units. If delt is less than this,
	!	exit with error because something has gone wrong.
	real,parameter	:: t_step_min = 1.e-6
	
	!	Scale lengths for plasma sheet population
	real,parameter	:: sheet_den=0.25
	real,parameter	:: alpha_s=4.
	
	!	Number of quantities to zero inside inner boundary
	integer,parameter :: num_zqt = 7
	!	Number of grid points to zero inside innermost boundary
	integer mzero
	!	Number of grid points in 'mid' region between inner boundary and
	!	zero region
	integer mmid
	!	Number of grid points on inner boundary surface
	integer msrf
	real,parameter	:: zero_bndry_offset = -1.5
	real,parameter	:: mid_bndry_offset = -0.5
	real,parameter	:: srf_bndry_offset = 0.6
	real zero_bndry
	real mid_bndry
	real srf_bndry
	real vlim
	integer,allocatable	:: ijsrf(:,:,:)
	integer,allocatable	:: ijmid(:,:,:)
	integer,allocatable	:: ijzero(:,:,:)
	integer	numsrf(mbndry)
	integer	nummid(mbndry)
	integer	numzero(mbndry)
	real,allocatable	:: parm_srf(:,:,:)
	real,allocatable	:: parm_mid(:,:,:)
	real,allocatable	:: parm_zero(:,:,:)
	
	!	********************
	!	Graphics parameters:
	!	********************

	!	Note: muvwp2=amax(mx,my,mz)+2,mz2=(mz-1)/2+1
	!	Deprecated; needed for NCAR graphics subroutines
	integer,parameter :: mx=61
	integer,parameter :: my=61
	integer,parameter :: mz=31
	integer,parameter :: muvwp2=63
	integer,parameter :: mz2=16
	
	!	*************************
	!	Coordinate work variables
	!	*************************

	!	Used for differences between points in various checks
	!	(stand-in for xspac)
	real dx,dy,dz
	!	Used in loops over grid points--x,y,z locations in re_equiv
	!	relative to origin
	real ax,ay,az
	!	Rotated coordinates
	real xr,yr,zr
	!	Dipole coordinates
	real xp,yp,zp
	!	Radial distance combining ax, ay, az
	real ar, ar2
	! 	?
	real ra
	!	Active loop values for velocity components
	real avx,avy,avz
	!	Distances relative to the origin
	real rx, ry, rz, rd
	!	Velocity variables derived from rotational motion
	real rvx,rvy,rvz
	!	Used as placeholders in various comparisons
	real apres, arho, astep, atilt, atime, old_tilt, told, ts1, ut2, tstep, &
		volume
	integer nsteps, nvx
	real hden, hflux_in, hflux_out, tot_h, &
		oden, oflux_in, oflux_out, ofrac, tot_o, &
		qden, qflux_in, qflux_out, tot_q
	!	Placeholders for input_fluid from upwind spacecraft
	real b_perp, bmag, delay, delt1, delt2, delt_old, dfrac, displace, &
		distance, dscale, dut, abtot, erho, dtilt, vut, rut, spd, wind_bnd
	
	!	Corotation speed for current distance
	real corotate
	!	Radial distance from origin to torus center in sim units
	real ar_tmid
	!	Distance from center of torus to current point
	real dr_tmid
	!	Corotation change for each time step
	real del_qp, del_hp, del_op
	real ar_iono
	real r_equat
	real rerg_sphere
	real rden_sphere
	real zheight
	real rho_iono
	real vt
	real scale, range, rscale, zscale

	!	Number of full rotations for planet and moon
	integer nrot_planet, nrot_moon, nrot2
	!	Number of hours into current local day/moon orbit
	real rot_hrs, rot_hrs_moon, rot_hr2
	real rot_angle, rot_angle_moon
	real sin_rot, cos_rot, sin_rot_moon, cos_rot_moon
	
	!	*********************
	!	Input file parameters
	!	*********************

	! group 'option':
	real tmax, stepsz, tsave
	integer ntgraph, ntinj
	logical start, isotropic, write_dat, diagnostics
	character*8 run_name
	! group 'planet'
	character*10 bodyname, moonname
	real xdip, ydip, zdip, r_inner, torus_rad, &
		tilt1, tilt2, &
		rmassq, rmassh, rmasso
	logical tilting
	! group 'speeds'
	real cs_inner, alf_inner1, alf_inner2, &
		alpha_e, denh_inner, denq_torus, denh_torus, deno_torus, &
		gravity, ti_te, gamma, reduct, t_torus, aniso_factor, aniso_limit, &
		qvlim, hvlim, ovlim, qclim, hclim, oclim, cslim
	logical ringo, update, reload, divb_lores, divb_hires
	! group 'windy'
	real re_wind, cs_wind, &
		vx_wind1, vx_wind2, vy_wind1, vy_wind2, vz_wind1, vz_wind2, &
		alfx_wind1, alfx_wind2, alfy_wind1, alfy_wind2, alfz_wind1, &
		alfz_wind2, den_wind1, den_wind2, &
		reynolds, resist, o_conc, rho_frac, bfrac, vfrac
	! group 'physical' 
	real re_equiv, v_equiv, rho_equiv
	real utstart
	logical spacecraft, input_fluid, warp, repeat_flybys
	! group 'smooth'
	real chirho, chipxyz, chierg, difrho, difpxyz, diferg
	! group 'crafthead'
	character*8	 :: cname = 'default'
	integer	:: num_vals = 0
	real	:: rot_closest = 0.0
	character*7 :: git_hash = 'pl_hold'
	
	!	*******************
	!	Spacecraft file I/O
	!	*******************

	!	File indices:
	!		1-9		Top-level I/O
	!		10-20	Fluid files
	!		21-29	Primary output data
	!		30-59	Spacecraft position/trajectory input files
	!		60-89	Spacecraft data recording
	!		100		Git hash file
	!		100-109	Dummy files
	!		110-150	Figure plotting

	integer,parameter :: input_f=1
	integer,parameter :: inp_o_f=2
	integer :: fluxs_f=21
	integer :: speed_f=22
	integer :: concs_f=23
	integer :: grdpt_f=24
	integer :: recdt_f=29
	integer :: fluid_f=10
	integer :: offset_f=-10
	integer,parameter :: scin=30
	integer,parameter :: scout=60
	integer,parameter :: git_f=100
	integer,parameter :: dummy_f=101
	integer,parameter :: dummy_fg=102
	
	character*5,parameter	:: fname_input='input'
	character*9,parameter	:: fname_inp_o='input_out'
	character*10,parameter	:: fname_fluxs='fluxes.dat'
	character*10,parameter	:: fname_speed='speeds.dat'
	character*8,parameter	:: fname_concs='conc.dat'
	character*7,parameter	:: fluid_pre='fluid'
	character*11,parameter	:: fname_restart='fluid_start'
	character*7		:: fname_reload=fluid_pre//'00'
	character*7		:: fname_fluid=''
	
	character*60,parameter	:: flux_header='main grid'//tab// &
		'UT'//tab//'box'//tab//'qflux'//tab//'hflux'//tab//'oflux'
	character*12,parameter	:: git_hash_file='git_hash.txt'
	
	!	************************************
	!	Physics plasma quantities: Main grid
	!	************************************

	real, dimension(nx,ny,nz,n_grids) :: bx,by,bz, qpx,qpy,qpz,qrho, &
		qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz, &
		hpx,hpy,hpz,hrho,hpresx,hpresy,hpresz,hpresxy,hpresxz,hpresyz, &
		opx,opy,opz,orho,opresx,opresy,opresz,opresxy,opresxz,opresyz, &
		epres
	
	!	****************************************************
	!	Work arrays for Runge-Kutta and smoothing: Main grid
	!	****************************************************

	real, dimension(nx,ny,nz,n_grids) :: &
		oldbx,oldby,oldbz, &
		oldqpx,oldqpy,oldqpz,oldqrho, &
		oldqpresx,oldqpresy,oldqpresz, &
		oldqpresxy,oldqpresxz,oldqpresyz, &
		oldhrho,oldhpx,oldhpy,oldhpz, &
		oldhpresx,oldhpresy,oldhpresz, &
		oldhpresxy,oldhpresxz,oldhpresyz, &
		oldorho,oldopx,oldopy,oldopz, &
		oldopresx,oldopresy,oldopresz, &
		oldopresxy,oldopresxz,oldopresyz, &
		oldepres
	
	real, dimension(nx,ny,nz,n_grids) :: &
		wrkbx,wrkby,wrkbz, &
		wrkqpx,wrkqpy,wrkqpz,wrkqrho, &
		wrkqpresx,wrkqpresy,wrkqpresz, &
		wrkqpresxy,wrkqpresxz,wrkqpresyz, &
		wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
		wrkhpresx,wrkhpresy,wrkhpresz, &
		wrkhpresxy,wrkhpresxz,wrkhpresyz, &
		wrkorho,wrkopx,wrkopy,wrkopz, &
		wrkopresx,wrkopresy,wrkopresz, &
		wrkopresxy,wrkopresxz,wrkopresyz, &
		wrkepres
	
	!	**********************
	!	Unperturbed quantities
	!	**********************

	real b_equiv
	real, dimension(nx,ny,nz,n_grids) :: bx0,by0,bz0
	
	real, dimension(nx,ny,nz) :: efldx,efldy,efldz, curx,cury,curz, &
		bsx,bsy,bsz,btot
	real resistive(nx,ny,nz,mbndry)
	
	!	*************************
	!	Variable time step arrays
	!	*************************
	
	!	We give t extra precision to make sure we don't easily hit an
	!	upper limit on runtime, and to have small steps for in-situ
	!	spacecraft measurements.
	real(dp) t, ut, utold
	real(dp) t_grid
	real(dp) t_old(n_grids)
	real(dp) t_new(n_grids)
	real(dp) t_step(n_grids)
	real(dp) t_stepnew(n_grids)
	real t_equiv
	!	For checking t and ut upon restart to warn of drastic changes
	real t_diff

	real, parameter :: t_step_max_def = 20.0
	real :: t_step_max = t_step_max_def
	
	!	*************************
	!	Boundary condition arrays
	!	*************************

	real, dimension(ny,nz) :: bxf,byf,bzf, rhof, svxf,svyf,svzf, &
		bxp,byp,bzp, rhop, svxp,svyp,svzp, future,past
	real bfld(ncts,4), rplas(ncts), svel(ncts,3)
	integer ncount(ny,nz), nc
	
	real, dimension(mx,my,mz) :: tx,ty,tz, tg1,tt
	real tg2(mx,my,mz2), work(muvwp2,muvwp2)
	real cross(ny,nz), along(nx,nz), flat(nx,ny)
	
	!	***************************
	!	Labels for graphics outputs
	!	***************************

	character*8 label
	character*15 title
	
	!	*************
	!	Miscellaneous
	!	*************

	!	Loop counters
	integer n, m, ms, nn, box
	!	Placeholders
	integer m_smallest, m_step, cbox
	!	For checking grid mating
	integer mating_index, next_box
	!	Keeps count of time ID for graphing data
	integer :: nplots=0
	!	For checking grid mating
	real grid_diff
	real smallest_step, fastest
	
	logical add_dip, save_dat, yes_step, yes_step_n, grid_reset
	
	!	************************
	!	Adjusted data parameters
	!	************************
	
	real gamma1
	real epress, eerg, spress, serg
	real rho_wind1, rho_wind2
	real srho, delrho
	real svelx, svely, svelz, zvelx, zvely, zvelz
	real spx, spy, spz
	
	real delvx_wind, delvy_wind, delvz_wind
	real delbx_wind, delby_wind, delbz_wind
	real sbx_wind, sby_wind, sbz_wind
	
	real deltg, deltinj, delt
	real b01, b02
	
	real	:: tgraph = 0.
	real	:: tinj = 0.
	real	:: tdiv = 10.
	real	:: del_tdiv = 5.
	
	!	spacecraft decides whether to include any spacecraft at all
	!	(input or output)
	!	(input_fluid is not currently implemented)
	
	!	ringo decides if you you want to plot 1 set of diagnostics
	!	with no time stepping.
	!	update decides if time setting is to be reset
	
	!	rho, pres, erg, px, py are the density, pressure, energy and
	!	momentum in the x and y directions, respectively, of the fluid.
	!	The indices 1, 2 are old, new values respectively.
	
	!	ijsrf gives position of ionosphere in grid units.
	!	Plasma parameters are stored in parm_srf.
	!	ijmid gives intermediate boundary, representing the	atmosphere.
	!	Plasma parameters are stored in parm_mid.
	!	ijzero gives position of all grid units interior to the inner
	!	boundary surface.
	
	!	d_min is the minimum allowable density
	!	stepsz is the size of the time step in terms of
	!	delta_x,y / (|umax|+c_s)
	
	!	Overall system dimensions are nx,ny,nz.
	!	Variable grid spacing enabled with rxyz > 1
	!	rx,ry,rz should not be set identically to zero

	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	@		SPACECRAFT INITIALIZATION		@
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	!	xcraft is the actual position of the spacecraft in body radii
	!	4th dimension is UT in hours
	!	zcraft is the future position of the spacecraft in body radii
	!	rcraft is the position of the spacecraft for which IMF is
	!	reference. No alteration from boundary conditions applied,
	!	and no data recorded. Not included in ncraft or ndef_craft.
	!	craft_info is directory for spacecraft position/times
	!	(relative to main directory)
	!	craft_datadir is directory for output data (relative to main
	!	directory)

	!	craftpos contains all (x,y,z,t) for all aux craft 
	!	sxyz is linearly interpolated spacecraft xyz for a given UT.
	!	gridpts holds two vertices of the cube containing the 8 grid
	!	coordinates closest to the spacecraft for a given UT.
	!	num_inst is how many 'instruments' the spacecraft has,
	!	i.e. how many scalar quantities are output to .dat file.
	!	When new measurements are to be added to spacecraft .dat files,
	!	num_inst must be incremented and 'SPACECRAFT DATA RECORDING'
	!	section also needs a new line for the new physical quantity:
	!	scdata(i) = physical_qty.
	!	This must be done in two places, one for default craft and one
	!	for aux craft.

	!	craftstat holds IOSTAT for reading in .craft files. Goes
	!	negative if EOF is reached before read() is done reading values.
	!	deflt_rec_skp is the number of time stepping loops to skip
	!	between default craft measurements.
	!	nskipped is the number of time stepping loops during which we
	!	skipped default craft recording.
	!	ntimes(:,1) holds the number of measurements an aux craft has
	!	made.
	!	ntimes(:,2) holds the number of (x,y,z,t) points we have for
	!	each aux craft.
	!	vals holds the max of ntimes(:,2) for use in subroutines.
	!	cgridpt(1:3) holds the xyz grid indices of the nearest grid
	!	point for a given craft. Default craft must be placed exactly
	!	on a grid point.
	!	cgridpt(4) holds the box number for the smallest box that fits
	!	our craft.
	!	recording is a boolean flag, true by default, which is set to
	!	false when a spacecraft has recorded data for its entire
	!	trajectory or when there is a problem.
	!	craft_rot_ca(:) holds values of rot_hrs for the orbiting moon
	!	for the time of closest approach for each trajectory craft.
	!	flyby_ref_time(:) is the sim time ut value to compare flyby
	!	relative times against.
	!	spam_reduct and spam_limit are used to reduce the amount of spam
	!	log messages printed for spacecraft data recording by a factor
	!	of spam_limit.

	!	Spacecraft data recording has a variable number of header lines.
	!	A dummy file is created to print a header and count the number
	!	of lines to skip when seeking to the end of these files. The
	!	number of header lines is stored in nheadlines.
	
	character*32, parameter	:: craft_info='spacecraft_info/'
	character*32, parameter	:: craft_datadir='data/spacecraft_data/'
	character*120	junkline
	character*8		numstring
	
	integer, parameter	:: deflt_rec_skp = 1000
	integer, parameter	:: num_inst = 3	! (Bx, By, Bz)
	character*8, parameter :: header_names(num_inst+4) = [ 'ut(h)', &
		'xpos', 'ypos', 'zpos', 'Bxval', 'Byval', 'Bzval' ]
	character*200	:: dat_header
	integer	:: nskipped = 0
	integer, parameter	:: spam_limit = 100
	integer	:: spam_reduct = spam_limit
	integer	ncraft
	integer	ndef_craft
	integer	naux_craft
	integer nheadlines
	integer n_vals
	
	!	For formatting spacecraft .dat file output:
	character*2 :: num_inst_char
	character*9 :: scfmt = '('
	character*6 :: scfmt_end = 'E15.7)'
	character*14 :: sc_headfmt = '('
	character*11 :: sc_headfmt_end = '(6X,A8,1X))'
	
	real	scdata(num_inst)
	real	rcraft(3)
	real	sxyz(3)
	real	sctime
	
	character*8, allocatable	:: craftnames(:)
	integer, allocatable		:: cgridpt(:,:)
	integer, allocatable		:: ntimes(:,:)
	logical, allocatable		:: recording(:)
	real, allocatable			:: xcraft(:,:)
	real, allocatable			:: zcraft(:,:)
	real, allocatable			:: craftpos(:,:,:)
	real, allocatable			:: craft_rot_ca(:)
	real, allocatable			:: flyby_ref_time(:)
	
	!	**********************
	!	Namelists for file I/O
	!	**********************

	namelist/option/tmax,ntgraph,stepsz,start,tsave,ntinj,isotropic, &
		run_name,write_dat,diagnostics
	namelist/planet/bodyname,moonname,xdip,ydip,zdip,r_inner,torus_rad, &
		tilt1,tilt2,tilting,rmassq,rmassh,rmasso
	namelist/speeds/cs_inner,alf_inner1,alf_inner2, &
		alpha_e,denh_inner,denq_torus,denh_torus,deno_torus, &
		gravity,ti_te,gamma,ringo,update,reload, &
		divb_lores,divb_hires,reduct,t_torus,aniso_factor,aniso_limit, &
		qvlim,hvlim,ovlim,qclim,hclim,oclim,cslim
	namelist/windy/re_wind,cs_wind,vx_wind1,vx_wind2, &
		vy_wind1,vy_wind2,vz_wind1,vz_wind2, &
		alfx_wind1,alfx_wind2, &
		alfy_wind1,alfy_wind2, &
		alfz_wind1,alfz_wind2, &
		den_wind1,den_wind2, &
		reynolds,resist,o_conc,rho_frac,bfrac,vfrac
	namelist/physical/re_equiv,v_equiv,rho_equiv, &
		spacecraft,input_fluid,warp,utstart,repeat_flybys
	namelist/smooth/chirho,chipxyz,chierg, &
		difrho,difpxyz,diferg
	namelist/crafthead/cname,num_vals,rot_closest,git_hash

	!	************************
	!	Planet & moon parameters
	!	************************
	
	!	Notes:	
	!		Calculations continued after input file is read in
	
	!		Adapted parameters for simulation use (boundaries etc.,
	!		see 'Planet & moon calculations')
	!		torus_inj_dist is used for making a correction to plasma
	!		torus injection
	real lunar_rad, torus_inj_dist, grav	
	real r_orbit, v_orbit, tilt
	real xmoon, ymoon, zmoon, rmoon, b0_moon
	
	!	*************
	!	Common blocks
	!	*************

	real, dimension(nx,ny,nz) :: vvx,vvy,vvz, tvx,tvy,tvz, evx,evy,evz
	real v_rot, r_rot, sin_tilt, cos_tilt, b0
	real planet_orbit_rad, planet_year, planet_rad, &
		planet_per, planet_mass, planet_obliq, planet_incl, &
		r_lim, torus_infall, planet_tilt, planet_init_long, &
		planet_xdip, planet_ydip, planet_zdip, torus_dist, &
		moon_orbit_rad, moon_per, moon_rad, moon_mass, &
		moon_incl, moon_init_rot
	
	common /space/vvx,vvy,vvz, tvx,tvy,tvz, evx,evy,evz
	common /rotation/v_rot,r_rot,rot_angle, xdip,ydip,zdip, &
		sin_tilt,cos_tilt,b0
	common /planetary/planet_orbit_rad, planet_year, planet_rad, &
		planet_per, planet_mass, planet_obliq, planet_incl, &
		r_lim, torus_infall, planet_tilt, planet_init_long, &
		planet_xdip, planet_ydip, planet_zdip, torus_dist, &
		moon_orbit_rad, moon_per, moon_rad, moon_mass, moon_incl, &
		moon_init_rot
	common /output_f/fluxs_f, speed_f, concs_f, grdpt_f, recdt_f
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!							EXECUTION SECTION
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	write(*,*) '-'
	call system('date +"START: %T-%Y-%m-%d"')
	write(*,*) '-'
	
	!	*************
	!	Grab git hash
	!	*************
		
	call system ('git rev-parse --short HEAD > '//trim(git_hash_file))
	open(git_f,file=trim(git_hash_file),status='unknown',form='formatted')
	read(git_f,*) git_hash
	write(*,*) 'Git hash: ', git_hash
	close(git_f)
	
	write(*,*) '-'
	
	!	********************************
	!	Open input and output data files
	!	********************************
	
	open(input_f,file='input',status='old',form='formatted')		
	open(fluxs_f,file='fluxes.dat',status='unknown',form='formatted')
	open(speed_f,file='speeds.dat',status='unknown',form='formatted')
	open(concs_f,file='conc.dat',status='unknown',form='formatted')
	!open(grdpt_f,file='grid.dat',status='unknown',form='formatted')
	!open(25,file='cur.dat',status='unknown',form='unformatted')
	!open(26,file='pot.dat',status='unknown',form='unformatted')
	if(write_dat) open(recdt_f,file='record.dat', &
		status='unknown',form='formatted')
	
	!	******************
	!	Open NCAR graphics
	!	******************
	
	!	ncar graphics disabled.
	!	call opngks
	!	call gsclip(0)

	!	********************
	!	Set NCAR color table
	!	********************

	!	call cpclrs

	!	call gselnt(0)
	!	call gsplci(1)
	!	call wtstr(.4,.975,'3d mutant code',2,0,0)
	
	!	*********************
	!	Read input parameters
	!	*********************
	
	read(input_f,option)
	read(input_f,planet)
	read(input_f,speeds)
	read(input_f,windy)
	read(input_f,physical)
	read(input_f,smooth)

	!	***********************************
	!	Write input data to stdout/log file
	!	***********************************
	
	write(*,option)
	write(*,planet)
	write(*,speeds)
	write(*,windy)
	write(*,physical)
	write(*,smooth)
	
	write(*,*) 'Grid xmin: ', 'Grid xmax: ', 'Grid ymin: ', &
		'Grid ymax: ', 'Grid zmin: ', 'Grid zmax: ', 'Grid spacing: '
	do i=1,n_grids
		xspac(i) = division**(i-1)
		grid_minvals(1,i) = -limit*xspac(i)
		grid_maxvals(1,i) = limit*xspac(i)
		grid_minvals(2,i) = -limit*xspac(i)
		grid_maxvals(2,i) = limit*xspac(i)
		grid_minvals(3,i) = -limit/2*xspac(i)
		grid_maxvals(3,i) = limit/2*xspac(i)
		
		if(i.eq.n_grids) then
			grid_diff = grid_maxvals(1,i) * wind_adjust - grid_maxvals(1,i)
			grid_maxvals(1,i) = grid_maxvals(1,i) + grid_diff
			grid_minvals(1,i) = grid_minvals(1,i) + grid_diff
		endif
		
		write(*,*) grid_minvals(1,i), grid_maxvals(1,i), &
			grid_minvals(2,i), grid_maxvals(2,i), &
			grid_minvals(3,i), grid_maxvals(3,i), xspac(i)
		ix = 1 + ( grid_maxvals(1,i) - grid_minvals(1,i) ) / xspac(i)
		iy = 1 + ( grid_maxvals(2,i) - grid_minvals(2,i) ) / xspac(i)
		iz = 1 + ( grid_maxvals(3,i) - grid_minvals(3,i) ) / xspac(i)

		if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz)) then
			write(*,*) 'ERROR: grid spacing does not match ', &
				'number of points.', ' Actual x,y,z: ', ix,iy,iz, &
				'Expected x,y,z: ', nx,ny,nz
			stop
		endif
	enddo
	
	!	**********************************************
	!	Check if subgrids mate properly to larger grid
	!	**********************************************
	
	do box=1,n_grids-1
		next_box = box + 1
	
		grid_diff = 1. + (grid_minvals(1,box) - grid_minvals(1,next_box)) &
			/ xspac(next_box)
		mating_index = int(grid_diff)
		dx = grid_diff - mating_index
		if(abs(dx).gt.0.001) then
			write(*,*)'Main grid: xmin dont match. ', &
				box, next_box, mating_index, grid_diff
			stop
		endif
	
		grid_diff = 1. + (grid_minvals(2,box) - grid_minvals(2,next_box)) &
			/ xspac(next_box)
		mating_index = int(grid_diff)
		dy = grid_diff - mating_index
		if(abs(dy).gt.0.001) then
			write(*,*)'Main grid: ymin dont match. ', &
				box, next_box, mating_index, grid_diff
			stop
		endif
	
		grid_diff = 1. + (grid_minvals(3,box) - grid_minvals(3,next_box)) &
			/ xspac(next_box)
		mating_index = int(grid_diff)
		dz = grid_diff - mating_index
		if(abs(dz).gt.0.001) then
			write(*,*)'Main grid: zmin dont match. ', &
				box, next_box, mating_index, grid_diff
			stop
		endif
	enddo
	
	!	**************************
	!	Set inner boundary regions
	!	**************************
	
	zero_bndry = r_inner + zero_bndry_offset
	mid_bndry = r_inner + mid_bndry_offset
	srf_bndry = r_inner + srf_bndry_offset
	
	!	**************************
	!	Planet & moon calculations
	!	**************************
	
	!	Uses astrometry module to set planetary parameters. Import list:
	!	planet_orbit_rad,	planet_year,
	!	planet_rad,			planet_mass,
	!	planet_obliq,		planet_incl,
	!	r_lim,				torus_infall,
	!	planet_tilt,		planet_init_long,
	!	planet_xdip,		planet_ydip,
	!	planet_zdip,		torus_dist,
	!	moon_orbit_rad,		moon_per,
	!	moon_rad,			moon_mass,
	!	moon_incl,			moon_init_rot
	call choose_system(bodyname,moonname)
	!	No rotation relative to plasma flow for tidally locked moons
	if(bodyname.eq."europa") planet_per = 0
	
	r_rot = r_lim
	!	Normalized rotation rate
	v_rot = 2.*pi*planet_rad/(planet_per*3600.)/v_equiv 
	torus_inj_dist = torus_dist + torus_infall
	!	Start exobase at 1.25 rt
	lunar_rad = 1.25 * moon_rad
	
	!	rmoon in # grid points
	rmoon = ( lunar_rad / planet_rad ) / re_equiv
	!	v_orbit is orbit velocity in sim units
	v_orbit = (moon_orbit_rad*2.*pi)/(moon_per*3600.)/v_equiv
	
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
	
	!	Gravity in m/s**2 at planet's surface 
	!	Need t_equiv in normalization
	t_equiv = planet_rad * re_equiv / v_equiv
	grav = gravity * (planet_rad*re_equiv*1000.) / (1000.*v_equiv)**2
	
	!	******************************
	!	Calculate inner boundary sizes
	!	******************************
	
	!	Calculates approx. number of grid point volume units in each
	!	region. 1.2 is fudge factor to make sure we include enough
	!	working points.
	msrf  = 1.2 * 1.33*pi*(srf_bndry**3 -mid_bndry**3)/xspac(mbndry)**3
	mmid  = 1.2 * 1.33*pi*(mid_bndry**3 -zero_bndry**3)/xspac(mbndry)**3
	mzero = 1.2 * 1.33*pi*(zero_bndry**3)/xspac(mbndry)**3
	
	!	*********************
	!	Untracked grid points
	!	*********************
	
	allocate( ijsrf(mbndry,3,msrf), ijmid(mbndry,3,mmid), &
		ijzero(mbndry,3,mzero), parm_srf(mbndry,num_zqt,msrf), &
		parm_mid(mbndry,num_zqt,mmid), parm_zero(mbndry,num_zqt,mzero) )
	
	
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	@		PARAMETER DECLARATIONS COMPLETE			@
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	!	Calculate effective magnetic field strength
	
	erho = denh_inner * rmassh
	b01 = alf_inner1 * sqrt(erho) * r_inner**3
	b02 = alf_inner2 * sqrt(erho) * r_inner**3
	b0 = b01
	write(*,*)'b0: ',b0
	
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!		@		RECONFIGURE INPUT DATA		@
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	!	The magnetic field, in dimensionless units, is in Alfven speed
	!	relative to the normalizing velocity
	!	The temperature is in sound speeds, relative to the normalizing
	!	speed, all squared.
	
	!	The magnetospheric plasma densities are assumed to vary as
	!	rho proportional to (r)**-alpha_e
	!	temperature proportional to (r)**alpha_e
	!	with total pressure constant, which is needed for equilibrium
	
	!	*****************
	!	Find equiv values
	!	*****************
	
	!	Sim units are in units of equiv values. These quantities are
	!	typically in SI times a power of 10, and are important for
	!	keeping simulation calculations as close to unity as possible.
	!	Numerical factors appear for the following conversions:
	!	rho_equiv is in cm^-3 of rmassq species
	!	v_equiv is in km/s
	
	! b_equiv is in nT
	b_equiv = v_equiv*1.e3 * sqrt(mu0 * rho_equiv*1.e6*m_prot) *1.e9

	!	Now, find the equivalent pressure of magnetosphere for the given
	!	sound speed
	
	gamma1 = gamma - 1.
	epress = cs_inner**2 * erho / gamma
	eerg = epress / gamma1
	
	!	Do the same for the solar wind
	
	rho_wind1 = den_wind1 * rmassq
	rho_wind2 = den_wind2 * rmassq
	srho = rho_wind1
	delrho = (rho_wind2 - rho_wind1) / tmax
	spress = (cs_wind**2 * srho / gamma) / gamma1
	svelx = vx_wind1
	svely = vy_wind1
	svelz = vz_wind1
	
	spx = srho * svelx
	spy = srho * svely
	spz = srho * svelz
	serg = 0.5 * (svelx**2 + svely**2 + svelz**2) * srho + spress
	delvx_wind = (vx_wind2 - vx_wind1) / tmax
	delvy_wind = (vy_wind2 - vy_wind1) / tmax
	delvz_wind = (vz_wind2 - vz_wind1) / tmax
	
	delbx_wind = ( alfx_wind2 * sqrt(rho_wind2) - alfx_wind1 &
		* sqrt(rho_wind1) ) / tmax
	delby_wind = ( alfy_wind2 * sqrt(rho_wind2) - alfy_wind1 &
		* sqrt(rho_wind1) ) / tmax
	delbz_wind = ( alfz_wind2 * sqrt(rho_wind2) - alfz_wind1 &
		* sqrt(rho_wind1) ) / tmax
	sbx_wind = alfx_wind1 * sqrt(rho_wind1)
	sby_wind = alfy_wind1 * sqrt(rho_wind1)
	sbz_wind = alfz_wind1 * sqrt(rho_wind1)
	
	deltg = tmax / float(ntgraph)
	deltinj = tmax / float(ntinj)
	delt = stepsz
	
	!	*****************
	!	Check for restart
	!	*****************
	
	if(.not.start) then
	
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!	@			RESTART FROM FLUID FILE				@
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
		if(reload) then
			open(fluid_f,file=fname_reload, &
				status='unknown',form='unformatted')
		else
			!	Restart fluid file typically named fluid_start
			open(fluid_f,file=fname_restart, &
				status='unknown',form='unformatted')
		endif
		
		!	*****************
		!	Read restart data
		!	*****************

		read(fluid_f) t
		read(fluid_f) qrho
		read(fluid_f) qpx
		read(fluid_f) qpy
		read(fluid_f) qpz
		read(fluid_f) qpresx
		read(fluid_f) qpresy
		read(fluid_f) qpresz
		read(fluid_f) qpresxy
		read(fluid_f) qpresxz
		read(fluid_f) qpresyz
		read(fluid_f) hrho
		read(fluid_f) hpx
		read(fluid_f) hpy
		read(fluid_f) hpz
		read(fluid_f) hpresx
		read(fluid_f) hpresy
		read(fluid_f) hpresz
		read(fluid_f) hpresxy
		read(fluid_f) hpresxz
		read(fluid_f) hpresyz
		read(fluid_f) orho
		read(fluid_f) opx
		read(fluid_f) opy
		read(fluid_f) opz
		read(fluid_f) opresx
		read(fluid_f) opresy
		read(fluid_f) opresz
		read(fluid_f) opresxy
		read(fluid_f) opresxz
		read(fluid_f) opresyz
		read(fluid_f) bx
		read(fluid_f) by
		read(fluid_f) bz
		read(fluid_f) epres
		read(fluid_f) bx0
		read(fluid_f) by0
		read(fluid_f) bz0
		read(fluid_f) parm_srf,parm_mid,parm_zero,ijzero,numzero, &
			ijmid,nummid,ijsrf,numsrf
		close(fluid_f)
		
		write(*,*) 'Restart from file: ', fname_restart
		fluid_f = fluid_f + 1
		
		!###########################################
		ut=utstart	!	utstart read from input file
		!###########################################
		if(utstart .gt. 0.00001) then
			t_diff = abs(ut*3600./t_equiv - t)
			if(t_diff .gt. 1.0) then
				if(t_diff .lt. 100.0) then
					write(*,*) '+++++++'
					write(*,*) 'WARNING: Non-zero utstart differs ', &
						'from t in restart fluid file.'
					write(*,*) 'utstart = ', utstart, ' t = ', t, &
						' t_equiv = ', t_equiv, &
						' ut(s)/t_equiv - t = ', t_diff
					write(*,*) '+++++++'
				else
					write(*,*) '+++++++'
					write(*,*) 'ERROR: utstart drastically differs ', &
						'from t in restart fluid file. ABORTING!'
					write(*,*) 'utstart = ', utstart, ' t = ', t, &
						' t_equiv = ', t_equiv, &
						' ut(s)/t_equiv - t = ', t_diff
					write(*,*) '+++++++'
					stop
				endif
			endif
			t = ut*3600./t_equiv
		endif
		
		!	***************
		!	Rotation timing
		!	***************

		nrot_planet = int( ut/planet_per )
		rot_hrs = ut - planet_per*nrot_planet
		rot_angle = 2.*pi*rot_hrs/planet_per
		sin_rot = sin(rot_angle)
		cos_rot = cos(rot_angle)
		
		nrot_moon = int( ut/moon_per )
		rot_hrs_moon = ut - moon_per*nrot_moon
		rot_angle_moon = 2.*pi*rot_hrs_moon/moon_per + moon_init_rot
		sin_rot_moon = sin(rot_angle_moon)
		cos_rot_moon = cos(rot_angle_moon)
		
		xmoon = moon_orbit_rad * cos_rot_moon
		ymoon = moon_orbit_rad * sin_rot_moon
		zmoon = 0.
		
		!	Check for div b errors
		
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
			
			!	Apply boundary conditions
			
			do box=n_grids-1,1,-1
				call flanks_synced(bx,nx,ny,nz,n_grids,box, &
					grid_minvals(1,:), grid_maxvals(1,:), &
					grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
				call flanks_synced(by,nx,ny,nz,n_grids,box, &
					grid_minvals(1,:), grid_maxvals(1,:), &
					grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
				call flanks_synced(bz,nx,ny,nz,n_grids,box, &
					grid_minvals(1,:), grid_maxvals(1,:), &
					grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
			enddo

			if(isotropic) then
				qpresy = qpresx
				qpresz = qpresx
				qpresxy = 0.
				qpresxz = 0.
				qpresyz = 0.

				opresy = opresx
				opresz = opresx
				opresxy = 0.
				opresxz = 0.
				opresyz = 0.

				hpresy = hpresx
				hpresz = hpresx
				hpresxy = 0.
				hpresxz = 0.
				hpresyz = 0.
			endif

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
					grid_minvals(1,:), grid_maxvals(1,:), &
					grid_minvals(2,:), grid_maxvals(2,:), &
					grid_minvals(3,:), grid_maxvals(3,:) )
			enddo  !end bndry_corer
		endif  ! end divb_lores
		
		if(spacecraft) write(*,*) 'Spacecraft positions not yet initialized.'

		!	NCAR graphics disabled.
		!	write(*,*) 'Entering lores visual'
		!	call visual(qrho,qpresx,qpresy,qpresz,qpresxy, &
		!		qpresxz,qpresyz,qpx,qpy,qpz,rmassq, &
		!		hrho,hpresx,hpresy,hpresz,hpresxy, &
		!		hpresxz,hpresyz,hpx,hpy,hpz,rmassh, &
		!		orho,opresx,opresy,opresz,opresxy, &
		!		opresxz,opresyz,opx,opy,opz,rmasso, &
		!		epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
		!		curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
		!		tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
		!		nx,ny,nz,n_grids,xspac, &
		!		cross,along,flat,xcraft,ncraft,re_equiv, &
		!		grid_minvals(1,:), grid_maxvals(1,:), &
		!		grid_minvals(2,:), grid_maxvals(2,:), &
		!		grid_minvals(3,:), grid_maxvals(3,:), &
		!		ut,b_equiv,ti_te,rho_equiv)

		ts1 = t + tsave
		tstep = tmax
		tmax = t + tmax
		tgraph = t + deltg
		tinj = t + deltinj
		tdiv = t
		
		!	Initialize plasma resistivity
		
		call set_resist(resistive,nx,ny,nz,mbndry,resist, &
			ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
			msrf,mmid,mzero,1.)
		
		write(*,*) 'Graphing data.'
		call write_graphing_data( &
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
			r_inner, run_name, dummy_fg, diagnostics, nplots, &
			isotropic, smallbit)
	else
		
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!	@			INITIALIZATION			@
		!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
		!###########################################
		ut=utstart	!	utstart read from input file
		!###########################################
		if(utstart .gt. 0.00001) t = ut*3600./t_equiv
		
		ts1 = t + tsave
		tmax = t + tmax
		tgraph = t + deltg
		tinj = t + deltinj
		tdiv = t
		fluid_f = fluid_f + 1
		nplots = 0
		
		!	***************
		!	Rotation timing
		!	***************

		nrot_planet = int( ut/planet_per )
		rot_hrs = ut - planet_per*nrot_planet
		rot_angle = 2.*pi*rot_hrs/planet_per
		sin_rot = sin(rot_angle)
		cos_rot = cos(rot_angle)
		
		nrot_moon = int( ut/moon_per )
		rot_hrs_moon = ut - moon_per*nrot_moon
		rot_angle_moon = 2.*pi*rot_hrs_moon/moon_per + moon_init_rot
		sin_rot_moon = sin(rot_angle_moon)
		cos_rot_moon = cos(rot_angle_moon)
		
		xmoon = moon_orbit_rad * cos_rot_moon
		ymoon = moon_orbit_rad * sin_rot_moon
		zmoon = 0.
		
		!	**************
		!	Ambient plasma
		!	**************

		!	Initialize indices and surface points
		
		numsrf=0
		nummid=0
		numzero=0
		
		write(*,*) 'Initial torus_inj_dist, rot_angle: ', &
			torus_inj_dist, rot_angle
		
		do box=1,n_grids
			dx = ( grid_maxvals(1,box) - grid_minvals(1,box) ) &
				/ ( nx - 1. )
			dy = ( grid_maxvals(2,box) - grid_minvals(2,box) ) &
				/ ( ny - 1. )
			dz = ( grid_maxvals(3,box) - grid_minvals(3,box) ) &
				/ ( nz - 1. )
		
			!	Create dipole magnetic field and load magnetospheric
			!	plasma
		
			do k=1,nz
				az = grid_minvals(3,box) + dz * (k-1) - zdip
		
				do j=1,ny
					ay = grid_minvals(2,box) + dy * (j-1) - ydip
		
					do i=1,nx
						ax = grid_minvals(1,box) + dx * (i-1) - xdip
		
						!	Rotate coordinates for planet motion
						
						xr=ax*cos_rot+ay*sin_rot
						yr=-ax*sin_rot+ay*cos_rot
						
						!	Tilt space to dipole space
						
						xp=xr*cos_tilt-az*sin_tilt
						yp=yr
						zp=xr*sin_tilt+az*cos_tilt
						ar=sqrt(xp**2+yp**2+zp**2)
						
						!	Determine magnetic dipole field
						
						call dipole(bx0(i,j,k,box),by0(i,j,k,box), &
							bz0(i,j,k,box),ax,ay,az)
						
						!	Zero interior magnetic field,
						!	so Alfven speed is small there
						
						if (ar.lt.zero_bndry) then
							bx0(i,j,k,box)=0.
							by0(i,j,k,box)=0.
							bz0(i,j,k,box)=0.
						endif
						
						!	Set up rotational properties
						
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
						
						!	For Jupiter top-hat distribution
						
						ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
						
						ar_iono=amax1(0.0001,ar_iono)
						ra=( (ar_iono+0.5*r_inner) / (1.5*r_inner) &
							)**(-alpha_e)
						zheight=amax1( 1., (zp**2+(1.5*r_inner)**2) &
							/(3.0*r_inner)**2 )
						ra=ra/zheight**2
						rho_iono=amax1(erho*ra,d_min)
						
						r_equat=(ar**3+0.001)/(xp**2+ay**2+0.001)
						r_equat=amax1(r_equat,r_inner)
						rerg_sphere=(0.001+(r_inner/r_equat)**4)
						
						if(r_equat.gt.5.*r_inner) then
							!	Constant pressure flux
							rerg_sphere =  0.001 &
								+ (r_inner/r_equat)**alpha_e
							rden_sphere=1.
						else
							rerg_sphere = 0.1 &
								+ 0.9*(r_equat-r_inner)/(4.*r_inner)
							!	Reduce O+ in plasmasphere
							rden_sphere = rerg_sphere**2	
						endif
						
						arho = amax1( rho_iono, d_min )
						vt = amax1( sqrt(rvy**2 + rvx**2), cs_inner )
						
						!	MAT - Corrected according to
						!	doi:10.1029/JA094iA12p17287
						!	and 
						!	doi:10.1029/2008GL035433
						
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
						
						!	Add small amount of heavies everywhere
						
						hrho(i,j,k,box)=0.05*qrho(i,j,k,box) &
							* rmassh/rmassq
						hpresx(i,j,k,box)=0.05*qpresx(i,j,k,box)
						hpresy(i,j,k,box)=hpresx(i,j,k,box)
						hpresz(i,j,k,box)=hpresx(i,j,k,box)
						hpresxy(i,j,k,box)=0.
						hpresxz(i,j,k,box)=0.
						hpresyz(i,j,k,box)=0.
						hpx(i,j,k,box)=hrho(i,j,k,box)*rvx
						hpy(i,j,k,box)=hrho(i,j,k,box)*rvy
						hpz(i,j,k,box)=0.
						
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
						
						epres(i,j,k,box)=(qpresx(i,j,k,box) &
							+ hpresx(i,j,k,box) &
							+ opresx(i,j,k,box) )/ti_te
						
						bx(i,j,k,box)=0.
						by(i,j,k,box)=0.
						bz(i,j,k,box)=0.
						
						!	Test for boundary point of planets surface
						!	or interior to planet
						
						if((ar.le.srf_bndry).and.(box.le.mbndry)) then
						
							if(ar.lt.zero_bndry) then
								numzero(box)=numzero(box)+1
								ijzero(box,1,numzero(box))=i
								ijzero(box,2,numzero(box))=j
								ijzero(box,3,numzero(box))=k
						
								parm_zero(box,1,numzero(box)) &
									= qrho(i,j,k,box)
								parm_zero(box,2,numzero(box)) &
									= hrho(i,j,k,box)
								parm_zero(box,3,numzero(box)) &
									= orho(i,j,k,box)
								parm_zero(box,4,numzero(box)) &
									= qpresx(i,j,k,box)
								parm_zero(box,5,numzero(box)) &
									= hpresx(i,j,k,box)
								parm_zero(box,6,numzero(box)) &
									= opresx(i,j,k,box)
								parm_zero(box,7,numzero(box)) &
									= epres(i,j,k,box)
							else if(ar.lt.mid_bndry) then
								nummid(box)=nummid(box)+1
								ijmid(box,1,nummid(box))=i
								ijmid(box,2,nummid(box))=j
								ijmid(box,3,nummid(box))=k
								ar2=sqrt(xp**2+ay**2)
								parm_mid(box,1,nummid(box)) &
									= qrho(i,j,k,box)
								parm_mid(box,2,nummid(box)) &
									= hrho(i,j,k,box)
								parm_mid(box,3,nummid(box)) &
									= orho(i,j,k,box)
								parm_mid(box,4,nummid(box)) &
									= qpresx(i,j,k,box)
								parm_mid(box,5,nummid(box)) &
									= hpresx(i,j,k,box)
								parm_mid(box,6,nummid(box)) &
									= opresx(i,j,k,box)
								parm_mid(box,7,nummid(box)) &
									= epres(i,j,k,box)
							else
								numsrf(box)=numsrf(box)+1
								ijsrf(box,1,numsrf(box))=i
								ijsrf(box,2,numsrf(box))=j
								ijsrf(box,3,numsrf(box))=k
								ar2=sqrt(xp**2+ay**2)
								parm_srf(box,1,numsrf(box)) &
									= qrho(i,j,k,box)
								parm_srf(box,2,numsrf(box)) &
									= hrho(i,j,k,box)
								parm_srf(box,3,numsrf(box)) &
									= orho(i,j,k,box)
								parm_srf(box,4,numsrf(box)) &
									= qpresx(i,j,k,box)
								parm_srf(box,5,numsrf(box)) &
									= hpresx(i,j,k,box)
								parm_srf(box,6,numsrf(box)) &
									= opresx(i,j,k,box)
								parm_srf(box,7,numsrf(box)) &
									= epres(i,j,k,box)
							endif
						endif  ! end bndry condition test
					enddo	! end x loop
				enddo	! end y loop
			enddo	! end z loop
		enddo	! end boxes loop
		
		write(*,*) 'Interior and zero points:'
		write(*,'(4(A12))') 'box','srf','mid','zero'
		do box=1,mbndry
			write(*,'(4(I12))') box,numsrf(box),nummid(box),numzero(box)
		enddo
		
		!	Initialize solar wind plasma. Inserted beyond
		!	the planet at a radius of re_wind.
		
		wind_bnd = r_rot / re_equiv
		ofrac = rho_frac
		do box=1,n_grids
		
			dx = ( grid_maxvals(1,box) - grid_minvals(1,box) ) / (nx-1.)
			dy = ( grid_maxvals(2,box) - grid_minvals(2,box) ) / (ny-1.)
			dz = ( grid_maxvals(3,box) - grid_minvals(3,box) ) / (nz-1.)
			
			do k=1,nz
				az = grid_minvals(3,box) + dz * (k-1) - zdip

				do j=1,ny
					ay = grid_minvals(2,box) + dy * (j-1) - ydip

					do i=1,nx
						ax = grid_minvals(1,box) + dx * (i-1) - xdip
						
						ar=sqrt(ax**2+ay**2+az**2)
						if( (ar.ge.1.5*wind_bnd) .and. (ax.lt.0.) ) then
							qrho(i,j,k,box)=srho+qrho(i,j,k,box)
							qpx(i,j,k,box)=spx+qpx(i,j,k,box)
							qpy(i,j,k,box)=spy+qpy(i,j,k,box)
							qpz(i,j,k,box)=spz+qpz(i,j,k,box)
							
							avx=qpx(i,j,k,box)/qrho(i,j,k,box)
							avy=qpy(i,j,k,box)/qrho(i,j,k,box)
							avz=qpz(i,j,k,box)/qrho(i,j,k,box)
							
							apres=cs_wind**2*qrho(i,j,k,box)/gamma
							qpresx(i,j,k,box) = 0.5*apres &
								+ qpresx(i,j,k,box)
							qpresy(i,j,k,box)=qpresx(i,j,k,box)
							qpresz(i,j,k,box)=qpresx(i,j,k,box)
							qpresxy(i,j,k,box)=0.
							qpresxz(i,j,k,box)=0.
							qpresyz(i,j,k,box)=0.
							
							epres(i,j,k,box) = 0.5*apres/ti_te &
								+ epres(i,j,k,box)
							
							hrho(i,j,k,box) = srho*rho_frac &
								* rmassh/rmassq + hrho(i,j,k,box)
							hpx(i,j,k,box)=avx*hrho(i,j,k,box)
							hpy(i,j,k,box)=avy*hrho(i,j,k,box)
							hpz(i,j,k,box)=avz*hrho(i,j,k,box)
							hpresx(i,j,k,box) = 0.5*spress*rho_frac &
								+ hpresx(i,j,k,box)
							hpresy(i,j,k,box)=hpresx(i,j,k,box)
							hpresz(i,j,k,box)=hpresx(i,j,k,box)
							hpresxy(i,j,k,box)=0.
							hpresxz(i,j,k,box)=0.
							hpresyz(i,j,k,box)=0.
							
							orho(i,j,k,box) = srho*ofrac*rmasso/rmassq &
								+ orho(i,j,k,box)
							opx(i,j,k,box)=avx*orho(i,j,k,box)
							opy(i,j,k,box)=avy*orho(i,j,k,box)
							opz(i,j,k,box)=avz*orho(i,j,k,box)
							opresx(i,j,k,box) = 0.5*spress*ofrac &
								+ opresx(i,j,k,box)
							opresy(i,j,k,box)=opresx(i,j,k,box)
							opresz(i,j,k,box)=opresx(i,j,k,box)
							opresxy(i,j,k,box)=0.
							opresxz(i,j,k,box)=0.
							opresyz(i,j,k,box)=0.
						endif
						
						if((ar.gt.wind_bnd) .and. (ar.lt.1.5*wind_bnd) &
							.and.(ax.lt.0.0)) then

							dfrac=(ar-wind_bnd)/(0.5*wind_bnd)
							qrho(i,j,k,box) = srho*dfrac+qrho(i,j,k,box)
							qpx(i,j,k,box)=spx*dfrac+qpx(i,j,k,box)
							qpy(i,j,k,box)=spy*dfrac+qpy(i,j,k,box)
							qpz(i,j,k,box)=spz*dfrac+qpz(i,j,k,box)
							avx=qpx(i,j,k,box)/qrho(i,j,k,box)
							avy=qpy(i,j,k,box)/qrho(i,j,k,box)
							avz=qpz(i,j,k,box)/qrho(i,j,k,box)
							
							apres=cs_wind**2*qrho(i,j,k,box)/gamma
							qpresx(i,j,k,box)=0.5*apres*dfrac &
								+ qpresx(i,j,k,box)
							qpresy(i,j,k,box)=qpresx(i,j,k,box)
							qpresz(i,j,k,box)=qpresx(i,j,k,box)
							qpresxy(i,j,k,box)=0.
							qpresxz(i,j,k,box)=0.
							qpresyz(i,j,k,box)=0.
							
							epres(i,j,k,box)=0.5*apres/ti_te*dfrac &
								+ epres(i,j,k,box)
							
							hrho(i,j,k,box)=srho*rho_frac &
								* rmassh/rmassq*dfrac + hrho(i,j,k,box)
							hpx(i,j,k,box)=avx*hrho(i,j,k,box)
							hpy(i,j,k,box)=avy*hrho(i,j,k,box)
							hpz(i,j,k,box)=avz*hrho(i,j,k,box)
							hpresx(i,j,k,box) = 0.5*spress*rho_frac &
								* dfrac + hpresx(i,j,k,box)
							hpresy(i,j,k,box)=hpresx(i,j,k,box)
							hpresz(i,j,k,box)=hpresx(i,j,k,box)
							hpresxy(i,j,k,box)=0.
							hpresxz(i,j,k,box)=0.
							hpresyz(i,j,k,box)=0.
							
							orho(i,j,k,box) = srho*ofrac*rmasso/rmassq &
								* dfrac + orho(i,j,k,box)
							opx(i,j,k,box)=avx*orho(i,j,k,box)
							opy(i,j,k,box)=avy*orho(i,j,k,box)
							opz(i,j,k,box)=avz*orho(i,j,k,box)
							opresx(i,j,k,box) = 0.5*spress*ofrac &
								* dfrac + opresx(i,j,k,box)
							opresy(i,j,k,box)=opresx(i,j,k,box)
							opresz(i,j,k,box)=opresx(i,j,k,box)
							opresxy(i,j,k,box)=0.
							opresxz(i,j,k,box)=0.
							opresyz(i,j,k,box)=0.
						endif
					enddo
				enddo
			enddo
		enddo
		
		!	Check speeds to see if they are ridiculous
		
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
						
						avx=hpx(i,j,k,box)/hrho(i,j,k,box)
						avy=hpy(i,j,k,box)/hrho(i,j,k,box)
						avz=hpz(i,j,k,box)/hrho(i,j,k,box)
						spd=sqrt(avx**2+avy**2+avz**2)
						if(spd.gt.1.) then
							hpx(i,j,k,box)=0.
							hpy(i,j,k,box)=0.
							hpz(i,j,k,box)=0.
						endif
						
						avx=opx(i,j,k,box)/orho(i,j,k,box)
						avy=opy(i,j,k,box)/orho(i,j,k,box)
						avz=opz(i,j,k,box)/orho(i,j,k,box)
						spd=sqrt(avx**2+avy**2+avz**2)
						if(spd.gt.1.) then
							opx(i,j,k,box)=0.
							opy(i,j,k,box)=0.
							opz(i,j,k,box)=0.
						endif
					enddo
				enddo
			enddo
		enddo
		
		!	Initialize plasma resistivity
		
		call set_resist(resistive,nx,ny,nz,mbndry,resist, &
			ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
			msrf,mmid,mzero,1.)
		
	endif ! if(.not.start) then
	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!						FINISHED INITIAL SETUP!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) 'Did we start from scratch?', start
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) 'UT = ', ut
	
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!		@		SPACECRAFT DATA INITIALIZATION		@
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	if(spacecraft) then
		write(*,*) 'Preparing spacecraft...'

		write(num_inst_char,'(I2.2)') 4+num_inst
		scfmt = trim(scfmt)//num_inst_char//scfmt_end
		sc_headfmt = trim(sc_headfmt)//num_inst_char//sc_headfmt_end

		write(dat_header,sc_headfmt) adjustr(header_names(:))
		dat_header = trim(dat_header)
		write(*,*) 'WARNING: For plotting .dat files in Matlab, ', &
			'line break + extra space in header row must be ', &
			'deleted manually.'

		!	Count crafts named in craft_info directory.
		call countcraft(craft_info,scin,ncraft,ndef_craft,naux_craft)
		
		allocate( xcraft(4,ncraft), zcraft(4,ncraft) )
		allocate( craftnames(ncraft), ntimes(ncraft,2), &
			recording(ncraft) )
		allocate( cgridpt(4,ncraft), craft_rot_ca(ncraft), &
			flyby_ref_time(ncraft) )
		
		call count_hlines( dat_header, dummy_f, nheadlines)
		call names_ntimes( craft_info, scin, ncraft, ndef_craft, &
			naux_craft, nheadlines, craftnames, ntimes )
		write(*,'(5(A,1X))') 'Default spacecraft found: ', &
			craftnames(1:ndef_craft)
		write(*,'(5(A,1X))') 'Trajectory spacecraft found: ', &
			craftnames(ndef_craft+1:ncraft)
		
		n_vals = maxval(ntimes(:,2))
		!	Allocate array which stores all spacecraft position/time
		!	input values
		allocate(craftpos(4, ncraft, n_vals))
		
		!	Fill in spacecraft info arrays and initialize/seek output
		!	files
		call initcraft( craft_info, craft_datadir, craftnames, scin, &
			scout, dat_header, ncraft, ndef_craft, naux_craft, &
			nheadlines, recording, ntimes, n_vals, cgridpt, craftpos, &
			craft_rot_ca, n_grids, re_equiv, xspac, grid_minvals, &
			grid_maxvals, git_hash )
		
		!	------------------------
		!	FINISHED SPACECRAFT INFO
		!	------------------------
		
		!	Reference craft is a copy of 'wind' default craft
		rcraft(1) = xcraft(1,1)
		rcraft(2) = xcraft(2,1)
		rcraft(3) = xcraft(3,1)
		
		!	Set current spacecraft location to init or last read in
		!	from .dat
		do n = 1, ndef_craft
			!	Set spacecraft actual position to first (and only)
			!	value in craftpos
			xcraft(:,n) = craftpos(:,n,1)
			!	Default spacecraft typically do not move; set future
			!	locations to current location.
			zcraft(:,n) = xcraft(:,n)
		enddo

		write(*,*) 'Reading aux craft:'

		do n = ndef_craft+1, ncraft
			flyby_ref_time(n) = craft_rot_ca(n) + moon_per
			!	xcraft is the location of the next measurement
			xcraft(:,n) = craftpos( :, n, ntimes(n,1)+1 )
			!	zcraft time used as a placeholder for the moment
			zcraft(4,n) = xcraft(4,n)
			xcraft(4,n) = zcraft(4,n) + flyby_ref_time(n)

			do while(xcraft(4,n) .lt. ut)
				flyby_ref_time(n) = flyby_ref_time(n) + moon_per
				xcraft(4,n) = zcraft(4,n) + flyby_ref_time(n)
			enddo
			write(*,*) craftnames(n), ' flyby_ref_time = ', &
				flyby_ref_time(n)

			if(ntimes(n,1) .lt. ntimes(n,2)-1) then
				!	Now zcraft holds the measurement after next...
				zcraft(:,n) = craftpos(:, n, ntimes(n,1)+2 )
				zcraft(4,n) = zcraft(4,n) + flyby_ref_time(n)
			else
				!	Except when the next measurement is the last one.
				zcraft(:,n) = xcraft(:,n)
			endif
		enddo
	else
			!	warp is deprecated, though it may make a comeback.
			warp = .false.	
	endif	! endif(spacecraft)
	

	if(input_fluid) then
		!	This function is included for historical reasons; including
		!	upstream input functionality is an intended future feature,
		!	but is not currently implemented.
		write(*,*) 'input_fluid is set to .true., but is ', &
			'not currently implemented. Aborting.'
		stop
		!	Calculate the distance between the reference spacecraft and
		!	wind boundary
		distance=grid_minvals(1,n_grids)-rcraft(1)/re_equiv
		write(*,*)'Wind displaced by: ', distance	
		
		do m=1,ncts
			read(29,*) bfld(m,4),bmag,bfld(m,1),bfld(m,2),bfld(m,3)
			read(27,*) rut,rplas(m),svel(m,1),svel(m,2),svel(m,3)
			read(27,*) rut,rplas(m)
			read(28,*) vut,svel(m,1),svel(m,2),svel(m,3)
			!	Keep bx in solar wind constant
			bfld(m,1)=-sbx_wind*b_equiv
		enddo
		
		!	Set timing
		do  k=1,nz
			do  j=1,ny
				ncount(j,k)=0
				future(j,k)=ut-0.01
			enddo
		enddo
		
		nvx=0
		vut=-999.
		do while(ut.gt.vut)
			svelx=zvelx
			nvx=nvx+1
			zvelx=-svel(nvx,1)/v_equiv
			vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
		enddo
		
		write(*,*)'UT=',ut,' Wind time: ',bfld(nvx,4)
		
		displace=0.
		dx = ( grid_maxvals(1,n_grids) - grid_minvals(1,n_grids) ) &
			/ (nx-1.)
		dy = ( grid_maxvals(2,n_grids) - grid_minvals(2,n_grids) ) &
			/ (ny-1.)
		dz = ( grid_maxvals(3,n_grids) - grid_minvals(3,n_grids) ) &
			/ (nz-1.)
		
		do k=1,nz
			do j=1,ny
				do while( (ut.gt.future(j,k)) &
					.and.(ncount(j,k)+1.le.ncts) )
					nc=ncount(j,k)+1
					bxp(j,k)=bxf(j,k)
					byp(j,k)=byf(j,k)
					bzp(j,k)=bzf(j,k)
					rhop(j,k)=rhof(j,k)
					svxp(j,k)=svxf(j,k)
					svyp(j,k)=svyf(j,k)
					svzp(j,k)=svzf(j,k)
					past(j,k)=future(j,k)
					
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
					
					if(warp) then
						write(*,*) 'warp is set to .true., but is ', &
							'not currently implemented. Aborting.'
						stop
						b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
						b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
						!	WARNING: This is not in a loop over 'box'.
						ay = grid_minvals(2,box) + dy * (j-1) &
							- rcraft(2) / re_equiv
						az = grid_minvals(3,box) + dz * (k-1) &
							- rcraft(3) / re_equiv
					
						!	Assuming Bz IMF is positive on average 
						!	and that we can ignore transients
						!	Also assuming By IMF is negative
						displace= -bxf(j,k) &
							* (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
					endif
					ar=distance+displace
					future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
				enddo
			enddo
		enddo
		
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
	endif	!	end input_fluid if
	
	!	Check initial conditions
	
	!	write(*,*) 'checking set speed'
	do box=n_grids,1,-1
	
		!	Sync time steps
	
		t_old(box)=0.
		t_new(box)=0.
	
		!	Check density
	
		call set_rho(qrho,qpresx,qpresy,qpresz, &
			qpresxy,qpresxz,qpresyz,rmassq, &
			hrho,hpresx,hpresy,hpresz, &
			hpresxy,hpresxz,hpresyz,rmassh, &
			orho,opresx,opresy,opresz, &
			opresxy,opresxz,opresyz,rmasso, &
			epres,nx,ny,nz,n_grids,box,o_conc)
	
		!	Check speeds of individual grids
	
		vlim=vlim_factor*box

		!	Find maximum velocity to determine time step
		call set_speed_agrd( &
			nx,ny,nz,n_grids, box, rmassq,rmassh,rmasso, &
			isotropic, smallbit, &
			aniso_limit, vlim, qvlim,hvlim,ovlim, &
			qclim,hclim,oclim, cslim, &
			bx,by,bz, bx0,by0,bz0, &
			qrho, qpx,qpy,qpz, &
			qpresx,qpresy,qpresz, &
			qpresxy,qpresxz,qpresyz, &
			hrho, hpx,hpy,hpz, &
			hpresx,hpresy,hpresz,&
			hpresxy,hpresxz,hpresyz, &
			orho, opx,opy,opz, &
			opresx,opresy,opresz, &
			opresxy,opresxz,opresyz, &
			epres, &
			btot, pxmax,pymax,pzmax, csmax, alfmax, fastest )
		write(*,speeds_fmt) box,csmax,alfmax,pxmax,pymax,pzmax
	
		t_stepnew(box)=amin1( stepsz*xspac(box)/fastest, t_step_max )
	enddo
	write(*,*)'Speeds checked.'
	
	delt_old=delt
	delt=stepsz/fastest
	delt=amin1(delt,1.25*delt_old, t_step_max)
	delt=amax1(3.e-3,delt)

	write(*,init_fmt) 'Init: t=', t, tab//'dt=', delt, tab//'ut=', ut

	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) '	INITIALIZATION COMPLETED!	'
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	
	
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!		@		********************************		@
	!		@			START MAIN TIME SEQUENCE			@
	!	  	@		********************************		@
	!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	
	
	do while(t.lt.tmax)
	
		write(*,*)'Start main loop, t=', t
		
		do box=n_grids,1,-1
			t_step(box)=t_stepnew(box)
			
			!	Sync time steps
			
			t_old(box)=0.
			t_new(box)=0.
			if(box.eq.n_grids) then
				smallest_step=t_step(box)
				m_smallest=box
			else
				!	Check to see if box grid doesnt outstep grid box+1
				if( t_step(box).gt.t_step(box+1) ) t_step(box)=t_step(box+1)
				if (smallest_step.gt.t_step(box)) then
					smallest_step=t_step(box)
					m_smallest=box
				endif
			endif
		enddo

		!	Sync time steps between boxes
		do box=1,n_grids
			if(box.le.m_smallest) then
				t_step(box)=smallest_step
			else
				!	Round up as needed
				astep=(t_step(box)/t_step(box-1)+.50)
				!	Nearest integer
				nsteps=astep
				if(nsteps.gt.2) nsteps=2
				if(nsteps.lt.1) nsteps=1
				t_step(box)=t_step(box-1)*nsteps
			endif
		enddo
		
		m_step = t_step(n_grids)/t_step(m_smallest)
		delt = t_step(n_grids)
		told = t
		utold = ut
		old_tilt = tilt
		
		if(delt .lt. t_step_min) then
			write(*,*) "ERROR: Main loop time step smaller than ", &
				t_step_min, ", exiting."
			stop
		endif
		t = t + delt
		ut = t*t_equiv/3600.
		tilt = tilt + dtilt*delt

		!	Increment planet rotation and moon orbit
		nrot_planet = int( ut/planet_per )
		rot_hrs = ut - planet_per*nrot_planet
		rot_angle = 2.*pi*rot_hrs/planet_per
		sin_rot = sin(rot_angle)
		cos_rot = cos(rot_angle)
		
		nrot_moon = int( ut/moon_per )
		rot_hrs_moon = ut - moon_per*nrot_moon
		rot_angle_moon = 2.*pi*rot_hrs_moon/moon_per + moon_init_rot
		sin_rot_moon = sin(rot_angle_moon)
		cos_rot_moon = cos(rot_angle_moon)
		
		xmoon = moon_orbit_rad * cos_rot_moon
		ymoon = moon_orbit_rad * sin_rot_moon
		zmoon = 0.
		
		!	Write to fluxes.dat file
		write(fluxs_f,*) flux_header
		do box=n_grids,1,-1
			scale=rho_equiv*v_equiv*1.e5 &
			* (xspac(box)*planet_rad*re_equiv*1.e5)**2

			!	Calculate flux from torus
			call flux_counter(qpx,qpy,qpz,hpx,hpy,hpz, &
				opx,opy,opz,nx,ny,nz,n_grids,box, &
				rmassq,rmassh,rmasso, &
				scale,qflux_in,qflux_out,hflux_in,hflux_out, &
				oflux_in,oflux_out)
			write(fluxs_f,*) ut, box, qflux_in, hflux_in, &
				oflux_in, 'grid_in'
			write(fluxs_f,*) ut, box, qflux_out, hflux_out, &
				oflux_out, 'grid_out'
		enddo
		
		if(input_fluid) then
		
			!	This block is not currently functional.			
			do while (ut.ge.vut)
				nvx=nvx+1
				zvelx=-svel(nvx,1)/v_equiv
				zvely=-svel(nvx,2)/v_equiv
				zvely=-svel(nvx,2)/v_equiv+0.03
				zvelz=svel(nvx,3)/v_equiv
				vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
			enddo

			displace=0.
			do k=1,nz
				do j=1,ny
					do while( (ut.ge.future(j,k)) &
						.and.(ncount(j,k)+1.le.ncts) )
						nc=ncount(j,k)+1
						future(j,k)=bfld(nc,4)
						bxf(j,k)=-bfld(nc,1)/b_equiv
						byf(j,k)=-bfld(nc,2)/b_equiv
						bzf(j,k)=bfld(nc,3)/b_equiv
						rhof(j,k)=rplas(nc)/rho_equiv
						svxf(j,k)=-svel(nc,1)/v_equiv
						svyf(j,k)=0.00
						svzf(j,k)=0.0
						svyf(j,k)=-svel(nc,2)/v_equiv
						svyf(j,k)=-svel(nc,2)/v_equiv+0.03
						svzf(j,k)=svel(nc,3)/v_equiv
						avx=svxf(j,k)
						ncount(j,k)=nc

						if(warp) then
							b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
							b_perp=amax1(b_perp,0.1*abs(bxf(j,k)))
							ay = (j-1.) * xspac(n_grids) &
								+ grid_minvals(2,n_grids) &
								- rcraft(2) / re_equiv
							az = (k-1.) * xspac(n_grids) &
								+ grid_minvals(3,n_grids) &
								- rcraft(3) / re_equiv

							!	Assuming Bz IMF is positive on average 
							!	and that we can ignore transients
							!	Also assuming By IMF is negative
							b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
							b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
							displace=-bxf(j,k)* &
							(ay*byf(j,k)+az*bzf(j,k))/b_perp**2
							ar=distance+displace
							future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
						endif
					enddo
				enddo
			enddo
			!	Update spacecraft position and magnetic field
			dut=(ut-xcraft(4,n))/(zcraft(4,n)-utold)
			xcraft(1,n)=xcraft(1,n)+dut*(zcraft(1,n)-xcraft(1,n))
			xcraft(2,n)=xcraft(2,n)+dut*(zcraft(2,n)-xcraft(2,n))
			xcraft(3,n)=xcraft(3,n)+dut*(zcraft(3,n)-xcraft(3,n))
			xcraft(4,n)=ut
		endif
		
		!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		!		@		SPACECRAFT DATA RECORDING		@
		!		@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		
		if(spacecraft) then
			!	WARNING: Arguments to craft_handling must be updated
			!	each time num_inst is incremented, to add the new
			!	quantities.

			call craft_handling( nx, ny, nz, n_grids, ncraft, n_vals, &
				ndef_craft, deflt_rec_skp, scout, t, ut, utold, tgraph, &
				t_equiv, t_step_max_def, xmoon, ymoon, zmoon, sxyz, &
				re_equiv, moon_per, repeat_flybys, craft_datadir, &
				grid_minvals, grid_maxvals, xspac, craftnames, scfmt, &
				dat_header, num_inst, spam_limit, scdata, craftpos, bx,by,bz, &
				nskipped, recording, ntimes, cgridpt, t_step_max, &
				flyby_ref_time, xcraft, zcraft )
			
		endif	!end if(spacecraft)
		
		if(input_fluid) then
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
			
			dut=(ut-utold)/(vut-utold)
			svelx=svelx+dut*(zvelx-svelx)
			svely=0.
			srho=srho/float(nz*ny)
			
			svelx=svelx+delvx_wind*delt
			svely=svely+delvy_wind*delt
			srho=srho+delrho*delt
			
			svelz=svelz+delvz_wind*delt
			
			sbx_wind=sbx_wind+delbx_wind*delt
			sby_wind=sby_wind+delby_wind*delt
			sbz_wind=sbz_wind+delbz_wind*delt
			
			spx=srho*svelx
			spy=srho*svely
			spz=srho*svelz
			
			spress=(cs_wind**2*srho/gamma)/gamma1
			serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
			
			delay=t_equiv*distance/svelx/3600.
		endif	!end if(input_fluid)
		
		!	*******************************
		!	Main grid loop over delt/m_step
		!	*******************************
		
		do ms=1,m_step
			do box=n_grids,1,-1
		
				!	Test if grid sector needs to be moved in time
		
				yes_step=.false.
		
				if(box.eq.1) then
					yes_step=.true.
				else
					if ( (t_old(box).eq.t_new(box) ) &
						.and. ( t_new(box).lt.t_step(n_grids) ) &
						.and. ( abs(t_new(box)-t_new(1)).le.0.005*t_step(1) ) &
						) then

						yes_step=.true.
					endif
				endif
		
				!	Time step grid
		
				if(yes_step) then
					t_old(box) = t_new(box)
					t_new(box) = t_old(box) + t_step(box)
		
					write(speed_f,*)'lores t step. box/t/old/new: ',box, t, &
						t_old(box),t_new(box)
		
					delt = t_step(box)
					delt2 = delt/2.
		
					if(tilting) then
						atilt=old_tilt+(t_old(box)+delt2)*dtilt
						sin_tilt=sin(atilt*pi/180.)
						cos_tilt=cos(atilt*pi/180.)
						ut2=utold+(t_old(box)+delt2)*t_equiv/3600.
						nrot2=ut2/planet_per
						rot_hr2=ut2-nrot2*planet_per
						rot_angle=2.*pi*rot_hr2/planet_per

						call mak_dip_grd(bx0,by0,bz0,nx,ny,nz, &
							n_grids,mbndry,box,ijzero,numzero,mzero,r_inner, &
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					
					!	*******************************
					!	Store initial plasma parameters
					!	*******************************
					
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
					
					call store_array(oldepres,epres,nx,ny,nz,n_grids,box)
					call store_array(oldbx,bx,nx,ny,nz,n_grids,box)
					call store_array(oldby,by,nx,ny,nz,n_grids,box)
					call store_array(oldbz,bz,nx,ny,nz,n_grids,box)
					
					!	*********************************
					!	2-step Lax-Wendroff: Step 1
					!	Estimate fluid quantites at n+1/2
					!	*********************************
					
					!	Calculate standard MHD current: j = curl B
					
					rx=xspac(box)
					ry=xspac(box)
					rz=xspac(box)
					
					!	write(*,*)'Entering subroutine: calcur'
					call calcur(bx,by,bz,nx,ny,nz,n_grids,box,curx,cury,curz, &
						rx,ry,rz)
					
					!	Find total magnetic field
					
					bsx(:,:,:) = bx0(:,:,:,box) + bx(:,:,:,box)
					bsy(:,:,:) = by0(:,:,:,box) + by(:,:,:,box)
					bsz(:,:,:) = bz0(:,:,:,box) + bz(:,:,:,box)
					
					!	Find magnitude of B
					
					btot = amax1( sqrt( bsx**2 + bsy**2 + bsz**2 ), smallbit )
					
					!	write(*,*)'Entering subroutine: fnd_evel'
					call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
						opx,opy,opz,orho,curx,cury,curz,evx,evy,evz, &
						tvx,tvy,tvz,nx,ny,nz,n_grids,box, &
						rmassq,rmassh,rmasso,reynolds)
					
					!	write(*,*)'Entering subroutine: bande'
					call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
						curx,cury,curz,evx,evy,evz,btot, &
						epres,qrho,hrho,orho,resistive,resist,reynolds, &
						nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso, &
						ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
						rx,ry,rz)
					
					!	write(*,*)'Entering subroutine: push_elec, first pass'
					call push_elec(wrkepres,oldepres,epres,evx,evy,evz, &
						gamma,gamma1,nx,ny,nz,n_grids,box,0.5*delt, &
						rx,ry,rz)
					
					!	write(*,*)'Entering subroutine: push_ion q, first pass'
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
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:), isotropic)
					
					!	write(*,*)'Entering subrout: push_ion h, first pass'
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
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:), isotropic)
					
					!	write(*,*)'Entering subroutine: push_ion o, first pass'
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
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:), isotropic)
					
					!	write(*,*)'Entering subroutine: push_bfld'
					call push_bfld(wrkbx,wrkby,wrkbz,oldbx,oldby,oldbz, &
						efldx,efldy,efldz,nx,ny,nz,n_grids,box,0.5*delt, &
						rx,ry,rz)
					
					!	****************************************
					!			Apply boundary conditions
					!	****************************************
					
					!	write(*,*)'Entering first main Lax-Wendroff loop.'
					
					if(box.eq.n_grids) then
						!	write(*,*)'Applying boundary conditions.'
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
						if(t_grid.gt.t_new(box+1)) t_grid=t_new(box+1)
						if (t_grid.lt.t_old(box+1)) t_grid=t_old(box+1)
		
						!	write(*,*) 'Entering subroutine: bndry_flanks'
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
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif

					if(box.le.mbndry) then
						!	write(*,*)'Entering subroutine: bndry_inner'
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
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					
					!	Check fluid parameters
					
					if(input_fluid) then
						call set_imf( wrkbx,wrkby,wrkbz,bx0,by0,bz0,&
						bxp,byp,bzp, &
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
						ti_te,rho_frac,nx,ny,nz,n_grids )
					endif
					
					call set_rho(wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
						wrkqpresxy,wrkqpresxz,wrkqpresyz,rmassq, &
						wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
						wrkhpresxy,wrkhpresxz,wrkhpresyz,rmassh, &
						wrkorho,wrkopresx,wrkopresy,wrkopresz, &
						wrkopresxy,wrkopresxz,wrkopresyz,rmasso, &
						wrkepres,nx,ny,nz,n_grids,box,o_conc)
					
					vlim=vlim_factor*box
					call set_speed_agrd( &
						nx,ny,nz,n_grids, box, rmassq,rmassh,rmasso, &
						isotropic, smallbit, &
						aniso_limit, vlim, qvlim,hvlim,ovlim, &
						qclim,hclim,oclim, cslim, &
						wrkbx,wrkby,wrkbz, bx0,by0,bz0, &
						wrkqrho, wrkqpx,wrkqpy,wrkqpz, &
						wrkqpresx,wrkqpresy,wrkqpresz, &
						wrkqpresxy,wrkqpresxz,wrkqpresyz, &
						wrkhrho, wrkhpx,wrkhpy,wrkhpz, &
						wrkhpresx,wrkhpresy,wrkhpresz, &
						wrkhpresxy,wrkhpresxz,wrkhpresyz, &
						wrkorho, wrkopx,wrkopy,wrkopz, &
						wrkopresx,wrkopresy,wrkopresz, &
						wrkopresxy,wrkopresxz,wrkopresyz, &
						wrkepres, &
						btot, pxmax,pymax,pzmax, csmax, alfmax, fastest )

					!	*******************************
					!	Lax-Wendroff: Step 2
					!	Use the predicted value to find
					!	corrected value for n+1
					!	*******************************

					if(tilting) then
						atilt=old_tilt+(t_old(box)+delt)*dtilt
						sin_tilt=sin(atilt*pi/180.)
						cos_tilt=cos(atilt*pi/180.)
						ut2=utold+(t_old(box)+delt)*t_equiv/3600.
						nrot2=ut2/planet_per
						rot_hr2=ut2-nrot2*planet_per
						rot_angle=2.*pi*rot_hr2/planet_per
						!	write(*,*)'mak_dip full with',t_old(box),delt
						call mak_dip_grd(bx0,by0,bz0,nx,ny,nz,n_grids, &
							mbndry,box,ijzero,numzero,mzero,r_inner, &
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					rx=xspac(box)
					ry=xspac(box)
					rz=xspac(box)
					
					!	Calculate standard MHD current: j = curl B
					
					call calcur(wrkbx,wrkby,wrkbz,nx,ny,nz,n_grids,box, &
						curx,cury,curz,rx,ry,rz)
					
					!	Find total magnetic field
					
					bsx(:,:,:) = bx0(:,:,:,box) + wrkbx(:,:,:,box)
					bsy(:,:,:) = by0(:,:,:,box) + wrkby(:,:,:,box)
					bsz(:,:,:) = bz0(:,:,:,box) + wrkbz(:,:,:,box)
					
					!	Find magnitude of B
					
					btot = amax1( sqrt( bsx**2 + bsy**2 + bsz**2 ), smallbit )
					
					!	Find the electric field from electron momentum eqn
					
					call fnd_evel(wrkqpx,wrkqpy,wrkqpz,wrkqrho, &
						wrkhpx,wrkhpy,wrkhpz,wrkhrho, &
						wrkopx,wrkopy,wrkopz,wrkorho, &
						curx,cury,curz,evx,evy,evz, &
						tvx,tvy,tvz,nx,ny,nz,n_grids,box, &
						rmassq,rmassh,rmasso,reynolds)
					
					call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
						curx,cury,curz,evx,evy,evz,btot, &
						wrkepres,wrkqrho,wrkhrho,wrkorho, &
						resistive,resist,reynolds, &
						nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso, &
						ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
						rx,ry,rz)
					
					!	write(*,*)'Entering subroutine: push_*, second pass'
					call push_elec(wrkepres,oldepres,epres,evx,evy,evz, &
						gamma,gamma1,nx,ny,nz,n_grids,box,delt,rx,ry,rz)
					
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
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:), isotropic)
					
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
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:), isotropic)
					
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
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:), isotropic)
					
					call push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
						efldx,efldy,efldz,nx,ny,nz,n_grids,box,delt, &
						rx,ry,rz)
					
					!	***************************************
					!			Apply boundary conditions
					!	***************************************
					
					!write(*,*)'Entering second main Lax-Wendroff loop.'
					
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
						if(t_grid.gt.t_new(box+1)) t_grid=t_new(box+1)
						if (t_grid.lt.t_old(box+1)) t_grid=t_old(box+1)
						
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
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					
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
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					
					if(input_fluid) then
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
					
					call set_rho(qrho,qpresx,qpresy,qpresz, &
						qpresxy,qpresxz,qpresyz,rmassq, &
						hrho,hpresx,hpresy,hpresz, &
						hpresxy,hpresxz,hpresyz,rmassh, &
						orho,opresx,opresy,opresz, &
						opresxy,opresxz,opresyz,rmasso, &
						epres,nx,ny,nz,n_grids,box,o_conc)
					
					vlim=vlim_factor*box
					call set_speed_agrd( &
						nx,ny,nz,n_grids, box, rmassq,rmassh,rmasso, &
						isotropic, smallbit, &
						aniso_limit, vlim, qvlim,hvlim,ovlim, &
						qclim,hclim,oclim, cslim, &
						bx,by,bz, bx0,by0,bz0, &
						qrho, qpx,qpy,qpz, &
						qpresx,qpresy,qpresz, &
						qpresxy,qpresxz,qpresyz, &
						hrho, hpx,hpy,hpz, &
						hpresx,hpresy,hpresz,&
						hpresxy,hpresxz,hpresyz, &
						orho, opx,opy,opz, &
						opresx,opresy,opresz, &
						opresxy,opresxz,opresyz, &
						epres, &
						btot, pxmax,pymax,pzmax, csmax, alfmax, fastest )
					
					!	......................
					!	Try Lapidus smoothing:
					!	smoothed results will
					!	appear in nt2
					!	......................
					
					rx=xspac(box)
					ry=xspac(box)
					rz=xspac(box)

					!	Species q
					
					curx(:,:,:) = qpx(:,:,:,box) &
						/ amax1( qrho(:,:,:,box), smallbit )
					cury(:,:,:) = qpy(:,:,:,box) &
						/ amax1( qrho(:,:,:,box), smallbit )
					curz(:,:,:) = qpz(:,:,:,box) &
						/ amax1( qrho(:,:,:,box), smallbit )
					call lap_plasma(qrho,qpx,qpy,qpz, &
						qpresx,qpresy,qpresz, &
						qpresxy,qpresxz,qpresyz, &
						wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
						wrkqpresx,wrkqpresy,wrkqpresz, &
						wrkqpresxy,wrkqpresxz,wrkqpresyz, &
						curx,cury,curz, &
						nx,ny,nz,n_grids,box, &
						chirho,chipxyz,chierg,delt,isotropic)
					
					!	Species h
					
					curx(:,:,:) = hpx(:,:,:,box) &
						/ amax1( hrho(:,:,:,box), smallbit )
					cury(:,:,:) = hpy(:,:,:,box) &
						/ amax1( hrho(:,:,:,box), smallbit )
					curz(:,:,:) = hpz(:,:,:,box) &
						/ amax1( hrho(:,:,:,box), smallbit )
					call lap_plasma(hrho,hpx,hpy,hpz, &
						hpresx,hpresy,hpresz, &
						hpresxy,hpresxz,hpresyz, &
						wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
						wrkhpresx,wrkhpresy,wrkhpresz, &
						wrkhpresxy,wrkhpresxz,wrkhpresyz, &
						curx,cury,curz, &
						nx,ny,nz,n_grids,box, &
						chirho,chipxyz,chierg,delt,isotropic)

					!	Species o
					
					curx(:,:,:) = opx(:,:,:,box) &
						/ amax1( orho(:,:,:,box), smallbit )
					cury(:,:,:) = opy(:,:,:,box) &
						/ amax1( orho(:,:,:,box), smallbit )
					curz(:,:,:) = opz(:,:,:,box) &
						/ amax1( orho(:,:,:,box), smallbit )
					call lap_plasma(orho,opx,opy,opz, &
						opresx,opresy,opresz, &
						opresxy,opresxz,opresyz, &
						wrkorho,wrkopx,wrkopy,wrkopz, &
						wrkopresx,wrkopresy,wrkopresz, &
						wrkopresxy,wrkopresxz,wrkopresyz, &
						curx,cury,curz, &
						nx,ny,nz,n_grids,box, &
						chirho,chipxyz,chierg,delt,isotropic)
					
					!	Electrons

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

					!	Magnetic fields
					
					call lap_bfld(bx,by,bz, &
						wrkbx,wrkby,wrkbz, &
						curx,cury,curz, &
						nx,ny,nz,n_grids,box, &
						chirho,chipxyz,chierg,delt)
					
					!	Reset boundary conditions
					
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
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					if(box.le.mbndry) then
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
							re_equiv,r_inner,sbx_wind,spress, &
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					
					if(input_fluid) then
						call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0, &
							bxp,byp,bzp, &
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
					
					vlim=vlim_factor*box
					call set_speed_agrd( &
						nx,ny,nz,n_grids, box, rmassq,rmassh,rmasso, &
						isotropic, smallbit, &
						aniso_limit, vlim, qvlim,hvlim,ovlim, &
						qclim,hclim,oclim, cslim, &
						wrkbx,wrkby,wrkbz, bx0,by0,bz0, &
						wrkqrho, wrkqpx,wrkqpy,wrkqpz, &
						wrkqpresx,wrkqpresy,wrkqpresz, &
						wrkqpresxy,wrkqpresxz,wrkqpresyz, &
						wrkhrho, wrkhpx,wrkhpy,wrkhpz, &
						wrkhpresx,wrkhpresy,wrkhpresz, &
						wrkhpresxy,wrkhpresxz,wrkhpresyz, &
						wrkorho, wrkopx,wrkopy,wrkopz, &
						wrkopresx,wrkopresy,wrkopresz, &
						wrkopresxy,wrkopresxz,wrkopresyz, &
						wrkepres, &
						btot, pxmax,pymax,pzmax, csmax, alfmax, fastest )
					
					!	write(*,*) 'Lapidus smoothing complete.'
					
					!	....................................
					!	Add a little bit of flux correction:
					!	....................................
					
					!	write(*,*)'Entering subroutine: flux_correct'
					call flux_correct(qrho,qpresx,qpresy,qpresz, &
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
					
					!	Set boundary conditions
					
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
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					
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
							grid_minvals(1,:), grid_maxvals(1,:), &
							grid_minvals(2,:), grid_maxvals(2,:), &
							grid_minvals(3,:), grid_maxvals(3,:) )
					endif
					
					if(input_fluid) then
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
					
					call set_rho(qrho,qpresx,qpresy,qpresz, &
						qpresxy,qpresxz,qpresyz,rmassq, &
						hrho,hpresx,hpresy,hpresz, &
						hpresxy,hpresxz,hpresyz,rmassh, &
						orho,opresx,opresy,opresz, &
						opresxy,opresxz,opresyz,rmasso, &
						epres,nx,ny,nz,n_grids,box,o_conc)
					
					vlim=vlim_factor*box
					call set_speed_agrd( &
						nx,ny,nz,n_grids, box, rmassq,rmassh,rmasso, &
						isotropic, smallbit, &
						aniso_limit, vlim, qvlim,hvlim,ovlim, &
						qclim,hclim,oclim, cslim, &
						bx,by,bz, bx0,by0,bz0, &
						qrho, qpx,qpy,qpz, &
						qpresx,qpresy,qpresz, &
						qpresxy,qpresxz,qpresyz, &
						hrho, hpx,hpy,hpz, &
						hpresx,hpresy,hpresz,&
						hpresxy,hpresxz,hpresyz, &
						orho, opx,opy,opz, &
						opresx,opresy,opresz, &
						opresxy,opresxz,opresyz, &
						epres, &
						btot, pxmax,pymax,pzmax, csmax, alfmax, fastest )
					
					t_stepnew(box)= amin1( stepsz*xspac(box)/fastest, &
						t_step_max )

					!	write(*,*)'needed step of',t_step(box),t_stepnew(box)
					
					!	Sync time steps if needed, and apply core conditions
					
					if(box.lt.n_grids) then
						if( abs(t_new(box)-t_new(box+1)) &
							.le. 0.005*t_step(box) ) then

							call bndry_corer( &
								qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
								hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
								orho,opresx,opresy,opresz,opx,opy,opz, &
								epres,bx,by,bz, &
								qpresxy,qpresxz,qpresyz, &
								hpresxy,hpresxz,hpresyz, &
								opresxy,opresxz,opresyz, &
								nx,ny,nz,n_grids,box, &
								grid_minvals(1,:), grid_maxvals(1,:), &
								grid_minvals(2,:), grid_maxvals(2,:), &
								grid_minvals(3,:), grid_maxvals(3,:) )
							t_old(box+1)=t_new(box+1)
						endif
					endif
					
				endif	! yes_step of main grid
			enddo	! lores box increment
		enddo	! box sweep of lores grid
		
		!	write(*,*) 'Final sync on boundary conditions.'
		t_old(1)=t_new(1)

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
				grid_minvals(1,:), grid_maxvals(1,:), &
				grid_minvals(2,:), grid_maxvals(2,:), &
				grid_minvals(3,:), grid_maxvals(3,:) )
			
			call set_rho(qrho,qpresx,qpresy,qpresz, &
				qpresxy,qpresxz,qpresyz,rmassq, &
				hrho,hpresx,hpresy,hpresz, &
				hpresxy,hpresxz,hpresyz,rmassh, &
				orho,opresx,opresy,opresz, &
				opresxy,opresxz,opresyz,rmasso, &
				epres,nx,ny,nz,n_grids,box,o_conc)
			t_old(box)=t_old(box-1)
			t_new(box)=t_new(box-1)
		enddo
		
		if(t.ge.tgraph) then
			
			!	Plot plasma properties
		
			write(*,*) 'Graphing data.'
			call write_graphing_data( &
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
				r_inner, run_name, dummy_fg, diagnostics, nplots, &
				isotropic, smallbit)
		
			!	NCAR graphics disabled.
			!	call visual(qrho,qpresx,qpresy,qpresz,qpresxy, &
			!		qpresxz,qpresyz,qpx,qpy,qpz,rmassq, &
			!		hrho,hpresx,hpresy,hpresz,hpresxy, &
			!		hpresxz,hpresyz,hpx,hpy,hpz,rmassh, &
			!		orho,opresx,opresy,opresz,opresxy, &
			!		opresxz,opresyz,opx,opy,opz,rmasso, &
			!		epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
			!		curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
			!		tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
			!		nx,ny,nz,n_grids,xspac, &
			!		cross,along,flat,xcraft,ncraft,re_equiv, &
			!		grid_minvals(1,:), grid_maxvals(1,:), &
			!		grid_minvals(2,:), grid_maxvals(2,:), &
			!		grid_minvals(3,:), grid_maxvals(3,:), &
			!		ut,b_equiv,ti_te,rho_equiv)
		
			tgraph=tgraph+deltg
			write(*,*) 'Graphics plotted. Next tgraph at t = ', tgraph
		endif !if(t.ge.tgraph)
		
		if(t.ge.tinj.and.ringo) then
		
			!	Inject torus plasma
		
			write(*,*) 'Injecting plasma torus, tinj:', tinj
			
			box=1
			dx = ( grid_maxvals(1,box) - grid_minvals(1,box) ) / (nx-1.)
			dy = ( grid_maxvals(2,box) - grid_minvals(2,box) ) / (ny-1.)
			dz = ( grid_maxvals(3,box) - grid_minvals(3,box) ) / (nz-1.)
			
			tot_o=0.
			tot_h=0.
			tot_q=0.
			
			bsx(:,:,:) = bx0(:,:,:,box) + bx(:,:,:,box)
			bsy(:,:,:) = by0(:,:,:,box) + by(:,:,:,box)
			bsz(:,:,:) = bz0(:,:,:,box) + bz(:,:,:,box)
			
			do k=1,nz
				az = grid_minvals(3,box) + dz * (k-1) - zdip

				do j=1,ny
					ay = grid_minvals(2,box) + dy * (j-1) - ydip

					do i=1,nx
						ax = grid_minvals(1,box) + dx * (i-1) - xdip
						
						rx=ax*re_equiv
						ry=ay*re_equiv
						rd=sqrt(rx**2+ry**2)
						!	Radial distance from planet center
						ar=sqrt( (ax+xdip)**2 + (ay+ydip)**2 &
							+ (az+zdip)**2 )
						
						rvy=rx*v_rot
						rvx=-ry*v_rot
						rvz=0.
						corotate=sqrt(rvx**2+rvy**2)
						
						ar_tmid=sqrt(ax**2+ay**2)*re_equiv
						dr_tmid=abs(ar_tmid - torus_inj_dist)
						
						!	MAT - Corrected, according to 
						!	doi:10.1029/2009JE003372
						!	and
						!	doi:10.1029/JA091iA08p08749
						
						if( (abs(dr_tmid).lt.2.*torus_rad) &
							.and. (ar.gt.srf_bndry) ) then
							
							! scale height in # planet_rad
							rscale=exp(-((dr_tmid)/(0.5*torus_rad))**2) 
							zscale=exp( -((az*re_equiv)&
								/ (0.5*torus_rad))**2 )
							abtot = ( bsx(i,j,k)**2+bsy(i,j,k)**2 ) &
								/ bsz(i,j,k)**2
							dscale=amax1(0.,1.-10.*abtot)
							
							! Inject h species
							
							hden=denh_torus*rmassh*rscale*zscale*dscale
							hrho(i,j,k,box)=hrho(i,j,k,box)+hden
							!	Temp. is proportional to v**2
							del_hp=hden*(corotate**2)*t_torus
							hpresx(i,j,k,box) = hpresx(i,j,k,box)+del_hp
							hpresy(i,j,k,box) = hpresy(i,j,k,box)+del_hp
							hpresz(i,j,k,box) = hpresz(i,j,k,box) &
								+ del_hp*aniso_factor
							hpresxy(i,j,k,box)=hpresxy(i,j,k,box)+del_hp
							
							hpx(i,j,k,box)=hpx(i,j,k,box) &
								+ reduct*hden*rvx
							hpy(i,j,k,box)=hpy(i,j,k,box) &
								+ reduct*hden*rvy
							
							! Inject q species
							
							qden=denq_torus*rmassq*rscale*zscale*dscale 
							qrho(i,j,k,box)=qrho(i,j,k,box)+qden
							!	Temp. is proportional to v**2
							del_qp=(qden/rmassq)*(corotate**2)*t_torus
							qpresx(i,j,k,box) = qpresx(i,j,k,box)+del_qp
							qpresy(i,j,k,box) = qpresy(i,j,k,box)+del_qp
							qpresz(i,j,k,box) = qpresz(i,j,k,box) &
								+ del_qp*aniso_factor
							qpresxy(i,j,k,box)=qpresxy(i,j,k,box)+del_qp
							
							qpx(i,j,k,box)=qpx(i,j,k,box) &
								+ reduct*qden*rvx
							qpy(i,j,k,box)=qpy(i,j,k,box) &
								+ reduct*qden*rvy
							
							!	Inject o species 
							
							oden=deno_torus*rmasso*rscale*zscale*dscale 
							orho(i,j,k,box)=orho(i,j,k,box)+oden
							!	Temp. is proportional to v**2
							del_op=(oden/rmasso)*(corotate**2)*t_torus 
							opresx(i,j,k,box) = opresx(i,j,k,box)+del_op
							opresy(i,j,k,box) = opresy(i,j,k,box)+del_op
							opresz(i,j,k,box) = opresz(i,j,k,box) &
								+ del_op*aniso_factor
							opresxy(i,j,k,box)=opresxy(i,j,k,box)+del_op
							
							opx(i,j,k,box)=opx(i,j,k,box) &
								+ reduct*oden*rvx
							opy(i,j,k,box)=opy(i,j,k,box) &
								+ reduct*oden*rvy

							!	Equal temp. for electrons
							epres(i,j,k,box)=epres(i,j,k,box) &
								+ (del_op+del_hp+del_qp)
							
							tot_q = tot_q + qden*dx*dy*dz
							tot_h = tot_h + hden*dx*dy*dz
							tot_o = tot_o + oden*dx*dy*dz
						endif
					enddo
				enddo
			enddo
			
			!	Scale factors to kg/s
			
			!	volume in cubic meters
			volume=(re_equiv*planet_rad*1.e3)**3
			atime=deltinj*t_equiv
			write(*,*)'ut,volume,t_equiv,atime',ut,volume,t_equiv,atime
			write(concs_f,*)'ut,volume,t_equiv,atime',ut,volume, &
				t_equiv,atime
			
			tot_q=tot_q*volume/atime*rho_equiv*1.e6/rmassq
			tot_h=tot_h*volume/atime*rho_equiv*1.e6/rmassh
			tot_o=tot_o*volume/atime*rho_equiv*1.e6/rmasso
			write(*,*)'tot torus ions/s q,h,o',ut,tot_q,tot_h,tot_o
			write(concs_f,*)'tot torus ions/s q,h,o',ut,tot_q,tot_h, &
				tot_o
			
			tot_q=tot_q*m_prot*rmassq
			tot_h=tot_h*m_prot*rmassh
			tot_o=tot_o*m_prot*rmasso
			write(*,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
			write(concs_f,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
			
			tinj=tinj+deltinj
			
		endif !if(t.ge.tinj.and.ringo)
		
		if(t.ge.ts1) then
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
			write(*,*) 'Writing to fluid file. t =', t
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
			write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
			
			write(fname_fluid,'(a5,i2.2)') fluid_pre, fluid_f+offset_f
			open(fluid_f,file=fname_fluid,status='unknown', &
				form='unformatted')
			
				!	Write restart data
			
				write(fluid_f)t
				write(fluid_f)qrho
				write(fluid_f)qpx
				write(fluid_f)qpy
				write(fluid_f)qpz
				write(fluid_f)qpresx
				write(fluid_f)qpresy
				write(fluid_f)qpresz
				write(fluid_f)qpresxy
				write(fluid_f)qpresxz
				write(fluid_f)qpresyz
				write(fluid_f)hrho
				write(fluid_f)hpx
				write(fluid_f)hpy
				write(fluid_f)hpz
				write(fluid_f)hpresx
				write(fluid_f)hpresy
				write(fluid_f)hpresz
				write(fluid_f)hpresxy
				write(fluid_f)hpresxz
				write(fluid_f)hpresyz
				write(fluid_f)orho
				write(fluid_f)opx
				write(fluid_f)opy
				write(fluid_f)opz
				write(fluid_f)opresx
				write(fluid_f)opresy
				write(fluid_f)opresz
				write(fluid_f)opresxy
				write(fluid_f)opresxz
				write(fluid_f)opresyz
				write(fluid_f)bx
				write(fluid_f)by
				write(fluid_f)bz
				write(fluid_f)epres
				write(fluid_f)bx0
				write(fluid_f)by0
				write(fluid_f)bz0
				write(fluid_f)parm_srf,parm_mid,parm_zero, &
					ijzero,numzero,ijmid,nummid,ijsrf,numsrf
			close(fluid_f)
			
			fluid_f = fluid_f + 1

			if(fluid_f.gt.20) then
				fluid_f = 11
				offset_f = offset_f + 10
				!	Prevent overwriting fluid files if we run long
				!	enough to write 100 files:
				if(offset_f.gt.89) tmax = t-1.
			endif
			ts1=ts1+tsave
			
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
				enddo	! end box loop
				
				!	Apply boundary conditions
				
				do box=n_grids-1,1,-1
					call flanks_synced(bx,nx,ny,nz,n_grids,box, &
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:) )
					call flanks_synced(by,nx,ny,nz,n_grids,box, &
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:) )
					call flanks_synced(bz,nx,ny,nz,n_grids,box, &
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
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
						grid_minvals(1,:), grid_maxvals(1,:), &
						grid_minvals(2,:), grid_maxvals(2,:), &
						grid_minvals(3,:), grid_maxvals(3,:) )
				enddo  !end bndry_corer
			endif  ! end divb_lores

		endif	!end if(t.ge.ts1)
	enddo	!end do while(t.lt.tmax)
	
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) 'Run complete! t =', t, 'ut = ', ut
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	
	tmax = int(t - utstart*3600./t_equiv) + 0.05
	utstart = ut
	start=.false.
	
	write(*,*) 'Writing input parameters to ',fname_inp_o
	open(inp_o_f,file=trim(fname_inp_o),status='unknown', &
		form='formatted')
		write(inp_o_f,option)
		write(inp_o_f,planet)
		write(inp_o_f,speeds)
		write(inp_o_f,windy)
		write(inp_o_f,physical)
		write(inp_o_f,smooth)
	close(inp_o_f)
	
	write(*,*) '-'
	call system('date +"FINISH: %T-%Y-%m-%d"')
	write(*,*) '-'

end program multifluid
