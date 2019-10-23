!	For generating G3M³ input files from relatable units
!	Author: Marshall J. Styczinski
!	

program make_input
	implicit none

	integer, parameter :: n_boxes = 6
	integer, parameter :: pts_per_box = 50
	integer, parameter :: box_scale = 2
	real, parameter :: z_fraction = 0.48
	!character*20, parameter :: box_nums = '()'

	integer, parameter :: dp = selected_real_kind(17,300)
	real(dp),parameter :: pi = 4.0*atan(1.0)
	real(dp),parameter :: c = 2.99792458e8_dp
	real(dp),parameter :: qe = 1.602e-19_dp
	real(dp),parameter :: melec = 9.1e-31_dp
	real(dp),parameter :: mprot = 1.67e-27_dp
	real(dp),parameter :: k_B = 1.38e-23_dp
	real(dp),parameter :: mu_o = pi*4e-7
	real(dp),parameter :: eps_o = 8.85e-12_dp
	character,parameter :: tab=char(9)

	integer :: input_f = 20
	character*10 :: input_fname = 'fmpd3din'

	real B_inner, B_wind
	real T_iono, m_iono, eden_iono
	real T_wind, m_wind, amuden_wind
	real Te

	integer box
	real grid_max

	!	Collision rates involve o and h species, (ions/neutral) and electrons.
	!	Letters indicate o species/h species or ion/neutral/electron involved in collisions.
	!	Expressions all come from Schunk and Nagy 2009.
	real cond_ions, cond_neut
	real rate_oi_hi1, rate_oi_hi2
	real rate_oi_hn
	real rate_oi_e1, rate_hi_e1
	real rate_oi_e2, rate_hi_e2
	real rate_hn_e

	real den_hn, h_frac, o_frac
	real Te_3halvs
	real :: hard_coded_hrho_factor = 10.	! This number is found in saturn3d_var.f90 in the definition of hrho0, at about line 874.
	! Sets ionospheric h number density to be this number times the inner boundary number density (den_earth, which sets q species number density).
	! In constrast, the o number density is set by multiplying o_conc by den_earth.

	! group 'option':
	real tmax, stepsz, tsave
	integer ntgraf
	logical start
	! group 'earth'
	real xdip, ydip, zdip, rearth, &
		tilt1, tilt2, &
		rmassq, rmassh, rmasso
	logical tilting
	! group 'speeds'
	real cs_inner, alf_inner1, alf_inner2, &
		alpha_e, den_earth, o_conc, &
		gravity, ti_te, gamma
	logical ringo, update, divb_on
	! group 'windy'
	real re_wind, cs_wind, &
		vx_wind1, vx_wind2, vy_wind1, vy_wind2, vz_wind1, vz_wind2, &
		alfx_wind1, alfx_wind2, alfy_wind1, alfy_wind2, alfz_wind1, &
		alfz_wind2, den_wind1, den_wind2, &
		reynolds, resist, rho_frac, bfrac, vfrac
	! group 'physical' 
	real re_equiv, b_equiv, v_equiv, rho_equiv
	real utstart
	logical spacecraft, warp
	! group 'smooth'
	real chirho, chipxyz, chierg, difrho, difpxyz, diferg

	namelist/option/tmax,ntgraf,stepsz,start,tsave
	namelist/earth/xdip,ydip,zdip,rearth, &
		tilt1,tilt2,tilting,rmassq,rmassh,rmasso
	namelist/speeds/cs_inner,alf_inner1,alf_inner2, &
		alpha_e,den_earth,o_conc, &
		gravity,ti_te,gamma,ringo,update,divb_on
	namelist/windy/re_wind,cs_wind,vx_wind1,vx_wind2, &
		vy_wind1,vy_wind2,vz_wind1,vz_wind2, &
		alfx_wind1,alfx_wind2, &
		alfy_wind1,alfy_wind2, &
		alfz_wind1,alfz_wind2, &
		den_wind1,den_wind2, &
		reynolds,resist,rho_frac,bfrac,vfrac
	namelist/physical/re_equiv,b_equiv,v_equiv,rho_equiv, &
		spacecraft,warp,utstart
	namelist/smooth/chirho,chipxyz,chierg, &
		difrho,difpxyz,diferg

! Executions

	!option
	tmax=10000.01
	ntgraf=20
	stepsz=0.08
	start=.true.
	tsave=200.

	!earth
	xdip=0.00001
	ydip=0.00000
	zdip=0.000000
	rearth=20.0	! in grid points contained within a radius
	tilt1=0.00
	tilt2=0.00
	tilting =.false. 
	rmassq=1.
	rmassh=16.
	rmasso=32.

	!physical
	re_equiv=1./rearth
	b_equiv=65.75	! in nT
	v_equiv=1000.	! in km/s
	rho_equiv=2.
	spacecraft=.false.
	warp=.false.
	utstart=0.00

	T_iono = 600.	! T in kelvin of ionosphere at inner boundary
	m_iono = 30.2	! Average m/q in amu/e of ionosphere species. Estimated from Rubin et al. 2015 figures indicating ~5000 O2+/cc and ~600 O+/cc
	eden_iono = 5000.	! Free charge density at ionosphere boundary
	B_inner = 450.	! Magnitude of magnetic field in nT at ionosphere boundary

	!speeds
	cs_inner = sqrt( 3.*k_B*T_iono/(m_iono*mprot) ) / 1000./v_equiv
	alf_inner1 = B_inner/b_equiv / sqrt(m_iono*eden_iono/rho_equiv)	! mu_o is 1 in sim units
	alf_inner2 = alf_inner1 + 0.
	alpha_e=6.0
	o_conc=hard_coded_hrho_factor/5.	! How should o density compare to h density at the inner boundary? This value taken from the same figures mentioned about from Rubin et al. 2015.
	den_earth = eden_iono / (hard_coded_hrho_factor + o_conc + 1) / rho_equiv	! This sets qden at the inner boundary, so we must divide by the additional h and o relative to qden.
	gravity=1.3 / (rearth*re_equiv)**2	! Number is surface gravity in m/s^2
	ti_te=7.
	gamma=1.6666
	ringo=.false.
	update=.false.
	divb_on=.true.

	T_wind = 130.*qe/k_B	! Number is in eV
	m_wind = 18.
	amuden_wind = 2000.
	B_wind = 450.

	!windy
	re_wind=35.
	cs_wind=sqrt(3.*k_B*T_wind/(m_wind*mprot)) / 1000./v_equiv
	vx_wind1 = 100. / v_equiv	! Number is in km/s
	vx_wind2 = vx_wind1 + 0.
	vy_wind1=0.000
	vy_wind2=0.000
	vz_wind1=0.0
	vz_wind2=0.0
	alfx_wind1=0.0000
	alfx_wind2=0.0000
	alfy_wind1=0.00
	alfy_wind2=0.00
	alfz_wind1 = B_wind/b_equiv / sqrt(amuden_wind/rho_equiv)
	alfz_wind2 = alfz_wind1 + 0.
	rho_frac=17./12.	! This must be calculated such that (rmassq + rho_frac*rmassh + rho_frac*rmasso) / (2*rho_frac+1) gives the average wind ion mass, m_wind.
	den_wind1 = amuden_wind / m_wind / (2*rho_frac + 1) / rho_equiv	! This sets qden at the upstream boundary, so we must divide by 2*rho_frac+1 so that qden*rho_frac gives correct hden and oden values.
	den_wind2 = den_wind1 + 0.
	reynolds=64.0
	bfrac=1.0
	vfrac=1.0

	den_hn = 5.e8	! From Rubin et al 2015 section 2.6
	o_frac = o_conc / ( o_conc + hard_coded_hrho_factor + 1. )
	h_frac = hard_coded_hrho_factor / ( o_conc + hard_coded_hrho_factor + 1. )
	Te = T_iono/ti_te
	Te_3halvs = sqrt(Te)**3

	!	These calculations are derived from Rubin et al 2015, who copied them from Schunk and Nagy 2009. The number densities are in 1/m^3 and the masses are all in amu, including electron mass.
	rate_oi_e1 = 54.5e-6 * o_frac*eden_iono*1.e6 / Te_3halvs
	rate_hi_e1 = 54.5e-6 * h_frac*eden_iono*1.e6 / Te_3halvs
	rate_oi_e2 = 1.27e-6 * (sqrt(melec/mprot)*eden_iono*1.e6) &
		/ (rmasso * Te_3halvs)
	rate_hi_e2 = 1.27e-6 * (sqrt(melec/mprot)*eden_iono*1.e6) &
		/ (rmassh * Te_3halvs)
	rate_hn_e = 1.82e-16 * den_hn*1.e6 * sqrt(Te) * (1. + 3.6e-2*sqrt(Te))
	cond_ions = (qe**2 * eden_iono) &
		/ (melec * (rate_oi_e1 + rate_oi_e2 + rate_hi_e1 + rate_hi_e2) )
	cond_neut = (qe**2 * eden_iono) &
		/ (melec * (rate_hn_e) )
	resist = 1./cond_ions + 1./cond_neut

	!smooth
	chirho=2.0
	chipxyz=2.0
	chierg=1.0    
    difrho=0.01
	difpxyz=0.008
	diferg=0.0025

	open(input_f,file=trim(input_fname),status='unknown',form='formatted')
	
	write(input_f,option)
	write(input_f,earth)
	write(input_f,speeds)
	write(input_f,windy)
	write(input_f,physical)
	write(input_f,smooth)

	do box=1, n_boxes
		grid_max = 0.5*pts_per_box*(2**box)
		write(input_f,*) -grid_max, grid_max, -grid_max, grid_max, &
			-z_fraction*grid_max, z_fraction*grid_max, 2**(box-1), &
			mod(box+1,n_boxes+1), box-1
	enddo

	write(*,*) 'Input file written to fname: ', input_fname

end program make_input
