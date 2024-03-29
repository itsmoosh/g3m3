!
!	Measured physical parameters and references at J2000.0.
!	Source: NASA Horizons ephemeris tables:
!	https://ssd.jpl.nasa.gov/horizons.cgi
!
!	Except for dipole tilt/offsets.
!	Earth obtained from:
!	https://www.ngdc.noaa.gov/geomag/data/poles/NP.xy
!	Jupiter obtained from:
!	https://radiojove.gsfc.nasa.gov/education/educ/jupiter/basics/jfacts.htm
!
!	Distances are in km
!	Times are in hours
!	Everything else in SI
!	Unless otherwise specified
!
module astrometry

	!	General
	real,parameter :: pi		= 4.0*atan(1.0)
	real,parameter :: mu0		= 4.E-7*pi
	real,parameter :: m_prot	= 1.67262E-27
	real,parameter :: q_elec	= 1.60218E-19
	real,parameter :: k_B		= 1.38065E-23
	real,parameter :: J2000		= 2451545.0
	real,parameter :: vsmall	= 1E-5
	real,parameter :: zero		= 0.

	!	Jupiter
	real,parameter :: jupiter_orbit_rad	= 7.78557619E+08
	real,parameter :: jupiter_year		= 4332.820*24.
	real,parameter :: jupiter_rad		= 71492.
	real,parameter :: jupiter_per		= 9.92491245
	real,parameter :: jupiter_mass		= 1898.13E+24
	real,parameter :: jupiter_obliq		= 3.12*pi/180.
	real,parameter :: jupiter_incl		= 6.0910615*pi/180.
	!	jupiter_r_rot units: re_equiv, where corotation stops
	real,parameter :: jupiter_r_rot		= 30.0
	!	jupiter_torus_infall: Torus infall allowance
	real,parameter :: jupiter_torus_infall = 0.075
	real,parameter :: jupiter_tilt		= 9.6*pi/180.
	real,parameter :: jupiter_tilt_long	= 201.7*pi/180.
	real,parameter :: jupiter_dip_long	= 148.6*pi/180.
	real,parameter :: jupiter_dip_z		= 98.0*pi/180.
	!	jupiter_dip_off units: planetary radii
	real,parameter :: jupiter_dip_off	= 0.131
	real,parameter :: jupiter_init_long	= atan(2.93857583/4.0011774)
	real,parameter :: jupiter_xdip		= jupiter_dip_off * &
		sin(jupiter_dip_z) * cos(jupiter_tilt_long - jupiter_dip_long)
	real,parameter :: jupiter_ydip		= jupiter_dip_off * &
		sin(jupiter_dip_z) * sin(jupiter_tilt_long - jupiter_dip_long)
	real,parameter :: jupiter_zdip		= jupiter_dip_off * &
		cos(jupiter_dip_z)

	!	Moons: Io
		!	io_orbit_rad units: jupiter_rad
		real,parameter :: io_orbit_rad	= 421.769E+03/jupiter_rad
		real,parameter :: io_sid_per	= 1.769138*24.
		!	io_per units: solar synodic period because solar wind is
		!	from a set direction in simulation coordinates.
		real,parameter :: io_per		= io_sid_per * ( 1 + &
			io_sid_per/jupiter_year )
		real,parameter :: io_rad		= 1821.3
		real,parameter :: io_mass		= 893.3E+20
		real,parameter :: io_incl		= 0.036*pi/180.
		!	io_init_rot units: Number of hours into orbit, past midnight
		!	line of parent body, at J2000.0
		real,parameter :: io_init_rot	= ( &
			atan(8.6409473E-04/2.6719244E-03) - jupiter_init_long + &
			2.0*pi ) / (2.0*pi) * io_per


	!	Moons: Europa
		!	eur_orbit_rad units: jupiter_rad
		real,parameter :: eur_orbit_rad	= 671.079E+03/jupiter_rad
		real,parameter :: eur_sid_per	= 3.551810*24.
		real,parameter :: eur_per		= 1. / ( 1./eur_sid_per - 1./jupiter_year )
		!	eur_rad note: R_Eur defined as 1560. in Galileo spacecraft data.
		!	Accepted value: 1565.
		real,parameter :: eur_rad		= 1560.	
		real,parameter :: eur_mass		= 479.7E+20
		real,parameter :: eur_incl		= 0.464*pi/180.
		real,parameter :: eur_init_rot	= ( atan((-2.379801E-03)/ &
			(-3.751687E-03)) + pi - jupiter_init_long ) / (2.0*pi) * &
			eur_per
		real,parameter :: eur_synodic	= 1. / ( 1./jupiter_per - 1./eur_sid_per )

	!	Moons: Ganymede
		!	gany_orbit_rad units: jupiter_rad
		real,parameter :: gany_orbit_rad	= 1070.0428E+03/jupiter_rad
		real,parameter :: gany_sid_per		= 7.154553*24.
		real,parameter :: gany_per			= 1. / ( 1./jupiter_per - 1./gany_sid_per )
		real,parameter :: gany_rad			= 2634.
		real,parameter :: gany_mass			= 1482.0E+20
		real,parameter :: gany_incl			= 0.186*pi/180.
		real,parameter :: gany_init_rot		= ( atan((-4.5815417E-03)/ &
			(-5.49035247E-03)) + pi - jupiter_init_long ) / (2.0*pi) * &
			gany_per

	!	Moons: Callisto
		!	call_orbit_rad units: jupiter_rad
		real,parameter :: call_orbit_rad	= 1883.0E+03/jupiter_rad
		real,parameter :: call_sid_per		= 16.689018*24.
		real,parameter :: call_per			= 1. / ( 1./jupiter_per - 1./call_sid_per )
		real,parameter :: call_rad			= 2403.
		real,parameter :: call_mass			= 1076.0E+20
		real,parameter :: call_incl			= 0.281*pi/180.
		real,parameter :: call_init_rot		= ( &
			atan(1.238159E-02/2.1730233E-03) - jupiter_init_long ) / &
			(2.0*pi) * call_per


	!	Saturn
	real,parameter :: saturn_orbit_rad	= 1.433436206E+09
	real,parameter :: saturn_year		= 10755.698*24.
	real,parameter :: saturn_rad		= 60268.
	real,parameter :: saturn_per		= 10.6562214
	real,parameter :: saturn_mass		= 568.319E+24
	real,parameter :: saturn_obliq		= 26.73*pi/180.
	real,parameter :: saturn_incl		= 6.0910615*pi/180.
	!	saturn_r_rot units: re_equiv, where corotation stops
	real,parameter :: saturn_r_rot		= 20.0
	real,parameter :: saturn_torus_infall = 0.05
	real,parameter :: saturn_tilt		= 0.0
	real,parameter :: saturn_init_long	= atan(6.5699885/6.40641043)
	real,parameter :: saturn_xdip		= vsmall
	real,parameter :: saturn_ydip		= zero
	real,parameter :: saturn_zdip		= zero

	!	Moons: Enceladus
		!	encel_orbit_rad units: saturn_rad
		real,parameter :: encel_orbit_rad	= 238.04E+03/saturn_rad
		real,parameter :: encel_sid_per		= 1.370218*24.
		real,parameter :: encel_per			= 1. / ( 1./encel_sid_per - 1./saturn_year )
		real,parameter :: encel_rad			= 252.3
		real,parameter :: encel_mass		= 10.805E+19
		real,parameter :: encel_incl		= 0.009*pi/180.
		real,parameter :: encel_init_rot	= ( &
			atan((-1.0650265E-03)/1.0809599E-03) + 2.0*pi - &
			saturn_init_long ) / (2.0*pi) * encel_per

	!	Moons: Titan
		!	titan_orbit_rad units: saturn_rad
		real,parameter :: titan_orbit_rad	= 1221.87E+03/saturn_rad
		real,parameter :: titan_sid_per		= 15.945421*24.
		real,parameter :: titan_per			= 1. / ( 1./titan_sid_per - 1./saturn_year )
		real,parameter :: titan_rad			= 2575.5
		real,parameter :: titan_mass		= 13455.3E+19
		real,parameter :: titan_incl		= 0.28
		real,parameter :: titan_init_rot	= ( atan(5.126197E-03/ &
			(-6.328986E-03)) + pi - saturn_init_long ) / (2.0*pi) * &
			titan_per
	

	!	Earth
	real,parameter :: earth_orbit_rad	= 1.49665015E+08
	real,parameter :: earth_year		= 8766.
	real,parameter :: earth_rad			= 6378.14
	real,parameter :: earth_per			= 24.
	real,parameter :: earth_mass		= 5.97219E+24
	real,parameter :: earth_incl		= 7.2515*pi/180.
	real,parameter :: earth_obliq		= 23.45*pi/180.
	!	earth_r_rot units: re_equiv
	real,parameter :: earth_r_rot		= 10.0
	real,parameter :: earth_torus_infall = 0.
	real,parameter :: earth_tilt		= (90.0 + 80.972)*pi/180.
	real,parameter :: earth_init_long	= atan(9.6724169E-01/ &
		(-1.771351E-01)) + pi
	!	earth_xdip units: planetary radii
	real,parameter :: earth_xdip		= 500./earth_rad
	real,parameter :: earth_ydip		= zero
	real,parameter :: earth_zdip		= zero


	!	The Moon
		!	luna_orbit_rad units: earth_rad
		real,parameter :: luna_orbit_rad	= 384400.0/earth_rad
		real,parameter :: luna_sid_per		= 27.321582*24.
		real,parameter :: luna_per			= 1. / ( 1./luna_sid_per - 1./earth_year )
		real,parameter :: luna_rad			= 1737.4
		real,parameter :: luna_mass			= 734.9E+20
		real,parameter :: luna_incl			= 5.145*pi/180.
		real,parameter :: luna_init_rot		= ( atan((-1.838126E-03)/ &
			(-1.949282E-03)) + pi - earth_init_long ) / (2.0*pi) * &
			luna_per

contains

	subroutine choose_system(bodyname,moonname)

		character*10, intent(in) :: bodyname, moonname

		real planet_orbit_rad, planet_year, planet_rad, &
			planet_per, planet_mass, planet_obliq, planet_incl, &
			r_lim, torus_infall, planet_tilt, planet_init_long, &
			planet_xdip, planet_ydip, planet_zdip, torus_dist, &
			moon_orbit_rad, moon_per, moon_rad, moon_mass, &
			moon_incl, moon_init_rot, moon_synodic

		common /planetary/planet_orbit_rad, planet_year, planet_rad, &
			planet_per, planet_mass, planet_obliq, planet_incl, &
			r_lim, torus_infall, planet_tilt, planet_init_long, &
			planet_xdip, planet_ydip, planet_zdip, torus_dist, &
			moon_orbit_rad, moon_per, moon_rad, moon_mass, moon_incl, &
			moon_init_rot, moon_synodic

		write(*,*) 'Choosing for system:'
		write(*,*) bodyname, ';', moonname


		select case (trim(bodyname))

			case ("jupiter")
				planet_orbit_rad = jupiter_orbit_rad
				planet_year = jupiter_year
				planet_rad = jupiter_rad
				planet_per = jupiter_per
				planet_mass = jupiter_mass
				planet_obliq = jupiter_obliq
				planet_incl = jupiter_incl
				r_lim = jupiter_r_rot
				torus_infall = jupiter_torus_infall
				torus_dist = io_orbit_rad
				planet_tilt = jupiter_tilt
				planet_init_long = jupiter_init_long
				planet_xdip = jupiter_xdip
				planet_ydip = jupiter_ydip
				planet_zdip = jupiter_zdip

			case ("saturn")
				planet_orbit_rad = saturn_orbit_rad
				planet_year = saturn_year
				planet_rad = saturn_rad
				planet_per = saturn_per
				planet_mass = saturn_mass
				planet_obliq = saturn_obliq
				planet_incl = saturn_incl
				r_lim = saturn_r_rot
				torus_infall = saturn_torus_infall
				torus_dist = encel_orbit_rad
				planet_tilt = saturn_tilt
				planet_init_long = saturn_init_long
				planet_xdip = saturn_xdip
				planet_ydip = saturn_ydip
				planet_zdip = saturn_zdip

			case ("earth")
				planet_orbit_rad = earth_orbit_rad
				planet_year = earth_year
				planet_rad = earth_rad
				planet_per = earth_per
				planet_mass = earth_mass
				planet_obliq = earth_obliq
				planet_incl = earth_incl
				r_lim = earth_r_rot
				torus_infall = earth_torus_infall
				torus_dist = moon_orbit_rad
				planet_tilt = earth_tilt
				planet_init_long = earth_init_long
				planet_xdip = earth_xdip	!	Incorrect values
				planet_ydip = earth_ydip
				planet_zdip = earth_zdip
				write(*,*)	'Warning: Dipole offset values not ', &
					'rigorously determined.'
				write(*,*)	'Warning: torus_dist location unknown.'

			case ("europa")
				planet_orbit_rad = saturn_orbit_rad
				planet_year = jupiter_year
				planet_rad = eur_rad
				planet_per = eur_per
				planet_mass = eur_mass
				planet_obliq = 0.0
				planet_incl = eur_incl
				r_lim = eur_rad
				torus_infall = 0.0
				torus_dist = 5*eur_rad	!	Arbitrary distance
				planet_tilt = vsmall
				planet_init_long = eur_init_rot
				planet_xdip = 0.0
				planet_ydip = 0.0
				planet_zdip = vsmall
				moon_synodic = eur_synodic

			case default
				write(*,*) 'Body name did not match valid options:'
				write(*,*) 'jupiter, saturn, earth, or europa'
				write(*,*) 'Defaulting to saturn.'
				planet_orbit_rad = saturn_orbit_rad
				planet_year = saturn_year
				planet_rad = saturn_rad
				planet_per = saturn_per
				planet_mass = saturn_mass
				planet_obliq = saturn_obliq
				planet_incl = saturn_incl
				r_lim = saturn_r_rot
				torus_infall = saturn_torus_infall
				torus_dist = encel_orbit_rad
				planet_tilt = saturn_tilt
				planet_init_long = saturn_init_long
				planet_xdip = saturn_xdip
				planet_ydip = saturn_ydip
				planet_zdip = saturn_zdip
				moon_synodic = eur_synodic

		end select


		select case (trim(moonname))

			case ("io")
				moon_orbit_rad = io_orbit_rad
				moon_per = io_per
				moon_rad = io_rad
				moon_mass = io_mass
				moon_incl = io_incl
				moon_init_rot = io_init_rot

			case ("europa")
				moon_orbit_rad = eur_orbit_rad
				moon_per = eur_per
				moon_rad = eur_rad
				moon_mass = eur_mass
				moon_incl = eur_incl
				moon_init_rot = eur_init_rot

			case ("ganymede")
				moon_orbit_rad = gany_orbit_rad
				moon_per = gany_per
				moon_rad = gany_rad
				moon_mass = gany_mass
				moon_incl = gany_incl
				moon_init_rot = gany_init_rot

			case ("callisto")
				moon_orbit_rad = call_orbit_rad
				moon_per = call_per
				moon_rad = call_rad
				moon_mass = call_mass
				moon_incl = call_incl
				moon_init_rot = call_init_rot

			case ("enceladus")
				moon_orbit_rad = encel_orbit_rad
				moon_per = encel_per
				moon_rad = encel_rad
				moon_mass = encel_mass
				moon_incl = encel_incl
				moon_init_rot = encel_init_rot

			case ("titan")
				moon_orbit_rad = titan_orbit_rad
				moon_per = titan_per
				moon_rad = titan_rad
				moon_mass = titan_mass
				moon_incl = titan_incl
				moon_init_rot = titan_init_rot

			case ("luna")
				moon_orbit_rad = luna_orbit_rad
				moon_per = luna_per
				moon_rad = luna_rad
				moon_mass = luna_mass
				moon_incl = luna_incl
				moon_init_rot = luna_init_rot

			case ("none")
				write(*,*) 'No moon selected for this body. ', &
					'Tracking a surface point.'
				moon_orbit_rad = planet_rad
				moon_per = planet_per
				moon_rad = vsmall
				moon_mass = vsmall
				moon_incl = planet_obliq
				moon_init_rot = planet_init_long

			case default
				write(*,*) 'Moon name did not match valid options:'
				write(*,*) 'io, europa, ganymede, callisto, ', &
					'enceladus, titan, or luna'
				write(*,*) 'Defaulting to enceladus.'
				moon_orbit_rad = encel_orbit_rad
				moon_per = encel_per
				moon_rad = encel_rad
				moon_mass = encel_mass
				moon_incl = encel_incl
				moon_init_rot = encel_init_rot
			
		end select

	end subroutine choose_system

end module astrometry
