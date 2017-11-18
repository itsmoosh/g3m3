!
!	********************************************************************
!	This program reads E-Phi-O position and time data for Galileo flybys
!	and prints a .craft file with xyz values aligned to simulation
!	grid, relative to Jupiter's center, and UT values in hours relative
!	to the time of closest approach.
!	********************************************************************
!
!	Requires:
!		m_astrometry.f90
!		c_utc_to_jd.f90
!
program make_craft

	use astrometry
	implicit none

!	Declarations

	integer, parameter :: dp = selected_real_kind(17,300)
	integer, parameter :: dummy_f = 10
	character*11, parameter	:: dummy_craft = 'dummy.craft'
	character,parameter		:: tab=char(9)
	character*120,parameter	:: dat_header='UT_rel(hrs)'//tab//'x(R_E)'//tab//'y(R_E)'//tab//'z(R_E)'!//tab//'Bxval'//tab//'Byval'//tab//'Bzval'

	!	J2000_adj is a reference UTC in JD that is offset from J2000.0 by a whole
	!	number of Europa synodic periods, and occurs before the first Galileo flyby.
	!	Using this reference causes all flybys to have a positive number of past_orbits,
	!	and allows us to use the same value for eur_init_rot as J2000.0. This is
	!	a difference of 350.0 orbits relative to the Jupiter noon meridian.
	real(dp), parameter	:: J2000_adj	= 2450300.84745

	integer :: past_orbits
	integer :: orbit_num
	integer :: flybystat
	integer :: line_num

	character*10 :: bodyname = "jupiter"
	character*10 :: moonname = "europa"
 
	!	Except for E25, these times are found in the .lbl file for each flyby.
	character*23, parameter :: utc_closest(26) = [ 'NA', 'NA', 'NA', &
		'1996-12-19T06:52:58.000', &	!	E4
		'NA', 'NA', 'NA', 'NA', 'NA', 'NA', &
		'1997-11-06T20:31:44.000', &	!	E11
		'1997-12-16T12:03:20.000', &	!	E12
		'NA', &
		'1998-03-29T13:21:06.000', &	!	E14
		'1998-05-31T21:12:58.000', &	!	E15
		'NA', &
		'1998-09-26T03:54:20.000', &	!	E17
		'NA', &
		'1999-02-01T02:19:50.000', &	!	E19
		'NA', 'NA', 'NA', 'NA', 'NA', &
		'1999-11-25T16:29:00.000', &	!	E25
		'2000-01-03T17:59:44.000' ]		!	E26

	!	jd_closest are UTC values in days in the J2000.0 epoch.
	real(dp) :: jd_closest(26)

	character*23 :: ut_string
	real :: x_eur
	real :: y_eur
	real :: z_eur
	real :: bx
	real :: by
	real :: bz
	real :: bmag

	real(dp) :: ut_jd
	real :: ut_rel
	real :: x_sim
	real :: y_sim
	real :: z_sim
	real :: bx_sim
	real :: by_sim
	real :: bz_sim
	
	character*150	:: junkline
	character*42	:: filepath
	character*30	:: craftpath
	character*13, parameter	:: datafolder = 'data/flybys/e'
	character*14, parameter	:: extension  = '-mag-ephio.tab'
	character*22, parameter	:: craftfolder=	'spacecraft_info/gali-e'
	character*6, parameter	:: craftxtn	  =	'.craft'
	character*2		:: orbit_num_char

	! namelist 'crafthead'
	character*8	 :: cname = 'default'
	integer	:: num_vals = 0
	real	:: rot_closest = 0.0
	character*7 :: git_hash = 'pl_hold'

	namelist/crafthead/cname,num_vals,rot_closest,git_hash

!	Executions

	do orbit_num=4,26
		if (len(trim(utc_closest(orbit_num))) .lt. 20) cycle

		write(orbit_num_char,'(I2)') orbit_num
		write(*,*)         '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
		write(*,'(1X,A,A,A)') 'Writing craft file for orbit E',trim(adjustl(orbit_num_char)),'.'
		write(*,*)         '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'

		write(cname,'(A,I2.2)') 'gali-e',orbit_num

		write(filepath,*) orbit_num
		craftpath = craftfolder//trim(adjustl(orbit_num_char))//craftxtn
		filepath = datafolder//trim(adjustl(orbit_num_char))//extension

		call system ('git rev-parse --short HEAD > '//trim(dummy_craft))
		call system ('wc -l '//trim(filepath)//' >> '//trim(dummy_craft))

		open(dummy_f,file=trim(dummy_craft),status='unknown',form='formatted')
			read(dummy_f,*) git_hash
			read(dummy_f,*) num_vals, junkline
		close(dummy_f)
		call system ('rm '//trim(dummy_craft))

		call utc_to_jd( utc_closest(orbit_num), jd_closest(orbit_num) )

		!	Calculate where Europa is located at the time of closest approach
		past_orbits = int( (jd_closest(orbit_num)-J2000_adj)/eur_per + eur_init_rot/eur_per )
		rot_closest = ( (jd_closest(orbit_num)-J2000_adj)/eur_per + eur_init_rot/eur_per - past_orbits )*eur_per
		if(rot_closest.ge.eur_per) then
			write(*,*) "Something's wrong, rot_closest too big: ", rot_closest, ' Aborted.'
			cycle
		else if(rot_closest.lt.0) then
			write(*,*) "Something's wrong, rot_closest is negative: ", rot_closest, ' Aborted.'
			cycle
		endif

		open(1,file=trim(filepath),iostat=flybystat,form='formatted')
		open(2,file=trim(craftpath),form='formatted')
			write(2,crafthead)
			write(2,*) trim(dat_header)
			write(*,*) 'Header lines:'
			write(*,crafthead)
			write(*,*) trim(dat_header)

			do line_num = 1, num_vals

				read(1,*) ut_string,bx,by,bz,bmag,x_eur,y_eur,z_eur

				!	Step 1: Convert UT values from YYYY-MM-DDTHH:MM:SS.SSS format to JD.
				call utc_to_jd(ut_string, ut_jd)
				ut_rel = ut_jd - jd_closest(orbit_num)
				ut_rel = ut_rel * 24.0	!	UT of the measurement, in hours relative to time of closest approach

				!	Step 2: Convert spacecraft locations into simulation coordinates.
				!		Spacecraft locations in .tab file are in E-Phi-O units, relative to the body center of Europa at the time of closest approach.
				!		E-Phi-O is x along the corotation direction of plasma (phi around Jupiter's spin axis), z is parallel to Jovian spin axis, and y completes the right-handed set (toward Jupiter).
				!		Sim unit distances are in units of re_equiv, which has units of planetary radius.
				x_sim = x_eur
				y_sim = y_eur
				z_sim = z_eur

				!	Step 3: Convert B values into simulation coordinates. No change is needed when x is directed along the corotating plasma direction.
				bx_sim = bx
				by_sim = by
				bz_sim = bz				

				!	Step 4: Write everything to .craft file.
				write(2,*) ut_rel, x_sim, y_sim, z_sim!, bx_sim, by_sim, bz_sim, bmag
			enddo
		close(2)
		close(1)
		write(*,*) '-'
		write(*,*) 'Finished craft file for orbit ', orbit_num_char
		write(*,*) '-'
	enddo	!	End loop over orbit number

end program make_craft
