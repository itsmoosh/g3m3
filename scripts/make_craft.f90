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

	integer :: orbit_num
	integer :: flybystat

	character*10 :: bodyname = "jupiter"
	character*10 :: moonname = "europa"
 
	!	t_closest are UTC values in days in the J2000.0 epoch.
	!	Except for E25, these values are found in the .lbl file for each flyby.
	real(dp), parameter :: t_closest(26) = [ 0., 0., 0., &
		2450436.786782407, &	!	E4
		0., 0., 0., 0., 0., 0., &
		2450759.355370370, &	!	E11
		2450799.002314815, &	!	E12
		0., &
		2450902.056319444, &	!	E14
		2450965.384004630, &	!	E15
		0., &
		2451082.662731481, &	!	E17
		0., &
		2451210.597106481, &	!	E19
		0., 0., 0., 0., 0., &
		2451508.186805556, &	!	E25
		2451547.249814815 ]		!	E26

	!	xyz values of Europa at the time of closest approach
	real(dp), parameter :: x_closest(26)
	real(dp), parameter :: y_closest(26)
	real(dp), parameter :: z_closest(26)

	character*23 :: ut_string
	real :: x_eur
	real :: y_eur
	real :: z_eur
	real :: bx
	real :: by
	real :: bz
	real :: bmag

	real(dp) :: ut_jd
	real(dp) :: ut_rel
	real :: x_sim
	real :: y_sim
	real :: z_sim
	real :: bx_sim
	real :: by_sim
	real :: bz_sim
	real :: rot_hrs
	
	character*150	:: junkline
	character*42	:: filepath
	character*20	:: craftpath
	character*13, parameter	:: datafolder = 'data/flybys/e'
	character*14, parameter	:: extension  = '-mag-ephio.tab'
	character*22, parameter	:: craftfolder=	'spacecraft_info/gali-e'
	character*6, parameter	:: craftxtn	  =	'.craft'

!	Executions

	do orbit_num=4,26
		if (t_closest(orbit_num) .lt. 1.0) cycle

		write(filepath,*) orbit_num
		craftpath = craftfolder//trim(adjustl(filepath))//craftxtn
		filepath = datafolder//trim(adjustl(filepath))//extension

		open(1,file=trim(filepath),iostat=flybystat,form='formatted')

		!open(2,file=trim(craftpath),form='formatted')
			do while(flybystat.eq.0)

				read(1,*) ut_string,bx,by,bz,bmag,x_eur,y_eur,z_eur

				!	Step 1: Convert UT values from YYYY-MM-DDTHH:MM:SS.SSS format to JD.
				call utc_to_jd(ut_string, ut_jd)
				ut_rel = ut_jd - t_closest(orbit_num)

				ut_rel * 24.0	!	UT in hours relative to time of closest approach
				write(*,*) rot_hrs, ut_rel
				exit

				!	Step 2: Convert spacecraft locations into simulation coordinates.
				!		Spacecraft locations in .tab file are in E-Phi-O units, relative to the body center of Europa at the time of closest approach.
				!		E-Phi-O is x along the corotation direction of plasma (phi around Jupiter's spin axis), z is parallel to Jovian spin axis, and y completes the right-handed set.
				

				!	Step 3: Convert B values into simulation coordinates.
				bx_sim = bx
				by_sim = by
				bz_sim = bz				

				!	Step 4: Write everything to .craft file.
				write(*,*) ut_adj, bx_sim, by_sim, bz_sim, bmag, x_sim, y_sim, z_sim
			enddo
		close(1)
		close(2)
		exit
	enddo	!	End loop over orbit number

end program make_craft
