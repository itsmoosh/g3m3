!
!	Calculates the number of hours into the current orbit, relative to the midnight line of the body.
!
subroutine find_rot_hrs(ut_jd, init_rot, solar_per, rot_hrs)

	integer, parameter	:: dp = selected_real_kind(17,300)

	real(dp), parameter	:: J2000 = 2451545.0
	
	real(dp), intent(in)	:: ut_jd	!	UT in JD format of desired instant
	real(dp), intent(in)	:: init_rot	!	rot_hrs at J2000.0 for the body
	real(dp), intent(in)	:: solar_per	!	Solar period of the body

	real(dp), intent(out)	:: rot_hrs	!	Number of hours after the midnight line of the body

	real(dp) n_rot

	n_rot = ( (ut_jd - J2000) * 24.0 + init_rot ) / solar_per	!	Decimal number of orbits/rotations after last midnight crossing before J2000
	rot_hrs = ( n_rot - floor(n_rot) ) * solar_per

	return
end subroutine find_rot_hrs
