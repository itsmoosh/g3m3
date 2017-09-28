subroutine utc_to_jd(ut_string, jd)
	!	Converts a UTC string to JD number.
	!	Expects a string with format YYYY-MM-DDTHH:MM:SS.SSS
	!	Outputs a real number that requires at least double precision

	integer, parameter	:: dp = selected_real_kind(17,300)
	real(dp), parameter	:: J2000 = 2451545.0

	character*23, intent(in)	:: ut_string
	real(dp), intent(out)		:: jd

	integer :: year
	integer :: month
	integer :: day
	integer :: hours
	integer :: mins
	real	:: seconds

	integer	:: leapdays
	integer	:: monthdays
	logical	:: leapyear

	!	Fill date/time parameters from UT string
	read(ut_string(1:4),  *) year
	read(ut_string(6:7),  *) month
	read(ut_string(9:10), *) day
	read(ut_string(12:13),*) hours
	read(ut_string(15:16),*) mins
	read(ut_string(18:23),*) seconds

	!	Adds or subtracts extra days for each leap year between current year and J2000.0
	leapdays = int((year-2000.0)/4.0)
	!	Accounts for leap days during leap years
	if ( (mod(year,4).eq.0) .and. (month.gt.2) ) leapdays = leapdays + 1
	!	Accounts for the 2000 leap day on years after J2000.0
	if ( (mod(year,4).ne.0) .and. (year.gt.2000) ) leapdays = leapdays + 1

	!	Add days passed before current month started
	select case(month)
		case (1)
			monthdays = 0
		case (2)
			monthdays = 31
		case (3)
			monthdays = 31+28
		case (4)
			monthdays = 31+28+31
		case (5)
			monthdays = 31+28+31+30
		case (6)
			monthdays = 31+28+31+30+31
		case (7)
			monthdays = 31+28+31+30+31+30
		case (8)
			monthdays = 31+28+31+30+31+30+31
		case (9)
			monthdays = 31+28+31+30+31+30+31+31
		case (10)
			monthdays = 31+28+31+30+31+30+31+31+30
		case (11)
			monthdays = 31+28+31+30+31+30+31+31+30+31
		case (12)
			monthdays = 31+28+31+30+31+30+31+31+30+31+30
	end select
	
	jd = (year-2000.0)*365.0 + leapdays + monthdays + day + hours/24.0 + mins/24.0/60.0 + seconds/24.0/3600.0
	jd = jd - 1.5 + J2000

	return
end subroutine utc_to_jd
