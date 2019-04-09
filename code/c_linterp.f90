!
!   Linearly interpolates between the given spacetime coordinates
!	xyzt1 and xyzt2 to find the xyz coordinate for 
!	the given input time.
!	Be sure to match time units for xyzt1(4), xyzt2(4), and time.
!
subroutine linterp(xyzt1,xyzt2,time,xyz)

	implicit none

	real, intent(in) :: xyzt1(4)
	real, intent(in) :: xyzt2(4)
	real, intent(in) :: time
	real, intent(out) :: xyz(3)
	
	!	If we input the same start and end point, return the same xyz
	if(xyzt1(4) .eq. xyzt2(4))then
		xyz(:) = xyzt1(1:3)
	else
		!	x' = x_o + v_x * delta-t
		xyz(1) = xyzt1(1) + (xyzt2(1)-xyzt1(1))/(xyzt2(4)-xyzt1(4)) * &
			(time-xyzt1(4))
		xyz(2) = xyzt1(2) + (xyzt2(2)-xyzt1(2))/(xyzt2(4)-xyzt1(4)) * &
			(time-xyzt1(4))
		xyz(3) = xyzt1(3) + (xyzt2(3)-xyzt1(3))/(xyzt2(4)-xyzt1(4)) * &
			(time-xyzt1(4))
	endif
	
    return
end subroutine linterp
