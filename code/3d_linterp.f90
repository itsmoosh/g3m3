subroutine linterp(xyzt1,xyzt2,time,xyz)
    !
    !   Linearly interpolates between the given spacetime coordinates
	!		xyzt1 and xyzt2 to find the xyz coordinate for 
	!		the given input time.
    !
	real, intent(in) :: xyzt1(4), xyzt2(4), time
	real, intent(out) :: xyz(3)
	!
	!	x' = x_o + v_x * delta-t
	xyz(1) = xyzt1(1) + (xyzt2(1)-xyzt1(1))/(xyzt2(4)-xyzt1(4)) * (time-xyzt1(4))
	xyz(2) = xyzt1(2) + (xyzt2(2)-xyzt1(2))/(xyzt2(4)-xyzt1(4)) * (time-xyzt1(4))
	xyz(3) = xyzt1(3) + (xyzt2(3)-xyzt1(3))/(xyzt2(4)-xyzt1(4)) * (time-xyzt1(4))
	!
    return
end subroutine linterp
