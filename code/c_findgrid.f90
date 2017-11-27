subroutine findgrid(xyz,n_grids,grid_minvals,grid_maxvals,xspac,gridpt,re_equiv,cname)
    !
    !   Finds the nearest grid point and returns the xyz indices.
    !
	implicit none

	integer, intent(in) :: n_grids
    real, intent(in) :: xyz(3)
    real, intent(in) :: grid_minvals(3,n_grids)
    real, intent(in) :: grid_maxvals(3,n_grids)
	real, intent(in) :: xspac(n_grids)
	real, intent(in) :: re_equiv
	character*8, intent(in) :: cname
	!
	integer, intent(out) :: gridpt(4)	!	gridpt contains: x_index,y_index,z_index,box_num
    integer box
	integer axis
	real xyz_adj(3)	!	Spacecraft location in sim units
	real delta_pos
	real abit
	!	delta_pos is the distance above the minval for a given axis. 
	!	grid_spacing is the distance between adjacent grid points for a given axis.

	abit = 0.001
	gridpt(4) = 0
	xyz_adj(:) = xyz(:)/re_equiv

	!	Find smallest grid xyz fits within
	do box=1, n_grids
		if( ( (xyz_adj(1).gt.grid_minvals(1,box)) .and. (xyz_adj(2).gt.grid_minvals(2,box)) .and. (xyz_adj(3).gt.grid_minvals(3,box)) ) &
			.and. &
		( (xyz_adj(1).lt.grid_maxvals(1,box)) .and. (xyz_adj(2).lt.grid_maxvals(2,box)) .and. (xyz_adj(3).lt.grid_maxvals(3,box)) ) ) then
			gridpt(4) = box
			exit
		else
			cycle
		endif
	enddo
	
	if(gridpt(4).eq.0) then
		write(*,*) 'Craft outside of grid limits. Name, xyz: ', cname, xyz
		gridpt(4) = n_grids
		xyz_adj(:) = amax1( xyz_adj(:), (grid_minvals(:,n_grids)+abit) )
		xyz_adj(:) = amin1( xyz_adj(:), (grid_maxvals(:,n_grids)-abit) )
	endif
	
	do axis=1,3
		delta_pos = xyz_adj(axis) - grid_minvals(axis,gridpt(4))
		gridpt(axis) = 1 + int(delta_pos/xspac(gridpt(4)))
	enddo
	
    return
end subroutine findgrid
