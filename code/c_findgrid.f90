subroutine findgrid(xyz,n_grids,grid_minvals,grid_maxvals,num_pts,craft_gridpt)
    !
    !   Finds the nearest grid point and returns the xyz indices.
    !
	integer, intent(in) :: n_grids, num_pts(3)
    real, intent(in) :: xyz(3), grid_minvals(3,n_grids), grid_maxvals(3,n_grids)
	!
	integer, intent(out) :: craft_gridpt(4)	!	craft_gridpt contains: x_index,y_index,z_index,box_num
	integer axis
	real delta_pos, grid_spacing
	!	delta_pos is the distance above the minval for a given axis. 
	!	grid_spacing is the distance between adjacent grid points for a given axis.
	!
	craft_gridpt(4) = 0
	!	Find smallest grid xyz fits within
	do box=1, n_grids
		if( ( (xyz(1).gt.grid_minvals(1,box)) .and. (xyz(2).gt.grid_minvals(2,box)) .and. (xyz(3).gt.grid_minvals(3,box)) ) &
			.and. &
		( (xyz(1).lt.grid_maxvals(1,box)) .and. (xyz(2).lt.grid_maxvals(2,box)) .and. (xyz(3).lt.grid_maxvals(3,box)) ) ) then
			craft_gridpt(4) = box
			exit
		else
			cycle
		endif
	enddo
	!
	if(craft_grdpt(4).eq.0) then
		write(*,*) 'Craft outside of grid limits.'
		craft_grdpt(1) = 0
	else
		!
		do axis=1,3
			grid_spacing = ( grid_maxvals(axis,craft_grdpt(4)) - grid_minvals(axis,craft_grdpt(4)) ) / num_pts(axis)
			delta_pos = xyz(axis) - grid_minvals(axis,craft_gridpt(4))
			craft_gridpt(axis) = 1 + int(delta_pos/grid_spacing)
		enddo
		!
	endif
	!
    return
end subroutine findgrid
