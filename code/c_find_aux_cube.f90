!
!	Writes spacecraft measurements to output data file.
!
subroutine find_aux_cube( cname, i_craft, scout, n_grids, &
	grid_minvals, grid_maxvals, xspac, re_equiv, xcraftn, zcraftn, ut, &
	ntimesn, recordingn, gridpt, sxyz, cube_vertices )
	
	implicit none

	character*8, intent(in)	:: cname

	integer, intent(in)		:: i_craft
	integer, intent(in)		:: scout
	integer, intent(in)		:: n_grids

	real, intent(in)	:: grid_minvals(3,n_grids)
	real, intent(in)	:: grid_maxvals(3,n_grids)
	real, intent(in)	:: xspac(n_grids)
	real, intent(in)	:: re_equiv
	real, intent(in)	:: xcraftn(4)
	real, intent(in)	:: zcraftn(4)
	real, intent(in)	:: ut

	integer, intent(inout)	:: ntimesn(2)
	logical, intent(inout)	:: recordingn

	!	Grid index values for point just below sxyz in each dimension,
	!	and number of smallest box that fits the craft
	integer, intent(out)	:: gridpt(4)
	!	Interpolated location of spacecraft in simulation coordinates
	real, intent(out)		:: sxyz(3)
	real, intent(out)		:: cube_vertices(3,2)

	integer	cbox	!	Dummy for smallest box that fits craft point
	integer	axis	!	For looping over spatial dimensions

	if( (xcraftn(4).eq.zcraftn(4)) .or. (ut.ge.zcraftn(4)) ) then
		sxyz(:) = xcraftn(1:3)
	else
		call linterp(xcraftn(:),zcraftn(:),ut,sxyz)
	endif
	!	Find closest xyz values below
	call findgrid(sxyz,n_grids,grid_minvals,grid_maxvals,xspac,gridpt,re_equiv,cname)
	cbox = gridpt(4)
	!	gridpt now contains x,y,z,box indices: the indices of the
	!	closest xyz LESS THAN the craft location, and the smallest box
	!	the craft fits within.

	!	x index is set to zero if there is a problem.
	if(gridpt(1) .le. 0) then
		recordingn = .false.
		close(scout+i_craft)
		write(*,*) 'Gridding problem with craft ', cname, &
		'at UT = ', ut, ' Future points skipped.'
		return
	endif

	!	Indices in gridpt are used to identify the physical parameters
	!	at nearby points
	!	Position values in cube_vertices are used to interpolate
	!	measurement values
	do axis=1,3
		cube_vertices(axis,1) = ( grid_minvals(axis,cbox) + (gridpt(axis)-1.)*xspac(cbox) )*re_equiv
		cube_vertices(axis,2) = cube_vertices(axis,1) + xspac(cbox)*re_equiv
	enddo
		
	return
end subroutine find_aux_cube
