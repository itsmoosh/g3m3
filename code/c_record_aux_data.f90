!
!	Writes spacecraft measurements to output data file.
!
subroutine record_aux_data( i_craft, scout, num_inst, scfmt, ut, &
	sxyz, cube_vertices, meas_qty, ntimesn )
	
	implicit none

	integer, intent(in)		:: i_craft
	integer, intent(in)		:: scout
	integer, intent(in)		:: num_inst
	character*9, intent(in)	:: scfmt

	real, intent(in)	:: ut
	real, intent(in)	:: sxyz(3)
	real, intent(in)	:: cube_vertices(3,2)
	real, intent(in)	:: meas_qty(2,2,2,num_inst)

	integer, intent(inout)		:: ntimesn(2)

	integer nn		!	Loop counter
	real	scdata(num_inst)

	!	Separate trilin_interp calls needed for each scalar measurement
	!	and each component of vector quantities
	do nn=1,num_inst
		call trilin_interp( sxyz, cube_vertices, meas_qty(:,:,:,nn), &
			scdata(nn) )
	enddo

	!	Write all spacecraft instruments to .dat file
	write(scout+i_craft,scfmt) ut, sxyz, scdata

	ntimesn(1) = ntimesn(1) + 1

	return
end subroutine record_aux_data
