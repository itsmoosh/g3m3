!
!	Calculates transpose(banded_matrix).x_vec
!	Warning: symmetric matrix assumed here,
!	taking advantage of the banded structure of banded_matrix
!
subroutine banded_tpose(nbands,npts,divB_band_pos,banded_matrix,x_vec,tp_result)

	implicit none

	!	*********
	!	Arguments
	!	*********

	integer, intent(in)	:: nbands, npts
	integer, intent(in)	:: divB_band_pos(nbands)
	real, intent(in)	:: banded_matrix(nbands,npts), x_vec(npts)
	real, intent(out)	:: tp_result(npts)
	integer band, i, ii

	!	Initialize product
	do i = 1, npts
		tp_result(i) = 0.
		
		do band = 1, nbands
			ii = i + divB_band_pos(band)
			if( (ii.ge.1) .and. (ii.le.npts) ) then
				tp_result(i) = tp_result(i) + banded_matrix(band,ii)*x_vec(ii)
			endif
		enddo
	enddo

	return
end subroutine banded_tpose
