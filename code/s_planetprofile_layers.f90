subroutine planetprofile_layers( fname_planetprofile, planetprofile_f, &
	n_layers, planet_rad, depths, conductivities )
    !
    !   Calculates interior layer properties based on PlanetProfile
	!		radial profiles of conductivity.
	!
	!	Expects r (in km) to be the third column of the file
	!		corresponding to file key planetprofile_f and conductivity
	!		to be the last column (in S/m).
    !
	implicit none

	character*30, intent(in) :: planetprofile_fname
	integer, intent(in) :: planetprofile_f, n_layers
	real, intent(out) :: depths(n_layers+1), conductivities(n_layers)

	integer,parameter :: dummy = 595

	! Quantities as ordered in PlanetProfile output data
	real, allocatable :: P(:),T(:), r(:), rho(:),vp(:),vs(:), &
		qs_gamma(:),ks(:),gs(:), sigma(:)

	integer n_vals, n_line, n_layer, n_avgd
	integer i_start, i_stop, i_layer
	character*120 junkline
	logical ice_zero = .true.

	call system ('wc -l '//trim(planetprofile_fname)//' > dummy.txt')
	open(dummy_f,'dummy.txt')
		read(dummy_f,*) n_vals, junkline
	close(dummy_f)
	n_vals = n_vals - 1
	allocate( P(n_vals), T(n_vals), r(n_vals), rho(n_vals), &
		vp(n_vals), vs(n_vals), qs_gamma(n_vals), ks(n_vals), &
		gs(n_vals), sigma(n_vals) )

	open(planetprofile_f,file=trim(planetprofile_fname),status='unknown',form='formatted')
		read(planetprofile_f,*) junkline
		do n_line = 1, n_vals-1
			read(planetprofile_f,*) P(n_line), T(n_line), r(n_line), &
				rho(n_line), vp(n_line), vs(n_line), qs_gamma(n_line), &
				ks(n_line), gs(n_line), sigma(n_line)
			if( (sigma(n_line) .lt. 0.01) .and. ice_zero ) then
				! Reset depth of ice-ocean interface until we reach non-zero conductivity
				depths(1) = ( planet_rad - r(n_line) ) / planet_rad
				i_start = n_line
			else if( sigma(n_line-1) .lt. 0.01 .and. sigma(n_line) .ge. 0.01 ) then
				! We just passed the ice-ocean boundary and now have non-zero conductivity.
				! Stop overwriting the ice-ocean interface depth
				ice_zero = .false.
			else if( sigma(n_line-1) .ge. 0.01 .and. sigma(n_line) .lt. 0.01 ) then
				! Now we have just passed the ocean-mantle boundary. Record the depth.
				depths(n_layers+1) = ( planet_rad - r(n_line) ) / planet_rad
				i_stop - n_line
			endif
		enddo
	close(planetprofile_f)

	! Find number of radial values to be averaged for each layer
	n_avgd = (i_stop - i_start ) / n_layers
	! Initialize i_layer for the case that n_layers = 1
	i_layer = i_start

	do n_layer = 1, n_layers
		if(n_layer .lt. n_layers) then
			i_layer = i_start + (n_layer-1)*n_avgd
			depths(n_layer+1) = ( planet_rad - r(i_layer+n_avgd) ) / planet_rad
			conductivities(n_layer) = sum( sigma(i_layer:i_layer+n_avgd) ) / n_avgd
		else
			! Averages from i_layer in the past loop to the bottom of the ocean layer.
			! Accounts for n_avgd being floored as an integer.
			conductivities(n_layers) = sum( sigma(i_stop-i_layer:i_stop) ) / (i_stop - i_layer)			
		endif
	enddo
	
    return
end subroutine planetprofile_layers
