!
!	Handles spacecraft checks and data recording.
!	meas_qty(i) are scalar physical quantities to be measured by a
!	spacecraft
!	scdata(i) are the interpolated values corresponding to
!	meas_qty(i) that a spacecraft actually records.
!	
subroutine craft_handling( nx, ny, nz, n_grids, ncraft, n_vals, &
	ndef_craft, deflt_rec_skp, scout, t, ut, utold, tgraph, &
	t_equiv, t_step_max_def, xmoon, ymoon, zmoon, sxyz, &
	re_equiv, moon_per, repeat_flybys, craft_datadir, &
	grid_minvals, grid_maxvals, xspac, craftnames, scfmt, &
	dat_header, num_inst, spam_limit, craftpos, bx,by,bz, &
	nskipped, recording, ntimes, cgridpt, t_step_max, &
	flyby_ref_time, xcraft, zcraft )

	implicit none

	integer, parameter :: dp = kind(1.d0)
	!	Set sc_delt to less than the interval between spacecraft
	!	measurements. t_step is capped at this amount so as to
	!	hit each spacecraft measurement. This number is in seconds.
	real, parameter :: sc_delt = 0.30

	integer, intent(in)		:: nx, ny, nz, n_grids
	integer, intent(in)		:: ncraft
	integer, intent(in)		:: n_vals
	integer, intent(in)		:: ndef_craft
	integer, intent(in)		:: deflt_rec_skp
	integer, intent(in)		:: scout
	real(dp), intent(in)	:: t, ut, utold
	real, intent(in)		:: tgraph, t_equiv
	real, intent(in)		:: t_step_max_def
	real, intent(in)		:: xmoon, ymoon, zmoon
	real, intent(in)		:: sxyz(3)
	real, intent(in)		:: re_equiv, moon_per

	logical, intent(in)		:: repeat_flybys
	character*32, intent(in)	:: craft_datadir

	real, intent(in), dimension(3,n_grids)	:: grid_minvals, grid_maxvals
	real, intent(in)	:: xspac(n_grids)

	character*8, intent(in)		:: craftnames(ncraft)
	character*9, intent(in)		:: scfmt
	character*200, intent(in)	:: dat_header
	
	integer, intent(in) :: num_inst
	integer, intent(in) :: spam_limit
	real, intent(in)	:: craftpos(4,ncraft,n_vals)


	real, intent(in), dimension(nx,ny,nz,n_grids)	:: bx,by,bz

	integer, intent(inout)	:: nskipped
	logical, intent(inout)	:: recording(ncraft)
	integer, intent(inout)	:: ntimes(ncraft,2)
	integer, intent(out)	:: cgridpt(4,ncraft)
	real, intent(inout)		:: t_step_max
	real, intent(inout), dimension(ncraft)	:: flyby_ref_time
	real, intent(inout), dimension(4,ncraft)	:: xcraft, zcraft

	!	***************
	!	Dummy variables
	!	***************

	!	For handling skipped measurement times:
	logical :: sc_stepping = .false.
	logical :: sc_record = .false.
	logical :: sc_fast_forward = .false.
	integer :: sc_ff_count = 0

	real cube_vertices(3,2)
	real scdata(num_inst)
	real meas_qty(2,2,2,num_inst)
	real sctime
	
	integer spam_reduct
	integer n

	spam_reduct = spam_limit

	!	---------------------
	!		DEFAULT CRAFT
	!	---------------------

	!	Only record values for default craft every
	!	deflt_rec_skp number of loops
	if(nskipped.ge.deflt_rec_skp .or. t.ge.tgraph) then
		write(*,*) 'Recording default spacecraft measurements.'

		!	Update position of moon diagnostic craft
		if(trim(craftnames(2)) == 'moondgn') then
			xcraft(1,2)=xmoon*re_equiv*1.05
			xcraft(2,2)=ymoon*re_equiv*1.05
			xcraft(3,2)=zmoon*re_equiv*1.05
			xcraft(4,2)=ut
			call findgrid( xcraft(1:3,2),n_grids,grid_minvals, &
				grid_maxvals,xspac,cgridpt(:,2),re_equiv, &
				craftnames(2) )
		endif
	
		do n=1,ndef_craft
			if(recording(n)) then
				!	WARNING: A line must be added below each time
				!	num_inst is incremented, to add the new
				!	quantity to the array scdata.
				sctime = ut
				scdata(1) = bx( cgridpt(1,n),cgridpt(2,n), &
					cgridpt(3,n),cgridpt(4,n) )
				scdata(2) = by( cgridpt(1,n),cgridpt(2,n), &
					cgridpt(3,n),cgridpt(4,n) )
				scdata(3) = bz( cgridpt(1,n),cgridpt(2,n), &
					cgridpt(3,n),cgridpt(4,n) )
				!	Default craft xyz are in simulation
				!	coordinates, including measurements.
				write(scout+n, scfmt) sctime, xcraft(1:3,n), &
					scdata
			endif
		enddo
	
		ntimes(n,1) = ntimes(n,1) + 1
		nskipped = 0
	else
		nskipped = nskipped + 1
	endif

	!	-----------------
	!		AUX CRAFT
	!	-----------------

	do n=ndef_craft+1,ncraft
		if( recording(n) ) then
			if( utold.le.xcraft(4,n) &
				.and. ut.gt.xcraft(4,n) ) then
				sc_record = .true.
			else if( utold.gt.xcraft(4,n) ) then
				write(*,*) 'WARNING: t_step too large, ', &
					'fast-forwarding for craft: ', &
					craftnames(n), ' at ut = ', ut, &
					' with xcraft(t) = ', xcraft(4,n)
				sc_record = .true.
				sc_fast_forward = .true.
				sc_ff_count = 0
			endif

			do while(sc_record)
				if( (spam_reduct .ge. spam_limit) &
					.or. (ntimes(n,1).eq.ntimes(n,2)-1) ) then
					write(*,*) '-'
					write(*,*) 'Recording data for craft: ', &
						craftnames(n)
					write(*,*) '-'
					write(*,*) ntimes(n,1)+1, ' of ', &
						ntimes(n,2), 'ut = ', ut, &
						'flyby_ref_time = ', &
						flyby_ref_time(n), 'xcraft(t) = ', &
						xcraft(4,n), 'craftpos(t) = ', &
						craftpos(4,n,ntimes(n,1)+1)
					spam_reduct = 1
				else
					spam_reduct = spam_reduct + 1
				endif

				sctime = ut
				call find_aux_cube( craftnames(n), n, scout, &
					n_grids, grid_minvals, grid_maxvals, &
					xspac, re_equiv, xcraft(:,n), zcraft(:,n), &
					sctime, ntimes(n,:), recording(n), &
					cgridpt(:,n), sxyz, cube_vertices, scfmt )

				!	WARNING: Each time num_inst is incremented,
				!	a line must be added below to add the new
				!	quantity to the array meas_qty. The lines
				!	are identical except the final index of
				!	meas_qty and the name of the array we are
				!	pulling a measurement from.
				meas_qty(:,:,:,1) = bx( &
					cgridpt(1,n):cgridpt(1,n)+1, &
					cgridpt(2,n):cgridpt(2,n)+1, &
					cgridpt(3,n):cgridpt(3,n)+1, cgridpt(4,n) )
				meas_qty(:,:,:,2) = by( &
					cgridpt(1,n):cgridpt(1,n)+1, &
					cgridpt(2,n):cgridpt(2,n)+1, &
					cgridpt(3,n):cgridpt(3,n)+1, cgridpt(4,n) )
				meas_qty(:,:,:,3) = bz( &
					cgridpt(1,n):cgridpt(1,n)+1, &
					cgridpt(2,n):cgridpt(2,n)+1, &
					cgridpt(3,n):cgridpt(3,n)+1, cgridpt(4,n) )

				!	Writes data to disk and increments
				!	ntimes(n,1)
				call record_aux_data( n, scout, num_inst, &
					scfmt, sctime, sxyz, cube_vertices, &
					meas_qty, ntimes(n,:) )
			
				!	Check if spacecraft has finished its whole
				!	trajectory
				if(ntimes(n,1).ge.ntimes(n,2)) then
					if(repeat_flybys) then
						flyby_ref_time(n) = flyby_ref_time(n) &
							+ moon_per
						call new_trajec( craft_datadir, &
							dat_header, scout+n, recording(n) )
						ntimes(n,1) = 0
						xcraft(:,n) = craftpos(:,n,1)
						xcraft(4,n) = xcraft(4,n) &
							+ flyby_ref_time(n)
						zcraft(:,n) = craftpos(:,n,2)
						zcraft(4,n) = zcraft(4,n) &
							+ flyby_ref_time(n)
					else
						recording(n) = .false.
						!	Arbitrarily high number so it's
						!	always greater than ut
						xcraft(4,n) = 9.9e30
						write(*,*) 'Finished recording ', &
							'trajectory for craft: ', &
							craftnames(n)
					endif

					!	Speed time stepping back up if we're not
					!	near another flyby
					if( ut+(ut-utold)*2 .le. &
						minval(xcraft(4,ndef_craft+1:ncraft)) &
						) then

						sc_stepping = .false.
						write(*,*) 'Trajectory finished, ', &
							'speeding back up.'
					endif
				else if( ntimes(n,1).ge.(ntimes(n,2)-1) ) then
					xcraft(:,n) = zcraft(:,n)
				else
					xcraft(:,n) = zcraft(:,n)
					zcraft(:,n) = craftpos(:,n,ntimes(n,1)+1)
					zcraft(4,n) = zcraft(4,n) &
						+ flyby_ref_time(n)
				endif	!end if(finished trajectory?)

				sc_ff_count = sc_ff_count + 1
				if(.not. sc_fast_forward) then
					sc_record = .false.
				else if(xcraft(4,n) .ge. ut) then
					sc_record = .false.
					write(*,*) '-'
					write(*,*) 'Done fast-forwarding, ', &
						'recorded ', sc_ff_count, &
						' points at once.'
					write(*,*) '-'
					sc_ff_count = 0
					sc_fast_forward = .false.
				endif
			enddo	!end do while(sc_record)
		endif	!end if(recording)
	enddo

	if(.not. sc_stepping) then
		!	When we get close to spacecraft recording times,
		!	slow down the time stepping so we don't skip points.
		if( ut+(ut-utold)*2 .ge. &
			minval(xcraft(4,ndef_craft+1:ncraft)) ) then
			write(*,*) 'Slowing down at ut = ', ut, &
				' for measurement at ut = ', &
				minval(xcraft(4,ndef_craft+1:ncraft)), &
				' with ut-utold = ', ut-utold
			sc_stepping = .true.
			t_step_max = sc_delt/t_equiv
		else
			t_step_max = t_step_max_def
		endif
	endif

	return
end subroutine craft_handling
