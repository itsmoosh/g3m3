!
!	Initializes spacecraft data recording.
!	1. Default crafts, with fixed coordinates, are read in from spacecraft_info/defaults.pos
!	2. Trajectory craft are read in from .craft files in the same directory
!	3. Output data files are created for each spacecraft, if they do not exist already.
!
subroutine initcraft( craft_info, craft_data, craftnames, scin, scout, &
	dat_header, ncraft, ndef_craft, naux_craft, nheadlines, recording, &
	ntimes, vals, craft_gridpt, craftpos, craft_rot_ca, n_grids, &
	re_equiv, xspac, grid_minvals, grid_maxvals, git_hash_3d )

	implicit none

	character*32,	intent(in)	:: craft_info
	character*32,	intent(in)	:: craft_data
	character*8,	intent(in)	:: craftnames(ncraft)
	character*200,	intent(in)	:: dat_header
	
	integer, intent(in)		:: scin
	integer, intent(in)		:: scout
	integer, intent(in)		:: ncraft
	integer, intent(in)		:: ndef_craft
	integer, intent(in)		:: naux_craft
	integer, intent(in)		:: nheadlines
	integer, intent(in)		:: n_grids
	integer, intent(in)		:: vals

	real, intent(in)	:: re_equiv
	real, intent(in)	:: xspac(n_grids)
	real, intent(in)	:: grid_minvals(3,n_grids)
	real, intent(in)	:: grid_maxvals(3,n_grids)
	character*7, intent(in) :: git_hash_3d

	logical, intent(inout)		:: recording(ncraft)
	integer, intent(inout)		:: ntimes(ncraft,2)
	integer, intent(out)		:: craft_gridpt(4,ncraft)
	real, intent(out)			:: craftpos(4,ncraft,vals)
	real, intent(out)			:: craft_rot_ca(ncraft)

	logical	:: dat_exists
	logical :: dummy_exists
	integer	:: craftstat
	integer :: i_headnum	!	Dummy for skipping header lines
	integer :: m_recnum	!	Recording index
	integer :: m_stepcount
	integer :: n	!	Craft index
	integer :: cbox	!	Dummy for smallest box that fits craft point
	integer :: axis	!	For looping over spatial dimensions
	integer :: n_trajec	!	Number of full trajectories completed

	character*120 fname_scdat	!	fpath relative to main directory
	character*120 new_scdat
	character*120 junkline
	character*3 trajec_num

	! namelist 'crafthead'
	character*8 cname
	integer num_vals
	real rot_closest
	character*7 git_hash

	namelist/crafthead/cname,num_vals,rot_closest,git_hash

	! Initialization
	dat_exists = .false.
	dummy_exists = .false.
	craftstat = 0
	n_trajec = 0

	cname = 'default'
	num_vals = 0
	rot_closest = 0.0

	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	1. Read in default craft positions
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	open(scin+1,file=trim(craft_info)//'defaults.pos',status='unknown',form='formatted')
		!	Skip header line, only one line in default craft file
		read(scin+1,*) junkline
		
		do n=1,ndef_craft
			recording(n)=.true.
			read(scin+1,*) cname, craftpos(1,n,1), &
				craftpos(2,n,1), craftpos(3,n,1)
			craftpos(4,n,1) = 0
			!	Default craft locations in defaults.pos are in sim
			!	units, but craftpos needs to be in real units
			craftpos(1:3,n,1) = craftpos(1:3,n,1)*re_equiv
			!	We don't have a set number of recordings for default
			!	craft
			ntimes(n,2) = -1
			
			!	Default craft must be placed on grid points.
			!	If not, snap to nearest grid point.
			call findgrid(craftpos(1:3,n,1),n_grids,grid_minvals, &
				grid_maxvals,xspac,craft_gridpt(:,n),re_equiv,cname)
			cbox = craft_gridpt(4,n)
			do axis=1,3
				craftpos(axis,n,1) = grid_minvals(axis,cbox) &
					+ xspac(cbox)*re_equiv*(craft_gridpt(axis,n)-1)
			enddo
		enddo
	close(scin+1)

	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	2. Read in aux craft positions & times
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	do n=ndef_craft+1, ncraft
		recording=.true.
		craftstat=0
		open(scin+n,file=trim(craft_info)//trim(craftnames(n))// &
			'.craft',status='unknown',iostat=craftstat,form='formatted')
			read(scin+n,crafthead)
			read(scin+n,*) junkline		!	Skip column headers

			!	rot_closest is read in with crafthead
			craft_rot_ca(n) = rot_closest
			
			do m_recnum=1,ntimes(n,2)
				!	.craft files are formatted as t, x, y, z
				read(scin+n,*) craftpos(4,n,m_recnum), &
					craftpos(1,n,m_recnum), craftpos(2,n,m_recnum), &
					craftpos(3,n,m_recnum)
			enddo
			
			if (craftstat .ne. 0) then
				write(*,*) 'Error reading craft file: ', craftnames(n)
				ntimes(n,2) = 0
			elseif (ntimes(n,2).eq.1) then
				write(*,*) 'Craft file has only one xyzt coordinate:', &
					craftnames(n)
				ntimes(n,2) = 0
				ntimes(n,1) = 0
			endif
			
		close(scin+n)
	enddo

	call limcraft(craftpos,ntimes,vals,ndef_craft,ncraft,re_equiv, &
		n_grids,grid_minvals,grid_maxvals)

	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!	3. Initialize spacecraft output files
	!	@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    !	Open and go to end of .dat files for each craft	
	rot_closest = 0.0	!	For default craft, which are handled first
    do n=1,ncraft
		cname = craftnames(n)
		num_vals = ntimes(n,2)
		fname_scdat=trim(craft_data)//trim(cname)//'.dat'
		inquire(file=fname_scdat, exist=dat_exists)
		craftstat=0
		!	Open .dat files for each craft. Files will be closed
		!	automatically when process completes or dies.
    	if(.not.dat_exists) then
			open(scout+n, file=fname_scdat, iostat=craftstat, &
				status='unknown', form='formatted')
			!	git_hash in .dat files corresponds to the time of
			!	creation for the .dat file, not the .craft file as read
			!	in, in step 2 above.
			git_hash = git_hash_3d
			write(scout+n,crafthead)
			write(scout+n,*) trim(dat_header)
			ntimes(n,1) = 0
		else
			call system ('wc -l '//fname_scdat//' > '// &
				trim(craft_data)//'dat_working.txt')
			open(scout,file=trim(craft_data)//'dat_working.txt', &
				status='unknown',form='formatted')

				read(scout,*) ntimes(n,1), junkline
				!	Subtract header lines
				ntimes(n,1) = ntimes(n,1) - nheadlines
			close(scout)

			open(scout+n,file=fname_scdat,iostat=craftstat, &
				status='unknown',form='formatted')
			
			!	If true, this is an aux. craft which has finished
			!	recording its data
			if(ntimes(n,2).le.ntimes(n,1) .and. n.gt.ndef_craft) then
				call new_trajec( craft_data, dat_header, scout+n, &
					recording(n) )
				ntimes(n,1) = 0
			else
				!	Step through the file until the end, so we are ready
				!	to write new measurements. We have to skip past
				!	header lines, too.
				do m_stepcount=1, (ntimes(n,1) + nheadlines)
					read(scout+n,*) junkline
				enddo
			endif
				
			if(craftstat.ne.0) then
				write(*,*) 'Problem reading .dat file for craft: ', &
					craftnames(n)
				recording(n) = .false.
			endif
		endif
	enddo

	!	4. Cleanup
	inquire(file=trim(craft_data)//'dat_working.txt', exist=dat_exists)
	if(dat_exists) call system ('rm '//trim(craft_data)// &
		'dat_working.txt')

	return
end subroutine initcraft
