!
!	Renames spacecraft output file based on how many
!	full trajectories have been completed.
!
subroutine new_trajec( craft_data, dat_header, outf_num, recording )

	implicit none
	
	character*32, intent(in)	:: craft_data
	character*200, intent(in)	:: dat_header
	integer, intent(in)			:: outf_num
	logical, intent(inout)		:: recording

	character*120 fname_scdat	!	File path relative to main directory
	character*120 new_scdat
	character*3 trajec_num
	logical datf_exists
	integer craftstat
	integer n_trajec

	! namelist 'crafthead'
	character*8 cname
	integer num_vals
	real rot_closest
	character*7 git_hash

	namelist/crafthead/cname,num_vals,rot_closest,git_hash

	n_trajec = 0
	!	We only enter this subroutine if this begins true.
	datf_exists = .true.

	rewind(outf_num)
	read(outf_num,crafthead)
	close(outf_num)

	fname_scdat=trim(craft_data)//trim(cname)//'.dat'
	n_trajec = 0
	do while (datf_exists)
		n_trajec = n_trajec + 1
		write(trajec_num,'(I2.2)') n_trajec
		!	Use format 'craftnm_06.dat' for 6th trajectory of craftnm
		trajec_num = '_'//trim(trajec_num)

		new_scdat = trim(craft_data)//trim(cname)//trajec_num//'.dat'
		inquire(file=new_scdat, exist=datf_exists)
	enddo
	call system('mv '//trim(fname_scdat)//' '//trim(new_scdat))
	!
	open(outf_num,file=fname_scdat,iostat=craftstat, &
		status='unknown',form='formatted')
	write(outf_num,crafthead)
	write(outf_num,*) trim(dat_header)

	if(craftstat.ne.0) then
		write(*,*) 'Problem creating new file for craft: ', cname
		recording = .false.
	else
		write(*,*) '------------------------------------'
		write(*,*) 'New file created for craft: ', cname
		write(*,*) '------------------------------------'
	endif

	return
end subroutine new_trajec
