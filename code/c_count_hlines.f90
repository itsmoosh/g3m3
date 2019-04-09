!
!	Prints craft output headers to a file and checks how many lines were printed.
!
subroutine count_hlines(dat_header, dummy_f, nheadlines)
	
	implicit none

	character*200, intent(in)	:: dat_header
	integer, intent(in)		:: dummy_f
	integer, intent(inout)	:: nheadlines

	character*11 dummy_craft
	character*120 junkline
	logical dummy_exists

	character*7 git_hash
	character*8 cname
	integer num_vals
	real rot_closest

	namelist/crafthead/cname,num_vals,rot_closest,git_hash

	dummy_craft='dummy.craft'
	git_hash = 'pl_hold'
	cname = 'default'
	num_vals = 0
	rot_closest = 0.0

	inquire(file=trim(dummy_craft), exist=dummy_exists)
	if(dummy_exists) call system ('rm '//trim(dummy_craft))
	
	open(dummy_f,file=trim(dummy_craft),status='unknown',form='formatted')
		write(dummy_f,crafthead)
		write(dummy_f,*) trim(dat_header)
		call system ('wc -l '//trim(dummy_craft)//' >> '//trim(dummy_craft))
		read(dummy_f,*) nheadlines, junkline
	close(dummy_f)
	call system ('rm '//trim(dummy_craft))

	return
end subroutine count_hlines
