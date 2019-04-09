!
!   Counts the number of spacecraft listed in defaults.pos and
!	those with their own .craft trajectory files. 
!	Returns the number of default craft, specified craft trajectories,
!	and the total number of spacecraft.
!
subroutine countcraft(craft_info,scin,ncraft,ndef_craft,naux_craft)

	implicit none

	character*32, intent(in) :: craft_info
	integer, intent(in) :: scin
	integer, intent(out) :: ncraft, ndef_craft, naux_craft
	character*80 junkline
	
	call system('wc -l '//trim(craft_info)//'defaults.pos > '// &
		trim(craft_info)//'crafts.txt')
	call system('ls -1 '//trim(craft_info)//'*.craft | wc -l >> '// &
		trim(craft_info)//'crafts.txt')
	
	open(scin,file=trim(craft_info)//'crafts.txt', &
		status='unknown',form='formatted')
		!	First line of crafts.txt, written above, is output of wc -l
		!	from defaults.pos, the number of lines in this file.
		read(scin,*) ndef_craft, junkline
		!	Second line is number of .craft files
		read(scin,*) naux_craft
	close(scin)
	
	!	Corrects for including header line in line count
	ndef_craft = ndef_craft - 1
	ncraft = ndef_craft + naux_craft
	
	return
end subroutine countcraft
