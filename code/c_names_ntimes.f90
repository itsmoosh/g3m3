subroutine names_ntimes(craft_info, scin, ncraft, ndef_craft, naux_craft, nheadlines, craftnames, ntimes)
	!
	!	Reads the name and number of measurements in each spacecraft trajectory file.
	!
	implicit none

	character*32, intent(in)	:: craft_info
	integer, intent(in)			:: scin
	integer, intent(in)			:: ncraft
	integer, intent(in)			:: ndef_craft
	integer, intent(in)			:: naux_craft
	integer, intent(in)			:: nheadlines

	character*8, intent(inout)	:: craftnames(ncraft)
	integer, intent(inout)		:: ntimes(ncraft,2)

	character*120 junkline
	integer n

	character,parameter :: tab = char(9)


	!	Read auxiliary spacecraft names into craftnames array
	!	Step 1: ls .craft names into crafts.txt file. 
	!		This file will contain only .craft names, not default craft.
	!		We use wc -l *.craft so that we get numbers of lines as well as file names.
	call system ('wc -l '//trim(craft_info)//'*.craft > '//trim(craft_info)//'crafts.txt')
	!
	open(scin,file=trim(craft_info)//'crafts.txt',status='unknown',form='formatted')
		do n=ndef_craft+1,ncraft
			!	Step 2: Read .craft filenames into placeholder string
			read(scin,'(A60)') junkline
			read(junkline, *) ntimes(n,2)
			ntimes(n,2) = ntimes(n,2) - nheadlines	!	Excludes header line in craft file
			!	Step 3: Print the filename before the .craft extension into craftnames array
			craftnames(n) = junkline(index(junkline,'/')+1:index(junkline,'.')-1)
		enddo
	close(scin)

	!@@@@@@@@@@@@@
	!Default craft
	!@@@@@@@@@@@@@

	open(scin+1,file=trim(craft_info)//'defaults.pos',status='unknown',form='formatted')
	read(scin+1,*) junkline	!	Skip header line, only one line in default craft file
		!
		do n=1,ndef_craft
			read(scin+1,*) craftnames(n), junkline
			ntimes(n,2) = -1	!	We don't have a set number of recordings for default craft
		enddo
	close(scin+1)

	return
end subroutine names_ntimes
