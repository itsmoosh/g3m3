program stack_gifs
	implicit none
	character*20, parameter :: img_dir = "figures/images/"
	integer, parameter :: input_f = 299

	character*16, allocatable :: qty_name(:)
	character*80, allocatable :: gif_path(:)
	character*1 :: box_char(9) = ['1','2','3','4','5','6','7','8','9']
	character*3 start_num
	character*3 end_num
	character*8 run_name
	integer n_qtys
	integer n_grids
	integer n_imgs

	character*4000, allocatable :: gif_cmds(:)
	character*4000, allocatable :: rm_cmds(:)

	integer i_box, i_qty, i_img

	open(input_f,file=trim(img_dir)//"gif_cmds.txt",status='unknown',form='formatted')
		read(input_f,*) n_imgs
		allocate( gif_cmds(n_imgs), rm_cmds(n_imgs) )
		do i_img=1, n_imgs
			read(input_f,'(A)') gif_cmds(i_img)
			read(input_f,'(A)') rm_cmds(i_img)
		enddo
	close(input_f)

	!$omp parallel do
	do i_img=1, n_imgs
		call system( trim(gif_cmds(i_img)) )
		call system( trim(rm_cmds(i_img)) )
	enddo
end program stack_gifs
