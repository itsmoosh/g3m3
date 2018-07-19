subroutine zero_qty(bx,nx,ny,nz,n_grids)
    !
    !	zeroes out the quantity for all grid points
    !
    dimension qty(nx,ny,nz,n_grids)
    !
    !$omp  parallel do
	do box=1,n_grids
		do k=1,nz
		    do j=1,ny
		        do i=1,nx
		            qty(i,j,k,box)=0
		        enddo
		    enddo
		enddo
	enddo
    !
    return
end
