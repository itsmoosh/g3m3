!
!	Converts ion momenta into net ion velocity for plotting purposes.
!
subroutine fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
	opx,opy,opz,orho,vx,vy,vz,nx,ny,nz,n_grids,box, &
	rmassq,rmassh,rmasso)

	implicit none

	integer, intent(in) :: nx, ny, nz, n_grids, box
	real, intent(in) :: qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids), &
		qpz(nx,ny,nz,n_grids),qrho(nx,ny,nz,n_grids), &
		hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids), &
		hpz(nx,ny,nz,n_grids),hrho(nx,ny,nz,n_grids), &
		opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids), &
		opz(nx,ny,nz,n_grids),orho(nx,ny,nz,n_grids)
	real, intent(in) :: rmassq, rmassh, rmasso
	real, intent(out) :: vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)

	integer i,j,k
	real tden
	
	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				tden = qrho(i,j,k,box) + hrho(i,j,k,box) &
					+ orho(i,j,k,box)
				
				vx(i,j,k) = ( qpx(i,j,k,box) + hpx(i,j,k,box) &
					+ opx(i,j,k,box) ) / tden
				vy(i,j,k) = ( qpy(i,j,k,box) + hpy(i,j,k,box) &
					+ opy(i,j,k,box) ) / tden
				vz(i,j,k) = ( qpz(i,j,k,box) + hpz(i,j,k,box) &
					+ opz(i,j,k,box) ) / tden
			enddo
		enddo
	enddo
	
	return
end subroutine fnd_vtot
