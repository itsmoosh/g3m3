!
!	Converts momenta into electron velocity for plotting purposes.
!
subroutine fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
	opx,opy,opz,orho,curx,cury,curz, &
	evx,evy,evz,tvx,tvy,tvz, &
	nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso,reynolds)

	implicit none

	integer, intent(in) :: nx, ny, nz, n_grids, box
	real, intent(in) :: qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids), &
	qpz(nx,ny,nz,n_grids),qrho(nx,ny,nz,n_grids), &
	hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids), &
	hpz(nx,ny,nz,n_grids),hrho(nx,ny,nz,n_grids), &
	opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids), &
	opz(nx,ny,nz,n_grids),orho(nx,ny,nz,n_grids), &
	curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz)
	real, intent(in) :: rmassq, rmassh, rmasso, reynolds
	real, intent(out) :: evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz), &
	tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz)

	integer i,j,k
	real tot_numden, tden, qden, hden, oden
	
	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				qden=(qrho(i,j,k,box)+0.000001) / rmassq
				hden=(hrho(i,j,k,box)+0.000001) / rmassh
				oden=(orho(i,j,k,box)+0.000001) / rmasso
				tot_numden = qden + hden + oden
				tden = qrho(i,j,k,box) + hrho(i,j,k,box) + orho(i,j,k,box)
				
				!	Keep separate the ion and current components
				
				tvx(i,j,k) = ( qpx(i,j,k,box) + hpx(i,j,k,box) &
					+ opx(i,j,k,box) ) / tden

				evx(i,j,k) = tvx(i,j,k) - curx(i,j,k)/tot_numden/reynolds

				tvy(i,j,k) = ( qpy(i,j,k,box) + hpy(i,j,k,box) &
					+ opy(i,j,k,box) ) / tden

				evy(i,j,k) = tvy(i,j,k) - cury(i,j,k)/tot_numden/reynolds
	
				tvz(i,j,k) = ( qpz(i,j,k,box) + hpz(i,j,k,box) &
					+ opz(i,j,k,box) ) / tden

				evz(i,j,k) = tvz(i,j,k) - curz(i,j,k)/tot_numden/reynolds
			enddo
		enddo
	enddo
	
	return
end subroutine fnd_evel
