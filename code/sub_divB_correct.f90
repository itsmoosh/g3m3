!
!	Solves the matrix equation
!		a.x=b
!	using using conjugate gradient method from
!	numerical recipes
!	3: d solution of poisson's eqn :
!		npts=nx*ny*nz
!
!	mapping is such that for position i,j,k:
!	n =i+nx*(j-1) +nx*ny*(k-1)
!
!	original solutions done on small grid
!
!	store only non-zero bands of matrix a in
!	array banded_matrix, and the position of these bands is
!	stored in divB_band_pos
!	divB is data vector
!	soln_vec is the solution we are seeking
!	xi,xj,g,h are work arrays needed for inversion
!
!	nband = 7 for 3d and 5 for 2d
!
!	size of array abd hard wired to be nband,nx*ny*nz=npts
!
subroutine divB_correct( &
	nx,ny,nz,n_grids, box, nbands, divB_band_pos, xspac, divB_errlim, &
	bx,by,bz )

	implicit none

	integer, parameter :: dp = kind(1.d0)
	
	!	*********
	!	Arguments
	!	*********

	integer, intent(in)	:: nx, ny, nz, n_grids, box, nbands
	integer, intent(in)	:: divB_band_pos(nbands)
	real, intent(in)	:: xspac(n_grids), divB_errlim

	real, intent(inout), dimension(nx,ny,nz,n_grids)	:: bx, by, bz

	!	***********
	!	Work arrays
	!	***********

	real, dimension(nx,ny,nz)	:: dbx, dby, dbz, poten
	real, dimension(nx*ny*nz)	:: divB, soln_vec, g, h, xi, xj
	real, dimension(nbands,nx*ny*nz)	:: banded_matrix

	integer i,j,k,n,npts
	real rx,ry,rz, scale
	!	residual2 and redo_lim are placeholders, used only in the recursive
	!	mode for sparse_solver.
	real residual2
	integer redo_lim

	integer bd_low3, bd_low2, bd_low1, bd_diag, bd_upp1, bd_upp2, bd_upp3
		
	rx=xspac(box)
	ry=xspac(box)
	rz=xspac(box)
	npts = nx*ny*nz
	bd_low3 = 7
	bd_low2 = 6
	bd_low1 = 5
	bd_diag = 4
	bd_upp1 = 3
	bd_upp2 = 2
	bd_upp3 = 1

	!write(*,*)'calling divb',nband,npts,nx,ny,nz,n_grids,box	

	!	Make the banded matrix
	
	!	Initialize elements of banded matrix
	!	for diagonal or element n requires
	!	upper2 (mu2) at n+nx  -> limits nx+1,npts
	!	upper1 (mu1) at n+1 -> limits 2,npts
	!	lower2 (ml2) at n-nx -> limits 1,npts-nx
	!	lower1 (ml1) at n-1 -> limits 1,npts-1
	
	!	Other than values identified below, remaining values in banded_matrix
	!	are all zero:
	banded_matrix(:,:) = 0.

	!	Center band
	do  n = 1, npts
		banded_matrix(bd_diag,n) = -6.
		divB(n) = 0.
	enddo
	
	!	Off-diagonal elements
	do k = 1, nz
		do j = 1, ny
			do i = 2, nx
				n = i + nx*(j-1) + nx*ny*(k-1)
				banded_matrix(bd_upp1,n) = 1.
			enddo
			
			do i = 1, nx-1
				n = i + nx*(j-1) + nx*ny*(k-1)
				banded_matrix(bd_low1,n) = 1.
			enddo
		enddo

		do j = nx+1, nx*ny
			n = j+nx*ny*(k-1)
			banded_matrix(bd_upp2,n) = 1.
		enddo
		
		do j = 1, nx*ny - nx
			n=j+nx*ny*(k-1)
			banded_matrix(bd_low2,n) = 1.
		enddo
	enddo
	
	do j = ny*nx+1, npts
		banded_matrix(bd_upp3,j) = 1.
	enddo
	do j = 1, npts - nx*ny
		banded_matrix(bd_low3,j) = 1.
	enddo
	
	!	Set boundary conditions: symmetric matrix
	!   so we don't have to worry about boundary conditions
	
	!	Find divB and load into matrix
	
	do k=2,nz-1
		do j=2,ny-1
			do i=2,nx-1
			
			n = i+nx*(j-1) + nx*ny*(k-1)
			divB(n) = ( (bx(i+1,j,k,box)-bx(i-1,j,k,box)) &
			+ (by(i,j+1,k,box)-by(i,j-1,k,box)) &
			+ (bz(i,j,k+1,box)-bz(i,j,k-1,box)) ) / (2.*rx)
			
			enddo
		enddo
	enddo
	
	call sparse_solver( nbands, npts, divB_band_pos, divB_errlim, &
		.false., redo_lim, divB, soln_vec, residual2, banded_matrix )
	
	!	Convert potential onto grid
	poten(:,:,:) = 0.
	scale=rx**2
	do k=1,nz
		do j=1,ny
			do i=1,nx
				n=i+nx*(j-1)+nx*ny*(k-1)
				poten(i,j,k)=scale*soln_vec(n)
			enddo
		enddo
	enddo

	!	Calculate perturbed magnetic field
	do k=2,nz-1
		do j=2,ny-1
			do i=2,nx-1
				dbx(i,j,k)=(poten(i+1,j,k)-poten(i-1,j,k))/(2.*rx)
				dby(i,j,k)=(poten(i,j+1,k)-poten(i,j-1,k))/(2.*ry)
				dbz(i,j,k)=(poten(i,j,k+1)-poten(i,j,k-1))/(2.*rz)

				!	Modified magnetic field
				bx(i,j,k,box)=bx(i,j,k,box)-dbx(i,j,k)
				by(i,j,k,box)=by(i,j,k,box)-dby(i,j,k)
				bz(i,j,k,box)=bz(i,j,k,box)-dbz(i,j,k)
			enddo
		enddo
	enddo
	
	return
end subroutine divB_correct
