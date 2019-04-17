!
!	Sparse matrix solver:
!		banded_matrix is the banded matrix
!		nbands is the number of bands
!		divB_band_pos is the position of the bands relative to the diagonal
!
!		npts is the order of the original matrix
!		divB(npts) the real vector for the right hand side
!		soln_vec(npts) is the solution vector
!
!		residual2 is the sum of the squares of the components of the residual
!		vector banded_matrix.soln_vec - divB
!
!	If this is not small then the matrix is numerically singular and the 
!	solution represents a least	squares best approximation
!
!	divB_errlim is an error limit on the solution
!
recursive subroutine sparse_solver( nbands, npts, divB_band_pos, divB_errlim, &
	recurse, redo_lim, divB, soln_vec, residual2, banded_matrix )

	implicit none

	!	*********
	!	Arguments
	!	*********

	integer, intent(in)	:: nbands, npts
	integer, intent(in)	:: divB_band_pos(nbands)
	real, intent(in)	:: divB_errlim
	logical, intent(in) :: recurse
	integer, intent(in) :: redo_lim
	real, intent(in)	:: divB(npts)

	real, intent(out)	:: soln_vec(npts), residual2
	real, intent(inout)	:: banded_matrix(nbands,npts)

	!	***************
	!	Local variables
	!	***************	

	!	Work arrays for calculating matrix solutions
	real g(npts), h(npts), xi(npts), xj(npts)

	integer	:: redos = 0
	integer i,j,k,iter
	real errlim2, sum_Bsqd, divB_resd
	real anum, aden, gg, dgg, gam

	!	Criterion for sum-squared residuals
	!	and number of restarts attempted internally

	errlim2 = npts * divB_errlim**2
	if(recurse) then
		redos = redo_lim
	else
		redos = 0
	endif
	
	redos = redos + 1
	call banded_tpose(nbands,npts,divB_band_pos,banded_matrix,soln_vec,xi)
	
	!	Add the magnitude of the right side and find the residue
	
	divB_resd = 0.
	sum_Bsqd = 0.
	do j = 1, npts
		sum_Bsqd = sum_Bsqd + divB(j)**2
		xi(j) = xi(j) - divB(j)
		divB_resd = divB_resd + xi(j)**2
	enddo
	
	call banded_tpose(nbands,npts,divB_band_pos,banded_matrix,xi,g)
	g(:)=-g(:)
	h(:)=g(:)
	
	!	Main iterative loop
	
	!	10*npts is max # of iterations
	do iter=1, 10*npts

		call banded_tpose(nbands,npts,divB_band_pos,banded_matrix,h,xi)
		
		!	Calculate the gradient
		
		anum=0.
		aden=0.
		
		do j = 1, npts
			anum = anum + g(j)*h(j)
			aden = aden + xi(j)**2
		enddo
		if(aden.eq.0.) write(*,*) 'Very singular matrix in sparse.'
		
		anum = anum/aden
		do j = 1, npts
			xi(j)=soln_vec(j)
			soln_vec(j)=soln_vec(j)+anum*h(j)
		enddo
		
		call banded_tpose(nbands,npts,divB_band_pos,banded_matrix,soln_vec,xj)
		
		residual2 = 0.
		do j = 1, npts
			xj(j) = xj(j) - divB(j)
			residual2 = residual2 + xj(j)**2
		enddo
		
		!	Test for convergence and exit if okay
		
		if( (residual2.eq.divB_resd) .or. (residual2.le.sum_Bsqd*errlim2) ) then
			return
		endif
		
		!	Test to see whether solution is improving
		!	and restart if necessary
		
		if(residual2.gt.divB_resd) then
			do j = 1, npts
				soln_vec(j) = xi(j)
			enddo
			
			!	Return if too many restarts - hitting roundoff error
			
			if(redos.ge.3) return
			
			call sparse_solver( nbands, npts, divB_band_pos, divB_errlim, &
				.true., redos, divB, soln_vec, residual2, banded_matrix )
		endif
		
		!	Compute gradient for next iteration
		
		divB_resd = residual2
		call banded_tpose(nbands,npts,divB_band_pos,banded_matrix,xj,xi)

		gg=0.
		dgg=0.
		do j = 1, npts
			gg=gg+g(j)**2
			dgg=dgg+(xi(j)+g(j))*xi(j)
		enddo
		
		!	Test to see if you have a solution and return if okay
		
		if(gg.eq.0.) return
		
		gam=dgg/gg
		do j = 1, npts
			g(j)=-xi(j)
			h(j)=g(j)+gam*h(j)
		enddo
		
	enddo
	
	!	Never found solution if you get here
	
	write(*,*) 'Too many iterations. No solution for sparse.'
	return
end subroutine sparse_solver
