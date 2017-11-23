subroutine limcraft(craftpos,ntimes,vals,ndef_craft,ncraft,re_equiv,n_grids,grid_minvals,grid_maxvals)
    !
    !	Tests to see whether trajectory craft locations are always in the system and
    !	resets their positions if not.
    !
	implicit none

    integer n, m_recnum
	integer, intent(in) :: n_grids
	integer, intent(in) :: grid_minvals(3,n_grids)
	integer, intent(in) :: grid_maxvals(3,n_grids)
	integer, intent(in) :: ntimes(ncraft,2)
	integer, intent(in) :: ndef_craft
	integer, intent(in) :: ncraft
	integer, intent(in) :: vals

	real, intent(in) :: re_equiv

	real, intent(inout) :: craftpos(4,ncraft,vals)

	real abit

	abit = 0.001
	
    do n = ndef_craft+1,ncraft
		do m_recnum = 1, ntimes(n,2)
			craftpos(1:3,n,m_recnum) = amax1( craftpos(1:3,n,m_recnum), (grid_minvals(:,n_grids) +abit)*re_equiv )
			craftpos(1:3,n,m_recnum) = amin1( craftpos(1:3,n,m_recnum), (grid_maxvals(:,n_grids) -abit)*re_equiv )
		enddo
    enddo

    return
end subroutine limcraft
