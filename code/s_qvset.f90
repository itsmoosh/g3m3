!
!	Set an array v of total dimension n to value a
!	MJS 4/9/19: This seems unnecessary, as this can be done with
!	v(:) = a
!
subroutine qvset(a,v,n)

	dimension v(1)
	
	!$omp  parallel do
	do i=1,n
		v(i)=a
	enddo

	return
end subroutine qvset
