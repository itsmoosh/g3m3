!
!	Calculates transpose(a).x
!	Warning: symmetrical matrix assumed here,
!	taking advantage of the banded structure of a
!
subroutine atsub(abd,nbands,n,ipvt,x,xi)

	real abd(nbands,n),x(n),xi(n)
	integer ipvt(nbands)
	integer band
	!	Initialize product
	do i=1,n
		xi(i)=0.
		!
		do band=1, nbands
			ii=i+ipvt(band)
			if((ii.ge.1).and.(ii.le.n)) xi(i)=xi(i)+abd(band,ii)*x(ii)
		enddo
	enddo

	return
end subroutine atsub
