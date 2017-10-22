!
!	This file contains two subroutines:
!	asub
!	atsub
!
subroutine asub(abd,nbands,n,ipvt,x,xi)
    !
    !    calculates a.x taking advantage of the banded
    !       structure of a
    !
    real abd(nbands,n),x(n),xi(n)
    integer ipvt(nbands)
    integer band
    !
    !     initialize product
    !
    do i=1,n
        xi(i)=0.
        !
        do band=1,nbands
            ii=i+ipvt(band)
            if((ii.ge.1).and.(ii.le.n)) xi(i)=xi(i)+abd(band,ii)*x(ii)
        enddo
    enddo
    !
    return
end
!
!
!	****************************************
!
!
subroutine atsub(abd,nbands,n,ipvt,x,xi)
    !
    !    warning: symmetrical matrix assumed here
    !    calculates transpose(a).x
    !         taking advantage of the banded
    !          structure of a
    !
    real abd(nbands,n),x(n),xi(n)
    integer ipvt(nbands)
    integer band
    !
    !     initialize product
    !
    do i=1,n
        xi(i)=0.
        !
        do band=1,nbands
            ii=i+ipvt(band)
            if((ii.ge.1).and.(ii.le.n)) xi(i)=xi(i)+abd(band,ii)*x(ii)
        enddo
    enddo
    !
    !
    return
end
