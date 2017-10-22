!
!	This file contains two subroutines:
!	asub
!	atsub
!
subroutine asub(abd,nband,n,ipvt,x,xi)
    !
    !    calculates a.x taking advantage of the banded
    !       structure of a
    !
    real abd(nband,n),x(n),xi(n)
    integer ipvt(nband)
    integer box
    !
    !     initialize product
    !
    do i=1,n
        xi(i)=0.
        !
        do box=1,nband
            ii=i+ipvt(box)
            if((ii.ge.1).and.(ii.le.n)) xi(i)=xi(i)+abd(box,ii)*x(ii)
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
subroutine atsub(abd,nband,n,ipvt,x,xi)
    !
    !    warning: symmetrical matrix assumed here
    !    calculates transpose(a).x
    !         taking advantage of the banded
    !          structure of a
    !
    real abd(nband,n),x(n),xi(n)
    integer ipvt(nband)
    integer box
    !
    !     initialize product
    !
    do i=1,n
        xi(i)=0.
        !
        do box=1,nband
            ii=i+ipvt(box)
            if((ii.ge.1).and.(ii.le.n))xi(i)=xi(i)+abd(m,ii)*x(ii)
        enddo
    enddo
    !
    !
    return
end
