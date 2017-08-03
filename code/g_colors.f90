!
!	This file includes six subroutines for assigning color lists:
!	cpclrs
!	isoclrs
!	isoclrs_hot
!	isorb
!	dfclrs
!	colram
!
subroutine cpclrs
    !
    !      8 basic colors, 4 extras for curve tracing, 15 more for isoplots
    !
    real,dimension(3,12) :: rgbv = &
    (/0.,0.,0., &     !black 
    1.,0.,0.,  &    !red 
    0.,1.,0.,  &   !green 
    0.,0.,1.,  &   !blue 
    0.,1.,1.,  &   !cyan 
    1.,0.,1.,  &   !magenta 
    1.,1.,0.,  &   !yellow 
    0.5, 1.0, 0.5, &    !bright green 
    1.0, 0.5, 1.0, &    !bright pink 
    0.3, 1.0, 1.0, &    !bright blue 
    1.0, 1.0, 0.3, &    !bright yellow 
    1.,1.,1./)     !white
    !
    call gscr(1,0,1.,1.,1.)
    !
    do i = 1, 12
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    !
    return
end
!
!
!*********************************
!
!
subroutine isoclrs
    real,dimension(3,16) :: rgbv = &
    (/0.,0.,0., .7,.7,.7, .75,.5,1., &
    .5,0.,1., 0.,0.,1.,  0.,.5,1., &
    0.,1.,1., 0.,1.,.6, 0.,.85,0., &
    .7,1.,0., 1.,1.,0., 1.,.75,0., &
    1.,.38,.38,  1.,0.,.75, &
    1.,0.,0.,  1.,1.,1./)
    call gscr(1,0,1.,1.,1.)
    do i=1,16
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    return
end
!
!
!*********************************
!
!
subroutine isoclrs_hot
    !
    !     arrays to plot colorbar
    !
    integer lind(20)
    real tcon(20)
    character*4 llbs(20)
    !
    real,dimension(3,21)::  rgbv = &
    (/0.,0.,0., .7,.7,.7, .75,.5,1., &
    .5,0.,1.,.25,0.,1.,  0.,0.,1., 0.,.25,1., &
    0.,.5,1., 0.,1.,1., 0.,1.,.6, 0.,.85,0., &
    0.3,.85,0., .7,1.,0., 1.,1.,0., 1.,.75,0., &
    1.,.38,.38,  1.,0.2,.2, &
    1.,0.1,0.1, 1.,0.,0., 0.85,0.,0.,1.,1.,1./)
    call gscr(1,0,1.,1.,1.)
    do i=1,21
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    return
end
!
!
!*********************************
!
!
subroutine isorb
	!
	!	creates a color map from the following RGB vectors:
	!
    real,dimension(3,16) :: rgbv = &
    (/0.,0.,0., .7,.7,.7, .75,.5,1., &
    .5,0.,1., 0.,0.,1.,  0.,.5,1., &
    0.,1.,1., 0.,1.,.6, &
    1.,.38,.38,  1.,0.,.38, &
    1.,0.,0.,  0.,0.65,0., &
    0.,.85,0., .7,1.,0., &
    1.,1.,0., 1.,.75,0./)
    call gscr(1,0,1.,1.,1.)
    do  i=1,16
        call gscr(1,i,rgbv(1,i),rgbv(2,i),rgbv(3,i))
    enddo
    return
end
!
!
!*********************************
!
!
subroutine dfclrs(ncol)
    !
    !  reads from the default input ($in on the cray) a color table of
    !  maxcol colors in free format 8-bit integer rgb values (i.e. 0 to 255)
    !  and selects ncol colors evenly distributed among the maxcol colors;
    !  these are bound to color indices 2 to ncol+1
    !
    parameter (maxcol=250)
    dimension irgbv(3,maxcol)
    data icall/0/
    save icall,irgbv
    if (ncol.le.0.or.ncol.gt.maxcol) then
        print*,'stop in dfclrs - input parameter ncol out of range'
        print*,'ncol = ',ncol
        stop
    endif
    if (icall.ne.1) then
        icall=1
        read(9,*) irgbv
    endif
    ainc=real(maxcol-1)/real(ncol-1)
    do jj=2,ncol+1
        j=nint(ainc*(jj-2))+1
        r=irgbv(1,j)/255.
        g=irgbv(2,j)/255.
        b=irgbv(3,j)/255.
        call gscr(1,jj,r,g,b)
    enddo
     return
end
!
!
!*********************************
!
!
subroutine colram(xcra,ycra,ncra,iaia,igia,naia)
    dimension xcra(*),ycra(*),iaia(*),igia(*)
    if (iaia(1).le.0) return
    !
    !  ..area numbering starts with 1, but color index numbering starts
    !    with 2 because 0 and 1 are reserved for background and foreground,
    !    respectively
    !
    call gsfaci(iaia(1))
    call gfa(ncra-1,xcra,ycra)
    return
end
