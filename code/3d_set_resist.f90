subroutine set_resist(rst,nx,ny,nz,mbndry,resist, &
    ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
    msrf,mmid,mzero,b0)
    !
    !      set resistance around the earth :
    !        include dayside and auroral conductivities
    !      magntiude set by resist
    !      width is set by del_ang=3.3 degrees = 0.058 rads
    !      radial dropoff as alpha=-8
    !      shifted of dipole axis by 2.5 degress
    !      radius of 22.5 degrees
    !
    dimension rst(nx,ny,nz,mbndry)
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
        ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    !
    !     resistivity characteristics
    !        initialize
    !
    rst=0.
    !     write(6,*)'set_resist',mbndry,nx,ny,nz,b0
    !     write(6,*)numsrf,nummid,numzero
    !
    !     interior resistivity
    !
    do n=1,numzero(mbndry)
        do m=1,mbndry
            i=ijzero(m,1,n)
            j=ijzero(m,2,n)
            k=ijzero(m,3,n)
            rst(i,j,k,m)=b0/resist
        enddo
    enddo
    !
    !     lower ionosphere resistivity
    !
    do n=1,nummid(mbndry)
        do m=1,mbndry
            i=ijmid(m,1,n)
            j=ijmid(m,2,n)
            k=ijmid(m,3,n)
            rst(i,j,k,m)=0.5*b0/resist
        enddo
    enddo
    !
    !     upper ionosphere
    !
    do n=1,numsrf(mbndry)
        do m=1,mbndry
            i=ijsrf(m,1,n)
            j=ijsrf(m,2,n)
            k=ijsrf(m,3,n)
            rst(i,j,k,m)=0.125*b0/resist
        enddo
    enddo 
    return
end
