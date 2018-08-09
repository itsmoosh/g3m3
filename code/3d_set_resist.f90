subroutine set_resist(rst,nx,ny,nz,mbndry,resist, &
    ijzero,numzero,ijmid,nummid,ijsrf,numsrf, ijlayers,numlayerpts, &
    msrf,mmid,mzero,m_layers,b0,conductivities,layer_id)
    !
    !      set resistance around the earth :
    !        include dayside and auroral conductivities
    !      magnitude set by resist
    !      width is set by del_ang=3.3 degrees = 0.058 rads
    !      radial dropoff as alpha=-8
    !      shifted of dipole axis by 2.5 degress
    !      radius of 22.5 degrees
    !
	integer, intent(in) :: nx,ny,nz,mbndry,msrf,mmid,mzero,m_layers
    dimension rst(nx,ny,nz,mbndry)
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
        ijzero(mbndry,3,mzero), ijlayers(mbndry,3,m_layers)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry),numlayerpts(mbndry)
    integer box,n,layer
    !
    !     resistivity characteristics
    !        initialize
    !
    rst=0.
    !     write(*,*)'set_resist: mbndry, nx, ny, nz, b0'
	!		write(*,*) mbndry, nx, ny, nz, b0
    !     write(*,*) 'numsrf, nummid, numzero'
	!		write(*,*) numsrf, nummid, numzero
    !
    !     interior resistivity
    !
    !$omp  parallel do
    do n=1,numzero(mbndry)
        do box=1,mbndry
            i=ijzero(box,1,n)
            j=ijzero(box,2,n)
            k=ijzero(box,3,n)
            rst(i,j,k,box)=b0/resist
        enddo
    enddo
    !
    !     lower ionosphere resistivity
    !
    !$omp  parallel do
    do n=1,nummid(mbndry)
        do box=1,mbndry
            i=ijmid(box,1,n)
            j=ijmid(box,2,n)
            k=ijmid(box,3,n)
            rst(i,j,k,box)=0.5*b0/resist
        enddo
    enddo
    !
    !     upper ionosphere
    !
    !$omp  parallel do
    do n=1,numsrf(mbndry)
        do box=1,mbndry
            i=ijsrf(box,1,n)
            j=ijsrf(box,2,n)
            k=ijsrf(box,3,n)
            rst(i,j,k,box)=0.125*b0/resist
        enddo
    enddo 
	!
	!	planetary conducting layers
	!
    !$omp  parallel do
    do n=1,numlayerpts(mbndry)
        do box=1,mbndry
			layer=layer_id(n)
            i=ijlayers(box,1,n)
            j=ijlayers(box,2,n)
            k=ijlayers(box,3,n)
            rst(i,j,k,box)=1./conductivities(layer)
        enddo
    enddo 
    return
end
