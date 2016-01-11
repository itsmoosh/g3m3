subroutine divb_correct(bx,by,bz,dbx,dby,dbz,poten, &
    b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,ngrd,m, &
    xspac)
    !
    !     this program solves the matrix equation
    !            a.x=b
    !      using using conjuage gradient method from
    !          numerical recipes
    !     3: d solution of poisson's eqn :
    !          ntot=nx*ny*nz
    !
    !       mapping is such that for position i,j,k:
    !          n =i+nx*(j-1) +nx*ny*(k-1)
    !
    !     original solutions done on small grid
    !     parameter (ntot=893101,nx=121,ny=121,nz=61,nband=7)
    !
    !     store only non-zero bnads of matrix a in
    !       array abd, and the position of these bands is
    !       stored in ipvt
    !       abd gives values of bands in matrix
    !       b is data vector
    !       x is the solution we are seeking
    !       xi,xj,g,h are work arrays needed for inversion
    !
    !           nband = 7 for 3d and 5 for 2d
    !
    !     size of array abd hard wired to be nband,nx*ny*nz=ntot
    !
    !
    !     common /space/abd(7,893101)
    real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
    real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz), &
        poten(nx,ny,nz)
    real xspac(ngrd)
    !
    real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
    integer ipvt(nband)
    !
    real, allocatable, dimension(:,:) :: abd
    allocate(abd(7,nx*ny*nz))
    !
    !write(6,*)'calling divb',nband,ntot,nx,ny,nz,ngrd,m
    !
    rx=xspac(m)
    ry=xspac(m)
    rz=xspac(m)
    !
    !
    !      make the banded matrix
    !
    !            highest band position :  mu2, mu1
    !            lowest band position :  ml2, ml1
    !            diagonal band position :  m0
    !
    !             note in matrix notation top row is 1
    !                                    bottom row is nband
    !
    m0=4
    ml1=5
    ml2=6
    ml3=7
    !
    mu1=3
    mu2=2
    mu3=1
    !
    !
    !      position of the bands
    !
    ipvt(m0)=0
    ipvt(ml1)=-1
    ipvt(ml2)=-nx
    ipvt(ml3)=-nx*ny
    ipvt(mu1)=1
    ipvt(mu2)=nx
    ipvt(mu3)=nx*ny
    !
    !      initialize elements of banded matrix
    !         for diagonal or element n requires
    !            mu2 at n+nx  -> limits nx+1,ntot
    !            mu1 at n+1 -> limits 2,ntot
    !            ml2 at n-nx -> limits 1,ntot-nx
    !            ml1 at n-1 -> limits 1,ntot-1
    !
    !
    !       initial bands to zero
    !
    !
    abd=0.
    !
    !       do centerband
    !
    do  n = 1, ntot
        abd(m0,n)=-6.
        b(n)=0.
    enddo
    !
    !     do off diagonal elements
    !
    do k=1,nz
        do j = 1,ny
            do i = 2, nx
                n=i+nx*(j-1)+nx*ny*(k-1)
                abd(mu1,n)=1.
            enddo
            do i = 1, nx-1
                n=i+nx*(j-1)+nx*ny*(k-1)
                abd(ml1,n)=1.
            enddo
        enddo
    enddo
    !
    do k=1,nz
        do j = nx+1, nx*ny
            n=j+nx*ny*(k-1)
            abd(mu2,n)=1.
        enddo
        do j = 1, nx*ny-nx
            n=j+nx*ny*(k-1)
            abd(ml2,n)=1.
        enddo
    enddo
    !
    do j = ny*nx+1, ntot
        abd(mu3,j)=1.
    enddo
    do j = 1, ntot-nx*ny
        abd(ml3,j)=1.
    enddo
    !
    !      set boundary conditions : symetrix matrix
    !         so don't have to worry about boundary conditions
    !
    !        abd(ml,1)=0.
    !        abd(mu,ntot)=0.
    !
    !      find div b and load into matrix
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                !
                n=i+nx*(j-1)+nx*ny*(k-1)
                divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m)) &
                    +(by(i,j+1,k,m)-by(i,j-1,k,m)) &
                    +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
                !
                b(n)=divb
            enddo
        enddo
    enddo
    !
    call sparse(abd,nband,ntot,b,x,ipvt,rsq, &
        g,h,xi,xj,.false.)
    !
    !     convert potential onto grid
    !
    poten=0.
    scale=rx**2
    do k=1,nz
        do j=1,ny
            do i=1,nx
                n=i+nx*(j-1)+nx*ny*(k-1)
                poten(i,j,k)=scale*x(n)
            enddo
        enddo
    enddo
    !
    !    calculate  perturbed magnetic field
    !
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                dbx(i,j,k)=(poten(i+1,j,k)-poten(i-1,j,k))/(2.*rx)
                dby(i,j,k)=(poten(i,j+1,k)-poten(i,j-1,k))/(2.*ry)
                dbz(i,j,k)=(poten(i,j,k+1)-poten(i,j,k-1))/(2.*rz)
            enddo
        enddo
    enddo
    !
    !
    !     modified magnetic field
    !
    
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                bx(i,j,k,m)=bx(i,j,k,m)-dbx(i,j,k)
                by(i,j,k,m)=by(i,j,k,m)-dby(i,j,k)
                bz(i,j,k,m)=bz(i,j,k,m)-dbz(i,j,k)
            enddo
        enddo
    enddo
    !
    if(allocated(abd)) deallocate(abd)
    return
end

!
!    **********************************************************
!
subroutine divb_correct_n(bx,by,bz,dbx,dby,dbz,poten, &
    b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,ngrd,m, &
    xspac)
    !
    !     this program solves the matrix equation
    !            a.x=b
    !      using using conjuage gradient method from
    !          numerical recipes
    !     3: d solution of poisson's eqn :
    !          ntot=nx*ny*nz
    !
    !       mapping is such that for position i,j,k:
    !          n =i+nx*(j-1) +nx*ny*(k-1)
    !
    !     original solutions done on small grid
    !     parameter (ntot=893101,nx=121,ny=121,nz=61,nband=7)
    !
    !     store only non-zero bnads of matrix a in
    !       array abd, and the position of these bands is
    !       stored in ipvt
    !       abd gives values of bands in matrix
    !       b is data vector
    !       x is the solution we are seeking
    !       xi,xj,g,h are work arrays needed for inversion
    !
    !           nband = 7 for 3d and 5 for 2d
    !
    !     size of array abd hard wired to be nband,nx*ny*nz=ntot
    !
    !
    common /space_n/abd(7,117649)
    real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
    real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz), &
        poten(nx,ny,nz)
    real xspac(ngrd)
    !
    real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
    integer ipvt(nband)
    !
    !      write(6,*)'calling divb',nband,ntot,nx,ny,nz,ngrd,m
    !
    rx=xspac(m)
    ry=xspac(m)
    rz=xspac(m)
    !
    !
    !      make the banded matrix
    !
    !            highest band position :  mu2, mu1
    !            lowest band position :  ml2, ml1
    !            diagonal band position :  m0
    !
    !             note in matrix notation top row is 1
    !                                    bottom row is nband
    !
    m0=4
    ml1=5
    ml2=6
    ml3=7
    !
    mu1=3
    mu2=2
    mu3=1
    !
    !
    !      position of the bands
    !
    ipvt(m0)=0
    ipvt(ml1)=-1
    ipvt(ml2)=-nx
    ipvt(ml3)=-nx*ny
    ipvt(mu1)=1
    ipvt(mu2)=nx
    ipvt(mu3)=nx*ny
    !
    !      initialize elements of banded matrix
    !         for diagonal or element n requires
    !            mu2 at n+nx  -> limits nx+1,ntot
    !            mu1 at n+1 -> limits 2,ntot
    !            ml2 at n-nx -> limits 1,ntot-nx
    !            ml1 at n-1 -> limits 1,ntot-1
    !
    !
    !       initial bands to zero
    !
    !
    abd=0.
    !
    !       do centerband
    !
    do  n = 1, ntot
        abd(m0,n)=-6.
        b(n)=0.
    enddo
    !
    !     do off diagonal elements
    !
    do k=1,nz
        do j = 1,ny
            do i = 2, nx
                n=i+nx*(j-1)+nx*ny*(k-1)
                abd(mu1,n)=1.
            enddo
            do i = 1, nx-1
                n=i+nx*(j-1)+nx*ny*(k-1)
                abd(ml1,n)=1.
            enddo
        enddo
    enddo
    !
    do k=1,nz
        do j = nx+1, nx*ny
            n=j+nx*ny*(k-1)
            abd(mu2,n)=1.
        enddo
        do j = 1, nx*ny-nx
            n=j+nx*ny*(k-1)
            abd(ml2,n)=1.
        enddo
    enddo
    !
    do j = ny*nx+1, ntot
        abd(mu3,j)=1.
    enddo
    do j = 1, ntot-nx*ny
        abd(ml3,j)=1.
    enddo
    !
    !      set boundary conditions : symetrix matrix
    !         so don't have to worry about boundary conditions
    !
    !        abd(ml,1)=0.
    !        abd(mu,ntot)=0.
    !
    !      find div b and load into matrix
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                !
                n=i+nx*(j-1)+nx*ny*(k-1)
                n=i+nx*(j-1)+nx*ny*(k-1)
                divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m)) &
                    +(by(i,j+1,k,m)-by(i,j-1,k,m)) &
                    +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
                !
                b(n)=divb
            enddo
        enddo
    enddo
    !
    call sparse(abd,nband,ntot,b,x,ipvt,rsq, &
        g,h,xi,xj,.false.)
    !
    !     convert potential onto grid
    !
    poten=0.
    scale=rx**2
    do k=1,nz
        do j=1,ny
            do i=1,nx
                n=i+nx*(j-1)+nx*ny*(k-1)
                poten(i,j,k)=scale*x(n)
            enddo
        enddo
    enddo
    !
    !    calculate  perturbed magnetic field
    !
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                dbx(i,j,k)=(poten(i+1,j,k)-poten(i-1,j,k))/(2.*rx)
                dby(i,j,k)=(poten(i,j+1,k)-poten(i,j-1,k))/(2.*ry)
                dbz(i,j,k)=(poten(i,j,k+1)-poten(i,j,k-1))/(2.*rz)
            enddo
        enddo
    enddo
    !
    !
    !     modified magnetic field
    !
    
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                bx(i,j,k,m)=bx(i,j,k,m)-dbx(i,j,k)
                by(i,j,k,m)=by(i,j,k,m)-dby(i,j,k)
                bz(i,j,k,m)=bz(i,j,k,m)-dbz(i,j,k)
            enddo
        enddo
    enddo
    !
    return
end

!
!     *****************************************************
!
subroutine divb_correct_tst(bx,by,bz,dbx,dby,dbz,poten, &
    b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,ngrd,m, &
    tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,xspac,rearth)
    !
    !     this program solves the matrix equation
    !            a.x=b
    !      using using conjuage gradient method from
    !          numerical recipes
    !     3: d solution of poisson's eqn :
    !          ntot=nx*ny*nz
    !
    !       mapping is such that for position i,j,k:
    !          n =i+nx*(j-1) +nx*ny*(k-1)
    !
    !     original solutions done on small grid
    !     parameter (ntot=893101,nx=121,ny=121,nz=61,nband=7)
    !
    !     store only non-zero bnads of matrix a in
    !       array abd, and the position of these bands is
    !       stored in ipvt
    !       abd gives values of bands in matrix
    !       b is data vector
    !       x is the solution we are seeking
    !       xi,xj,g,h are work arrays needed for inversion
    !
    !           nband = 7 for 3d and 5 for 2d
    !
    !     size of array abd hard wired to be nband,nx*ny*nz=ntot
    !
    !
    common /space/abd(7,893101)
    real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
    real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz), &
        poten(nx,ny,nz)
    real grd_xmin(ngrd),grd_xmax(ngrd),grd_ymin(ngrd),grd_ymax(ngrd), &
        grd_zmin(ngrd),grd_zmax(ngrd),xspac(ngrd)
    real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz), &
        tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2)
    !
    real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
    integer ipvt(nband)
    !
    write(6,*)'calling divb',nband,ntot,nx,ny,nz,ngrd,m
    !
    rx=xspac(m)
    ry=xspac(m)
    rz=xspac(m)
    !
    ymin=grd_ymin(m)+ry
    ymax=grd_ymax(m)-ry
    zmin=grd_zmin(m)+rz
    zmax=grd_zmax(m)-rz
    xmin=grd_xmin(m)+rx
    xmax=grd_xmax(m)-rx
    !
    xcut=(xmin+xmax)/2.
    !
    write(6,*)xmin,xmax,ymin,ymax,zmin,zmax
    !
    !      make the banded matrix
    !
    !            highest band position :  mu2, mu1
    !            lowest band position :  ml2, ml1
    !            diagonal band position :  m0
    !
    !             note in matrix notation top row is 1
    !                                    bottom row is nband
    !
    m0=4
    ml1=5
    ml2=6
    ml3=7
    !
    mu1=3
    mu2=2
    mu3=1
    !
    !
    !      position of the bands
    !
    ipvt(m0)=0
    ipvt(ml1)=-1
    ipvt(ml2)=-nx
    ipvt(ml3)=-nx*ny
    ipvt(mu1)=1
    ipvt(mu2)=nx
    ipvt(mu3)=nx*ny
    !
    !      initialize elements of banded matrix
    !         for diagonal or element n requires
    !            mu2 at n+nx  -> limits nx+1,ntot
    !            mu1 at n+1 -> limits 2,ntot
    !            ml2 at n-nx -> limits 1,ntot-nx
    !            ml1 at n-1 -> limits 1,ntot-1
    !
    !
    !       initial bands to zero
    !
    !       do  mm=1,nband
    !       do  n = 1, ntot
    !        abd(mm,n)=0.
    !       enddo
    !       enddo
    !
    abd=0.
    !
    !       do centerband
    !
    do  n = 1, ntot
        abd(m0,n)=-6.
        b(n)=0.
    enddo
    !
    !     do off diagonal elements
    !
    do k=1,nz
        do j = 1,ny
            do i = 2, nx
                n=i+nx*(j-1)+nx*ny*(k-1)
                abd(mu1,n)=1.
            enddo
            do i = 1, nx-1
                n=i+nx*(j-1)+nx*ny*(k-1)
                abd(ml1,n)=1.
            enddo
        enddo
    enddo
    !
    do k=1,nz
        do j = nx+1, nx*ny
            n=j+nx*ny*(k-1)
            abd(mu2,n)=1.
        enddo
        do j = 1, nx*ny-nx
            n=j+nx*ny*(k-1)
            abd(ml2,n)=1.
        enddo
    enddo
    !
    do j = ny*nx+1, ntot
        abd(mu3,j)=1.
    enddo
    do j = 1, ntot-nx*ny
        abd(ml3,j)=1.
    enddo
    !
    !      set boundary conditions : symetrix matrix
    !         so don't have to worry about boundary conditions
    !
    !        abd(ml,1)=0.
    !        abd(mu,ntot)=0.
    !
    !      find div b and load into matrix
    !
    poten=0.
    do k=2,nz-1
        az=grd_zmin(m)+rz*(k-1)
        do j=2,ny-1
            ay=grd_ymin(m)+ry*(j-1)
            do i=2,nx-1
                ax=grd_xmin(m)+rx*(i-1)
                ar=sqrt(ax**2+ay**2+az**2)
                !
                n=i+nx*(j-1)+nx*ny*(k-1)
                divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m)) &
                    +(by(i,j+1,k,m)-by(i,j-1,k,m)) &
                    +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
                !
                !          if(m.gt.1)then
                !            if (ar.gt.rx*rearth)then
                !             b(n)=divb
                !             poten(i,j,k)=abs(divb)    ! needed only for plotting
                !            endif
                !          else
                !             b(n)=divb
                !             poten(i,j,k)=abs(divb)    ! needed only for plotting
                !          endif
                !
                b(n)=divb
                poten(i,j,k)=abs(divb)
            enddo
        enddo
    enddo
    !
    !     plot div b error : poten
    !
    write(6,*)'plot divb in',m,nband,ntot,poten(10,10,10)
    write(6,*)xmin,xmax,ymin,ymax,zmin,zmax
    dbx=0.
    dby=0.
    dbz=0.
    call conflow(poten,dbx,dby,dbz,nx,ny,nz,1,m,m, &
        xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
        ut,'divb in',3,11,1,2.0, &
        tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
    !
    call sparse(abd,nband,ntot,b,x,ipvt,rsq, &
        g,h,xi,xj,.false.)
    !
    scale=rx**2
    do k=1,nz
        do j=1,ny
            do i=1,nx
                n=i+nx*(j-1)+nx*ny*(k-1)
                poten(i,j,k)=scale*x(n)
            enddo
        enddo
    enddo
    !
    !    calculate  perturbed magnetic field
    !
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                dbx(i,j,k)=(poten(i+1,j,k)-poten(i-1,j,k))/(2.*rx)
                dby(i,j,k)=(poten(i,j+1,k)-poten(i,j-1,k))/(2.*ry)
                dbz(i,j,k)=(poten(i,j,k+1)-poten(i,j,k-1))/(2.*rz)
            enddo
        enddo
    enddo
    !
    !     recalculate divb to see if this is correct
    !
    poten=0.
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                n=i+nx*(j-1)+nx*ny*(k-1)
                divb=( (dbx(i+1,j,k)-dbx(i-1,j,k)) &
                    +(dby(i,j+1,k)-dby(i,j-1,k)) &
                    +(dbz(i,j,k+1)-dbz(i,j,k-1)) )/(2.*rx)
                    poten(i,j,k)=abs(divb)
            enddo
        enddo
    enddo
    !
    !      plot the calculated divb error
    !
    write(6,*)'plot divb out',m,nband,ntot,poten(10,10,10)
    h=0.
    xi=0.
    xj=0.
    call conflow(poten,h,xi,xj,nx,ny,nz,1,m,m, &
        xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
        ut,'divb out',3,11,1,2.0, &
        tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
    !
    !     modified magnetic field
    !
    
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                bx(i,j,k,m)=bx(i,j,k,m)-dbx(i,j,k)
                by(i,j,k,m)=by(i,j,k,m)-dby(i,j,k)
                bz(i,j,k,m)=bz(i,j,k,m)-dbz(i,j,k)
            enddo
        enddo
    enddo
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                n=i+nx*(j-1)+nx*ny*(k-1)
                    divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m)) &
                    +(by(i,j+1,k,m)-by(i,j-1,k,m)) &
                    +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
                poten(i,j,k)=abs(divb)    ! needed only for plotting
            enddo
        enddo
    enddo
    !
    !     plot div b error : poten
    !
    h=0.
    xi=0.
    xj=0.
    write(6,*)'plot divb final',m,nband,ntot,poten(10,10,10)
    call conflow(poten,h,xi,xj,nx,ny,nz,1,m,m, &
        xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
        ut,'divb fin',3,11,1,2.0, &
        tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
    !
    return
end

!
!      ***************************************
!
recursive subroutine sparse(abd,nband,n,b,x,ipvt,rsq, &
    g,h,xi,xj,recurse,rirst)
    !
    !     sparse matrix solver:
    !       abd is the banded matrix
    !       nband is the number of bands
    !       ivpt is the position of the bands relative
    !              to the diagonal
    !
    !       n is the order of the original matrix
    !       b(n) the real vector for the right hand side
    !       x(n) is the solution vector
    !
    !       rsq is the sum of the squares of the components
    !         of the residual vector a.x -b
    !       if this is not small then the matrix is numerically
    !       singular and the solution represents a least
    !        squares best approximation
    !
    !    nmax is the anticpated maximum value of n
    !       should be reset if you want a big system
    !     eps is an error limit onthe solution
    !
    parameter(eps=1.e-4)
    !
    real abd(nband,n),b(n),x(n)
    integer ipvt(nband)
    real g(n),h(n),xi(n),xj(n)
    logical,optional :: recurse
    integer, optional :: rirst
    
    !
    !      criterion for sum-squared residuals
    !       and number of restarts attempted internally
    !
    eps2=n*eps**2
    
    if(recurse) then
        irst = rirst
    else
        irst=0
    endif
    
    !
    irst=irst+1
    call asub(abd,nband,n,ipvt,x,xi)
    rp=0.
    !
    !      add the magntiude of the right side
    !       and find the residue
    !
    bsq=0.
    do j=1,n
        bsq=bsq+b(j)**2
        xi(j)=xi(j)-b(j)
        rp=rp+xi(j)**2
    enddo
    !
    call atsub(abd,nband,n,ipvt,xi,g)
    do j=1,n
        g(j)=-g(j)
        h(j)=g(j)
    enddo
    !
    !      main iterative loop
    !
    max_iter=10*n
    do iter=1,max_iter
        call asub(abd,nband,n,ipvt,h,xi)
        !
        !     calculate the gradient
        !
        anum=0.
        aden=0.
        !
        do j=1,n
            anum=anum+g(j)*h(j)
            aden=aden+xi(j)**2
        enddo
        if(aden.eq.0.)pause 'very singular matrix'
        !
        anum=anum/aden
        do j=1,n
            xi(j)=x(j)
            x(j)=x(j)+anum*h(j)
        enddo
        !
        call asub(abd,nband,n,ipvt,x,xj)
        !
        rsq=0.
        do j=1,n
            xj(j)=xj(j)-b(j)
            rsq=rsq+xj(j)**2
        enddo
        !
        !      test for convergence and exit if okay
        !
        if(rsq.eq.rp.or.rsq.le.bsq*eps2)return
        !
        !     test to see whether solution is improving
        !     and restart if necessary
        !
        if(rsq.gt.rp)then
            do j=1,n
                x(j)=xi(j)
            enddo
            !
            !      return if too many restarts - hitting roundoff error
            !
            if(irst.ge.3)return
            !
            call sparse(abd,nband,n,b,x,ipvt,rsq, &
                g,h,xi,xj,.true.,irst)
        endif
        !
        !      compute gradient for next iteration
        !
        rp=rsq
        call atsub(abd,nband,n,ipvt,xj,xi)
        gg=0.
        dgg=0.
        do j=1,n
            gg=gg+g(j)**2
            dgg=dgg+(xi(j)+g(j))*xi(j)
        enddo
        !
        !     test to see if you have a solution and return if okay
        !
        if(gg.eq.0.)return
        !
        gam=dgg/gg
        do j=1,n
            g(j)=-xi(j)
            h(j)=g(j)+gam*h(j)
        enddo
        !
    enddo
    !
    !     never found solution if you get here
    !
    pause 'too many iterations'
    return
end

!
!     **************************
!
subroutine asub(abd,nband,n,ipvt,x,xi)
    !
    !    calculates a.x taking advantage of the banded
    !       structure of a
    !
    real abd(nband,n),x(n),xi(n)
    integer ipvt(nband)
    
    !
    !     initialize product
    !
    do i=1,n
        xi(i)=0.
        !
        do m=1,nband
            ii=i+ipvt(m)
            if((ii.ge.1).and.(ii.le.n)) xi(i)=xi(i)+abd(m,ii)*x(ii)
        enddo
    enddo
    !
    return
end

!
!     **************************
!
subroutine atsub(abd,nband,n,ipvt,x,xi)
    !
    !    beware: symmetrical matrix assumed here
    !    calculates transpose(a).x
    !         taking advantage of the banded
    !          structure of a
    !
    real abd(nband,n),x(n),xi(n)
    integer ipvt(nband)
    !
    !     initialize product
    !
    do i=1,n
        xi(i)=0.
        !
        do m=1,nband
            ii=i+ipvt(m)
            if((ii.ge.1).and.(ii.le.n))xi(i)=xi(i)+abd(m,ii)*x(ii)
        enddo
    enddo
    !
    !
    return
end