!
!	This file contains three subroutines:
!	divb_correct
!	divb_correct_n
!	divb_correct_tst
!
subroutine divb_correct(bx,by,bz,dbx,dby,dbz,poten, &
    b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,n_grids,box, &
    xspac)
    !
    !     this program solves the matrix equation
    !            a.x=b
    !      using using conjugate gradient method from
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
    !     store only non-zero bands of matrix a in
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
    real bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids)
    real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz), &
        poten(nx,ny,nz)
    real xspac(n_grids)
    !
    real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
    integer ipvt(nband)
    integer box
    !
    real, allocatable, dimension(:,:) :: abd
    allocate(abd(7,nx*ny*nz))
    !
    !write(*,*)'calling divb',nband,ntot,nx,ny,nz,n_grids,box
    !
    rx=xspac(box)
    ry=xspac(box)
    rz=xspac(box)
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
                divb=( (bx(i+1,j,k,box)-bx(i-1,j,k,box)) &
                    +(by(i,j+1,k,box)-by(i,j-1,k,box)) &
                    +(bz(i,j,k+1,box)-bz(i,j,k-1,box)) )/(2.*rx)
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
    !    calculate perturbed magnetic field
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
                bx(i,j,k,box)=bx(i,j,k,box)-dbx(i,j,k)
                by(i,j,k,box)=by(i,j,k,box)-dby(i,j,k)
                bz(i,j,k,box)=bz(i,j,k,box)-dbz(i,j,k)
            enddo
        enddo
    enddo
    !
    if(allocated(abd)) deallocate(abd)
    return
end
!
!
!	****************************************
!
!
subroutine divb_correct_n(bx,by,bz,dbx,dby,dbz,poten, &
    b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,n_grids,box, &
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
    real bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids)
    real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz), &
        poten(nx,ny,nz)
    real xspac(n_grids)
    !
    real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
    integer ipvt(nband)
    integer box
    !
    !      write(*,*)'calling divb',nband,ntot,nx,ny,nz,n_grids,box
    !
    rx=xspac(box)
    ry=xspac(box)
    rz=xspac(box)
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
                divb=( (bx(i+1,j,k,box)-bx(i-1,j,k,box)) &
                    +(by(i,j+1,k,box)-by(i,j-1,k,box)) &
                    +(bz(i,j,k+1,box)-bz(i,j,k-1,box)) )/(2.*rx)
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
    !    calculate perturbed magnetic field
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
                bx(i,j,k,box)=bx(i,j,k,box)-dbx(i,j,k)
                by(i,j,k,box)=by(i,j,k,box)-dby(i,j,k)
                bz(i,j,k,box)=bz(i,j,k,box)-dbz(i,j,k)
            enddo
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
subroutine divb_correct_tst(bx,by,bz,dbx,dby,dbz,poten, &
    b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,n_grids,box, &
    tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,xspac,r_inner)
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
    real bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids)
    real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz), &
        poten(nx,ny,nz)
    real grd_xmin(n_grids),grd_xmax(n_grids),grd_ymin(n_grids),grd_ymax(n_grids), &
        grd_zmin(n_grids),grd_zmax(n_grids),xspac(n_grids)
    real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz), &
        tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2)
    !
    real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
    integer ipvt(nband)
    integer box
    !
    write(*,*)'Calling divb: nband, ntot, nx, ny, nz, n_grids, box'
	write(*,*) nband,ntot,nx,ny,nz,n_grids,box
    !
    rx=xspac(box)
    ry=xspac(box)
    rz=xspac(box)
    !
    ymin=grd_ymin(box)+ry
    ymax=grd_ymax(box)-ry
    zmin=grd_zmin(box)+rz
    zmax=grd_zmax(box)-rz
    xmin=grd_xmin(box)+rx
    xmax=grd_xmax(box)-rx
    !
    xcut=(xmin+xmax)/2.
    !
	write(*,*) 'xmin, xmax, ymin, ymax, zmin, zmax'
    write(*,*) xmin,xmax,ymin,ymax,zmin,zmax
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
        az=grd_zmin(box)+rz*(k-1)
        do j=2,ny-1
            ay=grd_ymin(box)+ry*(j-1)
            do i=2,nx-1
                ax=grd_xmin(box)+rx*(i-1)
                ar=sqrt(ax**2+ay**2+az**2)
                !
                n=i+nx*(j-1)+nx*ny*(k-1)
                divb=( (bx(i+1,j,k,box)-bx(i-1,j,k,box)) &
                    +(by(i,j+1,k,box)-by(i,j-1,k,box)) &
                    +(bz(i,j,k+1,box)-bz(i,j,k-1,box)) )/(2.*rx)
                !
                !          if(box.gt.1)then
                !            if (ar.gt.rx*r_inner)then
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
    write(*,*)'plot divb in: box, nband, ntot, poten'
	write(*,*) box,nband,ntot,poten(10,10,10)
    write(*,*) 'xmin, xmax, ymin, ymax, zmin, zmax'
    write(*,*) xmin,xmax,ymin,ymax,zmin,zmax
    dbx=0.
    dby=0.
    dbz=0.
    call conflow(poten,dbx,dby,dbz,nx,ny,nz,1,box,box, &
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
    !    calculate perturbed magnetic field
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
    write(*,*) 'plot divb out: box, nband, ntot, poten'
	write(*,*) box,nband,ntot,poten(10,10,10)
    h=0.
    xi=0.
    xj=0.
    call conflow(poten,h,xi,xj,nx,ny,nz,1,box,box, &
        xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
        ut,'divb out',3,11,1,2.0, &
        tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
    !
    !     modified magnetic field
    !
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                bx(i,j,k,box)=bx(i,j,k,box)-dbx(i,j,k)
                by(i,j,k,box)=by(i,j,k,box)-dby(i,j,k)
                bz(i,j,k,box)=bz(i,j,k,box)-dbz(i,j,k)
            enddo
        enddo
    enddo
    do k=2,nz-1
        do j=2,ny-1
            do i=2,nx-1
                n=i+nx*(j-1)+nx*ny*(k-1)
                    divb=( (bx(i+1,j,k,box)-bx(i-1,j,k,box)) &
                    +(by(i,j+1,k,box)-by(i,j-1,k,box)) &
                    +(bz(i,j,k+1,box)-bz(i,j,k-1,box)) )/(2.*rx)
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
    write(*,*) 'plot divb final: box, nband, ntot, poten'
	write(*,*) box,nband,ntot,poten(10,10,10)
    call conflow(poten,h,xi,xj,nx,ny,nz,1,box,box, &
        xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
        ut,'divb fin',3,11,1,2.0, &
        tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
    !
    return
end
