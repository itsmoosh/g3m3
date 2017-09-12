      subroutine divb_correct(bx,by,bz,dbx,dby,dbz,poten,
     +           b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,ngrd,m,
     +           xspac)
c
c     this program solves the matrix equation
c            A.x=b
c      using using Conjuage Gradient Method from
c          Numerical recipes
c     3: D solution of Poisson's eqn :
c          ntot=nx*ny*nz
c
c       mapping is such that for position i,j,k:
c          n =i+nx*(j-1) +nx*ny*(k-1)
c
c     original solutions done on small grid
c     parameter (ntot=893101,nx=121,ny=121,nz=61,nband=7)
c
c     store only non-zero bnads of matrix A in
c       array abd, and the position of these bands is
c       stored in ipvt
c       abd gives values of bands in matrix
c       b is data vector
c       x is the solution we are seeking
c       xi,xj,g,h are work arrays needed for inversion
c
c           nband = 7 for 3d and 5 for 2d
c
c     size of array abd hard wired to be nband,nx*ny*nz=ntot
c
c
c     common /space/abd(7,893101)
      real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
      real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz),
     +     poten(nx,ny,nz)
      real xspac(ngrd)
c
      real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
      integer ipvt(nband)
c
      REAL, ALLOCATABLE, DIMENSION(:,:) :: abd
      allocate(abd(7,nx*ny*nz))
c
       write(6,*)'calling divb',nband,ntot,nx,ny,nz,ngrd,m
c
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
c
c
c      Make the banded matrix
c
c            highest band position :  mu2, mu1
c            lowest band position :  ml2, ml1
c            diagonal band position :  m0
c
c             Note in matrix notation top row is 1
c                                    bottom row is nband
c
         m0=4
         ml1=5
         ml2=6
         ml3=7
c
         mu1=3
         mu2=2
         mu3=1
c
c
c      position of the bands
c
      ipvt(m0)=0
      ipvt(ml1)=-1
      ipvt(ml2)=-nx
      ipvt(ml3)=-nx*ny
      ipvt(mu1)=1
      ipvt(mu2)=nx
      ipvt(mu3)=nx*ny
c
c      initialize elements of banded matrix
c         for diagonal or element N requires
c            mu2 at N+nx  -> limits nx+1,ntot
c            mu1 at N+1 -> limits 2,ntot
c            ml2 at N-nx -> limits 1,ntot-nx
c            ml1 at N-1 -> limits 1,ntot-1
c
c
c       initial bands to zero
c
c
        abd=0.
c
c       do centerband
c        
        do  n = 1, ntot
         abd(m0,n)=-6.
         b(n)=0.
        enddo
c
c     do off diagonal elements
c
      do 30 k=1,nz
      do 30 j = 1,ny
        do 24 i = 2, nx
          n=i+nx*(j-1)+nx*ny*(k-1)
         abd(mu1,n)=1.
 24     continue
        do 26 i = 1, nx-1
          n=i+nx*(j-1)+nx*ny*(k-1)
         abd(ml1,n)=1.
 26     continue
 30   continue
c
      do 40 k=1,nz
        do 34 j = nx+1, nx*ny
         n=j+nx*ny*(k-1)
         abd(mu2,n)=1.
 34     continue
        do 36 j = 1, nx*ny-nx
         n=j+nx*ny*(k-1)
         abd(ml2,n)=1.
 36     continue
 40   continue
c
        do 50 j = ny*nx+1, ntot
         abd(mu3,j)=1.
 50     continue
        do 55 j = 1, ntot-nx*ny
         abd(ml3,j)=1.
 55   continue
c
c      set boundary conditions : symetrix matrix
c         so don't have to worry about boundary conditions
c
c        abd(ml,1)=0.
c        abd(mu,ntot)=0.
c
c      find div B and load into matrix
c   
       do k=2,nz-1
        do j=2,ny-1
         do i=2,nx-1
c
           n=i+nx*(j-1)+nx*ny*(k-1)
           divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m))
     +           +(by(i,j+1,k,m)-by(i,j-1,k,m))
     +           +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
c
              b(n)=divb
         enddo
       enddo
       enddo
c
      call sparse(abd,nband,ntot,b,x,ipvt,rsq,
     +           g,h,xi,xj)
c
c     convert potential onto grid
c
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
c
c    calculate  perturbed magnetic field
c 
c
      do k=2,nz-1
      do j=2,ny-1
       do i=2,nx-1
        dbx(i,j,k)=(poten(i+1,j,k)-poten(i-1,j,k))/(2.*rx)
        dby(i,j,k)=(poten(i,j+1,k)-poten(i,j-1,k))/(2.*ry)
        dbz(i,j,k)=(poten(i,j,k+1)-poten(i,j,k-1))/(2.*rz)
       enddo
      enddo
      enddo
c
c
c     modified magnetic field
c

        do k=2,nz-1
        do j=2,ny-1
        do i=2,nx-1
          bx(i,j,k,m)=bx(i,j,k,m)-dbx(i,j,k)
          by(i,j,k,m)=by(i,j,k,m)-dby(i,j,k)
          bz(i,j,k,m)=bz(i,j,k,m)-dbz(i,j,k)
        enddo
        enddo
        enddo
c
      deallocate(abd)
c
      return
      end
c
c    **********************************************************
c
      subroutine divb_correct_n(bx,by,bz,dbx,dby,dbz,poten,
     +           b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,ngrd,m,
     +           xspac)
c
c     this program solves the matrix equation
c            A.x=b
c      using using Conjuage Gradient Method from
c          Numerical recipes
c     3: D solution of Poisson's eqn :
c          ntot=nx*ny*nz
c
c       mapping is such that for position i,j,k:
c          n =i+nx*(j-1) +nx*ny*(k-1)
c
c     original solutions done on small grid
c     parameter (ntot=893101,nx=121,ny=121,nz=61,nband=7)
c
c     store only non-zero bnads of matrix A in
c       array abd, and the position of these bands is
c       stored in ipvt
c       abd gives values of bands in matrix
c       b is data vector
c       x is the solution we are seeking
c       xi,xj,g,h are work arrays needed for inversion
c
c           nband = 7 for 3d and 5 for 2d
c
c     size of array abd hard wired to be nband,nx*ny*nz=ntot
c
c
      common /space_n/abd(7,117649)
      real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
      real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz),
     +     poten(nx,ny,nz)
      real xspac(ngrd)
c
      real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
      integer ipvt(nband)
c
c      write(6,*)'calling divb',nband,ntot,nx,ny,nz,ngrd,m
c
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
c
c
c      Make the banded matrix
c
c            highest band position :  mu2, mu1
c            lowest band position :  ml2, ml1
c            diagonal band position :  m0
c
c             Note in matrix notation top row is 1
c                                    bottom row is nband
c
         m0=4
         ml1=5
         ml2=6
         ml3=7
c
         mu1=3
         mu2=2
         mu3=1
c
c
c      position of the bands
c
      ipvt(m0)=0
      ipvt(ml1)=-1
      ipvt(ml2)=-nx
      ipvt(ml3)=-nx*ny
      ipvt(mu1)=1
      ipvt(mu2)=nx
      ipvt(mu3)=nx*ny
c
c      initialize elements of banded matrix
c         for diagonal or element N requires
c            mu2 at N+nx  -> limits nx+1,ntot
c            mu1 at N+1 -> limits 2,ntot
c            ml2 at N-nx -> limits 1,ntot-nx
c            ml1 at N-1 -> limits 1,ntot-1
c
c
c       initial bands to zero
c
c
        abd=0.
c
c       do centerband
c        
        do  n = 1, ntot
         abd(m0,n)=-6.
         b(n)=0.
        enddo
c
c     do off diagonal elements
c
      do 30 k=1,nz
      do 30 j = 1,ny
        do 24 i = 2, nx
          n=i+nx*(j-1)+nx*ny*(k-1)
         abd(mu1,n)=1.
 24     continue
        do 26 i = 1, nx-1
          n=i+nx*(j-1)+nx*ny*(k-1)
         abd(ml1,n)=1.
 26     continue
 30   continue
c
      do 40 k=1,nz
        do 34 j = nx+1, nx*ny
         n=j+nx*ny*(k-1)
         abd(mu2,n)=1.
 34     continue
        do 36 j = 1, nx*ny-nx
         n=j+nx*ny*(k-1)
         abd(ml2,n)=1.
 36     continue
 40   continue
c
        do 50 j = ny*nx+1, ntot
         abd(mu3,j)=1.
 50     continue
        do 55 j = 1, ntot-nx*ny
         abd(ml3,j)=1.
 55   continue
c
c      set boundary conditions : symetrix matrix
c         so don't have to worry about boundary conditions
c
c        abd(ml,1)=0.
c        abd(mu,ntot)=0.
c
c      find div B and load into matrix
c   
       do k=2,nz-1
        do j=2,ny-1
         do i=2,nx-1
c
           n=i+nx*(j-1)+nx*ny*(k-1)
           n=i+nx*(j-1)+nx*ny*(k-1)
           divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m))
     +           +(by(i,j+1,k,m)-by(i,j-1,k,m))
     +           +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
c
              b(n)=divb
         enddo
       enddo
       enddo
c
      call sparse(abd,nband,ntot,b,x,ipvt,rsq,
     +           g,h,xi,xj)
c
c     convert potential onto grid
c
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
c
c    calculate  perturbed magnetic field
c 
c
      do k=2,nz-1
      do j=2,ny-1
       do i=2,nx-1
        dbx(i,j,k)=(poten(i+1,j,k)-poten(i-1,j,k))/(2.*rx)
        dby(i,j,k)=(poten(i,j+1,k)-poten(i,j-1,k))/(2.*ry)
        dbz(i,j,k)=(poten(i,j,k+1)-poten(i,j,k-1))/(2.*rz)
       enddo
      enddo
      enddo
c
c
c     modified magnetic field
c

        do k=2,nz-1
        do j=2,ny-1
        do i=2,nx-1
          bx(i,j,k,m)=bx(i,j,k,m)-dbx(i,j,k)
          by(i,j,k,m)=by(i,j,k,m)-dby(i,j,k)
          bz(i,j,k,m)=bz(i,j,k,m)-dbz(i,j,k)
        enddo
        enddo
        enddo
c
      return
      end
c
c     *****************************************************
c
      subroutine divb_correct_tst(bx,by,bz,dbx,dby,dbz,poten,
     +           b,x,g,h,xi,xj,nband,ntot,nx,ny,nz,ngrd,m,
     +        tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, 
     +        grd_xmin,grd_xmax,grd_ymin,grd_ymax,
     +        grd_zmin,grd_zmax,xspac,rearth)
c
c     this program solves the matrix equation
c            A.x=b
c      using using Conjuage Gradient Method from
c          Numerical recipes
c     3: D solution of Poisson's eqn :
c          ntot=nx*ny*nz
c
c       mapping is such that for position i,j,k:
c          n =i+nx*(j-1) +nx*ny*(k-1)
c
c     original solutions done on small grid
c     parameter (ntot=893101,nx=121,ny=121,nz=61,nband=7)
c
c     store only non-zero bnads of matrix A in
c       array abd, and the position of these bands is
c       stored in ipvt
c       abd gives values of bands in matrix
c       b is data vector
c       x is the solution we are seeking
c       xi,xj,g,h are work arrays needed for inversion
c
c           nband = 7 for 3d and 5 for 2d
c
c     size of array abd hard wired to be nband,nx*ny*nz=ntot
c
c
      common /space/abd(7,893101)
      real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd)
      real dbx(nx,ny,nz),dby(nx,ny,nz),dbz(nx,ny,nz),
     +     poten(nx,ny,nz)
      real grd_xmin(ngrd),grd_xmax(ngrd),grd_ymin(ngrd),grd_ymax(ngrd),
     +     grd_zmin(ngrd),grd_zmax(ngrd),xspac(ngrd)
      real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz),
     +     tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2)
c
      real b(ntot),x(ntot),g(ntot),h(ntot),xi(ntot),xj(ntot)
      integer ipvt(nband)
c
       write(6,*)'calling divb',nband,ntot,nx,ny,nz,ngrd,m
c
       rx=xspac(m)
       ry=xspac(m)
       rz=xspac(m)
c
       ymin=grd_ymin(m)+ry
       ymax=grd_ymax(m)-ry
       zmin=grd_zmin(m)+rz
       zmax=grd_zmax(m)-rz
       xmin=grd_xmin(m)+rx
       xmax=grd_xmax(m)-rx
c
       xcut=(xmin+xmax)/2.
c
       write(6,*)xmin,xmax,ymin,ymax,zmin,zmax
c
c      Make the banded matrix
c
c            highest band position :  mu2, mu1
c            lowest band position :  ml2, ml1
c            diagonal band position :  m0
c
c             Note in matrix notation top row is 1
c                                    bottom row is nband
c
         m0=4
         ml1=5
         ml2=6
         ml3=7
c
         mu1=3
         mu2=2
         mu3=1
c
c
c      position of the bands
c
      ipvt(m0)=0
      ipvt(ml1)=-1
      ipvt(ml2)=-nx
      ipvt(ml3)=-nx*ny
      ipvt(mu1)=1
      ipvt(mu2)=nx
      ipvt(mu3)=nx*ny
c
c      initialize elements of banded matrix
c         for diagonal or element N requires
c            mu2 at N+nx  -> limits nx+1,ntot
c            mu1 at N+1 -> limits 2,ntot
c            ml2 at N-nx -> limits 1,ntot-nx
c            ml1 at N-1 -> limits 1,ntot-1
c
c
c       initial bands to zero
c
c       do  mm=1,nband
c       do  n = 1, ntot
c        abd(mm,n)=0.
c       enddo
c       enddo
c
        abd=0.
c
c       do centerband
c        
        do  n = 1, ntot
         abd(m0,n)=-6.
         b(n)=0.
        enddo
c
c     do off diagonal elements
c
      do 30 k=1,nz
      do 30 j = 1,ny
        do 24 i = 2, nx
          n=i+nx*(j-1)+nx*ny*(k-1)
         abd(mu1,n)=1.
 24     continue
        do 26 i = 1, nx-1
          n=i+nx*(j-1)+nx*ny*(k-1)
         abd(ml1,n)=1.
 26     continue
 30   continue
c
      do 40 k=1,nz
        do 34 j = nx+1, nx*ny
         n=j+nx*ny*(k-1)
         abd(mu2,n)=1.
 34     continue
        do 36 j = 1, nx*ny-nx
         n=j+nx*ny*(k-1)
         abd(ml2,n)=1.
 36     continue
 40   continue
c
        do 50 j = ny*nx+1, ntot
         abd(mu3,j)=1.
 50     continue
        do 55 j = 1, ntot-nx*ny
         abd(ml3,j)=1.
 55   continue
c
c      set boundary conditions : symetrix matrix
c         so don't have to worry about boundary conditions
c
c        abd(ml,1)=0.
c        abd(mu,ntot)=0.
c
c      find div B and load into matrix
c   
       poten=0.
       do k=2,nz-1
        az=grd_zmin(m)+rz*(k-1)
        do j=2,ny-1
         ay=grd_ymin(m)+ry*(j-1)
         do i=2,nx-1
           ax=grd_xmin(m)+rx*(i-1)
           ar=sqrt(ax**2+ay**2+az**2)
c
           n=i+nx*(j-1)+nx*ny*(k-1)
           divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m))
     +           +(by(i,j+1,k,m)-by(i,j-1,k,m))
     +           +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
c
c          if(m.gt.1)then
c            if (ar.gt.rx*rearth)then
c             b(n)=divb
c             poten(i,j,k)=abs(divb)    ! needed only for plotting
c            endif
c          else
c             b(n)=divb
c             poten(i,j,k)=abs(divb)    ! needed only for plotting
c          endif
c
              b(n)=divb
              poten(i,j,k)=abs(divb) 
         enddo
       enddo
       enddo
c
c     plot div B error : poten
c
      write(6,*)'plot divb IN',m,nband,ntot,poten(10,10,10)
      write(6,*)xmin,xmax,ymin,ymax,zmin,zmax
      dbx=0.
      dby=0.
      dbz=0.
      call conflow(poten,dbx,dby,dbz,nx,ny,nz,1,m,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +             ut,'divB in',3,11,1,2.0,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
      call sparse(abd,nband,ntot,b,x,ipvt,rsq,
     +           g,h,xi,xj)
c
      scale=rx**2
      do k=1,nz
      do j=1,ny
       do i=1,nx
        n=i+nx*(j-1)+nx*ny*(k-1)
        poten(i,j,k)=scale*x(n)
       enddo
      enddo
      enddo
c
c    calculate  perturbed magnetic field
c 
c
      do k=2,nz-1
      do j=2,ny-1
       do i=2,nx-1
        dbx(i,j,k)=(poten(i+1,j,k)-poten(i-1,j,k))/(2.*rx)
        dby(i,j,k)=(poten(i,j+1,k)-poten(i,j-1,k))/(2.*ry)
        dbz(i,j,k)=(poten(i,j,k+1)-poten(i,j,k-1))/(2.*rz)
       enddo
      enddo
      enddo
c
c     recalculate divb to see if this is correct
c
      poten=0.
c
       do k=2,nz-1
       do j=2,ny-1
         do i=2,nx-1
           n=i+nx*(j-1)+nx*ny*(k-1)
           divb=( (dbx(i+1,j,k)-dbx(i-1,j,k))
     +           +(dby(i,j+1,k)-dby(i,j-1,k))
     +           +(dbz(i,j,k+1)-dbz(i,j,k-1)) )/(2.*rx)
           poten(i,j,k)=abs(divb)
         enddo
       enddo
       enddo
c
c      plot the calculated divb error
c
      write(6,*)'plot divb OUT',m,nband,ntot,poten(10,10,10)
      h=0.
      xi=0.
      xj=0.
      call conflow(poten,h,xi,xj,nx,ny,nz,1,m,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +             ut,'divB out',3,11,1,2.0,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
c     modified magnetic field
c

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
           divb=( (bx(i+1,j,k,m)-bx(i-1,j,k,m))
     +           +(by(i,j+1,k,m)-by(i,j-1,k,m))
     +           +(bz(i,j,k+1,m)-bz(i,j,k-1,m)) )/(2.*rx)
           poten(i,j,k)=abs(divb)    ! needed only for plotting
         enddo
       enddo
       enddo
c
c     plot div B error : poten
c
      h=0.
      xi=0.
      xj=0.
      write(6,*)'plot divb FINAL',m,nband,ntot,poten(10,10,10)
      call conflow(poten,h,xi,xj,nx,ny,nz,1,m,m,
     +            xmin,xmax,ymin,ymax,zmin,zmax,xcut,
     +             ut,'divB fin',3,11,1,2.0,
     +             tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2)
c
      return
      end
c
c
c      ***************************************
c
      subroutine sparse(abd,nband,n,b,x,ipvt,rsq,
     +           g,h,xi,xj)
c
c     sparse matrix solver:
c       abd is the banded matrix
c       nband is the number of bands
c       ivpt is the position of the bands relative
c              to the diagonal
c
c       n is the order of the original matrix
c       b(n) the real vector for the right hand side
c       x(n) is the solution vector
c
c       rsq is the sum of the squares of the components
c         of the residual vector A.X -B
c       if this is not small then the matrix is numerically
c       singular and the solution represents a least
c        squares best approximation
c
c    nmax is the anticpated maximum value of n
c       SHOULD be reset if you want a big system
c     eps is an error limit onthe solution
c
      parameter(eps=1.e-4)
c
      real abd(nband,n),b(n),x(n)
      integer ipvt(nband)
      real g(n),h(n),xi(n),xj(n)
c
c      criterion for sum-squared residuals
c       and number of restarts attempted internally
c
      eps2=n*eps**2
      irst=0
c
    1 irst=irst+1
      call asub(abd,nband,n,ipvt,x,xi)
      rp=0.
c
c      add the magntiude of the right side
c       and find the residue
c
      bsq=0.
      do 11 j=1,n
        bsq=bsq+b(j)**2
        xi(j)=xi(j)-b(j)
        rp=rp+xi(j)**2
   11 continue
c
      call atsub(abd,nband,n,ipvt,xi,g)
      do 12 j=1,n
        g(j)=-g(j)
        h(j)=g(j)
   12 continue
c
c      main iterative loop
c
      max_iter=10*n
      do 100 iter=1,max_iter
      call asub(abd,nband,n,ipvt,h,xi)
c
c     calculate the gradient
c
      anum=0.
      aden=0.
c
      do 13 j=1,n
        anum=anum+g(j)*h(j)
        aden=aden+xi(j)**2
   13 continue
      if(aden.eq.0.)pause 'very singular matrix'
c
      anum=anum/aden
      do 14 j=1,n
         xi(j)=x(j)
         x(j)=x(j)+anum*h(j)
   14 continue
c
      call asub(abd,nband,n,ipvt,x,xj)
c
      rsq=0.
      do 15 j=1,n
        xj(j)=xj(j)-b(j)
        rsq=rsq+xj(j)**2
   15 continue
c
c      test for convergence and exit if okay
c
      if(rsq.eq.rp.or.rsq.le.bsq*eps2)return
c
c     test to see whether solution is improving
c     and restart if necessary
c
      if(rsq.gt.rp)then
        do 16 j=1,n
          x(j)=xi(j)
   16   continue
c
c      return if too many restarts - hitting roundoff error
c
        if(irst.ge.3)return
c
        goto 1
      endif
c
c      compute gradient for next iteration
c
      rp=rsq
      call atsub(abd,nband,n,ipvt,xj,xi)
      gg=0.
      dgg=0.
      do 17 j=1,n
        gg=gg+g(j)**2
       dgg=dgg+(xi(j)+g(j))*xi(j)
   17 continue
c
c     test to see if you have a solution and return if okay
c
      if(gg.eq.0.)return
c
      gam=dgg/gg
      do 18 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)+gam*h(j)
   18 continue
c
  100 continue
c
c     never found solution if you get here
c
      pause 'too many iterations'
      return
      end

c
c     **************************
c
      subroutine asub(abd,nband,n,ipvt,x,xi)
c
c    calculates A.X taking advantage of the banded
c       structure of A
c
      real abd(nband,n),x(n),xi(n)
      integer ipvt(nband)

c
c     initialize product
c
      do 20 i=1,n
        xi(i)=0.
c
        do 10 m=1,nband
          ii=i+ipvt(m)
          if((ii.lt.1).or.(ii.gt.n))goto 10
           xi(i)=xi(i)+abd(m,ii)*x(ii)
   10   continue
   20 continue
c
      return 
      end
c
c     **************************
c
      subroutine atsub(abd,nband,n,ipvt,x,xi)
c
c    BEWARE: symmetrical matrix assumed here
c    calculates transpose(A).X 
c         taking advantage of the banded
c          structure of A
c
      real abd(nband,n),x(n),xi(n)
      integer ipvt(nband)
c
c     initialize product
c
      do 20 i=1,n
        xi(i)=0.
c
        do 10 m=1,nband
          ii=i+ipvt(m)
          if((ii.lt.1).or.(ii.gt.n))goto 10
           xi(i)=xi(i)+abd(m,ii)*x(ii)
   10   continue
   20 continue
c
c
      return 
      end
c
