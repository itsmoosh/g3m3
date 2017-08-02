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
    !     eps is an error limit on the solution
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
