subroutine qvset(a,v,n)
    !
    !     set an array v of total dimension n to value a
    !
    dimension v(1)
    !
    ! parallelizes loop. RW, oct. 23, 2002
    !$omp  parallel do
    do i=1,n
        v(i)=a
    enddo
    return
end
