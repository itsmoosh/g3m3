subroutine limcraft(xcraft,ncraft,re_equiv,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !      tests to see whether spacecraft is in the system and
    !           resets their position if not
    !
    integer box
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
        grd_ymin(n_grids),grd_ymax(n_grids), &
        grd_zmin(n_grids),grd_zmax(n_grids)
    dimension xcraft(4,ncraft)
    !
    abit=0.001
    box=n_grids
    do n=1,ncraft
        xcraft(1,n)=amax1(xcraft(1,n),(grd_xmin(box)+abit)*re_equiv)
        xcraft(1,n)=amin1(xcraft(1,n),(grd_xmax(box)-abit)*re_equiv)
        xcraft(2,n)=amax1(xcraft(2,n),(grd_ymin(box)+abit)*re_equiv)
        xcraft(2,n)=amin1(xcraft(2,n),(grd_ymax(box)-abit)*re_equiv)
        xcraft(3,n)=amax1(xcraft(3,n),(grd_zmin(box)+abit)*re_equiv)
        xcraft(3,n)=amin1(xcraft(3,n),(grd_zmax(box)-abit)*re_equiv)
    enddo
    !
    return
end
