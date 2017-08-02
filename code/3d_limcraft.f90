subroutine limcraft(xcraft,ncraft,re_equiv,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !
    !      tests to see whether spacecraft is in the system and
    !           resets their position if not
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
        grd_ymin(ngrd),grd_ymax(ngrd), &
        grd_zmin(ngrd),grd_zmax(ngrd)
    dimension xcraft(4,ncraft)
    !
    abit=0.001
    m=ngrd
    do n=1,ncraft
        xcraft(1,n)=amax1(xcraft(1,n),(grd_xmin(m)+abit)*re_equiv)
        xcraft(1,n)=amin1(xcraft(1,n),(grd_xmax(m)-abit)*re_equiv)
        xcraft(2,n)=amax1(xcraft(2,n),(grd_ymin(m)+abit)*re_equiv)
        xcraft(2,n)=amin1(xcraft(2,n),(grd_ymax(m)-abit)*re_equiv)
        xcraft(3,n)=amax1(xcraft(3,n),(grd_zmin(m)+abit)*re_equiv)
        xcraft(3,n)=amin1(xcraft(3,n),(grd_zmax(m)-abit)*re_equiv)
    enddo
    !
    return
end
