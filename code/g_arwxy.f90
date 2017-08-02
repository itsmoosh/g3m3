subroutine arwxy(dx,dy,al,beta,adx1,ady1,adx2,ady2)
    !
    !     calculate relative positions of two points to draw arrows
    !        input: dx,dy  -- vector
    !           al - length  of arrow
    !           beta - angle of arrow
    !     output: (adx1,ady1), (adx2,ady2)  -- two relative positions
    !
    data pih,pi,pi32,pi2/1.57079,3.14159,4.71239,6.28318/
    !
    !       determine angle of the vector
    !
    if(abs(dx).le.0.00001) then
        theta=pih
        if(dy.lt.0.)  theta=pi32
    else if (abs(dy).lt.0.00001) then
        theta=0.
        if(dx.lt.0.) theta=pi
    else
        theta=atan(abs(dy/dx))
        if(dx.lt.0..and.dy.gt.0.) then
            theta=pi-theta
        else if(dx.lt.0..and.dy.lt.0.) then
            theta=pi+theta
        else if(dx.gt.0..and.dy.lt.0.) then
            theta=pi2-theta
        endif
    endif
    !
    !      calculate the relative position of two points
    !
    alfa1=theta-beta
    alfa2=pih-theta-beta
    adx1=-al*cos(alfa1)
    ady1=-al*sin(alfa1)
    adx2=-al*sin(alfa2)
    ady2=-al*cos(alfa2)
    !
    return
end
