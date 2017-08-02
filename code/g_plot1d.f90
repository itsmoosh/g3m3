subroutine plot1d(xary,yary,npt,time,label)
    !
    dimension xary(npt),yary(npt),xmax(2),ymax(2)
    character*12 label
    character*20 title
    !
    fmax=yary(1)
    fmin=fmax
    do i=2,npt
        fmax=amax1(yary(i),fmax)
        fmin=amin1(yary(i),fmin)
    enddo
    !
    if((fmax.gt.0.0).and.(fmin.gt.0.0))fmin=0.0
    if((fmax.lt.0.0).and.(fmin.lt.0.0))fmax=0.0
    !
    xmax(2)=xary(npt)
    xmax(1)=xary(1)
    ymax(2)=fmax
    ymax(1)=fmin
    !
    call frame
    call gselnt(0)
    call gsplci(1)
    call agseti('frame.',2)
    call agseti('set.',-1)
    write(title,'(f7.3)')time
    title='t= '//title
    call wtstr(.75,.97,title,1,0,0)
    call ezxy(xmax,ymax,2,label)
    call curve(xary,yary,npt)
    call agseti('set.',1)
    call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
    !
    return
end
