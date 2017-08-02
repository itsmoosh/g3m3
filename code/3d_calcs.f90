subroutine store_array(bx,bx0,nx,ny,nz,ngrd,m)
    !
    !
    !     calculates the total magnetic field from the perturbed and
    !        stationary magnetic field
    !
    dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                bx(i,j,k,m)=bx0(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine totfld(bx,bx0,btx,nx,ny,nz,ngrd,m)
    !
    !
    !     calculates the total magnetic field from the perturbed and
    !        stationary magnetic field
    !
    dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd),btx(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                btx(i,j,k)=bx0(i,j,k,m)+bx(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine fnd_fld(bx,btx,nx,ny,nz,ngrd,m)
    !
    !
    !     calculates the total dynamic magnetic field only
    !        dipole field added separately by by setting add_dip=.true.
    !
    dimension bx(nx,ny,nz,ngrd),btx(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                btx(i,j,k)=btx(i,j,k)+bx(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine fnd_vel(px,py,pz,rho,vx,vy,vz,nx,ny,nz,ngrd,m)
    !
    !     converts momentum into velocity for graphics
    !
    dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd), &
        pz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd), &
        vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                arho=amax1(rho(i,j,k,m),0.0001)
                vx(i,j,k)=px(i,j,k,m)/arho
                vy(i,j,k)=py(i,j,k,m)/arho
                vz(i,j,k)=pz(i,j,k,m)/arho
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
    opx,opy,opz,orho,vx,vy,vz,nx,ny,nz,ngrd,m, &
    rmassq,rmassh,rmasso)
    !
    !     converts momentum into velocity for graphics
    !
    dimension qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd), &
        qpz(nx,ny,nz,ngrd),qrho(nx,ny,nz,ngrd), &
        hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd), &
        hpz(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd), &
        opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd), &
        opz(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd), &
        vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                qden=(qrho(i,j,k,m)+0.000001)/rmassq
                hden=(hrho(i,j,k,m)+0.000001)/rmassh
                oden=(orho(i,j,k,m)+0.000001)/rmasso
                tden=qden+hden+oden
                !
                vx(i,j,k)=(qpx(i,j,k,m)/rmassq+hpx(i,j,k,m)/rmassh &
                    +opx(i,j,k,m)/rmasso)/tden
                vy(i,j,k)=(qpy(i,j,k,m)/rmassq+hpy(i,j,k,m)/rmassh &
                    +opy(i,j,k,m)/rmasso)/tden
                vz(i,j,k)=(qpz(i,j,k,m)/rmassq+hpz(i,j,k,m)/rmassh &
                    +opz(i,j,k,m)/rmasso)/tden
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
    opx,opy,opz,orho,curx,cury,curz, &
    evx,evy,evz,tvx,tvy,tvz, &
    nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,reynolds)
    !
    !     converts momentum into velocity for graphics
    !
    dimension qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd), &
        qpz(nx,ny,nz,ngrd),qrho(nx,ny,nz,ngrd), &
        hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd), &
        hpz(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd), &
        opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd), &
        opz(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd), &
        curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
        tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
        evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
    !
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                !      eden=amax1(qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m),
                !    +                 0.0001)
                qden=(qrho(i,j,k,m)+0.000001)/rmassq
                hden=(hrho(i,j,k,m)+0.000001)/rmassh
                oden=(orho(i,j,k,m)+0.000001)/rmasso
                tden=qden+hden+oden
                !
                !      keep sepearate the ion and current components
                !
                tvx(i,j,k)=(qpx(i,j,k,m)/rmassq+hpx(i,j,k,m)/rmassh &
                    +opx(i,j,k,m)/rmasso)/tden
                evx(i,j,k)= tvx(i,j,k) - curx(i,j,k)/tden/reynolds
                !
                tvy(i,j,k)=(qpy(i,j,k,m)/rmassq+hpy(i,j,k,m)/rmassh &
                    +opy(i,j,k,m)/rmasso)/tden
                evy(i,j,k)= tvy(i,j,k) - cury(i,j,k)/tden/reynolds
                !
                tvz(i,j,k)=(qpz(i,j,k,m)/rmassq+hpz(i,j,k,m)/rmassh &
                    +opz(i,j,k,m)/rmasso)/tden
                evz(i,j,k)= tvz(i,j,k)  -curz(i,j,k)/tden/reynolds
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine fnd_prss(px,py,pz,rho,erg,afld,nx,ny,nz,ngrd, &
    m,gamma1,igo)
    !
    !     now using pressure equation rather than an erg equation
    !
    dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd), &
        pz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd), &
        erg(nx,ny,nz,ngrd),afld(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                afld(i,j,k)=erg(i,j,k,m)
            enddo
        enddo
    enddo
    !
    !     if graphics refine plots
    !
    if(igo.le.0)return
    
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                afld(i,j,k)=sqrt(abs(afld(i,j,k)))
                afld(i,j,k)=amin1(afld(i,j,k),1.2)
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!
subroutine tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
    !
    !      initialize static magnetic field along entire grid
    !
    !
    dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
        btot(nx,ny,nz)
    !
    do k=1,nz
        do j=1,ny
            do i=1,nx
                atot=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2+bsz(i,j,k)**2)
                btot(i,j,k)=amax1(atot,1.e-5)
            enddo
        enddo
    enddo
    !
    return
end
!
!
!	******************************************
!
!

