subroutine set_rho(qrho,qpresx,qpresy,qpresz, &
    qpresxy,qpresxz,qpresyz,rmassq, &
    hrho,hpresx,hpresy,hpresz, &
    hpresxy,hpresxz,hpresyz,rmassh, &
    orho,opresx,opresy,opresz, &
    opresxy,opresxz,opresyz,rmasso, &
    epres,nx,ny,nz,ngrd,m,o_conc)
    !
    !    checks for minimum rho and negative pressure
    !     and resets value if necessary
    !
    dimension qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    !
    !     d_min is the minimum allowable density
    !
    q_min=2.e-4      ! min q aboslute density
    h_min=1.e-4      ! min h aboslute density
    o_min=2.e-5      ! min o aboslute density
    dstep=20.
    !
    !$omp  parallel do
    do  k=1,nz
        do  j=1,ny
            do  i=1,nx
                !
                tden=abs(qrho(i,j,k,m))/rmassq+abs(hrho(i,j,k,m))/rmassh &
                +abs(orho(i,j,k,m))/rmasso+1.e-6
                rhdens=abs(hrho(i,j,k,m))/rmassh/tden
                rodens=abs(orho(i,j,k,m))/rmasso/tden
                !
                if ((qrho(i,j,k,m).lt.q_min).or.(qpresx(i,j,k,m).lt.0.).or. &
                    (qpresy(i,j,k,m).lt.0.).or.(qpresz(i,j,k,m).lt.0.))then
                    !
                    qrho(i,j,k,m)=amax1(abs(qrho(i,j,k,m)),1.05*q_min)
                    apres=(abs(qpresx(i,j,k,m))+abs(qpresy(i,j,k,m))+ &
                    abs(qpresz(i,j,k,m)))/3.
                    qpresx(i,j,k,m)=apres
                    qpresy(i,j,k,m)=apres
                    qpresz(i,j,k,m)=apres
                    qpresxy(i,j,k,m)=0.
                    qpresxz(i,j,k,m)=0.
                    qpresyz(i,j,k,m)=0.
                endif
                !
                if ((hrho(i,j,k,m).lt.h_min).or.(hpresx(i,j,k,m).lt.0.).or. &
                    (hpresy(i,j,k,m).lt.0.).or.(hpresz(i,j,k,m).lt.0.))then
            
                    hrho(i,j,k,m)=amax1(abs(hrho(i,j,k,m)),1.05*h_min)
                    apres=(abs(hpresx(i,j,k,m))+abs(hpresy(i,j,k,m))+ &
                    abs(hpresz(i,j,k,m)))/3.
                    hpresx(i,j,k,m)=apres
                    hpresy(i,j,k,m)=apres
                    hpresz(i,j,k,m)=apres
                    hpresxy(i,j,k,m)=0.
                    hpresxz(i,j,k,m)=0.
                    hpresyz(i,j,k,m)=0.
                endif
                !
                if ((orho(i,j,k,m).lt.o_min).or.(opresx(i,j,k,m).lt.0.).or. &
                    (opresy(i,j,k,m).lt.0.).or.(opresz(i,j,k,m).lt.0.))then
                
                    orho(i,j,k,m)=amax1(abs(orho(i,j,k,m)),1.05*o_min)
                    apres=(abs(opresx(i,j,k,m))+abs(opresy(i,j,k,m))+ &
                    abs(opresz(i,j,k,m)))/3.
                    opresx(i,j,k,m)=apres
                    opresy(i,j,k,m)=apres
                    opresz(i,j,k,m)=apres
                    opresxy(i,j,k,m)=0.
                    opresxz(i,j,k,m)=0.
                    opresyz(i,j,k,m)=0.
                endif
                !
                epres(i,j,k,m)=abs(epres(i,j,k,m))
                !
            enddo
        enddo
    enddo
    !
    return
end
