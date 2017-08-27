subroutine set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
        orho,opresx,opresy,opresz,opx,opy,opz, &
        rmassq,rmassh,rmasso,epres, &
        qpresxy,qpresxz,qpresyz, &
        hpresxy,hpresxz,hpresyz, &
        opresxy,opresxz,opresyz, &
        rhop,svxp,svyp,svzp,svelx,spress, &
        ti_te,rho_frac,nx,ny,nz,ngrd)
    !
    !     Set IMF boundary conditions. This subroutine is 
	!		used only when spacecraft and craft_info are both true,
	!		as we are trying to use spacecraft data as physical
	!		inputs to set the IMF.
    !
    dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd),bxp(ny,nz), &
        by(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),byp(ny,nz), &
        bz(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd),bzp(ny,nz), &
        rhop(ny,nz),svxp(ny,nz),svyp(ny,nz),svzp(ny,nz), &
		!
        qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
        qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
        qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
        qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
        qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
		!
        hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
        hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
        hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
        hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
        hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
		!
        orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
        opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
        opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
        opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
        opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd)
    !
    m=ngrd
    !
    i=1
    do k=1,nz
        do j=1,ny
            bx(i,j,k,m)=bxp(j,k)-bx0(i,j,k,m)
            by(i,j,k,m)=byp(j,k)-by0(i,j,k,m)
            bz(i,j,k,m)=bzp(j,k)-bz0(i,j,k,m)
            !
            qrho(i,j,k,m)=rhop(j,k)
            qpx(i,j,k,m)=qrho(i,j,k,m)*svxp(j,k)
            qpy(i,j,k,m)=qrho(i,j,k,m)*svyp(j,k)
            qpz(i,j,k,m)=qrho(i,j,k,m)*svzp(j,k)
            qpresx(i,j,k,m)=0.5*spress
            qpresy(i,j,k,m)=qpresx(i,j,k,m)
            qpresz(i,j,k,m)=qpresx(i,j,k,m)
            qpresxy(i,j,k,m)=0.
            qpresxz(i,j,k,m)=0.
            qpresyz(i,j,k,m)=0.
            !
            hrho(i,j,k,m)=hfrac*frac_h*rhop(j,k)
            hpx(i,j,k,m)=hrho(i,j,k,m)*svxp(j,k)
            hpy(i,j,k,m)=hrho(i,j,k,m)*svyp(j,k)
            hpz(i,j,k,m)=hrho(i,j,k,m)*svzp(j,k)
            hpresx(i,j,k,m)=0.5*spress*hfrac
            hpresy(i,j,k,m)=hpresx(i,j,k,m)
            hpresz(i,j,k,m)=hpresx(i,j,k,m)
            hpresxy(i,j,k,m)=0.
            hpresxz(i,j,k,m)=0.
            hpresyz(i,j,k,m)=0.
            !
            orho(i,j,k,m)=ofrac*frac_o*rhop(j,k)
            opx(i,j,k,m)=orho(i,j,k,m)*svxp(j,k)
            opy(i,j,k,m)=orho(i,j,k,m)*svyp(j,k)
            opz(i,j,k,m)=orho(i,j,k,m)*svzp(j,k)
            opresx(i,j,k,m)=0.5*spress*ofrac
            opresy(i,j,k,m)=opresx(i,j,k,m)
            opresz(i,j,k,m)=opresx(i,j,k,m)
            opresxy(i,j,k,m)=0.
            opresxz(i,j,k,m)=0.
            opresyz(i,j,k,m)=0.
            !
            epres(i,j,k,m)=qpresx(i,j,k,m)/ti_te
            !
            avz=avz+svzp(j,k)
        enddo
    enddo
    !
    avz=avz/float(ny*nz)
    !      write(6,*)'set_imf avz',avz,ny,nz
    return
end
