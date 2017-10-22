subroutine set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
        orho,opresx,opresy,opresz,opx,opy,opz, &
        rmassq,rmassh,rmasso,epres, &
        qpresxy,qpresxz,qpresyz, &
        hpresxy,hpresxz,hpresyz, &
        opresxy,opresxz,opresyz, &
        rhop,svxp,svyp,svzp,svelx,spress, &
        ti_te,rho_frac,nx,ny,nz,n_grids)
    !
    !     Set IMF boundary conditions. This subroutine is 
	!		used only when spacecraft and craft_info are both true,
	!		as we are trying to use spacecraft data as physical
	!		inputs to set the IMF.
    !
    dimension bx(nx,ny,nz,n_grids),bx0(nx,ny,nz,n_grids),bxp(ny,nz), &
        by(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids),byp(ny,nz), &
        bz(nx,ny,nz,n_grids),bz0(nx,ny,nz,n_grids),bzp(ny,nz), &
        rhop(ny,nz),svxp(ny,nz),svyp(ny,nz),svzp(ny,nz), &
		!
        qrho(nx,ny,nz,n_grids),qpx(nx,ny,nz,n_grids), &
        qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
        qpresx(nx,ny,nz,n_grids),qpresy(nx,ny,nz,n_grids), &
        qpresz(nx,ny,nz,n_grids), qpresxy(nx,ny,nz,n_grids), &
        qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
		!
        hrho(nx,ny,nz,n_grids),hpx(nx,ny,nz,n_grids), &
        hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
        hpresx(nx,ny,nz,n_grids),hpresy(nx,ny,nz,n_grids), &
        hpresz(nx,ny,nz,n_grids), hpresxy(nx,ny,nz,n_grids), &
        hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
		!
        orho(nx,ny,nz,n_grids),opx(nx,ny,nz,n_grids), &
        opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
        opresx(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
        opresy(nx,ny,nz,n_grids), opresxy(nx,ny,nz,n_grids), &
        opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
    !
    epres(nx,ny,nz,n_grids)
    !
    box=n_grids
    !
    i=1
    do k=1,nz
        do j=1,ny
            bx(i,j,k,box)=bxp(j,k)-bx0(i,j,k,box)
            by(i,j,k,box)=byp(j,k)-by0(i,j,k,box)
            bz(i,j,k,box)=bzp(j,k)-bz0(i,j,k,box)
            !
            qrho(i,j,k,box)=rhop(j,k)
            qpx(i,j,k,box)=qrho(i,j,k,box)*svxp(j,k)
            qpy(i,j,k,box)=qrho(i,j,k,box)*svyp(j,k)
            qpz(i,j,k,box)=qrho(i,j,k,box)*svzp(j,k)
            qpresx(i,j,k,box)=0.5*spress
            qpresy(i,j,k,box)=qpresx(i,j,k,box)
            qpresz(i,j,k,box)=qpresx(i,j,k,box)
            qpresxy(i,j,k,box)=0.
            qpresxz(i,j,k,box)=0.
            qpresyz(i,j,k,box)=0.
            !
            hrho(i,j,k,box)=hfrac*frac_h*rhop(j,k)
            hpx(i,j,k,box)=hrho(i,j,k,box)*svxp(j,k)
            hpy(i,j,k,box)=hrho(i,j,k,box)*svyp(j,k)
            hpz(i,j,k,box)=hrho(i,j,k,box)*svzp(j,k)
            hpresx(i,j,k,box)=0.5*spress*hfrac
            hpresy(i,j,k,box)=hpresx(i,j,k,box)
            hpresz(i,j,k,box)=hpresx(i,j,k,box)
            hpresxy(i,j,k,box)=0.
            hpresxz(i,j,k,box)=0.
            hpresyz(i,j,k,box)=0.
            !
            orho(i,j,k,box)=ofrac*frac_o*rhop(j,k)
            opx(i,j,k,box)=orho(i,j,k,box)*svxp(j,k)
            opy(i,j,k,box)=orho(i,j,k,box)*svyp(j,k)
            opz(i,j,k,box)=orho(i,j,k,box)*svzp(j,k)
            opresx(i,j,k,box)=0.5*spress*ofrac
            opresy(i,j,k,box)=opresx(i,j,k,box)
            opresz(i,j,k,box)=opresx(i,j,k,box)
            opresxy(i,j,k,box)=0.
            opresxz(i,j,k,box)=0.
            opresyz(i,j,k,box)=0.
            !
            epres(i,j,k,box)=qpresx(i,j,k,box)/ti_te
            !
            avz=avz+svzp(j,k)
        enddo
    enddo
    !
    avz=avz/float(ny*nz)
    !      write(6,*)'set_imf avz',avz,ny,nz
    return
end
