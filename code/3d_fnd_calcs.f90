!
!	This file contains six subroutines:
!	fnd_fld
!	fnd_vel
!	fnd_vtot
!	fnd_evel
!	fnd_prss
!	fnd_pres
!
subroutine fnd_fld(bx,btx,nx,ny,nz,n_grids,box)
    !
    !
    !     calculates the total dynamic magnetic field only
    !        dipole field added separately by by setting add_dip=.true.
    !
    integer box
    dimension bx(nx,ny,nz,n_grids),btx(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                btx(i,j,k)=btx(i,j,k)+bx(i,j,k,box)
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
subroutine fnd_vel(px,py,pz,rho,vx,vy,vz,nx,ny,nz,n_grids,box)
    !
    !     converts momentum into velocity for graphics
    !
    integer box
    dimension px(nx,ny,nz,n_grids),py(nx,ny,nz,n_grids), &
        pz(nx,ny,nz,n_grids),rho(nx,ny,nz,n_grids), &
        vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                arho=amax1(rho(i,j,k,box),0.0001)
                vx(i,j,k)=px(i,j,k,box)/arho
                vy(i,j,k)=py(i,j,k,box)/arho
                vz(i,j,k)=pz(i,j,k,box)/arho
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
    opx,opy,opz,orho,vx,vy,vz,nx,ny,nz,n_grids,box, &
    rmassq,rmassh,rmasso)
    !
    !     converts momentum into velocity for graphics
    !
    integer box
	real tden
    dimension qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids), &
        qpz(nx,ny,nz,n_grids),qrho(nx,ny,nz,n_grids), &
        hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids), &
        hpz(nx,ny,nz,n_grids),hrho(nx,ny,nz,n_grids), &
        opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids), &
        opz(nx,ny,nz,n_grids),orho(nx,ny,nz,n_grids), &
        vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                tden = qrho(i,j,k,box) + hrho(i,j,k,box) + orho(i,j,k,box)
                !
                vx(i,j,k)=( qpx(i,j,k,box)+hpx(i,j,k,box) &
                    +opx(i,j,k,box) ) / tden
                vy(i,j,k)=( qpy(i,j,k,box)+hpy(i,j,k,box) &
                    +opy(i,j,k,box) ) / tden
                vz(i,j,k)=( qpz(i,j,k,box)+hpz(i,j,k,box) &
                    +opz(i,j,k,box) ) / tden
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
    nx,ny,nz,n_grids,box,rmassq,rmassh,rmasso,reynolds)
    !
    !     converts momentum into velocity for graphics
    !
    integer box
    dimension qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids), &
        qpz(nx,ny,nz,n_grids),qrho(nx,ny,nz,n_grids), &
        hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids), &
        hpz(nx,ny,nz,n_grids),hrho(nx,ny,nz,n_grids), &
        opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids), &
        opz(nx,ny,nz,n_grids),orho(nx,ny,nz,n_grids), &
        curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
        tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
        evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
    !
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                !      eden=amax1(qrho(i,j,k,box)+hrho(i,j,k,box)+orho(i,j,k,box),
                !    +                 0.0001)
                qden=(qrho(i,j,k,box)+0.000001)/rmassq
                hden=(hrho(i,j,k,box)+0.000001)/rmassh
                oden=(orho(i,j,k,box)+0.000001)/rmasso
                tden=qden+hden+oden
                !
                !      keep sepearate the ion and current components
                !
                tvx(i,j,k)=(qpx(i,j,k,box)/rmassq+hpx(i,j,k,box)/rmassh &
                    +opx(i,j,k,box)/rmasso)/tden
                evx(i,j,k)= tvx(i,j,k) - curx(i,j,k)/tden/reynolds
                !
                tvy(i,j,k)=(qpy(i,j,k,box)/rmassq+hpy(i,j,k,box)/rmassh &
                    +opy(i,j,k,box)/rmasso)/tden
                evy(i,j,k)= tvy(i,j,k) - cury(i,j,k)/tden/reynolds
                !
                tvz(i,j,k)=(qpz(i,j,k,box)/rmassq+hpz(i,j,k,box)/rmassh &
                    +opz(i,j,k,box)/rmasso)/tden
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
subroutine fnd_prss(px,py,pz,rho,erg,afld,nx,ny,nz,n_grids, &
    box,gamma1,igo)
    !
    !     now using pressure equation rather than an erg equation
    !
    integer box
    dimension px(nx,ny,nz,n_grids),py(nx,ny,nz,n_grids), &
        pz(nx,ny,nz,n_grids),rho(nx,ny,nz,n_grids), &
        erg(nx,ny,nz,n_grids),afld(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                afld(i,j,k)=erg(i,j,k,box)
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
!     *************************************************
!
!
subroutine fnd_pres(px,py,pz,pxy,pxz,pyz,p_para,p_perp,p_cross, &
	vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,n_grids,box)

	! This subroutine is only used in write_graphing_data.

	dimension px(nx,ny,nz,n_grids),py(nx,ny,nz,n_grids),pz(nx,ny,nz,n_grids), &
		pxy(nx,ny,nz,n_grids),pxz(nx,ny,nz,n_grids),pyz(nx,ny,nz,n_grids), &
		p_para(nx,ny,nz,n_grids),p_perp(nx,ny,nz,n_grids),p_cross(nx,ny,nz,n_grids), &
		vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz),bxt(nx,ny,nz), &
		byt(nx,ny,nz),bzt(nx,ny,nz)

	integer, intent(in) :: box

	real abx,aby,abz,bmag,avx,avy,avz,vmag
	real vcbx,vcby,vcbz,vcmag
	real p_para,p_perp,p_cross
	real presx,presy,presz,presmag
	real presxy,presxz,presyz

	!     Find anisotropy values for species 1
	!
	do k=1,nz
		do j=1,ny
			do i=1,nx
			!
			presx=px(i,j,k,box)
			presy=py(i,j,k,box)
			presz=pz(i,j,k,box)
			presxy=pxy(i,j,k,box)
			presxz=pxz(i,j,k,box)
			presyz=pyz(i,j,k,box)

			presmag=sqrt(presx**2+presy**2+presz**2)+1.e-11
			!
			abx=bxt(i,j,k)
			aby=byt(i,j,k)
			abz=bzt(i,j,k)
			bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
			!
			avx=vvx(i,j,k)
			avy=vvy(i,j,k)
			avz=vvz(i,j,k)
			vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
			!
			vcbx=(avy*abz-avz*aby)/bmag
			vcby=-(avx*abz-avz*abx)/bmag
			vcbz=(avx*aby-avy*abx)/bmag

			vcmag=sqrt(vcbx**2+vcby**2+vcbz**2) &
				+1.e-12

			p_para(i,j,k,box)=sqrt((abx*presx)**2 + (aby*presy)**2 + (abz*presz)**2) / bmag
			p_cross(i,j,k,box)=sqrt((vcbx*presx)**2 + (vcby*presy)**2 + (vcbz*presz)**2) / vcmag
			p_perp(i,j,k,box)=sqrt( abs( presmag**2-( p_para(i,j,k,box)**2+p_cross(i,j,k,box)**2 ) ) )

			enddo
		enddo
	enddo

end subroutine fnd_pres
