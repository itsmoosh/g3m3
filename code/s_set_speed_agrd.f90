!
!	Checks for minimum rho and negative pressure
!	and resets value if necessary
!
subroutine set_speed_agrd( &
	qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
	hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
	orho,opresx,opresy,opresz,opx,opy,opz, &
	epres,qpresxy,qpresxz,qpresyz, &
	hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
	bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
	rho,presx,presy,presz,presxy,presxz,presyz,px,py,pz, &
	rmassq,rmassh,rmasso,nx,ny,nz,n_grids,box, &
	pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
	vlim,alf_lim,o_conc,fastest,isotropic)

	integer box
	dimension qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
	qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
	qpresxy(nx,ny,nz,n_grids), &
	qpresxz(nx,ny,nz,n_grids),qpresyz(nx,ny,nz,n_grids), &
	qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
	
	hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
	hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
	hpresxy(nx,ny,nz,n_grids), &
	hpresxz(nx,ny,nz,n_grids),hpresyz(nx,ny,nz,n_grids), &
	hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
	
	orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids), &
	opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids), &
	opresxy(nx,ny,nz,n_grids), &
	opresxz(nx,ny,nz,n_grids),opresyz(nx,ny,nz,n_grids), &
	opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
	epres(nx,ny,nz,n_grids)
	
	dimension rho(nx,ny,nz),presx(nx,ny,nz),presy(nx,ny,nz), &
	presz(nx,ny,nz),px(nx,ny,nz),py(nx,ny,nz),pz(nx,ny,nz), &
	presxy(nx,ny,nz),presxz(nx,ny,nz),presyz(nx,ny,nz)
	dimension bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
	bx0(nx,ny,nz,n_grids),by0(nx,ny,nz,n_grids),bz0(nx,ny,nz,n_grids), &
	bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz)
	
	real, intent(inout), dimension(nx,ny,nz,n_grids) :: &
		qrho, qpx,qpy,qpz, qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz, &
		hrho, hpx,hpy,hpz, hpresx,hpresy,hpresz,hpresxy,hpresxz,hpresyz, &
		orho, opx,opy,opz, opresx,opresy,opresz,opresxy,opresxz,opresyz, &
		bx,by,bz, epres, bx0,by0,bz0

	!	0 or 1 to efficiently reset pressure tensor cross terms when isotropic
	real :: iso_flag = 1.
	
	!	**********
	!	Executions
	!	**********	

	!	Multiply off-diagonal pressure tensor components by 0 if isotropic
	if(isotropic) iso_flag = 0.

	!	Determine speeds on the code and changes accordingly	

	fastest=0.
	
	pxmax=0.
	pymax=0.
	pzmax=0.
	pmax=0.
	csmax=0.
	alfmax=0.
	
	vqlim=vlim
	vhlim=1.4*vlim
	volim=2.*vlim
	
	cqlim=1.0*vqlim
	chlim=1.0*vhlim
	colim=4.0*volim
	cslim=vlim
	
	ani=0.9995
	
	!	Find bulk velocities: species 1
	
	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				rho(i,j,k)=qrho(i,j,k,box)
				px(i,j,k)=qpx(i,j,k,box)
				py(i,j,k)=qpy(i,j,k,box)
				pz(i,j,k)=qpz(i,j,k,box)
				
				!	Limit anisotropy
				
				apres=(qpresx(i,j,k,box) + qpresy(i,j,k,box) &
					+ qpresz(i,j,k,box))/3.0
				presx(i,j,k)=ani*qpresx(i,j,k,box)+(1.-ani)*apres
				presy(i,j,k)=ani*qpresy(i,j,k,box)+(1.-ani)*apres
				presz(i,j,k)=ani*qpresz(i,j,k,box)+(1.-ani)*apres
				presxy(i,j,k)=ani*qpresxy(i,j,k,box)
				presxz(i,j,k)=ani*qpresxz(i,j,k,box)
				presyz(i,j,k)=ani*qpresyz(i,j,k,box)
			enddo
		enddo
	enddo

	if(isotropic) then
		presx(:,:,:) = qpresx(:,:,:,box)
		presy(:,:,:) = qpresx(:,:,:,box)
		presz(:,:,:) = qpresx(:,:,:,box)
		presxy(:,:,:) = 0.
		presxz(:,:,:) = 0.
		presyz(:,:,:) = 0.
	endif
	
	d_min=0.0333
	
	do k=2,nz-1
		kp=k+1
		km=k-1

		do j=2,ny-1
			jp=j+1
			jm=j-1

			do i=2,nx-1
				ip=i+1
				im=i-1
				
				arho=qrho(i,j,k,box)+1.e-5
				aqvx=abs(qpx(i,j,k,box)/arho)
				aqvy=abs(qpy(i,j,k,box)/arho)
				aqvz=abs(qpz(i,j,k,box)/arho)
				vq=sqrt(aqvx**2+aqvy**2+aqvz**2)
				cq=sqrt(0.33*(qpresx(i,j,k,box)+qpresy(i,j,k,box) &
					+ qpresz(i,j,k,box))/arho)
				
				if ((vq.gt.vqlim).or.(cq.gt.cqlim))then
					qrho(i,j,k,box) = ( rho(ip,j,k)+rho(im,j,k) &
						+rho(i,jp,k)+rho(i,jm,k)+rho(i,j,kp) &
						+rho(i,j,km) )/6.
					qpx(i,j,k,box) = ( px(ip,j,k)+px(im,j,k) &
						+px(i,jp,k)+px(i,jm,k)+px(i,j,kp) &
						+px(i,j,km) )/6.
					qpy(i,j,k,box) = ( py(ip,j,k)+py(im,j,k) &
						+py(i,jp,k)+py(i,jm,k)+py(i,j,kp) &
						+py(i,j,km) )/6.
					qpz(i,j,k,box) = ( pz(ip,j,k)+pz(im,j,k) &
						+pz(i,jp,k)+pz(i,jm,k)+pz(i,j,kp) &
						+pz(i,j,km) )/6.
					qpresx(i,j,k,box) = ( presx(ip,j,k)+presx(im,j,k) &
						+presx(i,jp,k)+presx(i,jm,k)+presx(i,j,kp) &
						+presx(i,j,km) )/6.
					qpresy(i,j,k,box) = ( presy(ip,j,k)+presy(im,j,k) &
						+presy(i,jp,k)+presy(i,jm,k)+presy(i,j,kp) &
						+presy(i,j,km) )/6.
					qpresz(i,j,k,box) = ( presz(ip,j,k)+presz(im,j,k) &
						+presz(i,jp,k)+presz(i,jm,k)+presz(i,j,kp) &
						+presz(i,j,km) )/6.
					qpresxy(i,j,k,box) = iso_flag * ( presxy(ip,j,k) &
						+presxy(im,j,k)+presxy(i,jp,k)+presxy(i,jm,k) &
						+presxy(i,j,kp)+presxy(i,j,km) )/6.
					qpresxz(i,j,k,box) = iso_flag * ( presxz(ip,j,k) &
						+presxz(im,j,k)+presxz(i,jp,k)+presxz(i,jm,k) &
						+presxz(i,j,kp)+presxz(i,j,km) )/6.
					qpresyz(i,j,k,box) = iso_flag * ( presyz(ip,j,k) &
						+presyz(im,j,k)+presyz(i,jp,k)+presyz(i,jm,k) &
						+presyz(i,j,kp)+presyz(i,j,km) )/6.
				endif
				
				pxmax=amax1(pxmax,aqvx)
				pymax=amax1(pymax,aqvy)
				pzmax=amax1(pzmax,aqvz)
				csmax=amax1(csmax,cq)
			enddo
		enddo
	enddo
	
	!	Find bulk velocities: species 2
	
	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				rho(i,j,k)=hrho(i,j,k,box)
				px(i,j,k)=hpx(i,j,k,box)
				py(i,j,k)=hpy(i,j,k,box)
				pz(i,j,k)=hpz(i,j,k,box)
				
				!      limit anisotropy
				
				apres = ( hpresx(i,j,k,box)+hpresy(i,j,k,box) &
					+hpresz(i,j,k,box) )/3.0
				presx(i,j,k)=ani*hpresx(i,j,k,box)+(1.-ani)*apres
				presy(i,j,k)=ani*hpresy(i,j,k,box)+(1.-ani)*apres
				presz(i,j,k)=ani*hpresz(i,j,k,box)+(1.-ani)*apres
				presxy(i,j,k)=ani*hpresxy(i,j,k,box)
				presxz(i,j,k)=ani*hpresxz(i,j,k,box)
				presyz(i,j,k)=ani*hpresyz(i,j,k,box)
			enddo
		enddo
	enddo

	if(isotropic) then
		presx(:,:,:) = hpresx(:,:,:,box)
		presy(:,:,:) = hpresx(:,:,:,box)
		presz(:,:,:) = hpresx(:,:,:,box)
		presxy(:,:,:) = 0.
		presxz(:,:,:) = 0.
		presyz(:,:,:) = 0.
	endif
	
	d_min=0.001
	
	do k=2,nz-1
		kp=k+1
		km=k-1

		do j=2,ny-1
			jp=j+1
			jm=j-1

			do i=2,nx-1
				ip=i+1
				im=i-1
				arho=hrho(i,j,k,box)+1.e-5
				ahvx=abs(hpx(i,j,k,box)/arho)
				ahvy=abs(hpy(i,j,k,box)/arho)
				ahvz=abs(hpz(i,j,k,box)/arho)
				vh=sqrt(ahvx**2+ahvy**2+ahvz**2)
				ch=sqrt(0.33*(hpresx(i,j,k,box)+hpresy(i,j,k,box) &
					+hpresz(i,j,k,box))/arho)
				
				if ((vh.gt.vhlim).or.(ch.gt.chlim)) then
					hrho(i,j,k,box) = ( rho(ip,j,k)+rho(im,j,k) &
						+rho(i,jp,k)+rho(i,jm,k)+rho(i,j,kp) &
						+rho(i,j,km) )/6.
					hpx(i,j,k,box) = ( px(ip,j,k)+px(im,j,k) &
						+px(i,jp,k)+px(i,jm,k)+px(i,j,kp) &
						+px(i,j,km) )/6.
					hpy(i,j,k,box) = ( py(ip,j,k)+py(im,j,k) &
						+py(i,jp,k)+py(i,jm,k)+py(i,j,kp) &
						+py(i,j,km) )/6.
					hpz(i,j,k,box) = ( pz(ip,j,k)+pz(im,j,k) &
						+pz(i,jp,k)+pz(i,jm,k)+pz(i,j,kp) &
						+pz(i,j,km) )/6.
					hpresx(i,j,k,box) = ( presx(ip,j,k)+presx(im,j,k) &
						+presx(i,jp,k)+presx(i,jm,k)+presx(i,j,kp) &
						+presx(i,j,km) )/6.
					hpresy(i,j,k,box) = ( presy(ip,j,k)+presy(im,j,k) &
						+presy(i,jp,k)+presy(i,jm,k)+presy(i,j,kp) &
						+presy(i,j,km) )/6.
					hpresz(i,j,k,box) = ( presz(ip,j,k)+presz(im,j,k) &
						+presz(i,jp,k)+presz(i,jm,k)+presz(i,j,kp) &
						+presz(i,j,km) )/6.
					hpresxy(i,j,k,box) = iso_flag * ( presxy(ip,j,k) &
						+presxy(im,j,k)+presxy(i,jp,k)+presxy(i,jm,k) &
						+presxy(i,j,kp)+presxy(i,j,km) )/6.
					hpresxz(i,j,k,box) = iso_flag * ( presxz(ip,j,k) &
						+presxz(im,j,k)+presxz(i,jp,k)+presxz(i,jm,k) &
						+presxz(i,j,kp)+presxz(i,j,km) )/6.
					hpresyz(i,j,k,box) = iso_flag * ( presyz(ip,j,k) &
						+presyz(im,j,k)+presyz(i,jp,k)+presyz(i,jm,k) &
						+presyz(i,j,kp)+presyz(i,j,km) )/6.
				endif
				
				pxmax=amax1(pxmax,ahvx)
				pymax=amax1(pymax,ahvy)
				pzmax=amax1(pzmax,ahvz)
				csmax=amax1(csmax,ch)
			enddo
		enddo
	enddo
	
	!	Find bulk velocities: species 3
	
	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				rho(i,j,k)=orho(i,j,k,box)
				px(i,j,k)=opx(i,j,k,box)
				py(i,j,k)=opy(i,j,k,box)
				pz(i,j,k)=opz(i,j,k,box)
				
				!      limit anisotropy
				
				apres = ( opresx(i,j,k,box)+opresy(i,j,k,box) &
					+opresz(i,j,k,box) )/3.0
				presx(i,j,k)=ani*opresx(i,j,k,box)+(1.-ani)*apres
				presy(i,j,k)=ani*opresy(i,j,k,box)+(1.-ani)*apres
				presz(i,j,k)=ani*opresz(i,j,k,box)+(1.-ani)*apres
				presxy(i,j,k)=ani*opresxy(i,j,k,box)
				presxz(i,j,k)=ani*opresxz(i,j,k,box)
				presyz(i,j,k)=ani*opresyz(i,j,k,box)
			enddo
		enddo
	enddo
	
	if(isotropic) then
		presx(:,:,:) = opresx(:,:,:,box)
		presy(:,:,:) = opresx(:,:,:,box)
		presz(:,:,:) = opresx(:,:,:,box)
		presxy(:,:,:) = 0.
		presxz(:,:,:) = 0.
		presyz(:,:,:) = 0.
	endif

	d_min=0.001
	
	do k=2,nz-1
		kp=k+1
		km=k-1

		do j=2,ny-1
			jp=j+1
			jm=j-1

			do i=2,nx-1
				ip=i+1
				im=i-1
				arho=orho(i,j,k,box)+1.e-5
				aovx=abs(opx(i,j,k,box)/arho)
				aovy=abs(opy(i,j,k,box)/arho)
				aovz=abs(opz(i,j,k,box)/arho)
				vo=sqrt(aovx**2+aovy**2+aovz**2)
				co=sqrt(0.33*(opresx(i,j,k,box)+opresy(i,j,k,box) &
					+opresz(i,j,k,box))/arho)
				
				if ((vo.gt.volim).or.(co.gt.colim))then
					orho(i,j,k,box) = ( rho(ip,j,k)+rho(im,j,k) &
						+rho(i,jp,k)+rho(i,jm,k)+rho(i,j,kp) &
						+rho(i,j,km) )/6.
					opx(i,j,k,box) = ( px(ip,j,k)+px(im,j,k) &
						+px(i,jp,k)+px(i,jm,k)+px(i,j,kp) &
						+px(i,j,km) )/6.
					opy(i,j,k,box) = ( py(ip,j,k)+py(im,j,k) &
						+py(i,jp,k)+py(i,jm,k)+py(i,j,kp) &
						+py(i,j,km) )/6.
					opz(i,j,k,box) = ( pz(ip,j,k)+pz(im,j,k) &
						+pz(i,jp,k)+pz(i,jm,k)+pz(i,j,kp) &
						+pz(i,j,km) )/6.
					opresx(i,j,k,box) = ( presx(ip,j,k)+presx(im,j,k) &
						+presx(i,jp,k)+presx(i,jm,k)+presx(i,j,kp) &
						+presx(i,j,km) )/6.
					opresy(i,j,k,box) = ( presy(ip,j,k)+presy(im,j,k) &
						+presy(i,jp,k)+presy(i,jm,k)+presy(i,j,kp) &
						+presy(i,j,km) )/6.
					opresz(i,j,k,box) = ( presz(ip,j,k)+presz(im,j,k) &
						+presz(i,jp,k)+presz(i,jm,k)+presz(i,j,kp) &
						+presz(i,j,km) )/6.
					opresxy(i,j,k,box) = iso_flag * ( presxy(ip,j,k) &
						+presxy(im,j,k)+presxy(i,jp,k)+presxy(i,jm,k) &
						+presxy(i,j,kp)+presxy(i,j,km) )/6.
					opresxz(i,j,k,box) = iso_flag * ( presxz(ip,j,k) &
						+presxz(im,j,k)+presxz(i,jp,k)+presxz(i,jm,k) &
						+presxz(i,j,kp)+presxz(i,j,km) )/6.
					opresyz(i,j,k,box) = iso_flag * ( presyz(ip,j,k) &
						+presyz(im,j,k)+presyz(i,jp,k)+presyz(i,jm,k) &
						+presyz(i,j,kp)+presyz(i,j,km) )/6.
				endif
				
				pxmax=amax1(pxmax,aovx)
				pymax=amax1(pymax,aovy)
				pzmax=amax1(pzmax,aovz)
				csmax=amax1(csmax,co)
			enddo
		enddo
	enddo
	
	!	Electron pressure and alfven speed
	
	bsx(:,:,:) = bx0(:,:,:,box) + bx(:,:,:,box)
	bsy(:,:,:) = by0(:,:,:,box) + by(:,:,:,box)
	bsz(:,:,:) = bz0(:,:,:,box) + bz(:,:,:,box)

	!	Find magnitude of b
	
	call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
	
	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				px(i,j,k)=bx(i,j,k,box)
				py(i,j,k)=by(i,j,k,box)
				pz(i,j,k)=bz(i,j,k,box)
			enddo
		enddo
	enddo
	
	do k=2,nz-1
		kp=k+1
		km=k-1

		do j=2,ny-1
			jp=j+1
			jm=j-1

			do i=2,nx-1
				ip=i+1
				im=i-1
				
				!	Electron pressure
				
				arho=qrho(i,j,k,box)+hrho(i,j,k,box)+orho(i,j,k,box) &
					+ 1.e-5
				cs=sqrt(epres(i,j,k,box)/arho)
				if(cs.gt.cslim) then
					epres(i,j,k,box)=epres(i,j,k,box)*cslim/cs
				endif
				
				!	Find alfven speed
				
				abfld=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
					+bsz(i,j,k)**2)
				alf=abfld/sqrt(arho)
				alfmax=amax1(alfmax,alf)
			enddo
		enddo
	enddo
	
	pmax=sqrt(pxmax**2+pymax**2+pzmax**2)
	
	fastest=amax1(fastest,sqrt(pxmax**2+pymax**2+pzmax**2 &
		+csmax**2+alfmax**2))
	
	return
end subroutine set_speed_agrd
