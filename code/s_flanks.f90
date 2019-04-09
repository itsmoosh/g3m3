!
!	Sets the outer boundaries of fine grid at position ns
!	using data from coarse grid at position nb
!
subroutine flanks(wrkrho,rho,oldrho,nx,ny,nz,n_grids,box, &
	t_new,t_old,t,work, &
	grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
	
	integer, parameter :: dp = kind(1.d0)
	integer box
	real(dp) t_old(n_grids), t_new(n_grids), t
	dimension rho(nx,ny,nz,n_grids),oldrho(nx,ny,nz,n_grids), &
	wrkrho(nx,ny,nz,n_grids),work(nx,ny,nz)
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)

	mb=box+1
	ms=box
	
	t1=(t_new(mb)-t)/(t_new(mb)-t_old(mb))
	t2=1.-t1
	
	delx=(grd_xmax(mb)-grd_xmin(mb))/(nx-1.)
	dely=(grd_ymax(mb)-grd_ymin(mb))/(ny-1.)
	delz=(grd_zmax(mb)-grd_zmin(mb))/(nz-1.)
	
	sx=(grd_xmax(ms)-grd_xmin(ms))/(nx-1.)
	sy=(grd_ymax(ms)-grd_ymin(ms))/(ny-1.)
	sz=(grd_zmax(ms)-grd_zmin(ms))/(nz-1.)
	
	!	Time interpolation
	
	!$omp  parallel do
	do k=1,nz
		do j=1,ny
			do i=1,nx
				work(i,j,k)=oldrho(i,j,k,mb)*t1+rho(i,j,k,mb)*t2
			enddo
		enddo
	enddo
	
	!	Interpolate onto grid
	
	!	Bottom panels
	
	!$omp  parallel do
	do kk=2,1,-1
		ks=kk
		az=grd_zmin(ms)+sz*(ks-1)
		ak=1.+(az-grd_zmin(mb))/delz
		kb=ak
		kb1=kb+1
		dz=ak-kb
		ddz=1.-dz
		frac=(kk-1)/2.
		frac1=1.-frac
		
		do js=1,ny
			ay=grd_ymin(ms)+sy*(js-1)
			aj=1.+(ay-grd_ymin(mb))/dely
			jb=aj
			jb1=jb+1
			dy=aj-jb
			ddy=1.-dy

			do is=1,nx
				ax=grd_xmin(ms)+sx*(is-1)
				ai=1.+(ax-grd_xmin(mb))/delx
				ib=ai
				ib1=ib+1
				dx=ai-ib
				ddx=1.-dx
				
				wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
					+frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
					+work(ib,jb,kb1)*ddx*ddy*dz &
					+work(ib,jb1,kb)*ddx*dy*ddz &
					+work(ib,jb1,kb1)*ddx*dy*dz &
					+work(ib1,jb,kb)*dx*ddy*ddz &
					+work(ib1,jb,kb1)*dx*ddy*dz &
					+work(ib1,jb1,kb)*dx*dy*ddz &
					+work(ib1,jb1,kb1)*dx*dy*dz)
			enddo
		enddo
	enddo
	
	!	Top panels
	
	!$omp  parallel do
	do kk=2,1,-1
		ks=nz-kk+1
		az=grd_zmin(ms)+sz*(ks-1)
		ak=1.+(az-grd_zmin(mb))/delz
		kb=ak
		kb1=kb+1
		dz=ak-kb
		ddz=1.-dz
		frac=(kk-1)/2.
		frac1=1.-frac
		
		do js=1,ny
			ay=grd_ymin(ms)+sy*(js-1)
			aj=1.+(ay-grd_ymin(mb))/dely
			jb=aj
			jb1=jb+1
			dy=aj-jb
			ddy=1.-dy

			do is=1,nx
				ax=grd_xmin(ms)+sx*(is-1)
				ai=1.+(ax-grd_xmin(mb))/delx
				ib=ai
				ib1=ib+1
				dx=ai-ib
				ddx=1.-dx
				
				wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
					+frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
					+work(ib,jb,kb1)*ddx*ddy*dz &
					+work(ib,jb1,kb)*ddx*dy*ddz &
					+work(ib,jb1,kb1)*ddx*dy*dz &
					+work(ib1,jb,kb)*dx*ddy*ddz &
					+work(ib1,jb,kb1)*dx*ddy*dz &
					+work(ib1,jb1,kb)*dx*dy*ddz &
					+work(ib1,jb1,kb1)*dx*dy*dz)
			enddo
		enddo
	enddo
	
	!	Front panels
	
	!$omp  parallel do
	do ks=1,nz
		az=grd_zmin(ms)+sz*(ks-1)
		ak=1.+(az-grd_zmin(mb))/delz
		kb=ak
		kb1=kb+1
		dz=ak-kb
		ddz=1.-dz

		do jj=2,1,-1
			js=jj
			ay=grd_ymin(ms)+sy*(js-1)
			aj=1.+(ay-grd_ymin(mb))/dely
			jb=aj
			jb1=jb+1
			dy=aj-jb
			ddy=1.-dy
			frac=(jj-1)/2.
			frac1=1.-frac

			do is=1,nx
				ax=grd_xmin(ms)+sx*(is-1)
				ai=1.+(ax-grd_xmin(mb))/delx
				ib=ai
				ib1=ib+1
				dx=ai-ib
				ddx=1.-dx
				
				wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
					+frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
					+work(ib,jb,kb1)*ddx*ddy*dz &
					+work(ib,jb1,kb)*ddx*dy*ddz &
					+work(ib,jb1,kb1)*ddx*dy*dz &
					+work(ib1,jb,kb)*dx*ddy*ddz &
					+work(ib1,jb,kb1)*dx*ddy*dz &
					+work(ib1,jb1,kb)*dx*dy*ddz &
					+work(ib1,jb1,kb1)*dx*dy*dz)
			enddo
		enddo
	enddo
	
	!	Back panels
	
	!$omp  parallel do
	do ks=1,nz
		az=grd_zmin(ms)+sz*(ks-1)
		ak=1.+(az-grd_zmin(mb))/delz
		kb=ak
		kb1=kb+1
		dz=ak-kb
		ddz=1.-dz

		do jj=2,1,-1
			js=ny-jj+1
			ay=grd_ymin(ms)+sy*(js-1)
			aj=1.+(ay-grd_ymin(mb))/dely
			jb=aj
			jb1=jb+1
			dy=aj-jb
			ddy=1.-dy
			frac=(jj-1)/2.
			frac1=1.-frac

			do is=1,nx
				ax=grd_xmin(ms)+sx*(is-1)
				ai=1.+(ax-grd_xmin(mb))/delx
				ib=ai
				ib1=ib+1
				dx=ai-ib
				ddx=1.-dx
				
				wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
					+frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
					+work(ib,jb,kb1)*ddx*ddy*dz &
					+work(ib,jb1,kb)*ddx*dy*ddz &
					+work(ib,jb1,kb1)*ddx*dy*dz &
					+work(ib1,jb,kb)*dx*ddy*ddz &
					+work(ib1,jb,kb1)*dx*ddy*dz &
					+work(ib1,jb1,kb)*dx*dy*ddz &
					+work(ib1,jb1,kb1)*dx*dy*dz)
			
			enddo
		enddo
	enddo
	
	!	Front-left panels
	
	!$omp  parallel do
	do ks=1,nz
		az=grd_zmin(ms)+sz*(ks-1)
		ak=1.+(az-grd_zmin(mb))/delz
		kb=ak
		kb1=kb+1
		dz=ak-kb
		ddz=1.-dz

		do js=1,ny
			ay=grd_ymin(ms)+sy*(js-1)
			aj=1.+(ay-grd_ymin(mb))/dely
			jb=aj
			jb1=jb+1
			dy=aj-jb
			ddy=1.-dy

			do ii=2,1,-1
				is=ii
				ax=grd_xmin(ms)+sx*(is-1)
				ai=1.+(ax-grd_xmin(mb))/delx
				ib=ai
				ib1=ib+1
				dx=ai-ib
				ddx=1.-dx
				frac=(ii-1)/2.
				frac1=1.-frac
				
				wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
					+frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
					+work(ib,jb,kb1)*ddx*ddy*dz &
					+work(ib,jb1,kb)*ddx*dy*ddz &
					+work(ib,jb1,kb1)*ddx*dy*dz &
					+work(ib1,jb,kb)*dx*ddy*ddz &
					+work(ib1,jb,kb1)*dx*ddy*dz &
					+work(ib1,jb1,kb)*dx*dy*ddz &
					+work(ib1,jb1,kb1)*dx*dy*dz)
			enddo
		enddo
	enddo
	
	!$omp  parallel do
	do ks=1,nz
		az=grd_zmin(ms)+sz*(ks-1)
		ak=1.+(az-grd_zmin(mb))/delz
		kb=ak
		kb1=kb+1
		dz=ak-kb
		ddz=1.-dz

		do js=1,ny
			ay=grd_ymin(ms)+sy*(js-1)
			aj=1.+(ay-grd_ymin(mb))/dely
			jb=aj
			jb1=jb+1
			dy=aj-jb
			ddy=1.-dy

			do ii=2,1,-1
				is=nx-ii+1
				ax=grd_xmin(ms)+sx*(is-1)
				ai=1.+(ax-grd_xmin(mb))/delx
				ib=ai
				ib1=ib+1
				dx=ai-ib
				ddx=1.-dx
				frac=(ii-1)/2.
				frac1=1.-frac
				
				wrkrho(is,js,ks,ms)=frac*wrkrho(is,js,ks,ms) &
					+frac1*(work(ib,jb,kb)*ddx*ddy*ddz &
					+work(ib,jb,kb1)*ddx*ddy*dz &
					+work(ib,jb1,kb)*ddx*dy*ddz &
					+work(ib,jb1,kb1)*ddx*dy*dz &
					+work(ib1,jb,kb)*dx*ddy*ddz &
					+work(ib1,jb,kb1)*dx*ddy*dz &
					+work(ib1,jb1,kb)*dx*dy*ddz &
					+work(ib1,jb1,kb1)*dx*dy*dz)
			enddo
		enddo
	enddo
	
	return
end subroutine flanks
