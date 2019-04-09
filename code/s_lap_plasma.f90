!
!	Apply the lapidus smoothing technique for particular ion component
!
subroutine lap_plasma(rho,px,py,pz, &
	presx,presy,presz,presxy,presxz,presyz, &
	wrkrho,wrkpx,wrkpy,wrkpz, &
	wrkpresx,wrkpresy,wrkpresz, &
	wrkpresxy,wrkpresxz,wrkpresyz, &
	vx,vy,vz, &
	nx,ny,nz,n_grids,box,chirho,chipxyz,chierg,delt,isotropic)

	integer box
	dimension rho(nx,ny,nz,n_grids),presx(nx,ny,nz,n_grids), &
	presy(nx,ny,nz,n_grids),presz(nx,ny,nz,n_grids), &
	px(nx,ny,nz,n_grids),py(nx,ny,nz,n_grids),pz(nx,ny,nz,n_grids), &
	wrkrho(nx,ny,nz,n_grids),wrkpresx(nx,ny,nz,n_grids), &
	wrkpresy(nx,ny,nz,n_grids),wrkpresz(nx,ny,nz,n_grids), &
	wrkpx(nx,ny,nz,n_grids),wrkpy(nx,ny,nz,n_grids), &
	wrkpz(nx,ny,nz,n_grids), &
	vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
	
	dimension presxy(nx,ny,nz,n_grids),presxz(nx,ny,nz,n_grids), &
	presyz(nx,ny,nz,n_grids),wrkpresxy(nx,ny,nz,n_grids), &
	wrkpresxz(nx,ny,nz,n_grids),wrkpresyz(nx,ny,nz,n_grids)
	
	common /gridding/grd_xmin(9),grd_xmax(9),grd_ymin(9),grd_ymax(9), &
	grd_zmin(9),grd_zmax(9),xspac(9),ncore(9),nbndry(9), &
	rx,ry,rz,xdip,ydip,zdip,r_inner,b0, &
	sin_tilt,cos_tilt
	
	!	Equal timing irrespective of grid size
	
	deltz=delt
	delty=delt
	deltx=delt
	
	!$omp  parallel do
	do  k=2,nz-1
		km=k-1
		kp=k+1

		do  j=2,ny-1
			jm=j-1
			jp=j+1

			do i=2,nx-1
				im=i-1
				ip=i+1
				
				uxp1=abs(vx(ip,j,k)-vx(i,j,k))
				uxm1=abs(vx(i,j,k)-vx(im,j,k))
				
				uyp1=abs(vy(i,jp,k)-vy(i,j,k))
				uym1=abs(vy(i,j,k)-vy(i,jm,k))
				
				uzp1=abs(vz(i,j,kp)-vz(i,j,k))
				uzm1=abs(vz(i,j,k)-vz(i,j,km))
				
				wrkrho(i,j,k,box)=rho(i,j,k,box)+chirho*( &
					deltx*(uxp1*(rho(ip,j,k,box)-rho(i,j,k,box)) &
					-uxm1*(rho(i,j,k,box)-rho(im,j,k,box))) &
					+delty*(uyp1*(rho(i,jp,k,box)-rho(i,j,k,box)) &
					-uym1*(rho(i,j,k,box)-rho(i,jm,k,box))) &
					+deltz*(uzp1*(rho(i,j,kp,box)-rho(i,j,k,box)) &
					-uzm1*(rho(i,j,k,box)-rho(i,j,km,box))) )
				
				wrkpx(i,j,k,box)=px(i,j,k,box)+chipxyz*( &
					+deltx*(uxp1*(px(ip,j,k,box)-px(i,j,k,box)) &
					-uxm1*(px(i,j,k,box)-px(im,j,k,box))) &
					+delty*(uyp1*(px(i,jp,k,box)-px(i,j,k,box)) &
					-uym1*(px(i,j,k,box)-px(i,jm,k,box))) &
					+deltz*(uzp1*(px(i,j,kp,box)-px(i,j,k,box)) &
					-uzm1*(px(i,j,k,box)-px(i,j,km,box))) )
				
				wrkpy(i,j,k,box)=py(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(py(ip,j,k,box)-py(i,j,k,box)) &
					-uxm1*(py(i,j,k,box)-py(im,j,k,box))) &
					+delty*(uyp1*(py(i,jp,k,box)-py(i,j,k,box)) &
					-uym1*(py(i,j,k,box)-py(i,jm,k,box))) &
					+deltz*(uzp1*(py(i,j,kp,box)-py(i,j,k,box)) &
					-uzm1*(py(i,j,k,box)-py(i,j,km,box))) )
				
				wrkpz(i,j,k,box)=pz(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(pz(ip,j,k,box)-pz(i,j,k,box)) &
					-uxm1*(pz(i,j,k,box)-pz(im,j,k,box))) &
					+delty*(uyp1*(pz(i,jp,k,box)-pz(i,j,k,box)) &
					-uym1*(pz(i,j,k,box)-pz(i,jm,k,box))) &
					+deltz*(uzp1*(pz(i,j,kp,box)-pz(i,j,k,box)) &
					-uzm1*(pz(i,j,k,box)-pz(i,j,km,box))) )
				
				!	Diagonal elements
				
				wrkpresx(i,j,k,box)=presx(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(presx(ip,j,k,box)-presx(i,j,k,box)) &
					-uxm1*(presx(i,j,k,box)- presx(im,j,k,box))) &
					+delty*(uyp1*(presx(i,jp,k,box)-presx(i,j,k,box)) &
					-uym1*(presx(i,j,k,box)- presx(i,jm,k,box))) &
					+deltz*(uzp1*(presx(i,j,kp,box)-presx(i,j,k,box)) &
					-uzm1*(presx(i,j,k,box)- presx(i,j,km,box))) )
				
				wrkpresy(i,j,k,box)=presy(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(presy(ip,j,k,box)-presy(i,j,k,box)) &
					-uxm1*(presy(i,j,k,box)- presy(im,j,k,box))) &
					+delty*(uyp1*(presy(i,jp,k,box)-presy(i,j,k,box)) &
					-uym1*(presy(i,j,k,box)- presy(i,jm,k,box))) &
					+deltz*(uzp1*(presy(i,j,kp,box)-presy(i,j,k,box)) &
					-uzm1*(presy(i,j,k,box)- presy(i,j,km,box))) )
				
				wrkpresz(i,j,k,box)=presz(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(presz(ip,j,k,box)-presz(i,j,k,box)) &
					-uxm1*(presz(i,j,k,box)- presz(im,j,k,box))) &
					+delty*(uyp1*(presz(i,jp,k,box)-presz(i,j,k,box)) &
					-uym1*(presz(i,j,k,box)- presz(i,jm,k,box))) &
					+deltz*(uzp1*(presz(i,j,kp,box)-presz(i,j,k,box)) &
					-uzm1*(presz(i,j,k,box)- presz(i,j,km,box))) )
				
				!	Off-diagonal elements
				
				wrkpresxy(i,j,k,box)=presxy(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(presxy(ip,j,k,box)-presxy(i,j,k,box)) &
					-uxm1*(presxy(i,j,k,box)- presxy(im,j,k,box))) &
					+delty*(uyp1*(presxy(i,jp,k,box)-presxy(i,j,k,box) &
					) -uym1*(presxy(i,j,k,box)- presxy(i,jm,k,box))) &
					+deltz*(uzp1*(presxy(i,j,kp,box)-presxy(i,j,k,box) &
					) -uzm1*(presxy(i,j,k,box)- presxy(i,j,km,box))) )
				
				wrkpresxz(i,j,k,box)=presxz(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(presxz(ip,j,k,box)-presxz(i,j,k,box)) &
					-uxm1*(presxz(i,j,k,box)- presxz(im,j,k,box))) &
					+delty*(uyp1*(presxz(i,jp,k,box)-presxz(i,j,k,box) &
					) -uym1*(presxz(i,j,k,box)- presxz(i,jm,k,box))) &
					+deltz*(uzp1*(presxz(i,j,kp,box)-presxz(i,j,k,box) &
					) -uzm1*(presxz(i,j,k,box)- presxz(i,j,km,box))) )
				
				wrkpresyz(i,j,k,box)=presyz(i,j,k,box)+chipxyz*( &
					deltx*(uxp1*(presyz(ip,j,k,box)-presyz(i,j,k,box)) &
					-uxm1*(presyz(i,j,k,box)- presyz(im,j,k,box))) &
					+delty*(uyp1*(presyz(i,jp,k,box)-presyz(i,j,k,box) &
					) -uym1*(presyz(i,j,k,box)- presyz(i,jm,k,box))) &
					+deltz*(uzp1*(presyz(i,j,kp,box)-presyz(i,j,k,box) &
					) -uzm1*(presyz(i,j,k,box)- presyz(i,j,km,box))) )
			enddo
		enddo
	enddo

	return
end subroutine lap_plasma
