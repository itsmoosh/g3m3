!
!	Calculates para/perp/cross terms of the pressure tensor, for
!	plotting purposes.
!
subroutine fnd_pres(px,py,pz,pxy,pxz,pyz,p_para,p_perp,p_cross, &
	vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,n_grids,box)

	implicit none

	integer, intent(in) :: nx, ny, nz, n_grids, box
	real, intent(in) :: px(nx,ny,nz,n_grids), py(nx,ny,nz,n_grids), &
		pz(nx,ny,nz,n_grids), pxy(nx,ny,nz,n_grids), &
		pxz(nx,ny,nz,n_grids), pyz(nx,ny,nz,n_grids), &
		vvx(nx,ny,nz), vvy(nx,ny,nz), vvz(nx,ny,nz), &
		bxt(nx,ny,nz), byt(nx,ny,nz), bzt(nx,ny,nz)
	real, intent(out) :: p_para(nx,ny,nz,n_grids), &
		p_perp(nx,ny,nz,n_grids), p_cross(nx,ny,nz,n_grids)

	integer i,j,k
	real abx,aby,abz,bmag,avx,avy,avz,vmag
	real vcbx,vcby,vcbz,vcmag
	real presx,presy,presz,presmag
	real presxy,presxz,presyz

	!     Find anisotropy values for species 1
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

			presmag=sqrt(presx**2+presy**2+presz**2) + 1.e-11
			
			abx=bxt(i,j,k)
			aby=byt(i,j,k)
			abz=bzt(i,j,k)
			bmag=sqrt(abx**2+aby**2+abz**2) + 1.e-12
			
			avx=vvx(i,j,k)
			avy=vvy(i,j,k)
			avz=vvz(i,j,k)
			vmag=sqrt(avx**2+avy**2+avz**2) + 1.e-8
			
			vcbx=(avy*abz-avz*aby)/bmag
			vcby=-(avx*abz-avz*abx)/bmag
			vcbz=(avx*aby-avy*abx)/bmag

			vcmag=sqrt(vcbx**2+vcby**2+vcbz**2) + 1.e-12

			p_para(i,j,k,box)=sqrt((abx*presx)**2 + (aby*presy)**2 &
				+ (abz*presz)**2) / bmag
			p_cross(i,j,k,box)=sqrt((vcbx*presx)**2 + (vcby*presy)**2 &
				+ (vcbz*presz)**2) / vcmag
			p_perp(i,j,k,box)=sqrt( abs( presmag**2 &
				- ( p_para(i,j,k,box)**2+p_cross(i,j,k,box)**2 ) ) )

			enddo
		enddo
	enddo

	return
end subroutine fnd_pres
