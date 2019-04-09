!
!	Apply straight artifical diffusion.
!
subroutine diffuse(qrho,qpres,qpx,qpy,qpz, &
	hrho,hpres,hpx,hpy,hpz, &
	orho,opres,opx,opy,opz,epres,bx,by,bz, &
	wrkqrho,wrkqpres,wrkqpx,wrkqpy,wrkqpz, &
	wrkhrho,wrkhpres,wrkhpx,wrkhpy,wrkhpz,wrkorho, &
	wrkopres,wrkopx,wrkopy,wrkopz, &
	wrkepres,wrkbx,wrkby,wrkbz,nx,ny,nz,n_grids,box, &
	difrho,difpxyz,diferg,delt,xspac,yspac,zspac)

	integer box
	dimension qrho(nx,ny,nz,n_grids),qpres(nx,ny,nz,n_grids), &
	qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
	hrho(nx,ny,nz,n_grids),hpres(nx,ny,nz,n_grids), &
	hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
	orho(nx,ny,nz,n_grids),opres(nx,ny,nz,n_grids), &
	opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
	bx(nx,ny,nz,n_grids),by(nx,ny,nz,n_grids),bz(nx,ny,nz,n_grids), &
	epres(nx,ny,nz,n_grids)
	
	dimension wrkqrho(nx,ny,nz,n_grids),wrkqpres(nx,ny,nz,n_grids), &
	wrkqpx(nx,ny,nz,n_grids),wrkqpy(nx,ny,nz,n_grids), &
	wrkqpz(nx,ny,nz,n_grids),wrkhrho(nx,ny,nz,n_grids), &
	wrkhpres(nx,ny,nz,n_grids),wrkhpx(nx,ny,nz,n_grids), &
	wrkhpy(nx,ny,nz,n_grids),wrkhpz(nx,ny,nz,n_grids), &
	wrkorho(nx,ny,nz,n_grids),wrkopres(nx,ny,nz,n_grids), &
	wrkopx(nx,ny,nz,n_grids),wrkopy(nx,ny,nz,n_grids), &
	wrkopz(nx,ny,nz,n_grids),wrkbx(nx,ny,nz,n_grids), &
	wrkby(nx,ny,nz,n_grids),wrkbz(nx,ny,nz,n_grids), &
	wrkepres(nx,ny,nz,n_grids)
	
	difb=diferg*delt/xspac
	dife=difrho*delt/xspac
	difp=difpxyz*delt/xspac
	difr=difrho*delt/xspac
	difo=0.25*difr
	
	!$omp  parallel do
	do k=2,nz-1
		km=k-1
		kp=k+1

		do j=2,ny-1
			jm=j-1
			jp=j+1

			do i=2,nx-1
				im=i-1
				ip=i+1
				
				!	Species 1
				
				qrho(i,j,k,box)=wrkqrho(i,j,k,box) +difr*( &
					+wrkqrho(ip,j,k,box)+wrkqrho(im,j,k,box) &
					+wrkqrho(i,jp,k,box)+wrkqrho(i,jm,k,box) &
					+wrkqrho(i,j,kp,box)+wrkqrho(i,j,km,box) &
					-6.*wrkqrho(i,j,k,box) )
				
				qpx(i,j,k,box)=wrkqpx(i,j,k,box) +difp*( &
					+wrkqpx(ip,j,k,box)+wrkqpx(im,j,k,box) &
					+wrkqpx(i,jp,k,box)+wrkqpx(i,jm,k,box) &
					+wrkqpx(i,j,kp,box)+wrkqpx(i,j,km,box) &
					-6.*wrkqpx(i,j,k,box) )
				
				qpy(i,j,k,box)=wrkqpy(i,j,k,box) +difp*( &
					+wrkqpy(ip,j,k,box)+wrkqpy(im,j,k,box) &
					+wrkqpy(i,jp,k,box)+wrkqpy(i,jm,k,box) &
					+wrkqpy(i,j,kp,box)+wrkqpy(i,j,km,box) &
					-6.*wrkqpy(i,j,k,box) )
				
				qpz(i,j,k,box)=wrkqpz(i,j,k,box) +difp*( &
					+wrkqpz(ip,j,k,box)+wrkqpz(im,j,k,box) &
					+wrkqpz(i,jp,k,box)+wrkqpz(i,jm,k,box) &
					+wrkqpz(i,j,kp,box)+wrkqpz(i,j,km,box) &
					-6.*wrkqpz(i,j,k,box) )
				
				qpres(i,j,k,box)=wrkqpres(i,j,k,box) +dife*( &
					+ wrkqpres(ip,j,k,box)+wrkqpres(im,j,k,box) &
					+ wrkqpres(i,jp,k,box)+wrkqpres(i,jm,k,box) &
					+ wrkqpres(i,j,kp,box)+wrkqpres(i,j,km,box) &
					-6.*wrkqpres(i,j,k,box) )
				
				!	Species 2
				
				hrho(i,j,k,box)=wrkhrho(i,j,k,box) +difr*( &
					+ wrkhrho(ip,j,k,box)+wrkhrho(im,j,k,box) &
					+ wrkhrho(i,jp,k,box)+wrkhrho(i,jm,k,box) &
					+ wrkhrho(i,j,kp,box)+wrkhrho(i,j,km,box) &
					-6.*wrkhrho(i,j,k,box) )
				
				hpx(i,j,k,box)=wrkhpx(i,j,k,box) +difp*( &
					+wrkhpx(ip,j,k,box)+wrkhpx(im,j,k,box) &
					+wrkhpx(i,jp,k,box)+wrkhpx(i,jm,k,box) &
					+wrkhpx(i,j,kp,box)+wrkhpx(i,j,km,box) &
					-6.*wrkhpx(i,j,k,box) )
				
				hpy(i,j,k,box)=wrkhpy(i,j,k,box) +difp*( &
					+wrkhpy(ip,j,k,box)+wrkhpy(im,j,k,box) &
					+wrkhpy(i,jp,k,box)+wrkhpy(i,jm,k,box) &
					+wrkhpy(i,j,kp,box)+wrkhpy(i,j,km,box) &
					-6.*wrkhpy(i,j,k,box) )
				
				hpz(i,j,k,box)=wrkhpz(i,j,k,box) +difp*( &
					+wrkhpz(ip,j,k,box)+wrkhpz(im,j,k,box) &
					+wrkhpz(i,jp,k,box)+wrkhpz(i,jm,k,box) &
					+wrkhpz(i,j,kp,box)+wrkhpz(i,j,km,box) &
					-6.*wrkhpz(i,j,k,box) )
				
				hpres(i,j,k,box)=wrkhpres(i,j,k,box) +dife*( &
					+ wrkhpres(ip,j,k,box)+wrkhpres(im,j,k,box) &
					+ wrkhpres(i,jp,k,box)+wrkhpres(i,jm,k,box) &
					+ wrkhpres(i,j,kp,box)+wrkhpres(i,j,km,box) &
					-6.*wrkhpres(i,j,k,box) )
				
				!	Species 3
				
				orho(i,j,k,box)=wrkorho(i,j,k,box) +difo*( &
					+ wrkorho(ip,j,k,box)+wrkorho(im,j,k,box) &
					+ wrkorho(i,jp,k,box)+wrkorho(i,jm,k,box) &
					+ wrkorho(i,j,kp,box)+wrkorho(i,j,km,box) &
					-6.*wrkorho(i,j,k,box) )
				
				opx(i,j,k,box)=wrkopx(i,j,k,box) +difo*( &
					+wrkopx(ip,j,k,box)+wrkopx(im,j,k,box) &
					+wrkopx(i,jp,k,box)+wrkopx(i,jm,k,box) &
					+wrkopx(i,j,kp,box)+wrkopx(i,j,km,box) &
					-6.*wrkopx(i,j,k,box) )
				
				opy(i,j,k,box)=wrkopy(i,j,k,box) +difo*( &
					+wrkopy(ip,j,k,box)+wrkopy(im,j,k,box) &
					+wrkopy(i,jp,k,box)+wrkopy(i,jm,k,box) &
					+wrkopy(i,j,kp,box)+wrkopy(i,j,km,box) &
					-6.*wrkopy(i,j,k,box) )
				
				opz(i,j,k,box)=wrkopz(i,j,k,box) +difo*( &
					+wrkopz(ip,j,k,box)+wrkopz(im,j,k,box) &
					+wrkopz(i,jp,k,box)+wrkopz(i,jm,k,box) &
					+wrkopz(i,j,kp,box)+wrkopz(i,j,km,box) &
					-6.*wrkopz(i,j,k,box) )
				
				opres(i,j,k,box)=wrkopres(i,j,k,box) +dife*( &
					+ wrkopres(ip,j,k,box)+wrkopres(im,j,k,box) &
					+ wrkopres(i,jp,k,box)+wrkopres(i,jm,k,box) &
					+ wrkopres(i,j,kp,box)+wrkopres(i,j,km,box) &
					-6.*wrkopres(i,j,k,box) )
				
				epres(i,j,k,box)=wrkepres(i,j,k,box) +dife*( &
					+ wrkepres(ip,j,k,box)+wrkepres(im,j,k,box) &
					+ wrkepres(i,jp,k,box)+wrkepres(i,jm,k,box) &
					+ wrkepres(i,j,kp,box)+wrkepres(i,j,km,box) &
					-6.*wrkepres(i,j,k,box) )
				
				bx(i,j,k,box)=wrkbx(i,j,k,box)+difb*( &
					+ wrkbx(ip,j,k,box)+wrkbx(im,j,k,box) &
					+ wrkbx(i,jp,k,box)+wrkbx(i,jm,k,box) &
					+ wrkbx(i,j,kp,box)+wrkbx(i,j,km,box) &
					-6.*wrkbx(i,j,k,box) )
				
				by(i,j,k,box)=wrkby(i,j,k,box)+difb*( &
					+ wrkby(ip,j,k,box)+wrkby(im,j,k,box) &
					+ wrkby(i,jp,k,box)+wrkby(i,jm,k,box) &
					+ wrkby(i,j,kp,box)+wrkby(i,j,km,box) &
					-6.*wrkby(i,j,k,box) )
				
				bz(i,j,k,box)=wrkbz(i,j,k,box)+difb*( &
					+ wrkbz(ip,j,k,box)+wrkbz(im,j,k,box) &
					+ wrkbz(i,jp,k,box)+wrkbz(i,jm,k,box) &
					+ wrkbz(i,j,kp,box)+wrkbz(i,j,km,box) &
					-6.*wrkbz(i,j,k,box) )
			enddo
		enddo
	enddo

	return
end subroutine diffuse
