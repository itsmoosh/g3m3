!
!	Cuts up array into a regular space sub-grid, then writes to binary
!	file 13
!
subroutine crafdatv(bx,by,bz, &
	qpx,qpy,qpz,qrho,qpresx,qpresy,qpresz,rmassq, &
	hpx,hpy,hpz,hrho,hpresx,hpresy,hpresz,rmassh, &
	opx,opy,opz,orho,opresx,opresy,opresz,rmasso, &
	epres,nx,ny,nz,n_grids,box,craft,ncraft,n,ut, &
	re_equiv,b_equiv,v_equiv,rho_equiv,gamma, &
	grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)

	integer box
	real hpressx, hpressy, hpressz
	dimension bx(nx,ny,nz),by(nx,ny,nz), &
	bz(nx,ny,nz),epres(nx,ny,nz,n_grids), &
	qpx(nx,ny,nz,n_grids),qpy(nx,ny,nz,n_grids),qpz(nx,ny,nz,n_grids), &
	qrho(nx,ny,nz,n_grids),qpresx(nx,ny,nz,n_grids), &
	qpresy(nx,ny,nz,n_grids),qpresz(nx,ny,nz,n_grids), &
	hpx(nx,ny,nz,n_grids),hpy(nx,ny,nz,n_grids),hpz(nx,ny,nz,n_grids), &
	hrho(nx,ny,nz,n_grids),hpresx(nx,ny,nz,n_grids), &
	hpresy(nx,ny,nz,n_grids),hpresz(nx,ny,nz,n_grids), &
	opx(nx,ny,nz,n_grids),opy(nx,ny,nz,n_grids),opz(nx,ny,nz,n_grids), &
	orho(nx,ny,nz,n_grids),opresx(nx,ny,nz,n_grids) , &
	opresy(nx,ny,nz,n_grids),opresz(nx,ny,nz,n_grids)
	dimension grd_xmin(n_grids),grd_xmax(n_grids), &
	grd_ymin(n_grids),grd_ymax(n_grids), &
	grd_zmin(n_grids),grd_zmax(n_grids)
	
	dimension craft(4,ncraft)
	
	!	Set gridding
	
	delx=(grd_xmax(box)-grd_xmin(box))/(nx-1.)
	dely=(grd_ymax(box)-grd_ymin(box))/(ny-1.)
	delz=(grd_zmax(box)-grd_zmin(box))/(nz-1.)
	
	!	Load t stuff
	
	az=craft(3,n)/re_equiv
	z1=1+(az-grd_zmin(box))/delz
	k1=z1
	k2=k1+1
	dz=(z1-k1)
	
	ay=craft(2,n)/re_equiv
	y1=1+(ay-grd_ymin(box))/dely
	j1=y1
	j2=j1+1
	dy=(y1-j1)
	
	ax=craft(1,n)/re_equiv
	x1=1+(ax-grd_xmin(box))/delx
	i1=x1
	i2=i1+1
	dx=(x1-i1)
	
	sbx=bx(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
		+bx(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
		+bx(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
		+bx(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
		+bx(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
		+bx(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
		+bx(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
		+bx(i2,j2,k2)*(dx)*(dy)*(dz)
	sby=by(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
		+by(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
		+by(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
		+by(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
		+by(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
		+by(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
		+by(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
		+by(i2,j2,k2)*(dx)*(dy)*(dz)
	sbz=bz(i1,j1,k1)*(1.-dx)*(1.-dy)*(1.-dz) &
		+bz(i1,j1,k2)*(1.-dx)*(1.-dy)*(dz) &
		+bz(i1,j2,k1)*(1.-dx)*(dy)*(1.-dz) &
		+bz(i1,j2,k2)*(1.-dx)*(dy)*(dz) &
		+bz(i2,j1,k1)*(dx)*(1.-dy)*(1.-dz) &
		+bz(i2,j1,k2)*(dx)*(1.-dy)*(dz) &
		+bz(i2,j2,k1)*(dx)*(dy)*(1.-dz) &
		+bz(i2,j2,k2)*(dx)*(dy)*(dz)
	
	sqpx=qpx(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+qpx(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+qpx(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+qpx(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+qpx(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+qpx(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+qpx(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+qpx(i2,j2,k2,box)*(dx)*(dy)*(dz)
	sqpy=qpy(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+qpy(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+qpy(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+qpy(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+qpy(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+qpy(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+qpy(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+qpy(i2,j2,k2,box)*(dx)*(dy)*(dz)
	sqpz=qpz(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+qpz(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+qpz(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+qpz(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+qpz(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+qpz(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+qpz(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+qpz(i2,j2,k2,box)*(dx)*(dy)*(dz)
	aqrho=qrho(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+qrho(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+qrho(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+qrho(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+qrho(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+qrho(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+qrho(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+qrho(i2,j2,k2,box)*(dx)*(dy)*(dz)
	qden=(aqrho+0.00001)/rmassq+0.0000001
	qpres1=qpresx(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+qpresx(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+qpresx(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+qpresx(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+qpresx(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+qpresx(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+qpresx(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+qpresx(i2,j2,k2,box)*(dx)*(dy)*(dz)
	qpres2=qpresy(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+qpresy(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+qpresy(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+qpresy(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+qpresy(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+qpresy(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+qpresy(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+qpresy(i2,j2,k2,box)*(dx)*(dy)*(dz)
	qpres3=qpresz(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+qpresz(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+qpresz(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+qpresz(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+qpresz(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+qpresz(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+qpresz(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+qpresz(i2,j2,k2,box)*(dx)*(dy)*(dz)
	rtq=sqrt( ( (qpres1+qpres2+qpres3)/3.)/qden)
	
	shpx=hpx(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+hpx(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+hpx(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+hpx(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+hpx(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+hpx(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+hpx(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+hpx(i2,j2,k2,box)*(dx)*(dy)*(dz)
	shpy=hpy(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+hpy(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+hpy(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+hpy(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+hpy(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+hpy(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+hpy(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+hpy(i2,j2,k2,box)*(dx)*(dy)*(dz)
	shpz=hpz(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+hpz(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+hpz(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+hpz(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+hpz(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+hpz(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+hpz(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+hpz(i2,j2,k2,box)*(dx)*(dy)*(dz)
	ahrho=hrho(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+hrho(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+hrho(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+hrho(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+hrho(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+hrho(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+hrho(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+hrho(i2,j2,k2,box)*(dx)*(dy)*(dz)
	hden=(ahrho+0.00001)/rmassh+0.0000001
	hpres1=hpresx(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+hpresx(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+hpresx(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+hpresx(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+hpresx(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+hpresx(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+hpresx(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+hpresx(i2,j2,k2,box)*(dx)*(dy)*(dz)
	hpres2=hpresy(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+hpresy(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+hpresy(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+hpresy(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+hpresy(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+hpresy(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+hpresy(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+hpresy(i2,j2,k2,box)*(dx)*(dy)*(dz)
	hpres3=hpresz(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+hpresz(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+hpresz(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+hpresz(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+hpresz(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+hpresz(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+hpresz(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+hpresz(i2,j2,k2,box)*(dx)*(dy)*(dz)
	rth=sqrt( ( (hpres1+hpres2+hpres3)/3.)/hden)
	
	sopx=opx(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+opx(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+opx(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+opx(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+opx(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+opx(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+opx(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+opx(i2,j2,k2,box)*(dx)*(dy)*(dz)
	sopy=opy(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+opy(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+opy(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+opy(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+opy(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+opy(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+opy(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+opy(i2,j2,k2,box)*(dx)*(dy)*(dz)
	sopz=opz(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+opz(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+opz(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+opz(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+opz(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+opz(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+opz(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+opz(i2,j2,k2,box)*(dx)*(dy)*(dz)
	aorho=orho(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+orho(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+orho(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+orho(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+orho(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+orho(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+orho(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+orho(i2,j2,k2,box)*(dx)*(dy)*(dz)
	oden=(aorho+0.00001)/rmasso +0.0000001
	opres1=opresx(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+opresx(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+opresx(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+opresx(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+opresx(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+opresx(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+opresx(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+opresx(i2,j2,k2,box)*(dx)*(dy)*(dz)
	opres2=opresy(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+opresy(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+opresy(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+opresy(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+opresy(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+opresy(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+opresy(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+opresy(i2,j2,k2,box)*(dx)*(dy)*(dz)
	opres3=opresz(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+opresz(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+opresz(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+opresz(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+opresz(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+opresz(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+opresz(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+opresz(i2,j2,k2,box)*(dx)*(dy)*(dz)
	rto=sqrt( ( (opres1+opres2+opres3)/3.)/oden)
	
	epress=epres(i1,j1,k1,box)*(1.-dx)*(1.-dy)*(1.-dz) &
		+epres(i1,j1,k2,box)*(1.-dx)*(1.-dy)*(dz) &
		+epres(i1,j2,k1,box)*(1.-dx)*(dy)*(1.-dz) &
		+epres(i1,j2,k2,box)*(1.-dx)*(dy)*(dz) &
		+epres(i2,j1,k1,box)*(dx)*(1.-dy)*(1.-dz) &
		+epres(i2,j1,k2,box)*(dx)*(1.-dy)*(dz) &
		+epres(i2,j2,k1,box)*(dx)*(dy)*(1.-dz) &
		+epres(i2,j2,k2,box)*(dx)*(dy)*(dz)
	rte=sqrt(epress/(qden+oden+hden))
	vte=42.8*rte*v_equiv
	
	!	Put back in real units and GSM
	
	ax=-ax*re_equiv
	ay=-ay*re_equiv
	az=az*re_equiv
	sbx=-sbx*b_equiv
	sby=-sby*b_equiv
	sbz=sbz*b_equiv
	btot=sqrt(sbx**2+sby**2+sbz**2)
	
	vxq=-sqpx*v_equiv/(aqrho+0.00000001)
	vyq=-sqpy*v_equiv/(aqrho+0.00000001)
	vzq=sqpz*v_equiv/(aqrho+0.00000001)
	vtq=rtq/sqrt(rmassq)*v_equiv
	vtotq=sqrt(vxq**2+vyq**2+vzq**2)
	
	vxh=-shpx*v_equiv/(ahrho+0.00000001)
	vyh=-shpy*v_equiv/(ahrho+0.00000001)
	vzh=shpz*v_equiv/(ahrho+0.00000001)
	vth=rth/sqrt(rmassh)*v_equiv
	vtoth=sqrt(vxh**2+vyh**2+vzh**2)
	
	vxo=-sopx*v_equiv/(aorho+0.00000001)
	vyo=-sopy*v_equiv/(aorho+0.00000001)
	vzo=sopz*v_equiv/(aorho+0.00000001)
	vto=rto/sqrt(rmasso)*v_equiv
	vtoto=sqrt(vxo**2+vyo**2+vzo**2)
	
	qden=qden*rho_equiv
	hden=hden*rho_equiv
	oden=oden*rho_equiv
	rden=qden/(qden+oden+hden)
	xden=oden/(qden+oden+hden)
	
	rtx=hpres1/(hpres1+hpres2+hpres3+1.e-8)
	rty=hpres2/(hpres1+hpres2+hpres3+1.e-8)
	rtz=hpres3/(hpres1+hpres2+hpres3+1.e-8)
	
	mout=50+n
	write(mout,*)ut,ax,ay,az
	write(mout,*)btot,sbx,sby,sbz
	write(mout,*)qden,vxq,vyq,vzq
	write(mout,*)hden,vxh,vyh,vzh
	write(mout,*)oden,vxo,vyo,vzo
	write(mout,*)rden,vte,vtq,vto
	write(mout,*)xden,rtx,rty,rtz
	
	return
end subroutine crafdatv
