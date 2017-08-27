subroutine trilin_interp(sxyz,gridpts,qty,sdata)
    !
    !   Performs trilinear interpolation to approximate the value
	!		of the measurable quantity qty at the spatial coordinate
	!		of the spacecraft sxyz using the nearest xyz above and below
	!		the spacecraft location. Grid is assumed regular, so gridpts
	!		is only two ordered triplets: (x0,y0,z0) and (x1,y1,z1).
	!	Only accepts scalar input. Call for each component of measured
	!		vectors.
    !
	real, intent(in) :: sxyz(3), gridpts(3,2), qty(2,2,2)
	real, intent(out) :: sdata
	!
	real x0, x1, y0, y1, z0, z1
	!
	!	We interpolate between corners to approximate x-parallel edge values
	real edge00, edge01, edge10, edge11
	!	Then interpolate between edge values to find xy-plane face values
	real face0, face1
	!
	real xd, yd, zd
	!
	!	**********************************
	!
	x0=gridpts(1,1)
	x1=gridpts(1,2)
	y0=gridpts(2,1)
	y1=gridpts(2,2)
 	z0=gridpts(3,1)
	z1=gridpts(3,2)
	!
	!	Find % of the way from pt. 0 to pt. 1 for each axis
	xd = (sxyz(1) - x0)/(x1 - x0)
	yd = (sxyz(2) - y0)/(y1 - y0)
	zd = (sxyz(3) - z0)/(z1 - z0)
	!
	!	Find interpolated value at the right % from low-x face to high-x face
	!		along each of the x-parallel edges
	edge00 = qty(1,1,1)*(1-xd) + qty(2,1,1)*xd
	edge01 = qty(1,1,2)*(1-xd) + qty(2,1,2)*xd
	edge10 = qty(1,2,1)*(1-xd) + qty(2,2,1)*xd
	edge11 = qty(1,2,2)*(1-xd) + qty(2,2,2)*xd
	!
	!	Use interpolated edge values to interpolate along y-parallel
	!		direction to find values above and below the desired point in z
	face0 = edge00*(1-yd) + edge10*yd
	face1 = edge01*(1-yd) + edge11*yd
	!
	!	Interpolate between face values to get trilinearly interpolated
	!		value for qty at sxyz location
	sdata = face0*(1-zd) + face1*zd
	!
    return
end subroutine trilin_interp
