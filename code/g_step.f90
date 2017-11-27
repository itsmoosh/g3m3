subroutine step(bx,by,bz,n1,n2,n3,box,add_dip,x,y,z,ds,errin,n_grids, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !     re-calculates coords x,y,z for one step along field line.
    !     ds is step size,
    !     errin is permissible error value.
    !
    integer box
    dimension bx(n1,n2,n3),by(n1,n2,n3),bz(n1,n2,n3)
    dimension grd_xmin(n_grids),grd_xmax(n_grids), &
    grd_ymin(n_grids),grd_ymax(n_grids), &
    grd_zmin(n_grids),grd_zmax(n_grids)
    logical roc,add_dip
    !
    do 
        ds3=-ds/3.
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x,y,z,r11,r12,r13,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        if(roc) return
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+r11,y+r12,z+r13,r21,r22,r23,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+.5*(r11+r21),y+.5*(r12+r22),z+.5* &
        (r13+r23),r31,r32,r33,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+.375*(r11+3.*r31),y+.375*(r12+3.*r32 &
        ),z+.375*(r13+3.*r33),r41,r42,r43,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call intpol(bx,by,bz,n1,n2,n3,box,add_dip, &
        x+1.5*(r11-3.*r31+4.*r41),y+1.5*(r12- &
        3.*r32+4.*r42),z+1.5*(r13-3.*r33+4.*r43), &
        r51,r52,r53,ds3,roc,n_grids, &
        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        errcur=abs(r11-4.5*r31+4.*r41-.5*r51)+abs(r12-4.5*r32+4.*r42-.5* &
        r52)+abs(r13-4.5*r33+4.*r43-.5*r53)
        if (errcur.lt.errin) exit
        ds=ds*.5
        if (ds.lt.0.2) exit
    enddo
    x=x+.5*(r11+4.*r41+r51)
    y=y+.5*(r12+4.*r42+r52)
    z=z+.5*(r13+4.*r43+r53)
    if(errcur.lt.errin*.04.and.abs(ds).lt.1.33) ds=ds*1.5
    return
end
