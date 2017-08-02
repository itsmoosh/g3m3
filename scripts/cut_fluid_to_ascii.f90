
!     this is a 3-d modified three fluid simulation using the
!           electrons :  arrays starting with e
!           solar wind : arrays starting with q
!           ionospheric: arrays starting with o oxygen,
!                                         h hydrogen
!
!      and change set33(0.1,1.,0.1,1.0  to
!                 set3(0.,1.,0.,1.
!     warning - make sure space is compatible with graphics
!               routines
!               the arrays ijzero,ijmid  and ijsrf have to
!               be modified in bsubs.f if modified in main
!     add in mirror dipole  so bx is zero at wind boundary
!     graphics - contour and flows - have to be manually set
!                for right aspect ratios
!
!      warming plasma and magnetic field data must be aligned in time
!
program cut_fluid_to_ascii

      implicit none
      integer, parameter :: nx=121,ny=121,nz=61,ngrd=5, &
                  mbndry=1,msrf=2000,mmid=1500,mzero=5000, &
                  ncraft=30,ncts=281

!     --------------------------------------------------------------
      integer, parameter :: nlines=479 !synthetic trajectory length
!      integer, parameter :: nlines2=6302
      integer, parameter :: skip=1 !sample skip for flythrough data
!      integer, parameter :: nlines=404 !mex trajectory length
!      integer, parameter :: skip=5 !sample skip for flythrough data
      integer, parameter :: vbins=100 !number of energy/velocity bins on y-axis
      integer, parameter :: phi_res=10 !number of angle bins on y-axis

       real, parameter :: ev = 1.602e-19, &
                          mu_o = 1.2566e-06, &
                          mp = 1.673e-27, &
                          radius = 6.0268e06, &
                          mr = 1836.0, &
                          pi = 3.141592654, &
                          norm = 609082663.5

!      norm  = (1/(2pi)^3/2)*(10^6)*(1.602x10^-19)*(10^-4)/mp
!                            -----   --------------------
!                               |          |
!      convert num dens from cm^-3 to m^-3 |
!      convert 1/s/m^2/j to 1/s/cm^2/ev ---|
!
!     radius = planetary radius in meters
!     mr = mp/me
!
      real rx,ry,rz
!
!      grid limits now set by grd_min grd_max arrays
!      ncore denotes couser grid to be hollowed out by fine grid
!      nbndry denotes finer grid to which coaser grid sets flanks
!      xspac is the relative grid spacing relative to inner grid system
!      main_ngrd gives the box number from which the fine gridding is
!                   interpolated
       real grd_xmin(ngrd),grd_xmax(ngrd),grd_ymin(ngrd), &
           grd_ymax(ngrd), grd_zmin(ngrd),grd_zmax(ngrd), &
           xspac(ngrd)
       integer ncore(ngrd),nbndry(ngrd)
!

!      xcraft is the actual position of the spacecraft in re
!          4th dimension of the actual time
!      zcraft is the future position of the spacecraft in re
!      rcraft is the position of the spacecraft for which
!           imf is reference. no alteration from boundary conditions applied
!
      real xcraft(4,ncraft),zcraft(4,ncraft),rcraft(3)

!     physics plasma quantities
      real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
           qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
           qrho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
           opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
           orho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
           hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
           hrho(nx,ny,nz,ngrd),epres(nx,ny,nz,ngrd)
!
      real qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
           qpresz(nx,ny,nz,ngrd),qpresxy(nx,ny,nz,ngrd), &
           qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
           qpara(nx,ny,nz,ngrd),qperp(nx,ny,nz,ngrd), &
           qcross(nx,ny,nz,ngrd)

      real hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
           hpresz(nx,ny,nz,ngrd),hpresxy(nx,ny,nz,ngrd), &
           hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
           hpara(nx,ny,nz,ngrd),hperp(nx,ny,nz,ngrd), &
           hcross(nx,ny,nz,ngrd)

      real opresx(nx,ny,nz,ngrd),opresy(nx,ny,nz,ngrd), &
           opresz(nx,ny,nz,ngrd),opresxy(nx,ny,nz,ngrd), &
           opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
           opara(nx,ny,nz,ngrd),operp(nx,ny,nz,ngrd), &
           ocross(nx,ny,nz,ngrd)

      real qvx(nx,ny,nz,ngrd), qvy(nx,ny,nz,ngrd), &
           qvz(nx,ny,nz,ngrd), ovx(nx,ny,nz,ngrd), &
           ovy(nx,ny,nz,ngrd), ovz(nx,ny,nz,ngrd), &
           hvx(nx,ny,nz,ngrd), hvy(nx,ny,nz,ngrd), &
           hvz(nx,ny,nz,ngrd)


      real evelx(nx,ny,nz,ngrd), evely(nx,ny,nz,ngrd), &
           evelz(nx,ny,nz,ngrd)

      real evx(nx,ny,nz), evy(nx,ny,nz), evz(nx,ny,nz)
      real vvx(nx,ny,nz), vvy(nx,ny,nz), vvz(nx,ny,nz)

      real efldx(nx,ny,nz,ngrd), efldy(nx,ny,nz,ngrd), &
           efldz(nx,ny,nz,ngrd)

      real ex(nx,ny,nz), ey(nx,ny,nz), ez(nx,ny,nz)
      real tvx(nx,ny,nz), tvy(nx,ny,nz), tvz(nx,ny,nz)

!     unperturbed quantities
!
      real bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
      real bxs(nx,ny,nz,ngrd),bys(nx,ny,nz,ngrd),bzs(nx,ny,nz,ngrd)
      real br0(nx,ny,nz,ngrd),bt0(nx,ny,nz,ngrd),bp0(nx,ny,nz,ngrd)
      real bxt(nx,ny,nz),byt(nx,ny,nz),bzt(nx,ny,nz),btot(nx,ny,nz)


      real curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
           resistive(nx,ny,nz),curx_all(nx,ny,nz,ngrd), &
           cury_all(nx,ny,nz,ngrd),curz_all(nx,ny,nz,ngrd)



!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
      real nq(nx,ny,nz,ngrd), no(nx,ny,nz,ngrd), &
           nh(nx,ny,nz,ngrd), ne(nx,ny,nz,ngrd), &
           vte(nx,ny,nz,ngrd), vbe(nx,ny,nz,ngrd), &
           vtq(nx,ny,nz,ngrd), vto(nx,ny,nz,ngrd), &
           vth(nx,ny,nz,ngrd), vbq(nx,ny,nz,ngrd), &
           vbo(nx,ny,nz,ngrd), vbh(nx,ny,nz,ngrd), &
           vbq_para(nx,ny,nz,ngrd),vbq_perp(nx,ny,nz,ngrd), &
           vbh_para(nx,ny,nz,ngrd),vbh_perp(nx,ny,nz,ngrd), &
           vbo_para(nx,ny,nz,ngrd),vbo_perp(nx,ny,nz,ngrd)

      real bxtmp, bytmp, bztmp, bmag, dotq, cosq, &
           doth, cosh, doto, coso, para, perp, &
           vq_flat, vo_flat, vh_flat, ve_flat, &
           vbq_para_intp, vbq_perp_intp, vbh_para_intp, &
           vbh_perp_intp, vbo_para_intp, vbo_perp_intp, &
           vtq_intp, vth_intp, vto_intp, &
           numq_intp, numh_intp, numo_intp, &
           vbe_intp, vte_intp, nume_intp, &
           thetaq,thetah,thetao,drftnorm, &
           qeratio, oeratio, heratio

      real theta,q_costheta,o_costheta,h_costheta,sintheta
      real s_xmin(ngrd),s_xmax(ngrd),s_ymin(ngrd), &
           s_ymax(ngrd),s_zmin(ngrd),s_zmax(ngrd), &
           r_xmin(ngrd), r_xmax(ngrd), r_ymin(ngrd), &
           r_ymax(ngrd), r_zmin(ngrd), r_zmax(ngrd)

      real nrgy,mq,mo,mh,lognrgy,kreal,vbinsreal

      real ml(nlines),mlt(nlines),rad_sat(nlines), &
           xsat(nlines),ysat(nlines),zsat(nlines), &
           xplot_f, yplot_f


      real,allocatable :: gx(:),gy(:),gz(:),gr(:), &
           gbx(:),gby(:),gbz(:),xpos(:),ypos(:),zpos(:), &
           inbox(:,:)

      real xtemp, ytemp, ztemp

      real phi, vqpar, vqprp, vopar, voprp, vhpar, vhprp

      integer, allocatable :: nbox(:), grd_sys(:)

      integer a,b,i,j,k,m,n,posit, boxtemp, count1, count2, &
           count3, count4, count5, count6, count7, count8, &
           count9, moon
      integer ix, iy, iz, i_xmin, i_xmax, i_ymin, i_ymax, &
           i_zmin, i_zmax, nchf

      real aqpres, ahpres, aopres, aepres, qden, hden, oden, &
           eden, qtemp, htemp, otemp, etemp, qpress, hpress, &
           opress, epress, qdens, hdens, odens, edens, &
           qtemps, htemps, otemps, etemps, abx, aby, abz, &
           ri, rj, rk, rad, bsurmag, time, abx2, aby2, abz2, &
           ri_moon, rj_moon, rk_moon,aoden,aqden,ahden,aeden, &
           aotemp,ahtemp,aqtemp,aetemp

      integer nx1, nx2, ny1, ny2, nz1, nz2

      character*1 cut,box

      character*32 wd1,wd2,wd3,wd4,wd5,wd6,wd7,wd8,wd9

      real t, planet_rad, moon_rad, pl_ratio, t_equiv, ut, &
           ofrac, frac_o,dist,tot_v,rot_mach,j_rad,j_phi, &
           cur_norm, planet_per

      real,allocatable :: yplot(:,:),flux_p(:,:), &
           flux_q(:,:),flux_o(:,:),flux_h(:,:), &
           flux_e(:,:), bxg(:), byg(:), bzg(:), &
           brg(:), btg(:), bpg(:),cxg(:), cyg(:), czg(:)

!~*~*~*~*~*~**~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~

      integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
              ijzero(mbndry,3,mzero)
      integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
      real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
              parm_zero(mbndry,7,mzero)
!

      logical start,add_dip,ringo,update,save_dat,write_dat, &
              spacecraft,tilting,warp,reload, &
              divb_lores,divb_hires,divb_on

      real tmax,stepsz,tsave,xdip,ydip,zdip,rearth, &
           tilt1,tilt2,rmassq,rmassh,rmasso,cs_inner, &
           alf_inner1,alf_inner2,alpha_e,den_earth, &
           o_conc,gravity,ti_te,gamma, &
           re_wind,cs_wind,vx_wind1,vx_wind2, &
           vy_wind1,vy_wind2,vz_wind1,vz_wind2, &
           alfx_wind1,alfx_wind2, &
           alfy_wind1,alfy_wind2, &
           alfz_wind1,alfz_wind2, &
           den_wind1,den_wind2, &
           reynolds,resist,rho_frac,bfrac,vfrac, &
           re_equiv,b_equiv,v_equiv,rho_equiv, &
           utstart,chirho,chipxyz,chierg, &
           difrho,difpxyz,diferg

      real xmoon_re,ymoon_re,zmoon_re,cs_moon, &
           qden_moon,hden_moon,oden_moon,ti_te_moon, &
           alf_moon,xdip_moon,ydip_moon,zdip_moon,offset, &
           orbit_moon,theta_moon,ani,resist_moon,den_lunar, &
           reduct,t_torus,aniso_factor,ani_q,ani_h,ani_o,uday

      integer no_shock,ntgraf

      logical isotropic,flow_reset

!
!      ringo decides if you you want to plot 1 set of diagnostics
!              with no time stepping
!      update decides if time setting is to be reset
!
!     rho,pres,erg,px,py are the density, pressure, energy and momentum
!              in the x and y directions, respectively, of the fluid
!       the indices 1, 2 is required to store old and new values
!
!     ijsrf give position of ionosphere in grid units - plasma
!           parameters stored in parm_srf
!     ijmid gives intermediate boundary - say representing the
!            atmosphere - plasma parameters stored in parm_mid
!     ijzero gives position of all grid units interior to surface
!
!     frho,ferg and fpx,fpy are the estimates of the fluid quantities
!           at n+1/2 as defined by the lax-wendroff scheme
!       the index of 3 is needed to store adjacent x values for these fns
!
!     d_min is the minimum allowable density
!     stepsz is the size of the time step in terms of delta_x,y/(|umax|+c_s)
!
!     system dimensions are nx,ny,nz
!     variable grid spacing enabled with rxyz >1
!           rx,ry,rz should not be set identically to zero
!
        namelist/option/tmax,ntgraf,stepsz,start,tsave,isotropic
        namelist/earth/xdip,ydip,zdip,rearth, &
            tilt1,tilt2,tilting,rmassq,rmassh,rmasso
        namelist/speeds/cs_inner,alf_inner1,alf_inner2, &
            alpha_e,den_earth,den_lunar,o_conc,gravity, &
            ti_te,gamma,ringo,update,reload, &
            divb_lores,divb_hires,reduct,t_torus,aniso_factor, &
            ani_q,ani_h,ani_o
        namelist/windy/re_wind,cs_wind,vx_wind1,vx_wind2, &
            vy_wind1,vy_wind2,vz_wind1,vz_wind2, &
            alfx_wind1,alfx_wind2, &
            alfy_wind1,alfy_wind2, &
            alfz_wind1,alfz_wind2, &
            den_wind1,den_wind2, &
            reynolds,resist,rho_frac,bfrac,vfrac
        namelist/lunar/orbit_moon,theta_moon,cs_moon, &
            qden_moon,hden_moon,oden_moon, &
            alf_moon,ti_te_moon, &
            xdip_moon,ydip_moon,zdip_moon,offset
        namelist/physical/re_equiv,b_equiv,v_equiv,rho_equiv, &
            spacecraft,warp,uday,utstart
        namelist/smooth/chirho,chipxyz,chierg, &
            difrho,difpxyz,diferg

      cut='y'


!      open input data file. Must be in same directory as executable.

!*******************************************************************

      open(5,file='cutter_input',status='old',form='formatted')
      open(15,file='timebox.dat',status='unknown',form='formatted')

!     read input parameters
!
      read(5,option)
      read(5,earth)
      read(5,speeds)
      read(5,windy)
      read(5,lunar)
      read(5,physical)
      read(5,smooth)
!
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
!     large grid system
!        ncore -> grid gets all the information from this grid no.
!        nbdry -> which grid data will be applied to this grid no.for bndry
!
      do m=1,ngrd
         read(5,*)grd_xmin(m),grd_xmax(m),grd_ymin(m),grd_ymax(m), &
              grd_zmin(m),grd_zmax(m),xspac(m)
         write(6,*)grd_xmin(m),grd_xmax(m),grd_ymin(m),grd_ymax(m), &
              grd_zmin(m),grd_zmax(m),xspac(m)
         ix=1+int((grd_xmax(m)-grd_xmin(m))/xspac(m))
         iy=1+int((grd_ymax(m)-grd_ymin(m))/xspac(m))
         iz=1+int((grd_zmax(m)-grd_zmin(m))/xspac(m))
         if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz))then
            write(6,*)' warning: sizes',ix,iy,iz,nx,ny,nz
            stop
         endif
      enddo
!
!
!!!!

!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
! read in data file to cut 
!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

      print *,'begin reading data files'

      ! input fluid file 'cutter_fluid'
      open(12,file='cutter_fluid',status='unknown',form='unformatted')
      nchf=12
!
       read(nchf)t

       read(nchf)qrho
       read(nchf)qpx
       read(nchf)qpy
       read(nchf)qpz
       read(nchf)qpresx
       read(nchf)qpresy
       read(nchf)qpresz
       read(nchf)qpresxy
       read(nchf)qpresxz
       read(nchf)qpresyz

       read(nchf)hrho
       read(nchf)hpx
       read(nchf)hpy
       read(nchf)hpz
       read(nchf)hpresx
       read(nchf)hpresy
       read(nchf)hpresz
       read(nchf)hpresxy
       read(nchf)hpresxz
       read(nchf)hpresyz

       read(nchf)orho
       read(nchf)opx
       read(nchf)opy
       read(nchf)opz
       read(nchf)opresx
       read(nchf)opresy
       read(nchf)opresz
       read(nchf)opresxy
       read(nchf)opresxz
       read(nchf)opresyz

       read(nchf)bx
       read(nchf)by
       read(nchf)bz
       read(nchf)epres
       read(nchf)bx0
       read(nchf)by0
       read(nchf)bz0
       read(nchf)parm_srf,parm_mid,parm_zero, &
                ijzero,numzero,ijmid,nummid,ijsrf,numsrf
       close(12)

!
      moon_rad = 2575.   !km
      planet_rad = 60268.   !km
      planet_per = 10.7  !hours
      t_equiv=planet_rad*re_equiv/v_equiv

      cur_norm = 1.37e-9
!      t_equiv=(radius/1000.0)*re_equiv/v_equiv


  170 ut=utstart+t*t_equiv/3600.
      print*,'time = ',ut,'ut'
      write(15,*)ut
      close(15)

      close(11)

!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
! calculate fields/currents/velocities from data file
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
      do m = 1, 3
         rx=xspac(m)
         ry=xspac(m)
         rz=xspac(m)

!     calculate current

         call calcur(bx,by,bz,nx,ny,nz,ngrd,m, &
              curx,cury,curz,rx,ry,rz)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  curx_all(i,j,k,m) = curx(i,j,k)
                  cury_all(i,j,k,m) = cury(i,j,k)
                  curz_all(i,j,k,m) = curz(i,j,k)

               end do
            end do
         end do


!     write(6,*)' totbfld ing now'
      call totfld(bx,bx0,bxt,nx,ny,nz,ngrd,m)
      call totfld(by,by0,byt,nx,ny,nz,ngrd,m)
      call totfld(bz,bz0,bzt,nx,ny,nz,ngrd,m)
!
!     find magnitude of b
!
      call tot_b(btot,bxt,byt,bzt,nx,ny,nz)

      call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz, &
            evx,evy,evz,tvx,tvy,tvz, &
            nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,reynolds)


!     write(6,*)' bande ing now' 
         call bande(ex,ey,ez,bxt,byt,bzt, &
              curx,cury,curz,evx,evy,evz,btot, &
              epres,qrho,hrho,orho,resistive,resist,reynolds, &
              nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
              ijmid,nummid,ijzero,numzero,mbndry,mmid,mzero, &
              rx,ry,rz)

         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  efldx(i,j,k,m) = ex(i,j,k)/9.26
                  efldy(i,j,k,m) = ey(i,j,k)/9.26 ! drift speed normalization
                  efldz(i,j,k,m) = ez(i,j,k)/9.26

               end do
            end do
         end do


!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~
! below here, you can calculate specifics 
! you want out, like pressure anisotropy 
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~

!   EXAMPLE: 
!     find anisotropy values for species 1
!
      call fnd_vel(qpx,qpy,qpz,qrho,vvx,vvy,vvz,nx,ny,nz,ngrd,m)

        do k=1,nz
           do j=1,ny
              do i=1,nx
                 qvx(i,j,k,m) = vvx(i,j,k)
                 qvy(i,j,k,m) = vvy(i,j,k)
                 qvz(i,j,k,m) = vvz(i,j,k)
              end do
           end do
        end do
!
!
!      write(*,*) "finding qpres-s now!"
      call fnd_pres(qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz,&
                    qpara,qcross,qperp, &
                    vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,ngrd,m)


      end do !loop over m

!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
      mq=rmassq
      mo=rmasso
      mh=rmassh
!~*~*~*~*~*~**~*~*~*~*~*~**~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*

      if(cut == "y") then
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
! here, write out files for DATA explorer (ascii data)
! you can also just output as needed for matplotlib, etc.
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
      print *,'making ascii data files'

!     output data file names
!     ------------

       do m=1,3
           write(box,'(i1)') m
           wd1 = "plasma.dat"//trim(adjustl(box))
           wd2 = "flow.dat"//trim(adjustl(box))
           wd3 = "field.dat"//trim(adjustl(box))
           wd4 = "rad.dat"//trim(adjustl(box))
           wd5 = "model.dat"//trim(adjustl(box))
           wd6 = "efield.dat"//trim(adjustl(box))
           wd7 = "pressure.dat"//trim(adjustl(box))
           wd8 = "rotational.dat"//trim(adjustl(box))
           wd9 = "currents.dat"//trim(adjustl(box))
           open(7,file=wd1,status='replace',form='formatted')
           open(8,file=wd2,status='replace',form='formatted')
           open(9,file=wd3,status='replace',form='formatted')
           open(10,file=wd4,status='replace',form='formatted')
           open(11,file=wd5,status='replace',form='formatted')
           open(12,file=wd6,status='replace',form='formatted')
           open(13,file=wd7,status='replace',form='formatted')
           open(14,file=wd8,status='replace',form='formatted')
           open(15,file=wd9,status='replace',form='formatted')


          ! output loop, put data in physical units, 
          ! and essentially write out what you are interested
          ! in plotting
          do k=1,nz
             do j=1,ny
                do i=1,nx
                   ! normalized  pressures
                   aqpres=amax1(0.0000001,0.333*(qpresx(i,j,k,m)+ &
                            qpresy(i,j,k,m)+qpresz(i,j,k,m)))
                   ahpres=amax1(0.0000001,0.333*(hpresx(i,j,k,m)+ &
                            hpresy(i,j,k,m)+hpresz(i,j,k,m)))
                   aopres=amax1(0.0000001,0.333*(opresx(i,j,k,m)+ &
                            opresy(i,j,k,m)+opresz(i,j,k,m)))
                   aepres=amax1(0.0000001,epres(i,j,k,m))
    
                   ! normalized densities
                   qden=amax1(qrho(i,j,k,m),0.000001)
                   hden=amax1(0.0001,hrho(i,j,k,m))
                   oden=amax1(0.0001,orho(i,j,k,m))
                   qden=qden/rmassq
                   hden=hden/rmassh
                   oden=oden/rmasso
                   eden=oden+hden+qden+0.0001

                   ! get temps from pressures/densities
                   qtemp=amax1(0.00001,aqpres/qden)
                   htemp=ahpres/hden
                   otemp=aopres/oden
                   etemp=aepres/eden

                   qdens=alog10(qden)
                   hdens=alog10(hden)
                   odens=alog10(oden)
                   edens=alog10(eden)

                   ! total b field
                   abx=b_equiv*(bx(i,j,k,m) + bx0(i,j,k,m))
                   aby=b_equiv*(by(i,j,k,m) + by0(i,j,k,m))
                   abz=b_equiv*(bz(i,j,k,m) + bz0(i,j,k,m))

                   ! radial distance component
                   ri = grd_xmin(m) + (xspac(m)*real(i-1))
                   ri = ri*re_equiv
                   rj = grd_ymin(m) + (xspac(m)*real(j-1))
                   rj = rj*re_equiv
                   rk = grd_zmin(m) + (xspac(m)*real(k-1))
                   rk = rk*re_equiv

                   rad = sqrt((ri*ri) + (rj*rj) + (rk*rk))

                   ! k == 31 implied equatorial planet in my sims

                   if(k == 31) then
                       j_rad = cur_norm*(curx(i,j,k)*ri + &
                           cury(i,j,k)*rj)/rad

                       j_phi = cur_norm*sqrt(curx(i,j,k)**2 + &
                           cury(i,j,k)**2) - j_rad

                       dist = sqrt(ri**2+rj**2)
                       tot_v = (sqrt(qvx(i,j,k,m)**2+qvy(i,j,k,m)**2)-&
                               (qvx(i,j,k,m)*ri+qvy(i,j,k,m)*rj)/dist) + &
                               (sqrt(hvx(i,j,k,m)**2+hvy(i,j,k,m)**2)-&
                               (hvx(i,j,k,m)*ri+hvy(i,j,k,m)*rj)/dist) + &
                               (sqrt(ovx(i,j,k,m)**2+ovy(i,j,k,m)**2)-&
                               (ovx(i,j,k,m)*ri+ovy(i,j,k,m)*rj)/dist)

                       rot_mach = tot_v / &
                           ( 2*pi*rad*planet_rad/(planet_per*3600) )

                       if( rad == 0.0 ) rot_mach=0.0
                   
                       write(14,'(3(f9.2),10(es14.6))')ri,rj,rot_mach
                       write(15,'(2(f9.2),2(es14.6))')ri,rj,j_rad,j_phi
                   endif


                   bsurmag = sqrt((bx0(i,j,k,m)**2) &
                          + (by0(i,j,k,m)**2) &
                          + (bz0(i,j,k,m)**2))
                   bsurmag = bsurmag*b_equiv

                   write(7,'(3(f9.2),8(es14.6))')ri,rj,rk,&
                            qdens,qtemp,hdens,htemp, &
                            odens,otemp,edens,etemp

                   write(8,'(3(f9.2),9(es14.6))')ri,rj,rk, &
                            qvx(i,j,k,m),qvy(i,j,k,m), &
                            qvz(i,j,k,m),ovx(i,j,k,m),ovy(i,j,k,m), &
                            ovz(i,j,k,m),hvx(i,j,k,m),hvy(i,j,k,m), &
                            hvz(i,j,k,m)

                   write(9,'(3(f9.2),6(es14.6))')ri,rj,rk, &
                            abx,aby,abz,curx_all(i,j,k,m), &
                            cury_all(i,j,k,m),curz_all(i,j,k,m)

                   write(10,'(4(es14.6))')rad, bsurmag, ri, rj, rk

                   abx=b_equiv*(bx(i,j,k,m))
                   abx2=b_equiv*(bx0(i,j,k,m))
                   aby=b_equiv*(by(i,j,k,m))
                   aby2=b_equiv*(by0(i,j,k,m))
                   abz=b_equiv*(bz(i,j,k,m))
                   abz2=b_equiv*(bz0(i,j,k,m))

                   write(11,'(3(f9.2),6(es14.6))')ri,rj,rk,&
                            abx,aby,abz,abx2,aby2,abz2

                   write(12,'(3(f9.2),3(es14.6))')ri,rj,rk, &
                            efldx(i,j,k,m),efldy(i,j,k,m),efldz(i,j,k,m)

                   write(13,'(3(f9.2),10(es14.6))')ri,rj,rk, &
                            qpara(i,j,k,m),qperp(i,j,k,m), &
                            qcross(i,j,k,m)

                enddo
             enddo
          enddo

        close(7)
        close(8)
        close(9)
        close(10)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)


      end do ! loop over m

      endif
!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~**~*~
      end


!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~**~*~
!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~**~*~
!     *************************************************
      subroutine interpolate(plas_int,plas,box,xpos,ypos,zpos, &
          nbx,nx,ny,nz)

      integer, intent (in) :: nbx, nx, ny, nz, box

      real, intent (in) :: plas(nbx,nx,ny,nz), xpos, ypos, zpos
      real, intent (inout) :: plas_int

      integer i, j, k, m, ii, jj, kk
      real dx, dy, dz


      m = box

      i = int(xpos)
      dx = xpos - real(i)

      j = int(ypos)
      dy = ypos - real(j)

      k = int(zpos)
      dz = zpos - real(k)

      ii = i + 1
      jj = j + 1
      kk = k + 1

      plas_int = plas(i,j,k,m)*(1.0 - dx)*(1.0 - dy)*(1.0 - dz) &
               + plas(ii,jj,kk,m)*dx*dy*dz &
               + plas(ii,j,k,m)*dx*(1.0 - dy)*(1.0 - dz) &
               + plas(i,jj,k,m)*(1.0 - dx)*dy*(1.0 - dz) &
               + plas(i,j,k,mk)*(1.0 - dx)*(1.0 - dy)*dz &
               + plas(ii,jj,k,m)*dx*dy*(1.0 - dz) &
               + plas(ii,j,kk,m)*dx*(1.0 - dy)*dz &
               + plas(i,jj,kk,m)*(1.0 - dx)*dy*dz

      end subroutine interpolate
!     ************************************************
      subroutine set_gridpt(xpos,ypos,zpos, &
          xsat,ysat,zsat,r_xmin,r_ymin,r_zmin,xspac,re_equiv, &
          moon,xmoon_re,ymoon_re,zmoon_re)

      real, intent (inout) :: xpos,ypos,zpos

      real, intent (in) :: xsat,ysat,zsat,r_xmin, &
           r_ymin,r_zmin,xspac,re_equiv, &
           xmoon_re,ymoon_re,zmoon_re

      integer, intent (in) :: moon

      real xtmp, ytmp, ztmp, df1, df2

         df2 = xspac*re_equiv

         df1 = r_xmin - xsat
         xtmp = df1/df2
         xpos = abs(xtmp) + 1

         df1 = r_ymin - ysat
         ytmp = df1/df2
         ypos = abs(ytmp) + 1

         df1 = r_zmin - zsat
         ztmp = df1/df2
         zpos = abs(ztmp) + 1


      end subroutine set_gridpt

!     *************************************************
      subroutine calcur(bx,by,bz,nx,ny,nz,ngrd,m,curx,cury,curz, &
                    rx,ry,rz)
!
!     this calculates the current associated with the perturbed b field
!          i.e. curz= dby/dx - dbx/dy
!
      dimension bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
                curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz)
!
      dxt=2.*rx
      dyt=2.*ry
      dzt=2.*rz
!
!     parallelizes loop. RW, oct. 23, 2002
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
!
               curx(i,j,k)=(bz(i,jp,k,m)-bz(i,jm,k,m))/dyt &
                    - (by(i,j,kp,m)-by(i,j,km,m))/dzt
!
               cury(i,j,k)=(bx(i,j,kp,m)-bx(i,j,km,m))/dzt &
                    - (bz(ip,j,k,m)-bz(im,j,k,m))/dxt
!
               curz(i,j,k)=(by(ip,j,k,m)-by(im,j,k,m))/dxt &
                    - (bx(i,jp,k,m)-bx(i,jm,k,m))/dyt
            enddo
         enddo
      enddo
!
!     following boundary conditions are set so that no forward
!     communication , i.e. same x required
!     and symmetry between k=1 and k=2
!
!     set boundary regions - bottom and top panels
!
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
!     parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do j=2,ny1
         do i=2,nx1
!
!           regular boundary conditions
!
            curx(i,j,1)=curx(i,j,2)
            curx(i,j,nz)=curx(i,j,nz1)
!
            cury(i,j,1)=cury(i,j,2)
            cury(i,j,nz)=cury(i,j,nz1)
!
            curz(i,j,1)=curz(i,j,2)
            curz(i,j,nz)=curz(i,j,nz1)
         enddo
      enddo
!
!     set boundary regions - front and back
!
!     parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do k=2,nz1
         do j=2,ny1
            curx(1,j,k)=curx(2,j,k)
            curx(nx,j,k)=curx(nx1,j,k)
!
            cury(1,j,k)=cury(2,j,k)
            cury(nx,j,k)=cury(nx1,j,k)
!
            curz(1,j,k)=curz(2,j,k)
            curz(nx,j,k)=curz(nx1,j,k)
         enddo
      enddo
!
!     set boundary regions - left and right
!
!     parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do j=2,nz1
         do i=2,nx1
            curx(i,1,j)=curx(i,2,j)
            curx(i,ny,j)=curx(i,ny1,j)
!
            cury(i,1,j)=cury(i,2,j)
            cury(i,ny,j)=cury(i,ny1,j)
!
            curz(i,1,j)=curz(i,2,j)
            curz(i,ny,j)=curz(i,ny1,j)
         enddo
      enddo
!
!     set corner lines
!
!     parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do i=2,nx1
         curx(i,1,1)=curx(i,1,2)
         curx(i,1,nz)=curx(i,1,nz1)
         cury(i,1,1)=cury(i,1,2)
         cury(i,1,nz)=cury(i,1,nz1)
         curz(i,1,1)=curz(i,1,2)
         curz(i,1,nz)=curz(i,1,nz1)
!
         curx(i,ny,1)=curx(i,ny,2)
         curx(i,ny,nz)=curx(i,ny1,nz)
         cury(i,ny,1)=cury(i,ny,2)
         cury(i,ny,nz)=cury(i,ny1,nz)
         curz(i,ny,1)=curz(i,ny,2)
         curz(i,ny,nz)=curz(i,ny1,nz)
      enddo
!
!     parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do j=2,ny1
         curx(1,j,1)=curx(1,j,2)
         curx(1,j,nz)=curx(1,j,nz1)
         cury(1,j,1)=cury(1,j,2)
         cury(1,j,nz)=cury(1,j,nz1)
         curz(1,j,1)=curz(1,j,2)
         curz(1,j,nz)=curz(1,j,nz1)
!
         curx(nx,j,1)=curx(nx,j,2)
         curx(nx,j,nz)=curx(nx,j,nz1)
         cury(nx,j,1)=cury(nx,j,2)
         cury(nx,j,nz)=cury(nx,j,nz1)
         curz(nx,j,1)=curz(nx,j,2)
         curz(nx,j,nz)=curz(nx,j,nz1)
      enddo
!
!     parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do k=2,nz1
         curx(1,1,k)=curx(1,2,k)
         curx(nx,1,k)=curx(nx,2,k)
         cury(1,1,k)=cury(1,2,k)
         cury(nx,1,k)=cury(nx,2,k)
         curz(1,1,k)=curz(1,2,k)
         curz(nx,1,k)=curz(nx,2,k)
!
         curx(1,ny,k)=curx(1,ny1,k)
         curx(nx,ny,k)=curx(nx,ny1,k)
         cury(1,ny,k)=cury(1,ny1,k)
         cury(nx,ny,k)=cury(nx,ny1,k)
         curz(1,ny,k)=curz(1,ny1,k)
         curz(nx,ny,k)=curz(nx,ny1,k)
!
      enddo
!
      curx(1,1,1)=curx(1,1,2)
      curx(1,ny,1)=curx(1,ny,2)
      curx(1,1,nz)=(curx(1,1,nz1)+curx(1,2,nz))/2.
      curx(1,ny,nz)=(curx(1,ny,nz1)+curx(1,ny1,nz))/2.
      curx(nx,1,1)=(curx(nx,1,2)+curx(nx,2,1)+curx(nx1,1,1))/3.
      curx(nx,ny,1)=(curx(nx,ny,2)+curx(nx,ny1,1)+curx(nx1,ny,1))/3.
      curx(nx,1,nz)=(curx(nx,1,nz1)+curx(nx,2,nz)+curx(nx1,1,nz))/3.
      curx(nx,ny,nz)=(curx(nx,ny,nz1)+curx(nx,ny1,nz) &
           + curx(nx1,ny,nz))/3.
!
      cury(1,1,1)=cury(1,1,2)
      cury(1,ny,1)=cury(1,ny,2)
      cury(1,1,nz)=(cury(1,1,nz1)+cury(1,2,nz))/2.
      cury(1,ny,nz)=(cury(1,ny,nz1)+cury(1,ny1,nz))/2.
      cury(nx,1,1)=(cury(nx,1,2)+cury(nx,2,1)+cury(nx1,1,1))/3.
      cury(nx,ny,1)=(cury(nx,ny,2)+cury(nx,ny1,1)+cury(nx1,ny,1))/3.
      cury(nx,1,nz)=(cury(nx,1,nz1)+cury(nx,2,nz)+cury(nx1,1,nz))/3.
      cury(nx,ny,nz)=(cury(nx,ny,nz1)+cury(nx,ny1,nz) &
           + cury(nx1,ny,nz))/3.
!
      curz(1,1,1)=curz(1,1,2)
      curz(1,ny,1)=curz(1,ny,2)
      curz(1,1,nz)=(curz(1,1,nz1)+curz(1,2,nz))/2.
      curz(1,ny,nz)=(curz(1,ny,nz1)+curz(1,ny1,nz))/2.
      curz(nx,1,1)=(curz(nx,1,2)+curz(nx,2,1)+curz(nx1,1,1))/3.
      curz(nx,ny,1)=(curz(nx,ny,2)+curz(nx,ny1,1)+curz(nx1,ny,1))/3.
      curz(nx,1,nz)=(curz(nx,1,nz1)+curz(nx,2,nz)+curz(nx1,1,nz))/3.
      curz(nx,ny,nz)=(curz(nx,ny,nz1)+curz(nx,ny1,nz) &
           + curz(nx1,ny,nz))/3.
!
!
      return
      end
!
!     *************************************************
!
      subroutine bande(efldx,efldy,efldz,bsx,bsy,bsz, &
             curx,cury,curz,evx,evy,evz,btot, &
             epres,qrho,hrho,orho,rst,resist,reynolds, &
             nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
             ijmid,nummid,ijzero,numzero,mbndry,mmid,mzero, &
             rx,ry,rz)
!
!      calculates the surface magnetic field bs
!      and the body electric field eb
!
      dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
            efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
            curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
            evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz), &
            qrho(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd), &
            epres(nx,ny,nz,ngrd),rst(mbndry,nx,ny,nz),btot(nx,ny,nz)
      integer ijmid(mbndry,3,mmid),ijzero(mbndry,3,mzero)
!
!      ohm's law: ve = vi-j/ne
!
!      print*, 'in bande subroutine'

! parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do k=1,nz
      do j=1,ny
      do i=1,nx
!
      avx=evx(i,j,k)
      avy=evy(i,j,k)
      avz=evz(i,j,k)
!
      abx=bsx(i,j,k)
      aby=bsy(i,j,k)
      abz=bsz(i,j,k)
!
      efldx(i,j,k)=-(avy*abz-avz*aby)
      efldy(i,j,k)=-(avz*abx-avx*abz)
      efldz(i,j,k)=-(avx*aby-avy*abx)
!
      enddo
      enddo
      enddo
!
!     add in grad p term
!
      dxt=2.*rx
      dyt=2.*ry
      dzt=2.*rz
! parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do k=2,nz-1
      kp=k+1
      km=k-1
      do j=2,ny-1
      jp=j+1
      jm=j-1
      do i=2,nx-1
      ip=i+1
      im=i-1
!
      arho=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
             orho(i,j,k,m)/rmasso)*reynolds
!
      efldx(i,j,k)=efldx(i,j,k)- &
              ((epres(ip,j,k,m)-epres(im,j,k,m))/dxt)/arho
      efldy(i,j,k)=efldy(i,j,k)- &
              ((epres(i,jp,k,m)-epres(i,jm,k,m))/dyt)/arho
      efldz(i,j,k)=efldz(i,j,k)- &
              ((epres(i,j,k,mp)-epres(i,j,k,mm))/dzt)/arho
!
!     with grad b corrections if proven - comment above if implemented
!
!     abtot=btot(i,j,k)+0.125*(btot(ip,j,k)+btot(im,j,k)+
!    +     btot(i,jp,k)+btot(i,jm,k)+btot(i,j,kp)+btot(i,j,km))/6.
!     epress=1.5*epres(i,j,k,m)
!     efldx(i,j,k)=efldx(i,j,k)-(1./dxt/arho)*
!    +        ( (epres(m,ip,j,k)-epres(m,im,j,k))
!    +      +epress*(btot(ip,j,k)-btot(im,j,k))/abtot)
!     efldy(i,j,k)=efldy(i,j,k)-(1./dyt/arho)*
!    +        ( (epres(m,i,jp,k)-epres(m,i,jm,k))
!    +      +epress*(btot(ip,j,k)-btot(im,j,k))/abtot)
!     efldz(i,j,k)=efldz(i,j,k)-(1./dzt/arho)*
!    +        ( (epres(i,j,k,mp)-epres(i,j,k,mm))
!    +      +epress*(btot(ip,j,k)-btot(im,j,k))/abtot)
      enddo
      enddo
      enddo
!

!      print*,'boundary conditions'
!      add in ionospheric resistance
!
      if((m.le.mbndry).and.(resist.lt.5000.))then
! parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
       do n=1,nummid
         i=ijmid(m,1,n)
         j=ijmid(m,2,n)
         k=ijmid(m,3,n)
!
         eden=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
             orho(i,j,k,m)/rmasso)
!
         efldx(i,j,k)=efldx(i,j,k)+rst(i,j,k,m)*curx(i,j,k)/eden
         efldy(i,j,k)=efldy(i,j,k)+rst(i,j,k,m)*cury(i,j,k)/eden
         efldz(i,j,k)=efldz(i,j,k)+rst(i,j,k,m)*curz(i,j,k)/eden
       enddo

!       print*, numzero
       do n=1,numzero
         i=ijzero(m,1,n)
         j=ijzero(m,2,n)
         k=ijzero(m,3,n)
!
         eden=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
             orho(i,j,k,m)/rmasso)

         efldx(i,j,k)=efldx(i,j,k)+rst(i,j,k,m)*curx(i,j,k)/eden
         efldy(i,j,k)=efldy(i,j,k)+rst(i,j,k,m)*cury(i,j,k)/eden
         efldz(i,j,k)=efldz(i,j,k)+rst(i,j,k,m)*curz(i,j,k)/eden

       enddo
      endif
!
!      boundary condtions for terrestrial surface currents
!
!     if(m.eq.1)then
!      do n=1,nummid
!        i=ijmid(1,n)
!        j=ijmid(2,n)
!        k=ijmid(3,n)
!        efldx(i,j,k)=0.
!        efldy(i,j,k)=0.
!        efldz(i,j,k)=0.
!       enddo
!      endif
!
!      set flank bounday conditions
!
!      print*,'flanks conditions'
      nx1=nx-1
      ny1=ny-1
      nz1=nz-1
!
! parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do k=1,nz
      do j=1,ny
       efldx(1,j,k)=efldx(2,j,k)
       efldx(nx,j,k)=efldx(nx1,j,k)
       efldy(1,j,k)=efldy(2,j,k)
       efldy(nx,j,k)=efldy(nx1,j,k)
       efldz(1,j,k)=efldz(2,j,k)
       efldz(nx,j,k)=efldz(nx1,j,k)
      enddo
      enddo
!
! parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do j=1,ny
      do i=1,nx
       efldx(i,j,1)=efldx(i,j,2)
       efldx(i,j,nz)=efldx(i,j,nz1)
       efldy(i,j,1)=efldy(i,j,2)
       efldy(i,j,nz)=efldy(i,j,nz1)
       efldz(i,j,1)=efldz(i,j,2)
       efldz(i,j,nz)=efldz(i,j,nz1)
      enddo
      enddo
!
! parallelizes loop. RW, oct. 23, 2002
!$omp  parallel do
      do k=1,nz
      do i=1,ny
       efldx(i,1,k)=efldx(i,2,k)
       efldx(i,ny,k)=efldx(i,ny1,k)
       efldy(i,1,k)=efldy(i,2,k)
       efldy(i,ny,k)=efldy(i,ny1,k)
       efldz(i,1,k)=efldz(i,2,k)
       efldz(i,ny,k)=efldz(i,ny1,k)
      enddo

      enddo
!
      return
      end
!
!     ***********************************************
!
!     ************************************************
!
      subroutine fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz, &
            evx,evy,evz,tvx,tvy,tvz, &
            nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,reynolds)
!
!     converts momentum into velocity for graphics
!
      dimension qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd), &
                qpz(nx,ny,nz,ngrd),qrho(nx,ny,nz,ngrd), &
                hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd), &
                hpz(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd), &
                opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd), &
                opz(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd), &
                curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
                tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
                evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
!
!
!$omp  parallel do
      do k=1,nz
      do j=1,ny
      do i=1,nx
!      eden=amax1(qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m),
!    +                 0.0001)
       qden=(qrho(i,j,k,m)+0.000001)/rmassq
       hden=(hrho(i,j,k,m)+0.000001)/rmassh
       oden=(orho(i,j,k,m)+0.000001)/rmasso
       tden=qden+hden+oden
!      keep sepearate the ion and current components
!
       tvx(i,j,k)=(qpx(i,j,k,m)/rmassq+hpx(i,j,k,m)/rmassh &
                     +opx(i,j,k,m)/rmasso)/tden
       evx(i,j,k)= tvx(i,j,k) - curx(i,j,k)/tden/reynolds
!
       tvy(i,j,k)=(qpy(i,j,k,m)/rmassq+hpy(i,j,k,m)/rmassh &
                     +opy(i,j,k,m)/rmasso)/tden
       evy(i,j,k)= tvy(i,j,k) - cury(i,j,k)/tden/reynolds
!
       tvz(i,j,k)=(qpz(i,j,k,m)/rmassq+hpz(i,j,k,m)/rmassh &
                     +opz(i,j,k,m)/rmasso)/tden
       evz(i,j,k)= tvz(i,j,k)  -curz(i,j,k)/tden/reynolds
      enddo
      enddo
      enddo
!
      return
      end

!
!     ************************************************
      subroutine fnd_vel(px,py,pz,rho,vx,vy,vz,nx,ny,nz,ngrd,m)
!
!     converts momentum into velocity for graphics
!
      dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd), &
                pz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd), &
                vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
!
!$omp  parallel do
      do k=1,nz
      do j=1,ny
      do i=1,nx
       arho=amax1(rho(i,j,k,m),0.0001)
       vx(i,j,k)=px(i,j,k,m)/arho
       vy(i,j,k)=py(i,j,k,m)/arho
       vz(i,j,k)=pz(i,j,k,m)/arho
      enddo
      enddo
      enddo
!
      return
      end



!     ************************************************
!
      subroutine fnd_pres(px,py,pz,pxy,pxz,pyz,p_para,p_perp,p_cross, &
                          vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,ngrd,m)

      dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd), &
        pxy(nx,ny,nz,ngrd),pxz(nx,ny,nz,ngrd),pyz(nx,ny,nz,ngrd), &
        p_para(nx,ny,nz,ngrd),p_perp(nx,ny,nz,ngrd),p_cross(nx,ny,nz,ngrd), &
        vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz),bxt(nx,ny,nz), &
        byt(nx,ny,nz),bzt(nx,ny,nz)

      real abx,aby,abz,bmag,avx,avy,avz,vmag
      real vcbx,vcby,vcbz,vcmag
      real p_para,p_perp,p_cross
      real presx,presy,presz,presmag
      real presxy,presxz,presyz

!     find anisotropy values for species 
!
         do k=1,nz
            do j=1,ny
               do i=1,nx
!
                  presx=px(i,j,k,m)
                  presy=py(i,j,k,m)
                  presz=pz(i,j,k,m)
                  presxy=pxy(i,j,k,m)
                  presxz=pxz(i,j,k,m)
                  presyz=pyz(i,j,k,m)
                  presmag=sqrt(presx**2+presy**2+presz**2 + &
                    2*(presxy**2) +2*(presxz**2) + 2*(presyz**2))+1.e-11
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
                  vcbx=(avy*abz-avz*aby)
                  vcby=-(avx*abz-avz*abx)
                  vcbz=(avx*aby-avy*abx)
                  vcmag=sqrt(vcbx**2+vcby**2+vcbz**2) &
                       +1.e-12
!
!     find vparallel
!
                  p_para(i,j,k,m)=sqrt((abs(abx)*presx+abs(aby)*presxy+abs(abz)*presxz)**2+ &
                              (abs(abx)*presxy+abs(aby)*presy+abs(abz)*presyz)**2+ &
                              (abs(abx)*presxz+abs(aby)*presyz+abs(abz)*presz)**2)/bmag
                  p_cross(i,j,k,m)=sqrt((abs(vcbx)*presx+abs(vcby)*presxy+abs(vcbz)*presxz)**2+ &
                               (abs(vcbx)*presxy+abs(vcby)*presy+abs(vcbz)*presyz)**2+ &
                               (abs(vcbx)*presxz+abs(vcby)*presyz+abs(vcbz)*presz)**2)/vcmag
                  p_perp(i,j,k,m)=sqrt(abs(presmag**2-(p_para(i,j,k,m)**2+p_cross(i,j,k,m)**2)))

               enddo
            enddo
         enddo

      end
!     ************************************************
!
      subroutine totfld(bx,bx0,btx,nx,ny,nz,ngrd,m)
!
!
!     calculates the total magnetic field from the perturbed and
!        stationary magnetic field
!
      dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd),btx(nx,ny,nz)
!
!$omp  parallel do
      do k=1,nz
      do j=1,ny
      do i=1,nx
        btx(i,j,k)=bx0(i,j,k,m)+bx(i,j,k,m)
      enddo
      enddo
      enddo
!
      return
      end
!
!     *************************************************
!
      subroutine tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
!
!      initialize static magnetic field along entire grid
!

!
      dimension bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz), &
            btot(nx,ny,nz)
!
      do k=1,nz
      do j=1,ny
      do i=1,nx
        atot=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2+bsz(i,j,k)**2)
        btot(i,j,k)=amax1(atot,1.e-5)
      enddo
      enddo
      enddo
!
      return
      end

