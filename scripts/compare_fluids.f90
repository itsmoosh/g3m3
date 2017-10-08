
!   This program slices binary fluid files from multifluid
!	simulations into several ASCII text files.
!
!	Required inputs:
!	cutter_input	input file from the run generating the fluid file
!	cutter_fluid	binary fluid file to be sliced
!
program cut_fluid_to_ascii

      implicit none
      integer, parameter :: nx=121,ny=121,nz=61,ngrd=5, &
                  mbndry=1,msrf=2000,mmid=1500,mzero=5000, &
                  ncraft=30,ncts=281, division=2

      integer, parameter :: dp = kind(1.d0)  !	Double precision

		integer :: percent = 0

      real t, t2
		real t_equiv, ut, ut2, utstart

	  real,parameter :: wind_adjust=4./3., limit=60.

	character*12 :: wd_pre = 'diff_'
	character*30 :: timebox
	integer, parameter :: m_max = 2
	  
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
           xspac(ngrd), grd_diff
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

      real planet_rad, moon_rad, pl_ratio, &
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real bx2(nx,ny,nz,ngrd),by2(nx,ny,nz,ngrd),bz2(nx,ny,nz,ngrd), &
           qpx2(nx,ny,nz,ngrd),qpy2(nx,ny,nz,ngrd),qpz2(nx,ny,nz,ngrd), &
           qrho2(nx,ny,nz,ngrd),opx2(nx,ny,nz,ngrd), &
           opy2(nx,ny,nz,ngrd),opz2(nx,ny,nz,ngrd), &
           orho2(nx,ny,nz,ngrd),hpx2(nx,ny,nz,ngrd), &
           hpy2(nx,ny,nz,ngrd),hpz2(nx,ny,nz,ngrd), &
           hrho2(nx,ny,nz,ngrd),epres2(nx,ny,nz,ngrd)
!
      real qpresx2(nx,ny,nz,ngrd),qpresy2(nx,ny,nz,ngrd), &
           qpresz2(nx,ny,nz,ngrd),qpresxy2(nx,ny,nz,ngrd), &
           qpresxz2(nx,ny,nz,ngrd),qpresyz2(nx,ny,nz,ngrd), &
           qpara2(nx,ny,nz,ngrd),qperp2(nx,ny,nz,ngrd), &
           qcross2(nx,ny,nz,ngrd)

      real hpresx2(nx,ny,nz,ngrd),hpresy2(nx,ny,nz,ngrd), &
           hpresz2(nx,ny,nz,ngrd),hpresxy2(nx,ny,nz,ngrd), &
           hpresxz2(nx,ny,nz,ngrd),hpresyz2(nx,ny,nz,ngrd), &
           hpara2(nx,ny,nz,ngrd),hperp2(nx,ny,nz,ngrd), &
           hcross2(nx,ny,nz,ngrd)

      real opresx2(nx,ny,nz,ngrd),opresy2(nx,ny,nz,ngrd), &
           opresz2(nx,ny,nz,ngrd),opresxy2(nx,ny,nz,ngrd), &
           opresxz2(nx,ny,nz,ngrd),opresyz2(nx,ny,nz,ngrd), &
           opara2(nx,ny,nz,ngrd),operp2(nx,ny,nz,ngrd), &
           ocross2(nx,ny,nz,ngrd)

      real qvx2(nx,ny,nz,ngrd), qvy2(nx,ny,nz,ngrd), &
           qvz2(nx,ny,nz,ngrd), ovx2(nx,ny,nz,ngrd), &
           ovy2(nx,ny,nz,ngrd), ovz2(nx,ny,nz,ngrd), &
           hvx2(nx,ny,nz,ngrd), hvy2(nx,ny,nz,ngrd), &
           hvz2(nx,ny,nz,ngrd)


      real evelx2(nx,ny,nz,ngrd), evely2(nx,ny,nz,ngrd), &
           evelz2(nx,ny,nz,ngrd)

      real evx2(nx,ny,nz), evy2(nx,ny,nz), evz2(nx,ny,nz)
      real vvx2(nx,ny,nz), vvy2(nx,ny,nz), vvz2(nx,ny,nz)

      real efldx2(nx,ny,nz,ngrd), efldy2(nx,ny,nz,ngrd), &
           efldz2(nx,ny,nz,ngrd)

      real ex2(nx,ny,nz), ey2(nx,ny,nz), ez2(nx,ny,nz)
      real tvx2(nx,ny,nz), tvy2(nx,ny,nz), tvz2(nx,ny,nz)

!     unperturbed quantities
!
      real bx02(nx,ny,nz,ngrd),by02(nx,ny,nz,ngrd),bz02(nx,ny,nz,ngrd)
      real bxs2(nx,ny,nz,ngrd),bys2(nx,ny,nz,ngrd),bzs2(nx,ny,nz,ngrd)
      real br02(nx,ny,nz,ngrd),bt02(nx,ny,nz,ngrd),bp02(nx,ny,nz,ngrd)
      real bxt2(nx,ny,nz),byt2(nx,ny,nz),bzt2(nx,ny,nz),btot2(nx,ny,nz)


      real curx2(nx,ny,nz),cury2(nx,ny,nz),curz2(nx,ny,nz), &
           resistive2(nx,ny,nz),curx_all2(nx,ny,nz,ngrd), &
           cury_all2(nx,ny,nz,ngrd),curz_all2(nx,ny,nz,ngrd)



!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
      real nq2(nx,ny,nz,ngrd), no2(nx,ny,nz,ngrd), &
           nh2(nx,ny,nz,ngrd), ne2(nx,ny,nz,ngrd), &
           vte2(nx,ny,nz,ngrd), vbe2(nx,ny,nz,ngrd), &
           vtq2(nx,ny,nz,ngrd), vto2(nx,ny,nz,ngrd), &
           vth2(nx,ny,nz,ngrd), vbq2(nx,ny,nz,ngrd), &
           vbo2(nx,ny,nz,ngrd), vbh2(nx,ny,nz,ngrd), &
           vbq_para2(nx,ny,nz,ngrd),vbq_perp2(nx,ny,nz,ngrd), &
           vbh_para2(nx,ny,nz,ngrd),vbh_perp2(nx,ny,nz,ngrd), &
           vbo_para2(nx,ny,nz,ngrd),vbo_perp2(nx,ny,nz,ngrd)

      real bxtmp2, bytmp2, bztmp2, bmag2, dotq2, cosq2, &
           doth2, cosh2, doto2, coso2, para2, perp2, &
           vq_flat2, vo_flat2, vh_flat2, ve_flat2, &
           vbq_para_intp2, vbq_perp_intp2, vbh_para_intp2, &
           vbh_perp_intp2, vbo_para_intp2, vbo_perp_intp2, &
           vtq_intp2, vth_intp2, vto_intp2, &
           numq_intp2, numh_intp2, numo_intp2, &
           vbe_intp2, vte_intp2, nume_intp2, &
           thetaq2,thetah2,thetao2,drftnorm2, &
           qeratio2, oeratio2, heratio2

      real theta2,q_costheta2,o_costheta2,h_costheta2,sintheta2
      real s_xmin2(ngrd),s_xmax2(ngrd),s_ymin2(ngrd), &
           s_ymax2(ngrd),s_zmin2(ngrd),s_zmax2(ngrd), &
           r_xmin2(ngrd), r_xmax2(ngrd), r_ymin2(ngrd), &
           r_ymax2(ngrd), r_zmin2(ngrd), r_zmax2(ngrd)

      real nrgy2,mq2,mo2,mh2,lognrgy2,kreal2,vbinsreal2

      real ml2(nlines),mlt2(nlines),rad_sat2(nlines), &
           xsat2(nlines),ysat2(nlines),zsat2(nlines), &
           xplot_f2, yplot_f2

      real,allocatable :: gx2(:),gy2(:),gz2(:),gr2(:), &
           gbx2(:),gby2(:),gbz2(:),xpos2(:),ypos2(:),zpos2(:), &
           inbox2(:,:)

      real xtemp2, ytemp2, ztemp2

      real phi2, vqpar2, vqprp2, vopar2, voprp2, vhpar2, vhprp2

      integer, allocatable :: nbox2(:), grd_sys2(:)

!      integer a,b,i,j,k,m,n,posit, boxtemp, count1, count2, &
!           count3, count4, count5, count6, count7, count8, &
!           count9, moon
!      integer ix, iy, iz, i_xmin, i_xmax, i_ymin, i_ymax, &
!           i_zmin, i_zmax, nchf

      real aqpres2, ahpres2, aopres2, aepres2, qden2, hden2, oden2, &
           eden2, qtemp2, htemp2, otemp2, etemp2, qpress2, hpress2, &
           opress2, epress2, qdens2, hdens2, odens2, edens2, &
           qtemps2, htemps2, otemps2, etemps2, abx_2, aby_2, abz_2, &
           ri2, rj2, rk2, rad2, bsurmag2, time2, abx22, aby22, abz22, &
           ri_moon2, rj_moon2, rk_moon2,aoden2,aqden2,ahden2,aeden2, &
           aotemp2,ahtemp2,aqtemp2,aetemp2

!      integer nx1, nx2, ny1, ny2, nz1, nz2

!      character*1 cut,box

!      character*32 wd1,wd2,wd3,wd4,wd5,wd6,wd7,wd8,wd9

      real planet_rad2, moon_rad2, pl_ratio2, &
           ofrac2, frac_o2,dist2,tot_v2,rot_mach2,j_rad2,j_phi2, &
           cur_norm2, planet_per2

      real,allocatable :: yplot2(:,:),flux_p2(:,:), &
           flux_q2(:,:),flux_o2(:,:),flux_h2(:,:), &
           flux_e2(:,:), bxg2(:), byg2(:), bzg2(:), &
           brg2(:), btg2(:), bpg2(:),cxg2(:), cyg2(:), czg2(:)

!~*~*~*~*~*~**~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~

      integer ijsrf2(mbndry,3,msrf),ijmid2(mbndry,3,mmid), &
              ijzero2(mbndry,3,mzero)
      integer numsrf2(mbndry),nummid2(mbndry),numzero2(mbndry)
      real    parm_srf2(mbndry,7,msrf),parm_mid2(mbndry,7,mmid), &
              parm_zero2(mbndry,7,mzero)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical start,add_dip,ringo,update,save_dat,write_dat, &
              spacecraft,tilting,warp,reload, &
              divb_lores,divb_hires,divb_on

	character*10 bodyname,moonname

      real tmax,stepsz,tsave,xdip,ydip,zdip,r_inner, &
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
           chirho,chipxyz,chierg, &
           difrho,difpxyz,diferg

      real xmoon_re,ymoon_re,zmoon_re,cs_moon, &
           qden_moon,hden_moon,oden_moon,ti_te_moon, &
           alf_moon,xdip_moon,ydip_moon,zdip_moon,offset, &
           orbit_moon,theta_moon,ani,resist_moon,den_lunar, &
           reduct,t_torus,aniso_factor,ani_q,ani_h,ani_o,uday

      integer no_shock,ntgraph

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
        namelist/option/tmax,ntgraph,stepsz,start,tsave,isotropic
        namelist/planet/bodyname,moonname,xdip,ydip,zdip,r_inner, &
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

	timebox = trim(wd_pre)//'timebox.dat'
      open(2,file='cutter_input',status='old',form='formatted')
      open(3,file=trim(timebox),status='unknown',form='formatted')

!     read input parameters
!
      read(2,option)
      read(2,planet)
      read(2,speeds)
      read(2,windy)
      read(2,lunar)
      read(2,physical)
      read(2,smooth)

	close(2)
!
!~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
!     large grid system
!        ncore -> grid gets all the information from this grid no.
!        nbdry -> which grid data will be applied to this grid no.for bndry
!
      write(*,*) 'Grid xmin: ','Grid xmax: ','Grid ymin: ','Grid ymax: ', &
		'Grid zmin: ','Grid zmax: ','Grid spacing: '
	  do i=1,ngrd
		xspac(i) = division**(i-1)
		grd_xmin(i) = -limit*xspac(i)
		grd_xmax(i) = limit*xspac(i)
	    grd_ymin(i) = -limit*xspac(i)
		grd_ymax(i) = limit*xspac(i)        
		grd_zmin(i) = -limit/2*xspac(i)
		grd_zmax(i) = limit/2*xspac(i)
		!
		if(i.eq.ngrd) then
			grd_diff = grd_xmax(i) * wind_adjust - grd_xmax(i)
			grd_xmax(i) = grd_xmax(i) + grd_diff
			grd_xmin(i) = grd_xmin(i) + grd_diff
		endif
		!
		write(*,*) grd_xmin(i), grd_xmax(i), &
		grd_ymin(i), grd_ymax(i), &
		grd_zmin(i), grd_zmax(i), xspac(i)
		ix = 1 + ( grd_xmax(i) - grd_xmin(i) ) / xspac(i)
		iy = 1 + ( grd_ymax(i) - grd_ymin(i) ) / xspac(i)
	    iz = 1 + ( grd_zmax(i) - grd_zmin(i) ) / xspac(i)
	    if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz)) then
	        write(*,*)'ERROR: grid spacing does not match number of points.', &
				' Actual x,y,z: ',ix,iy,iz,'Expected x,y,z: ',nx,ny,nz
	        stop
	    endif
	  enddo
!
!
!!!!

!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~
! read in data file to cut 
!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~

      write(*,*) 'Begin reading data files:'

      ! input fluid file 'cutter_fluid1'
      nchf=1
      open(nchf,file='cutter_fluid1',status='unknown',form='unformatted')

!
       read(nchf)t
		write(*,*) 't1 = ', t
		write(3,*) 't1 = ', t
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
       close(nchf)

      ! input fluid file 'cutter_fluid1'
		nchf=10
      open(nchf,file='cutter_fluid2',status='unknown',form='unformatted')
!
       read(nchf)t2
		write(*,*) 't2 = ', t2
		write(3,*) 't2 = ', t2
       read(nchf)qrho2
       read(nchf)qpx2
       read(nchf)qpy2
       read(nchf)qpz2
       read(nchf)qpresx2
       read(nchf)qpresy2
       read(nchf)qpresz2
       read(nchf)qpresxy2
       read(nchf)qpresxz2
       read(nchf)qpresyz2

       read(nchf)hrho2
       read(nchf)hpx2
       read(nchf)hpy2
       read(nchf)hpz2
       read(nchf)hpresx2
       read(nchf)hpresy2
       read(nchf)hpresz2
       read(nchf)hpresxy2
       read(nchf)hpresxz2
       read(nchf)hpresyz2

       read(nchf)orho2
       read(nchf)opx2
       read(nchf)opy2
       read(nchf)opz2
       read(nchf)opresx2
       read(nchf)opresy2
       read(nchf)opresz2
       read(nchf)opresxy2
       read(nchf)opresxz2
       read(nchf)opresyz2

       read(nchf)bx2
       read(nchf)by2
       read(nchf)bz2
       read(nchf)epres2
       read(nchf)bx02
       read(nchf)by02
       read(nchf)bz02
       read(nchf)parm_srf2,parm_mid2,parm_zero2, &
                ijzero2,numzero2,ijmid2,nummid2,ijsrf2,numsrf2
       close(nchf)

!
      moon_rad = 2575.   !km
      planet_rad = 60268.   !km
      planet_per = 10.7  !hours
      moon_rad2 = 2575.   !km
      planet_rad2 = 60268.   !km
      planet_per2 = 10.7  !hours
      t_equiv=planet_rad*re_equiv/v_equiv

      cur_norm = 1.37e-9
      cur_norm2 = 1.37e-9
!      t_equiv=(radius/1000.0)*re_equiv/v_equiv


      ut=utstart+t*t_equiv/3600.
      ut2=utstart+t2*t_equiv/3600.
      write(*,*) 'ut1 = ',ut
      write(3,*) 'ut1 = ',ut
      write(*,*) 'ut2 = ',ut2
      write(3,*) 'ut2 = ',ut2
      close(3)


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
		 call calcur(bx2,by2,bz2,nx,ny,nz,ngrd,m, &
              curx2,cury2,curz2,rx,ry,rz)
         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  curx_all(i,j,k,m) = curx(i,j,k)
                  cury_all(i,j,k,m) = cury(i,j,k)
                  curz_all(i,j,k,m) = curz(i,j,k)

                  curx_all2(i,j,k,m) = curx2(i,j,k)
                  cury_all2(i,j,k,m) = cury2(i,j,k)
                  curz_all2(i,j,k,m) = curz2(i,j,k)
               end do
            end do
         end do


!     write(6,*)' totbfld ing now'
      call totfld(bx,bx0,bxt,nx,ny,nz,ngrd,m)
      call totfld(by,by0,byt,nx,ny,nz,ngrd,m)
      call totfld(bz,bz0,bzt,nx,ny,nz,ngrd,m)

      call totfld(bx2,bx02,bxt2,nx,ny,nz,ngrd,m)
      call totfld(by2,by02,byt2,nx,ny,nz,ngrd,m)
      call totfld(bz2,bz02,bzt2,nx,ny,nz,ngrd,m)
!
!     find magnitude of b
!
      call tot_b(btot,bxt,byt,bzt,nx,ny,nz)
      call tot_b(btot2,bxt2,byt2,bzt2,nx,ny,nz)

      call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz, &
            evx,evy,evz,tvx,tvy,tvz, &
            nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,reynolds)
      call fnd_evel(qpx2,qpy2,qpz2,qrho2,hpx2,hpy2,hpz2,hrho2, &
            opx2,opy2,opz2,orho2,curx2,cury2,curz2, &
            evx2,evy2,evz2,tvx2,tvy2,tvz2, &
            nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso,reynolds)


!     write(6,*)' bande ing now' 
         call bande(ex,ey,ez,bxt,byt,bzt, &
              curx,cury,curz,evx,evy,evz,btot, &
              epres,qrho,hrho,orho,resistive,resist,reynolds, &
              nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
              ijmid,nummid,ijzero,numzero,mbndry,mmid,mzero, &
              rx,ry,rz)
		 call bande(ex2,ey2,ez2,bxt2,byt2,bzt2, &
              curx2,cury2,curz2,evx2,evy2,evz2,btot2, &
              epres2,qrho2,hrho2,orho2,resistive2,resist,reynolds, &
              nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
              ijmid2,nummid2,ijzero2,numzero2,mbndry,mmid,mzero, &
              rx,ry,rz)

         do k = 1, nz
            do j = 1, ny
               do i = 1, nx
                  efldx(i,j,k,m) = ex(i,j,k)/9.26
                  efldy(i,j,k,m) = ey(i,j,k)/9.26 ! drift speed normalization
                  efldz(i,j,k,m) = ez(i,j,k)/9.26

                  efldx2(i,j,k,m) = ex2(i,j,k)/9.26
                  efldy2(i,j,k,m) = ey2(i,j,k)/9.26 ! drift speed normalization
                  efldz2(i,j,k,m) = ez2(i,j,k)/9.26
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
      call fnd_vel(qpx2,qpy2,qpz2,qrho2,vvx2,vvy2,vvz2,nx,ny,nz,ngrd,m)

        do k=1,nz
           do j=1,ny
              do i=1,nx
                 qvx(i,j,k,m) = vvx(i,j,k)
                 qvy(i,j,k,m) = vvy(i,j,k)
                 qvz(i,j,k,m) = vvz(i,j,k)

                 qvx2(i,j,k,m) = vvx2(i,j,k)
                 qvy2(i,j,k,m) = vvy2(i,j,k)
                 qvz2(i,j,k,m) = vvz2(i,j,k)
              enddo
           enddo
        enddo
		!
		!
		!      write(*,*) "finding qpres-s now!"
		call fnd_pres(qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz,&
		    qpara,qcross,qperp, &
		    vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,ngrd,m)
		call fnd_pres(qpresx2,qpresy2,qpresz2,qpresxy2,qpresxz2,qpresyz2,&
		    qpara2,qcross2,qperp2, &
		    vvx2,vvy2,vvz2,bxt2,byt2,bzt2,nx,ny,nz,ngrd,m)
	enddo !loop over m
	!
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
		write(*,*) 'Making ASCII data files.'
		!     output data file names
		!     ------------
		!
		do m=1,m_max
			percent = (m-1)*20
			call loadbar(percent)
			write(box,'(i1)') m
			wd1 = trim(wd_pre)//"plasma-box"//trim(adjustl(box)//".dat")
			wd2 = trim(wd_pre)//"flow-box"//trim(adjustl(box)//".dat")
			wd3 = trim(wd_pre)//"field-box"//trim(adjustl(box)//".dat")
			wd4 = trim(wd_pre)//"rad-box"//trim(adjustl(box)//".dat")
			wd5 = trim(wd_pre)//"model-box"//trim(adjustl(box)//".dat")
			wd6 = trim(wd_pre)//"efield-box"//trim(adjustl(box)//".dat")
			wd7 = trim(wd_pre)//"pressure-box"//trim(adjustl(box)//".dat")
			wd8 = trim(wd_pre)//"rotational-box"//trim(adjustl(box)//".dat")
			wd9 = trim(wd_pre)//"currents-box"//trim(adjustl(box)//".dat")
			open(7,file=wd1,status='replace',form='formatted')
			open(8,file=wd2,status='replace',form='formatted')
			open(9,file=wd3,status='replace',form='formatted')
			open(10,file=wd4,status='replace',form='formatted')
			open(11,file=wd5,status='replace',form='formatted')
			open(12,file=wd6,status='replace',form='formatted')
			open(13,file=wd7,status='replace',form='formatted')
			open(14,file=wd8,status='replace',form='formatted')
			open(15,file=wd9,status='replace',form='formatted')
			!	Write header lines to cut data files
			write(14,'(3(A9))') 'ri','rj','rot_mach'
			write(15,'(2(A9),2(A14))') 'ri','rj','j_rad','j_phi'
			write(7,'(3(A9),8(A14))') 'ri','rj','rk',	'qdens','qtemp','hdens','htemp', 'odens','otemp','edens','etemp'
			write(8,'(3(A9),9(A14))') 'ri','rj','rk', 'qvx','qvy','qvz', 'ovx','ovy','ovz', 'hvx','hvy','hvz'
			write(9,'(3(A9),6(A14))') 'ri','rj','rk', 'abx','aby','abz','curx_all', 'cury_all','curz_all'
			write(10,'(5(A14))') 'rad', 'bsurmag', 'ri', 'rj', 'rk'
			write(11,'(3(A9),6(A14))') 'ri','rj','rk', 'abx','aby','abz','abx2','aby2','abz2'
			write(12,'(3(A9),3(A14))') 'ri','rj','rk',	'efldx','efldy','efldz'
			write(13,'(3(A9),10(A14))') 'ri','rj','rk', 'qpara','qperp','qcross'
			!
			! output loop, put data in physical units, 
			! and essentially write out what you are interested
			! in plotting
			!
			do k=1,nz
				percent = 100.*( (nz*(m-1)*1.0 + k*1.0)/(nz*m_max*1.0) )
				call loadbar(percent)

		    	do j=1,ny
		    		do i=1,nx
						!
		    			! normalized  pressures
						aqpres=amax1(0.0000001,0.333*(qpresx(i,j,k,m)+ &
							qpresy(i,j,k,m)+qpresz(i,j,k,m)))
						ahpres=amax1(0.0000001,0.333*(hpresx(i,j,k,m)+ &
							hpresy(i,j,k,m)+hpresz(i,j,k,m)))
						aopres=amax1(0.0000001,0.333*(opresx(i,j,k,m)+ &
							opresy(i,j,k,m)+opresz(i,j,k,m)))
						aepres=amax1(0.0000001,epres(i,j,k,m))
						!
						! normalized densities
						qden=amax1(qrho(i,j,k,m),0.000001)
						hden=amax1(0.0001,hrho(i,j,k,m))
						oden=amax1(0.0001,orho(i,j,k,m))
						qden=qden/rmassq
						hden=hden/rmassh
						oden=oden/rmasso
						eden=oden+hden+qden+0.0001
						!
						! get temps from pressures/densities
						qtemp=amax1(0.00001,aqpres/qden)
						htemp=ahpres/hden
						otemp=aopres/oden
						etemp=aepres/eden
						!
						qdens=alog10(qden)
						hdens=alog10(hden)
						odens=alog10(oden)
						edens=alog10(eden)
						!
						! total b field
						abx=b_equiv*(bx(i,j,k,m) + bx0(i,j,k,m))
						aby=b_equiv*(by(i,j,k,m) + by0(i,j,k,m))
						abz=b_equiv*(bz(i,j,k,m) + bz0(i,j,k,m))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						! normalized  pressures
						aqpres2=amax1(0.0000001,0.333*(qpresx2(i,j,k,m)+ &
							qpresy2(i,j,k,m)+qpresz2(i,j,k,m)))
						ahpres2=amax1(0.0000001,0.333*(hpresx2(i,j,k,m)+ &
							hpresy2(i,j,k,m)+hpresz2(i,j,k,m)))
						aopres2=amax1(0.0000001,0.333*(opresx2(i,j,k,m)+ &
							opresy2(i,j,k,m)+opresz2(i,j,k,m)))
						aepres2=amax1(0.0000001,epres2(i,j,k,m))
						!
						! normalized densities
						qden2=amax1(qrho2(i,j,k,m),0.000001)
						hden2=amax1(0.0001,hrho2(i,j,k,m))
						oden2=amax1(0.0001,orho2(i,j,k,m))
						qden2=qden2/rmassq
						hden2=hden2/rmassh
						oden2=oden2/rmasso
						eden2=oden2+hden2+qden2+0.0001
						!
						! get temps from pressures/densities
						qtemp2=amax1(0.00001,aqpres2/qden2)
						htemp2=ahpres2/hden2
						otemp2=aopres2/oden2
						etemp2=aepres2/eden2
						!
						qdens2=alog10(qden2)
						hdens2=alog10(hden2)
						odens2=alog10(oden2)
						edens2=alog10(eden2)
						!
						! total b field
						abx_2=b_equiv*(bx2(i,j,k,m) + bx02(i,j,k,m))
						aby_2=b_equiv*(by2(i,j,k,m) + by02(i,j,k,m))
						abz_2=b_equiv*(bz2(i,j,k,m) + bz02(i,j,k,m))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						!
						! radial distance component
						ri = grd_xmin(m) + (xspac(m)*real(i-1))
						ri = ri*re_equiv
						rj = grd_ymin(m) + (xspac(m)*real(j-1))
						rj = rj*re_equiv
						rk = grd_zmin(m) + (xspac(m)*real(k-1))
						rk = rk*re_equiv
						!
						rad = sqrt((ri*ri) + (rj*rj) + (rk*rk))
						!
						! k == 31 implied equatorial planet in my sims (whose?)
						if(k == 31) then
							j_rad = cur_norm*(curx(i,j,k)*ri + &
								cury(i,j,k)*rj)/rad
							j_rad2 = cur_norm2*(curx2(i,j,k)*ri + &
								cury2(i,j,k)*rj)/rad
							!
							j_phi = cur_norm*sqrt(curx(i,j,k)**2 + &
								cury(i,j,k)**2) - j_rad
							j_phi2 = cur_norm2*sqrt(curx2(i,j,k)**2 + &
								cury2(i,j,k)**2) - j_rad2
							!
							dist = sqrt(ri**2+rj**2)
							tot_v = (sqrt(qvx(i,j,k,m)**2+qvy(i,j,k,m)**2)-&
								(qvx(i,j,k,m)*ri+qvy(i,j,k,m)*rj)/dist) + &
								(sqrt(hvx(i,j,k,m)**2+hvy(i,j,k,m)**2)-&
								(hvx(i,j,k,m)*ri+hvy(i,j,k,m)*rj)/dist) + &
								(sqrt(ovx(i,j,k,m)**2+ovy(i,j,k,m)**2)-&
								(ovx(i,j,k,m)*ri+ovy(i,j,k,m)*rj)/dist)
							tot_v2 = (sqrt(qvx2(i,j,k,m)**2+qvy2(i,j,k,m)**2)-&
								(qvx2(i,j,k,m)*ri+qvy2(i,j,k,m)*rj)/dist) + &
								(sqrt(hvx2(i,j,k,m)**2+hvy2(i,j,k,m)**2)-&
								(hvx2(i,j,k,m)*ri+hvy2(i,j,k,m)*rj)/dist) + &
								(sqrt(ovx2(i,j,k,m)**2+ovy2(i,j,k,m)**2)-&
								(ovx2(i,j,k,m)*ri+ovy2(i,j,k,m)*rj)/dist)
							!
							rot_mach = tot_v / &
								( 2*pi*rad*planet_rad/(planet_per*3600) )
							rot_mach2 = tot_v2 / &
								( 2*pi*rad*planet_rad/(planet_per*3600) )
							!
							if( rad == 0.0 ) then
								rot_mach=0.0
								rot_mach2=0.0
							endif
							!
							rot_mach = (rot_mach - rot_mach2)/rot_mach
							write(14,'(3(f9.2),10(es14.6))')ri,rj,rot_mach
							j_rad = (j_rad - j_rad2)/j_rad
							j_phi = (j_phi - j_phi2)/j_phi
							write(15,'(2(f9.2),2(es14.6))')ri,rj,j_rad,j_phi
							!
						endif
						!
						bsurmag = sqrt((bx0(i,j,k,m)**2) &
							+ (by0(i,j,k,m)**2) &
							+ (bz0(i,j,k,m)**2))
							bsurmag = bsurmag*b_equiv
						bsurmag2 = sqrt((bx02(i,j,k,m)**2) &
							+ (by02(i,j,k,m)**2) &
							+ (bz02(i,j,k,m)**2))
							bsurmag2 = bsurmag2*b_equiv
						!
						qdens = (qdens - qdens2)/qdens
						qtemp = (qtemp - qtemp2)/qtemp
						hdens = (hdens - hdens2)/hdens
						htemp = (htemp - htemp2)/htemp
						odens = (odens - odens2)/odens
						otemp = (otemp - otemp2)/otemp
						edens = (edens - edens2)/edens
						etemp = (etemp - etemp2)/etemp
						write(7,'(3(f9.2),8(es14.6))')ri,rj,rk,&
							qdens,qtemp,hdens,htemp, &
							odens,otemp,edens,etemp
						!
						qvx(i,j,k,m) = (qvx(i,j,k,m) - qvx2(i,j,k,m))/qvx(i,j,k,m)
						qvy(i,j,k,m) = (qvy(i,j,k,m) - qvy2(i,j,k,m))/qvy(i,j,k,m)
						qvz(i,j,k,m) = (qvz(i,j,k,m) - qvz2(i,j,k,m))/qvz(i,j,k,m)
						ovx(i,j,k,m) = (ovx(i,j,k,m) - ovx2(i,j,k,m))/ovx(i,j,k,m)
						ovy(i,j,k,m) = (ovy(i,j,k,m) - ovy2(i,j,k,m))/ovy(i,j,k,m)
						ovz(i,j,k,m) = (ovz(i,j,k,m) - ovz2(i,j,k,m))/ovz(i,j,k,m)
						hvx(i,j,k,m) = (hvx(i,j,k,m) - hvx2(i,j,k,m))/hvx(i,j,k,m)
						hvy(i,j,k,m) = (hvy(i,j,k,m) - hvy2(i,j,k,m))/hvy(i,j,k,m)
						hvz(i,j,k,m) = (hvz(i,j,k,m) - hvz2(i,j,k,m))/hvz(i,j,k,m)
						write(8,'(3(f9.2),9(es14.6))')ri,rj,rk, &
							qvx(i,j,k,m),qvy(i,j,k,m), &
							qvz(i,j,k,m),ovx(i,j,k,m),ovy(i,j,k,m), &
							ovz(i,j,k,m),hvx(i,j,k,m),hvy(i,j,k,m), &
							hvz(i,j,k,m)
						!
						abx = (abx - abx_2)/abx
						aby = (aby - aby_2)/aby
						abz = (abz - abz_2)/abz
						curx_all(i,j,k,m) = (curx_all(i,j,k,m) - curx_all2(i,j,k,m))/curx_all(i,j,k,m)
						cury_all(i,j,k,m) = (cury_all(i,j,k,m) - cury_all2(i,j,k,m))/cury_all(i,j,k,m)
						curz_all(i,j,k,m) = (curz_all(i,j,k,m) - curz_all2(i,j,k,m))/curz_all(i,j,k,m)
						write(9,'(3(f9.2),6(es14.6))')ri,rj,rk, &
							abx,aby,abz,curx_all(i,j,k,m), &
							cury_all(i,j,k,m),curz_all(i,j,k,m)
						!
						bsurmag = (bsurmag - bsurmag2)/bsurmag
						write(10,'(5(es14.6))')rad, bsurmag, ri, rj, rk
						!
						abx=b_equiv*(bx(i,j,k,m))
						abx2=b_equiv*(bx0(i,j,k,m))
						aby=b_equiv*(by(i,j,k,m))
						aby2=b_equiv*(by0(i,j,k,m))
						abz=b_equiv*(bz(i,j,k,m))
						abz2=b_equiv*(bz0(i,j,k,m))

						abx_2=b_equiv*(bx2(i,j,k,m))
						abx22=b_equiv*(bx02(i,j,k,m))
						aby_2=b_equiv*(by2(i,j,k,m))
						aby22=b_equiv*(by02(i,j,k,m))
						abz_2=b_equiv*(bz2(i,j,k,m))
						abz22=b_equiv*(bz02(i,j,k,m))
						!
						abx = (abx - abx_2)/abx
						aby = (aby - aby_2)/aby
						abz = (abz - abz_2)/abz
						abx2 = (abx2 - abx22)/abx2
						aby2 = (aby2 - aby22)/aby2
						abz2 = (abz2 - abz22)/abz2
						write(11,'(3(f9.2),6(es14.6))')ri,rj,rk,&
							abx,aby,abz,abx2,aby2,abz2
						!
						efldx(i,j,k,m) = (efldx(i,j,k,m) - efldx2(i,j,k,m))/efldx(i,j,k,m)
						efldy(i,j,k,m) = (efldy(i,j,k,m) - efldy2(i,j,k,m))/efldy(i,j,k,m)
						efldz(i,j,k,m) = (efldz(i,j,k,m) - efldz2(i,j,k,m))/efldz(i,j,k,m)
						write(12,'(3(f9.2),3(es14.6))')ri,rj,rk, &
							efldx(i,j,k,m),efldy(i,j,k,m),efldz(i,j,k,m)
						!
						qpara(i,j,k,m) = (qpara(i,j,k,m) - qpara2(i,j,k,m))/qpara(i,j,k,m)
						qperp(i,j,k,m) = (qperp(i,j,k,m) - qperp2(i,j,k,m))/qperp(i,j,k,m)
						qcross(i,j,k,m) = (qcross(i,j,k,m) - qcross2(i,j,k,m))/qcross(i,j,k,m)
						write(13,'(3(f9.2),10(es14.6))')ri,rj,rk, &
							qpara(i,j,k,m),qperp(i,j,k,m), &
							qcross(i,j,k,m)
!7  = plasma-box#.dat
!8  = flow-box#.dat
!9  = field-box#.dat
!10 = rad-box#.dat
!11 = model-box#.dat
!12 = efield-box#.dat
!13 = pressure-box#.dat
!14 = rotational-box#.dat
!15 = currents-box#.dat
						!
		        	enddo
		        enddo
		    enddo
			!
		    close(7)
		    close(8)
		    close(9)
		    close(10)
		    close(11)
		    close(12)
		    close(13)
		    close(14)
		    close(15)
			!
		enddo ! loop over m
	!
	endif
	if (percent.lt.100)	call loadbar(100)
	!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~**~*~
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
	write(*,*) '                Done writing ASCII files!'
	write(*,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
end
!
!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~**~*~
!~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~*~*~*~*~*~*~*~**~*~*~*~*~*~*~*~**~*~
!
!
!	*************************************************
!
!
subroutine interpolate(plas_int,plas,box,xpos,ypos,zpos, &
    nbx,nx,ny,nz)
	!
    integer, intent (in) :: nbx, nx, ny, nz, box
	!
    real, intent (in) :: plas(nbx,nx,ny,nz), xpos, ypos, zpos
    real, intent (inout) :: plas_int
	!
    integer i, j, k, m, ii, jj, kk
    real dx, dy, dz
	!
    m = box
	!
    i = int(xpos)
    dx = xpos - real(i)
	!
    j = int(ypos)
    dy = ypos - real(j)
	!
    k = int(zpos)
    dz = zpos - real(k)
	!
    ii = i + 1
    jj = j + 1
    kk = k + 1
	!
	plas_int = plas(i,j,k,m)*(1.0 - dx)*(1.0 - dy)*(1.0 - dz) &
		+ plas(ii,jj,kk,m)*dx*dy*dz &
		+ plas(ii,j,k,m)*dx*(1.0 - dy)*(1.0 - dz) &
		+ plas(i,jj,k,m)*(1.0 - dx)*dy*(1.0 - dz) &
		+ plas(i,j,k,mk)*(1.0 - dx)*(1.0 - dy)*dz &
		+ plas(ii,jj,k,m)*dx*dy*(1.0 - dz) &
		+ plas(ii,j,kk,m)*dx*(1.0 - dy)*dz &
		+ plas(i,jj,kk,m)*(1.0 - dx)*dy*dz
	!
end subroutine interpolate
!
!
!     ************************************************
!
!
subroutine set_gridpt(xpos,ypos,zpos, &
	xsat,ysat,zsat,r_xmin,r_ymin,r_zmin,xspac,re_equiv, &
	moon,xmoon_re,ymoon_re,zmoon_re)
	!
	real, intent (inout) :: xpos,ypos,zpos
	!
	real, intent (in) :: xsat,ysat,zsat,r_xmin, &
		r_ymin,r_zmin,xspac,re_equiv, &
		xmoon_re,ymoon_re,zmoon_re
	!
	integer, intent (in) :: moon
	!
	real xtmp, ytmp, ztmp, df1, df2
	!
	df2 = xspac*re_equiv
	!
	df1 = r_xmin - xsat
	xtmp = df1/df2
	xpos = abs(xtmp) + 1
	!
	df1 = r_ymin - ysat
	ytmp = df1/df2
	ypos = abs(ytmp) + 1
	!
    df1 = r_zmin - zsat
    ztmp = df1/df2
    zpos = abs(ztmp) + 1
	!
end subroutine set_gridpt
!
!
!     *************************************************
!
!
subroutine calcur(bx,by,bz,nx,ny,nz,ngrd,m,curx,cury,curz, &
	rx,ry,rz)
	!
	!     this subroutine calculates the current associated with the perturbed b field
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
				!    regular boundary conditions
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
!
!     *************************************************
!
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
			!
    		do j=2,ny-1
				jp=j+1
				jm=j-1
				!
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
	!
	!	print*, numzero
	do n=1,numzero
        i=ijzero(m,1,n)
        j=ijzero(m,2,n)
        k=ijzero(m,3,n)
		!
        eden=(qrho(i,j,k,m)/rmassq+hrho(i,j,k,m)/rmassh+ &
        	orho(i,j,k,m)/rmasso)
		!
        efldx(i,j,k)=efldx(i,j,k)+rst(i,j,k,m)*curx(i,j,k)/eden
        efldy(i,j,k)=efldy(i,j,k)+rst(i,j,k,m)*cury(i,j,k)/eden
        efldz(i,j,k)=efldz(i,j,k)+rst(i,j,k,m)*curz(i,j,k)/eden
		!
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
!
!     ************************************************
!
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
					!	eden=amax1(qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m),
					!   +                 0.0001)
		 			qden=(qrho(i,j,k,m)+0.000001)/rmassq
				    hden=(hrho(i,j,k,m)+0.000001)/rmassh
					oden=(orho(i,j,k,m)+0.000001)/rmasso
       				tden=qden+hden+oden
					!
					!	keep separate ion and current components
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
!
!     ************************************************
!
!
subroutine fnd_vel(px,py,pz,rho,vx,vy,vz,nx,ny,nz,ngrd,m)
	!
	!     converts momentum into velocity for graphics
	!
    dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd), &
	    pz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd), &
        vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
	!
	!$omp parallel do
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
!
!
!     ************************************************
!
!
subroutine fnd_pres(px,py,pz,pxy,pxz,pyz,p_para,p_perp,p_cross, &
    vvx,vvy,vvz,bxt,byt,bzt,nx,ny,nz,ngrd,m)
	!
    dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd),pz(nx,ny,nz,ngrd), &
        pxy(nx,ny,nz,ngrd),pxz(nx,ny,nz,ngrd),pyz(nx,ny,nz,ngrd), &
        p_para(nx,ny,nz,ngrd),p_perp(nx,ny,nz,ngrd),p_cross(nx,ny,nz,ngrd), &
        vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz),bxt(nx,ny,nz), &
        byt(nx,ny,nz),bzt(nx,ny,nz)
	!
    real abx,aby,abz,bmag,avx,avy,avz,vmag
    real vcbx,vcby,vcbz,vcmag
    real p_para,p_perp,p_cross
    real presx,presy,presz,presmag
    real presxy,presxz,presyz
	!
	!    find anisotropy values for species 
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
                vcmag=sqrt(vcbx**2+vcby**2+vcbz**2)+1.e-12
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
				!
            enddo
		enddo
	enddo
end
!
!
!     ************************************************
!
!
subroutine totfld(bx,bx0,btx,nx,ny,nz,ngrd,m)
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
!
!     *************************************************
!
!
subroutine tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
	!
	!      initialize static magnetic field along entire grid
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
!
!
!     *************************************************
!
!
subroutine loadbar(percent)
	integer, intent(in)	:: percent
	integer, parameter	:: mult = 6	! (mult*10)+2 is number of characters used
	character,parameter	:: esc  = char(27)

	integer :: loop		= 1
	integer :: barperc
    integer :: totalbar	= mult*10
	integer :: progress
	integer :: remain

	character*10 :: full10 = '==========', empty10 = '          '
	character*60 :: full = '', empty = ''

	character*6 :: fmtbeg
	character*6 :: fmtend
	character*3 :: progfmt
	character*3 :: remfmt

	barperc = percent

	if(barperc.lt.0) barperc = 0

	if(barperc.ge.100) then
		barperc = 100
		progress = totalbar
		remain = 0
	else
		progress = barperc*totalbar/100
		remain = totalbar - progress - 1
	endif

	do n=1, mult
		full	= trim(full)//full10
		empty	= trim(empty)//empty10
	enddo

	write(progfmt,'(I3)') progress
	fmtbeg   = '(A'//trim(adjustl(progfmt))//'$)'
	write(remfmt,'(I3)') remain
	fmtend = '(A'//trim(adjustl(remfmt))//'$)'

	write(*,'(A1$)') '['
	write(*,fmtbeg) full
	if(progress.lt.totalbar)	write(*,'(A1$)') '>'
	write(*,fmtend) empty
	write(*,'(A1)') ']'
	write(*,'(A,I3,A)') 'Working,', barperc, '% done.'

	if(barperc.lt.100) write(*,'(A,A$)') esc, '[2A'

	return
end subroutine
