!
!     this is a 3-d modified three fluid simulation using the
!           electrons : arrays starting with e
!           solar wind : arrays starting wit q
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
!       grid within grid size nt = 2*ngrd
!       ncts is the size of the data array for the imf data file
!
program multifluid
    parameter (nx=121,ny=121,nz=61,ngrd=7, &
    mbndry=1,msrf=2000,mmid=1500,mzero=5000, &
    nx_n=49,ny_n=49,nz_n=49,ngrd_n=6,mbndry_n=1, &
    msrf_n=1400,mmid_n=1000,mzero_n=2800, &
    ncraft=30,ncts=281)
    !
    !      graphics parameters:muvwp2=amax(mx,my,mz)+2,mz2=(mz-1)/2+1
    !
    parameter (mx=61,my=61,mz=31,muvwp2=63,mz2=16)
    parameter (mx_n=25,my_n=25,mz_n=25,muvwp2_n=28,mz2_n=13)
    !
    common /space/vvx(nx,ny,nz),vvy(nx,ny,nz),vvz(nx,ny,nz), &
    tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
    evx(nx,ny,nz),evy(nx,ny,nz),evz(nx,ny,nz)
    common /space_n/vvx_n(nx_n,ny_n,nz_n),vvy_n(nx_n,ny_n,nz_n), &
    vvz_n(nx_n,ny_n,nz_n),tvx_n(nx_n,ny_n,nz_n), &
    tvy_n(nx_n,ny_n,nz_n),tvz_n(nx_n,ny_n,nz_n), &
    evx_n(nx_n,ny_n,nz_n),evy_n(nx_n,ny_n,nz_n), &
    evz_n(nx_n,ny_n,nz_n)
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    common /moon/xmoon,ymoon,zmoon,rmoon,b0_moon, &
    xdip_moon,ydip_moon,zdip_moon,offset
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
    grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd), &
    xspac(ngrd)
    dimension grd_xmin_n(ngrd_n),grd_xmax_n(ngrd_n), &
    grd_ymin_n(ngrd_n),grd_ymax_n(ngrd_n), &
    grd_zmin_n(ngrd_n),grd_zmax_n(ngrd_n), &
    xspac_n(ngrd_n)
    !
    dimension grd_time_n(ngrd_n),grd_dx_n(ngrd_n),grd_dy_n(ngrd_n), &
    grd_vx_n(ngrd_n),grd_vy_n(ngrd_n)
    !
    !      grid limits now set by grd_min grd_max arrays
    !      xspac is the relative grid spacing relative to inner grid system
    !      xcraft is the actual position of the spacecraft in re
    !          4th dimension of the actual time
    !      zcraft is the future position of the spacecraft in re
    !      rcraft is the position of the spacecraft for which
    !           imf is reference. no alteration from boundary conditions applied
    !
    real xcraft(4,ncraft),zcraft(4,ncraft),rcraft(3)
    !
    !     physics plasma quantities :main grid
    !
    real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    qpresxy(nx,ny,nz,ngrd),qpresxz(nx,ny,nz,ngrd), &
    qpresyz(nx,ny,nz,ngrd), &
    !
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    hpresxy(nx,ny,nz,ngrd),hpresxz(nx,ny,nz,ngrd), &
    hpresyz(nx,ny,nz,ngrd), &
    !
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    opresxy(nx,ny,nz,ngrd),opresxz(nx,ny,nz,ngrd), &
    opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd)
    !
    !     physics plasma quantities:hires grid
    !
    real bx_n(nx_n,ny_n,nz_n,ngrd_n),by_n(nx_n,ny_n,nz_n,ngrd_n), &
    bz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    qpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    hpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    hpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    opx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opz_n(nx_n,ny_n,nz_n,ngrd_n), &
    orho_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    opresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    epres_n(nx_n,ny_n,nz_n,ngrd_n)
    !
    !
    !     work arrays for runge-kutta and smothing: main grid
    !
    real, allocatable, dimension(:,:,:,:) :: &
    oldbx,oldby,oldbz, &
    oldqrho,oldqpx,oldqpy,oldqpz, &
    oldqpresx,oldqpresy,oldqpresz, &
    oldqpresxy,oldqpresxz,oldqpresyz, &
    oldhrho,oldhpx,oldhpy,oldhpz, &
    oldhpresx,oldhpresy,oldhpresz, &
    oldhpresxy,oldhpresxz,oldhpresyz, &
    oldorho,oldopx,oldopy,oldopz, &
    oldopresx,oldopresy,oldopresz, &
    oldopresxy,oldopresxz,oldopresyz, &
    oldepres
    !
    real, allocatable, dimension(:,:,:,:) :: &
    wrkbx,wrkby,wrkbz, &
    wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
    wrkqpresx,wrkqpresy,wrkqpresz, &
    wrkqpresxy,wrkqpresxz,wrkqpresyz, &
    wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
    wrkhpresx,wrkhpresy,wrkhpresz, &
    wrkhpresxy,wrkhpresxz,wrkhpresyz, &
    wrkorho,wrkopx,wrkopy,wrkopz, &
    wrkopresx,wrkopresy,wrkopresz, &
    wrkopresxy,wrkopresxz,wrkopresyz, &
    wrkepres
    !
    !     work arrays for runge-kutta and smothing: hires grid
    !
    real, allocatable, dimension(:,:,:,:) :: &
    oldbx_n,oldby_n,oldbz_n, &
    oldqrho_n,oldqpx_n,oldqpy_n,oldqpz_n, &
    oldqpresx_n,oldqpresy_n,oldqpresz_n, &
    oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
    oldhrho_n,oldhpx_n,oldhpy_n,oldhpz_n, &
    oldhpresx_n,oldhpresy_n,oldhpresz_n, &
    oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
    oldorho_n,oldopx_n,oldopy_n,oldopz_n, &
    oldopresx_n,oldopresy_n,oldopresz_n, &
    oldopresxy_n,oldopresxz_n,oldopresyz_n, &
    oldepres_n
    !
    real, allocatable, dimension(:,:,:,:) :: &
    wrkbx_n,wrkby_n,wrkbz_n, &
    wrkqrho_n,wrkqpx_n,wrkqpy_n,wrkqpz_n, &
    wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
    wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
    wrkhrho_n,wrkhpx_n,wrkhpy_n,wrkhpz_n, &
    wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
    wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
    wrkorho_n,wrkopx_n,wrkopy_n,wrkopz_n, &
    wrkopresx_n,wrkopresy_n,wrkopresz_n, &
    wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
    wrkepres_n
    !
    !
    !     unperturbed quantities
    !
    real bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
    real bx0_n(nx_n,ny_n,nz_n,ngrd_n),by0_n(nx_n,ny_n,nz_n,ngrd_n), &
    bz0_n(nx_n,ny_n,nz_n,ngrd_n)
    
    real efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz),btot(nx,ny,nz), &
    resistive(nx,ny,nz,mbndry)
    real efldx_n(nx_n,ny_n,nz_n),efldy_n(nx_n,ny_n,nz_n), &
    efldz_n(nx_n,ny_n,nz_n),curx_n(nx_n,ny_n,nz_n), &
    cury_n(nx_n,ny_n,nz_n),curz_n(nx_n,ny_n,nz_n), &
    bsx_n(nx_n,ny_n,nz_n),bsy_n(nx_n,ny_n,nz_n), &
    bsz_n(nx_n,ny_n,nz_n),btot_n(nx_n,ny_n,nz_n), &
    resistive_n(nx_n,ny_n,nz_n,mbndry_n)
    !
    real lunar_rad,lunar_dist
    !
    !      variable time step arrays
    !
    real t_old(ngrd),t_new(ngrd),t_step(ngrd),t_stepnew(ngrd)
    real t_old_n(ngrd_n),t_new_n(ngrd_n),t_step_n(ngrd_n), &
    t_stepnew_n(ngrd)
    !
    !     boundary condition arrays
    !
    dimension bxf(ny,nz),byf(ny,nz),bzf(ny,nz), &
    rhof(ny,nz),svxf(ny,nz),svyf(ny,nz),svzf(ny,nz)
    dimension bxp(ny,nz),byp(ny,nz),bzp(ny,nz), &
    rhop(ny,nz),svxp(ny,nz),svyp(ny,nz),svzp(ny,nz)
    dimension future(ny,nz),past(ny,nz), &
    bfld(ncts,4),rplas(ncts),svel(ncts,3)
    integer ncount(ny,nz)
    !
    real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz), &
    tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2), &
    cross(ny,nz),along(nx,nz),flat(nx,ny)
    !
    real tx_n(mx_n,my_n,mz_n),ty_n(mx_n,my_n,mz_n), &
    tz_n(mx_n,my_n,mz_n),tg1_n(mx_n,my_n,mz_n), &
    tg2_n(mx_n,my_n,mz2_n),tt_n(mx_n,my_n,mz_n), &
    work_n(muvwp2_n,muvwp2_n), &
    cross_n(ny_n,nz_n),along_n(nx_n,nz_n), &
    flat_n(nx_n,ny_n)
    !
    character*8 wd1,wd2,wd3,wd4
    character*8 label
    character*15 title
    !
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
    ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    real    parm_srf(mbndry,7,msrf),parm_mid(mbndry,7,mmid), &
    parm_zero(mbndry,7,mzero)
    !
    integer ijsrf_n(mbndry_n,3,msrf_n),ijmid_n(mbndry_n,3,mmid_n), &
    ijzero_n(mbndry_n,3,mzero_n)
    integer numsrf_n(mbndry_n),nummid_n(mbndry_n),numzero_n(mbndry_n)
    real   parm_srf_n(mbndry_n,7,msrf_n), &
    parm_mid_n(mbndry_n,7,mmid_n), &
    parm_zero_n(mbndry_n,7,mzero_n)
    !
    logical start,add_dip,ringo,update,save_dat,write_dat, &
    spacecraft,tilting,warp,reload,divb_lores,divb_hires, &
    yes_step,yes_step_n,grid_reset,isotropic
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
    difrho,difpxyz,diferg, &
    chirho_n,chipxyz_n,chierg_n, &
    difrho_n,difpxyz_n,diferg_n
    namelist/subgrid/main_grd_n
    !
    !     allocate arrays
    !
    !     work arrays for runge-kutta and smoothing
    
    allocate(oldbx(nx,ny,nz,ngrd),oldby(nx,ny,nz,ngrd), &
    oldbz(nx,ny,nz,ngrd), &
    !
    oldqpx(nx,ny,nz,ngrd),oldqpy(nx,ny,nz,ngrd), &
    oldqpz(nx,ny,nz,ngrd),oldqrho(nx,ny,nz,ngrd), &
    oldqpresx(nx,ny,nz,ngrd),oldqpresy(nx,ny,nz,ngrd), &
    oldqpresz(nx,ny,nz,ngrd), oldqpresxy(nx,ny,nz,ngrd), &
    oldqpresxz(nx,ny,nz,ngrd),oldqpresyz(nx,ny,nz,ngrd), &
    !
    oldhpx(nx,ny,nz,ngrd),oldhpy(nx,ny,nz,ngrd), &
    oldhpz(nx,ny,nz,ngrd),oldhrho(nx,ny,nz,ngrd), &
    oldhpresx(nx,ny,nz,ngrd),oldhpresy(nx,ny,nz,ngrd), &
    oldhpresz(nx,ny,nz,ngrd), oldhpresxy(nx,ny,nz,ngrd), &
    oldhpresxz(nx,ny,nz,ngrd),oldhpresyz(nx,ny,nz,ngrd), &
    !
    oldopx(nx,ny,nz,ngrd),oldopy(nx,ny,nz,ngrd), &
    oldopz(nx,ny,nz,ngrd),oldorho(nx,ny,nz,ngrd), &
    oldopresx(nx,ny,nz,ngrd),oldopresy(nx,ny,nz,ngrd), &
    oldopresz(nx,ny,nz,ngrd), oldopresxy(nx,ny,nz,ngrd), &
    oldopresxz(nx,ny,nz,ngrd),oldopresyz(nx,ny,nz,ngrd), &
    !
    oldepres(nx,ny,nz,ngrd))
    !
    allocate(wrkbx(nx,ny,nz,ngrd),wrkby(nx,ny,nz,ngrd), &
    wrkbz(nx,ny,nz,ngrd), &
    !
    wrkqpx(nx,ny,nz,ngrd),wrkqpy(nx,ny,nz,ngrd), &
    wrkqpz(nx,ny,nz,ngrd),wrkqrho(nx,ny,nz,ngrd), &
    wrkqpresx(nx,ny,nz,ngrd),wrkqpresy(nx,ny,nz,ngrd), &
    wrkqpresz(nx,ny,nz,ngrd), wrkqpresxy(nx,ny,nz,ngrd), &
    wrkqpresxz(nx,ny,nz,ngrd),wrkqpresyz(nx,ny,nz,ngrd), &
    !
    wrkhpx(nx,ny,nz,ngrd),wrkhpy(nx,ny,nz,ngrd), &
    wrkhpz(nx,ny,nz,ngrd),wrkhrho(nx,ny,nz,ngrd), &
    wrkhpresx(nx,ny,nz,ngrd),wrkhpresy(nx,ny,nz,ngrd), &
    wrkhpresz(nx,ny,nz,ngrd), wrkhpresxy(nx,ny,nz,ngrd), &
    wrkhpresxz(nx,ny,nz,ngrd),wrkhpresyz(nx,ny,nz,ngrd), &
    !
    wrkopx(nx,ny,nz,ngrd),wrkopy(nx,ny,nz,ngrd), &
    wrkopz(nx,ny,nz,ngrd),wrkorho(nx,ny,nz,ngrd), &
    wrkopresx(nx,ny,nz,ngrd),wrkopresy(nx,ny,nz,ngrd), &
    wrkopresz(nx,ny,nz,ngrd), wrkopresxy(nx,ny,nz,ngrd), &
    wrkopresxz(nx,ny,nz,ngrd),wrkopresyz(nx,ny,nz,ngrd), &
    !
    wrkepres(nx,ny,nz,ngrd))
    !
    !     allocate arrays for hires section
    !
    !     work arrays for runge-kutta and smoothing
    
    allocate(oldbx_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldby_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldbz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldqpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    oldhpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldhpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    oldopx_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldorho_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    oldopresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    oldepres_n(nx_n,ny_n,nz_n,ngrd_n))
    !
    allocate(wrkbx_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkby_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkbz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkqpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    wrkhpx_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhrho_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkhpresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    wrkopx_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkorho_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopresx_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopresy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopresz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopresxy_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopresxz_n(nx_n,ny_n,nz_n,ngrd_n), &
    wrkopresyz_n(nx_n,ny_n,nz_n,ngrd_n), &
    !
    wrkepres_n(nx_n,ny_n,nz_n,ngrd_n))
    !
    !      open input data file
    !
    open(3,file='fluxes.dat',status='unknown',form='formatted')
    open(5,file='hiresin',status='old',form='formatted')
    !     open(6,file='rmpd3dout',status='unknown',form='formatted')
    open(7,file='grid.dat',status='unknown',form='formatted')
    open(8,file='cur.dat',status='unknown',form='unformatted')
    open(9,file='pot.dat',status='unknown',form='unformatted')
    open(10,file='conc.dat',status='unknown',form='formatted')
    !      open ncargraphics
    !
    !
    call opngks
    call gsclip(0)
    !
    !     set color table
    !
    call cpclrs
    !
    call gselnt(0)
    call gsplci(1)
    call wtstr(.4,.975,'3d mutant code',2,0,0)
    !
    !
    !     read input parameters
    !
    read(5,option)
    write(6,option)
    read(5,earth)
    write(6,earth)
    read(5,speeds)
    write(6,speeds)
    read(5,windy)
    write(6,windy)
    read(5,lunar)
    write(6,lunar)
    read(5,physical)
    write(6,physical)
    read(5,smooth)
    write(6,smooth)
    read(5,subgrid)
    write(6,subgrid)
    !
    !     output to test whether the parameters are the actual ones you want
    !
    write(6,option)
    write(6,earth)
    write(6,speeds)
    write(6,windy)
    write(6,physical)
    write(6,smooth)
    !
    do i=1,ngrd
        read(5,*)grd_xmin(i),grd_xmax(i),grd_ymin(i),grd_ymax(i), &
        grd_zmin(i),grd_zmax(i),xspac(i)
        write(6,*)grd_xmin(i),grd_xmax(i),grd_ymin(i),grd_ymax(i), &
        grd_zmin(i),grd_zmax(i),xspac(i)
        ix=1+(grd_xmax(i)-grd_xmin(i))/xspac(i)
        iy=1+(grd_ymax(i)-grd_ymin(i))/xspac(i)
        iz=1+(grd_zmax(i)-grd_zmin(i))/xspac(i)
        if((ix.ne.nx).or.(iy.ne.ny).or.(iz.ne.nz))then
            write(6,*)' warning: sizes',ix,iy,iz,nx,ny,nz
            stop
        endif
    enddo
    !
    !
    !      check if subgrids mate properly to larger grid
    !
    do m=1,ngrd-1
        mm=m+1
        ai=1.+(grd_xmin(m)-grd_xmin(mm))/xspac(mm)
        i=ai
        dx=ai-i
        if(abs(dx).gt.0.001)then
            write(6,*)'main grd: xmin dont match',m,mm,i,ai
            stop
        endif
        !
        aj=1.+(grd_ymin(m)-grd_ymin(mm))/xspac(mm)
        j=aj
        dy=aj-j
        if(abs(dy).gt.0.001)then
            write(6,*)'main grd: ymin dont match',m,mm,j,aj
            stop
        endif
        !
        ak=1.+(grd_zmin(m)-grd_zmin(mm))/xspac(mm)
        k=ak
        dz=ak-k
        if(abs(dz).gt.0.001)then
            write(6,*)'main grd: zmin dont match',m,mm,k,ak
            stop
        endif
        !
    enddo
    !
    !      _n grid system
    !
    do m=1,ngrd_n
        read(5,*)grd_xmin_n(m),grd_xmax_n(m),grd_ymin_n(m),grd_ymax_n(m), &
        grd_zmin_n(m),grd_zmax_n(m),xspac_n(m)
        write(6,*)grd_xmin_n(m),grd_xmax_n(m),grd_ymin_n(m), &
        grd_ymax_n(m),grd_zmin_n(m),grd_zmax_n(m),xspac_n(m)
        ix=1+(grd_xmax_n(m)-grd_xmin_n(m))/xspac_n(m)
        iy=1+(grd_ymax_n(m)-grd_ymin_n(m))/xspac_n(m)
        iz=1+(grd_zmax_n(m)-grd_zmin_n(m))/xspac_n(m)
        if((ix.ne.nx_n).or.(iy.ne.ny_n).or.(iz.ne.nz_n))then
            write(6,*)' warning: _n sizes',m,ix,iy,iz,nx_n,ny_n,nz_n
            stop
        endif
    enddo
    !
    !     check to determine whether the two systems mate
    !
    m=main_grd_n
    m_n=ngrd_n
    !
    i=1+(grd_xmin_n(m_n)-grd_xmin(m))/xspac(m)
    ax= grd_xmin(m)+(i-1)*xspac(m)
    if(ax.ne.grd_xmin_n(m_n))then
        write(6,*)'warning: twos grid dont match: xmin', &
        grd_xmin_n(m_n),ax,m,m_n
        stop
    endif
    !
    i=1+(grd_xmax_n(m_n)-grd_xmin(m))/xspac(m)
    ax= grd_xmin(m)+(i-1)*xspac(m)
    if(ax.ne.grd_xmax_n(m_n))then
        write(6,*)'warning: twos grid dont match: xmax', &
        grd_xmax_n(m_n),ax,m,m_n
        stop
    endif
    !
    j=1+(grd_ymin_n(m_n)-grd_ymin(m))/xspac(m)
    ay= grd_ymin(m)+(j-1)*xspac(m)
    if(ay.ne.grd_ymin_n(m_n))then
        write(6,*)'warning: twos grid dont match: ymin', &
        grd_ymin_n(m_n),ay,m,m_n
        stop
    endif
    !
    j=1+(grd_ymax_n(m_n)-grd_ymin(m))/xspac(m)
    ay= grd_ymin(m)+(j-1)*xspac(m)
    if(ay.ne.grd_ymax_n(m_n))then
        write(6,*)'warning: twos grid dont match: ymax', &
        grd_ymax_n(m_n),ay,m,m_n
        stop
    endif
    !
    k=1+(grd_zmin_n(m_n)-grd_zmin(m))/xspac(m)
    az= grd_zmin(m)+(k-1)*xspac(m)
    if(az.ne.grd_zmin_n(m_n))then
        write(6,*)'warning: twos grid dont match: zmin', &
        grd_zmin_n(m_n),az,m,m_n
        stop
    endif
    !
    k=1+(grd_zmax_n(m_n)-grd_zmin(m))/xspac(m)
    az= grd_zmin(m)+(k-1)*xspac(m)
    if(az.ne.grd_zmax_n(m_n))then
        write(6,*)'warning: twos grid dont match: zmax', &
        grd_zmax_n(m_n),az,m,m_n
        stop
    endif
    !
    do m=1,ngrd_n-1    ! start m loop
        mm=m+1
        ai=1.+(grd_xmin_n(m)-grd_xmin_n(mm))/xspac_n(mm)
        i=ai
        dx=ai-i
        if(abs(dx).gt.0.001)then
            write(6,*)'sub grd: xmin dont match',m,mm,i,ai
            stop
        endif
        !
        aj=1.+(grd_ymin_n(m)-grd_ymin_n(mm))/xspac_n(mm)
        j=aj
        dy=aj-j
        if(abs(dy).gt.0.001)then
            write(6,*)'sub grd: ymin dont match',m,mm,j,aj
            stop
        endif
        !
        ak=1.+(grd_zmin_n(m)-grd_zmin_n(mm))/xspac_n(mm)
        k=ak
        dz=ak-k
        if(abs(dz).gt.0.001)then
            write(6,*)'sub grd: zmin dont match',m,mm,k,ak
            stop
        endif
        !
    enddo  ! end m loop
    !
    !
    !     write important data to graphics file
    !
    write(wd4,'(f6.3)')stepsz
    !
    title='stepsz = '//wd4
    call wtstr(.75,.85,title,1,0,0)
    !
    write(wd1,'(f6.3)')xdip
    write(wd2,'(f6.3)')ydip
    write(wd3,'(f6.3)')zdip
    !
    title='xdip = '//wd1
    call wtstr(.15,.82,title,1,0,0)
    title='ydip = '//wd2
    call wtstr(.35,.82,title,1,0,0)
    title='zdip = '//wd3
    call wtstr(.55,.82,title,1,0,0)
    !
    
    write(wd1,'(f5.1)')rearth
    !
    title='rearth = '//wd1
    call wtstr(.15,.79,title,1,0,0)
    !
    !     calculate effective magnetic field strength
    !
    erho=den_earth*rmassh
    b01=alf_inner1*sqrt(erho)*rearth**3
    b02=alf_inner2*sqrt(erho)*rearth**3
    alf_lim=6.00*alf_inner1
    b0=b01
    write(6,*)'b0',b0
    bmax=b02
    delb0=(b02-b01)/tmax
    !
    write(wd1,'(f6.3)')cs_inner
    write(wd2,'(f6.3)')alf_inner1
    write(wd3,'(f6.3)')alf_inner2
    write(wd4,'(f5.1)')den_earth
    !
    title='cs_inner = '//wd1
    call wtstr(.15,.76,title,1,0,0)
    title='alf_inner1= '//wd2
    call wtstr(.35,.76,title,1,0,0)
    title='alf_inner2= '//wd3
    call wtstr(.55,.76,title,1,0,0)
    title='den_earth = '//wd4
    call wtstr(.75,.76,title,1,0,0)
    !
    
    write(wd1,'(f6.3)')o_conc
    write(wd2,'(f6.3)')gravity
    write(wd3,'(f6.3)')rmasso
    !
    title='o_conc = '//wd1
    call wtstr(.15,.73,title,1,0,0)
    title='gravity= '//wd2
    call wtstr(.35,.73,title,1,0,0)
    title='rmasso= '//wd2
    call wtstr(.35,.73,title,1,0,0)
    !
    write(wd1,'(f6.3)')rho_wind1
    write(wd2,'(f6.3)')rho_wind2
    write(wd3,'(f6.3)')vx_wind1
    write(wd4,'(f6.3)')vx_wind2
    !
    title='rho_wind1= '//wd1
    call wtstr(.15,.7,title,1,0,0)
    title='rho_wind2= '//wd2
    call wtstr(.35,.7,title,1,0,0)
    title='vx_wind1 = '//wd3
    call wtstr(.55,.7,title,1,0,0)
    title='vx_wind2 = '//wd4
    call wtstr(.75,.7,title,1,0,0)
    
    !
    write(wd1,'(f6.3)')vy_wind1
    write(wd2,'(f6.3)')vy_wind2
    write(wd3,'(f6.3)')vz_wind1
    write(wd4,'(f6.3)')vz_wind2
    !
    title='vy_wind1= '//wd1
    call wtstr(.15,.67,title,1,0,0)
    title='vy_wind2= '//wd2
    call wtstr(.35,.67,title,1,0,0)
    title='vz_wind1 = '//wd3
    call wtstr(.55,.67,title,1,0,0)
    title='vz_wind2 = '//wd4
    call wtstr(.75,.67,title,1,0,0)
    !
    write(wd1,'(f6.3)')alfx_wind1
    write(wd2,'(f6.3)')alfx_wind2
    write(wd3,'(f6.3)')alfy_wind1
    write(wd4,'(f6.3)')alfy_wind2
    
    !
    title='alfx1 = '//wd1
    call wtstr(.15,.64,title,1,0,0)
    title='alfx2 = '//wd2
    call wtstr(.35,.64,title,1,0,0)
    title='alfy1 = '//wd3
    call wtstr(.55,.64,title,1,0,0)
    title='alfy2 = '//wd4
    call wtstr(.75,.64,title,1,0,0)
    !
    write(wd3,'(f6.3)')alfz_wind1
    write(wd4,'(f6.3)')alfz_wind2
    title='alfz1 = '//wd3
    call wtstr(.55,.61,title,1,0,0)
    title='alfz2 = '//wd4
    call wtstr(.75,.61,title,1,0,0)
    
    write(wd1,'(f5.1)')re_wind
    write(wd2,'(f5.0)')reynolds
    write(wd3,'(f5.0)')resist
    !     re_wind sets raduius from earth where initial wind placed
    !     reynolds coefficient for surface currents
    !     resist equivalent if you wish to run anomalous resistivity
    !     bfrac determines the percentage of the tangential magnetic
    !     field allowed at earth's surface
    !
    title='re_wind = '//wd1
    call wtstr(.15,.58,title,1,0,0)
    title='reynolds = '//wd2
    call wtstr(.35,.58,title,1,0,0)
    title='resist = '//wd3
    call wtstr(.55,.58,title,1,0,0)
    !
    write(wd1,'(f6.3)')bfrac
    write(wd2,'(f6.3)')vfrac
    title='bfrac = '//wd1
    call wtstr(.35,.55,title,1,0,0)
    title='vfrac = '//wd2
    call wtstr(.55,.55,title,1,0,0)
    !
    write(wd1,'(f6.3)')chirho
    write(wd2,'(f6.3)')chipxyz
    write(wd3,'(f6.3)')chierg
    !
    title='chirho = '//wd1
    call wtstr(.15,.52,title,1,0,0)
    title='chipxyz = '//wd2
    call wtstr(.35,.52,title,1,0,0)
    title='chierg = '//wd3
    call wtstr(.55,.52,title,1,0,0)
    !
    write(wd1,'(f6.3)')difrho
    write(wd2,'(f6.3)')difpxyz
    write(wd3,'(f6.3)')diferg
    !
    title='chirho = '//wd1
    call wtstr(.15,.49,title,1,0,0)
    title='chipxyz = '//wd2
    call wtstr(.35,.49,title,1,0,0)
    title='chierg = '//wd3
    call wtstr(.55,.49,title,1,0,0)
    !
    write(wd1,'(f4.1)')tilt1
    title='tilt1 = '//wd1
    call wtstr(.15,.40,title,1,0,0)
    write(wd1,'(f4.1)')tilt2
    title='tilt2 = '//wd1
    call wtstr(.30,.40,title,1,0,0)
    tilt=tilt1
    sin_tilt=sin(tilt*.0174533)
    cos_tilt=cos(tilt*.0174533)
    dtilt=(tilt2-tilt1)/tmax
    !
    !     jupiter parameters
    !
    planet_rad=71000.   !km
    planet_per=9.7     !hr
    lunar_dist=5.9 + 0.075     !orbital radii + torus infall allowance
    torus_rad=1.0
    v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
    r_rot=30.0   !re where corotation stops
    !
    !     earth parameters
    !     planet_rad=6371.   !km
    !     planet_per=24.     !hr
    !     lunar_rad=60.      !orbital radii
    !     v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
    !     r_rot=10.0   !re where corotation stops
    !
    !     saturn parameters
    !
    !     planet_rad=60268.   !km
    !     planet_per=10.6     !hr
    !     lunar_rad=60.0      !orbital radii
    !     v_rot=6.2832*planet_rad/(planet_per*3600.)/v_equiv  ! normalized units
    !     r_rot=7.5   !re where corotation stops
    !
    !     lunar stuff: moon radius 1738. titan radius 2575
    !
    lunar_rad=1560.*1.25  !km start at exobase at 1.25 rt
    rmoon=(lunar_rad/planet_rad)/re_equiv   ! in grid points
    !
    r_orbit=orbit_moon/re_equiv ! grid pts
    ut_orbit=42.46       !hrs set articifically fast for titan
    v_orbit=(orbit_moon*planet_rad*2.*3.1414)/(ut_orbit*3600.) &
    /v_equiv    !sim units
    !
    write(6,*)'moon orbit parms',ut_orbit,v_orbit
    !
    !     spacecraft stuff
    !
    !     gravity in m/s**2 at earths surface !need t_equiv in normalization
    t=0.
    t_equiv=planet_rad*re_equiv/v_equiv
    grav=gravity*(planet_rad*re_equiv*1000.)/(1000.*v_equiv)**2
    ut=utstart
    !
    nrot=ut/planet_per
    rot_hrs=ut-nrot*planet_per
    rot_angle=6.2832*rot_hrs/planet_per
    d_min=0.1
    !
    !      ionospheric parameters
    !
    old_o_conc=o_conc
    o_conc_min=0.01
    o_conc_max=1.0
    cur_min=0.75
    cur_max=20.0
    !
    !
    !     initial position of spacecraft in re but simulation directions
    !        wind :
    xcraft(1,1)=-126.07
    xcraft(2,1)=-0.6
    xcraft(3,1)=-0.00
    !
    !      reference craft
    !
    rcraft(1)=xcraft(1,1)
    rcraft(2)=xcraft(2,1)
    rcraft(3)=xcraft(3,1)
    !
    !      io1
    !
    xcraft(1,2)=-0.00
    xcraft(2,2)=-6.10
    xcraft(3,2)=0.0003
    !
    !      io2
    !
    xcraft(1,3)=-0.00
    xcraft(2,3)=6.00
    xcraft(3,3)=0.0003
    !
    !      io3
    !
    xcraft(1,4)=5.95
    xcraft(2,4)=0.0
    xcraft(3,4)=0.0003
    !
    !      io4
    !
    xcraft(1,5)=-5.95
    xcraft(2,5)=-0.00
    xcraft(3,5)=0.0003
    !
    !      europa1
    !
    xcraft(1,6)=0.05
    xcraft(2,6)=9.4
    xcraft(3,6)=-0.0
    !
    !      europa2
    !
    xcraft(1,7)=0.0001
    xcraft(2,7)=-9.4
    xcraft(3,7)=0.0
    !
    !      europa3
    !
    xcraft(1,8)=9.4
    xcraft(2,8)=0.0
    xcraft(3,8)=-0.0
    !
    !      europa4
    !
    xcraft(1,9)=-9.4
    xcraft(2,9)=0.0
    xcraft(3,9)=0.0
    !
    !      ganymede1
    !
    xcraft(1,10)=10.6
    xcraft(2,10)=10.6
    xcraft(3,10)=0.
    !
    !      ganymede2
    !
    xcraft(1,11)=10.6
    xcraft(2,11)=-10.6
    xcraft(3,11)=0.
    !
    !      ganymede3
    !
    xcraft(1,12)=-10.6
    xcraft(2,12)=10.6
    xcraft(3,12)=0.
    !
    !      ganymede4
    !
    xcraft(1,13)=-10.6
    xcraft(2,13)=-10.6
    xcraft(3,13)=0.
    !
    !      calisto
    !
    xcraft(1,14)=26.7
    xcraft(2,14)=00.0
    xcraft(3,14)=0.
    !
    !      tail1
    !
    xcraft(1,15)=50.
    xcraft(2,15)=00.0
    xcraft(3,15)=0.
    !
    !      tail2
    !
    xcraft(1,16)=100.
    xcraft(2,16)=0.0
    xcraft(3,16)=0.
    !
    !      tail 3
    !
    xcraft(1,17)=200.
    xcraft(2,17)=00.0
    xcraft(3,17)=0.
    !
    !      tail 4
    !
    xcraft(1,18)=400.
    xcraft(2,18)=00.0
    xcraft(3,18)=0.
    !
    !
    !      tail5
    !
    xcraft(1,19)=50.
    xcraft(2,19)=10.0
    xcraft(3,19)=0.
    !
    !      tail 6
    !
    xcraft(1,20)=50.
    xcraft(2,20)=20.0
    xcraft(3,20)=0.
    !
    !      tail 7
    !
    xcraft(1,21)=50.
    xcraft(2,21)=15.0
    xcraft(3,21)=10.
    !
    !      tail 8
    !
    xcraft(1,22)=50.
    xcraft(2,22)=15.0
    xcraft(3,22)=-10.
    !
    !
    !      tail 9
    !
    xcraft(1,23)=100.
    xcraft(2,23)=20.0
    xcraft(3,23)=0.
    !
    !      tail 10
    !
    xcraft(1,24)=100.
    xcraft(2,24)=40.0
    xcraft(3,24)=0.
    !
    !      tail 11
    !
    xcraft(1,25)=100.
    xcraft(2,25)=20.0
    xcraft(3,25)=20.
    !
    !      tail 12
    !
    xcraft(1,26)=100.
    xcraft(2,26)=20.0
    xcraft(3,26)=-20.
    !
    !      tail 13
    !
    xcraft(1,27)=200.
    xcraft(2,27)=40.0
    xcraft(3,27)=0.
    !
    !      tail 14
    !
    xcraft(1,28)=200.
    xcraft(2,28)=80.0
    xcraft(3,28)=0.
    !
    !      tail 15
    !
    xcraft(1,29)=200.
    xcraft(2,29)=80.0
    xcraft(3,29)=40.
    !
    !      tail 16
    !
    xcraft(1,30)=200.
    xcraft(2,30)=80.0
    xcraft(3,30)=-40.
    !
    zut=ut-.01
    rut=zut
    vut=rut
    do n=1,ncraft
        xcraft(4,n)=zut
    enddo
    !
    call limcraft(xcraft,ncraft,re_equiv,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    do n=1,ncraft
        do i=1,4
            zcraft(i,n)=xcraft(i,n)
        enddo
    enddo
    !
    open(51,file='wind.dat',status='unknown',form='formatted')
    open(52,file='io1.dat',status='unknown',form='formatted')
    open(53,file='io2.dat',status='unknown',form='formatted')
    open(54,file='io3.dat',status='unknown',form='formatted')
    open(55,file='io4.dat',status='unknown',form='formatted')
    open(56,file='europa1.dat',status='unknown',form='formatted')
    open(57,file='europa2.dat',status='unknown',form='formatted')
    open(58,file='europa3.dat',status='unknown',form='formatted')
    open(59,file='europa4.dat',status='unknown',form='formatted')
    open(60,file='ganymede1.dat',status='unknown',form='formatted')
    open(61,file='ganymede2.dat',status='unknown',form='formatted')
    open(62,file='ganymede3.dat',status='unknown',form='formatted')
    open(63,file='ganymede4.dat',status='unknown',form='formatted')
    open(64,file='callisto.dat',status='unknown',form='formatted')
    open(65,file='tail01.dat',status='unknown',form='formatted')
    open(66,file='tail02.dat',status='unknown',form='formatted')
    open(67,file='tail03.dat',status='unknown',form='formatted')
    open(68,file='tail04.dat',status='unknown',form='formatted')
    open(69,file='tail05.dat',status='unknown',form='formatted')
    open(70,file='tail06.dat',status='unknown',form='formatted')
    open(71,file='tail07.dat',status='unknown',form='formatted')
    open(72,file='tail08.dat',status='unknown',form='formatted')
    open(73,file='tail09.dat',status='unknown',form='formatted')
    open(74,file='tail10.dat',status='unknown',form='formatted')
    open(75,file='tail11.dat',status='unknown',form='formatted')
    open(76,file='tail12.dat',status='unknown',form='formatted')
    open(77,file='tail13.dat',status='unknown',form='formatted')
    open(78,file='tail14.dat',status='unknown',form='formatted')
    open(79,file='tail15.dat',status='unknown',form='formatted')
    open(80,file='tail16.dat',status='unknown',form='formatted')
    if(spacecraft) then
        open(41,file='wind.pos',status='unknown',form='formatted')
        !       open(42,file='polar.pos',status='unknown',form='formatted')
        !       open(43,file='equators.pos',status=unknown'',form='formatted')
        !       open(44,file='geotail.pos',status='unknown',form='formatted')
        !       open(47,file='wind.den',status='unknown',form='formatted')
        !       open(48,file='wind.vel',status='unknown',form='formatted')
        open(47,file='wind.plas',status='unknown',form='formatted')
        open(49,file='wind.mag',status='unknown',form='formatted')
    endif
    !
    !     the magnetic field in dimensionless units is actually in alfven speeds
    !             relative to the normalizing velocity
    !     the temperature is in sound speeds relative to the normalizing speed
    !             all squared
    !
    !     the magnetospheric plasma density are assumed to vary as
    !             rho proportional to (r)**-alpha_e
    !             temperatue proportional to (r)**alpha_e
    !         with total pressure constant which is needed for equilibrium
    !
    !     now find the equivalent pressure of magnetosphere for the given
    !         sound speed
    !
    gamma1=gamma-1.
    epress=cs_inner**2*erho/gamma
    eerg=epress/gamma1
    !
    !     do the same for the solar wind
    !
    rho_wind1=den_wind1*rmassq
    rho_wind2=den_wind2*rmassq
    srho=rho_wind1
    delrho=(rho_wind2-rho_wind1)/tmax
    spress=(cs_wind**2*srho/gamma)/gamma1
    svelx=vx_wind1
    svely=vy_wind1
    svelz=vz_wind1
    !
    spx=srho*svelx
    spy=srho*svely
    spz=srho*svelz
    serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
    delvx_wind=(vx_wind2-vx_wind1)/tmax
    delvy_wind=(vy_wind2-vy_wind1)/tmax
    delvz_wind=(vz_wind2-vz_wind1)/tmax
    !
    delbx_wind=(alfx_wind2*sqrt(rho_wind2) &
    -alfx_wind1*sqrt(rho_wind1))/tmax
    delby_wind=(alfy_wind2*sqrt(rho_wind2) &
    -alfy_wind1*sqrt(rho_wind1))/tmax
    delbz_wind=(alfz_wind2*sqrt(rho_wind2) &
    -alfz_wind1*sqrt(rho_wind1))/tmax
    sbx_wind=alfx_wind1*sqrt(rho_wind1)
    sby_wind=alfy_wind1*sqrt(rho_wind1)
    sbz_wind=alfz_wind1*sqrt(rho_wind1)
    !
    deltg=tmax/float(ntgraf)
    delt=stepsz
    tgraf=0.
    ts1=tsave
    tdiv=10.
    del_tdiv=5.
    write_dat=.true.
    !
    !     ************************************************
    !     check for restart
    !     ************************************************
    !
    nchf=11
    if(.not.start) then
        !
        !     main grid restart files
        !
        if(reload)then
            nchf=15
            open(nchf,file='fluid00',status='unknown',form='unformatted')
        else
            open(11,file='fluid11',status='unknown',form='unformatted')
            open(12,file='fluid12',status='unknown',form='unformatted')
            read(nchf) t
            rewind nchf
            inchf=23-nchf
            read(inchf) t2
            rewind inchf
            !
            if(t.lt.t2) then
                close(nchf)
                nchf=inchf
            else
                close(inchf)
            endif
        endif
        !
        !     read restart data
        !
        read(nchf)t
        read(nchf)qrho
        read(nchf)qpx
        read(nchf)qpy
        read(nchf)qpz
        read(nchf)qpresx
        if(isotropic)then
            qpresy=qpresx
            qpresz=qpresx
            qpresxy=0.
            qpresxz=0.
            qpresyz=0.
        else
            write(6,*)'anistropic pressure read'
            read(nchf)qpresy
            read(nchf)qpresz
            read(nchf)qpresxy
            read(nchf)qpresxz
            read(nchf)qpresyz
        endif
        read(nchf)hrho
        read(nchf)hpx
        read(nchf)hpy
        read(nchf)hpz
        read(nchf)hpresx
        if(isotropic)then
            hpresy=hpresx
            hpresz=hpresx
            hpresxy=0.
            hpresxz=0.
            hpresyz=0.
        else
            write(6,*)'anistropic pressure read'
            read(nchf)hpresy
            read(nchf)hpresz
            read(nchf)hpresxy
            read(nchf)hpresxz
            read(nchf)hpresyz
        endif
        read(nchf)orho
        read(nchf)opx
        read(nchf)opy
        read(nchf)opz
        read(nchf)opresx
        if(isotropic)then
            opresy=opresx
            opresz=opresx
            opresxy=0.
            opresxz=0.
            opresyz=0.
        else
            write(6,*)'anistropic pressure read'
            read(nchf)opresy
            read(nchf)opresz
            read(nchf)opresxy
            read(nchf)opresxz
            read(nchf)opresyz
        endif
        read(nchf)bx
        read(nchf)by
        read(nchf)bz
        read(nchf)epres
        read(nchf)bx0
        read(nchf)by0
        read(nchf)bz0
        !      read(nchf)qrho0
        !      read(nchf)hrho0
        !      read(nchf)orho0
        !      read(nchf)qpres0
        !      read(nchf)hpres0
        !      read(nchf)opres0
        !      read(nchf)epres0
        read(nchf)parm_srf,parm_mid,parm_zero, &
        ijzero,numzero,ijmid,nummid,ijsrf,numsrf
        close(nchf)
        !
        !
        !      check for div b errors
        !
        if(divb_lores)then
            range=1.33*lunar_dist/re_equiv
            write(6,*)'range for divb taper',lunar_dist,range
            do m=ngrd,1,-1
                write(6,*)'divb on box',m
                !
                call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                    curx,cury,curz,efldx,efldy,efldz, &
                    7,nx*ny*nz,nx,ny,nz,ngrd,m,xspac)
                !
                call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                    curx,cury,curz,efldx,efldy,efldz, &
                    7,nx*ny*nz,nx,ny,nz,ngrd,m,xspac)
                write(6,*)'completed divb on box',m
                !
                !
            enddo !end m loop
            !
            !        apply bndry conditions
            !
            do m=ngrd-1,1,-1
                call flanks_synced(bx,nx,ny,nz,ngrd,m, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                call flanks_synced(by,nx,ny,nz,ngrd,m, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                call flanks_synced(bz,nx,ny,nz,ngrd,m, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            enddo
            do m=1,ngrd-1
                call bndry_corer( &
                    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                    orho,opresx,opresy,opresz,opx,opy,opz, &
                    epres,bx,by,bz, &
                    qpresxy,qpresxz,qpresyz, &
                    hpresxy,hpresxz,hpresyz, &
                    opresxy,opresxz,opresyz, &
                    nx,ny,nz,ngrd,m, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
            enddo  !end bndry_corer
        endif  ! end divb_lores
        !
        !     if(update)t=0.
        write(6,*)'entering lores visual'
        ut=utstart+t*t_equiv/3600.
        call visual(qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
            orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
            epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
            curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
            tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
            nx,ny,nz,ngrd,xspac, &
            cross,along,flat,xcraft,ncraft,re_equiv, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
            grd_zmin,grd_zmax,ut,b_equiv,ti_te,rho_equiv)

        !
        !      read fine grid if .not.reload
        !
        if(.not.reload)then
            mchf=nchf+10
            if(mchf.eq.21)then
                open(21,file='fluid21',status='unknown',form='unformatted')
            else
                open(22,file='fluid22',status='unknown',form='unformatted')
            endif
            write(6,*)'reading refined grid'
            !
            !      read restart data
            !
            read(mchf)t_n,ut_insert
            if(t.ne.t_n)then
                write(6,*)'time starts don"t match',t,t_n,nchf,mchf
                stop
            endif
            read(mchf)qrho_n
            read(mchf)qpx_n
            read(mchf)qpy_n
            read(mchf)qpz_n
            read(mchf)qpresx_n
            if(isotropic)then
                qpresy_n=qpresx_n
                qpresz_n=qpresx_n
                qpresxy_n=0.
                qpresxz_n=0.
                qpresyz_n=0.
            else
                write(6,*)'anistropic pressure read'
                read(mchf)qpresy_n
                read(mchf)qpresz_n
                read(mchf)qpresxy_n
                read(mchf)qpresxz_n
                read(mchf)qpresyz_n
            endif
            read(mchf)hrho_n
            read(mchf)hpx_n
            read(mchf)hpy_n
            read(mchf)hpz_n
            read(mchf)hpresx_n
            if(isotropic)then
                hpresy_n=hpresx_n
                hpresz_n=hpresx_n
                hpresxy_n=0.
                hpresxz_n=0.
                hpresyz_n=0.
            else
                write(6,*)'anistropic pressure read'
                read(mchf)hpresy_n
                read(mchf)hpresz_n
                read(mchf)hpresxy_n
                read(mchf)hpresxz_n
                read(mchf)hpresyz_n
            endif
            read(mchf)orho_n
            read(mchf)opx_n
            read(mchf)opy_n
            read(mchf)opz_n
            read(mchf)opresx_n
            if(isotropic)then
                opresy_n=opresx_n
                opresz_n=opresx_n
                opresxy_n=0.
                opresxz_n=0.
                opresyz_n=0.
            else
                write(6,*)'anistropic pressure read'
                read(mchf)opresy_n
                read(mchf)opresz_n
                read(mchf)opresxy_n
                read(mchf)opresxz_n
                read(mchf)opresyz_n
            endif
            read(mchf)bx_n
            read(mchf)by_n
            read(mchf)bz_n
            read(mchf)epres_n
            read(mchf)bx0_n
            read(mchf)by0_n
            read(mchf)bz0_n
            read(mchf)parm_srf_n,parm_mid_n,parm_zero_n, &
                ijzero_n,numzero_n,ijmid_n,nummid_n,ijsrf_n,numsrf_n
            read(mchf)grd_xmin_n,grd_xmax_n, &
                grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n, &
                grd_time_n,grd_dx_n,grd_dy_n,grd_vx_n,grd_vy_n
            close(mchf)
            !
    
            !
            !      initialize position of the moon
            !
            !      recalibrate ut_insert  due to wrong initial orbital period 85.5
            !
            if(update)then
                write(6,*)'modifying ut_inset at ut=',ut
                theta=((ut-ut_insert)/85.5)*2.*3.1414
                write(6,*)'ut_insert was',ut_insert,theta
                !
                ut_insert=ut-theta*ut_orbit/(2.*3.1414)
                theta=((ut-ut_insert)/ut_orbit)*2.*3.1414
                write(6,*)'ut_insert changed to',ut_insert,theta
            endif
            !
    
            tempi=cs_moon**2/gamma    ! units in proton masses
            ut=utstart+t*t_equiv/3600.
            theta=((ut-ut_insert)/ut_orbit)*2.*3.1414
            xmoon=r_orbit*sin(theta)
            ymoon=-r_orbit*cos(theta)
            zmoon=0.
            vx_moon=v_orbit*cos(theta)
            vy_moon=v_orbit*sin(theta)
            vz_moon=0.
            write(6,*)'restart moon at'
            write(6,*)ut,ut_insert,xmoon,ymoon,zmoon
        else
            !
            !      initialize quantities on refined grid from main grid
            !
            ut=utstart+t*t_equiv/3600.
    
            ut_insert=ut-(theta_moon/360.) *ut_orbit   !  save time of refinement
    
            theta=((ut-ut_insert)/ut_orbit)*2.*3.1414
            xmoon=r_orbit*sin(theta)
            ymoon=-r_orbit*cos(theta)
            zmoon=0.
            vx_moon=v_orbit*cos(theta)
            vy_moon=v_orbit*sin(theta)
            vz_moon=0.
            write(6,*)'refinement moon at',ut,theta
            !
            !         species 1
            !
    
            write(6,*)'starting refinement'
            call refinement(qrho,nx,ny,nz,ngrd,main_grd_n, &
                qrho_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpx,nx,ny,nz,ngrd,main_grd_n, &
                qpx_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpy,nx,ny,nz,ngrd,main_grd_n, &
                qpy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpz,nx,ny,nz,ngrd,main_grd_n, &
                qpz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpresx,nx,ny,nz,ngrd,main_grd_n, &
                qpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpresy,nx,ny,nz,ngrd,main_grd_n, &
                qpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpresz,nx,ny,nz,ngrd,main_grd_n, &
                qpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpresxy,nx,ny,nz,ngrd,main_grd_n, &
                qpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpresxz,nx,ny,nz,ngrd,main_grd_n, &
                qpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(qpresyz,nx,ny,nz,ngrd,main_grd_n, &
                qpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            !
            !               species 2
            !
            call refinement(hrho,nx,ny,nz,ngrd,main_grd_n, &
                hrho_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpx,nx,ny,nz,ngrd,main_grd_n, &
                hpx_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpy,nx,ny,nz,ngrd,main_grd_n, &
                hpy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpz,nx,ny,nz,ngrd,main_grd_n, &
                hpz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpresx,nx,ny,nz,ngrd,main_grd_n, &
                hpresx_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpresy,nx,ny,nz,ngrd,main_grd_n, &
                hpresy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpresz,nx,ny,nz,ngrd,main_grd_n, &
                hpresz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpresxy,nx,ny,nz,ngrd,main_grd_n, &
                hpresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpresxz,nx,ny,nz,ngrd,main_grd_n, &
                hpresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(hpresyz,nx,ny,nz,ngrd,main_grd_n, &
                hpresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            !
            !              species 3
            !
            call refinement(orho,nx,ny,nz,ngrd,main_grd_n, &
                orho_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opx,nx,ny,nz,ngrd,main_grd_n, &
                opx_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opy,nx,ny,nz,ngrd,main_grd_n, &
                opy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opz,nx,ny,nz,ngrd,main_grd_n, &
                opz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opresx,nx,ny,nz,ngrd,main_grd_n, &
                opresx_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opresy,nx,ny,nz,ngrd,main_grd_n, &
                opresy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opresz,nx,ny,nz,ngrd,main_grd_n, &
                opresz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opresxy,nx,ny,nz,ngrd,main_grd_n, &
                opresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opresxz,nx,ny,nz,ngrd,main_grd_n, &
                opresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(opresyz,nx,ny,nz,ngrd,main_grd_n, &
                opresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            !
            !              electrons
            !
            call refinement(epres,nx,ny,nz,ngrd,main_grd_n, &
                epres_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            !
            !              magnetic field
            !
            call refinement(bx,nx,ny,nz,ngrd,main_grd_n, &
                bx_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(by,nx,ny,nz,ngrd,main_grd_n, &
                by_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(bz,nx,ny,nz,ngrd,main_grd_n, &
                bz_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            !
            !
            !       set background magnetic field
            !
            !
            call refinement(bx0,nx,ny,nz,ngrd,main_grd_n, &
                bx0_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(by0,nx,ny,nz,ngrd,main_grd_n, &
                by0_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            call refinement(bz0,nx,ny,nz,ngrd,main_grd_n, &
                bz0_n,nx_n,ny_n,nz_n,ngrd_n, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n)
            !
            !
            !      determine any regions in contact with inner bndry
            !
            write(6,*)'moon grd pts',xmoon,ymoon,zmoon,rmoon
            !
            tempi=cs_moon**2/gamma    ! units in proton masses
            alpha_m=alpha_e
            range_moon=12.*rmoon       ! distance from the moon where dens applied
            !
            !      determine position and speed of the moon
            !
            theta=((ut-ut_insert)/ut_orbit)*2.*3.1414
            xmoon=r_orbit*sin(theta)
            ymoon=-r_orbit*cos(theta)
            zmoon=0.
            vx_moon=v_orbit*cos(theta)
            vy_moon=v_orbit*sin(theta)
            vz_moon=0.
            !
            tempi=cs_moon**2/gamma    ! units in proton masses
            alpha_m=alpha_e
            range_moon=12.*rmoon       ! distance from the moon where dens applied
            !
            !      determine position and speed of the moon
            !
            theta=((ut-ut_insert)/ut_orbit)*2.*3.1414
            xmoon=r_orbit*sin(theta)
            ymoon=-r_orbit*cos(theta)
            zmoon=0.
            vx_moon=v_orbit*cos(theta)
            vy_moon=v_orbit*sin(theta)
            vz_moon=0.
            !
            !       add density of the moon
            !
            do m=1,ngrd_n
                !
                !        initialize indices
                !
                if (m.le.mbndry_n) then
                    numsrf_n(m)=0
                    nummid_n(m)=0
                    numzero_n(m)=0
                endif
                !
                !        initialize plasma
                !
                dx=(grd_xmax_n(m)-grd_xmin_n(m))/(nx_n-1.)
                dy=(grd_ymax_n(m)-grd_ymin_n(m))/(ny_n-1.)
                dz=(grd_zmax_n(m)-grd_zmin_n(m))/(nz_n-1.)
                !         write(6,*)'grid spacing at',m,grd_xmin_n(m),grd_xmax_n(m),
                !    +        grd_ymin_n(m),grd_ymax_n(m),grd_zmin_n(m),grd_zmax_n(m)
                !         write(6,*)'rmoon',rmoon,dx
                !
                !
                !        speed for ram induced ionosphere
                !
                vtot=sqrt(vx_moon**2+vy_moon**2)+1.e-11
                !
                do k=1,nz_n
                    az=grd_zmin_n(m)+dz*(k-1)
                    zm=az-zmoon
                    do j=1,ny_n
                        ay=grd_ymin_n(m)+dy*(j-1)
                        ym=ay-ymoon
                        do i=1,nx_n
                            ax=grd_xmin_n(m)+dx*(i-1)
                            xm=ax-xmoon
                            ar_moon=sqrt(xm**2+ym**2+zm**2)
                            if(ar_moon.le.range_moon)then
                                !
                                ra_moon=((ar_moon+0.3*rmoon)/(1.3*rmoon))**(-alpha_m)
                                ra_moon=amin1(1.,ra_moon)
                                !
                                !               allow for ram density profile
                                !
                                !
                                rm=(vx_moon*xm+vy_moon*ym)/vtot
                                !
                                if(rm.lt.0.0)then
                                    xscale=1.
                                else
                                    xr=sqrt(rm**2+(offset*rmoon)**2)-(offset*rmoon)
                                    xscale=exp(-xr/(offset*rmoon))
                                endif
                                !
                                qden=qden_moon*ra_moon*xscale
                                hden=hden_moon*ra_moon*xscale
                                oden=oden_moon*ra_moon*xscale
                                !
                                !               qden=amax1(qrho_n(i,j,k,m)/rmassq,qden)
                                !               hden=amax1(hrho_n(i,j,k,m)/rmassh,hden)
                                !               oden=amax1(orho_n(i,j,k,m)/rmasso,oden)
                                !
                                qrho_n(i,j,k,m)=qrho_n(i,j,k,m)+qden*rmassq
                                hrho_n(i,j,k,m)=hrho_n(i,j,k,m)+hden*rmassh
                                orho_n(i,j,k,m)=orho_n(i,j,k,m)+oden*rmasso
                                !
                                qpx_n(i,j,k,m)=qpx_n(i,j,k,m)+qden*rmassq*vx_moon
                                qpy_n(i,j,k,m)=qpy_n(i,j,k,m)+qden*rmassq*vy_moon
                                !
                                hpx_n(i,j,k,m)=hpx_n(i,j,k,m)+hden*rmassh*vx_moon
                                hpy_n(i,j,k,m)=hpy_n(i,j,k,m)+hden*rmassh*vy_moon
                                !
                                opx_n(i,j,k,m)=opx_n(i,j,k,m)+oden*rmasso*vx_moon
                                opy_n(i,j,k,m)=opy_n(i,j,k,m)+oden*rmasso*vy_moon
                                !
                                qpresx_n(i,j,k,m)=qpresx_n(i,j,k,m)+tempi*qden
                                qpresy_n(i,j,k,m)=qpresy_n(i,j,k,m)+tempi*qden
                                qpresz_n(i,j,k,m)=qpresz_n(i,j,k,m)+tempi*qden
                                !
                                hpresx_n(i,j,k,m)=hpresx_n(i,j,k,m)+tempi*hden
                                hpresy_n(i,j,k,m)=hpresy_n(i,j,k,m)+tempi*hden
                                hpresz_n(i,j,k,m)=hpresz_n(i,j,k,m)+tempi*hden*0.9
                                !
                                opresx_n(i,j,k,m)=opresx_n(i,j,k,m)+tempi*oden
                                opresy_n(i,j,k,m)=opresy_n(i,j,k,m)+tempi*oden
                                opresz_n(i,j,k,m)=opresz_n(i,j,k,m)+tempi*oden
                                !
                                epres_n(i,j,k,m)=epres_n(i,j,k,m)+ &
                                (tempi*(qden+hden+oden)) &
                                /ti_te_moon
                            endif
    
                        enddo
                    enddo
                enddo
            enddo  !end moon plasma conditions
            !
            do m_n=1,ngrd_n
                call set_rho(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                qpresxy_n,qpresxz_n,qpresyz_n,rmassq, &
                hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                hpresxy_n,hpresxz_n,hpresyz_n,rmassh, &
                orho_n,opresx_n,opresy_n,opresz_n, &
                opresxy_n,opresxz_n,opresyz_n,rmasso, &
                epres_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                o_conc)
            enddo
            !
            !
            !      paramters for moving grid  : all in sim units
            !
            do m=1,ngrd_n
                grd_time_n(m)=t
                grd_dx_n(m)=0.
                grd_dy_n(m)=0.
                grd_vx_n(m)=vx_moon
                grd_vy_n(m)=vy_moon
            enddo
            !
            write(7,*)'time',ut,'moon pos',xmoon,ymoon
            do m_n=1,ngrd_n
                write(7,*)m_n,grd_xmin_n(m_n),grd_xmax_n(m_n), &
                grd_ymin_n(m_n),grd_ymax_n(m_n)
            enddo
            !
            !
            write(6,*)'moon in grid points at',xmoon,ymoon,zmoon
            write(6,*)'moon in planetary r at',xmoon*re_equiv, &
            ymoon*re_equiv, &
            zmoon*re_equiv
            !
            !
        endif   !if(.not.reload) (restart/refinement)

        !
        !      apply boundary conditions to test if code is working and replot
        !
        !    first set inner boundary conditions
        !
        write(6,*)'calling bndry_moon'
        do m_n=1,mbndry_n
            !
            call set_bndry_moon_ram(rmassq,rmassh,rmasso, &
            m_n,nx_n,ny_n,nz_n,ngrd_n, &
            parm_srf_n,parm_mid_n,parm_zero_n, &
            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
            qden_moon,hden_moon,oden_moon, &
            vx_moon,vy_moon,tempi,gamma, &
            ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_e, &
            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
            grd_zmin_n,grd_zmax_n)
            !
            call bndry_moon(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
            qpresxy_n,qpresxz_n,qpresyz_n, &
            qpx_n,qpy_n,qpz_n,rmassq, &
            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
            hpresxy_n,hpresxz_n,hpresyz_n, &
            hpx_n,hpy_n,hpz_n,rmassh, &
            orho_n,opresx_n,opresy_n,opresz_n, &
            opresxy_n,opresxz_n,opresyz_n, &
            opx_n,opy_n,opz_n,rmasso, &
            epres_n,bx_n,by_n,bz_n,m_n, &
            nx_n,ny_n,nz_n,ngrd_n, &
            parm_srf_n,parm_mid_n,parm_zero_n, &
            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
            qden_moon,hden_moon,oden_moon,cs_moon,gamma, &
            ti_te_moon,vx_moon,vy_moon,vz_moon, &
            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
            grd_zmin_n,grd_zmax_n)
        enddo
        !
        write(6,*)mbndry_n,numsrf_n,nummid_n,numzero_n
        !
        do m_n=1,ngrd_n
            call set_rho(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
            qpresxy_n,qpresxz_n,qpresyz_n,rmassq, &
            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
            hpresxy_n,hpresxz_n,hpresyz_n,rmassh, &
            orho_n,opresx_n,opresy_n,opresz_n, &
            opresxy_n,opresxz_n,opresyz_n,rmasso, &
            epres_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
            o_conc)
            !
            vlim=0.6
            call set_speed_agrd( &
            qrho_n,qpresx_n,qpresy_n,qpresz_n,qpx_n,qpy_n,qpz_n, &
            hrho_n,hpresx_n,hpresy_n,hpresz_n,hpx_n,hpy_n,hpz_n, &
            orho_n,opresx_n,opresy_n,opresz_n,opx_n,opy_n,opz_n, &
            epres_n,qpresxy_n,qpresxz_n,qpresyz_n, &
            hpresxy_n,hpresxz_n,hpresyz_n, &
            opresxy_n,opresxz_n,opresyz_n, &
            bx_n,by_n,bz_n,bx0_n,by0_n,bz0_n, &
            bsx_n,bsy_n,bsz_n,btot_n, &
            vvx_n,tvx_n,tvy_n,tvz_n,evx_n,evy_n,evz_n, &
            curx_n,cury_n,curz_n, &
            rmassq,rmassh,rmasso,nx_n,ny_n,nz_n,ngrd_n,m_n, &
            pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
            vlim,alf_lim,o_conc,fastest)
        enddo
        !
        !      apply bndry conditions to other grids
        !
        do m_n=1,ngrd_n-1
            write(6,*)'bndry-corer m_n',m_n
            call bndry_corer(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
            qpx_n,qpy_n,qpz_n, &
            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
            hpx_n,hpy_n,hpz_n, &
            orho_n,opresx_n,opresy_n,opresz_n, &
            opx_n,opy_n,opz_n, &
            epres_n,bx_n,by_n,bz_n, &
            qpresxy_n,qpresxz_n,qpresyz_n, &
            hpresxy_n,hpresxz_n,hpresyz_n, &
            opresxy_n,opresxz_n,opresyz_n, &
            nx_n,ny_n,nz_n,ngrd_n,m_n, &
            grd_xmin_n,grd_xmax_n, &
            grd_ymin_n,grd_ymax_n,grd_zmin_n,grd_zmax_n)
        enddo
        !
        write(6,*)'calling hires visual'
    
        call visual_hires(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
        qpx_n,qpy_n,qpz_n,rmassq, &
        hrho_n,hpresx_n,hpresy_n,hpresz_n, &
        hpx_n,hpy_n,hpz_n,rmassh, &
        orho_n,opresx_n,opresy_n,opresz_n, &
        opx_n,opy_n,opz_n,rmasso, &
        epres_n,bx_n,by_n,bz_n,bx0_n,by0_n,bz0_n, &
        bsx_n,bsy_n,bsz_n, &
        curx_n,cury_n,curz_n,efldx_n,efldy_n,efldz_n, &
        tvx_n,tvy_n,tvz_n, &
        tx_n,ty_n,tz_n,tg1_n,tg2_n,tt_n,work_n, &
        mx_n,my_n,mz_n,mz2_n,muvwp2_n, &
        nx_n,ny_n,nz_n,ngrd_n,xspac_n, &
        cross_n,along_n,flat_n,xcraft,ncraft,re_equiv, &
        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
        grd_zmin_n,grd_zmax_n,ut,b_equiv, &
        ti_te_moon,rho_equiv)
        !
        ts1=t+tsave
        tstep=tmax
        tmax=t+tmax
        tgraf=t+deltg
        tdiv=t
        nchf=11
        ut=utstart+t*t_equiv/3600.
        !
        !     rescale oxygen density at inner boundary by oxygen scale
        !
        if(ringo)then
            !
            !          check for io plasma torus addtions
            !
            m=1
            dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
            dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
            dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
            !
            tot_o=0.
            tot_h=0.
            tot_q=0.
            !
            write(6,*)'ringo',lunar_dist,reduct,v_rot
            !
            do k=1,nz
                az=(grd_zmin(m)+dz*(k-1)-zdip)
                do j=1,ny
                    ay=grd_ymin(m)+dy*(j-1)-ydip
                    do i=1,nx
                        ax=(grd_xmin(m)+dx*(i-1)-xdip)
                        !
                        rx=ax*re_equiv
                        ry=ay*re_equiv
                        rd=sqrt(rx**2+ry**2)
                        !
                        rvy=rx*v_rot
                        rvx=-ry*v_rot
                        rvz=0.
                        corotate=sqrt(rvx**2+rvy**2)
                        !
                        !
                        !
                        !          ar=sqrt(ax**2+ay**2)
                        !          if((ar.lt.1.66*rearth).and.
                        !    +        (abs(az).le.0.5*rearth) )then
                        !           qpx(i,j,k,m)=rvx*qrho(i,j,k,m)
                        !           qpy(i,j,k,m)=rvy*qrho(i,j,k,m)
                        !           hpx(i,j,k,m)=rvx*hrho(i,j,k,m)
                        !           hpy(i,j,k,m)=rvy*hrho(i,j,k,m)
                        !           opx(i,j,k,m)=rvx*orho(i,j,k,m)
                        !           opy(i,j,k,m)=rvy*orho(i,j,k,m)
                        !          endif
                        !
                        !
                        ar_io=sqrt(ax**2+ay**2)*re_equiv
                        dr_io=abs(ar_io -lunar_dist)
                        !
                        !          check for torus
                        !
                        if(abs(dr_io.lt.2.*torus_rad)) then
                            rscale=exp(-((dr_io)/(0.5*torus_rad))**2)  ! scale height in rjs
                            zscale=exp(-((az*re_equiv)/(0.5*torus_rad))**2)
                            !
                            oden=den_lunar*rmasso*rscale*zscale
                            !            write(6,*)i,j,k,m,oden
                            orho(i,j,k,m)=orho(i,j,k,m) +oden
                            del_op=(oden/rmasso)*(corotate**2)*t_torus!temp goes as v**2
                            opresx(i,j,k,m)=opresx(i,j,k,m)+del_op
                            opresy(i,j,k,m)=opresy(i,j,k,m)+del_op
                            opresz(i,j,k,m)=opresz(i,j,k,m)+del_op*aniso_factor
                            opx(i,j,k,m)=opx(i,j,k,m)+reduct*oden*rvx
                            opy(i,j,k,m)=opy(i,j,k,m)+reduct*oden*rvy
                            !
                            hden=oden*rmassh/rmasso/o_conc
                            hrho(i,j,k,m)=hrho(i,j,k,m) +hden
                            del_hp=(hden/rmassh)*(corotate**2)*t_torus  !temp goes as v**2
                            hpresx(i,j,k,m)=hpresx(i,j,k,m)+del_hp
                            hpresy(i,j,k,m)=hpresy(i,j,k,m)+del_hp
                            hpresz(i,j,k,m)=hpresz(i,j,k,m)+del_hp*aniso_factor
                            hpx(i,j,k,m)=hpx(i,j,k,m)+reduct*hden*rvx
                            hpy(i,j,k,m)=hpy(i,j,k,m)+reduct*hden*rvy
                            !
                            qden=0.5*hden*rmassq/rmassh
                            qrho(i,j,k,m)=qrho(i,j,k,m) +qden
                            del_qp=(qden/rmassq)*(corotate**2)*t_torus    !temp goes as v**2
                            qpresx(i,j,k,m)=qpresx(i,j,k,m)+del_qp
                            qpresy(i,j,k,m)=qpresy(i,j,k,m)+del_qp
                            qpresz(i,j,k,m)=qpresz(i,j,k,m)+del_qp*aniso_factor
                            qpx(i,j,k,m)=qpx(i,j,k,m)+reduct*qden*rvx
                            qpy(i,j,k,m)=qpy(i,j,k,m)+reduct*qden*rvy
                            !
                            epres(i,j,k,m)=epres(i,j,k,m)+(del_op+del_hp+del_qp) ! equal temps
                            !
                            tot_q=tot_q+qden*dx*dy*dz
                            tot_h=tot_h+hden*dx*dy*dz
                            tot_o=tot_o+oden*dx*dy*dz
                            !
                        endif
                    enddo
                enddo
            enddo
            !
            !     scale factors to kg/s
            !
            volume=(re_equiv*planet_rad*1.e3)**3  !(cubic meters)
            atime=tsave*t_equiv
            injections=tstep/tsave
            proton_mass=1.67e-27
            write(6,*)'volume,t_equiv,atime',volume,t_equiv,atime
            !
            tot_q=tot_q*volume/atime*rho_equiv*1.e6/rmassq
            tot_h=tot_h*volume/atime*rho_equiv*1.e6/rmassh
            tot_o=tot_o*volume/atime*rho_equiv*1.e6/rmasso
            write(6,*)'tot torus ions/s',tot_q,tot_h,tot_o
            write(10,*)'tot torus ions/s',tot_q,tot_h,tot_o
            !
            tot_q=tot_q*proton_mass*rmassq
            tot_h=tot_h*proton_mass*rmassh
            tot_o=tot_o*proton_mass*rmasso
            write(6,*)'injections',injections, '  single at'
            write(6,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
            write(10,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
            !
            !
            !          add into high res grid
            !
            do m=ngrd_n,1,-1
                dx=(grd_xmax_n(m)-grd_xmin_n(m))/(nx_n-1.)
                dy=(grd_ymax_n(m)-grd_ymin_n(m))/(ny_n-1.)
                dz=(grd_zmax_n(m)-grd_zmin_n(m))/(nz_n-1.)
                !
                tot_o=0.
                tot_h=0.
                tot_q=0.
                !
                do k=1,nz_n
                    az=grd_zmin_n(m)+dz*(k-1)
                    do j=1,ny_n
                        ay=grd_ymin_n(m)+dy*(j-1)
                        do i=1,nx_n
                            ax=grd_xmin_n(m)+dx*(i-1)
                            !
                            rx=ax*re_equiv
                            ry=ay*re_equiv
                            rd=sqrt(rx**2+ry**2)
                            !
                            rvy=rx*v_rot
                            rvx=-ry*v_rot
                            rvz=0.
                            corotate=sqrt(rvx**2+rvy**2)
                            !
                            !          enforce co-rotation
                            !
                            !          ar=sqrt(ax**2+ay**2)
                            !          if((ar.lt.1.66*rearth).and.
                            !    +        (abs(az).le.0.5*rearth) )then
                            !           qpx_n(i,j,k,m)=rvx*qrho_n(i,j,k,m)
                            !           qpy_n(i,j,k,m)=rvy*qrho_n(i,j,k,m)
                            !           hpx_n(i,j,k,m)=rvx*hrho_n(i,j,k,m)
                            !           hpy_n(i,j,k,m)=rvy*hrho_n(i,j,k,m)
                            !           opx_n(i,j,k,m)=rvx*orho_n(i,j,k,m)
                            !           opy_n(i,j,k,m)=rvy*orho_n(i,j,k,m)
                            !          endif
                            !
                            ar_io=sqrt(ax**2+ay**2)*re_equiv
                            dr_io=abs(ar_io -lunar_dist)
                            !
                            if(abs(dr_io.lt.2.*torus_rad)) then
                                rscale=exp(-((dr_io)/(0.5*torus_rad))**2) ! scale height in rjs
                                zscale=exp(-((az*re_equiv)/(0.5*torus_rad))**2)
                                !
                                oden=den_lunar*rmasso*rscale*zscale
                                !            write(6,*)i,j,k,m,oden
                                orho_n(i,j,k,m)=orho_n(i,j,k,m) +oden
                                del_op=(oden/rmasso)*(corotate**2)*t_torus !temp goes as v**2
                                opresx_n(i,j,k,m)=opresx_n(i,j,k,m)+del_op
                                opresy_n(i,j,k,m)=opresy_n(i,j,k,m)+del_op
                                opresz_n(i,j,k,m)=opresz_n(i,j,k,m)+del_op*aniso_factor
                                opx_n(i,j,k,m)=opx_n(i,j,k,m)+reduct*oden*rvx
                                opy_n(i,j,k,m)=opy_n(i,j,k,m)+reduct*oden*rvy
                                !
                                hden=oden*rmassh/rmasso/o_conc
                                hrho_n(i,j,k,m)=hrho_n(i,j,k,m) +hden
                                del_hp=(hden/rmassh)*(corotate**2)*t_torus  !temp goes as v**2
                                hpresx_n(i,j,k,m)=hpresx_n(i,j,k,m)+del_hp
                                hpresy_n(i,j,k,m)=hpresy_n(i,j,k,m)+del_hp
                                hpresz_n(i,j,k,m)=hpresz_n(i,j,k,m)+del_hp*aniso_factor
                                hpx_n(i,j,k,m)=hpx_n(i,j,k,m)+reduct*hden*rvx
                                hpy_n(i,j,k,m)=hpy_n(i,j,k,m)+reduct*hden*rvy
                                !
                                qden=0.5*hden*rmassq/rmassh
                                qrho_n(i,j,k,m)=qrho_n(i,j,k,m) +qden
                                del_qp=(qden/rmassq)*(corotate**2) *t_torus  !temp goes as v**2
                                qpresx_n(i,j,k,m)=qpresx_n(i,j,k,m)+del_qp
                                qpresy_n(i,j,k,m)=qpresy_n(i,j,k,m)+del_qp
                                qpresz_n(i,j,k,m)=qpresz_n(i,j,k,m)+del_qp*aniso_factor
                                qpx_n(i,j,k,m)=qpx_n(i,j,k,m)+reduct*qden*rvx
                                qpy_n(i,j,k,m)=qpy_n(i,j,k,m)+reduct*qden*rvy
                                !
                                epres_n(i,j,k,m)=epres_n(i,j,k,m) &
                                +(del_op+del_hp+del_qp) ! equal temps
                                !
                                tot_q=tot_q+qden*dx*dy*dz
                                tot_h=tot_h+hden*dx*dy*dz
                                tot_o=tot_o+oden*dx*dy*dz
                                !
                            endif
                        enddo
                    enddo
                enddo
                !
                !     scale factors to kg/s
                !
                volume=(re_equiv*planet_rad*1.e3)**3  !(cubic meters)
                atime=tsave*t_equiv
                proton_mass=1.67e-27
                !      write(6,*)'volume,t_equiv,atime',volume,t_equiv,atime
                !
                tot_q=tot_q*volume/atime*rho_equiv*1.e6/rmassq
                tot_h=tot_h*volume/atime*rho_equiv*1.e6/rmassh
                tot_o=tot_o*volume/atime*rho_equiv*1.e6/rmasso
                write(6,*)'hires tot ions/s',m,tot_q,tot_h,tot_o
                !      write(10,*)'hires tot ions/s',m,tot_q,tot_h,tot_o
                !
                tot_q=tot_q*proton_mass*rmassq
                tot_h=tot_h*proton_mass*rmassh
                tot_o=tot_o*proton_mass*rmasso
                write(6,*)'hires tot  kg/s',m,tot_q,tot_h,tot_o
                !      write(10,*)'hires tot  kg/s',m,tot_q,tot_h,tot_o
            enddo
            !
        endif  ! end ringo if
        !
        !     initialize plasma resistivity
        !
        call set_resist(resistive,nx,ny,nz,mbndry,resist, &
            ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
            msrf,mmid,mzero,1.)
        !
        call set_resist(resistive_n,nx_n,ny_n,nz_n,mbndry_n,resist, &
            ijzero_n,numzero_n,ijmid_n,nummid_n,ijsrf_n, &
            numsrf_n,msrf_n,mmid_n,mzero_n,0.01)
    
    
        write(6,*)'entering lores visual'
        ut=utstart+t*t_equiv/3600.
        call visual(qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
            orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
            epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
            curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
            tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
            nx,ny,nz,ngrd,xspac, &
            cross,along,flat,xcraft,ncraft,re_equiv, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
            grd_zmin,grd_zmax,ut,b_equiv,ti_te,rho_equiv)
        !
        write(6,*)'calling hires visual'
        !
        call visual_hires(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
            qpx_n,qpy_n,qpz_n,rmassq, &
            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
            hpx_n,hpy_n,hpz_n,rmassh, &
            orho_n,opresx_n,opresy_n,opresz_n, &
            opx_n,opy_n,opz_n,rmasso, &
            epres_n,bx_n,by_n,bz_n,bx0_n,by0_n,bz0_n, &
            bsx_n,bsy_n,bsz_n, &
            curx_n,cury_n,curz_n,efldx_n,efldy_n,efldz_n, &
            tvx_n,tvy_n,tvz_n, &
            tx_n,ty_n,tz_n,tg1_n,tg2_n,tt_n,work_n, &
            mx_n,my_n,mz_n,mz2_n,muvwp2_n, &
            nx_n,ny_n,nz_n,ngrd_n,xspac_n, &
            cross_n,along_n,flat_n,xcraft,ncraft,re_equiv, &
            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
            grd_zmin_n,grd_zmax_n,ut,b_equiv, &
            ti_te_moon,rho_equiv)
        !
        write(6,79) nchf
        79 format('  restart from   mpd3d',i2)
        rewind nchf
        nchf=11
    else
    
        !
        !     ******************************
        !            initialization
        !     ******************************
        !
        !     ambient plasma
        !
        !      initialize indices and surface points
        !
        numsrf=0
        nummid=0
        numzero=0
        !
        !       rotation parameters
        !
        sin_rot=sin(rot_angle)
        cos_rot=cos(rot_angle)
        write(6,*)'init',lunar_dist,rot_angle
        !
        !      scale lengths for plasma sheet population
        !
        sheet_den=0.25
        alpha_s=4.
        !
        do m=1,ngrd
            dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
            dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
            dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
            !
            !      create dipole magnetic field and load magnetospheric plasma
            !
            do k=1,nz
                az=(grd_zmin(m)+dz*(k-1)-zdip)
                !
                do j=1,ny
                    ay=grd_ymin(m)+dy*(j-1)-ydip
                    !
                    do i=1,nx
                        ax=(grd_xmin(m)+dx*(i-1)-xdip)
                        !
                        !       rotate corrdinates for planet motion
                        !
                        xr=ax*cos_rot+ay*sin_rot
                        yr=-ax*sin_rot+ay*cos_rot
                        !
                        !       tilt space to dipole space
                        !
                        xp=xr*cos_tilt-az*sin_tilt
                        yp=yr
                        zp=xr*sin_tilt+az*cos_tilt
                        ar=sqrt(xp**2+yp**2+zp**2)
                        radius=ar*re_equiv
    
                        !
                        !        determine magnetic dipole field
                        !
                        call dipole(bx0(i,j,k,m),by0(i,j,k,m), &
                            bz0(i,j,k,m),ax,ay,az)
                        !
                        !         zero interior magnetic field so alfven speed small
                        !
                        if (ar.lt.rearth-1.5)then
                            bx0(i,j,k,m)=0.
                            by0(i,j,k,m)=0.
                            bz0(i,j,k,m)=0.
                        endif
                        !
                        !        set up rotational properties
                        !
                        rx=ax*re_equiv
                        ry=ay*re_equiv
                        rd=sqrt(rx**2+ry**2)
                        if(rd.lt.r_rot)then
                            vfrac=1.
                        else
                            vfrac=((2.*r_rot**2)/(r_rot**2+ rd**2))**2
                        endif
                        rvy=rx*v_rot*vfrac
                        rvx=-ry*v_rot*vfrac
                        !
                        !        for jupiter top hat distribution
    
                        ar_iono=sqrt(xp**2+ay**2+1.25*zp**2)
                        !        isotropic
                        !        ar_iono=sqrt(xp**2+ay**2+zp**2)
                        !
                        ar_iono=amax1(0.0001,ar_iono)
                        ra=((ar_iono+0.5*rearth)/(1.5*rearth))**(-alpha_e)
                        zheight=amax1(1.,(zp**2+(1.5*rearth)**2)/(3.0*rearth)**2)
                        ra=ra/zheight**2
                        rho_iono=amax1(erho*ra,0.001)
                        !
                        r_equat=(ar**3+0.001)/(xp**2+ay**2+0.001)
                        r_equat=amax1(r_equat,rearth)
                        rerg_sphere=(0.001+(rearth/r_equat)**4)
                        !
                        if(r_equat.gt.5.*rearth)then
                            rerg_sphere=(0.001+(rearth/r_equat)**alpha_e)!constant pressure flux
                            rden_sphere=1.
                        else
                            rerg_sphere=0.1+0.9*(r_equat-rearth)/(4.*rearth)
                            rden_sphere=(rerg_sphere**2)          !reduce o+ in plasmasphere
                        endif
                        !
                        arho=amax1(rho_iono,d_min)
                        vt=amax1(sqrt(rvy**2+rvx**2),cs_inner)
                        !        vt=amin1(vt,20.*cs_inner)
                        !
                        !
                        !
                        qrho(i,j,k,m)=arho*rmassq/zheight
                        qpresx(i,j,k,m)=arho*vt
                        qpresy(i,j,k,m)=qpresx(i,j,k,m)
                        qpresz(i,j,k,m)=qpresx(i,j,k,m)
                        qpresxy(i,j,k,m)=0.
                        qpresxz(i,j,k,m)=0.
                        qpresyz(i,j,k,m)=0.
                        qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
                        qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
                        qpz(i,j,k,m)=0.
                        !
                        !          add 1% heavies everywhere
                        !
    
                        hrho(i,j,k,m)=0.005*qrho(i,j,k,m)*rmassh/rmassq
                        hpresx(i,j,k,m)=0.005*qpresx(i,j,k,m)
                        hpresy(i,j,k,m)=hpresx(i,j,k,m)
                        hpresz(i,j,k,m)=hpresx(i,j,k,m)
                        hpresxy(i,j,k,m)=0.
                        hpresxz(i,j,k,m)=0.
                        hpresyz(i,j,k,m)=0.
                        hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
                        hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
                        hpz(i,j,k,m)=0.
                        !
                        orho(i,j,k,m)=0.005*o_conc*qrho(i,j,k,m)*rmasso/rmassq
                        opresx(i,j,k,m)=0.005*o_conc*qpresx(i,j,k,m)
                        opresy(i,j,k,m)=opresx(i,j,k,m)
                        opresz(i,j,k,m)=opresx(i,j,k,m)
                        opresxy(i,j,k,m)=0.
                        opresxz(i,j,k,m)=0.
                        opresyz(i,j,k,m)=0.
                        opx(i,j,k,m)=orho(i,j,k,m)*rvx
                        opy(i,j,k,m)=orho(i,j,k,m)*rvy
                        opz(i,j,k,m)=0.
                        !
                        epres(i,j,k,m)=(qpresx(i,j,k,m)+hpresx(i,j,k,m)+ &
                            opresx(i,j,k,m))/ti_te
                        !
                        !          check for io - thermal mach number set at 1/4
                        amach=4.
                        !
                        ar_io=sqrt(ax**2+ay**2)*re_equiv
                        dr_io=abs(ar_io -lunar_dist)
    
                        if(abs(dr_io.lt.2.*torus_rad)) then
                            rscale=exp(-((dr_io)/(0.5*torus_rad))**2) ! scale height in rjs
                            zscale=exp(-((zp*re_equiv)/(0.5*torus_rad))**2)
                            !
                            oden=den_lunar*rmasso*rscale*zscale
                            !            write(6,*)i,j,k,m,oden
                            orho(i,j,k,m)=orho(i,j,k,m) +oden
                            opresx(i,j,k,m)=opresx(i,j,k,m) &
                                +amach**2*oden*(v_rot**2)/rmasso  !temp goes as v**2
                            opresy(i,j,k,m)=opresx(i,j,k,m)
                            opresz(i,j,k,m)=opresx(i,j,k,m)
                            opresxy(i,j,k,m)=0.
                            opresxz(i,j,k,m)=0.
                            opresyz(i,j,k,m)=0.
                            opx(i,j,k,m)=opx(i,j,k,m)+reduct**oden*rvx
                            opy(i,j,k,m)=opy(i,j,k,m)+reduct*oden*rvy
                            !
                            hden=oden*rmassh/rmasso/o_conc
                            hrho(i,j,k,m)=hrho(i,j,k,m) +hden
                            hpresx(i,j,k,m)=hpresx(i,j,k,m) &
                                +amach**2*hden*(v_rot**2)/rmassh  !temp goes as v**2
                            hpresy(i,j,k,m)=hpresx(i,j,k,m)
                            hpresz(i,j,k,m)=hpresx(i,j,k,m)
                            hpresxy(i,j,k,m)=0.
                            hpresxz(i,j,k,m)=0.
                            hpresyz(i,j,k,m)=0.
                            hpx(i,j,k,m)=hpx(i,j,k,m)+reduct*hden*rvx
                            hpy(i,j,k,m)=hpy(i,j,k,m)+reduct*hden*rvy
                            !
                            epres(i,j,k,m)=(qpresx(i,j,k,m)+hpresx(i,j,k,m)+ &
                                opresx(i,j,k,m))/ti_te
                        endif
                        !
                        bx(i,j,k,m)=0.
                        by(i,j,k,m)=0.
                        bz(i,j,k,m)=0.
                        !
                        !        test for boundary point of planets surface or
                        !        interior to planet
                        !
                        if((ar.le.rearth+.6).and.(m.le.mbndry))then
                            !
                            if(ar.lt.rearth-1.5) then
                                numzero(m)=numzero(m)+1
                                ijzero(m,1,numzero(m))=i
                                ijzero(m,2,numzero(m))=j
                                ijzero(m,3,numzero(m))=k
                                !
                                !
                                parm_zero(m,1,numzero(m))=qrho(i,j,k,m)
                                parm_zero(m,2,numzero(m))=hrho(i,j,k,m)
                                parm_zero(m,3,numzero(m))=orho(i,j,k,m)
                                parm_zero(m,4,numzero(m))=qpresx(i,j,k,m)
                                parm_zero(m,5,numzero(m))=hpresx(i,j,k,m)
                                parm_zero(m,6,numzero(m))=opresx(i,j,k,m)
                                parm_zero(m,7,numzero(m))=epres(i,j,k,m)
                                !
                                !         qpx(i,j,k,m)=qrho(i,j,k,m)*rvx
                                !         qpy(i,j,k,m)=qrho(i,j,k,m)*rvy
                                !         qpz(i,j,k,m)=0.
                                !         hpx(i,j,k,m)=hrho(i,j,k,m)*rvx
                                !         hpy(i,j,k,m)=hrho(i,j,k,m)*rvy
                                !         hpz(i,j,k,m)=0.
                                !         opx(i,j,k,m)=orho(i,j,k,m)*rvx
                                !         opy(i,j,k,m)=orho(i,j,k,m)*rvy
                                !         opz(i,j,k,m)=0.
                                !
                            else  if(ar.lt.rearth-.5) then
                                nummid(m)=nummid(m)+1
                                ijmid(m,1,nummid(m))=i
                                ijmid(m,2,nummid(m))=j
                                ijmid(m,3,nummid(m))=k
                                ar2=sqrt(xp**2+ay**2)
                                parm_mid(m,1,nummid(m))=qrho(i,j,k,m)
                                parm_mid(m,2,nummid(m))=hrho(i,j,k,m)
                                parm_mid(m,3,nummid(m))=orho(i,j,k,m)
                                parm_mid(m,4,nummid(m))=qpresx(i,j,k,m)
                                parm_mid(m,5,nummid(m))=hpresx(i,j,k,m)
                                parm_mid(m,6,nummid(m))=opresx(i,j,k,m)
                                parm_mid(m,7,nummid(m))=epres(i,j,k,m)
                                !
                            else
                                numsrf(m)=numsrf(m)+1
                                ijsrf(m,1,numsrf(m))=i
                                ijsrf(m,2,numsrf(m))=j
                                ijsrf(m,3,numsrf(m))=k
                                ar2=sqrt(xp**2+ay**2)
                                parm_srf(m,1,numsrf(m))=qrho(i,j,k,m)
                                parm_srf(m,2,numsrf(m))=hrho(i,j,k,m)
                                parm_srf(m,3,numsrf(m))=orho(i,j,k,m)
                                parm_srf(m,4,numsrf(m))=qpresx(i,j,k,m)
                                parm_srf(m,5,numsrf(m))=hpresx(i,j,k,m)
                                parm_srf(m,6,numsrf(m))=opresx(i,j,k,m)
                                parm_srf(m,7,numsrf(m))=epres(i,j,k,m)
                            endif
                        endif  ! end bndry condition test
                        !
                    enddo
                enddo
            enddo
        enddo
        !
        write(6,*)'interior and zero points'
        do m=1,mbndry
            write(6,*)m,numsrf(m),nummid(m),numzero(m)
        enddo
        !
        !     initialize solar wind plasa can be placed beyond
        !       the earth at a radius of re_wind
        !
        wind_bnd=r_rot/re_equiv
        ofrac=rho_frac
        do m=1,ngrd
            !
            dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
            dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
            dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
            !
            do k=1,nz
                az=(grd_zmin(m)+dz*(k-1)-zdip)
                !
                do j=1,ny
                    ay=grd_ymin(m)+dy*(j-1)-ydip
                    do i=1,nx
                        ax=(grd_xmin(m)+dx*(i-1)-xdip)
    
                        ar=sqrt(ax**2+ay**2+az**2)
                        if((ar.ge.1.5*wind_bnd).and.(ax.lt.0.))then
                            qrho(i,j,k,m)=srho+qrho(i,j,k,m)
                            qpx(i,j,k,m)=spx+qpx(i,j,k,m)
                            qpy(i,j,k,m)=spy+qpy(i,j,k,m)
                            qpz(i,j,k,m)=spz+qpz(i,j,k,m)
                            !
                            avx=qpx(i,j,k,m)/qrho(i,j,k,m)
                            avy=qpy(i,j,k,m)/qrho(i,j,k,m)
                            avz=qpz(i,j,k,m)/qrho(i,j,k,m)
                            !
                            apres=cs_wind**2*qrho(i,j,k,m)/gamma
                            qpresx(i,j,k,m)=0.5*apres+qpresx(i,j,k,m)
                            qpresy(i,j,k,m)=qpresx(i,j,k,m)
                            qpresz(i,j,k,m)=qpresx(i,j,k,m)
                            qpresxy(i,j,k,m)=0.
                            qpresxz(i,j,k,m)=0.
                            qpresyz(i,j,k,m)=0.
                            !
                            epres(i,j,k,m)=0.5*apres/ti_te+epres(i,j,k,m)
                            !
                            hrho(i,j,k,m)=srho*rho_frac*rmassh/rmassq+ &
                                hrho(i,j,k,m)
                            hpx(i,j,k,m)=avx*hrho(i,j,k,m)
                            hpy(i,j,k,m)=avy*hrho(i,j,k,m)
                            hpz(i,j,k,m)=avz*hrho(i,j,k,m)
                            hpresx(i,j,k,m)=0.5*spress*rho_frac &
                                +hpresx(i,j,k,m)
                            hpresy(i,j,k,m)=hpresx(i,j,k,m)
                            hpresz(i,j,k,m)=hpresx(i,j,k,m)
                            hpresxy(i,j,k,m)=0.
                            hpresxz(i,j,k,m)=0.
                            hpresyz(i,j,k,m)=0.
                            !
                            orho(i,j,k,m)=srho*ofrac*rmasso/rmassq+ &
                                orho(i,j,k,m)
                            opx(i,j,k,m)=avx*orho(i,j,k,m)
                            opy(i,j,k,m)=avy*orho(i,j,k,m)
                            opz(i,j,k,m)=avz*orho(i,j,k,m)
                            opresx(i,j,k,m)=0.5*spress*ofrac &
                                +opresx(i,j,k,m)
                            opresy(i,j,k,m)=opresx(i,j,k,m)
                            opresz(i,j,k,m)=opresx(i,j,k,m)
                            opresxy(i,j,k,m)=0.
                            opresxz(i,j,k,m)=0.
                            opresyz(i,j,k,m)=0.
                        endif

                        if((ar.gt.wind_bnd).and.(ar.lt.1.5*wind_bnd) &
                            .and.(ax.lt.0.0))then
                            dfrac=(ar-wind_bnd)/(0.5*wind_bnd)
                            qrho(i,j,k,m)=srho*dfrac+qrho(i,j,k,m)
                            qpx(i,j,k,m)=spx*dfrac+qpx(i,j,k,m)
                            qpy(i,j,k,m)=spy*dfrac+qpy(i,j,k,m)
                            qpz(i,j,k,m)=spz*dfrac+qpz(i,j,k,m)
                            avx=qpx(i,j,k,m)/qrho(i,j,k,m)
                            avy=qpy(i,j,k,m)/qrho(i,j,k,m)
                            avz=qpz(i,j,k,m)/qrho(i,j,k,m)
                            !
                            apres=cs_wind**2*qrho(i,j,k,m)/gamma
                            qpresx(i,j,k,m)=0.5*apres*dfrac+qpresx(i,j,k,m)
                            qpresy(i,j,k,m)=qpresx(i,j,k,m)
                            qpresz(i,j,k,m)=qpresx(i,j,k,m)
                            qpresxy(i,j,k,m)=0.
                            qpresxz(i,j,k,m)=0.
                            qpresyz(i,j,k,m)=0.
!
                            epres(i,j,k,m)=0.5*apres/ti_te*dfrac+epres(i,j,k,m)
!
                            hrho(i,j,k,m)=srho*rho_frac*rmassh/rmassq*dfrac+ &
                                hrho(i,j,k,m)
                            hpx(i,j,k,m)=avx*hrho(i,j,k,m)
                            hpy(i,j,k,m)=avy*hrho(i,j,k,m)
                            hpz(i,j,k,m)=avz*hrho(i,j,k,m)
                            hpresx(i,j,k,m)=0.5*spress*rho_frac*dfrac &
                                +hpresx(i,j,k,m)
                            hpresy(i,j,k,m)=hpresx(i,j,k,m)
                            hpresz(i,j,k,m)=hpresx(i,j,k,m)
                            hpresxy(i,j,k,m)=0.
                            hpresxz(i,j,k,m)=0.
                            hpresyz(i,j,k,m)=0.
!
                            orho(i,j,k,m)=srho*ofrac*rmasso/rmassq*dfrac+ &
                                orho(i,j,k,m)
                            opx(i,j,k,m)=avx*orho(i,j,k,m)
                            opy(i,j,k,m)=avy*orho(i,j,k,m)
                            opz(i,j,k,m)=avz*orho(i,j,k,m)
                            opresx(i,j,k,m)=0.5*spress*ofrac*dfrac &
                                +opresx(i,j,k,m)
                            opresy(i,j,k,m)=opresx(i,j,k,m)
                            opresz(i,j,k,m)=opresx(i,j,k,m)
                            opresxy(i,j,k,m)=0.
                            opresxz(i,j,k,m)=0.
                            opresyz(i,j,k,m)=0.
                        endif
                        !
                    enddo
                enddo
            enddo
        enddo
        !
        !
        !     check speeds
        !
        do  m=ngrd,1,-1
            do  i=1,nx
                do  j=1,ny
                    do  k=1,nz
                        avx=qpx(i,j,k,m)/qrho(i,j,k,m)
                        avy=qpy(i,j,k,m)/qrho(i,j,k,m)
                        avz=qpz(i,j,k,m)/qrho(i,j,k,m)
                        spd=sqrt(avx**2+avy**2+avz**2)
                        if(spd.gt.1.)then
                            qpx(i,j,k,m)=0.
                            qpy(i,j,k,m)=0.
                            qpz(i,j,k,m)=0.
                        endif
                        !
                        avx=hpx(i,j,k,m)/hrho(i,j,k,m)
                        avy=hpy(i,j,k,m)/hrho(i,j,k,m)
                        avz=hpz(i,j,k,m)/hrho(i,j,k,m)
                        spd=sqrt(avx**2+avy**2+avz**2)
                        if(spd.gt.1.)then
                            hpx(i,j,k,m)=0.
                            hpy(i,j,k,m)=0.
                            hpz(i,j,k,m)=0.
                        endif
                        !
                        avx=opx(i,j,k,m)/orho(i,j,k,m)
                        avy=opy(i,j,k,m)/orho(i,j,k,m)
                        avz=opz(i,j,k,m)/orho(i,j,k,m)
                        spd=sqrt(avx**2+avy**2+avz**2)
                        if(spd.gt.1.)then
                            opx(i,j,k,m)=0.
                            opy(i,j,k,m)=0.
                            opz(i,j,k,m)=0.
                        endif
                        !
                    enddo
                enddo
            enddo
        enddo
    endif ! if(.not.start) then
    
    !
    !     initialized other important stuff
    !
    ut=utstart+t*t_equiv/3600.
    nrot=ut/planet_per
    rot_hrs=ut-nrot*planet_per
    rot_angle=6.2832*rot_hrs/planet_per
    !
    !     initialize plasma resistivity
    !
    call set_resist(resistive,nx,ny,nz,mbndry,resist, &
        ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
        msrf,mmid,mzero,1.)
    call set_resist(resistive_n,nx_n,ny_n,nz_n,mbndry_n,resist, &
        ijzero_n,numzero_n,ijmid_n,nummid_n,ijsrf_n, &
        numsrf_n,msrf_n,mmid_n,mzero_n,0.01)
    !
    !
    !     read down relevant data list to find correct pointer
    !
    if(spacecraft)then
        !
        !      calculate the distance between the reference spacecraft and
        !          solar wind boundary
        !
        m=ngrd
        distance=grd_xmin(m)-rcraft(1)/re_equiv
        !
        !      read down data file until correct time in data file
        !
        !      do  n=1,ncraft
        !      do  n=1,2
        n=1
        mout=40+n
        do while(ut.ge.zcraft(4,n))
            read(mout,*)zcraft(4,n),zcraft(1,n), &
                zcraft(2,n),zcraft(3,n)
            !
            !          change direction to get from gsm to simulation coords
            !
            zcraft(1,n)=-zcraft(1,n)
            zcraft(2,n)=-zcraft(2,n)
            if(n.eq.1)then
                rcraft(1)=zcraft(1,1)
                rcraft(2)=zcraft(2,1)
                rcraft(3)=zcraft(3,1)
            endif
            !
        end do
        !      end do
        !
        !
        call limcraft(zcraft,ncraft,re_equiv,ngrd, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
            grd_zmin,grd_zmax)
        !
        !      calculate the distance between the reference spacecraft and
        !          solar wind boundary
        !
        distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
        write(6,*)'wind displaced by',distance
        do n=1,ncraft
            do nn=1,4
                xcraft(nn,n)=zcraft(nn,n)
            end do
        end do
        !
        !      do while(ut.gt.vut)
        !           svelx=zvelx
        !           svely=zvely
        !           svelz=zvelz
        !           read(28,*)vut,zvelx,zvely,zvelz
        !            zvelx=-zvelx/v_equiv
        !            zvely=-zvely/v_equiv+0.03
        !            zvelz=zvelz/v_equiv
        !            vut=vut+t_equiv*distance/zvelx/3600.
        !      end do
        !
        !      do while (ut.gt.rut)
        !           srho=zrho
        !           read(27,*)rut,zrho
        !           rut=rut+t_equiv*distance/zvelx/3600.
        !           zrho=zrho/rho_equiv
        !      end do
        !
        !
        !     initialize counting array
        !
        do  j=1,ny
            do  k=1,nz
                ncount(j,k)=0
                future(j,k)=ut-0.01
            enddo
        enddo
        !
        !      read all the magnetic field data to minize data sorting
        !
        do  m=1,ncts
            read(29,*)bfld(m,4),bmag,bfld(m,1),bfld(m,2),bfld(m,3)
            read(27,*)rut,rplas(m),svel(m,1),svel(m,2),svel(m,3)
            !      read(27,*)rut,rplas(m)
            !      read(28,*)vut,svel(m,1),svel(m,2),svel(m,3)
            !      warning recalibration
            !         keep bx in solar wind constant
            !      bfld(m,1)=-sbx_wind*b_equiv
        enddo
        !
        !      set timing
        !
        nvx=0
        vut=-999.
        do while(ut.gt.vut)
            svelx=zvelx
            nvx=nvx+1
            zvelx=-svel(nvx,1)/v_equiv
            vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
        end do
        !
        write(6,*)'ut=',ut,'wind time',bfld(nvx,4)
        !
        displace=0.
        m=ngrd
        dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
        dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
        dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
        !
        do k=1,nz
            do j=1,ny
                do while((ut.gt.future(j,k)) &
                    .and.(ncount(j,k)+1.le.ncts))
                    nc=ncount(j,k)+1
                    bxp(j,k)=bxf(j,k)
                    byp(j,k)=byf(j,k)
                    bzp(j,k)=bzf(j,k)
                    rhop(j,k)=rhof(j,k)
                    svxp(j,k)=svxf(j,k)
                    svyp(j,k)=svyf(j,k)
                    svzp(j,k)=svzf(j,k)
                    past(j,k)=future(j,k)
                    !
                    future(j,k)=bfld(nc,4)
                    bxf(j,k)=-bfld(nc,1)/b_equiv
                    byf(j,k)=-bfld(nc,2)/b_equiv
                    bzf(j,k)=bfld(nc,3)/b_equiv
                    rhof(j,k)=rplas(nc)/rho_equiv
                    svxf(j,k)=-svel(nc,1)/v_equiv
                    svyf(j,k)=0.
                    svzf(j,k)=0.
                    !       svyf(j,k)=-svel(nc,2)/v_equiv +0.03
                    !       svyf(j,k)=-svel(nc,2)/v_equiv
                    !       svzf(j,k)=svel(nc,3)/v_equiv
                    ncount(j,k)=nc
                    avx=svxf(j,k)
                    !
                    !      calculate delay
                    !
                    if(warp)then
                        b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
                        b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
                        ay=grd_ymin(m)+dy*(j-1)-rcraft(2)/re_equiv
                        az=grd_zmin(m)+dz*(k-1)-rcraft(3)/re_equiv
                        !
                        !       going to assume bz imf on average pos and ignore transients
                        !               and by imf is negative
                        displace=-bxf(j,k)* &
                            (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
                    endif
                    ar=distance+displace
                    future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
                enddo
            enddo
        enddo
        !
        call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
            orho,opresx,opresy,opresz,opx,opy,opz, &
            rmassq,rmassh,rmasso,epres, &
            qpresxy,qpresxz,qpresyz, &
            hpresxy,hpresxz,hpresyz, &
            opresxy,opresxz,opresyz, &
            rhop,svxp,svyp,svzp,svelx,spress, &
            ti_te,rho_frac,nx,ny,nz,ngrd)
    endif   ! end spacecraft if
    !
    !     check initial conditions
    !
    !      write(6,*) 'checking set speed'
    do m=ngrd,1,-1
        !
        !       sync time steps
        !
        t_old(m)=0.
        t_new(m)=0.
        !
        !       check density
        !
        call set_rho(qrho,qpresx,qpresy,qpresz, &
            qpresxy,qpresxz,qpresyz,rmassq, &
            hrho,hpresx,hpresy,hpresz, &
            hpresxy,hpresxz,hpresyz,rmassh, &
            orho,opresx,opresy,opresz, &
            opresxy,opresxz,opresyz,rmasso, &
            epres,nx,ny,nz,ngrd,m,o_conc)
        !
        !       check speeds of individual grids
        !
        vlim=0.6*m
        call set_speed_agrd( &
            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
            orho,opresx,opresy,opresz,opx,opy,opz, &
            epres,qpresxy,qpresxz,qpresyz, &
            hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
            bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
            vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
            rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
            pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
            vlim,alf_lim,o_conc,fastest)
        write(6,195)m,csmax,alfmax,pxmax,pymax,pzmax
        195   format(1x,i2,5(1x,1pe12.5))
        !
        t_stepnew(m)=stepsz*xspac(m)/fastest
    enddo
    write(6,*)'speeds checked'
    !
    delt_old=delt
    delt=stepsz/fastest
    delt=amin1(delt,1.25*delt_old)
    delt=amax1(3.e-3,delt)
    write(6,163)t,delt,ut,bfld(nvx,4)
    163 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=', &
        1pe12.5,' wind time',1pe12.5)
    !
    write(6,161)
    161 format(' initialization completed')
    !
    !
    !     ********************************************
    !     start time sequence
    !     ********************************************
    !
    do while(t.lt.tmax)
        !
        !
        !     find maximum velocity to determine time step
        !
        write(6,*)'start main loop'
        !
        do m=ngrd,1,-1
            t_step(m)=t_stepnew(m)
            !
            !       sync time steps
            !
            t_old(m)=0.
            t_new(m)=0.
  
            if(m.eq.ngrd)then
                smallest_step=t_step(m)
                mallest_step=m
                write(6,*)'syncing',m,t_step(m)
            else
                !
                !          check to see if m grid doesnt outstep grid m+1
                write(6,*)'unsync',m,t_step(m)
                if(t_step(m).gt.t_step(m+1))t_step(m)=t_step(m+1)
                !
                if (smallest_step.gt.t_step(m))then
                    smallest_step=t_step(m)
                    mallest_step=m
                endif
            endif
        enddo
        !
        !     set variable steps for each grid box
        !       round_step size off
        !
        !       i_step=smallest_step*10000
        !       smallest_step=i_step/10000.
        !
        !     write(6,*)'unsync steps',t_step,mallest_step,smalslest_step
        !
        do m=1,ngrd
            if(m.le.mallest_step) then
                t_step(m)=smallest_step
                write(6,*)'sync step',m,t_step(m)
            else
                astep=(t_step(m)/t_step(m-1)+.50)  !round up as needed
                nsteps=astep                          !nearest integer
                !           write(6,*)astep,nsteps
                if(nsteps.gt.2)nsteps=2
                if(nsteps.lt.1)nsteps=1
                t_step(m)=t_step(m-1)*nsteps
                write(6,*)'sync step',m,t_step(m)
            endif
        enddo
        !
        m_step=t_step(ngrd)/t_step(mallest_step)
        write(6,*)'lores steps ',mallest_step,m_step,t_step
        !
        !
        mallest_step_n=4
        do m=ngrd_n,1,-1
            t_old_n(m)=0.
            t_new_n(m)=0.
            !
            if(m.eq.ngrd_n)then
                t_step_n(m)=t_step(main_grd_n)
                !    +              *xspac_n(m)/xspac(main_grd_n)
            else if(m.le.mallest_step_n)then
                t_step_n(m)=t_step_n(ngrd_n)*xspac_n(mallest_step_n) &
                    /xspac_n(ngrd_n)
            else
                t_step_n(m)=t_step_n(ngrd_n)*xspac_n(m)/xspac_n(ngrd_n)
            endif
        enddo
        m_step_n=t_step(main_grd_n)/t_step_n(mallest_step_n)
        !
        write(6,*)'hires steps ',mallest_step_n,m_step_n,t_step_n
        !
        delt=t_step(ngrd)
        told=t
        utold=ut
        old_tilt=tilt
        !
        t=t+delt
        ut=utstart+t*t_equiv/3600.
        tilt=tilt+dtilt*delt
        delay=t_equiv*distance/svelx/3600.
        !
        write(6,201)t,delt,ut,bfld(nvx,4)
        201 format(1x,'t=',1pe12.5,' dt=',1pe12.5,' ut=', &
            1pe12.5,' wind time',1pe12.5)
        !
        nrot=ut/planet_per
        rot_hrs=ut-nrot*planet_per
        rot_angle=6.2832*rot_hrs/planet_per
        !
        !     write out data if necessary - only using course gridding
        !        the momemt
        !
        !     if(spacecraft) then
        !
        !         update position of moon diagnostics
        !
        xcraft(1,2)=xmoon*re_equiv*1.05
        xcraft(2,2)=ymoon*re_equiv*1.05
        xcraft(3,2)=zmoon*re_equiv*1.05
        !
        do n=1,ncraft
            call qvset(0.,bsx,nx*ny*nz)
            call qvset(0.,bsy,nx*ny*nz)
            call qvset(0.,bsz,nx*ny*nz)
            !
            !
            m=1
            do while ((xcraft(1,n).gt.grd_xmax(m)*re_equiv).or. &
                (xcraft(1,n).lt.grd_xmin(m)*re_equiv).or. &
                (xcraft(2,n).gt.grd_ymax(m)*re_equiv).or. &
                (xcraft(2,n).lt.grd_ymin(m)*re_equiv).or. &
                (xcraft(3,n).gt.grd_zmax(m)*re_equiv).or. &
                (xcraft(3,n).lt.grd_zmin(m)*re_equiv).and. &
                (m+1.le.ngrd))
                m=m+1
            enddo
            rx=xspac(m)
            ry=rx
            rz=rz
            !
            call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
            call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
            call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
            !
            add_dip=.false.
            call crafdatv(bsx,bsy,bsz, &
                qpx,qpy,qpz,qrho,qpresx,qpresy,qpresz,rmassq, &
                hpx,hpy,hpz,hrho,hpresx,hpresy,hpresz,rmassh, &
                opx,opy,opz,orho,opresx,opresy,opresz,rmasso, &
                epres,nx,ny,nz,ngrd,m,xcraft,ncraft,n,ut, &
                re_equiv,b_equiv,v_equiv,rho_equiv,gamma, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax)
        enddo
        !
        !      calculate flux from torus
        !
        write(3,*)'main grid'
        do m=ngrd,1,-1
            scale=rho_equiv*v_equiv*1.e3* &
            (xspac(m)*planet_rad*re_equiv*1.e5)**2
            call flux_counter(qpx,qpy,qpz,hpx,hpy,hpz, &
                opx,opy,opz,nx,ny,nz,ngrd,m, &
                rmassq,rmassh,rmasso, &
                scale,qflux_in,qflux_out,hflux_in,hflux_out, &
                oflux_in,oflux_out)
            write(3,*)ut,m,qflux_in,hflux_in,oflux_in,'grd_in'
            write(3,*)ut,m,qflux_out,hflux_out,oflux_out,'grd_out'
        enddo
        !
        !      calculate flux from the moon
        !
        write(3,*)'grid n'
        do m=ngrd_n,1,-1
            scale=rho_equiv*v_equiv*1.e3* &
                (xspac_n(m)*planet_rad*re_equiv*1.e5)**2
            call flux_counter(qpx_n,qpy_n,qpz_n,hpx_n,hpy_n,hpz_n, &
                opx_n,opy_n,opz_n,nx_n,ny_n,nz_n,ngrd_n,m, &
                rmassq,rmassh,rmasso, &
                scale,qflux_in,qflux_out,hflux_in,hflux_out, &
                oflux_in,oflux_out)
            write(3,*)ut,m,qflux_in,hflux_in,oflux_in,'moon in'
            write(3,*)ut,m,qflux_out,hflux_out,oflux_out,'moon out'
        enddo
  
        !
        !      endif
        !
        !     test to see whether scraft positions need to be updated
        !
        if(spacecraft)then
            !      do 210 n=1,ncraft
            !      do 210 n=1,2
            n=1
            if(ut.ge.zcraft(4,n))then
                mout=40+n
                read(mout,*)zcraft(4,n),zcraft(1,n), &
                    zcraft(2,n),zcraft(3,n)
                !          if(n.eq.1)zcraft(4,n)=zcraft(4,n)/3600.
                !
                !          change direction to get from gsm to simulation coords
                !
                zcraft(1,n)=-zcraft(1,n)
                zcraft(2,n)=-zcraft(2,n)
                if(n.eq.1)then
                    rcraft(1)=zcraft(1,1)
                    rcraft(2)=zcraft(2,1)
                    rcraft(3)=zcraft(3,1)
                endif
                !
            endif
            !
            zcraft(4,3)=zcraft(4,2)
            zcraft(4,4)=zcraft(4,2)
            !
            !     zcraft(1,2)=zcraft(1,2)
            !     zcraft(2,2)=zcraft(2,2)-1.2
            !     zcraft(3,2)=zcraft(3,2)
            !
            !         set refernce spacecraft position
            !                   and spacecraft limits
            !
            call limcraft(zcraft,ncraft,re_equiv,ngrd, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax)
            !
            !         set density and velocity
            !
            distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
            do while (ut.ge.vut)
                nvx=nvx+1
                zvelx=-svel(nvx,1)/v_equiv
                zvely=-svel(nvx,2)/v_equiv
                !         zvely=-svel(nvx,2)/v_equiv+0.03
                zvelz=svel(nvx,3)/v_equiv
                vut=bfld(nvx,4)+t_equiv*distance/zvelx/3600.
                !          read(28,*)vut,zvelx,zvely,zvelz
                !            zvelx=-zvelx/v_equiv
                !            zvely=-zvely/v_equiv+0.03
                !            zvelz=zvelz/v_equiv
                !            vut=vut+t_equiv*distance/zvelx/3600.
            end do
            !         do while (ut.ge.rut)
            !          read(27,*)rut,zrho
            !            rut=rut+t_equiv*distance/zvelx/3600.
            !            zrho=zrho/rho_equiv
            !         end do
            !
            !        fix up magnetic field
            !
            displace=0.
            do k=1,nz
                do j=1,ny
                    do while((ut.ge.future(j,k)) &
                        .and.(ncount(j,k)+1.le.ncts))
                        nc=ncount(j,k)+1
                        future(j,k)=bfld(nc,4)
                        bxf(j,k)=-bfld(nc,1)/b_equiv
                        byf(j,k)=-bfld(nc,2)/b_equiv
                        bzf(j,k)=bfld(nc,3)/b_equiv
                        rhof(j,k)=rplas(nc)/rho_equiv
                        svxf(j,k)=-svel(nc,1)/v_equiv
                        svyf(j,k)=0.00
                        svzf(j,k)=0.0
                        !          svyf(j,k)=-svel(nc,2)/v_equiv
                        !          svyf(j,k)=-svel(nc,2)/v_equiv+0.03
                        !          svzf(j,k)=svel(nc,3)/v_equiv
                        avx=svxf(j,k)
                        ncount(j,k)=nc
                        !
                        !      calculate delay
                        !
                        if(warp)then
                            !           b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
                            !           b_perp=amax1(b_perp,0.1*abs(bxf(j,k)))
                            ay=(j-1.)*xspac(ngrd)+grd_ymin(ngrd)-rcraft(2)/re_equiv
                            az=(k-1.)*xspac(ngrd)+grd_zmin(ngrd)-rcraft(3)/re_equiv
                            !
                            !       going to assume bz imf on average pos and ignore transients
                            !               and by imf is negative
                            !
                            b_perp=sqrt(bzf(j,k)**2+byf(j,k)**2)
                            b_perp=amax1(b_perp,0.33*abs(bxf(j,k)))
                            displace=-bxf(j,k)* &
                                (ay*byf(j,k)+az*bzf(j,k))/b_perp**2
                        endif
                        ar=distance+displace
                        future(j,k)=future(j,k)+t_equiv*ar/avx/3600.
                    enddo
                enddo
            enddo
            !
        endif
        !
        if(spacecraft)then
            !
            !     update spacecraft position and magnetic field
            !
            do n=1,ncraft
                dut=(ut-xcraft(4,n))/(zcraft(4,n)-utold)
                xcraft(1,n)=xcraft(1,n)+dut*(zcraft(1,n)-xcraft(1,n))
                xcraft(2,n)=xcraft(2,n)+dut*(zcraft(2,n)-xcraft(2,n))
                xcraft(3,n)=xcraft(3,n)+dut*(zcraft(3,n)-xcraft(3,n))
                xcraft(4,n)=ut
            enddo
            distance=grd_xmin(ngrd)-rcraft(1)/re_equiv
            call limcraft(xcraft,ncraft,re_equiv,ngrd, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax)
            !
            srho=0.
            do j=1,ny
                do k=1,nz
                    dut=(ut-utold)/(future(j,k)-utold)
                    bxp(j,k)=bxp(j,k)+dut*(bxf(j,k)-bxp(j,k))
                    byp(j,k)=byp(j,k)+dut*(byf(j,k)-byp(j,k))
                    bzp(j,k)=bzp(j,k)+dut*(bzf(j,k)-bzp(j,k))
                    rhop(j,k)=rhop(j,k)+dut*(rhof(j,k)-rhop(j,k))
                    svxp(j,k)=svxp(j,k)+dut*(svxf(j,k)-svxp(j,k))
                    svyp(j,k)=svyp(j,k)+dut*(svyf(j,k)-svyp(j,k))
                    svzp(j,k)=svzp(j,k)+dut*(svzf(j,k)-svzp(j,k))
                    srho=srho+rhop(j,k)
                end do
            end do
            !
            dut=(ut-utold)/(vut-utold)
            svelx=svelx+dut*(zvelx-svelx)
            svely=0.
            svelz=svelz+delvz_wind*delt
            srho=srho/float(nz*ny)
            !
            sbx_wind=sbx_wind+delbx_wind*delt
            sby_wind=sby_wind+delby_wind*delt
            sbz_wind=sbz_wind+delbz_wind*delt
            spx=srho*svelx
            spy=srho*svely
            spz=srho*svelz
            spress=(cs_wind**2*srho/gamma)/gamma1
            serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
        else
            sbx_wind=sbx_wind+delbx_wind*delt
            sby_wind=sby_wind+delby_wind*delt
            sbz_wind=sbz_wind+delbz_wind*delt
            !
            svelx=svelx+delvx_wind*delt
            svely=svely+delvy_wind*delt
            svelz=svelz+delvz_wind*delt
            srho=srho+delrho*delt
            spx=srho*svelx
            spy=srho*svely
            spz=srho*svelz
            spress=(cs_wind**2*srho/gamma)/gamma1
            serg=0.5*(svelx**2+svely**2+svelz**2)*srho+spress
        endif

        delay=t_equiv*distance/svelx/3600.
        !
        !
        !     determine if it is time to move hires grid
        !
        grid_reset=.false.
        do m_n=ngrd_n,1,-1
            dx=(grd_xmax_n(m_n)-grd_xmin_n(m_n))/(nx_n-1.)
            dy=(grd_ymax_n(m_n)-grd_ymin_n(m_n))/(ny_n-1.)
            avx=0.5*(vx_moon+grd_vx_n(m_n))
            avy=0.5*(vy_moon+grd_vy_n(m_n))
            step_t=t-grd_time_n(m_n)
            delx=step_t*avx+grd_dx_n(m_n)
            dely=step_t*avy+grd_dy_n(m_n)
            !
            !       grids must move by integrals of 2
            !
            alx=delx/dx
            lx=alx
            lx=lx/2
            aly=dely/dy
            ly=aly
            ly=ly/2
            !       write(6,*)'shftgrd #',m_n,step_t,t,grd_time_n(m_n)
            !       write(6,*)'shftgrd ##',delx,dely,alx,aly,dx,dy,lx,ly
            if((lx.ne.0).or.(ly.ne.0))then
                grid_reset=.true.
                !
                !       reset grid parameters
                !
                grd_time_n(m_n)=t
                grd_dx_n(m_n)=delx-dx*lx*2
                grd_dy_n(m_n)=dely-dy*ly*2
                grd_vx_n(m_n)=vx_moon
                grd_vy_n(m_n)=vy_moon
                !
                !        reset grid
                !
                agrd_xmin=grd_xmin_n(m_n)+dx*lx*2
                agrd_xmax=grd_xmax_n(m_n)+dx*lx*2
                agrd_ymin=grd_ymin_n(m_n)+dy*ly*2
                agrd_ymax=grd_ymax_n(m_n)+dy*ly*2
                !
                !        reset grid
                !
                if(m_n.lt.ngrd_n)then
                    if( (agrd_xmin.gt.grd_xmin_n(m_n+1)).and. &
                        (agrd_xmax.lt.grd_xmax_n(m_n+1)) )then
                        grd_xmin_n(m_n)=agrd_xmin
                        grd_xmax_n(m_n)=agrd_xmax
                    endif

                    if( (agrd_ymin.gt.grd_ymin_n(m_n+1)).and. &
                        (agrd_ymax.lt.grd_ymax_n(m_n+1)) )then
                        grd_ymin_n(m_n)=agrd_ymin
                        grd_ymax_n(m_n)=agrd_ymax
                    endif
                else
                    if( (agrd_xmin.gt.grd_xmin(main_grd_n)).and. &
                        (agrd_xmax.lt.grd_xmax(main_grd_n)) )then
                        grd_xmin_n(m_n)=agrd_xmin
                        grd_xmax_n(m_n)=agrd_xmax
                    endif
                    if( (agrd_ymin.gt.grd_ymin(main_grd_n)).and. &
                        (agrd_ymax.lt.grd_ymax(main_grd_n)) )then
                        grd_ymin_n(m_n)=agrd_ymin
                        grd_ymax_n(m_n)=agrd_ymax
                    endif
                endif
                !
                write(6,*)'grid shft',grd_xmin_n(m_n),grd_xmax_n(m_n), &
                    grd_ymin_n(m_n),grd_ymax_n(m_n)
                !
                !         species 1
                !
                call shift_grd(qrho,nx,ny,nz,ngrd,main_grd_n, &
                    qrho_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpresx,nx,ny,nz,ngrd,main_grd_n, &
                    qpresx_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpresy,nx,ny,nz,ngrd,main_grd_n, &
                    qpresy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpresz,nx,ny,nz,ngrd,main_grd_n, &
                    qpresz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpresxy,nx,ny,nz,ngrd,main_grd_n, &
                    qpresxy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpresxz,nx,ny,nz,ngrd,main_grd_n, &
                    qpresxz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpresyz,nx,ny,nz,ngrd,main_grd_n, &
                    qpresyz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpx,nx,ny,nz,ngrd,main_grd_n, &
                    qpx_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpy,nx,ny,nz,ngrd,main_grd_n, &
                    qpy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(qpz,nx,ny,nz,ngrd,main_grd_n, &
                    qpz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                !
                !         species 2
                !
                call shift_grd(hrho,nx,ny,nz,ngrd,main_grd_n, &
                    hrho_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpresx,nx,ny,nz,ngrd,main_grd_n, &
                    hpresx_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpresy,nx,ny,nz,ngrd,main_grd_n, &
                    hpresy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpresz,nx,ny,nz,ngrd,main_grd_n, &
                    hpresz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpresxy,nx,ny,nz,ngrd,main_grd_n, &
                    hpresxy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpresxz,nx,ny,nz,ngrd,main_grd_n, &
                    hpresxz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpresyz,nx,ny,nz,ngrd,main_grd_n, &
                    hpresyz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpx,nx,ny,nz,ngrd,main_grd_n, &
                    hpx_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpy,nx,ny,nz,ngrd,main_grd_n, &
                    hpy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(hpz,nx,ny,nz,ngrd,main_grd_n, &
                    hpz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                !
                !         species 3
                !
                call shift_grd(orho,nx,ny,nz,ngrd,main_grd_n, &
                    orho_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opresx,nx,ny,nz,ngrd,main_grd_n, &
                    opresx_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opresy,nx,ny,nz,ngrd,main_grd_n, &
                    opresy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opresz,nx,ny,nz,ngrd,main_grd_n, &
                    opresz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opresxy,nx,ny,nz,ngrd,main_grd_n, &
                    opresxy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opresxz,nx,ny,nz,ngrd,main_grd_n, &
                    opresxz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opresyz,nx,ny,nz,ngrd,main_grd_n, &
                    opresyz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opx,nx,ny,nz,ngrd,main_grd_n, &
                    opx_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opy,nx,ny,nz,ngrd,main_grd_n, &
                    opy_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(opz,nx,ny,nz,ngrd,main_grd_n, &
                    opz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                !
                !       electron pressure
                !
                call shift_grd(epres,nx,ny,nz,ngrd,main_grd_n, &
                    epres_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                !
                !        magnetic field
                !
                call shift_grd(bx,nx,ny,nz,ngrd,main_grd_n, &
                    bx_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(by,nx,ny,nz,ngrd,main_grd_n, &
                    by_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(bz,nx,ny,nz,ngrd,main_grd_n, &
                    bz_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                !
                call shift_grd(bx0,nx,ny,nz,ngrd,main_grd_n, &
                    bx0_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(by0,nx,ny,nz,ngrd,main_grd_n, &
                    by0_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                call shift_grd(bz0,nx,ny,nz,ngrd,main_grd_n, &
                    bz0_n,nx_n,ny_n,nz_n,ngrd_n,m_n,lx,ly, &
                    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                    grd_zmin,grd_zmax, &
                    grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                    grd_zmin_n,grd_zmax_n)
                !
                if(m_n.le.mbndry_n) then
                    theta=((ut-ut_insert)/ut_orbit)*2.*3.1414
                    xmoon=r_orbit*sin(theta)
                    ymoon=-r_orbit*cos(theta)
                    zmoon=0.
                    vx_moon=v_orbit*cos(theta)
                    vy_moon=v_orbit*sin(theta)
                    vz_moon=0.
                    !        write(66,*)'new moon',ut,ut_insert,theta
                    !        write(66,*)' re,spd ',xmoon*re_requiv,ymoon*re_equiv,
                    !    +                      vx_moon*v_equiv,vy_moon*v_equiv
                    !
                    call set_bndry_moon_ram(rmassq,rmassh,rmasso, &
                        m_n,nx_n,ny_n,nz_n,ngrd_n, &
                        parm_srf_n,parm_mid_n,parm_zero_n, &
                        ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                        numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                        qden_moon,hden_moon,oden_moon, &
                        vx_moon,vy_moon,tempi,gamma, &
                        ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_e, &
                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                        grd_zmin_n,grd_zmax_n)
                    write(6,*)'shifted moon'
                    call bndry_moon(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                        qpresxy_n,qpresxz_n,qpresyz_n, &
                        qpx_n,qpy_n,qpz_n,rmassq, &
                        hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                        hpresxy_n,hpresxz_n,hpresyz_n, &
                        hpx_n,hpy_n,hpz_n,rmassh, &
                        orho_n,opresx_n,opresy_n,opresz_n, &
                        opresxy_n,opresxz_n,opresyz_n, &
                        opx_n,opy_n,opz_n,rmasso, &
                        epres_n,bx_n,by_n,bz_n,m_n, &
                        nx_n,ny_n,nz_n,ngrd_n, &
                        parm_srf_n,parm_mid_n,parm_zero_n, &
                        ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                        numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                        qden_moon,hden_moon,oden_moon,cs_moon,gamma, &
                        ti_te_moon,vx_moon,vy_moon,vz_moon, &
                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                        grd_zmin_n,grd_zmax_n)
                endif  ! end bndry reset
            endif  ! end grid move
        enddo
        if(grid_reset)then
            write(7,*)'time', ut, 'moon pos',xmoon,ymoon
            do m_n=1,ngrd_n
                write(7,*)m_n,grd_xmin_n(m_n),grd_xmax_n(m_n), &
                    grd_ymin_n(m_n),grd_ymax_n(m_n)
            enddo
        endif
        !
       
        !
        !     ***********************************************
        !     main grid loop over delt/m_step
        !     **********************************************
        !
        do ms=1,m_step
            do m=ngrd,1,-1
                !
                !     test if grid sector needs to be moved in time
                !
                yes_step=.false.
                !
                if(m.eq.1)then
                    yes_step=.true.
                else
                    if ((t_old(m).eq.t_new(m)).and.(t_new(m).lt.t_step(ngrd)) &
                        .and.(abs(t_new(m)-t_new(1)).le.0.005*t_step(1)) )then
                        yes_step=.true.
                    endif
                endif
                !
                !      time step grid
                !
                if(yes_step) then
                    t_old(m)=t_new(m)
                    t_new(m)=t_old(m)+t_step(m)
                    !
                    write(6,*)'lores time stepping',m, t_old(m),t_new(m)
                    !
                    delt= t_step(m)
                    delt2=delt/2.
                    !
                    if(tilting)then
                        atilt=old_tilt+(t_old(m)+delt2)*dtilt
                        sin_tilt=sin(atilt*.0174533)
                        cos_tilt=cos(atilt*.0174533)
                        ut2=utold+(t_old(m)+delt2)*t_equiv/3600.
                        nrot2=ut2/planet_per
                        rot_hr2=ut2-nrot2*planet_per
                        rot_angle=6.2832*rot_hr2/planet_per
                        !        write(6,*)'mak_dip half with', t_old(m),delt2
                        call mak_dip_grd(bx0,by0,bz0,nx,ny,nz, &
                            ngrd,mbndry,m,ijzero,numzero,mzero,rearth, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
                    !
                    !     *******************************************************
                    !     store initial plasma parameters
                    !     ******************************************************
                    !c
                    !     store initial plasma parameters
                    !
                    call store_array(oldqrho,qrho,nx,ny,nz,ngrd,m)
                    call store_array(oldqpresx,qpresx,nx,ny,nz,ngrd,m)
                    call store_array(oldqpresy,qpresy,nx,ny,nz,ngrd,m)
                    call store_array(oldqpresz,qpresz,nx,ny,nz,ngrd,m)
                    call store_array(oldqpresxy,qpresxy,nx,ny,nz,ngrd,m)
                    call store_array(oldqpresxz,qpresxz,nx,ny,nz,ngrd,m)
                    call store_array(oldqpresyz,qpresyz,nx,ny,nz,ngrd,m)
                    call store_array(oldqpx,qpx,nx,ny,nz,ngrd,m)
                    call store_array(oldqpy,qpy,nx,ny,nz,ngrd,m)
                    call store_array(oldqpz,qpz,nx,ny,nz,ngrd,m)
                    !
                    call store_array(oldhrho,hrho,nx,ny,nz,ngrd,m)
                    call store_array(oldhpresx,hpresx,nx,ny,nz,ngrd,m)
                    call store_array(oldhpresy,hpresy,nx,ny,nz,ngrd,m)
                    call store_array(oldhpresz,hpresz,nx,ny,nz,ngrd,m)
                    call store_array(oldhpresxy,hpresxy,nx,ny,nz,ngrd,m)
                    call store_array(oldhpresxz,hpresxz,nx,ny,nz,ngrd,m)
                    call store_array(oldhpresyz,hpresyz,nx,ny,nz,ngrd,m)
                    call store_array(oldhpx,hpx,nx,ny,nz,ngrd,m)
                    call store_array(oldhpy,hpy,nx,ny,nz,ngrd,m)
                    call store_array(oldhpz,hpz,nx,ny,nz,ngrd,m)
                    !
                    call store_array(oldorho,orho,nx,ny,nz,ngrd,m)
                    call store_array(oldopresx,opresx,nx,ny,nz,ngrd,m)
                    call store_array(oldopresy,opresy,nx,ny,nz,ngrd,m)
                    call store_array(oldopresz,opresz,nx,ny,nz,ngrd,m)
                    call store_array(oldopresxy,opresxy,nx,ny,nz,ngrd,m)
                    call store_array(oldopresxz,opresxz,nx,ny,nz,ngrd,m)
                    call store_array(oldopresyz,opresyz,nx,ny,nz,ngrd,m)
                    call store_array(oldopx,opx,nx,ny,nz,ngrd,m)
                    call store_array(oldopy,opy,nx,ny,nz,ngrd,m)
                    call store_array(oldopz,opz,nx,ny,nz,ngrd,m)
                    !
                    call store_array(oldepres,epres,nx,ny,nz,ngrd,m)
                    call store_array(oldbx,bx,nx,ny,nz,ngrd,m)
                    call store_array(oldby,by,nx,ny,nz,ngrd,m)
                    call store_array(oldbz,bz,nx,ny,nz,ngrd,m)
                    !
                    !
                    !     **********************************************************
                    !     two step lax-wendroff : step 1
                    !          estimate of the fluid quantites at n+1/2
                    !     **********************************************************
                    !
                    !     store initial plasma parameters
                    !
                    !     calculate standard mhd current j = curl b
                    !
                    !
                    rx=xspac(m)
                    ry=xspac(m)
                    rz=xspac(m)
                    !
                    !      write(6,*)' calcur ing now'
                    call calcur(bx,by,bz,nx,ny,nz,ngrd,m,curx,cury,curz, &
                        rx,ry,rz)
                    !
                    !     find total magnetic field
                    !
                    !     write(6,*)' totbfld ing now'
                    call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
                    call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
                    call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
                    !
                    !     find magnitude of b
                    !
                    call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
                    !
                    !
                    !      write(6,*)' fnd_evel ing now'
                    call fnd_evel(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
                        opx,opy,opz,orho,curx,cury,curz,evx,evy,evz, &
                        tvx,tvy,tvz,nx,ny,nz,ngrd,m, &
                        rmassq,rmassh,rmasso,reynolds)
                    !
                    !     write(6,*)' bande ing now'
                    call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
                        curx,cury,curz,evx,evy,evz,btot, &
                        epres,qrho,hrho,orho,resistive,resist,reynolds, &
                        nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
                        ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
                        rx,ry,rz)
                    !
                    !      write(6,*)' push elec ing now'
                    call push_elec(wrkepres,oldepres,epres,evx,evy,evz, &
                        gamma,gamma1,nx,ny,nz,ngrd,m,0.5*delt, &
                        rx,ry,rz)
                    !
                    !      write(6,*)' push ion 1 ing now'
                    call push_ion(wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                        oldqpx,oldqpy,oldqpz, &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        oldqpresxy,oldqpresxz,oldqpresyz, &
                        qpresxy,qpresxz,qpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,ngrd,m,0.5*delt,grav,re_equiv,reynolds, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                        grd_zmin,grd_zmax,ani_q)
           
                    !      write(6,*)' push ion 2 ing now'
                    call push_ion(wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                        oldhpx,oldhpy,oldhpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        oldhpresxy,oldhpresxz,oldhpresyz, &
                        hpresxy,hpresxz,hpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,ngrd,m,0.5*delt,grav,re_equiv,reynolds, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                        grd_zmin,grd_zmax,ani_h)
                    !
                    !      write(6,*)' push ion 3 ing now'
                    call push_ion(wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        oldorho,oldopresx,oldopresy,oldopresz, &
                        oldopx,oldopy,oldopz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        oldopresxy,oldopresxz,oldopresyz, &
                        opresxy,opresxz,opresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,ngrd,m,0.5*delt,grav,re_equiv,reynolds, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                        grd_zmin,grd_zmax,ani_o)
                    !
                    call push_bfld(wrkbx,wrkby,wrkbz,oldbx,oldby,oldbz, &
                        efldx,efldy,efldz,nx,ny,nz,ngrd,m,0.5*delt, &
                        rx,ry,rz)
                    !
                    !     write(6,988)
                    ! 988 format(' main lax loop now')
                    !
                    !     *************************************************************
                    !     apply boundary conditions
                    !     *************************************************************
                    !
                    !     write(6,989)
                    ! 989 format(' doing bndry conditions')
                    !
                    if(m.eq.ngrd) then
                        call bndry_outer( &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz,wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rmassq,rmassh,rmasso,wrkbx,wrkby,wrkbz, &
                            bx0,by0,bz0,nx,ny,nz,ngrd, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te)
                    else
                        t_grd=t_old(m)+delt2
                        if(t_grd.gt.t_new(m+1))then
                            t_grd=t_new(m+1)
                        endif
                        if (t_grd.lt.t_old(m+1))then
                            t_grd=t_old(m+1)
                        endif
                        !      write(6,*) 'flanks1', m,t_grd, t_old(m+1),t_new(m+1)
                        call bndry_flanks( &
                            wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
                            wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
                            wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkorho,wrkopx,wrkopy,wrkopz, &
                            wrkopresx,wrkopresy,wrkopresz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            wrkepres,wrkbx,wrkby,wrkbz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres,oldbx,oldby,oldbz,vvx, &
                            nx,ny,nz,ngrd,m,t_old,t_new,t_grd, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
                    if(m.le.mbndry)then
                        !      write(6,*)'calling bndry inner'
                        call bndry_inner( &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz,rmassq, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz,rmassh, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz,rmasso, &
                            wrkepres,wrkbx,wrkby,wrkbz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
                            ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
                            mbndry,msrf,mmid,mzero, &
                            erho,epress,alpha_e,ti_te,o_conc, &
                            re_equiv,rearth,sbx_wind,spress, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
                    !
                    !      write(6,*)'lax 1 speeds'
                    !
                    !     check fluid parameters
                    !
                    if(spacecraft)then
                        call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp, &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz, &
                            rmassq,rmassh,rmasso,wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,ngrd)
                    endif
                    call set_rho(wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz,rmassq, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz,rmassh, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopresxy,wrkopresxz,wrkopresyz,rmasso, &
                        wrkepres,nx,ny,nz,ngrd,m,o_conc)
                    !
                    vlim=0.6*m
                    call set_speed_agrd( &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        wrkepres,wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        wrkbx,wrkby,wrkbz,bx0,by0,bz0, &
                        bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest)
                    !
                    !
                    !      ***********************************************************
                    !      lax-wendroff: step 2
                    !            use the predicted value to find corrected value for n+1
                    !      ***********************************************************
                    !
                    if(tilting)then
                        atilt=old_tilt+(t_old(m)+delt)*dtilt
                        sin_tilt=sin(atilt*.0174533)
                        cos_tilt=cos(atilt*.0174533)
                        ut2=utold+(t_old(m)+delt)*t_equiv/3600.
                        nrot2=ut2/planet_per
                        rot_hr2=ut2-nrot2*planet_per
                        rot_angle=6.2832*rot_hr2/planet_per
                        !         write(6,*)'mak_dip full with',t_old(m),delt
                        call mak_dip_grd(bx0,by0,bz0,nx,ny,nz,ngrd, &
                            mbndry,m,ijzero,numzero,mzero,rearth, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
           
                    endif
                    rx=xspac(m)
                    ry=xspac(m)
                    rz=xspac(m)
                    !
                    !     calculate standard mhd current j = curl b
                    !
                    call calcur(wrkbx,wrkby,wrkbz,nx,ny,nz,ngrd,m,curx,cury,curz, &
                        rx,ry,rz)
                    !
                    !     find total magnetic field
                    !
                    call totfld(wrkbx,bx0,bsx,nx,ny,nz,ngrd,m)
                    call totfld(wrkby,by0,bsy,nx,ny,nz,ngrd,m)
                    call totfld(wrkbz,bz0,bsz,nx,ny,nz,ngrd,m)
                    !
                    !     find magnitude of b
                    !
                    call tot_b(btot,bsx,bsy,bsz,nx,ny,nz)
                    !
                    !
                    !     find the  electric field from electron momentum eqn
                    !
                    call fnd_evel(wrkqpx,wrkqpy,wrkqpz,wrkqrho, &
                        wrkhpx,wrkhpy,wrkhpz,wrkhrho, &
                        wrkopx,wrkopy,wrkopz,wrkorho, &
                        curx,cury,curz,evx,evy,evz, &
                        tvx,tvy,tvz,nx,ny,nz,ngrd,m, &
                        rmassq,rmassh,rmasso,reynolds)
                    !c
                    call bande(efldx,efldy,efldz,bsx,bsy,bsz, &
                        curx,cury,curz,evx,evy,evz,btot, &
                        wrkepres,wrkqrho,wrkhrho,wrkorho, &
                        resistive,resist,reynolds, &
                        nx,ny,nz,ngrd,m,rmassq,rmassh,rmasso, &
                        ijmid,nummid,ijzero,mbndry,numzero,mmid,mzero, &
                        rx,ry,rz)
                    !
                    call push_elec(epres,oldepres,wrkepres,evx,evy,evz, &
                        gamma,gamma1,nx,ny,nz,ngrd,m,delt, &
                        rx,ry,rz)
                    !
                    !       write(6,*)' push ion 1 ing now'
                    call push_ion(qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                        oldqpx,oldqpy,oldqpz, &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        qpresxy,qpresxz,qpresyz, &
                        oldqpresxy,oldqpresxz,oldqpresyz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassq, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,ngrd,m,delt,grav,re_equiv,reynolds, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                        grd_zmin,grd_zmax,ani_q)
                    !
                    !      write(6,*)' push ion 2 ing now'
                    call push_ion(hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                        oldhpx,oldhpy,oldhpz, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        hpresxy,hpresxz,hpresyz, &
                        oldhpresxy,oldhpresxz,oldhpresyz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmassh, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,ngrd,m,delt,grav,re_equiv,reynolds, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                        grd_zmin,grd_zmax,ani_h)
                    !
                    !       write(6,*)' push ion 3 ing now'
                    call push_ion(orho,opresx,opresy,opresz,opx,opy,opz, &
                        oldorho,oldopresx,oldopresy,oldopresz, &
                        oldopx,oldopy,oldopz, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        opresxy,opresxz,opresyz, &
                        oldopresxy,oldopresxz,oldopresyz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        bsx,bsy,bsz,btot,efldx,efldy,efldz,rmasso, &
                        vvx,vvy,vvz,tvx,tvy,tvz,gamma,gamma1, &
                        nx,ny,nz,ngrd,m,delt,grav,re_equiv,reynolds, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                        grd_zmin,grd_zmax,ani_o)
                    !
                    call push_bfld(bx,by,bz,oldbx,oldby,oldbz, &
                        efldx,efldy,efldz,nx,ny,nz,ngrd,m,delt, &
                        rx,ry,rz)
                    !
                    !
                    !     write(6,992)
                    ! 992 format(' main 2nd lax loop now')
                    !
                    !     *************************************************************
                    !     apply boundary conditions
                    !     *************************************************************
                    !
                    !
                    if(m.eq.ngrd) then
                        call bndry_outer( &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rmassq,rmassh,rmasso,bx,by,bz, &
                            bx0,by0,bz0,nx,ny,nz,ngrd, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te)
                    else
                        t_grd=t_old(m)+delt
                        if(t_grd.gt.t_new(m+1))then
                            !        write(6,*)'warning on lax2 max',m,t_grd,t_new(m+1),
                            !    +                         t_grd-t_new(m+1)
                            t_grd=t_new(m+1)
                        endif
                        if (t_grd.lt.t_old(m+1))then
                            t_grd=t_old(m+1)
                        endif
                        !       write(6,*)'flanks2', m,t_grd, t_old(m+1),t_new(m+1)
                        call bndry_flanks( &
                            qrho,qpx,qpy,qpz, &
                            qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz, &
                            hpresx,hpresy,hpresz,hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz, &
                            opresx,opresy,opresz,opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres,oldbx,oldby,oldbz,vvx, &
                            nx,ny,nz,ngrd,m,t_old,t_new,t_grd, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
                    !
                    if(m.le.mbndry)then
                        !       write(6,*)'calling bndry inner step 2'
                        call bndry_inner( &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
                            orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
                            epres,bx,by,bz, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
                            ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
                            mbndry,msrf,mmid,mzero, &
                            erho,epress,alpha_e,ti_te,o_conc, &
                            re_equiv,rearth,sbx_wind,spress, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
                    !
                    if(spacecraft)then
                        call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            rmassq,rmassh,rmasso,epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,ngrd)
                    endif
                    !
                    !      write(6,*)'lax 2 speeds'
                    !
                    call set_rho(qrho,qpresx,qpresy,qpresz, &
                        qpresxy,qpresxz,qpresyz,rmassq, &
                        hrho,hpresx,hpresy,hpresz, &
                        hpresxy,hpresxz,hpresyz,rmassh, &
                        orho,opresx,opresy,opresz, &
                        opresxy,opresxz,opresyz,rmasso, &
                        epres,nx,ny,nz,ngrd,m,o_conc)
                    !
                    vlim=0.6*m
                    call set_speed_agrd( &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        epres,qpresxy,qpresxz,qpresyz, &
                        hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
                        bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest)
                    !     .......................................................
                    !     try lapdius smoothing - smoothed results will appear in nt2
                    !     .......................................................
                    !
                    !     write(6,*)' doing smoothing okay'
                    rx=xspac(m)
                    ry=xspac(m)
                    rz=xspac(m)
                    !      write(6,*)'calling lapidus'
                    !
                    !
                    !     lapidus smoothing
                    !
                    !      species 1
                    !
                    !     write(*,*)'lap ion1'
                    call fnd_vel(qpx,qpy,qpz,qrho, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m)
                    call lap_plasma(qrho,qpx,qpy,qpz, &
                        qpresx,qpresy,qpresz, &
                        qpresxy,qpresxz,qpresyz, &
                        wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
                        wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m, &
                        chirho,chipxyz,chierg,delt)
                    !
                    !     species 2
                    !
                    !     write(*,*)'lap ion2'
                    call fnd_vel(hpx,hpy,hpz,hrho, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m)
                    call lap_plasma(hrho,hpx,hpy,hpz, &
                        hpresx,hpresy,hpresz, &
                        hpresxy,hpresxz,hpresyz, &
                        wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
                        wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m, &
                        chirho,chipxyz,chierg,delt)
                    !
                    !     species 3
                    !
                    !     write(*,*)'lap ion3'
                    call fnd_vel(opx,opy,opz,orho, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m)
                    call lap_plasma(orho,opx,opy,opz, &
                        opresx,opresy,opresz, &
                        opresxy,opresxz,opresyz, &
                        wrkorho,wrkopx,wrkopy,wrkopz, &
                        wrkopresx,wrkopresy,wrkopresz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m, &
                        chirho,chipxyz,chierg,delt)
                    !
                    !     electrons
                    !
                    !     write(*,*)'lap electrons'
                    call fnd_vtot(qpx,qpy,qpz,qrho, &
                        hpx,hpy,hpz,hrho, &
                        opx,opy,opz,orho, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m, &
                        rmassq,rmassh,rmasso)
                    call lap_elec(epres,wrkepres, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m, &
                        chirho,chipxyz,chierg,delt)
                    !
                    !     amgnetic field fields
                    !
                    !     write(*,*)'lap bfld'
                    call lap_bfld(bx,by,bz, &
                    wrkbx,wrkby,wrkbz, &
                        curx,cury,curz, &
                        nx,ny,nz,ngrd,m, &
                        chirho,chipxyz,chierg,delt)
                        !
                    !
                    !
                    !     reset boundary conditions
                    !
                    if(m.eq.ngrd) then
                        call bndry_outer( &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz, &
                            wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rmassq,rmassh,rmasso,wrkbx,wrkby,wrkbz, &
                            bx0,by0,bz0,nx,ny,nz,ngrd, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te)
                    else
                        call bndry_flanks( &
                            wrkqrho,wrkqpx,wrkqpy,wrkqpz, &
                            wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhrho,wrkhpx,wrkhpy,wrkhpz, &
                            wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkorho,wrkopx,wrkopy,wrkopz, &
                            wrkopresx,wrkopresy,wrkopresz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            wrkepres,wrkbx,wrkby,wrkbz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres,bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres,oldbx,oldby,oldbz, vvx, &
                            nx,ny,nz,ngrd,m,t_old,t_new,t_grd, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
                    if(m.le.mbndry)then
                        !       write(6,*)'calling bndry inner lap'
                        call bndry_inner( &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz,rmassq, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz,rmassh, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz,rmasso, &
                            wrkepres,  wrkbx, wrkby, wrkbz, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
                            ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
                            mbndry,msrf,mmid,mzero, &
                            erho,epress,alpha_e,ti_te,o_conc, &
                            sbx_wind,spress, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
                    !
                    if(spacecraft)then
                        call set_imf(wrkbx,wrkby,wrkbz,bx0,by0,bz0,bxp,byp,bzp, &
                            wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                            wrkqpx,wrkqpy,wrkqpz, &
                            wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                            wrkhpx,wrkhpy,wrkhpz, &
                            wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                            wrkopx,wrkopy,wrkopz, &
                            rmassq,rmassh,rmasso,wrkepres, &
                            wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                            wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                            wrkopresxy,wrkopresxz,wrkopresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,ngrd)
                    endif
                    call set_rho(wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz,rmassq, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz,rmassh, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopresxy,wrkopresxz,wrkopresyz,rmasso, &
                        wrkepres,nx,ny,nz,ngrd,m,o_conc)
                    !
                    !      write(6,*)'lapidus speeds'
                    !
                    vlim=0.6*m
                    call set_speed_agrd( &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopx,wrkopy,wrkopz, &
                        wrkepres,wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        wrkbx,wrkby,wrkbz,bx0,by0,bz0, &
                        bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest)
                    !
                    !     write(6,994)
                    ! 994 format(' lapidus done')
                    !
                    !
                    !     .......................................................
                    !     add a little bit of flux correction  :
                    !     ........................................................
                    !
                    !
                    !       write(6,*)'flux correct starting'
                    call flux_correct(qrho,qpresx,qpresy,qpresz, &
                    qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
                        wrkqrho,wrkqpresx,wrkqpresy,wrkqpresz, &
                        wrkqpresxy,wrkqpresxz,wrkqpresyz, &
                        wrkqpx,wrkqpy,wrkqpz, &
                        oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                        oldqpresxy,oldqpresxz,oldqpresyz, &
                        oldqpx,oldqpy,oldqpz, &
                    !
                        hrho,hpresx,hpresy,hpresz, &
                        hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
                        wrkhrho,wrkhpresx,wrkhpresy,wrkhpresz, &
                        wrkhpresxy,wrkhpresxz,wrkhpresyz, &
                        wrkhpx,wrkhpy,wrkhpz, &
                        oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                        oldhpresxy,oldhpresxz,oldhpresyz, &
                        oldhpx,oldhpy,oldhpz, &
                    !
                        orho,opresx,opresy,opresz, &
                        opresxy,opresxz,opresyz,opx,opy,opz, &
                        wrkorho,wrkopresx,wrkopresy,wrkopresz, &
                        wrkopresxy,wrkopresxz,wrkopresyz, &
                        wrkopx,wrkopy,wrkopz, &
                        oldorho,oldopresx,oldopresy,oldopresz, &
                        oldopresxy,oldopresxz,oldopresyz, &
                        oldopx,oldopy,oldopz, &
                    !
                        epres,wrkepres,oldepres, &
                        bx,by,bz,wrkbx,wrkby,wrkbz, &
                        oldbx,oldby,oldbz,vvx,vvy,vvz, &
                        nx,ny,nz,ngrd,m,difrho,diferg,xspac)
                    !
                    !      set bndry conditions
                    !
                    if(m.eq.ngrd) then
                        call bndry_outer( &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rmassq,rmassh,rmasso,bx,by,bz, &
                            bx0,by0,bz0,nx,ny,nz,ngrd, &
                            srho,rho_frac,o_conc,spress,spx,spy,spz, &
                            sbx_wind,sby_wind,sbz_wind,ti_te)
                    else
                        call bndry_flanks( &
                            qrho,qpx,qpy,qpz, &
                            qpresx,qpresy,qpresz,qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz, &
                            hpresx,hpresy,hpresz,hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz, &
                            opresx,opresy,opresz,opresxy,opresxz,opresyz, &
                            epres, &
                            bx,by,bz, &
                            qrho,qpx,qpy,qpz,qpresx,qpresy,qpresz, &
                            qpresxy,qpresxz,qpresyz, &
                            hrho,hpx,hpy,hpz,hpresx,hpresy,hpresz, &
                            hpresxy,hpresxz,hpresyz, &
                            orho,opx,opy,opz,opresx,opresy,opresz, &
                            opresxy,opresxz,opresyz, &
                            epres, &
                            bx,by,bz, &
                            oldqrho,oldqpx,oldqpy,oldqpz, &
                            oldqpresx,oldqpresy,oldqpresz, &
                            oldqpresxy,oldqpresxz,oldqpresyz, &
                            oldhrho,oldhpx,oldhpy,oldhpz, &
                            oldhpresx,oldhpresy,oldhpresz, &
                            oldhpresxy,oldhpresxz,oldhpresyz, &
                            oldorho,oldopx,oldopy,oldopz, &
                            oldopresx,oldopresy,oldopresz, &
                            oldopresxy,oldopresxz,oldopresyz, &
                            oldepres, &
                            oldbx,oldby,oldbz, vvx, &
                            nx,ny,nz,ngrd,m,t_old,t_new,t_grd, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                    endif
           
                    !
                    if(m.le.mbndry)then
                        !       write(6,*)'calling bndry inner'
                        call bndry_inner( &
                            qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
                            orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
                            epres,bx,by,bz, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            nx,ny,nz,ngrd,parm_srf,parm_mid,parm_zero, &
                            ijsrf,numsrf,ijmid,nummid,ijzero,numzero, &
                            mbndry,msrf,mmid,mzero, &
                            erho,epress,alpha_e,ti_te,o_conc, &
                            re_equiv,rearth,sbx_wind,spress, &
                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                            grd_zmin,grd_zmax)
                        !
                    endif   ! end mbndry
           
                    if(spacecraft)then
                        call set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                            hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                            orho,opresx,opresy,opresz,opx,opy,opz, &
                            rmassq,rmassh,rmasso,epres, &
                            qpresxy,qpresxz,qpresyz, &
                            hpresxy,hpresxz,hpresyz, &
                            opresxy,opresxz,opresyz, &
                            rhop,svxp,svyp,svzp,svelx,spress, &
                            ti_te,rho_frac,nx,ny,nz,ngrd)
                    endif
                    !
                    !      write(6,*)'fcsmooth speeds'
                    !
                    call set_rho(qrho,qpresx,qpresy,qpresz, &
                        qpresxy,qpresxz,qpresyz,rmassq, &
                        hrho,hpresx,hpresy,hpresz, &
                        hpresxy,hpresxz,hpresyz,rmassh, &
                        orho,opresx,opresy,opresz, &
                        opresxy,opresxz,opresyz,rmasso, &
                        epres,nx,ny,nz,ngrd,m,o_conc)
                    !
                    vlim=0.6*m
                    call set_speed_agrd( &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        epres,qpresxy,qpresxz,qpresyz, &
                        hpresxy,hpresxz,hpresyz,opresxy,opresxz,opresyz, &
                        bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz,btot, &
                        vvx,tvx,tvy,tvz,evx,evy,evz,curx,cury,curz, &
                        rmassq,rmassh,rmasso,nx,ny,nz,ngrd,m, &
                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                        vlim,alf_lim,o_conc,fastest)
                        t_stepnew(m)=stepsz*xspac(m)/fastest
                    !     write(6,*)'needed step of',t_step(m),t_stepnew(m)
                    !
                    !
                    !     sync time steps if needed and apply core conditions
                    !
                    if(m.lt.ngrd)then
                        if(abs(t_new(m)-t_new(m+1)).le.0.005*t_step(m))then
                            call bndry_corer( &
                                qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                                hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                                orho,opresx,opresy,opresz,opx,opy,opz, &
                                epres,bx,by,bz, &
                                qpresxy,qpresxz,qpresyz, &
                                hpresxy,hpresxz,hpresyz, &
                                opresxy,opresxz,opresyz, &
                                nx,ny,nz,ngrd,m, &
                                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                                grd_zmin,grd_zmax)
                            !         write(6,*)'corer',m,t_new(m),m+1,t_new(m+1)
                            t_old(m+1)=t_new(m+1)
                        endif
                    endif
                    !
                    !
                    !
           
                    !     if (m.eq.0.) then
                    if (m.eq.main_grd_n) then
                        !     ***********************************************
                        !     start hires grid variable time loop over delt/m_step
                        !     **********************************************
                        !
                        !     write(6,*)'main grd loop times',t_old(main_grd_n),
                        !    +                     t_new(main_grd_n)
                        do ms_n=1,m_step_n      !box sweep of hires
                            do m_n=ngrd_n,1,-1      !hires box increment
                                !
                                !     test if grid sector needs to be moved in time
                                !
                                yes_step_n=.false.
                                !
                                !     write(6,*)'hires loop',m_n,t_old_n(m_n),t_new_n(m_n)
                                if(m_n.eq.1)then
                                    yes_step_n=.true.
                                else
                                    if((t_old_n(m_n).eq.t_new_n(m_n)) &
                                        .and.(abs(t_new_n(m_n)-t_new_n(1)).le. &
                                        0.005*t_step_n(1)))then
                                        yes_step_n=.true.
                                    endif
                                endif
                                !
                                !      time step grid
                                !
                                if(yes_step_n) then
                                    t_old_n(m_n)=t_new_n(m_n)
                                    t_new_n(m_n)=t_old_n(m_n)+t_step_n(m_n)
                                    write(6,*)'hires loop times',m_n,t_old_n(m_n),t_new_n(m_n)
                                    !
                                    delt_n= t_step_n(m_n)
                                    delt2_n=delt_n/2.
                                    !
                                    if(tilting)then
                                        atilt=old_tilt+(t_old_n(m_n)+delt2_n)*dtilt
                                        sin_tilt=sin(atilt*.0174533)
                                        cos_tilt=cos(atilt*.0174533)
                                        ut2=utold+(t_old_n(m_n)+delt2_n)*t_equiv/3600.
                                        nrot2=ut2/planet_per
                                        rot_hr2=ut2-nrot2*planet_per
                                        rot_angle=6.2832*rot_hr2/planet_per
                                        !       write(6,*)'mak_dip moon', t_old_n(m_n),delt2_n
                                        call mak_dip_moon(bx0_n,by0_n,bz0_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            mbndry_n,m_n,ijzero_n,numzero_n,mzero_n,rearth, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    !     move the moon
                                    !
                                    !      write(6,*)'fixing to move moon'
                                    if(m_n.le.mbndry_n)then
                                        ut2=utold+(t_old_n(m_n)+delt2_n)*t_equiv/3600.
                                        theta=((ut2-ut_insert)/ut_orbit)*2.*3.1414
                                        xmoon=r_orbit*sin(theta)
                                        ymoon=-r_orbit*cos(theta)
                                        zmoon=0.
                                        vx_moon=v_orbit*cos(theta)
                                        vy_moon=v_orbit*sin(theta)
                                        vz_moon=0.
                                        !      write(6,*)'move moon'
                                        call set_bndry_moon_ram(rmassq,rmassh,rmasso, &
                                            m_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            parm_srf_n,parm_mid_n,parm_zero_n, &
                                            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                                            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                                            qden_moon,hden_moon,oden_moon, &
                                            vx_moon,vy_moon,tempi,gamma, &
                                            ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_e, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    !     **************************************************************
                                    !     sub grid loop with variable time stepping
                                    !     **************************************************************
                                    !
                                    !     store old parameters
                                    !
                                    !      write(6,*)'saving old data'
                                    call store_array(oldqrho_n,qrho_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldqpresx_n,qpresx_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldqpresy_n,qpresy_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldqpresz_n,qpresz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldqpresxy_n,qpresxy_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldqpresxz_n,qpresxz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldqpresyz_n,qpresyz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldqpx_n,qpx_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldqpy_n,qpy_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldqpz_n,qpz_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
           !
                                    call store_array(oldhrho_n,hrho_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldhpresx_n,hpresx_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldhpresy_n,hpresy_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldhpresz_n,hpresz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldhpresxy_n,hpresxy_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldhpresxz_n,hpresxz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldhpresyz_n,hpresyz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldhpx_n,hpx_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldhpy_n,hpy_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldhpz_n,hpz_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
           !
                                    call store_array(oldorho_n,orho_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldopresx_n,opresx_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldopresy_n,opresy_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldopresz_n,opresz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldopresxy_n,opresxy_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldopresxz_n,opresxz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldopresyz_n,opresyz_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldopx_n,opx_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldopy_n,opy_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldopz_n,opz_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    !
                                    call store_array(oldepres_n,epres_n,nx_n,ny_n,nz_n, &
                                        ngrd_n,m_n)
                                    call store_array(oldbx_n,bx_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldby_n,by_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call store_array(oldbz_n,bz_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    !
                                    !      write(6,*)'hires lax1'
                                    !
                                    !     **********************************************************
                                    !     lax step 1  - hires
                                    !     **********************************************************
                                    !
                                    rx=xspac_n(m_n)
                                    ry=xspac_n(m_n)
                                    rz=xspac_n(m_n)
                                    !
                                    call calcur(bx_n,by_n,bz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        curx_n,cury_n,curz_n,rx,ry,rz)
                                    !
                                    !     find total magnetic field
                                    !
                                    call totfld(bx_n,bx0_n,bsx_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call totfld(by_n,by0_n,bsy_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call totfld(bz_n,bz0_n,bsz_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    !
                                    !     find magnitude of b
                                    !
                                    call tot_b(btot_n,bsx_n,bsy_n,bsz_n,nx_n,ny_n,nz_n)
                                    !
                                    !
                                    call fnd_evel(qpx_n,qpy_n,qpz_n,qrho_n, &
                                        hpx_n,hpy_n,hpz_n,hrho_n, &
                                        opx_n,opy_n,opz_n,orho_n, &
                                        curx_n,cury_n,curz_n,evx_n,evy_n,evz_n, &
                                        tvx_n,tvy_n,tvz_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        rmassq,rmassh,rmasso,reynolds)
                                    !
                                    call bande(efldx_n,efldy_n,efldz_n,bsx_n,bsy_n,bsz_n, &
                                        curx_n,cury_n,curz_n,evx_n,evy_n,evz_n,btot_n, &
                                        epres_n,qrho_n,hrho_n,orho_n, &
                                        resistive_n,resist,reynolds, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,rmassq,rmassh,rmasso, &
                                        ijmid_n,nummid_n,ijzero_n,mbndry_n,numzero_n, &
                                        mmid_n,mzero_n,rx,ry,rz)
                                    !
                                    call push_elec(wrkepres_n,oldepres_n,epres_n, &
                                        evx_n,evy_n,evz_n, &
                                        gamma,gamma1,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        0.5*delt_n,rx,ry,rz)
                                    !
                                    call push_ion(wrkqrho_n,wrkqpresx_n,wrkqpresy_n, &
                                        wrkqpresz_n, wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                        oldqrho_n,oldqpresx_n,oldqpresy_n,oldqpresz_n, &
                                        oldqpx_n,oldqpy_n,oldqpz_n, &
                                        qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                        qpx_n,qpy_n,qpz_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                        oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
                                        qpresxy_n,qpresxz_n,qpresyz_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n,efldx_n,efldy_n,efldz_n, &
                                        rmassq, &
                                        vvx_n,vvy_n,vvz_n,tvx_n,tvy_n,tvz_n,gamma,gamma1, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,0.5*delt_n,grav, &
                                        re_equiv,reynolds, &
                                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                        grd_zmin_n,grd_zmax_n,ani_q)
                                    !
                                    call push_ion(wrkhrho_n,wrkhpresx_n,wrkhpresy_n, &
                                        wrkhpresz_n, wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                        oldhrho_n,oldhpresx_n,oldhpresy_n,oldhpresz_n, &
                                        oldhpx_n,oldhpy_n,oldhpz_n, &
                                        hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                        hpx_n,hpy_n,hpz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                        oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n,efldx_n,efldy_n,efldz_n, &
                                        rmassh, &
                                        vvx_n,vvy_n,vvz_n,tvx_n,tvy_n,tvz_n,gamma,gamma1, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,0.5*delt_n,grav, &
                                        re_equiv,reynolds, &
                                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                        grd_zmin_n,grd_zmax_n,ani_h)
                                    !
                                    !
                                    call push_ion(wrkorho_n,wrkopresx_n,wrkopresy_n, &
                                        wrkopresz_n, wrkopx_n,wrkopy_n,wrkopz_n, &
                                        oldorho_n,oldopresx_n,oldopresy_n,oldopresz_n, &
                                        oldopx_n,oldopy_n,oldopz_n, &
                                        orho_n,opresx_n,opresy_n,opresz_n, &
                                        opx_n,opy_n,opz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                        oldopresxy_n,oldopresxz_n,oldopresyz_n, &
                                        opresxy_n,opresxz_n,opresyz_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n,efldx_n,efldy_n,efldz_n, &
                                        rmasso, &
                                        vvx_n,vvy_n,vvz_n,tvx_n,tvy_n,tvz_n,gamma,gamma1, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,0.5*delt_n,grav, &
                                        re_equiv,reynolds, &
                                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                        grd_zmin_n,grd_zmax_n,ani_o)
                                    !
                                    call push_bfld(wrkbx_n,wrkby_n,wrkbz_n, &
                                        oldbx_n,oldby_n,oldbz_n, &
                                        efldx_n,efldy_n,efldz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n, &
                                        m_n,0.5*delt_n,rx,ry,rz)
                                    !
                                    !      write(6,*)'bndry lax 1 hires'
                                    !
                                    !     *************************************************************
                                    !     apply boundary conditions
                                    !     *************************************************************
                                    !
                                    !
                                    !     fine grid bndry conditions
                                    !
                                    !
                                    if(m_n.eq.ngrd_n)then
                                        t_grd_n=t_old_n(m_n)+delt2_n
                                        if(t_grd_n.gt.t_new(main_grd_n))then
                                            t_grd_n=t_new(main_grd_n)
                                        endif
                                        if (t_grd_n.lt.t_old(main_grd_n))then
                                            t_grd_n=t_old(main_grd_n)
                                        endif
                                        call bndry_grds( &
                                            wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                            wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                            wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                            wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                            wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                            wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                            wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                            wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                            wrkopx_n,wrkopy_n,wrkopz_n, &
                                            wrkepres_n,wrkbx_n,wrkby_n,wrkbz_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,main_grd_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n, &
                                            qrho,qpresx,qpresy,qpresz, &
                                            qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
                                            hrho,hpresx,hpresy,hpresz, &
                                            hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
                                            orho,opresx,opresy,opresz, &
                                            opresxy,opresxz,opresyz,opx,opy,opz, &
                                            epres,bx,by,bz, &
                                            oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                                            oldqpresxy,oldqpresxz,oldqpresyz, &
                                            oldqpx,oldqpy,oldqpz, &
                                            oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                                            oldhpresxy,oldhpresxz,oldhpresyz, &
                                            oldhpx,oldhpy,oldhpz, &
                                            oldorho,oldopresx,oldopresy,oldopresz, &
                                            oldopresxy,oldopresxz,oldopresyz, &
                                            oldopx,oldopy,oldopz, &
                                            oldepres,oldbx,oldby,oldbz,vvx, &
                                            nx,ny,nz,ngrd,t_old,t_new, &
                                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                                            grd_zmin,grd_zmax)
                                    else
                                        t_grd_n=t_old_n(m_n)+delt2_n
                                        if(t_grd_n.gt.t_new_n(m_n+1))then
                                            t_grd_n=t_new(m_n+1)
                                        endif
                                        if (t_grd_n.lt.t_old_n(m_n+1))then
                                            t_grd_n=t_old_n(m_n+1)
                                        endif
                                        !      write(6,*)'calling bndry flanks',t_old_n,t_new_n,t_grd_n
                                        call bndry_flanks( &
                                            wrkqrho_n,wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                            wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                            wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                            wrkhrho_n,wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                            wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                            wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                            wrkorho_n,wrkopx_n,wrkopy_n,wrkopz_n, &
                                            wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                            wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                            wrkepres_n,wrkbx_n,wrkby_n,wrkbz_n, &
                                            qrho_n,qpx_n,qpy_n,qpz_n, &
                                            qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            hrho_n,hpx_n,hpy_n,hpz_n, &
                                            hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            orho_n,opx_n,opy_n,opz_n, &
                                            opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            oldqrho_n,oldqpx_n,oldqpy_n,oldqpz_n, &
                                            oldqpresx_n,oldqpresy_n,oldqpresz_n, &
                                            oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
                                            oldhrho_n,oldhpx_n,oldhpy_n,oldhpz_n, &
                                            oldhpresx_n,oldhpresy_n,oldhpresz_n, &
                                            oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
                                            oldorho_n,oldopx_n,oldopy_n,oldopz_n, &
                                            oldopresx_n,oldopresy_n,oldopresz_n, &
                                            oldopresxy_n,oldopresxz_n,oldopresyz_n, &
                                            oldepres_n,oldbx_n,oldby_n,oldbz_n,vvx_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,m_n,t_old_n,t_new_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    if(m_n.le.mbndry_n)then
                                        call bndry_moon(wrkqrho_n, &
                                            wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                            wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                            wrkqpx_n,wrkqpy_n,wrkqpz_n,rmassq, &
                                            wrkhrho_n, &
                                            wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                            wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                            wrkhpx_n,wrkhpy_n,wrkhpz_n,rmassh, &
                                            wrkorho_n, &
                                            wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                            wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                            wrkopx_n,wrkopy_n,wrkopz_n,rmasso, &
                                            wrkepres_n,wrkbx_n,wrkby_n,wrkbz_n,m_n, &
                                            nx_n,ny_n,nz_n,ngrd_n, &
                                            parm_srf_n,parm_mid_n,parm_zero_n, &
                                            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                                            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                                            qden_moon,hden_moon,oden_moon,cs_moon,gamma, &
                                            ti_te_moon,vx_moon,vy_moon,vz_moon, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    call set_rho(wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n,rmassq, &
                                        wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n,rmassh, &
                                        wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n,rmasso, &
                                        wrkepres_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        o_conc)
               
                                    !
                                    vlim=0.6
                                    call set_speed_agrd( &
                                        wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                        wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                        wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                        wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                        wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                        wrkopx_n,wrkopy_n,wrkopz_n,wrkepres_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                        wrkbx_n,wrkby_n,wrkbz_n,bx0_n,by0_n,bz0_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n, &
                                        vvx_n,tvx_n,tvy_n,tvz_n,evx_n,evy_n,evz_n, &
                                        curx_n,cury_n,curz_n, &
                                        rmassq,rmassh,rmasso,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                                        vlim,alf_lim,o_conc,fastest)
                                    !
               
                                    !
                                    !     ***********************************************************
                                    !
                                    !     lax step 2
                                    !
                                    !     *************************************************************
                                    !
                                    rx=xspac_n(m_n)
                                    ry=xspac_n(m_n)
                                    rz=xspac_n(m_n)
                                    !
                                    if(tilting)then
                                        atilt=old_tilt+(t_old_n(m_n)+delt_n)*dtilt
                                        sin_tilt=sin(atilt*.0174533)
                                        cos_tilt=cos(atilt*.0174533)
                                        ut2=utold+(t_old_n(m_n)+delt_n)*t_equiv/3600.
                                        nrot2=ut2/planet_per
                                        rot_hr2=ut2-nrot2*planet_per
                                        rot_angle=6.2832*rot_hr2/planet_per
                                        !       write(6,*)'mak_dip with moon',t_old_n(m_n),delt_n
                                        call mak_dip_moon(bx0_n,by0_n,bz0_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            mbndry_n,m_n,ijzero_n,numzero_n,mzero_n,rearth, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    !     move the moon
                                    !
                                    if(m_n.le.mbndry_n)then
                                        ut2=utold+(t_old_n(m_n)+delt_n)*t_equiv/3600.
                                        theta=((ut2-ut_insert)/ut_orbit)*2.*3.1414
                                        xmoon=r_orbit*sin(theta)
                                        ymoon=-r_orbit*cos(theta)
                                        zmoon=0.
                                        vx_moon=v_orbit*cos(theta)
                                        vy_moon=v_orbit*sin(theta)
                                        vz_moon=0.
                                        !       write(6,*)' set_bndry_moon'
                                        call set_bndry_moon_ram(rmassq,rmassh,rmasso, &
                                            m_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            parm_srf_n,parm_mid_n,parm_zero_n, &
                                            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                                            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                                            qden_moon,hden_moon,oden_moon, &
                                            vx_moon,vy_moon,tempi,gamma, &
                                            ti_te_moon,xmoon,ymoon,zmoon,rmoon,offset,alpha_e, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !      write(6,*)'strt lax 2'
                                    !
                                    !     calculate standard mhd current j = curl b
                                    !
                                    call calcur(wrkbx_n,wrkby_n,wrkbz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        curx_n,cury_n,curz_n,rx,ry,rz)
                                    !
                                    !     find total magnetic field
                                    !
                                    call totfld(wrkbx_n,bx0_n,bsx_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call totfld(wrkby_n,by0_n,bsy_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call totfld(wrkbz_n,bz0_n,bsz_n,nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    !
                                    !     find magnitude of b
                                    !
                                    call tot_b(btot_n,bsx_n,bsy_n,bsz_n,nx_n,ny_n,nz_n)
                                    !
                                    !     find the  electric field from electron momentum eqn
                                    !
                                    call fnd_evel(wrkqpx_n,wrkqpy_n,wrkqpz_n,wrkqrho_n, &
                                        wrkhpx_n,wrkhpy_n,wrkhpz_n,wrkhrho_n, &
                                        wrkopx_n,wrkopy_n,wrkopz_n,wrkorho_n, &
                                        curx_n,cury_n,curz_n,evx_n,evy_n,evz_n, &
                                        tvx_n,tvy_n,tvz_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        rmassq,rmassh,rmasso,reynolds)
                                    !
                                    !
                                    call bande(efldx_n,efldy_n,efldz_n,bsx_n,bsy_n,bsz_n, &
                                        curx_n,cury_n,curz_n,evx_n,evy_n,evz_n,btot_n, &
                                        wrkepres_n,wrkqrho_n,wrkhrho_n,wrkorho_n, &
                                        resistive_n,resist,reynolds, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,rmassq,rmassh,rmasso, &
                                        ijmid_n,nummid_n,ijzero_n,mbndry_n,numzero_n, &
                                        mmid_n,mzero_n,rx,ry,rz)
                                    !
                                    call push_elec(epres_n,oldepres_n,wrkepres_n, &
                                        evx_n,evy_n,evz_n, &
                                    gamma,gamma1,nx_n,ny_n,nz_n,ngrd_n,m_n,delt,rx,ry,rz)
                                    !
                                    !
                                    call push_ion(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                        qpx_n,qpy_n,qpz_n, &
                                        oldqrho_n,oldqpresx_n,oldqpresy_n,oldqpresz_n, &
                                        oldqpx_n,oldqpy_n,oldqpz_n, &
                                        wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                        wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                        qpresxy_n,qpresxz_n,qpresyz_n, &
                                        oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n,efldx_n,efldy_n,efldz_n,rmassq, &
                                        vvx_n,vvy_n,vvz_n,tvx_n,tvy_n,tvz_n,gamma,gamma1, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,delt_n,grav,re_equiv,reynolds, &
                                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                        grd_zmin_n,grd_zmax_n,ani_q)
                                    !
                                    call push_ion(hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                        hpx_n,hpy_n,hpz_n, &
                                        oldhrho_n,oldhpresx_n,oldhpresy_n,oldhpresz_n, &
                                        oldhpx_n,oldhpy_n,oldhpz_n, &
                                        wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                        wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n, &
                                        oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n,efldx_n,efldy_n,efldz_n,rmassh, &
                                        vvx_n,vvy_n,vvz_n,tvx_n,tvy_n,tvz_n,gamma,gamma1, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,delt_n,grav,re_equiv,reynolds, &
                                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                        grd_zmin_n,grd_zmax_n,ani_h)
                                    !
                                    !
                                    call push_ion(orho_n,opresx_n,opresy_n,opresz_n, &
                                        opx_n,opy_n,opz_n, &
                                        oldorho_n,oldopresx_n,oldopresy_n,oldopresz_n, &
                                        oldopx_n,oldopy_n,oldopz_n, &
                                        wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                        wrkopx_n,wrkopy_n,wrkopz_n, &
                                        opresxy_n,opresxz_n,opresyz_n, &
                                        oldopresxy_n,oldopresxz_n,oldopresyz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n,efldx_n,efldy_n,efldz_n,rmasso, &
                                        vvx_n,vvy_n,vvz_n,tvx_n,tvy_n,tvz_n,gamma,gamma1, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,delt_n,grav,re_equiv,reynolds, &
                                        grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                        grd_zmin_n,grd_zmax_n,ani_o)
                                    !
                                    call push_bfld(bx_n,by_n,bz_n,oldbx_n,oldby_n,oldbz_n, &
                                        efldx_n,efldy_n,efldz_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        delt_n,rx,ry,rz)
                                    !      write(6,*)'bndry lax 2 hires'
                                    !
                                    !     fine grid bndry conditions
                                    !
                                    if(m_n.eq.ngrd_n)then
                                        t_grd_n=t_old_n(m_n)+delt_n
                                        if(t_grd_n.gt.t_new(main_grd_n))then
                                            t_grd_n=t_new(main_grd_n)
                                        endif
                                        if (t_grd_n.lt.t_old(main_grd_n))then
                                            t_grd_n=t_old(main_grd_n)
                                        endif
               
                                        call bndry_grds( &
                                            qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            qpx_n,qpy_n,qpz_n, &
                                            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            hpx_n,hpy_n,hpz_n, &
                                            orho_n,opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            opx_n,opy_n,opz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,main_grd_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n, &
                                            qrho,qpresx,qpresy,qpresz, &
                                            qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
                                            hrho,hpresx,hpresy,hpresz, &
                                            hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
                                            orho,opresx,opresy,opresz, &
                                            opresxy,opresxz,opresyz,opx,opy,opz, &
                                            epres,bx,by,bz, &
                                            oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                                            oldqpresxy,oldqpresxz,oldqpresyz, &
                                            oldqpx,oldqpy,oldqpz, &
                                            oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                                            oldhpresxy,oldhpresxz,oldhpresyz, &
                                            oldhpx,oldhpy,oldhpz, &
                                            oldorho,oldopresx,oldopresy,oldopresz, &
                                            oldopresxy,oldopresxz,oldopresyz, &
                                            oldopx,oldopy,oldopz, &
                                            oldepres,oldbx,oldby,oldbz,vvx, &
                                            nx,ny,nz,ngrd,t_old,t_new, &
                                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                                            grd_zmin,grd_zmax)
                                    else
                                        t_grd_n=t_old_n(m_n)+delt_n
                                        if(t_grd_n.gt.t_new_n(m_n+1))then
                                            t_grd_n=t_new(m_n+1)
                                        endif
                                        if (t_grd_n.lt.t_old_n(m_n+1))then
                                            t_grd_n=t_old_n(m_n+1)
                                        endif
                                        call bndry_flanks( &
                                            qrho_n,qpx_n,qpy_n,qpz_n, &
                                            qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            hrho_n,hpx_n,hpy_n,hpz_n, &
                                            hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            orho_n,opx_n,opy_n,opz_n, &
                                            opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            qrho_n,qpx_n,qpy_n,qpz_n, &
                                            qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            hrho_n,hpx_n,hpy_n,hpz_n, &
                                            hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            orho_n,opx_n,opy_n,opz_n, &
                                            opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            oldqrho_n,oldqpx_n,oldqpy_n,oldqpz_n, &
                                            oldqpresx_n,oldqpresy_n,oldqpresz_n, &
                                            oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
                                            oldhrho_n,oldhpx_n,oldhpy_n,oldhpz_n, &
                                            oldhpresx_n,oldhpresy_n,oldhpresz_n, &
                                            oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
                                            oldorho_n,oldopx_n,oldopy_n,oldopz_n, &
                                            oldopresx_n,oldopresy_n,oldopresz_n, &
                                            oldopresxy_n,oldopresxz_n,oldopresyz_n, &
                                            oldepres_n,oldbx_n,oldby_n,oldbz_n,vvx_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,m_n,t_old_n,t_new_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    if(m_n.le.mbndry_n)then
                                        call bndry_moon(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            qpx_n,qpy_n,qpz_n,rmassq, &
                                            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            hpx_n,hpy_n,hpz_n,rmassh, &
                                            orho_n,opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            opx_n,opy_n,opz_n,rmasso, &
                                            epres_n,bx_n,by_n,bz_n,m_n, &
                                            nx_n,ny_n,nz_n,ngrd_n, &
                                            parm_srf_n,parm_mid_n,parm_zero_n, &
                                            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                                            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                                            qden_moon,hden_moon,oden_moon,cs_moon,gamma, &
                                            ti_te_moon,vx_moon,vy_moon,vz_moon, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    !      write(6,*)'calling set_rho'
                                    call set_rho(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                        qpresxy_n,qpresxz_n,qpresyz_n,rmassq, &
                                        hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n,rmassh, &
                                        orho_n,opresx_n,opresy_n,opresz_n, &
                                        opresxy_n,opresxz_n,opresyz_n,rmasso, &
                                        epres_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        o_conc)
                                    !
                                    !      write(6,*)'calling set_speed'
                                    vlim=0.6
                                    call set_speed_agrd( &
                                        qrho_n,qpresx_n,qpresy_n,qpresz_n,qpx_n,qpy_n,qpz_n, &
                                        hrho_n,hpresx_n,hpresy_n,hpresz_n,hpx_n,hpy_n,hpz_n, &
                                        orho_n,opresx_n,opresy_n,opresz_n,opx_n,opy_n,opz_n, &
                                        epres_n,qpresxy_n,qpresxz_n,qpresyz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n, &
                                        opresxy_n,opresxz_n,opresyz_n, &
                                        bx_n,by_n,bz_n,bx0_n,by0_n,bz0_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n, &
                                        vvx_n,tvx_n,tvy_n,tvz_n,evx_n,evy_n,evz_n, &
                                        curx_n,cury_n,curz_n, &
                                        rmassq,rmassh,rmasso,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                                        vlim,alf_lim,o_conc,fastest)
                                    !
                                    !     lapidus smoothing
                                    !
                                    !      species 1
                                    !
                                    !     write(*,*)'lap ion1'
                                    call fnd_vel(qpx_n,qpy_n,qpz_n,qrho_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call lap_plasma(qrho_n,qpx_n,qpy_n,qpz_n, &
                                        qpresx_n,qpresy_n,qpresz_n, &
                                        qpresxy_n,qpresxz_n,qpresyz_n, &
                                        wrkqrho_n,wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                        wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        chirho_n,chipxyz_n,chierg_n,delt_n)
                                    !
                                    !     species 2
                                    !
                                    !     write(*,*)'lap ion2'
                                    call fnd_vel(hpx_n,hpy_n,hpz_n,hrho_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call lap_plasma(hrho_n,hpx_n,hpy_n,hpz_n, &
                                        hpresx_n,hpresy_n,hpresz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n, &
                                        wrkhrho_n,wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                        wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        chirho_n,chipxyz_n,chierg_n,delt_n)
                                    !
                                    !     species 3
                                    !
                                    !     write(*,*)'lap ion3'
                                    call fnd_vel(opx_n,opy_n,opz_n,orho_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n)
                                    call lap_plasma(orho_n,opx_n,opy_n,opz_n, &
                                        opresx_n,opresy_n,opresz_n, &
                                        opresxy_n,opresxz_n,opresyz_n, &
                                        wrkorho_n,wrkopx_n,wrkopy_n,wrkopz_n, &
                                        wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        chirho_n,chipxyz_n,chierg_n,delt_n)
                                    !
                                    !     electrons
                                    !
                                    call fnd_vtot(qpx_n,qpy_n,qpz_n,qrho_n, &
                                        hpx_n,hpy_n,hpz_n,hrho_n, &
                                        opx_n,opy_n,opz_n,orho_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        rmassq,rmassh,rmasso)
                                    call lap_elec(epres_n,wrkepres_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        chirho_n,chipxyz_n,chierg_n,delt_n)
                                    !
                                    !     amgnetic field fields
                                    !
                                    call lap_bfld(bx_n,by_n,bz_n, &
                                        wrkbx_n,wrkby_n,wrkbz_n, &
                                        curx_n,cury_n,curz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        chirho_n,chipxyz_n,chierg_n,delt_n)
               
                                    !      write(6,*)'bndry lapidius'
                                    !
                                    !     fine grid bndry conditions
                                    !
                                    !
                                    if(m_n.eq.ngrd_n)then
                                        t_grd_n=t_old_n(m_n)+delt_n
                                        if(t_grd_n.gt.t_new(main_grd_n))then
                                            t_grd_n=t_new(main_grd_n)
                                        endif
                                        if (t_grd_n.lt.t_old(main_grd_n))then
                                            t_grd_n=t_old(main_grd_n)
                                        endif
                                        call bndry_grds( &
                                            wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                            wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                            wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                            wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                            wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                            wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                            wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                            wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                            wrkopx_n,wrkopy_n,wrkopz_n, &
                                            wrkepres_n,wrkbx_n,wrkby_n,wrkbz_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,main_grd_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n, &
                                            qrho,qpresx,qpresy,qpresz, &
                                            qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
                                            hrho,hpresx,hpresy,hpresz, &
                                            hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
                                            orho,opresx,opresy,opresz, &
                                            opresxy,opresxz,opresyz,opx,opy,opz, &
                                            epres,bx,by,bz, &
                                            oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                                            oldqpresxy,oldqpresxz,oldqpresyz, &
                                            oldqpx,oldqpy,oldqpz, &
                                            oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                                            oldhpresxy,oldhpresxz,oldhpresyz, &
                                            oldhpx,oldhpy,oldhpz, &
                                            oldorho,oldopresx,oldopresy,oldopresz, &
                                            oldopresxy,oldopresxz,oldopresyz, &
                                            oldopx,oldopy,oldopz, &
                                            oldepres,oldbx,oldby,oldbz,vvx, &
                                            nx,ny,nz,ngrd,t_old,t_new, &
                                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                                            grd_zmin,grd_zmax)
                                    else
                                        t_grd_n=t_old_n(m_n)+delt_n
                                        if(t_grd_n.gt.t_new_n(m_n+1))then
                                            t_grd_n=t_new(m_n+1)
                                        endif
                                        if (t_grd_n.lt.t_old_n(m_n+1))then
                                            t_grd_n=t_old_n(m_n+1)
                                        endif
                                        call bndry_flanks( &
                                            wrkqrho_n,wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                            wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                            wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                            wrkhrho_n,wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                            wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                            wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                            wrkorho_n,wrkopx_n,wrkopy_n,wrkopz_n, &
                                            wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                            wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                            wrkepres_n,wrkbx_n,wrkby_n,wrkbz_n, &
                                            qrho_n,qpx_n,qpy_n,qpz_n, &
                                            qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            hrho_n,hpx_n,hpy_n,hpz_n, &
                                            hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            orho_n,opx_n,opy_n,opz_n, &
                                            opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            oldqrho_n,oldqpx_n,oldqpy_n,oldqpz_n, &
                                            oldqpresx_n,oldqpresy_n,oldqpresz_n, &
                                            oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
                                            oldhrho_n,oldhpx_n,oldhpy_n,oldhpz_n, &
                                            oldhpresx_n,oldhpresy_n,oldhpresz_n, &
                                            oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
                                            oldorho_n,oldopx_n,oldopy_n,oldopz_n, &
                                            oldopresx_n,oldopresy_n,oldopresz_n, &
                                            oldopresxy_n,oldopresxz_n,oldopresyz_n, &
                                            oldepres_n,oldbx_n,oldby_n,oldbz_n,vvx_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,m_n,t_old_n,t_new_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
           
                                    if(m_n.le.mbndry_n)then
                                        call bndry_moon(wrkqrho_n, &
                                            wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                            wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                            wrkqpx_n,wrkqpy_n,wrkqpz_n,rmassq, &
                                            wrkhrho_n, &
                                            wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                            wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                            wrkhpx_n,wrkhpy_n,wrkhpz_n,rmassh, &
                                            wrkorho_n, &
                                            wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                            wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                            wrkopx_n,wrkopy_n,wrkopz_n,rmasso, &
                                            wrkepres_n,wrkbx_n,wrkby_n,wrkbz_n,m_n, &
                                            nx_n,ny_n,nz_n,ngrd_n, &
                                            parm_srf_n,parm_mid_n,parm_zero_n, &
                                            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                                            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                                            qden_moon,hden_moon,oden_moon,cs_moon,gamma, &
                                            ti_te_moon,vx_moon,vy_moon,vz_moon, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    call set_rho(wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n,rmassq, &
                                        wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n,rmassh, &
                                        wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n,rmasso, &
                                        wrkepres_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        o_conc)
                                    !
                                    vlim=0.6
                                    call set_speed_agrd( &
                                        wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                        wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                        wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                        wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                        wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                        wrkopx_n,wrkopy_n,wrkopz_n,wrkepres_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                        wrkbx_n,wrkby_n,wrkbz_n,bx0_n,by0_n,bz0_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n, &
                                        vvx_n,tvx_n,tvy_n,tvz_n,evx_n,evy_n,evz_n, &
                                        curx_n,cury_n,curz_n, &
                                        rmassq,rmassh,rmasso,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                                        vlim,alf_lim,o_conc,fastest)
               
                                    !
                                    !    flux correction smoothing
                                    !
                                    !
                                    !       write(6,*)'calling hires flux_correct'
                                    call flux_correct(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                        qpresxy_n,qpresxz_n,qpresyz_n, &
                                        qpx_n,qpy_n,qpz_n, &
                                        wrkqrho_n,wrkqpresx_n,wrkqpresy_n,wrkqpresz_n, &
                                        wrkqpresxy_n,wrkqpresxz_n,wrkqpresyz_n, &
                                        wrkqpx_n,wrkqpy_n,wrkqpz_n, &
                                        oldqrho_n,oldqpresx_n,oldqpresy_n,oldqpresz_n, &
                                        oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
                                        oldqpx_n,oldqpy_n,oldqpz_n, &
       !  
                                        hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n, &
                                        hpx_n,hpy_n,hpz_n, &
                                        wrkhrho_n,wrkhpresx_n,wrkhpresy_n,wrkhpresz_n, &
                                        wrkhpresxy_n,wrkhpresxz_n,wrkhpresyz_n, &
                                        wrkhpx_n,wrkhpy_n,wrkhpz_n, &
                                        oldhrho_n,oldhpresx_n,oldhpresy_n,oldhpresz_n, &
                                        oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
                                        oldhpx_n,oldhpy_n,oldhpz_n, &
       !   
                                        orho_n,opresx_n,opresy_n,opresz_n, &
                                        opresxy_n,opresxz_n,opresyz_n, &
                                        opx_n,opy_n,opz_n, &
                                        wrkorho_n,wrkopresx_n,wrkopresy_n,wrkopresz_n, &
                                        wrkopresxy_n,wrkopresxz_n,wrkopresyz_n, &
                                        wrkopx_n,wrkopy_n,wrkopz_n, &
                                        oldorho_n,oldopresx_n,oldopresy_n,oldopresz_n, &
                                        oldopresxy_n,oldopresxz_n,oldopresyz_n, &
                                        oldopx_n,oldopy_n,oldopz_n, &
       !    
                                        epres_n,wrkepres_n,oldepres_n, &
                                        bx_n,by_n,bz_n,wrkbx_n,wrkby_n,wrkbz_n, &
                                        oldbx_n,oldby_n,oldbz_n,vvx_n,vvy_n,vvz_n, &
                                        nx_n,ny_n,nz_n,ngrd_n,m_n,difrho_n,diferg_n, &
                                        xspac_n)
                                    !
                                    !      write(6,*)'hires flux-correct complete'
                                    !
                                    !     smooth epres to remove whislter wave turbulence
                                    !
                                    if(m_n.le.2)then
                                        !       write(6,*)'ssmooth'
                                        call ssmooth(wrkepres_n,epres_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            m_n,0.0025)
                                    endif
                                    if(m_n.le.3)then
                                        call ssmooth(wrkopresx_n,opresx_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            m_n,0.001)
                                        call ssmooth(wrkopresy_n,opresy_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            m_n,0.001)
                                        call ssmooth(wrkopresz_n,opresz_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            m_n,0.001)
                                        call ssmooth(wrkopresxy_n,opresxy_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            m_n,0.001)
                                        call ssmooth(wrkopresxz_n,opresxz_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            m_n,0.001)
                                        call ssmooth(wrkopresyz_n,opresyz_n,nx_n,ny_n,nz_n,ngrd_n, &
                                            m_n,0.001)
                                    endif
               
                                    !
                                    !     fine grid bndry conditions
                                    !
                                    !     write(6,*)'bndry flux_correct'
                                    if(m_n.eq.ngrd_n)then
                                        t_grd_n=t_old_n(m_n)+delt_n
                                        if(t_grd_n.gt.t_new(main_grd_n))then
                                            t_grd_n=t_new(main_grd_n)
                                        endif
                                        if (t_grd_n.lt.t_old(main_grd_n))then
                                            t_grd_n=t_old(main_grd_n)
                                        endif
                                        call bndry_grds( &
                                            qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            qpx_n,qpy_n,qpz_n, &
                                            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            hpx_n,hpy_n,hpz_n, &
                                            orho_n,opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            opx_n,opy_n,opz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,main_grd_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n, &
                                            qrho,qpresx,qpresy,qpresz, &
                                            qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
                                            hrho,hpresx,hpresy,hpresz, &
                                            hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
                                            orho,opresx,opresy,opresz, &
                                            opresxy,opresxz,opresyz,opx,opy,opz, &
                                            epres,bx,by,bz, &
                                            oldqrho,oldqpresx,oldqpresy,oldqpresz, &
                                            oldqpresxy,oldqpresxz,oldqpresyz, &
                                            oldqpx,oldqpy,oldqpz, &
                                            oldhrho,oldhpresx,oldhpresy,oldhpresz, &
                                            oldhpresxy,oldhpresxz,oldhpresyz, &
                                            oldhpx,oldhpy,oldhpz, &
                                            oldorho,oldopresx,oldopresy,oldopresz, &
                                            oldopresxy,oldopresxz,oldopresyz, &
                                            oldopx,oldopy,oldopz, &
                                            oldepres,oldbx,oldby,oldbz,vvx, &
                                            nx,ny,nz,ngrd,t_old,t_new, &
                                            grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                                            grd_zmin,grd_zmax)
                                    else
                                        t_grd_n=t_old_n(m_n)+delt_n
                                        if(t_grd_n.gt.t_new_n(m_n+1))then
                                            t_grd_n=t_new(m_n+1)
                                        endif
                                        if (t_grd_n.lt.t_old_n(m_n+1))then
                                            t_grd_n=t_old_n(m_n+1)
                                        endif
                                        call bndry_flanks( &
                                            qrho_n,qpx_n,qpy_n,qpz_n, &
                                            qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            hrho_n,hpx_n,hpy_n,hpz_n, &
                                            hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            orho_n,opx_n,opy_n,opz_n, &
                                            opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            qrho_n,qpx_n,qpy_n,qpz_n, &
                                            qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            hrho_n,hpx_n,hpy_n,hpz_n, &
                                            hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            orho_n,opx_n,opy_n,opz_n, &
                                            opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            epres_n,bx_n,by_n,bz_n, &
                                            oldqrho_n,oldqpx_n,oldqpy_n,oldqpz_n, &
                                            oldqpresx_n,oldqpresy_n,oldqpresz_n, &
                                            oldqpresxy_n,oldqpresxz_n,oldqpresyz_n, &
                                            oldhrho_n,oldhpx_n,oldhpy_n,oldhpz_n, &
                                            oldhpresx_n,oldhpresy_n,oldhpresz_n, &
                                            oldhpresxy_n,oldhpresxz_n,oldhpresyz_n, &
                                            oldorho_n,oldopx_n,oldopy_n,oldopz_n, &
                                            oldopresx_n,oldopresy_n,oldopresz_n, &
                                            oldopresxy_n,oldopresxz_n,oldopresyz_n, &
                                            oldepres_n,oldbx_n,oldby_n,oldbz_n,vvx_n, &
                                            nx_n,ny_n,nz_n,ngrd_n,m_n,t_old_n,t_new_n,t_grd_n, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    if(m_n.le.mbndry_n)then
                                        call bndry_moon(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                            qpresxy_n,qpresxz_n,qpresyz_n, &
                                            qpx_n,qpy_n,qpz_n,rmassq, &
                                            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                            hpresxy_n,hpresxz_n,hpresyz_n, &
                                            hpx_n,hpy_n,hpz_n,rmassh, &
                                            orho_n,opresx_n,opresy_n,opresz_n, &
                                            opresxy_n,opresxz_n,opresyz_n, &
                                            opx_n,opy_n,opz_n,rmasso, &
                                            epres_n,bx_n,by_n,bz_n,m_n, &
                                            nx_n,ny_n,nz_n,ngrd_n, &
                                            parm_srf_n,parm_mid_n,parm_zero_n, &
                                            ijsrf_n,numsrf_n,ijmid_n,nummid_n,ijzero_n, &
                                            numzero_n,mbndry_n,msrf_n,mmid_n,mzero_n, &
                                            qden_moon,hden_moon,oden_moon,cs_moon,gamma, &
                                            ti_te_moon,vx_moon,vy_moon,vz_moon, &
                                            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                            grd_zmin_n,grd_zmax_n)
                                    endif
                                    !
                                    call set_rho(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                        qpresxy_n,qpresxz_n,qpresyz_n,rmassq, &
                                        hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n,rmassh, &
                                        orho_n,opresx_n,opresy_n,opresz_n, &
                                        opresxy_n,opresxz_n,opresyz_n,rmasso, &
                                        epres_n,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        o_conc)
                                    !
                                    vlim=0.6
                                    call set_speed_agrd( &
                                        qrho_n,qpresx_n,qpresy_n,qpresz_n,qpx_n,qpy_n,qpz_n, &
                                        hrho_n,hpresx_n,hpresy_n,hpresz_n,hpx_n,hpy_n,hpz_n, &
                                        orho_n,opresx_n,opresy_n,opresz_n,opx_n,opy_n,opz_n, &
                                        epres_n,qpresxy_n,qpresxz_n,qpresyz_n, &
                                        hpresxy_n,hpresxz_n,hpresyz_n, &
                                        opresxy_n,opresxz_n,opresyz_n, &
                                        bx_n,by_n,bz_n,bx0_n,by0_n,bz0_n, &
                                        bsx_n,bsy_n,bsz_n,btot_n, &
                                        vvx_n,tvx_n,tvy_n,tvz_n,evx_n,evy_n,evz_n, &
                                        curx_n,cury_n,curz_n, &
                                        rmassq,rmassh,rmasso,nx_n,ny_n,nz_n,ngrd_n,m_n, &
                                        pxmax,pymax,pzmax,pmax,csmax,alfmax,gamma, &
                                        vlim,alf_lim,o_conc,fastest)
                                    !
                                    !       check for divb errors on high resoultion grid
                                    !
                                    !        call divb_moon(bx_n,by_n,bz_n,qpx_n,qpy_n,qpz_n,qrho_n,
                                    !    +              qpres_n,nx_n,ny_n,nz_n,ngrd_n,m_n,srho,rearth,
                                    !    +              grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n,
                                    !    +              grd_zmin_n,grd_zmax_n)
                                    !
                                    !
                                    !     apply internal/core boundary conditions if needed
                                    !
                                    if(m_n.eq.ngrd_n)then
                                        dift=abs(t_new_n(m_n)-t_new(main_grd_n))
                                        errort=0.003*t_step_n(m_n)
                                        !       write(6,*)'checking bndry_grd_core',dift,errort
                                        if(dift.le.errort)then
                                            write(6,*)'calling bndry_grd_core'
                                            call bndry_grd_core(qrho,qpresx,qpresy,qpresz, &
                                                qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
                                                hrho,hpresx,hpresy,hpresz, &
                                                hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
                                                orho,opresx,opresy,opresz, &
                                                opresxy,opresxz,opresyz,opx,opy,opz, &
                                                epres,bx,by,bz,nx,ny,nz,ngrd, &
                                                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
                                                qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                                qpresxy_n,qpresxz_n,qpresyz_n,qpx_n,qpy_n,qpz_n, &
                                                hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                                hpresxy_n,hpresxz_n,hpresyz_n,hpx_n,hpy_n,hpz_n, &
                                                orho_n,opresx_n,opresy_n,opresz_n, &
                                                opresxy_n,opresxz_n,opresyz_n,opx_n,opy_n,opz_n, &
                                                epres_n,bx_n,by_n,bz_n, &
                                                nx_n,ny_n,nz_n,ngrd_n,main_grd_n, &
                                                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                                                grd_zmin_n,grd_zmax_n)
                                        endif
                                    else
                                        dift=abs(t_new_n(m_n)-t_new_n(m_n+1))
                                        errort=0.003*t_step_n(m_n)
                                        !       write(6,*)'checking bndry_corer',dift,errort
                                        if(dift.le.errort)then
                                            call bndry_corer(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                                                qpx_n,qpy_n,qpz_n, &
                                                hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                                                hpx_n,hpy_n,hpz_n, &
                                                orho_n,opresx_n,opresy_n,opresz_n, &
                                                opx_n,opy_n,opz_n, &
                                                epres_n,bx_n,by_n,bz_n, &
                                                qpresxy_n,qpresxz_n,qpresyz_n, &
                                                hpresxy_n,hpresxz_n,hpresyz_n, &
                                                opresxy_n,opresxz_n,opresyz_n, &
                                                nx_n,ny_n,nz_n,ngrd_n,m_n,grd_xmin_n,grd_xmax_n, &
                                                grd_ymin_n,grd_ymax_n,grd_zmin_n,grd_zmax_n)
                                            !
                                            !         this tells code that all has been completed
                                            t_old_n(m_n+1)=t_new_n(m_n+1)
                                            !         write(6,*)'corer',m_n,t_new_n(m_n),t_old_n(m_n+1)
                                        endif
                                    endif
                                    !
                                    !
                                endif   !yes_step of hires grid
                                !
                            enddo    !hires box increment
                        enddo    !box sweep of hires grid
                        !
                        !
                    endif  ! end of test for hires grid time step                        !
                    endif   !yes_step of main gri
        
            enddo   ! lores box increment
        enddo    ! box sweep of lores grid            
        !
        !
        !    final sync on boundary conditions
        !
        !     write(6,*)'time sync',t,t_new
        t_new_n(1)=t_new(1)
        t_old_n(1)=t_new_n(1)
        do m_n=1,ngrd_n-1
            !         write(6,*)'m_n',m_n
            call bndry_corer(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                qpx_n,qpy_n,qpz_n, &
                hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                hpx_n,hpy_n,hpz_n, &
                orho_n,opresx_n,opresy_n,opresz_n, &
                opx_n,opy_n,opz_n, &
                epres_n,bx_n,by_n,bz_n, &
                qpresxy_n,qpresxz_n,qpresyz_n, &
                hpresxy_n,hpresxz_n,hpresyz_n, &
                opresxy_n,opresxz_n,opresyz_n, &
                nx_n,ny_n,nz_n,ngrd_n,m_n, &
                grd_xmin_n,grd_xmax_n, &
                grd_ymin_n,grd_ymax_n,grd_zmin_n,grd_zmax_n)
            !         write(6,*)m_n-1,t_old_n(m_n-1),t_new_n(m_n-1)
            t_old_n(m_n)=t_old_n(m_n-1)
            t_new_n(m_n)=t_new_n(m_n-1)
        enddo
        !
        !      write(6,*)'final bndry_grd_core'
        call bndry_grd_core(qrho,qpresx,qpresy,qpresz, &
            qpresxy,qpresxz,qpresyz,qpx,qpy,qpz, &
            hrho,hpresx,hpresy,hpresz, &
            hpresxy,hpresxz,hpresyz,hpx,hpy,hpz, &
            orho,opresx,opresy,opresz, &
            opresxy,opresxz,opresyz,opx,opy,opz, &
            epres,bx,by,bz,nx,ny,nz,ngrd, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax, &
            qrho_n,qpresx_n,qpresy_n,qpresz_n, &
            qpresxy_n,qpresxz_n,qpresyz_n,qpx_n,qpy_n,qpz_n, &
            hrho_n,hpresx_n,hpresy_n,hpresz_n, &
            hpresxy_n,hpresxz_n,hpresyz_n,hpx_n,hpy_n,hpz_n, &
            orho_n,opresx_n,opresy_n,opresz_n, &
            opresxy_n,opresxz_n,opresyz_n,opx_n,opy_n,opz_n, &
            epres_n,bx_n,by_n,bz_n, &
            nx_n,ny_n,nz_n,ngrd_n,main_grd_n, &
            grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
            grd_zmin_n,grd_zmax_n)
        !
        t_old(1)=t_new(1)
        !
        !     write(6,*)'final bndry__core'
        do m=1,ngrd-1
            call bndry_corer( &
                qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                orho,opresx,opresy,opresz,opx,opy,opz, &
                epres,bx,by,bz, &
                qpresxy,qpresxz,qpresyz, &
                hpresxy,hpresxz,hpresyz, &
                opresxy,opresxz,opresyz, &
                nx,ny,nz,ngrd,m, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax)
            !
            call set_rho(qrho,qpresx,qpresy,qpresz, &
                qpresxy,qpresxz,qpresyz,rmassq, &
                hrho,hpresx,hpresy,hpresz, &
                hpresxy,hpresxz,hpresyz,rmassh, &
                orho,opresx,opresy,opresz, &
                opresxy,opresxz,opresyz,rmasso, &
                epres,nx,ny,nz,ngrd,m,o_conc)
                t_old(m)=t_old(m-1)
                t_new(m)=t_new(m-1)
            !
        enddo
      
        !     write(6,*)'end loop 999'
        !
        !     write(6,999)t
        ! 999 format(' step 2 complete at t= ',1pe12.5)
        !
        if(t.ge.tgraf) then
            !
            !     plot plasma propeties
            !
            !      calculate size of plotting stuff and ensure no distortions
            !         over desired scales
            !
            519  write(6,*)'graphics plotted'
            call visual(qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
                hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
                orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
                epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
                curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
                tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
                nx,ny,nz,ngrd,xspac, &
                cross,along,flat,xcraft,ncraft,re_equiv, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
                grd_zmin,grd_zmax,ut,b_equiv,ti_te,rho_equiv)
            !
            call visual_hires(qrho_n,qpresx_n,qpresy_n,qpresz_n, &
                qpx_n,qpy_n,qpz_n,rmassq, &
                hrho_n,hpresx_n,hpresy_n,hpresz_n, &
                hpx_n,hpy_n,hpz_n,rmassh, &
                orho_n,opresx_n,opresy_n,opresz_n, &
                opx_n,opy_n,opz_n,rmasso, &
                epres_n,bx_n,by_n,bz_n,bx0_n,by0_n,bz0_n, &
                bsx_n,bsy_n,bsz_n, &
                curx_n,cury_n,curz_n,efldx_n,efldy_n,efldz_n, &
                tvx_n,tvy_n,tvz_n, &
                tx_n,ty_n,tz_n,tg1_n,tg2_n,tt_n,work_n, &
                mx_n,my_n,mz_n,mz2_n,muvwp2_n, &
                nx_n,ny_n,nz_n,ngrd_n,xspac_n, &
                cross_n,along_n,flat_n,xcraft,ncraft,re_equiv, &
                grd_xmin_n,grd_xmax_n,grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n,ut,b_equiv, &
                ti_te_moon,rho_equiv)
            !
            tgraf=tgraf+deltg
        endif !if(t.ge.tgraf)
      
        if(t.ge.ts1) then
            if(nchf.eq.11) &
                open(11,file='fluid11',status='unknown',form='unformatted')
            if(nchf.eq.12) &
                open(12,file='fluid12',status='unknown',form='unformatted')
            if(nchf.eq.13) &
                open(13,file='fluid13',status='unknown',form='unformatted')
            if(nchf.eq.14) &
                open(14,file='fluid14',status='unknown',form='unformatted')
            if(nchf.eq.15) &
                open(15,file='fluid15',status='unknown',form='unformatted')
            !
            !     write restart data
            !
            write(nchf)t
            write(nchf)qrho
            write(nchf)qpx
            write(nchf)qpy
            write(nchf)qpz
            write(nchf)qpresx
            write(nchf)qpresy
            write(nchf)qpresz
            write(nchf)qpresxy
            write(nchf)qpresxz
            write(nchf)qpresyz
            write(nchf)hrho
            write(nchf)hpx
            write(nchf)hpy
            write(nchf)hpz
            write(nchf)hpresx
            write(nchf)hpresy
            write(nchf)hpresz
            write(nchf)hpresxy
            write(nchf)hpresxz
            write(nchf)hpresyz
            write(nchf)orho
            write(nchf)opx
            write(nchf)opy
            write(nchf)opz
            write(nchf)opresx
            write(nchf)opresy
            write(nchf)opresz
            write(nchf)opresxy
            write(nchf)opresxz
            write(nchf)opresyz
            write(nchf)bx
            write(nchf)by
            write(nchf)bz
            write(nchf)epres
            write(nchf)bx0
            write(nchf)by0
            write(nchf)bz0
            !     write(nchf)qrho0
            !     write(nchf)hrho0
            !     write(nchf)orho0
            !     write(nchf)qpres0
            !     write(nchf)hpres0
            !     write(nchf)opres0
            !     write(nchf)epres0
            write(nchf)parm_srf,parm_mid,parm_zero, &
                ijzero,numzero,ijmid,nummid,ijsrf,numsrf
            !     write(nchf)abx0,aby0,abz0
            !     write(nchf)apx0,apy0,apz0
            !     write(nchf)aerg0
            close(nchf)
            !
            mchf=nchf+10
            if(mchf.eq.21) &
                open(21,file='fluid21',status='unknown',form='unformatted')
            if(mchf.eq.22) &
                open(22,file='fluid22',status='unknown',form='unformatted')
            if(mchf.eq.23) &
                open(23,file='fluid23',status='unknown',form='unformatted')
            if(mchf.eq.24) &
                open(24,file='fluid24',status='unknown',form='unformatted')
            if(mchf.eq.25) &
                open(25,file='fluid25',status='unknown',form='unformatted')
            !
            !     write restart data
            !
            write(mchf)t,ut_insert
            write(mchf)qrho_n
            write(mchf)qpx_n
            write(mchf)qpy_n
            write(mchf)qpz_n
            write(mchf)qpresx_n
            write(mchf)qpresy_n
            write(mchf)qpresz_n
            write(mchf)qpresxy_n
            write(mchf)qpresxz_n
            write(mchf)qpresyz_n
            write(mchf)hrho_n
            write(mchf)hpx_n
            write(mchf)hpy_n
            write(mchf)hpz_n
            write(mchf)hpresx_n
            write(mchf)hpresy_n
            write(mchf)hpresz_n
            write(mchf)hpresxy_n
            write(mchf)hpresxz_n
            write(mchf)hpresyz_n
            write(mchf)orho_n
            write(mchf)opx_n
            write(mchf)opy_n
            write(mchf)opz_n
            write(mchf)opresx_n
            write(mchf)opresy_n
            write(mchf)opresz_n
            write(mchf)opresxy_n
            write(mchf)opresxz_n
            write(mchf)opresyz_n
            write(mchf)bx_n
            write(mchf)by_n
            write(mchf)bz_n
            write(mchf)epres_n
            write(mchf)bx0_n
            write(mchf)by0_n
            write(mchf)bz0_n
            write(mchf)parm_srf_n,parm_mid_n,parm_zero_n, &
                ijzero_n,numzero_n,ijmid_n,nummid_n,ijsrf_n,numsrf_n
            write(mchf)grd_xmin_n,grd_xmax_n, &
                grd_ymin_n,grd_ymax_n, &
                grd_zmin_n,grd_zmax_n, &
                grd_time_n,grd_dx_n,grd_dy_n,grd_vx_n,grd_vy_n
            close(mchf)
      
            !     nchf=23-nchf
            nchf=nchf+1
            if(nchf.gt.15)nchf=11
            ts1=ts1+tsave
            !
            !
            !      check for div b errors
            !
            if(divb_hires)then
                range=1.33*lunar_dist/re_equiv
                write(6,*)'range for divb taper',lunar_dist,range
                do m=ngrd,1,-1
                    write(6,*)'divb on box',m
                    call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                        curx,cury,curz,efldx,efldy,efldz, &
                        7,nx*ny*nz,nx,ny,nz,ngrd,m,xspac)
                    call divb_correct(bx,by,bz,bsx,bsy,bsz,btot, &
                        curx,cury,curz,efldx,efldy,efldz, &
                        7,nx*ny*nz,nx,ny,nz,ngrd,m,xspac)
                    write(6,*)'completed divb on box',m
                    !
                    !
                enddo   !end m loop
                !
                !        apply bndry conditions
                !
                do m=ngrd-1,1,-1
                    call flanks_synced(bx,nx,ny,nz,ngrd,m, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                    call flanks_synced(by,nx,ny,nz,ngrd,m, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                    call flanks_synced(bz,nx,ny,nz,ngrd,m, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                enddo
                do m=1,ngrd-1
                    call bndry_corer( &
                        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
                        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
                        orho,opresx,opresy,opresz,opx,opy,opz, &
                        epres,bx,by,bz, &
                        qpresxy,qpresxz,qpresyz, &
                        hpresxy,hpresxz,hpresyz, &
                        opresxy,opresxz,opresyz, &
                        nx,ny,nz,ngrd,m, &
                        grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
                enddo !end bndry_corer
            endif  ! end divb_lores
            !
            if(ringo)then
                !
                !          check for io - thermal mach number set at 1/4
                !
                m=1
                dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
                dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
                dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
                !
                tot_o=0.
                tot_h=0.
                tot_q=0.
                !
                do k=1,nz
                    az=(grd_zmin(m)+dz*(k-1)-zdip)
                    do j=1,ny
                        ay=grd_ymin(m)+dy*(j-1)-ydip
                        do i=1,nx
                            ax=(grd_xmin(m)+dx*(i-1)-xdip)
                            !
                            rx=ax*re_equiv
                            ry=ay*re_equiv
                            rd=sqrt(rx**2+ry**2)
                            !
                            rvy=rx*v_rot
                            rvx=-ry*v_rot
                            rvz=0.
                            corotate=sqrt(rvx**2+rvy**2)
                            !
                            ar_io=sqrt(ax**2+ay**2)*re_equiv
                            dr_io=abs(ar_io -lunar_dist)
      
                            if(abs(dr_io.lt.2.*torus_rad)) then
                                rscale=exp(-((dr_io)/(0.5*torus_rad))**2) ! scale height in rjs
                                zscale=exp(-((az*re_equiv)/(0.5*torus_rad))**2)
                                !
                                oden=den_lunar*rmasso*rscale*zscale
                                !            write(6,*)i,j,k,m,oden
                                orho(i,j,k,m)=orho(i,j,k,m) +oden
                                del_op=(oden/rmasso)*(corotate**2)*t_torus!temp goes as v**2
                                opresx(i,j,k,m)=opresx(i,j,k,m)+del_op
                                opresy(i,j,k,m)=opresy(i,j,k,m)+del_op
                                opresz(i,j,k,m)=opresz(i,j,k,m)+del_op*aniso_factor
                                opx(i,j,k,m)=opx(i,j,k,m)+reduct*oden*rvx
                                opy(i,j,k,m)=opy(i,j,k,m)+reduct*oden*rvy
                                !
                                hden=oden*rmassh/rmasso/o_conc
                                hrho(i,j,k,m)=hrho(i,j,k,m) +hden
                                del_hp=(hden/rmassh)*(corotate**2)*t_torus  !temp goes as v**2
                                hpresx(i,j,k,m)=hpresx(i,j,k,m)+del_hp
                                hpresy(i,j,k,m)=hpresy(i,j,k,m)+del_hp
                                hpresz(i,j,k,m)=hpresz(i,j,k,m)+del_hp*aniso_factor
                                hpx(i,j,k,m)=hpx(i,j,k,m)+reduct*hden*rvx
                                hpy(i,j,k,m)=hpy(i,j,k,m)+reduct*hden*rvy
                                !
                                qden=0.5*hden*rmassq/rmassh
                                qrho(i,j,k,m)=qrho(i,j,k,m) +qden
                                del_qp=(qden/rmassq)*(corotate**2)*t_torus   !temp goes as v**2
                                qpresx(i,j,k,m)=qpresx(i,j,k,m)+del_qp
                                qpresy(i,j,k,m)=qpresy(i,j,k,m)+del_qp
                                qpresz(i,j,k,m)=qpresz(i,j,k,m)+del_qp*aniso_factor
                                qpx(i,j,k,m)=qpx(i,j,k,m)+reduct*qden*rvx
                                qpy(i,j,k,m)=qpy(i,j,k,m)+reduct*qden*rvy
                                !
                                epres(i,j,k,m)=epres(i,j,k,m)+(del_op+del_hp+del_qp) ! equal temps
                                !
                                tot_q=tot_q+qden*dx*dy*dz
                                tot_h=tot_h+hden*dx*dy*dz
                                tot_o=tot_o+oden*dx*dy*dz
                                !
                            endif
                        enddo
                    enddo
                enddo
                !
                !     scale factors to kg/s
                !
                volume=(re_equiv*planet_rad*1.e3)**3  !(cubic meters)
                atime=tsave*t_equiv
                injections=tstep/tsave
                proton_mass=1.67e-27
                write(6,*)'volume,t_equiv,atime',ut,volume,t_equiv,atime
                !
                tot_q=tot_q*volume/atime*rho_equiv*1.e6/rmassq
                tot_h=tot_h*volume/atime*rho_equiv*1.e6/rmassh
                tot_o=tot_o*volume/atime*rho_equiv*1.e6/rmasso
                write(6,*)'tot torus ions/s',ut,tot_q,tot_h,tot_o
                write(10,*)'tot torus ions/s',ut,tot_q,tot_h,tot_o
                !
                tot_q=tot_q*proton_mass*rmassq
                tot_h=tot_h*proton_mass*rmassh
                tot_o=tot_o*proton_mass*rmasso
                write(6,*)'injections',injections, '  single at'
                write(6,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
                write(10,*)'tot torus kg/s',ut,tot_q,tot_h,tot_o
                !
                !
                !
                !          add into high res grid
                !
                do m=ngrd_n,1,-1
                    dx=(grd_xmax_n(m)-grd_xmin_n(m))/(nx_n-1.)
                    dy=(grd_ymax_n(m)-grd_ymin_n(m))/(ny_n-1.)
                    dz=(grd_zmax_n(m)-grd_zmin_n(m))/(nz_n-1.)
                    !
                    tot_o=0.
                    tot_h=0.
                    tot_q=0.
                    !
                    do k=1,nz_n
                        az=grd_zmin_n(m)+dz*(k-1)
                        do j=1,ny_n
                            ay=grd_ymin_n(m)+dy*(j-1)
                            do i=1,nx_n
                                ax=grd_xmin_n(m)+dx*(i-1)
                                !
                                rx=ax*re_equiv
                                ry=ay*re_equiv
                                rd=sqrt(rx**2+ry**2)
                                !
                                rvy=rx*v_rot
                                rvx=-ry*v_rot
                                rvz=0.
                                corotate=sqrt(rvx**2+rvy**2)
                                !
                                ar=sqrt(ax**2+ay**2)
                                if((ar.lt.1.66*rearth).and. &
                                    (abs(az).le.0.5*rearth) )then
                                    qpx_n(i,j,k,m)=rvx*qrho_n(i,j,k,m)
                                    qpy_n(i,j,k,m)=rvy*qrho_n(i,j,k,m)
                                    hpx_n(i,j,k,m)=rvx*hrho_n(i,j,k,m)
                                    hpy_n(i,j,k,m)=rvy*hrho_n(i,j,k,m)
                                    opx_n(i,j,k,m)=rvx*orho_n(i,j,k,m)
                                    opy_n(i,j,k,m)=rvy*orho_n(i,j,k,m)
                                endif
                                !
                                !
                                ar_io=sqrt(ax**2+ay**2)*re_equiv
                                dr_io=abs(ar_io -lunar_dist)
                                !
                                if(abs(dr_io.lt.2.*torus_rad)) then
                                    rscale=exp(-((dr_io)/(0.5*torus_rad))**2) ! scale height in rjs
                                    zscale=exp(-((az*re_equiv)/(0.5*torus_rad))**2)
                                    !
                                    oden=den_lunar*rmasso*rscale*zscale
                                    !            write(6,*)i,j,k,m,oden
                                    orho_n(i,j,k,m)=orho_n(i,j,k,m) +oden
                                    del_op=(oden/rmasso)*(corotate**2)*t_torus!temp goes as v**2
                                    opresx_n(i,j,k,m)=opresx_n(i,j,k,m)+del_op
                                    opresy_n(i,j,k,m)=opresy_n(i,j,k,m)+del_op
                                    opresz_n(i,j,k,m)=opresz_n(i,j,k,m)+del_op*aniso_factor
                                    opx_n(i,j,k,m)=opx_n(i,j,k,m)+reduct*oden*rvx
                                    opy_n(i,j,k,m)=opy_n(i,j,k,m)+reduct*oden*rvy
                                    !
                                    hden=oden*rmassh/rmasso/o_conc
                                    hrho_n(i,j,k,m)=hrho_n(i,j,k,m) +hden
                                    del_hp=(hden/rmassh)*(corotate**2)*t_torus  !temp goes as v**2
                                    hpresx_n(i,j,k,m)=hpresx_n(i,j,k,m)+del_hp
                                    hpresy_n(i,j,k,m)=hpresy_n(i,j,k,m)+del_hp
                                    hpresz_n(i,j,k,m)=hpresz_n(i,j,k,m)+del_hp*aniso_factor
                                    hpx_n(i,j,k,m)=hpx_n(i,j,k,m)+reduct*hden*rvx
                                    hpy_n(i,j,k,m)=hpy_n(i,j,k,m)+reduct*hden*rvy
                                    !
                                    qden=0.5*hden*rmassq/rmassh
                                    qrho_n(i,j,k,m)=qrho_n(i,j,k,m) +qden
                                    del_qp=(qden/rmassq)*(corotate**2)*t_torus    !temp goes as v**2
                                    qpresx_n(i,j,k,m)=qpresx_n(i,j,k,m)+del_qp
                                    qpresy_n(i,j,k,m)=qpresy_n(i,j,k,m)+del_qp
                                    qpresz_n(i,j,k,m)=qpresz_n(i,j,k,m)+del_qp*aniso_factor
                                    qpx_n(i,j,k,m)=qpx_n(i,j,k,m)+reduct*qden*rvx
                                    qpy_n(i,j,k,m)=qpy_n(i,j,k,m)+reduct*qden*rvy
                                    !
                                    epres_n(i,j,k,m)=epres_n(i,j,k,m) &
                                        +(del_op+del_hp+del_qp) ! equal temps
                                    !
                                    tot_q=tot_q+qden*dx*dy*dz
                                    tot_h=tot_h+hden*dx*dy*dz
                                    tot_o=tot_o+oden*dx*dy*dz
                                    !
                                endif
                            enddo
                        enddo
                    enddo
                    !
                    !     scale factors to kg/s
                    !
                    volume=(re_equiv*planet_rad*1.e3)**3  !(cubic meters)
                    atime=tsave*t_equiv
                    proton_mass=1.67e-27
                    !      write(6,*)'volume,t_equiv,atime',volume,t_equiv,atime
                    !
                    tot_q=tot_q*volume/atime*rho_equiv*1.e6/rmassq
                    tot_h=tot_h*volume/atime*rho_equiv*1.e6/rmassh
                    tot_o=tot_o*volume/atime*rho_equiv*1.e6/rmasso
                    write(6,*)'hires tot ions/s',m,tot_q,tot_h,tot_o
                    write(10,*)'hires tot ions/s',m,tot_q,tot_h,tot_o
                    !
                    tot_q=tot_q*proton_mass*rmassq
                    tot_h=tot_h*proton_mass*rmassh
                    tot_o=tot_o*proton_mass*rmasso
                    write(6,*)'hires tot  kg/s',m,tot_q,tot_h,tot_o
                    write(10,*)'hires tot  kg/s',m,tot_q,tot_h,tot_o
                enddo
                !
            endif  ! end ringo if
            !
        endif ! if(t.ge.ts1)
    enddo ! do while(t.lt.tmax)
          
    call clsgks
end

!
!     **************************************************
!
subroutine set_imf(bx,by,bz,bx0,by0,bz0,bxp,byp,bzp, &
        qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz, &
        hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz, &
        orho,opresx,opresy,opresz,opx,opy,opz, &
        rmassq,rmassh,rmasso,epres, &
        qpresxy,qpresxz,qpresyz, &
        hpresxy,hpresxz,hpresyz, &
        opresxy,opresxz,opresyz, &
        rhop,svxp,svyp,svzp,svelx,spress, &
        ti_te,rho_frac,nx,ny,nz,ngrd)
    !
    !     set imf boundary conditions
    !
    !
    dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd),bxp(ny,nz), &
        by(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),byp(ny,nz), &
        bz(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd),bzp(ny,nz), &
        rhop(ny,nz),svxp(ny,nz),svyp(ny,nz),svzp(ny,nz), &
!
        qrho(nx,ny,nz,ngrd),qpx(nx,ny,nz,ngrd), &
        qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
        qpresx(nx,ny,nz,ngrd),qpresy(nx,ny,nz,ngrd), &
        qpresz(nx,ny,nz,ngrd), qpresxy(nx,ny,nz,ngrd), &
        qpresxz(nx,ny,nz,ngrd),qpresyz(nx,ny,nz,ngrd), &
!
        hrho(nx,ny,nz,ngrd),hpx(nx,ny,nz,ngrd), &
        hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
        hpresx(nx,ny,nz,ngrd),hpresy(nx,ny,nz,ngrd), &
        hpresz(nx,ny,nz,ngrd), hpresxy(nx,ny,nz,ngrd), &
        hpresxz(nx,ny,nz,ngrd),hpresyz(nx,ny,nz,ngrd), &
! 
        orho(nx,ny,nz,ngrd),opx(nx,ny,nz,ngrd), &
        opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
        opresx(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
        opresy(nx,ny,nz,ngrd), opresxy(nx,ny,nz,ngrd), &
        opresxz(nx,ny,nz,ngrd),opresyz(nx,ny,nz,ngrd), &
    !
    epres(nx,ny,nz,ngrd)
    !
    m=ngrd
    !
    i=1
    do k=1,nz
        do j=1,ny
            bx(i,j,k,m)=bxp(j,k)-bx0(i,j,k,m)
            by(i,j,k,m)=byp(j,k)-by0(i,j,k,m)
            bz(i,j,k,m)=bzp(j,k)-bz0(i,j,k,m)
            !
            !
            qrho(i,j,k,m)=rhop(j,k)
            qpx(i,j,k,m)=qrho(i,j,k,m)*svxp(j,k)
            qpy(i,j,k,m)=qrho(i,j,k,m)*svyp(j,k)
            qpz(i,j,k,m)=qrho(i,j,k,m)*svzp(j,k)
            qpresx(i,j,k,m)=0.5*spress
            qpresy(i,j,k,m)=qpresx(i,j,k,m)
            qpresz(i,j,k,m)=qpresx(i,j,k,m)
            qpresxy(i,j,k,m)=0.
            qpresxz(i,j,k,m)=0.
            qpresyz(i,j,k,m)=0.
            !
            hrho(i,j,k,m)=hfrac*frac_h*rhop(j,k)
            hpx(i,j,k,m)=hrho(i,j,k,m)*svxp(j,k)
            hpy(i,j,k,m)=hrho(i,j,k,m)*svyp(j,k)
            hpz(i,j,k,m)=hrho(i,j,k,m)*svzp(j,k)
            hpresx(i,j,k,m)=0.5*spress*hfrac
            hpresy(i,j,k,m)=hpresx(i,j,k,m)
            hpresz(i,j,k,m)=hpresx(i,j,k,m)
            hpresxy(i,j,k,m)=0.
            hpresxz(i,j,k,m)=0.
            hpresyz(i,j,k,m)=0.
            !
            orho(i,j,k,m)=ofrac*frac_o*rhop(j,k)
            opx(i,j,k,m)=orho(i,j,k,m)*svxp(j,k)
            opy(i,j,k,m)=orho(i,j,k,m)*svyp(j,k)
            opz(i,j,k,m)=orho(i,j,k,m)*svzp(j,k)
            opresx(i,j,k,m)=0.5*spress*ofrac
            opresy(i,j,k,m)=opresx(i,j,k,m)
            opresz(i,j,k,m)=opresx(i,j,k,m)
            opresxy(i,j,k,m)=0.
            opresxz(i,j,k,m)=0.
            opresyz(i,j,k,m)=0.
            !
            epres(i,j,k,m)=qpresx(i,j,k,m)/ti_te
            !
            avz=avz+svzp(j,k)
        enddo
    enddo
    !
    avz=avz/float(ny*nz)
    !      write(6,*)'set_imf avz',avz,ny,nz
    return
end

!
!     **************************************************
!
subroutine set_resist(rst,nx,ny,nz,mbndry,resist, &
    ijzero,numzero,ijmid,nummid,ijsrf,numsrf, &
    msrf,mmid,mzero,b0)
    !
    !      set resistance around the earth :
    !        include dayside and auroral conductivities
    !      magntiude set by resist
    !      width is set by del_ang=3.3 degrees = 0.058 rads
    !      radial dropoff as alpha=-8
    !      shifted of dipole axis by 2.5 degress
    !      radius of 22.5 degrees
    !
    dimension rst(nx,ny,nz,mbndry)
    integer ijsrf(mbndry,3,msrf),ijmid(mbndry,3,mmid), &
        ijzero(mbndry,3,mzero)
    integer numsrf(mbndry),nummid(mbndry),numzero(mbndry)
    !
    !     resistivity characteristics
    !        initialize
    !
    rst=0.
    !     write(6,*)'set_resist',mbndry,nx,ny,nz,b0
    !     write(6,*)numsrf,nummid,numzero
    !
    !     interior resistivity
    !
    do m=1,mbndry
        do n=1,numzero(m)
            i=ijzero(m,1,n)
            j=ijzero(m,2,n)
            k=ijzero(m,3,n)
            rst(i,j,k,m)=b0/resist
        enddo
        !
        !     lower ionosphere resistivity
        !
        do n=1,nummid(m)
            i=ijmid(m,1,n)
            j=ijmid(m,2,n)
            k=ijmid(m,3,n)
            rst(i,j,k,m)=0.5*b0/resist
        enddo
        !
        !     upper ionosphere
        !
        do n=1,numsrf(m)
            i=ijsrf(m,1,n)
            j=ijsrf(m,2,n)
            k=ijsrf(m,3,n)
            rst(i,j,k,m)=0.125*b0/resist
        enddo
        !
    enddo ! end m loop
    return
end

!
!     ************************************************
!
subroutine store_array(bx,bx0,nx,ny,nz,ngrd,m)
    !
    !
    !     calculates the total magnetic field from the perturbed and
    !        stationary magnetic field
    !
    dimension bx(nx,ny,nz,ngrd),bx0(nx,ny,nz,ngrd)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                bx(i,j,k,m)=bx0(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end

!
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
!     ************************************************
!
subroutine fnd_fld(bx,btx,nx,ny,nz,ngrd,m)
    !
    !
    !     calculates the total dynamic magnetic field only
    !        dipole field added separately by by setting add_dip=.true.
    !
    dimension bx(nx,ny,nz,ngrd),btx(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                btx(i,j,k)=btx(i,j,k)+bx(i,j,k,m)
            enddo
        enddo
    enddo
    !
    return
end

!
!     ************************************************
!
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

!
!     ************************************************
!
subroutine fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
    opx,opy,opz,orho,vx,vy,vz,nx,ny,nz,ngrd,m, &
    rmassq,rmassh,rmasso)
    !
    !     converts momentum into velocity for graphics
    !
    dimension qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd), &
        qpz(nx,ny,nz,ngrd),qrho(nx,ny,nz,ngrd), &
        hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd), &
        hpz(nx,ny,nz,ngrd),hrho(nx,ny,nz,ngrd), &
        opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd), &
        opz(nx,ny,nz,ngrd),orho(nx,ny,nz,ngrd), &
        vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)
    !
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                qden=(qrho(i,j,k,m)+0.000001)/rmassq
                hden=(hrho(i,j,k,m)+0.000001)/rmassh
                oden=(orho(i,j,k,m)+0.000001)/rmasso
                tden=qden+hden+oden
                !
                vx(i,j,k)=(qpx(i,j,k,m)/rmassq+hpx(i,j,k,m)/rmassh &
                    +opx(i,j,k,m)/rmasso)/tden
                vy(i,j,k)=(qpy(i,j,k,m)/rmassq+hpy(i,j,k,m)/rmassh &
                    +opy(i,j,k,m)/rmasso)/tden
                vz(i,j,k)=(qpz(i,j,k,m)/rmassq+hpz(i,j,k,m)/rmassh &
                    +opz(i,j,k,m)/rmasso)/tden
            enddo
        enddo
    enddo
    !
    return
end

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
                !
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
!
subroutine fnd_prss(px,py,pz,rho,erg,afld,nx,ny,nz,ngrd, &
    m,gamma1,igo)
    !
    !     now using pressure equation rather than an erg equation
    !
    dimension px(nx,ny,nz,ngrd),py(nx,ny,nz,ngrd), &
        pz(nx,ny,nz,ngrd),rho(nx,ny,nz,ngrd), &
        erg(nx,ny,nz,ngrd),afld(nx,ny,nz)
    !
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                afld(i,j,k)=erg(i,j,k,m)
            enddo
        enddo
    enddo
    !
    !     if graphics refine plots
    !
    if(igo.le.0)return
    
    !$omp  parallel do
    do k=1,nz
        do j=1,ny
            do i=1,nx
                afld(i,j,k)=sqrt(abs(afld(i,j,k)))
                afld(i,j,k)=amin1(afld(i,j,k),1.2)
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

!
!     *************************************************
!
subroutine mak_dip_all(bx0,by0,bz0,nx,ny,nz,ngrd, &
    mbndry,ijzero,numzero,mzero,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      initialize static magnetic field along entire grid
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
        grd_ymin(ngrd),grd_ymax(ngrd), &
        grd_zmin(ngrd),grd_zmax(ngrd)
    dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd), &
        bz0(nx,ny,nz,ngrd)
    integer ijzero(mbndry,3,mzero),numzero(mbndry)
    !
    !
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    !$omp  parallel do
    do m=1,ngrd
        dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
        dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
        dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
        !
        do k=1,nz
            az=grd_zmin(m)+dz*(k-1)
            z1=(az-zdip)
            !
            do  j=1,ny
                ay=grd_ymin(m)+dy*(j-1)
                y1=(ay-ydip)
                !
                do  i=1,nx
                    ax=grd_xmin(m)+dx*(i-1)
                    x1=(ax-xdip)
                    !
                    !       rotate corrdinates for planet motion
                    !
                    xr=x1*cos_rot+y1*sin_rot
                    yr=-x1*sin_rot+y1*cos_rot
                    !
                    !
                    !       tilt space to dipole space
                    !
                    xp=xr*cos_tilt-z1*sin_tilt
                    yp=yr
                    zp=xr*sin_tilt+z1*cos_tilt
                    !
                    x2=xp**2
                    y2=yp**2
                    z2=zp**2
                    ar=sqrt(x2+y2+z2)
                    !
                    !        cartesian equivalent
                    !
                    bmag=-b0/ar**5
                    dbx=-3.*bmag*xp*zp
                    dby=-3.*bmag*yp*zp
    
                    dbz=bmag*(x2+y2-2.*z2)
                    !
                    !        tilt b field back to coordinate space
                    !
                    rbx=dbx*cos_tilt+dbz*sin_tilt
                    rby=dby
                    rbz=-dbx*sin_tilt+dbz*cos_tilt
                    !
                    !         rotate b
                    !
                    abx=rbx*cos_rot-rby*sin_rot
                    aby=rbx*sin_rot+rby*cos_rot
                    abz=rbz
                    !
                    if(ar.gt.rearth-1.5) then
                        bx0(i,j,k,m)=abx
                        by0(i,j,k,m)=aby
                        bz0(i,j,k,m)=abz
                    else
                        bx0(i,j,k,m)=0.
                        by0(i,j,k,m)=0.
                        bz0(i,j,k,m)=0.
                    endif
                    !
                enddo
            enddo
        enddo
        !
    enddo
    !
    !     boundary conditions
    !
    do m=1,mbndry
        !$omp  parallel do
        do n=1,numzero(m)
            !
            !        get coords of point
            !
            i=ijzero(m,1,n)
            j=ijzero(m,2,n)
            k=ijzero(m,3,n)
            !
            bx0(i,j,k,m)=0.
            by0(i,j,k,m)=0.
            bz0(i,j,k,m)=0.
        enddo
    enddo
    !
    return
end

!
!     *************************************************
!
subroutine mak_dip_grd(bx0,by0,bz0,nx,ny,nz,ngrd, &
    mbndry,m,ijzero,numzero,mzero,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      initialize static magnetic field along entire grid
    !
    
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
        grd_ymin(ngrd),grd_ymax(ngrd), &
        grd_zmin(ngrd),grd_zmax(ngrd)
    dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd), &
        bz0(nx,ny,nz,ngrd)
    integer ijzero(mbndry,3,mzero),numzero(mbndry)
    !
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    !     write(6,*)'in mak dip',m,rearth,b0
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !$omp  parallel do
    do k=1,nz
        az=grd_zmin(m)+dz*(k-1)
        z1=(az-zdip)
        !
        do  j=1,ny
            ay=grd_ymin(m)+dy*(j-1)
            y1=(ay-ydip)
            !
            do  i=1,nx
                ax=grd_xmin(m)+dx*(i-1)
                x1=(ax-xdip)
                !
                !        determine magnetic dipole field
                !
                !
                !        determine magnetic dipole field
                !
                !       rotate corrdinates for planet motion
                !
                xr=x1*cos_rot+y1*sin_rot
                yr=-x1*sin_rot+y1*cos_rot
                !
                !
                !       tilt space to dipole space
                !
                xp=xr*cos_tilt-z1*sin_tilt
                yp=yr
                zp=xr*sin_tilt+z1*cos_tilt
                !
                x2=xp**2
                y2=yp**2
                z2=zp**2
                ar=sqrt(x2+y2+z2)
                !
                !        cartesian equivalent
                !
                bmag=-b0/ar**5
                dbx=-3.*bmag*xp*zp
                dby=-3.*bmag*yp*zp
    
                dbz=bmag*(x2+y2-2.*z2)
                !
                !        tilt b field back to coordinate space
                !
                rbx=dbx*cos_tilt+dbz*sin_tilt
                rby=dby
                rbz=-dbx*sin_tilt+dbz*cos_tilt
                !
                !         rotate b
                !
                abx=rbx*cos_rot-rby*sin_rot
                aby=rbx*sin_rot+rby*cos_rot
                abz=rbz
                !
                if(ar.gt.rearth-1.5) then
                    bx0(i,j,k,m)=abx
                    by0(i,j,k,m)=aby
                    bz0(i,j,k,m)=abz
                else
                    bx0(i,j,k,m)=0.
                    by0(i,j,k,m)=0.
                    bz0(i,j,k,m)=0.
                endif
                !
            enddo
        enddo
    enddo
    !
    !
    !     boundary conditions
    !
    if(m.le.mbndry)then
        do n=1,numzero(m)
            !
            !        get coords of point
            !
            i=ijzero(m,1,n)
            j=ijzero(m,2,n)
            k=ijzero(m,3,n)
            !
            bx0(i,j,k,m)=0.
            by0(i,j,k,m)=0.
            bz0(i,j,k,m)=0.
        enddo
    endif
    !
    !
    return
end

!
!     *************************************************
!
subroutine mak_dip_moon(bx0,by0,bz0,nx,ny,nz,ngrd, &
    mbndry,m,ijzero,numzero,mzero,rearth, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    !
    !      initialize static magnetic field along entire grid
    !
    
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
        grd_ymin(ngrd),grd_ymax(ngrd), &
        grd_zmin(ngrd),grd_zmax(ngrd)
    dimension bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd), &
        bz0(nx,ny,nz,ngrd)
    integer ijzero(mbndry,3,mzero),numzero(mbndry)
    !
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    !     write(6,*)'in mak dip moon',m,nx,ny,nz,ngrd,mbndry
    !
    dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
    dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
    dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
    !
    !$omp  parallel do
    do k=1,nz
        az=grd_zmin(m)+dz*(k-1)
        z1=(az-zdip)
        !
        do  j=1,ny
            ay=grd_ymin(m)+dy*(j-1)
            y1=(ay-ydip)
            !
            do  i=1,nx
                ax=grd_xmin(m)+dx*(i-1)
                x1=(ax-xdip)
                !
                !        determine magnetic dipole field
                !
                !
                !        determine magnetic dipole field
                !
                !       rotate corrdinates for planet motion
                !
                xr=x1*cos_rot+y1*sin_rot
                yr=-x1*sin_rot+y1*cos_rot
                !
                !
                !       tilt space to dipole space
                !
                xp=xr*cos_tilt-z1*sin_tilt
                yp=yr
                zp=xr*sin_tilt+z1*cos_tilt
                !
                x2=xp**2
                y2=yp**2
                z2=zp**2
                ar=sqrt(x2+y2+z2)+1.e-8
                !        write(6,*)x2,y2,z2,ar
                !
                !        cartesian equivalent
                !
                bmag=-b0/ar**5
                dbx=-3.*bmag*xp*zp
                dby=-3.*bmag*yp*zp
    
                dbz=bmag*(x2+y2-2.*z2)
                !
                !        tilt b field back to coordinate space
                !
                rbx=dbx*cos_tilt+dbz*sin_tilt
                rby=dby
                rbz=-dbx*sin_tilt+dbz*cos_tilt
                !
                !         rotate b
                !
                abx=rbx*cos_rot-rby*sin_rot
                aby=rbx*sin_rot+rby*cos_rot
                abz=rbz
                !
                if(ar.gt.rearth-1.5) then
                    bx0(i,j,k,m)=abx
                    by0(i,j,k,m)=aby
                    bz0(i,j,k,m)=abz
                else
                    bx0(i,j,k,m)=0.
                    by0(i,j,k,m)=0.
                    bz0(i,j,k,m)=0.
                endif
                !
            enddo
        enddo
    enddo
    !
    !
    !
    return
end

!
!     *************************************************
!
subroutine dipole(abx,aby,abz,ax,ay,az)
    !
    !     calculates magnetic field for a dipole
    !
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
        sin_tilt,cos_tilt,b0
    !
    sin_rot=sin(rot_angle)
    cos_rot=cos(rot_angle)
    !
    x1=(ax-xdip)
    y1=(ay-ydip)
    z1=(az-zdip)
    !
    !       rotate corrdinates for planet motion
    !
    xr=x1*cos_rot+y1*sin_rot
    yr=-x1*sin_rot+y1*cos_rot
    !
    !
    !       tilt space to dipole space
    !
    xp=xr*cos_tilt-z1*sin_tilt
    yp=yr
    zp=xr*sin_tilt+z1*cos_tilt
    !
    x2=xp**2
    y2=yp**2
    z2=zp**2
    ar=sqrt(x2+y2+z2)
    !
    !        cartesian equivalent
    !
    bmag=-b0/ar**5
    dbx=-3.*bmag*xp*zp
    dby=-3.*bmag*yp*zp
    dbz=bmag*(x2+y2-2.*z2)
    !
    !        tilt b field back to coordinate space
    !
    rbx=dbx*cos_tilt+dbz*sin_tilt
    rby=dby
    rbz=-dbx*sin_tilt+dbz*cos_tilt
    !
    !         rotate b
    !
    abx=rbx*cos_rot-rby*sin_rot
    aby=rbx*sin_rot+rby*cos_rot
    abz=rbz
    !
    return
end

!
!     **************************************************************
!
subroutine limcraft(xcraft,ncraft,re_equiv,ngrd, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax)
    !
    !
    !      tests to see whether spacecraft is in the system and
    !           resets their position in not
    !
    dimension grd_xmin(ngrd),grd_xmax(ngrd), &
        grd_ymin(ngrd),grd_ymax(ngrd), &
        grd_zmin(ngrd),grd_zmax(ngrd)
    dimension xcraft(4,ncraft)
    !
    abit=0.001
    m=ngrd
    do n=1,ncraft
        xcraft(1,n)=amax1(xcraft(1,n),(grd_xmin(m)+abit)*re_equiv)
        xcraft(1,n)=amin1(xcraft(1,n),(grd_xmax(m)-abit)*re_equiv)
        xcraft(2,n)=amax1(xcraft(2,n),(grd_ymin(m)+abit)*re_equiv)
        xcraft(2,n)=amin1(xcraft(2,n),(grd_ymax(m)-abit)*re_equiv)
        xcraft(3,n)=amax1(xcraft(3,n),(grd_zmin(m)+abit)*re_equiv)
        xcraft(3,n)=amin1(xcraft(3,n),(grd_zmax(m)-abit)*re_equiv)
    enddo
    !
    return
end
    
!
!     **************************************************
!
subroutine visual( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
    orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
    epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
    curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
    tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
    nx,ny,nz,ngrd,xspac, &
    cross,along,flat,xcraft,ncraft,re_equiv, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,ut,b_equiv,ti_te,rho_equiv)
    !
    real grd_xmin(ngrd),grd_xmax(ngrd),grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd),xspac(ngrd)
    real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    real bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
    real efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz)
    real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz), &
    tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2), &
    cross(ny,nz),along(nx,nz),flat(nx,ny)
    real xcraft(4,ncraft)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    character*5 wd1
    character*8 label
    !
    logical add_dip
    !
    do m=ngrd,1,-1
        rx=xspac(m)
        ry=xspac(m)
        rz=xspac(m)
        !
        ymin=grd_ymin(m)
        ymax=grd_ymax(m)
        zmin=grd_zmin(m)
        zmax=grd_zmax(m)
        xmin=grd_xmin(m)
        xmax=grd_xmax(m)
        !
        xcut=((xmin+xmax))/2.
        !
        add_dip=.false.
        !
        call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
        call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
        call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
        if(m.le.3)then
            preslim=16.0/float(m)
            po=preslim*0.33
        else
            preslim=16.0/(m-1.)
            po=preslim*0.33
        endif
        !
        call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=qpresx(i,j,k,m)
                    presy=qpresy(i,j,k,m)
                    presz=qpresz(i,j,k,m)
                    presmag=sqrt(presx**2+presy**2+presz**2)+1.e-11
                    apres=(presx+presy+presz)/3.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    vbx=avx-vcrossb_x
                    vby=avy-vcrossb_y
                    vbz=avz-vcrossb_z
                    vbmag=sqrt((vbx**2+vby**2+vbz**2))+1.e-8
                    !
                    p_para=sqrt((presx*abx)**2+(presy*aby)**2+(presz*abz)**2) &
                        /(bmag)
                    p_cross=sqrt((presx*vcrossb_x)**2+(presy*vcrossb_y)**2+ &
                    (presz*vcrossb_z)**2) &
                        /(vcmag)
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
                    !
    
                    tvx(i,j,k)=sqrt(presx)
                    tvy(i,j,k)=sqrt(presy)
                    tvz(i,j,k)=sqrt(presz)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=presx/apres+0.0001
                    efldy(i,j,k)=presy/apres+0.0001
                    efldz(i,j,k)=presz/apres+0.0001
                enddo
            enddo
        enddo
        !
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'qpresx',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'qpresy',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'qpresz',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqprx',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqpry',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqprz',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
    
        call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=hpresx(i,j,k,m)
                    presy=hpresy(i,j,k,m)
                    presz=hpresz(i,j,k,m)
                    presmag=sqrt(presx**2+presy**2+presz**2)+1.e-11
                    apres=(presx+presy+presz)/3.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    vbx=avx-vcrossb_x
                    vby=avy-vcrossb_y
                    vbz=avz-vcrossb_z
                    vbmag=sqrt((vbx**2+vby**2+vbz**2))+1.e-8
                    !
                    p_para=sqrt((presx*abx)**2+(presy*aby)**2+(presz*abz)**2) &
                        /(bmag)
                    p_cross=sqrt((presx*vcrossb_x)**2+(presy*vcrossb_y)**2+ &
                    (presz*vcrossb_z)**2) &
                        /(vcmag)
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
                    !
    
                    tvx(i,j,k)=sqrt(presx)
                    tvy(i,j,k)=sqrt(presy)
                    tvz(i,j,k)=sqrt(presz)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=presx/apres+0.0001
                    efldy(i,j,k)=presy/apres+0.0001
                    efldz(i,j,k)=presz/apres+0.0001
                enddo
            enddo
        enddo
        !
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'hpresx',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'hpresy',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'hpresz',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhprx',3,18,1,2.0,0.75,1.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhpry',3,18,1,2.0,0.75,1.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhprz',3,18,1,2.0,0.75,1.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
    
        call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=opresx(i,j,k,m)
                    presy=opresy(i,j,k,m)
                    presz=opresz(i,j,k,m)
                    presmag=sqrt(presx**2+presy**2+presz**2)+1.e-11
                    apres=(presx+presy+presz)/3.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    vbx=avx-vcrossb_x
                    vby=avy-vcrossb_y
                    vbz=avz-vcrossb_z
                    vbmag=sqrt((vbx**2+vby**2+vbz**2))+1.e-8
                    !
                    p_para=sqrt((presx*abx)**2+(presy*aby)**2+(presz*abz)**2) &
                        /(bmag)
                    p_cross=sqrt((presx*vcrossb_x)**2+(presy*vcrossb_y)**2+ &
                    (presz*vcrossb_z)**2) &
                        /(vcmag)
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
                    !
    
                    tvx(i,j,k)=sqrt(presx)
                    tvy(i,j,k)=sqrt(presy)
                    tvz(i,j,k)=sqrt(presz)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=presx/apres+0.0001
                    efldy(i,j,k)=presy/apres+0.0001
                    efldz(i,j,k)=presz/apres+0.0001
                enddo
            enddo
        enddo
        !
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'opresx',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'opresy',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'opresz',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'roprx',3,18,1,2.0,0.75,1.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'ropry',3,18,1,2.0,0.75,1.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'roprz',3,18,1,2.0,0.75,1.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=sqrt(epres(i,j,k,m))
                enddo
            enddo
        enddo
        !
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'epres',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=abs(qpresx(i,j,k,m)+hpresx(i,j,k,m) &
                        +opresx(i,j,k,m))/3.
                    efldy(i,j,k)=abs(qpresy(i,j,k,m)+hpresy(i,j,k,m) &
                        +opresy(i,j,k,m))/3.
                    efldz(i,j,k)=abs(qpresz(i,j,k,m)+hpresz(i,j,k,m) &
                        +opresz(i,j,k,m))/3.
                    tvx(i,j,k)=sqrt(efldx(i,j,k)+efldy(i,j,k)+efldz(i,j,k) &
                        +epres(i,j,k,m))
                enddo
            enddo
        enddo
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'tpres',3,18,1,2.0,preslim*3., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !      trace velocity streams
        !
        call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m, &
            rmassq,rmassh,rmasso)
        !
        call conflow(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'pres-vel',3,11,1,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !     find total magnetic field
        !
        call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
        call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
        call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
    
        write(wd1,'(i3)')m
        label='box '//wd1
        call concross(tvx,curx,cury,curz,bsx,bsy,bsz, &
            nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,rearth, &
            xmin,xmax,ymin,ymax,zmin,zmax, &
            ut,label,3,11,2.0,add_dip,1,-2,start, &
            tx,ty,tz,tt,tg1,tg2,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        call contop(tvx,curx,cury,curz,bsx,bsy,bsz, &
            nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,rearth, &
            xmin,xmax,ymin,ymax,zmin,zmax, &
            ut,label,3,11,2.0,add_dip,1,0,start, &
            tx,ty,tz,tt,tg1,tg2,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       magnitude of b
        !
        blim=0.
    
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldx(i,j,k)=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                    +bsz(i,j,k)**2)
                    efldx(i,j,k)=alog10(b_equiv*efldx(i,j,k)+1.e-10)
                    blim=amax1(blim,efldx(i,j,k))
                enddo
            enddo
        enddo
        !
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'bmag',3,18,1,2.0,0.1,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        blim=0.001
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldx(i,j,k)=bsz(i,j,k)
                    if (efldx(i,j,k).gt.blim)efldx(i,j,k)=blim
                    if (efldx(i,j,k).lt.-blim)efldx(i,j,k)=-blim
                    efldx(i,j,k)=efldx(i,j,k)+1.001*blim
                enddo
            enddo
        enddo
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,' bz ',3,12,1,2.0,2.*blim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !      calculate alfven mach number
        !
        call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m, &
            rmassq,rmassh,rmasso)
        !
        dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
        dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
        dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    aden=qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m)+1.e-8
                    efldx(i,j,k)=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                        +bsz(i,j,k)**2)
                    alf_spd=(efldx(i,j,k)+1.e-2)/sqrt(aden)
                    apres=tvx(i,j,k)+1.e-10
                    cs_spd=sqrt(apres/aden)
                    aspd=sqrt(curx(i,j,k)**2+cury(i,j,k)**2 &
                        +curz(i,j,k)**2)
                    efldy(i,j,k)=aspd/alf_spd
                    efldz(i,j,k)=aspd/cs_spd
                    !
                    !       calculate rotation velocity
                    !
    
                    ay=grd_ymin(m)+dy*(j-1)-ydip
                    ax=grd_xmin(m)+dx*(i-1)-xdip
                    ar=sqrt(ax**2+ay**2)+1.e-8
    
                    rspd=abs( (curx(i,j,k)*ay-cury(i,j,k)*ax) / ar)
                    radius=ar*re_equiv
                    rot_spd=radius*v_rot+.001
                    efldx(i,j,k)=rspd/rot_spd
                enddo
            enddo
        enddo
        !
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'alf_mach',3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'cs_mach',3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rot_mach',3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !      plot individual temperatures
        !
        tempx=0.50
        temph=2.*tempx
        tempo=4.*tempx
        if(m.gt.1) then
            rho_lim=10.
        else
            rho_lim=5.0
        endif
        rfrac=0.5
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    qden=qrho(i,j,k,m)/rmassq
                    if(qden.gt.0.001)then
                        apres=(qpresx(i,j,k,m)+qpresy(i,j,k,m) &
                            +qpresz(i,j,k,m))/3.
                        bsx(i,j,k)=amin1(tempx,sqrt(apres/qden))
                    else
                        bsx(i,j,k)=0.
                    endif
                    hden=hrho(i,j,k,m)/rmassh
                    if(hden.gt.0.0005)then
                        apres=(hpresx(i,j,k,m)+hpresy(i,j,k,m) &
                            +hpresz(i,j,k,m))/3.
                        bsy(i,j,k)=amin1(temph,sqrt(apres/hden))
                    else
                        bsy(i,j,k)=0.
                    endif
                    !
                    oden=orho(i,j,k,m)/rmasso
                    if(oden.gt.0.00002)then
                        apres=(opresx(i,j,k,m)+opresy(i,j,k,m) &
                            +opresz(i,j,k,m))/3.
                        bsz(i,j,k)=amin1(tempo,sqrt(apres/oden))
                    else
                        bsz(i,j,k)=0.
                    endif
                    !
                    efldx(i,j,k)=sqrt(epres(i,j,k,m)/(qden+hden+oden))
                    !
                enddo
            enddo
        enddo
        !
        call conhot(bsx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'q temp',3,18,1,2.0,tempx, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(bsy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'h temp',3,18,1,2.0,temph, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(bsz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'o temp',3,18,1,2.0,tempo, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'e temp',3,18,1,2.0,tempx/sqrt(ti_te), &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : solar wind
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !       call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)
                    efldx(i,j,k)=alog10(efldx(i,j,k)*rho_equiv)+6. ! per m**3
                enddo
            enddo
        enddo
        !
        if(m.le.3)then
            plot_min=5.
            plot_max=9.
        else
            plot_min=4.-0.5*(m-3)
            plot_max=8.-0.5*(m-3)
        endif
        !
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'q den v',3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : ionospheric h
        !
        !      call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldy(i,j,k)=hrho(i,j,k,m)/rmassh
                    efldy(i,j,k)=alog10(efldy(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m,&
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'h den v',3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : ionospheric o
        !
        !       call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldz(i,j,k)=orho(i,j,k,m)/rmasso
                    efldz(i,j,k)=alog10(efldz(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'o den v',3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : total
        !
        !      call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
        !    +       opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m,
        !    +      rmassq,rmassh,rmasso)
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,m)/rmasso+ &
                        hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
                    efldx(i,j,k)=alog10(tden*rho_equiv)+6.
                    if(efldx(i,j,k).lt.0.1)then
                        curx(i,j,k)=0.
                        cury(i,j,k)=0.
                        curz(i,j,k)=0.
                    endif
                enddo
            enddo
        enddo
        !
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'tden v',3,18,1,2.0,plot_min,plot_max, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,m)/rmasso+ &
                    hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
                    efldy(i,j,k)=(hrho(i,j,k,m)/rmassh)/(tden+1.e-6)
                    efldz(i,j,k)=(orho(i,j,k,m)/rmasso)/(tden+1.e-6)
                    efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)/(tden+1.e-6)
                enddo
            enddo
        enddo
        !
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqdens',3,18,1,2.0,0.5, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhdens',3,18,1,2.0,1.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rodens',3,18,1,2.0,0.5, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    enddo
    !
    return
end
    
!
!     **************************************************
!
subroutine visual_hires( &
    qrho,qpresx,qpresy,qpresz,qpx,qpy,qpz,rmassq, &
    hrho,hpresx,hpresy,hpresz,hpx,hpy,hpz,rmassh, &
    orho,opresx,opresy,opresz,opx,opy,opz,rmasso, &
    epres,bx,by,bz,bx0,by0,bz0,bsx,bsy,bsz, &
    curx,cury,curz,efldx,efldy,efldz,tvx,tvy,tvz, &
    tx,ty,tz,tg1,tg2,tt,work,mx,my,mz,mz2,muvwp2, &
    nx,ny,nz,ngrd,xspac, &
    cross,along,flat,xcraft,ncraft,re_equiv, &
    grd_xmin,grd_xmax,grd_ymin,grd_ymax, &
    grd_zmin,grd_zmax,ut,b_equiv,ti_te,rho_equiv)
    !
    real grd_xmin(ngrd),grd_xmax(ngrd),grd_ymin(ngrd),grd_ymax(ngrd), &
    grd_zmin(ngrd),grd_zmax(ngrd),xspac(ngrd)
    real bx(nx,ny,nz,ngrd),by(nx,ny,nz,ngrd),bz(nx,ny,nz,ngrd), &
    qpx(nx,ny,nz,ngrd),qpy(nx,ny,nz,ngrd),qpz(nx,ny,nz,ngrd), &
    qrho(nx,ny,nz,ngrd),qpresx(nx,ny,nz,ngrd), &
    qpresy(nx,ny,nz,ngrd),qpresz(nx,ny,nz,ngrd), &
    opx(nx,ny,nz,ngrd),opy(nx,ny,nz,ngrd),opz(nx,ny,nz,ngrd), &
    orho(nx,ny,nz,ngrd),opresx(nx,ny,nz,ngrd), &
    opresy(nx,ny,nz,ngrd),opresz(nx,ny,nz,ngrd), &
    hpx(nx,ny,nz,ngrd),hpy(nx,ny,nz,ngrd),hpz(nx,ny,nz,ngrd), &
    hrho(nx,ny,nz,ngrd),hpresx(nx,ny,nz,ngrd), &
    hpresy(nx,ny,nz,ngrd),hpresz(nx,ny,nz,ngrd), &
    epres(nx,ny,nz,ngrd)
    real bx0(nx,ny,nz,ngrd),by0(nx,ny,nz,ngrd),bz0(nx,ny,nz,ngrd)
    real efldx(nx,ny,nz),efldy(nx,ny,nz),efldz(nx,ny,nz), &
    tvx(nx,ny,nz),tvy(nx,ny,nz),tvz(nx,ny,nz), &
    curx(nx,ny,nz),cury(nx,ny,nz),curz(nx,ny,nz), &
    bsx(nx,ny,nz),bsy(nx,ny,nz),bsz(nx,ny,nz)
    real tx(mx,my,mz),ty(mx,my,mz),tz(mx,my,mz),tg1(mx,my,mz), &
    tg2(mx,my,mz2),tt(mx,my,mz),work(muvwp2,muvwp2), &
    cross(ny,nz),along(nx,nz),flat(nx,ny)
    real xcraft(4,ncraft)
    common /rotation/v_rot,r_rot,rot_angle,xdip,ydip,zdip, &
    sin_tilt,cos_tilt,b0
    !
    character*5 wd1
    character*8 label
    !
    logical add_dip
    !
    do m=ngrd,1,-1
        rx=xspac(m)
        ry=xspac(m)
        rz=xspac(m)
        !
        ymin=grd_ymin(m)+ry
        ymax=grd_ymax(m)-ry
        zmin=grd_zmin(m)+rz
        zmax=grd_zmax(m)-rz
        xmin=grd_xmin(m)+rx
        xmax=grd_xmax(m)-rx
        !
        xcut=((xmin+xmax)/2.+xmax)/2.
        !
        add_dip=.false.
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
        preslim=16.
        po=preslim
        call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=qpresx(i,j,k,m)
                    presy=qpresy(i,j,k,m)
                    presz=qpresz(i,j,k,m)
                    presmag=sqrt(presx**2+presy**2+presz**2)+1.e-11
                    apres=(presx+presy+presz)/3.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    vbx=avx-vcrossb_x
                    vby=avy-vcrossb_y
                    vbz=avz-vcrossb_z
                    vbmag=sqrt((vbx**2+vby**2+vbz**2))+1.e-8
                    !
                    p_para=sqrt((presx*abx)**2+(presy*aby)**2+(presz*abz)**2) &
                        /(bmag)
                    p_cross=sqrt((presx*vcrossb_x)**2+(presy*vcrossb_y)**2+ &
                    (presz*vcrossb_z)**2) &
                        /(vcmag)
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
                    !
    
                    tvx(i,j,k)=sqrt(presx)
                    tvy(i,j,k)=sqrt(presy)
                    tvz(i,j,k)=sqrt(presz)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=presx/apres+0.0001
                    efldy(i,j,k)=presy/apres+0.0001
                    efldz(i,j,k)=presz/apres+0.0001
                enddo
            enddo
        enddo
        !
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'qpresx',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'qpresy',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'qpresz',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqprx',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqpry',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqprz',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
    
        call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=hpresx(i,j,k,m)
                    presy=hpresy(i,j,k,m)
                    presz=hpresz(i,j,k,m)
                    presmag=sqrt(presx**2+presy**2+presz**2)+1.e-11
                    apres=(presx+presy+presz)/3.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    vbx=avx-vcrossb_x
                    vby=avy-vcrossb_y
                    vbz=avz-vcrossb_z
                    vbmag=sqrt((vbx**2+vby**2+vbz**2))+1.e-8
                    !
                    p_para=sqrt((presx*abx)**2+(presy*aby)**2+(presz*abz)**2) &
                        /(bmag)
                    p_cross=sqrt((presx*vcrossb_x)**2+(presy*vcrossb_y)**2+ &
                    (presz*vcrossb_z)**2) &
                        /(vcmag)
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
                    !
    
                    tvx(i,j,k)=sqrt(presx)
                    tvy(i,j,k)=sqrt(presy)
                    tvz(i,j,k)=sqrt(presz)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=presx/apres+0.0001
                    efldy(i,j,k)=presy/apres+0.0001
                    efldz(i,j,k)=presz/apres+0.0001
                enddo
            enddo
        enddo
        !
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'hpresx',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'hpresy',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'hpresz',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhprx',3,18,1,2.0,0.9,1.1, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhpry',3,18,1,2.0,0.9,1.1, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhprz',3,18,1,2.0,0.9,1.1, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        !
        !
    
        call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    !
                    presx=opresx(i,j,k,m)
                    presy=opresy(i,j,k,m)
                    presz=opresz(i,j,k,m)
                    presmag=sqrt(presx**2+presy**2+presz**2)+1.e-11
                    apres=(presx+presy+presz)/3.+1.e-11
                    !
                    abx=bsx(i,j,k)
                    aby=bsy(i,j,k)
                    abz=bsz(i,j,k)
                    bmag=sqrt(abx**2+aby**2+abz**2)+1.e-12
                    !
                    avx=curx(i,j,k)
                    avy=cury(i,j,k)
                    avz=curz(i,j,k)
                    vmag=sqrt(avx**2+avy**2+avz**2)+1.e-8
                    !
                    vcrossb_x=(avy*abz-avz*aby)/bmag
                    vcrossb_y=-(avx*abz-avz*abx)/bmag
                    vcrossb_z=(avx*aby-avy*abx)/bmag
                    vcmag=sqrt(vcrossb_x**2+vcrossb_y**2+vcrossb_z**2)+1.e-12
                    !
                    !       find vparallel
                    !
                    vbx=avx-vcrossb_x
                    vby=avy-vcrossb_y
                    vbz=avz-vcrossb_z
                    vbmag=sqrt((vbx**2+vby**2+vbz**2))+1.e-8
                    !
                    p_para=sqrt((presx*abx)**2+(presy*aby)**2+(presz*abz)**2) &
                        /(bmag)
                    p_cross=sqrt((presx*vcrossb_x)**2+(presy*vcrossb_y)**2+ &
                    (presz*vcrossb_z)**2) &
                        /(vcmag)
                    p_perp=sqrt(abs(presmag**2-(p_para**2+p_cross**2)))
                    !
    
                    tvx(i,j,k)=sqrt(presx)
                    tvy(i,j,k)=sqrt(presy)
                    tvz(i,j,k)=sqrt(presz)
                    !
                    !       efldx(i,j,k)=p_para/((presmag+1.e-13)/sqrt(3.))
                    !       efldy(i,j,k)=p_cross/((presmag+1.e-13)/sqrt(3.))
                    !       efldz(i,j,k)=p_perp/((presmag+1.e-13)/sqrt(3.))
                    !
                    efldx(i,j,k)=presx/apres+0.0001
                    efldy(i,j,k)=presy/apres+0.0001
                    efldz(i,j,k)=presz/apres+0.0001
                enddo
            enddo
        enddo
        !
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'opresx',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'opresy',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(tvz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'opresz',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'roprx',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'ropry',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'roprz',3,18,1,2.0,0.5,2.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=sqrt(epres(i,j,k,m))
                enddo
            enddo
        enddo
        !
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'epres',3,18,1,2.0,preslim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=abs(qpresx(i,j,k,m)+hpresx(i,j,k,m) &
                    +opresx(i,j,k,m))/3.
                    efldy(i,j,k)=abs(qpresy(i,j,k,m)+hpresy(i,j,k,m) &
                    +opresy(i,j,k,m))/3.
                    efldz(i,j,k)=abs(qpresz(i,j,k,m)+hpresz(i,j,k,m) &
                    +opresz(i,j,k,m))/3.
                    tvx(i,j,k)=efldx(i,j,k)+efldy(i,j,k)+efldz(i,j,k) &
                    +epres(i,j,k,m)
                enddo
            enddo
        enddo
        call conhot(tvx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'tpres',3,18,1,2.0,preslim*3., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !
        !     find total magnetic field
        !
        call totfld(bx,bx0,bsx,nx,ny,nz,ngrd,m)
        call totfld(by,by0,bsy,nx,ny,nz,ngrd,m)
        call totfld(bz,bz0,bsz,nx,ny,nz,ngrd,m)
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !
        if(m.gt.3)then
            write(wd1,'(i3)')m
            label='box '//wd1
            call concross(tvx,curx,cury,curz,bsx,bsy,bsz, &
                nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,rearth, &
                xmin,xmax,ymin,ymax,zmin,zmax, &
                ut,label,3,11,2.0,add_dip,1,-2,start, &
                tx,ty,tz,tt,tg1,tg2,work,mx,my,mz,mz2,muvwp2, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    
            call contop(tvx,curx,cury,curz,bsx,bsy,bsz, &
                nx,ny,nz,1,1,m,xcraft,ncraft,re_equiv,rearth, &
                xmin,xmax,ymin,ymax,zmin,zmax, &
                ut,label,3,11,2.0,add_dip,1,0,start, &
                tx,ty,tz,tt,tg1,tg2,work,mx,my,mz,mz2,muvwp2, &
                grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        endif
        !
        !       magnitude of b
        !
        blim=0.
    
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldx(i,j,k)=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                    +bsz(i,j,k)**2)
                    efldx(i,j,k)=alog10(b_equiv*efldx(i,j,k)+1.e-10)
                    blim=amax1(blim,efldx(i,j,k))
                enddo
            enddo
        enddo
        !
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'bmag',3,18,1,2.0,2.7,3.7, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        blim=0.001
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldx(i,j,k)=bsz(i,j,k)
                    if (efldx(i,j,k).gt.blim)efldx(i,j,k)=blim
                    if (efldx(i,j,k).lt.-blim)efldx(i,j,k)=-blim
                    efldx(i,j,k)=efldx(i,j,k)+1.001*blim
                enddo
            enddo
        enddo
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,' bz ',3,18,1,2.0,2.*blim, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !      calculate alfven mach number
        !
        call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho, &
            opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m, &
            rmassq,rmassh,rmasso)
        !
        dx=(grd_xmax(m)-grd_xmin(m))/(nx-1.)
        dy=(grd_ymax(m)-grd_ymin(m))/(ny-1.)
        dz=(grd_zmax(m)-grd_zmin(m))/(nz-1.)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    aden=qrho(i,j,k,m)+hrho(i,j,k,m)+orho(i,j,k,m)+1.e-8
                    efldx(i,j,k)=sqrt(bsx(i,j,k)**2+bsy(i,j,k)**2 &
                    +bsz(i,j,k)**2)
                    alf_spd=(efldx(i,j,k)+1.e-2)/sqrt(aden)
                    apres=tvx(i,j,k)+1.e-10
                    cs_spd=sqrt(apres/aden)
                    aspd=sqrt(curx(i,j,k)**2+cury(i,j,k)**2 &
                    +curz(i,j,k)**2)
                    efldy(i,j,k)=aspd/alf_spd
                    efldz(i,j,k)=aspd/cs_spd
                    !
                    !       calculate rotation velocity
                    !
                    !
                    ay=grd_ymin(m)+dy*(j-1)-ydip
                    ax=grd_xmin(m)+dx*(i-1)-xdip
                    ar=sqrt(ax**2+ay**2)+1.e-8
                    !
                    !       rspd=sqrt(curx(i,j,k)**2+cury(i,j,k)**2)
                    !
                    rspd=abs( (curx(i,j,k)*ay-cury(i,j,k)*ax) / ar)
                    !
                    radius=ar*re_equiv
                    rot_spd=radius*v_rot+.001
                    efldx(i,j,k)=rspd/rot_spd
                enddo
            enddo
        enddo
        !
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'alf_mach',3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'cs_mach',3,18,1,2.0,4., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rot_mach',3,18,1,2.0,2.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !
        !      plot individual temperatures
        !
        tempx=0.40
        tempo=2.*tempx
        if(m.gt.1) then
            rho_lim=10.
        else
            rho_lim=5.0
        endif
        rfrac=0.5
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    qden=qrho(i,j,k,m)/rmassq
                    if(qden.gt.0.001)then
                        apres=(qpresx(i,j,k,m)+qpresy(i,j,k,m) &
                            +qpresz(i,j,k,m))/3.
                        bsx(i,j,k)=amin1(tempx,sqrt(apres/qden))
                    else
                        bsx(i,j,k)=0.
                    endif
                    hden=hrho(i,j,k,m)/rmassh
                    if(hden.gt.0.0005)then
                        apres=(hpresx(i,j,k,m)+hpresy(i,j,k,m) &
                            +hpresz(i,j,k,m))/3.
                        bsy(i,j,k)=amin1(tempo,sqrt(apres/hden))
                    else
                        bsy(i,j,k)=0.
                    endif
                    !
                    oden=orho(i,j,k,m)/rmasso
                    if(oden.gt.0.00002)then
                        apres=(opresx(i,j,k,m)+opresy(i,j,k,m) &
                            +opresz(i,j,k,m))/3.
                        bsz(i,j,k)=amin1(tempo,sqrt(apres/oden))
                    else
                        bsz(i,j,k)=0.
                    endif
                    !
                    efldx(i,j,k)=sqrt(epres(i,j,k,m)/(qden+hden+oden))
                    !
                enddo
            enddo
        enddo
        !
        !
        call conhot(bsx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'q temp',3,18,1,2.0,tempx, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(bsy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'h temp',3,18,1,2.0,tempx*2., &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(bsz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'o temp',3,18,1,2.0,tempo, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'e temp',3,18,1,2.0,tempx/sqrt(ti_te), &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : solar wind
        !
        call qvset(0.,curx,nx*ny*nz)
        call qvset(0.,cury,nx*ny*nz)
        call qvset(0.,curz,nx*ny*nz)
        !       call fnd_vel(qpx,qpy,qpz,qrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        do k=1,nz
            do j=1,ny
                do i=1,nx
                    efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)
                    efldx(i,j,k)=alog10(efldx(i,j,k)*rho_equiv)+6. ! per m**3
                enddo
            enddo
        enddo
        !
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'q den v',3,18,1,2.0,7.,9.99, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : ionospheric h
        !
        !      call fnd_vel(hpx,hpy,hpz,hrho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldy(i,j,k)=hrho(i,j,k,m)/rmassh
                    efldy(i,j,k)=alog10(efldy(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        call conlog(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'h den v',3,18,1,2.0,7.,9.99, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : ionospheric o
        !
        !       call fnd_vel(opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    efldz(i,j,k)=orho(i,j,k,m)/rmasso
                    efldz(i,j,k)=alog10(efldz(i,j,k)*rho_equiv)+6.
                enddo
            enddo
        enddo
        !
        call conlog(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'o den v',3,18,1,2.0,7.,9.99, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        !       plot mass flows : total
        !
        !      call fnd_vtot(qpx,qpy,qpz,qrho,hpx,hpy,hpz,hrho,
        !    +       opx,opy,opz,orho,curx,cury,curz,nx,ny,nz,ngrd,m,
        !    +      rmassq,rmassh,rmasso)
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,m)/rmasso+ &
                    hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
                    efldx(i,j,k)=alog10(tden*rho_equiv)+6.
                    if(efldx(i,j,k).lt.0.1)then
                        curx(i,j,k)=0.
                        cury(i,j,k)=0.
                        curz(i,j,k)=0.
                    endif
                enddo
            enddo
        enddo
        !
        call conlog(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'tden v',3,18,1,2.0,7.,9.99, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        !
        do  k=1,nz
            do  j=1,ny
                do  i=1,nx
                    tden=(orho(i,j,k,m)/rmasso+ &
                    hrho(i,j,k,m)/rmassh+qrho(i,j,k,m)/rmassq)
                    efldy(i,j,k)=(hrho(i,j,k,m)/rmassh)/(tden+1.e-6)
                    efldz(i,j,k)=(orho(i,j,k,m)/rmasso)/(tden+1.e-6)
                    efldx(i,j,k)=(qrho(i,j,k,m)/rmassq)/(tden+1.e-6)
                enddo
            enddo
        enddo
        !
        call conhot(efldx,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rqdens',3,18,1,2.0,0.5, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(efldy,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rhdens',3,18,1,2.0,1.0, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
        call conhot(efldz,curx,cury,curz,nx,ny,nz,1,1,m, &
            xmin,xmax,ymin,ymax,zmin,zmax,xcut, &
            ut,'rodens',3,18,1,2.0,0.25, &
            tx,ty,tz,tg1,tt,work,mx,my,mz,mz2,muvwp2, &
            grd_xmin,grd_xmax,grd_ymin,grd_ymax,grd_zmin,grd_zmax)
    enddo
    !
    return
end
