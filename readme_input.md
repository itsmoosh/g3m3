NOTES:
	@ Some quantities are set to 1 in sim units:
		#	mu_o, the vacuum permeability
		#	amu, in the values for rmassq, o, and h
	@ Many of these pre-date me (MAT), and I never updated the var names - confusion expected
	@ Some variable names have been changed. Previous names are indicated in [square brackets.] (MJS)
 
$ option
tmax - duration in units of t_equiv of physical time for single call to executable
ntgraph - number of gmeta frames for each var during tmax run
stepsz - upper limit to solution timestep (can drive instability, try lowering this until stable)
start - true for new run, false if continuing with input fluid files
tsave - duration in units of t_equiv per fluid file (number of fluid files * tsave should = tmax, above)
isotropic - deprecated after splitting code into two versions (can be removed/ignored)
 
 
$planet	[formerly earth]
xdip - x offset in sim units for dip moment (if @ center, enter small value in one avoid singularity)
ydip - y offset in sim units for dip moment (if @ center, enter small value in one avoid singularity)
zdip - z offset in sim units for dip moment (if @ center, enter small value in one avoid singularity)
r_inner - [formerly rearth] # grid points for the inner boundary radius (this * re_equiv = planetary radii of inner b)
torus_rad - minor radius of plasma torus in planetary radii. 1-sigma radius assuming Gaussian cross-sectional profile.
tilt1 - initial tilt of dipole moment for this execution call
tilt2 - final tilt of dipole moment for this execution call
tilting - turn on/off tilting
rmassq - mass in amu for solar wind species 
rmassh - mass in amu for ‘heavy’ species
rmasso - mass in amu for ‘oxygen’ species 
 
 
$speeds
cs_inner - sound speed of ‘heavy’ species at inner boundary (sets ion temperature)
alf_inner1 - initial alfven speed at equator of inner boundary for this exec (sets dipole strength). alfven speed in sim units = B/sqrt()
alf_inner2 - final alfven speed at equator of inner boundary  for this exec (sets dipole strength)
alpha_e - exponent for latitudinal density scaling at inner boundary (‘top hat’ (?) distribution)
denh_inner - number density in sim units for ‘heavy’ species at planetary inner boundary for this exec. denh_inner*rho_equiv = n_h (in cm^-3)
denh_torus - number density in sim units for 'heavy' species at moon inner boundary for this exec 
gravity - acceleration at planetary surface in SI
ti_te - ion to electron temperature ratio
gamma - polytropic index
ringo - control plasma torus conditional @Saturn or ring current outflow @Earth (see code)
update - ask robert … unsure, didn’t use
reload - ask robert … unsure, didn’t use
divb_lores - ensure that div B is zero in planetary grid
divb_hires - ensure that div B is zero in moon grid
reduct - multiplier for torus injected plasma momentum
t_torus - fraction of corotation energy to pickup ions in torus (\propto corotation vel^2)
aniso_factor - scale z pressure component of injected torus plasma (0 = purely perpendicular aniso, 1 = isotropic)
ani_q - deprecated (originally used to set definable aniso ratios per ion, I think?)
ani_h - deprecated (originally used to set definable aniso ratios per ion, I think?)
ani_o - deprecated (originally used to set definable aniso ratios per ion, I think?)
 
$windy
re_wind - deprecated, used to set distance scaling for solar wind
vx_wind1 - initial x velocity in sim units for stellar wind for this exec
vx_wind2 - final x velocity in sim units for stellar wind for this exec
vy_wind1 - initial y velocity in sim units for stellar wind for this exec
vy_wind2 - final y velocity in sim units for stellar wind for this exec
vz_wind1 - initial z velocity in sim units for stellar wind for this exec
vz_wind2 - final z velocity in sim units for stellar wind for this exec
alfx_wind1 - initial x alfven speed in sim units for stellar wind for this exec (sets IMF x)
alfx_wind2 - final x alfven speed in sim units for stellar wind for this exec (sets IMF x)
alfy_wind1 - initial y alfven speed in sim units for stellar wind for this exec (sets IMF y)
alfy_wind2 - final y alfven speed in sim units for stellar wind for this exec (sets IMF y)
alfz_wind1 - initial z alfven speed in sim units for stellar wind for this exec (sets IMF z)
alfz_wind2 - final z alfven speed in sim units for stellar wind for this exec (sets IMF z)
den_wind1 - initial  solar wind density in sim units for this exec
den_wind2 - final  solar wind density in sim units for this exec
reynolds - ion skin depths per grid point (start high = ideal MHD-ish, decr. 4 realism)
resist - ionospheric resistance scaling (the only ionospheric parameter we set)
o_conc - fraction of ‘oxygen’ species in solar wind (ratio of n_o to n_prot)
rho_frac - fraction of heavy ions in stellar wind
bfrac - fraction of tangential magnetic field allowed at planetary surface (not used)
vfrac - scale plasma corotation as function of distance during initialization
 
$lunar
details for the orbiting hi-res subgrid, extrapolate these from planetary vars or ask erika/robert, I didn’t use since I didn’t input physical moon in sub-grid
 
$physical 
re_equiv - spacing between grid points, for innermost grid, in units of planetary radii
b_equiv - magnetic field normalization factor (obtained from next two norms + alfven speed def)
v_equiv - velocity normalization factor (chosen to keep calculations near unity). This is 1 sim unit in km/s
rho_equiv - density normalization factor (chosen to keep calculations near unity). This is 1 sim unit in rmassq/cm^3
spacecraft - conditional to record spacecraft data (unused in my code MAT)
craft_input - conditional to use spacecraft data as physical inputs for sims
warp - conditional related to spacecraft data (unused in my code)
ut_init - time to set UT to for an initialization run (when start=.true.)
utstart - UT to resume at
 
$smooth
	these parameters used in introducing artificial viscosity to stabilize the simulation
	don’t change unless you are certain of what you intend
 
[Removed (MJS):]
$randomlinesofnumbers
grid sizes in grid points, each successive grid needs to be 2X the last, not necessarily centered

