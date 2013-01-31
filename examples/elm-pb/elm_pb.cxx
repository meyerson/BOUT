/*******************************************************************************
 * High-Beta Flute-Reduced MHD
 * see elm_reduced.pdf
 * Basically the same as Hazeltine-Meiss but different normalisations.
 * Can also include the Vpar compressional term
 *******************************************************************************/

#include "bout.hxx"
#include "initialprofiles.hxx"
#include "invert_laplace.hxx"
#include "interpolation.hxx"
#include "derivs.hxx"
#include <math.h>
#include "sourcex.hxx"
#include <boutmain.hxx>
#include <bout/constants.hxx>
#include <msg_stack.hxx>


//xia:       rf waves
BoutReal Prf; // power of rf waves 
Field2D Frf, Prf_h; // force of rf waves and the normalized power
BoutReal n_par; // parallel refractive index of rf
BoutReal rf_width_x,rf_width_y,rf_x,rf_y; // the width and the center coordinate of Gaussian rf waves 



// 2D inital profiles
Field2D  J0, P0; // Current and pressure
Vector2D b0xcv; // Curvature term
Field2D beta, gradparB;   // Used for Vpar terms
Field2D phi0;   // When diamagnetic terms used

// B field vectors
Vector2D B0vec; // B0 field vector

// V0 field vectors
Vector2D V0vec; // V0 field vector in convection
Vector2D V0eff; // effective V0 field vector in Ohm's law

// 3D evolving variables
Field3D U, Psi, P, Vpar;

// Derived 3D variables
Field3D Jpar, phi; // Parallel current, electric potential

Field3D Jpar2; //  Delp2 of Parallel current

// Constraint
Field3D C_phi;

// Parameters
BoutReal density; // Number density [m^-3]
BoutReal Bbar, Lbar, Tbar, Va; // Normalisation constants
BoutReal dnorm; // For diamagnetic terms: 1 / (2. * wci * Tbar)
BoutReal dia_fact; // Multiply diamagnetic term by this
BoutReal delta_i; // Normalized ion skin depth

BoutReal diffusion_par;  // Parallel pressure diffusion
BoutReal diffusion_Vpar;  //xia: Parallel velocitye diffusion
BoutReal diffusion_psiP;  //xia: Parallel magnetic flux diffusion
BoutReal diffusion_p4;   //M: 4th Parallel pressure diffusion
BoutReal diffusion_v4;   //M: 4th Parallel velocitye diffusion
BoutReal diffusion_U4;  //xia: 4th order Parall vorticity diffusion

BoutReal heating_P;  // heating power in pressure
BoutReal hp_width;  // heating profile radial width in pressure
BoutReal hp_length;  // heating radial domain in pressure
BoutReal sink_P;     // sink in pressure
BoutReal sp_width;   // sink profile radial width in pressure
BoutReal sp_length;  // sink radial domain in pressure

BoutReal sink_Ul;     // left edge sink in vorticity
BoutReal su_widthl;   // left edge sink profile radial width in vorticity
BoutReal su_lengthl;  // left edge sink radial domain in vorticity

BoutReal sink_Ur;     // right edge sink in vorticity
BoutReal su_widthr;   // right edge sink profile radial width in vorticity
BoutReal su_lengthr;  // right edge sink radial domain in vorticity

BoutReal viscos_par;  // Parallel viscosity
BoutReal viscos_perp; // Perpendicular viscosity
BoutReal hyperviscos; // Hyper-viscosity (radial)
Field3D hyper_mu_x; // Hyper-viscosity coefficient

// options
bool include_curvature, include_jpar0, compress0;
bool evolve_pressure;

BoutReal vacuum_pressure;
BoutReal vacuum_trans; // Transition width
Field3D vac_mask;

int phi_flags, apar_flags;
bool nonlinear;
bool evolve_jpar; 
BoutReal g; // Only if compressible
bool phi_curv;

bool diamag;
bool diamag_grad_t; // Grad_par(Te) term in Psi equation
bool diamag_phi0;   // Include the diamagnetic equilibrium phi0

bool eHall;
BoutReal AA; // ion mass in units of the proton mass; AA=Mi/Mp

BoutReal Vt0; // equilibrium toroidal flow normalized to Alfven velocity
BoutReal Vp0; // equilibrium poloidal flow normalized to Alfven velocity

bool nogradparj;
bool filter_z;
int filter_z_mode;
int low_pass_z;
int zonal_flow;
int zonal_field;
int zonal_bkgd;
bool relax_j_vac;
BoutReal relax_j_tconst; // Time-constant for j relax
Field3D Psitarget;   // The (moving) target to relax to

bool smooth_j_x;  // Smooth Jpar in the x direction

int jpar_bndry_width; // Zero jpar in a boundary region

bool parallel_lr_diff; // Use left and right shifted stencils for parallel differences

bool parallel_lagrange; // Use (semi-) Lagrangian method for parallel derivatives
bool parallel_project;  // Use Apar to project field-lines



//********************


Field3D Xip_x, Xip_z;     // Displacement of y+1 (in cell index space)

Field3D Xim_x, Xim_z;     // Displacement of y-1 (in cell index space)

bool phi_constraint; // Solver for phi using a solver constraint 

bool include_rmp; // Include RMP coil perturbation
bool simple_rmp;  // Just use a simple form for the perturbation
int rmp_n, rmp_m; // toroidal and poloidal mode numbers
BoutReal rmp_polwid;  // Poloidal width (-ve -> full, fraction of 2pi)
BoutReal rmp_polpeak; // Peak poloidal location (fraction of 2pi)
BoutReal rmp_factor;  // Multiply amplitude by this factor
BoutReal rmp_ramp;    // Ramp-up time for RMP [s]. negative -> instant
BoutReal rmp_freq;    // Amplitude oscillation frequency [Hz] (negative -> no oscillation)
BoutReal rmp_rotate;  // Rotation rate [Hz]
Field3D rmp_Psi0; // Parallel vector potential from Resonant Magnetic Perturbation (RMP) coils
Field3D rmp_Psi;  // Value used in calculations 
Field3D rmp_dApdt; // Time variation

BoutReal vac_lund, core_lund;       // Lundquist number S = (Tau_R / Tau_A). -ve -> infty
BoutReal vac_resist,  core_resist;  // The resistivities (just 1 / S)
Field3D eta;                    // Resistivity profile (1 / S)
bool spitzer_resist;  // Use Spitzer formula for resistivity
BoutReal Zeff;            // Z effective for resistivity formula

BoutReal hyperresist;    // Hyper-resistivity coefficient (in core only)
BoutReal ehyperviscos;   // electron Hyper-viscosity coefficient
Field3D hyper_eta_x; // Radial resistivity profile
Field3D hyper_eta_z; // Toroidal resistivity profile

int damp_width;     // Width of inner damped region
BoutReal damp_t_const;  // Timescale of damping

// Metric coefficients
Field2D Rxy, Bpxy, Btxy, B0, hthe;
Field2D I; // Shear factor

const BoutReal MU0 = 4.0e-7*PI;
const BoutReal Mi = 2.0*1.6726e-27; // Ion mass

// Communication objects
FieldGroup comms;

int jacobian(BoutReal t); // Jacobian-vector multiply

int precon_phi(BoutReal t, BoutReal cj, BoutReal delta); // Preconditioner with phi constraint

void advect_tracer(const Field3D &p,  // phi (input)
		   const Field3D &delta_x, const Field3D &delta_z, // Current location (input) 
		   Field3D &F_dx, Field3D &F_dz); // Time-derivative of location

const Field2D rfpower_exp(BoutReal rf_width_x, BoutReal rf_width_y, BoutReal rf_x,  BoutReal rf_y);
//xia: inital normalized profile of rf power with Gaussian form

const Field3D Grad2_par2new(const Field3D &f); //for 4th order diffusion


const Field3D Grad2_par2new(const Field3D &f)
{
  /*
   * This function implements d2/dy2 where y is the poloidal coordinate theta
   */


#ifdef CHECK
  int msg_pos = msg_stack.push("Grad2_par2new( Field3D )");
#endif



  Field3D result = D2DY2(f);
  

#ifdef TRACK
  result.name = "Grad2_par2new("+f.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}

const Field2D rfpower_exp(BoutReal rf_width_x, BoutReal rf_width_y, BoutReal rf_x,  BoutReal rf_y)
{
  Field2D result;
  result.allocate();

  if(rf_y < 0.0 ){
    for(int jx=0;jx<mesh->ngx;++jx)
      for(int jy=0;jy<mesh->ngy;++jy)
	//      for(int jz=0;jz<mesh->ngz;jz++)
	{
	  BoutReal lx = mesh->GlobalX(jx) - rf_x;
	  //BoutReal ly = mesh->GlobalX(jy) - rf_y;
	  //	BoutReal lz = mesh->GlobalX(jz) - rf_z;
	  //	BoutReal dampl = exp(-(lx*lx+ly*ly)/(rf_width*rf_width));
	  BoutReal dampl = exp(-(lx*lx)/(rf_width_x*rf_width_x));
	  result[jx][jy] = dampl/(sqrt(PI)*rf_width_x*rf_x+rf_width_x*rf_width_x*(1.0-1.0/2.71828)/2.0);
	}
  }
  else {
    for(int jx=0;jx<mesh->ngx;++jx)
      for(int jy=0;jy<mesh->ngy;++jy)
	//      for(int jz=0;jz<mesh->ngz;jz++)
	{
	  BoutReal lx = mesh->GlobalX(jx) - rf_x;
	  BoutReal ly = mesh->GlobalY(jy) - rf_y;
	  //	BoutReal lz = mesh->GlobalX(jz) - rf_z;
	  BoutReal dampl = exp(-(lx*lx)/(rf_width_x*rf_width_x)-(ly*ly)/(rf_width_y*rf_width_y));
	  //BoutReal dampl = exp(-(lx*lx)/(rf_width*rf_width));
	  result[jx][jy] = dampl/(PI*rf_width_x*rf_width_x*rf_x+sqrt(PI)*rf_width_x*rf_width_x*rf_width_y*(1.0-1.0/2.71828)/2.0);
	}
  }

  mesh->communicate(result);

  return result;
}


int physics_init(bool restarting)
{
  bool noshear;
  
  output.write("Solving high-beta flute reduced equations\n");
  output.write("\tFile    : %s\n", __FILE__);
  output.write("\tCompiled: %s at %s\n", __DATE__, __TIME__);

  //////////////////////////////////////////////////////////////
  // Load data from the grid

  // Load 2D profiles
  mesh->get(J0, "Jpar0");    // A / m^2
  mesh->get(P0, "pressure"); // Pascals
  
  // Load curvature term
  b0xcv.covariant = false; // Read contravariant components
  mesh->get(b0xcv, "bxcv"); // mixed units x: T y: m^-2 z: m^-2

  // Load metrics
  if(mesh->get(Rxy,  "Rxy")) { // m
    output.write("Error: Cannot read Rxy from grid\n");
    return 1;
  }
  if(mesh->get(Bpxy, "Bpxy")) { // T
    output.write("Error: Cannot read Bpxy from grid\n");
    return 1;
  }
  mesh->get(Btxy, "Btxy"); // T
  mesh->get(B0,   "Bxy");  // T
  mesh->get(hthe, "hthe"); // m
  mesh->get(I,    "sinty");// m^-2 T^-1

  //////////////////////////////////////////////////////////////
  // Read parameters from the options file
  // 
  // Options.get ( NAME,    VARIABLE,    DEFAULT VALUE)
  //
  // or if NAME = "VARIABLE" then just
  //
  // OPTION(VARIABLE, DEFAULT VALUE) 
  //
  // Prints out what values are assigned
  /////////////////////////////////////////////////////////////

  //  options.setSection("highbeta");    // The section header for all these settings

  Options *globalOptions = Options::getRoot();
  Options *options = globalOptions->getSection("highbeta");

  // xia: rf waves input
  OPTION(options, Prf,               0);      // Input rf power
  OPTION(options, n_par,             0);      // the parallel refractive index
  OPTION(options, rf_width_x,           0.1);    // the width of rf waves of x
  OPTION(options, rf_width_y,           0.1);    // the width of rf waves of y
  OPTION(options, rf_x,              0.8);   // The center of rf injection
  OPTION(options, rf_y,              -1.);

  BoutReal scaleP;
  OPTION(options, scaleP,            1.0);  // Multiply P0 by a number
  P0 *= scaleP;              // NOTE: Not self consistent with B field
  
  OPTION(options, density,           1.0e19); // Number density [m^-3]

  OPTION(options, evolve_jpar,       false);  // If true, evolve J raher than Psi
  OPTION(options, phi_constraint,    false);  // Use solver constraint for phi

  // Effects to include/exclude
  OPTION(options, include_curvature, true);
  OPTION(options, include_jpar0,     true);
  OPTION(options, evolve_pressure,   true);
  
  OPTION(options, compress0,          false);
  OPTION(options, nonlinear,         false);

  OPTION(options, eHall,            false);  // electron Hall or electron parallel pressue gradient effects? 
  OPTION(options, AA,               1.0);    // ion mass in units of proton mass

  OPTION(options, diamag,            false);  // Diamagnetic effects? 
  OPTION(options, diamag_grad_t,     diamag); // Grad_par(Te) term in Psi equation
  OPTION(options, diamag_phi0,       diamag); // Include equilibrium phi0
  OPTION(options, dia_fact,          1.0);    // Scale diamagnetic effects by this factor

  OPTION(options, Vt0,          0.0);    // equilibrium toroidal flow normalized to Alfven velocity
  OPTION(options, Vp0,          0.0);    // equilibrium poloidal flow normalized to Alfven velocity

  OPTION(options, noshear,           false);

  OPTION(options, relax_j_vac,       false); // Relax vacuum current to zero
  OPTION(options, relax_j_tconst,    0.1);
  
  // Toroidal filtering
  OPTION(options, filter_z,          false);  // Filter a single n
  OPTION(options, filter_z_mode,     1);
  OPTION(options, low_pass_z,       -1);      // Low-pass filter
  OPTION(options, zonal_flow,       -1);      // zonal flow filter
  OPTION(options, zonal_field,      -1);      // zonal field filter
  OPTION(options, zonal_bkgd,       -1);      // zonal background P filter

  // Radial smoothing
  OPTION(options, smooth_j_x,       false);  // Smooth Jpar in x

  // Jpar boundary region
  OPTION(options, jpar_bndry_width, -1);

  // Parallel differencing
  OPTION(options, parallel_lr_diff, false);
  OPTION(options, parallel_lagrange, false); // Use a (semi-) Lagrangian method for Grad_parP
  OPTION(options, parallel_project, false);

  // RMP-related options
  OPTION(options, include_rmp,       false);  // Read RMP data from grid
  
  OPTION(options, simple_rmp,        false);  // Include a simple RMP model
  OPTION(options, rmp_factor,         1.0);
  OPTION(options, rmp_ramp,          -1.0);
  OPTION(options, rmp_freq,          -1.0);
  OPTION(options, rmp_rotate,         0.0);

  // Vacuum region control
  OPTION(options, vacuum_pressure,   0.02);   // Fraction of peak pressure
  OPTION(options, vacuum_trans,      0.005);  // Transition width in pressure
  
  // Resistivity and hyper-resistivity options
  OPTION(options, vac_lund,          0.0);    // Lundquist number in vacuum region
  OPTION(options, core_lund,         0.0);    // Lundquist number in core region
  OPTION(options, hyperresist,       -1.0);
  OPTION(options, ehyperviscos,      -1.0);
  OPTION(options, spitzer_resist,    false);  // Use Spitzer resistivity
  OPTION(options, Zeff,              2.0);    // Z effective

  // Inner boundary damping
  OPTION(options, damp_width,        0);
  OPTION(options, damp_t_const,      0.1);

  // Viscosity and hyper-viscosity
  OPTION(options, viscos_par,        -1.0);  // Parallel viscosity
  OPTION(options, viscos_perp,       -1.0);  // Perpendicular viscosity
  OPTION(options, hyperviscos,       -1.0);  // Radial hyperviscosity
  
  // parallel pressure diffusion
  OPTION(options, diffusion_par,        -1.0);  // Parallel pressure diffusion
  OPTION(options, diffusion_Vpar,        -1.0);  //xia: Parallel velocity diffusion
  OPTION(options, diffusion_psiP,        -1.0);  //xia: Parallel magnetic flux diffusion
  OPTION(options, diffusion_p4,        -1.0);  // M: 4th Parallel pressure diffusion
  OPTION(options, diffusion_v4,        -1.0);  //M: 4th Parallel velocity diffusion
  OPTION(options, diffusion_U4,        -1.0);  //M: 4th Parallel vorticity diffusion


  // heating factor in pressure
  OPTION(options, heating_P,        -1.0);  //  heating power in pressure
  OPTION(options, hp_width,         0.1);  //  the percentage of radial grid points for heating profile radial width in pressure
  OPTION(options, hp_length,        0.04);  //  the percentage of radial grid points for heating profile radial domain in pressure

  // sink factor in pressure
  OPTION(options, sink_P,           -1.0);  //  sink in pressure
  OPTION(options, sp_width,         0.05);  //  the percentage of radial grid points for sink profile radial width in pressure
  OPTION(options, sp_length,        0.04);  //  the percentage of radial grid points for sink profile radial domain in pressure


  // left edge sink factor in vorticity
  OPTION(options, sink_Ul,           -1.0);  //  left edge sink in vorticity
  OPTION(options, su_widthl,         0.06);  //  the percentage of left edge radial grid points for sink profile radial width in vorticity
  OPTION(options, su_lengthl,        0.15);  //  the percentage of left edge radial grid points for sink profile radial domain in vorticity

  // right edge sink factor in vorticity
  OPTION(options, sink_Ur,           -1.0);  //  right edge sink in vorticity
  OPTION(options, su_widthr,         0.06);  //  the percentage of right edge radial grid points for sink profile radial width in vorticity
  OPTION(options, su_lengthr,        0.15);  //  the percentage of right edge radial grid points for sink profile radial domain in vorticity

  // Compressional terms
  OPTION(options, phi_curv,          true);
  options->get("gamma",             g,                 5.0/3.0);
  
  // Field inversion flags
  OPTION(options, phi_flags,         0);
  OPTION(options, apar_flags,        0);

  if(simple_rmp)
    include_rmp = true;

  if(include_rmp) {
    // Including external field coils. 
    if(simple_rmp) {
      // Use a fairly simple form for the perturbation

      Field2D pol_angle;
      if(mesh->get(pol_angle, "pol_angle")) {
	output.write("     ***WARNING: need poloidal angle for simple RMP\n");
	include_rmp = false;
      }else {
	OPTION(options, rmp_n,  3);
	OPTION(options, rmp_m,  9);
	OPTION(options, rmp_polwid, -1.0);
	OPTION(options, rmp_polpeak, 0.5);
	// Divide n by the size of the domain
        int zperiod;
	//     options->get(NULL, "zperiod", zperiod, 1);
	options->get("zperiod", zperiod, 1);
	if((rmp_n % zperiod) != 0)
	  output.write("     ***WARNING: rmp_n (%d) not a multiple of zperiod (%d)\n", rmp_n, zperiod);

	output.write("\tMagnetic perturbation: n = %d, m = %d, magnitude %e Tm\n", 
		     rmp_n, rmp_m, rmp_factor);
	
	rmp_Psi0 = 0.0;
	BoutReal ***d = rmp_Psi0.getData();
	// Set the outer boundary
	for(int jx=mesh->ngx-4;jx<mesh->ngx;jx++)
	  for(int jy=0;jy<mesh->ngy;jy++)
	    for(int jz=0;jz<mesh->ngz;jz++) {
	      
	      BoutReal angle = rmp_m * pol_angle[jx][jy] + rmp_n * ((BoutReal) jz) * mesh->dz;
	      d[jx][jy][jz] = (((BoutReal)(jx - 4)) / ((BoutReal)(mesh->ngx - 5))) * rmp_factor * cos(angle);
              if(rmp_polwid > 0.0) {
                // Multiply by a Gaussian in poloidal angle
                BoutReal gx = ((pol_angle[jx][jy] / (2.*PI)) - rmp_polpeak) / rmp_polwid;
                d[jx][jy][jz] *= exp(-gx*gx);
              }
	    }
	
	rmp_Psi0 = rmp_Psi0.shiftZ(false); // Shift into field-aligned coords

	// Now have a simple model for Psi due to coils at the outer boundary
	// Need to calculate Psi inside the domain, enforcing j = 0
	
	Jpar = 0.0;
	rmp_Psi0 = invert_laplace(Jpar, INVERT_4TH_ORDER | INVERT_OUT_SET, NULL);
	// Currently only implemented for 4th-order serial inversion (NXPE = 1)
      }
    }else {
      // Load perturbation from grid file.
      include_rmp = mesh->get(rmp_Psi0, "rmp_A"); // Only include if found
      
      // Multiply by factor
      rmp_Psi0 *= rmp_factor;
    }
  }

  if(!include_curvature)
    b0xcv = 0.0;
  
  if(!include_jpar0)
    J0 = 0.0;

  if(noshear) {
    if(include_curvature)
      b0xcv.z += I*b0xcv.x;
    mesh->ShiftXderivs = false;
    I = 0.0;
  }
  
  //////////////////////////////////////////////////////////////
  // SHIFTED RADIAL COORDINATES

  if(mesh->ShiftXderivs) {
    if(mesh->IncIntShear) {
      // BOUT-06 style, using d/dx = d/dpsi + I * d/dz
      mesh->IntShiftTorsion = I;
      
    }else {
      // Dimits style, using local coordinate system
      if(include_curvature)
	b0xcv.z += I*b0xcv.x;
      I = 0.0;  // I disappears from metric
    }
  }

  //////////////////////////////////////////////////////////////
  // NORMALISE QUANTITIES
  
  if(mesh->get(Bbar, "bmag")) // Typical magnetic field
    Bbar = 1.0;
  if(mesh->get(Lbar, "rmag")) // Typical length scale
    Lbar = 1.0;
  
  Va = sqrt(Bbar*Bbar / (MU0*density*Mi));

  Tbar = Lbar / Va;

  dnorm = dia_fact * Mi / (2.*1.602e-19*Bbar*Tbar);

  delta_i = AA*60.67*5.31e5/sqrt(density/1e6)/(Lbar*100.0);

  output.write("Normalisations: Bbar = %e T   Lbar = %e m\n", Bbar, Lbar);
  output.write("                Va = %e m/s   Tbar = %e s\n", Va, Tbar);
  output.write("                dnorm = %e\n", dnorm);
  output.write("    Resistivity\n");

  if(eHall)
    output.write("                delta_i = %e   AA = %e \n", delta_i, AA);

  if(vac_lund > 0.0) {
    output.write("        Vacuum  Tau_R = %e s   eta = %e Ohm m\n", vac_lund * Tbar, 
		 MU0 * Lbar * Lbar / (vac_lund * Tbar));
    vac_resist = 1. / vac_lund;
  }else {
    output.write("        Vacuum  - Zero resistivity -\n");
    vac_resist = 0.0;
  }
  if(core_lund > 0.0) {
    output.write("        Core    Tau_R = %e s   eta = %e Ohm m\n", core_lund * Tbar,
		 MU0 * Lbar * Lbar / (core_lund * Tbar));
    core_resist = 1. / core_lund;
  }else {
    output.write("        Core    - Zero resistivity -\n");
    core_resist = 0.0;
  }

  if(hyperresist > 0.0) {
    output.write("    Hyper-resistivity coefficient: %e\n", hyperresist);
    dump.add(hyper_eta_x, "hyper_eta_x", 1);
    dump.add(hyper_eta_z, "hyper_eta_z", 1);
  }

  if(ehyperviscos > 0.0) {
    output.write("    electron Hyper-viscosity coefficient: %e\n", ehyperviscos);
  }

  if(hyperviscos > 0.0) {
    output.write("    Hyper-viscosity coefficient: %e\n", hyperviscos);
    dump.add(hyper_mu_x, "hyper_mu_x", 1);
  }

  if(diffusion_par > 0.0) {
    output.write("    diffusion_par: %e\n", diffusion_par);
    dump.add(diffusion_par, "diffusion_par", 1);
  }

  //M: 4th order diffusion of p
  if(diffusion_p4 > 0.0) {
    output.write("    diffusion_p4: %e\n", diffusion_p4);
    dump.add(diffusion_p4, "diffusion_p4", 1);
  }

  // xia: add the parallel velocity
  if(diffusion_Vpar > 0.0) {
    output.write("    diffusion_Vpar: %e\n", diffusion_Vpar);
    dump.add(diffusion_Vpar, "diffusion_Vpar", 1);
  }

  // M: 4th order parallel velocity diffusion
  if(diffusion_v4 > 0.0) {
    output.write("    diffusion_v4: %e\n", diffusion_v4);
    dump.add(diffusion_v4, "diffusion_v4", 1);
  }

  // xia: 4th order vorticity diffusion
  if(diffusion_U4 > 0.0) {
    output.write("    diffusion_U4: %e\n", diffusion_U4);
    dump.add(diffusion_v4, "diffusion_U4", 1);
  }


  // xia: add the magnetic flux diffusion
  if(diffusion_psiP > 0.0) {
    output.write("    diffusion_psiP: %e\n", diffusion_psiP);
    dump.add(diffusion_psiP, "diffusion_psiP", 1);
  }

  if(heating_P > 0.0) {
    output.write("    heating_P(watts): %e\n", heating_P);
    dump.add(heating_P, "heating_P", 1);

    output.write("    hp_width(%): %e\n",hp_width);
    dump.add(hp_width, "hp_width", 1);

    output.write("    hp_length(%): %e\n",hp_length);
    dump.add(hp_length, "hp_length", 1);
  }

  if(sink_P > 0.0) {
    output.write("    sink_P(rate): %e\n", sink_P);
    dump.add(sink_P, "sink_P", 1);

    output.write("    sp_width(%): %e\n",sp_width);
    dump.add(sp_width, "sp_width", 1);

    output.write("    sp_length(%): %e\n",sp_length);
    dump.add(sp_length, "sp_length", 1);
  }

  //xia:  rf initial
  if(Prf>0.0)
    {
      output.write("    rf_power(Watts): %e\n ", Prf);
      dump.add(Prf, "Pinput", 1);

      output.write("    n_parallel: %e\n ", n_par);
      dump.add(n_par, "n_par", 1);

      output.write("    rf_width_x(%): %e\n ", rf_width_x);
      dump.add(rf_width_x, "rf_width_x", 1);

      output.write("    rf_width_y(%): %e\n ", rf_width_y);
      dump.add(rf_width_y, "rf_width_y", 1);

      output.write("    rf_x(%): %e\n ", rf_x);
      dump.add(rf_x, "rf_x", 1);

      output.write("    rf_y(%): %e\n ", rf_y);
      dump.add(rf_y, "rf_y", 1);

      Prf_h = Prf*rfpower_exp(rf_width_x, rf_width_y, rf_x, rf_y)* 2.*MU0*Tbar/(Bbar*Bbar);
      dump.add(Prf_h, "Prf", 1);
      Frf = Lbar/Tbar *2.*Prf_h *n_par/Va;
      dump.add(Frf, "Frf", 1);

    }


  Field2D Te;
  Te = P0 / (2.0*density * 1.602e-19); // Temperature in eV

  J0 = - MU0*Lbar * J0 / B0;
  P0 = 2.0*MU0 * P0 / (Bbar*Bbar);

  b0xcv.x /= Bbar;
  b0xcv.y *= Lbar*Lbar;
  b0xcv.z *= Lbar*Lbar;

  Rxy  /= Lbar;
  Bpxy /= Bbar;
  Btxy /= Bbar;
  B0   /= Bbar;
  hthe /= Lbar;
  mesh->dx   /= Lbar*Lbar*Bbar;
  I    *= Lbar*Lbar*Bbar;
  
  BoutReal pnorm = max(P0, true); // Maximum over all processors
  
  vacuum_pressure *= pnorm; // Get pressure from fraction
  vacuum_trans *= pnorm;

  // Transitions from 0 in core to 1 in vacuum
  vac_mask = (1.0 - tanh( (P0 - vacuum_pressure) / vacuum_trans )) / 2.0;

  if(spitzer_resist) {
    // Use Spitzer resistivity 
    output.write("\tTemperature: %e -> %e [eV]\n", min(Te), max(Te));
    eta = 0.51*1.03e-4*Zeff*20.*(Te^(-1.5)); // eta in Ohm-m. NOTE: ln(Lambda) = 20
    output.write("\tSpitzer resistivity: %e -> %e [Ohm m]\n", min(eta), max(eta));
    eta /= MU0 * Va * Lbar;
    output.write("\t -> Lundquist %e -> %e\n", 1.0/max(eta), 1.0/min(eta));
  }else {
    // transition from 0 for large P0 to resistivity for small P0
    eta = core_resist + (vac_resist - core_resist) * vac_mask;
  }

  dump.add(eta, "eta", 0);

  if(include_rmp) {
    // Normalise RMP quantities

    rmp_Psi0 /= Bbar * Lbar;
    
    rmp_ramp /= Tbar;
    rmp_freq *= Tbar;
    rmp_rotate *= Tbar;
    
    rmp_Psi = rmp_Psi0;
    rmp_dApdt = 0.0;

    bool apar_changing = false;

    output.write("Including magnetic perturbation\n");
    if(rmp_ramp > 0.0) {
      output.write("\tRamping up over period t = %e (%e ms)\n", rmp_ramp, rmp_ramp*Tbar*1000.);
      apar_changing = true;
    }
    if(rmp_freq > 0.0) {
      output.write("\tOscillating with frequency f = %e (%e kHz)\n", rmp_freq, rmp_freq/Tbar/1000.);
      apar_changing = true;
    }
    if(rmp_rotate != 0.0) {
      output.write("\tRotating with a frequency f = %e (%e kHz)\n", rmp_rotate, rmp_rotate/Tbar/1000.);
      apar_changing = true;
    }

    if(apar_changing) {
      dump.add(rmp_Psi, "rmp_Psi", 1);
      dump.add(rmp_dApdt, "rmp_dApdt", 1);
    }else {
      dump.add(rmp_Psi, "rmp_Psi", 0);
    }
  }else
    rmp_Psi = 0.0;
	
  /**************** CALCULATE METRICS ******************/

  mesh->g11 = (Rxy*Bpxy)^2;
  mesh->g22 = 1.0 / (hthe^2);
  mesh->g33 = (I^2)*mesh->g11 + (B0^2)/mesh->g11;
  mesh->g12 = 0.0;
  mesh->g13 = -I*mesh->g11;
  mesh->g23 = -Btxy/(hthe*Bpxy*Rxy);
  
  mesh->J = hthe / Bpxy;
  mesh->Bxy = B0;
  
  mesh->g_11 = 1.0/mesh->g11 + ((I*Rxy)^2);
  mesh->g_22 = (B0*hthe/Bpxy)^2;
  mesh->g_33 = Rxy*Rxy;
  mesh->g_12 = Btxy*hthe*I*Rxy/Bpxy;
  mesh->g_13 = I*Rxy*Rxy;
  mesh->g_23 = Btxy*hthe*Rxy/Bpxy;
  
  mesh->geometry(); // Calculate quantities from metric tensor

  // Set B field vector
  
  B0vec.covariant = false;
  B0vec.x = 0.;
  B0vec.y = Bpxy / hthe;
  B0vec.z = 0.;

  // Set V0vec field vector
  
  V0vec.covariant = false;
  V0vec.x = 0.;
  V0vec.y = Vp0 / hthe;
  V0vec.z = Vt0 / Rxy;

  // Set V0eff field vector

  V0eff.covariant = false;
  V0eff.x = 0.;
  V0eff.y = -(Btxy/(B0*B0))*(Vp0*Btxy-Vt0*Bpxy) / hthe;
  V0eff.z =  (Bpxy/(B0*B0))*(Vp0*Btxy-Vt0*Bpxy) / Rxy;

  /**************** SET VARIABLE LOCATIONS *************/

  P.setLocation(CELL_CENTRE);
  U.setLocation(CELL_CENTRE);
  phi.setLocation(CELL_CENTRE);
  Psi.setLocation(CELL_YLOW);
  Jpar.setLocation(CELL_YLOW);
  Vpar.setLocation(CELL_YLOW);

  /**************** SET EVOLVING VARIABLES *************/

  // Tell BOUT which variables to evolve
  SOLVE_FOR(U);
  SOLVE_FOR(P);

  if(evolve_jpar) {
    output.write("Solving for jpar: Inverting to get Psi\n");
    SOLVE_FOR(Jpar);
    dump.add(Psi, "Psi", 1);
  }else {
    output.write("Solving for Psi, Differentiating to get jpar\n");
    SOLVE_FOR(Psi);
    dump.add(Jpar, "jpar", 1);
  }
  
  if(parallel_lagrange) {
    // Evolving the distortion of the flux surfaces (Ideal-MHD only!)
    
    bout_solve(Xip_x, "Xip_x");
    bout_solve(Xip_z, "Xip_z");
    
    bout_solve(Xim_x, "Xim_x");
    bout_solve(Xim_z, "Xim_z");
  }
  
  if(parallel_project) {
    // Add Xi to the dump file
    dump.add(Xip_x, "Xip_x", 1);
    dump.add(Xip_z, "Xip_z", 1);
    
    dump.add(Xim_x, "Xim_x", 1);
    dump.add(Xim_z, "Xim_z", 1);
  }

  if(compress0) {
    output.write("Including compression (Vpar) effects\n");
    
    SOLVE_FOR(Vpar);

    beta = B0*B0 / ( 0.5 + (B0*B0 / (g*P0)));
    gradparB = Grad_par(B0) / B0;

    output.write("Beta in range %e -> %e\n",
      min(beta), max(beta));
  }

  if(phi_constraint) {
    // Implicit Phi solve using IDA
    
    if(!bout_constrain(phi, C_phi, "phi")) {
      output.write("ERROR: Cannot constrain. Run again with phi_constraint=false\n");
      bout_error("Aborting.\n");
    }
    
    // Set preconditioner
    solver->setPrecon(precon_phi);

  }else {
    // Phi solved in RHS (explicitly)
    dump.add(phi, "phi", 1);

    // Set Jacobian
    solver->setJacobian(jacobian);
  }

  // Diamagnetic phi0
  if(diamag && diamag_phi0) {
    // Stationary equilibrium plasma. ExB velocity balances diamagnetic drift
    phi0 = -0.5*dnorm*P0/B0;
    SAVE_ONCE(phi0);
  }

  // Add some equilibrium quantities and normalisations
  // everything needed to recover physical units
  SAVE_ONCE2(J0, P0);
  SAVE_ONCE4(density, Lbar, Bbar, Tbar);
  SAVE_ONCE2(Va, B0);

  /////////////// CHECK VACUUM ///////////////////////
  // In vacuum region, initial vorticity should equal zero
  
  if(!restarting) {
    // Only if not restarting: Check initial perturbation

    // Set U to zero where P0 < vacuum_pressure
    U = where(P0 - vacuum_pressure, U, 0.0);

    // Phi should be consistent with U
    phi = invert_laplace(U, phi_flags, NULL);
    
    if(diamag) {
      phi -= 0.5*dnorm * P / B0;
    }
  }

  /************** SETUP COMMUNICATIONS **************/
  
  comms.add(U);
  comms.add(Psi);
  comms.add(P);
  comms.add(phi);

  phi.setBoundary("phi"); // Set boundary conditions
  
  if(evolve_jpar) {
    comms.add(Jpar);
  }else {
    // otherwise Need to communicate Jpar separately
    Jpar.setBoundary("J");
  }
  Jpar2.setBoundary("J");

  return 0;
}


/*const Field3D Div_par_CtoL(const Field3D &var)
{
  return mesh->Bxy * Grad_par_CtoL(var / mesh->Bxy);
  }*/

// Parallel gradient along perturbed field-line
const Field3D Grad_parP(const Field3D &f, CELL_LOC loc = CELL_DEFAULT)
{
  Field3D result;
  
  if(parallel_lagrange || parallel_project) {
    // Moving stencil locations
    
    Field3D fp, fm; // Interpolated on + and - y locations
    
    fp = interpolate(f, Xip_x, Xip_z);
    fm = interpolate(f, Xim_x, Xim_z);
    
    result.allocate();
    for(int i=0;i<mesh->ngx;i++)
      for(int j=1;j<mesh->ngy-1;j++)
	for(int k=0;k<mesh->ngz-1;k++) {
	  result[i][j][k] = (fp[i][j+1][k] - fm[i][j-1][k])/(2.*mesh->dy[i][j]*sqrt(mesh->g_22[i][j]));
	}
  }else {
    if(parallel_lr_diff) {
      // Use left/right biased stencils. NOTE: First order only!
      if(loc == CELL_YLOW) {
	result = Grad_par_CtoL(f);
      }else
	result = Grad_par_LtoC(f);
    }else
      result = Grad_par(f, loc);
    
    if(nonlinear) {
      result -= b0xGrad_dot_Grad(Psi + rmp_Psi, f);
    }
  }

  return result;
}

bool first_run = true; // For printing out some diagnostics first time around

int physics_run(BoutReal t)
{
  ////////////////////////////////////////////
  // Transitions from 0 in core to 1 in vacuum
  if(nonlinear) {
    vac_mask = (1.0 - tanh( ((P0 + P) - vacuum_pressure) / vacuum_trans )) / 2.0;
    
    // Update resistivity
    if(spitzer_resist) {
      // Use Spitzer formula
      Field3D Te;
      Te = (P0+P)*Bbar*Bbar/(4.*MU0) / (density * 1.602e-19); // eV
      eta = 0.51*1.03e-4*Zeff*20.*(Te^(-1.5)); // eta in Ohm-m. ln(Lambda) = 20
      eta /= MU0 * Va * Lbar; // Normalised eta
    }else {
      // Use specified core and vacuum Lundquist numbers
      eta = core_resist + (vac_resist - core_resist) * vac_mask;
    }
  }

  ////////////////////////////////////////////
  // Resonant Magnetic Perturbation code

  if(include_rmp) {

    if( (rmp_ramp > 0.0) || (rmp_freq > 0.0) || (rmp_rotate != 0.0) ) {
      // Need to update the RMP terms
    
      if((rmp_ramp > 0.0) && (t < rmp_ramp)) {
	// Still in ramp phase
	
	rmp_Psi = (t / rmp_ramp) * rmp_Psi0 ; // Linear ramp
	
	rmp_dApdt = rmp_Psi0 / rmp_ramp;
      }else {
	rmp_Psi = rmp_Psi0;
	rmp_dApdt = 0.0;
      }
      
      if(rmp_freq > 0.0) {
	// Oscillating the amplitude
	
	rmp_dApdt = rmp_dApdt * sin(2.*PI * rmp_freq * t)
	  + rmp_Psi * (2.*PI * rmp_freq) * cos(2.*PI * rmp_freq * t);
	
	rmp_Psi *= sin(2.*PI * rmp_freq * t);
	
      }
      
      if(rmp_rotate != 0.0) {
	// Rotate toroidally at given frequency
	
	rmp_Psi = rmp_Psi.shiftZ(2*PI * rmp_rotate * t);
	
	rmp_dApdt = rmp_dApdt.shiftZ(2*PI * rmp_rotate * t);
	
	// Add toroidal rotation term. CHECK SIGN
	
	rmp_dApdt += DDZ(rmp_Psi) * 2*PI * rmp_rotate;
      }
      
      // Set to zero in the core
      rmp_Psi *= vac_mask;
    }else {
      // Set to zero in the core region
      rmp_Psi = rmp_Psi0 * vac_mask;  // Only in vacuum -> skin current -> diffuses inwards
    }
  }
  
  ////////////////////////////////////////////
  // Inversion

  if(evolve_jpar) {
    // Invert laplacian for Psi
    Psi = invert_laplace(Jpar, apar_flags, NULL);
  }

  if(phi_constraint) {
    // Phi being solved as a constraint
    
    Field3D Ctmp = phi;
    Ctmp.setBoundary("phi"); // Look up boundary conditions for phi
    Ctmp.applyBoundary();
    Ctmp -= phi; // Now contains error in the boundary
    
    C_phi = Delp2(phi) - U; // Error in the bulk
    C_phi.setBoundaryTo(Ctmp);

  }else {
    // Invert laplacian for phi
    phi = invert_laplace(U, phi_flags, NULL);
    
    if(diamag) {
      phi -= 0.5*dnorm * P / B0;
    }
    
    // Apply a boundary condition on phi for target plates
    phi.applyBoundary();
  }

  // Perform communications
  mesh->communicate(comms);

  if(!evolve_jpar) {
    // Get J from Psi
    Jpar = Delp2(Psi + rmp_Psi);

    Jpar.applyBoundary();
    mesh->communicate(Jpar);

    if(jpar_bndry_width > 0) {
      // Zero j in boundary regions. Prevents vorticity drive
      // at the boundary
      
      for(int i=0;i<jpar_bndry_width;i++)
	for(int j=0;j<mesh->ngy;j++)
	  for(int k=0;k<mesh->ngz-1;k++) {
	    if(mesh->firstX())
	      Jpar[i][j][k] = 0.0;
	    if(mesh->lastX())
	      Jpar[mesh->ngx-1-i][j][k] = 0.0;
	  }
    }

    // Smooth j in x
    if(smooth_j_x)
      Jpar = smooth_x(Jpar);

    //xqx begin
    // Get Delp2(J) from J
    Jpar2 = Delp2(Jpar);

    Jpar2.applyBoundary();
    mesh->communicate(Jpar2);

    if(jpar_bndry_width > 0) {
      // Zero jpar2 in boundary regions. Prevents vorticity drive
      // at the boundary
      
      for(int i=0;i<jpar_bndry_width;i++)
	for(int j=0;j<mesh->ngy;j++)
	  for(int k=0;k<mesh->ngz-1;k++) {
	    if(mesh->firstX())
	      Jpar2[i][j][k] = 0.0;
	    if(mesh->lastX())
	      Jpar2[mesh->ngx-1-i][j][k] = 0.0;
	  }
    }
    //xqx end
  }
  ////////////////////////////////////////////////////
  // Project perturbed field-lines

  if(parallel_project) {
    Vector3D Btilde; // Perturbed field
    Btilde = Grad(Psi) ^ B0vec;
    // Want contravariant components: B^x is change in psi, B^z is change in toroidal angle
    Btilde.toContravariant();
    
    // Calculate displacements in index space
    Xip_x = 0.;
    Xip_z = 0.;
    Xim_x = 0.;
    Xim_z = 0.;
    Field2D sgx = sqrt(mesh->g_11); // 1/(R*Bp)
    Field2D sgz = sqrt(mesh->g_33); // 1/R
    for(int i=0;i<mesh->ngx;i++)
      for(int j=1;j<mesh->ngy-1;j++)
	for(int k=0;k<mesh->ngz-1;k++) {
	  BoutReal bxratio = (Btilde.x[i][j][k]*sgx[i][j])/B0[i][j]; // |Bx| / B
	  bxratio *= (hthe[i][j] * B0[i][j] / Bpxy[i][j]); // Multiply by parallel length
	  Xip_x[i][j+1][k] = bxratio / (sgx[i][j+1] * mesh->dx[i][j+1]); // Convert to psi, then index
	  Xim_x[i][j-1][k] = bxratio / (sgx[i][j-1] * mesh->dx[i][j-1]); 
	  
	  BoutReal bzratio = (Btilde.z[i][j][k]*sgz[i][j])/B0[i][j]; // |Bz| / B
	  bzratio *= (hthe[i][j] * B0[i][j] / Bpxy[i][j]); // Multiply by parallel length
	  Xip_z[i][j+1][k] = bzratio / (sgz[i][j+1] * mesh->dz); // Convert to psi, then index
	  Xim_z[i][j-1][k] = bzratio / (sgz[i][j-1] * mesh->dz); 
	}
  }

  ////////////////////////////////////////////////////
  // Parallel electric field

  if(evolve_jpar) {
    // Jpar
    
    ddt(Jpar) = -Grad_parP(B0*U, CELL_YLOW) / B0 + eta*Delp2(Jpar);

    /*
      ddt(Psi) = -Grad_parP(B0*phi, CELL_YLOW) / B0 + eta*Jpar;
      // Apply boundary condition on Psi
      apply_boundary(ddt(Psi), "Psi");
      mesh->communicate(ddt(psi));
      ddt(Jpar) = Delp2(ddt(Psi));
    */

    if(relax_j_vac) {
      // Make ddt(Jpar) relax to zero.
      
      ddt(Jpar) -= vac_mask * Jpar / relax_j_tconst;
    }
  }else {
    // Vector potential
      
    ddt(Psi) = -Grad_parP(B0*phi, CELL_YLOW) / B0 + Jpar*eta;
    //xqx      ddt(Psi) = -Grad_parP(B0*phi, CELL_YLOW) / B0 + eta*Jpar;

    if(eHall) {
      ddt(Psi) +=  0.25*delta_i*(Grad_parP(B0*P, CELL_YLOW) / B0 
                                 +b0xGrad_dot_Grad(P0, Psi));   // electron parallel pressure
    }
    if(diamag && diamag_phi0)
      ddt(Psi) -= b0xGrad_dot_Grad(B0*phi0, Psi)/B0;   // Equilibrium flow

    if(Vt0!=0 || Vp0!=0)
      ddt(Psi) += 1.0*V_dot_Grad(V0eff, Psi);

    //F_Psi = -Grad_par_CtoL(B0*phi) / B0 + eta*Jpar;
    
    if(diamag_grad_t) {
      // grad_par(T_e) correction

      ddt(Psi) += 1.71 * dnorm * 0.5 * Grad_parP(P, CELL_YLOW) / B0;
    }

    // Hyper-resistivity
    if(hyperresist > 0.0) {
      ddt(Psi) -= eta*hyperresist * Delp2(Jpar);
    }

    // electron Hyper-viscosity coefficient
    if(ehyperviscos > 0.0) {
      ddt(Psi) -= eta*ehyperviscos * Delp2(Jpar2);
    }

    //xia: adding hyper diffusion of psi
    if(diffusion_psiP > 0.0) {
      ddt(Psi) -= diffusion_psiP * Grad2_par2new(Grad2_par2new(Psi));
    }



    // Vacuum solution
    if(relax_j_vac) {
      // Calculate the J and Psi profile we're aiming for
      Field3D Jtarget = Jpar * (1.0 - vac_mask); // Zero in vacuum
      
      // Invert laplacian for Psi
      Psitarget = invert_laplace(Jtarget, apar_flags, NULL);
      
      // Add a relaxation term in the vacuum
      ddt(Psi) = ddt(Psi)*(1. - vac_mask) - (Psi - Psitarget)*vac_mask / relax_j_tconst;
    }

    //xia: rf-force
    if(Prf>0.0)
      {
	ddt(Psi) -= Bbar*Tbar/(2.*1836*MU0*B0*Lbar*Lbar)*Frf;
      }
  }

  if(parallel_lagrange) {
    // Move tracers with the plasma
    //advect_tracer(phi, Xi_x, Xi_z, F_Xi_x, F_Xi_z);

    // Follow intersection of fieldlines with neighbouring
    // poloidal planes

    // Calculate the ExB velocity on the grid-points
    Field3D vx0, vz0, vxp, vzp, vxm, vzm;
    vx0 = DDZ(phi);
    vz0 = -DDX(phi);
    
    // Interpolate to get velocity at intersections
    vxp = interpolate(vx0, Xip_x, Xip_z);
    vzp = interpolate(vz0, Xip_x, Xip_z);
    
    vxm = interpolate(vx0, Xim_x, Xim_z);
    vzm = interpolate(vz0, Xim_x, Xim_z);
    
    // Convert into relative velocity of
    // grid points and intersections.
    // NOTE: Should be Christoffel terms here
  }

  ////////////////////////////////////////////////////
  // Vorticity equation

  ddt(U) = (B0^2) * b0xGrad_dot_Grad(Psi + rmp_Psi, J0, CELL_CENTRE); // Grad j term

  ddt(U) += b0xcv*Grad(P);  // curvature term

  if(!nogradparj) {
    // Parallel current term
    ddt(U) -= (B0^2)*Grad_parP(Jpar, CELL_CENTRE); // b dot grad j
    //ddt(U) -= (B0^2)*Grad_par_LtoC(Jpar);
  }
  
  if(diamag && diamag_phi0)
    ddt(U) -= b0xGrad_dot_Grad(phi0, U);   // Equilibrium flow

  if(Vt0!=0 || Vp0!=0)
    ddt(U) -= V_dot_Grad(V0vec, U);

  if(nonlinear) {
    ddt(U) -= b0xGrad_dot_Grad(phi, U);    // Advection
  }

  // Viscosity terms 
  if(viscos_par > 0.0)
    ddt(U) += viscos_par * Grad2_par2(U); // Parallel viscosity
  
    //xia: 4th order diffusion of velocity
  if(diffusion_U4 > 0.0)
    ddt(U) -= diffusion_U4 * Grad2_par2new(Grad2_par2new(U));

  if(viscos_perp > 0.0)
    ddt(U) += viscos_perp * Delp2(U);     // Perpendicular viscosity

  // Hyper-viscosity
  if(hyperviscos > 0.0) {
    // Calculate coefficient.
    
    hyper_mu_x = hyperviscos * mesh->g_11*SQ(mesh->dx) * abs(mesh->g11*D2DX2(U)) / (abs(U) + 1e-3);
    hyper_mu_x.applyBoundary("dirichlet"); // Set to zero on all boundaries
    
    ddt(U) += hyper_mu_x * mesh->g11*D2DX2(U);
    
    if(first_run) { // Print out maximum values of viscosity used on this processor
      output.write("   Hyper-viscosity values:\n");
      output.write("      Max mu_x = %e, Max_DC mu_x = %e\n", max(hyper_mu_x), max(hyper_mu_x.DC()));
    }
  }

  // left edge sink terms 
  if(sink_Ul > 0.0){
    ddt(U) -=  sink_Ul*sink_tanhxl(P0,U,su_widthl,su_lengthl); // core sink
  }

  // right edge sink terms 
  if(sink_Ur > 0.0){
    ddt(U) -=  sink_Ur*sink_tanhxr(P0,U,su_widthr,su_lengthr); //  sol sink
  }

  ////////////////////////////////////////////////////
  // Pressure equation

  ddt(P) = 0.0;
  if(evolve_pressure) {
    ddt(P) -= b0xGrad_dot_Grad(phi, P0);
    
    if(diamag && diamag_phi0)
      ddt(P) -= b0xGrad_dot_Grad(phi0, P);   // Equilibrium flow

    if(Vt0!=0 || Vp0!=0)
      ddt(P) -= V_dot_Grad(V0vec, P);

    if(nonlinear)
      ddt(P) -= b0xGrad_dot_Grad(phi, P);
  }

  // Parallel diffusion terms 
  if(diffusion_par > 0.0)
    ddt(P) += diffusion_par * Grad2_par2(P); // Parallel diffusion

  //another choice
  //  if(diffusion_par > 0.0)
  //  ddt(P) += diffusion_par * Grad2_par2new(P)

  //M: 4th order Parallel diffusion terms 
  if(diffusion_p4 > 0.0)
    ddt(P) -= diffusion_p4 * Grad2_par2new(Grad2_par2new(P));

  // xia: rf power
  if(Prf>0.0)
    ddt(P) += 2./3.*Prf_h; 

  // heating source terms 
  if(heating_P > 1.0){
    BoutReal pnorm = P0[0][0];
    //    ddt(P) += 2./3.*Prf_h;  // xia:  rf waves 
    ddt(P) += heating_P*source_expx2(P0,2.*hp_width,0.5*hp_length)*(Tbar/pnorm); // heat source
    //ddt(P) += (100.*source_tanhx(P0,hp_width,hp_length)+0.01) * mesh->g11 * D2DX2(P) * (Tbar/Lbar/Lbar) ;     // radial diffusion
  }

  // sink terms 
  if(sink_P > 0.0){
    ddt(P) -= sink_P*sink_tanhxr(P0,P,sp_width,sp_length)*Tbar; // sink
  }

  ////////////////////////////////////////////////////
  // Compressional effects
  
  if(compress0) {
    
    ddt(P) += beta*( - Grad_parP(Vpar, CELL_CENTRE) + Vpar*gradparB );
    //    ddt(P) -= beta*Div_par_CtoL(Vpar);
    
    if(phi_curv) {
      ddt(P) -= 2.*beta*b0xcv*Grad(phi);
    }

    // Vpar equation

    ddt(Vpar) = -0.5*Grad_parP(P + P0, CELL_CENTRE);
    //    ddt(Vpar) = -0.5*Grad_parP(P + P0, CELL_YLOW);
    //    ddt(Vpar) = -0.5*Grad_par_LtoC(P + P0);

    if(nonlinear)
      ddt(Vpar) -= b0xGrad_dot_Grad(phi, Vpar); // Advection

    //xia: diffusion of velocity
    if(diffusion_Vpar > 0.0)
      ddt(Vpar) += diffusion_Vpar * Grad2_par2(Vpar); // Parallel diffusion

    //another choice
    //    if(diffusion_Vpar > 0.0)
    //  ddt(Vpar) += diffusion_Vpar * Grad2_par2new(Vpar);

    //M: 4th order diffusion of velocity
    if(diffusion_v4 > 0.0)
      ddt(Vpar) -= diffusion_v4 * Grad2_par2new(Grad2_par2new(Vpar));

    //xqx include vpar convection
    //#if 0
    if(nonlinear) {
      if(evolve_jpar) {
	ddt(Jpar) -= Vpar_Grad_par(Vpar, Jpar);
      }//else
	//	ddt(Psi) -= Vpar_Grad_par(Vpar, Psi);

      ddt(U) -= Vpar_Grad_par(Vpar, U);
      ddt(P) -= Vpar_Grad_par(Vpar, P+P0);
      ddt(Vpar) -= Vpar_Grad_par(Vpar, Vpar);
    }
    //#endif

    //xia: rf force
    if(Prf>0.0)
      ddt(Vpar) += Frf;
  }
  
  if(filter_z) {
    // Filter out all except filter_z_mode
    
    if(evolve_jpar) {
      ddt(Jpar) = filter(ddt(Jpar), filter_z_mode);
    }else
      ddt(Psi) = filter(ddt(Psi), filter_z_mode);

    ddt(U) = filter(ddt(U), filter_z_mode);
    ddt(P) = filter(ddt(P), filter_z_mode);

    if(compress0) {
      ddt(Vpar) = filter(ddt(Vpar), filter_z_mode);
    }
  }

  if(low_pass_z > 0) {
    // Low-pass filter, keeping n up to low_pass_z
    if(evolve_jpar) {
      ddt(Jpar) = lowPass(ddt(Jpar), low_pass_z, zonal_field);
    }else
      ddt(Psi) = lowPass(ddt(Psi), low_pass_z, zonal_field);

    ddt(U) = lowPass(ddt(U), low_pass_z, zonal_flow);
    ddt(P) = lowPass(ddt(P), low_pass_z, zonal_bkgd);

    if(compress0) {
      ddt(Vpar) = lowPass(ddt(Vpar), low_pass_z, zonal_bkgd);
    }
  }

  if(damp_width > 0) {
    for(int i=0;i<damp_width;i++) {
      for(int j=0;j<mesh->ngy;j++)
	for(int k=0;k<mesh->ngz;k++) {
	  if(mesh->firstX())
	    ddt(U)[i][j][k] -= U[i][j][k] / damp_t_const;
	  if(mesh->lastX())
	    ddt(U)[mesh->ngx-1-i][j][k] -= U[mesh->ngx-1-i][j][k] / damp_t_const;
	}
    }
  }

  first_run = false;

  return 0;
}


/*******************************************************************************
 * Preconditioner described in elm_reduced.cpp
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in F_vars
 * 
 * o Return values should be in vars (overwriting system state)
 *
 * NOTE: EXPERIMENTAL
 * enable by setting solver / use_precon = true in BOUT.inp
 *******************************************************************************/

const Field3D Lp(const Field3D &f)
{
  return b0xcv*Grad(f);
}

const Field3D Lpsi(const Field3D &f)
{
  Field3D jpre = Delp2(f);
  jpre.setBoundary("J");
  jpre.applyBoundary();

  mesh->communicate(jpre);

  Field3D result = b0xGrad_dot_Grad(f, J0);
  
  return (B0^2)*(result - Grad_parP(jpre) );
}

const Field3D Pschur(const Field3D &f)
{
  Field3D phitmp = invert_laplace(f, phitmp, phi_flags, NULL);
  // Need to communicate phi
  mesh->communicate(phitmp);

  Field3D dP = b0xGrad_dot_Grad(phitmp, P0);
  Field3D dPsi = Grad_parP(B0*phitmp) / B0;
  dP.setBoundary("P"); dP.applyBoundary();
  dPsi.setBoundary("Psi"); dPsi.applyBoundary();
  
  Field3D dJ = Delp2(dPsi);
  dJ.setBoundary("J"); dJ.applyBoundary();

  mesh->communicate(dP, dPsi, dJ);
  
  Field3D result = b0xcv*Grad(dP) + (B0^2)*(Grad_par(dJ) - b0xGrad_dot_Grad(dPsi, J0));
  result.setBoundary("U"); result.applyBoundary();
  
  return result;
}

/*******************************************************************************
 * Jacobian-vector multiply
 *
 * Input
 *   System state is in (P, Psi, U)
 *   Vector v is in (F_P, F_Psi, F_U)
 * Output
 *   Jacobian-vector multiplied Jv should be in (P, Psi, U)
 * 
 * NOTE: EXPERIMENTAL
 * enable by setting solver / use_jacobian = true in BOUT.inp
 *******************************************************************************/

int jacobian(BoutReal t)
{
  // NOTE: LINEAR ONLY!
  
  // Communicate
  mesh->communicate(ddt(P), ddt(Psi), ddt(U));

  phi = invert_laplace(ddt(U), phi_flags, NULL);

  Jpar = Delp2(ddt(Psi));

  mesh->communicate(phi, Jpar);

  Field3D JP = -b0xGrad_dot_Grad(phi, P0);
  JP.setBoundary("P"); JP.applyBoundary();
  
  Field3D JPsi = -Grad_par(B0*phi, CELL_YLOW) / B0;
  JPsi.setBoundary("Psi"); JPsi.applyBoundary();

  Field3D JU = b0xcv*Grad(ddt(P))
    - (B0^2)*Grad_par(Jpar, CELL_CENTRE)
    + (B0^2) * b0xGrad_dot_Grad(ddt(Psi), J0, CELL_CENTRE);
  JU.setBoundary("U"); JU.applyBoundary();

  // Put result into vars

  P = JP;
  Psi = JPsi;
  U = JU;

  return 0;
}

/*******************************************************************************
 * Preconditioner for when phi solved as a constraint
 * Currently only possible with the IDA solver
 *
 * o System state in variables (as in rhs function)
 * o Values to be inverted in F_vars
 * 
 * o Return values should be in vars (overwriting system state)
 *******************************************************************************/

int precon_phi(BoutReal t, BoutReal cj, BoutReal delta)
{
  P = ddt(P);
  Psi = ddt(Psi);
  /*
  invert_laplace(F_U, phi, phi_flags, NULL);
  phi = C_phi + phi;
  */
  phi = invert_laplace(C_phi - ddt(U), phi_flags, NULL);
  U = ddt(U);
  return 0;
}

/*******************************************************************************
 * Tracer advection code
 * 
 * May be used to track field-lines and calculate Grad_par along
 * highly perturbed field-lines.
 *
 * Aug 2009
 *******************************************************************************/

void advect_tracer(const Field3D &p,  // phi (input)
		   const Field3D &delta_x, const Field3D &delta_z, // Current location (input) 
		   Field3D &F_dx, Field3D &F_dz)
{
  // Calculate the ExB velocity
  Field3D vx, vz;
  
  vx = DDZ(p);
  vz = -DDX(p);

  // Interpolate these velocities onto the current locations
  vx = interpolate(vx, delta_x, delta_z);
  vz = interpolate(vz, delta_x, delta_z);
  
  // Time derivative of the cell index:
  F_dx = vx / mesh->dx;
  F_dz = vz / mesh->dz;
}
