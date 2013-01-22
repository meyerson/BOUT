/*******************************************************************
 * Bare Bones Bout
 *
 * D. Meyerson 2013
 *******************************************************************/

#include <bout.hxx>
#include <boutmain.hxx>
#include <derivs.hxx>
#include <initialprofiles.hxx>
#include <invert_laplace.hxx>

// Evolving variables 
Field3D u, n; //vorticity, density

//derived variables
Field3D phi,brkt;
int phi_flags;


//other params
BoutReal alpha, nu, mu,gam, beta;

FieldGroup comms; // Group of variables for communications

const Field3D mybracket(const Field3D &phi, const Field3D &A);

int physics_init(bool restarting)
{
  // // 2D initial profiles
  // Field2D N0, P0;
  // Vector2D V0;
  // BoutReal v0_multiply;

  // // Read initial conditions

  // mesh->get(N0, "density");
  // mesh->get(P0, "pressure");
  // V0.covariant = false; // Read contravariant components
  // V.covariant = false; // Evolve contravariant components
  // mesh->get(V0, "v");
  // g.covariant = false;
  // mesh->get(g, "g");
  Options *options = Options::getRoot();
  options = options->getSection("physics");
  OPTION(options, phi_flags, 0);
  OPTION(options, alpha,3e-5);
  OPTION(options, nu, 2e-3);
  //OPTION(options, mu, 0.040);
  OPTION(options, mu, 2e-3);
  OPTION(options, gam, 1);
  OPTION(options, beta, 6e-4);

  // read options
  phi.setLocation(CELL_CENTRE);
  
  //overide
  //beta = 1e-5;
  //mu = 1e-2;
  //nu = 1e-2;

  bout_solve(u, "u");
  comms.add(u);
  phi = invert_laplace(u, phi_flags); 
  comms.add(phi);
  bout_solve(n, "n");

  comms.add(n);
   
  //brkt = b0xGrad_dot_Grad(phi, u);

  dump.add(phi,"phi",1);
  dump.add(brkt,"brkt",1);

  return 0;
}

#define bracket3D(f, g) ( b0xGrad_dot_Grad(f, g) )

int physics_run(BoutReal t)
{
  // Run communications
  mesh->communicate(comms);
  //output.write("mesh->dx): %g \n",beta);
  phi = invert_laplace(u, phi_flags);
  phi.applyBoundary("neumann");
  //phi.applyBoundary("dirichlet");
  // Density
  //f = lowPass(f,1);
  //f = lowPass(g,1);
  mesh->communicate(comms);
  //mesh->communicate(phi);
  ddt(u)=0;
  ddt(n)=0;
 
  // brkt = mybracket(phi,n);
  // //brkt.applyBoundary("neumann");
  // brkt.applyBoundary("dirichlet");

  ddt(u) -= mybracket(phi,u);
  ddt(u) += alpha * phi;
  ddt(u) += nu * Laplacian(u);
  //ddt(u) -= beta * DDY(n)/n; 
  //ddt(u) -= beta* Grad_par(n)/n; 
  ddt(u) -= Grad_par(n); 
  //ddt(u).applyBoundary("dirichlet");

  //mesh->communicate(comms); no don't do this here
  //.applyBoundary();
  //brkt = VDDY(DDY(phi), n) +  VDDZ(DDZ(phi), n) ;
 

  ddt(n) -= mybracket(phi,n);
  ddt(n) += mu * Laplacian(n);
  ddt(n) -= alpha* n;
  //ddt(n).applyBoundary("dirichlet");

  //ddt(n) -= VDDZ(n,n) + mu*VDDY(u,n);

  // if (driver){
    
  // }

  

  //ddt(f) -= 10*f;
  //ddt(f) = lowPass(ddt(f),5);

  
 
  return 0;
}


const Field3D mybracket(const Field3D &phi, const Field3D &A)
{
  Field3D dpdx, dpdy, dpdz;
  Field3D vx, vy, vz;
  Field3D result;

  //output.write("mesh->Bxy = %e\n", (mesh->J*sqrt(mesh->g_22))[2][2]);

#ifdef CHECK
  int msg_pos = msg_stack.push("b0xGrad_dot_Grad( Field3D , Field3D )");
#endif

  // Calculate phi derivatives
  #pragma omp parallel sections
  {
    #pragma omp section
    dpdx = DDX(phi); 
    
    #pragma omp section
    dpdy = DDY(phi);
    
    #pragma omp section
    dpdz = 0;
  }
  
  // Calculate advection velocity
  #pragma omp parallel sections
  {

    #pragma omp section
    vx = mesh->g_23*dpdz - mesh->g_33*dpdy; 
    //vx = mesh->g_22*dpdz - mesh->g_23*dpdy;
    
    #pragma omp section
    vy = mesh->g_33*dpdx - mesh->g_13*dpdz;
      //vy = mesh->g_23*dpdx - mesh->g_12*dpdz;
    
    #pragma omp section
    //vz = mesh->g_12*dpdy - mesh->g_22*dpdx;
    vz =0;
  }


  // Upwind A using these velocities
  
  Field3D ry, rz;
  #pragma omp parallel sections
  {
    #pragma omp section
    result = VDDX(vx, A);
    
    #pragma omp section
    ry = VDDY(vy, A);
    
    #pragma omp section
    rz = VDDZ(vz, A);
  }
  //output.write("mesh->g_22: %g  \n" ,vx[4][4]);
  //result = (ry + rz); //usually  convention
  result = (result + ry) / (mesh->J*sqrt(mesh->g_22));

#ifdef TRACK
  result.name = "b0xGrad_dot_Grad("+phi.name+","+A.name+")";
#endif
#ifdef CHECK
  msg_stack.pop(msg_pos);
#endif
  return result;
}
