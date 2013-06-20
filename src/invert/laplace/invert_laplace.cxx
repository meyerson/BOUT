/*!
 * \file invert_laplace.cxx
 *
 * \brief Perpendicular Laplacian inversion using FFT and Tridiagonal solver
 *
 * Equation solved is \f$d*\nabla^2_\perp x + (1./c)\nabla_perp c\cdot\nabla_\perp x + a x = b \f$, where
 * \f$x\f$ and \f$x\f$ are perpendicular (X-Z) or 3D fields, 
 * and \f$a\f$ and d are 2D fields. If d is not supplied then it is 1
 * 
 * Flags control the boundary conditions (see header file)
 *
 **************************************************************************
 * Copyright 2010 B.D.Dudson, S.Farley, M.V.Umansky, X.Q.Xu
 *
 * Contact: Ben Dudson, bd512@york.ac.uk
 * 
 * This file is part of BOUT++.
 *
 * BOUT++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * BOUT++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with BOUT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <globals.hxx>
#include <invert_laplace.hxx>
#include <bout_types.hxx>
#include <options.hxx>
#include <boutexception.hxx>
#include <utils.hxx>
#include <cmath>
#include <bout/sys/timer.hxx>
#include <output.hxx>
#include <msg_stack.hxx>
#include <bout/constants.hxx>

#include "laplacefactory.hxx"

/**********************************************************************************
 *                         INITIALISATION AND CREATION
 **********************************************************************************/

/// Laplacian inversion initialisation. Called once at the start to get settings
Laplacian::Laplacian(Options *options) {

  if(options == NULL) {
    // Use the default options
    options = Options::getRoot()->getSection("laplace");
  }
  
  output.write("Initialising Laplacian inversion routines\n");
  
  // Communication option. Controls if asyncronous sends are used
  options->get("async", async_send, true);
  
  BoutReal filter; ///< Fraction of Z modes to filter out. Between 0 and 1
  OPTION(options, filter, 0.2);
  int ncz = mesh->ngz-1;
 
  //int ncz = mesh->ngy;
  // convert filtering into an integer number of modes
  maxmode = ROUND((1.0 - filter) * ((double) (ncz / 2)));
  // Can be overriden by max_mode option
  OPTION(options, maxmode, maxmode);
  if(maxmode < 0) maxmode = 0;
  if(maxmode > ncz/2) maxmode = ncz/2;
  
  OPTION3(options, low_mem, all_terms, nonuniform, false);
  
  OPTION(options, flags, 0);
}

Laplacian* Laplacian::create(Options *opts) {
  return LaplaceFactory::getInstance()->createLaplacian(opts);
}

Laplacian* Laplacian::instance = NULL;

Laplacian* Laplacian::defaultInstance() {
  if(instance == NULL)
    instance = create();
  return instance;
}

/**********************************************************************************
 *                                 Solve routines
 **********************************************************************************/

const Field3D Laplacian::solve(const Field3D &b) {
#ifdef CHECK
  msg_stack.push("Laplacian::solve(Field3D)");
#endif

  int ys = mesh->ystart, ye = mesh->yend;

  if(mesh->hasBndryLowerY())
    ys = 0; // Mesh contains a lower boundary
  if(mesh->hasBndryUpperY())
    ye = mesh->ngy-1; // Contains upper boundary

  Field3D x;
  x.allocate();

  for(int jy=ys; jy <= ye; jy++) {
    x = solve(b.slice(jy));
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  x.setLocation(b.getLocation());

  return x;
}

// NB: Really inefficient, but functional
const Field2D Laplacian::solve(const Field2D &b) {
  Field3D f = b;
  f = solve(f);
  return f.DC();
}

const Field3D Laplacian::solve(const Field3D &b, const Field3D &x0) {
  #ifdef CHECK
  msg_stack.push("Laplacian::solve(Field3D, Field3D)");
#endif

  int ys = mesh->ystart, ye = mesh->yend;
  if(mesh->hasBndryLowerY())
    ys = 0; // Mesh contains a lower boundary
  if(mesh->hasBndryUpperY())
    ye = mesh->ngy-1; // Contains upper boundary
  
  Field3D x;
  x.allocate();

  for(int jy=ys; jy <= ye; jy++) {
    x = solve(b.slice(jy), x0.slice(jy));
  }
  
#ifdef CHECK
  msg_stack.pop();
#endif

  x.setLocation(b.getLocation());

  return x;
}

const Field2D Laplacian::solve(const Field2D &b, const Field2D &x0) {
  Field3D f = b, g = x0;
  f = solve(f, g);
  return f.DC();
}

/**********************************************************************************
 *                              MATRIX ELEMENTS
 **********************************************************************************/

void Laplacian::tridagCoefs(int jx, int jy, int jz, 
                            dcomplex &a, dcomplex &b, dcomplex &c, 
                            const Field2D *ccoef, const Field2D *d) {
  
  BoutReal coef1, coef2, coef3, coef4, coef5, kwave;
  
  kwave=jz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
  
  coef1=mesh->g11[jx][jy];     ///< X 2nd derivative coefficient
  coef2=mesh->g33[jx][jy];     ///< Z 2nd derivative coefficient
  coef3=2.*mesh->g13[jx][jy];  ///< X-Z mixed derivative coefficient

  coef4 = 0.0;
  coef5 = 0.0;
  if(all_terms) {
    coef4 = mesh->G1[jx][jy]; // X 1st derivative
    coef5 = mesh->G3[jx][jy]; // Z 1st derivative
  }

  if(d != (Field2D*) NULL) {
    // Multiply Delp2 component by a factor
    coef1 *= (*d)[jx][jy];
    coef2 *= (*d)[jx][jy];
    coef3 *= (*d)[jx][jy];
    coef4 *= (*d)[jx][jy];
    coef5 *= (*d)[jx][jy];
  }

  if(nonuniform) {
    // non-uniform mesh correction
    if((jx != 0) && (jx != (mesh->ngx-1))) {
      //coef4 += mesh->g11[jx][jy]*0.25*( (1.0/dx[jx+1][jy]) - (1.0/dx[jx-1][jy]) )/dx[jx][jy]; // SHOULD BE THIS (?)
      coef4 -= 0.5*((mesh->dx[jx+1][jy] - mesh->dx[jx-1][jy])/SQ(mesh->dx[jx][jy]))*coef1; // BOUT-06 term
    }
  }

  if(ccoef != NULL) {
    // A first order derivative term
    
    if((jx > 0) && (jx < (mesh->ngx-1)))
      coef4 += mesh->g11[jx][jy] * ((*ccoef)[jx+1][jy] - (*ccoef)[jx-1][jy]) / (2.*mesh->dx[jx][jy]*((*ccoef)[jx][jy]));
  }
  
  if(mesh->ShiftXderivs && mesh->IncIntShear) {
    // d2dz2 term
    coef2 += mesh->g11[jx][jy] * mesh->IntShiftTorsion[jx][jy] * mesh->IntShiftTorsion[jx][jy];
    // Mixed derivative
    coef3 = 0.0; // This cancels out
  }
  
  coef1 /= SQ(mesh->dx[jx][jy]);
  coef3 /= 2.*mesh->dx[jx][jy];
  coef4 /= 2.*mesh->dx[jx][jy];

  a = dcomplex(coef1 - coef4,-kwave*coef3);
  b = dcomplex(-2.0*coef1 - SQ(kwave)*coef2,kwave*coef5);
  c = dcomplex(coef1 + coef4,kwave*coef3);
}

/// Sets the coefficients for parallel tridiagonal matrix inversion
/*!
 * Uses the laplace_tridag_coefs routine above to fill a matrix [kz][ix] of coefficients
 */
void Laplacian::tridagMatrix(dcomplex **avec, dcomplex **bvec, dcomplex **cvec,
                             dcomplex **bk, int jy, int flags, 
                             const Field2D *a, const Field2D *ccoef, 
                             const Field2D *d) {
  int ix, kz;
  
  int ncx = mesh->ngx-1;

  int xbndry = 2;
  if(flags & INVERT_BNDRY_ONE)
    xbndry = 1;
  
  for(kz = 0; kz <= maxmode; kz++) {
    
    // Entire domain. Change to boundaries later

    for(ix=0;ix<=ncx;ix++) {

      tridagCoefs(ix, jy, kz, avec[kz][ix], bvec[kz][ix], cvec[kz][ix], ccoef, d);
      
      if(a != (Field2D*) NULL)
	bvec[kz][ix] += (*a)[ix][jy];
    }

    if(!mesh->periodicX) {
      // Boundary conditions
      
      if(mesh->firstX()) {
	// INNER BOUNDARY ON THIS PROCESSOR
	
	if(!(flags & (INVERT_IN_RHS | INVERT_IN_SET))) {
	  for(ix=0;ix<xbndry;ix++)
	    bk[kz][ix] = 0.;
	}
	
	if(kz == 0) {
	  // DC
	  
	  if(flags & INVERT_DC_IN_GRAD) {
	    // Zero gradient at inner boundary
	    for (ix=0;ix<xbndry;ix++){
	      avec[kz][ix] =  0.;
	      bvec[kz][ix] =  1.;
	      cvec[kz][ix] = -1.;
	    }
	  }else {
	    // Zero value at inner boundary or INVERT_IN_SET
	    for (ix=0;ix<xbndry;ix++){
	      avec[kz][ix] = 0.;
	      bvec[kz][ix] = 1.;
	      cvec[kz][ix] = 0.;
	    }
	  }
	  
	}else {
	  // AC
	
	  if(flags & INVERT_AC_IN_GRAD) {
	    // Zero gradient at inner boundary
	    for (ix=0;ix<xbndry;ix++){
	      avec[kz][ix]=dcomplex(0.,0.);
	      bvec[kz][ix]=dcomplex(1.,0.);
	      cvec[kz][ix]=dcomplex(-1.,0.);
	    }
	  }else if(flags & INVERT_AC_IN_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=kz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (ix=0;ix<xbndry;ix++) {
	      avec[kz][ix] = 0.0;
	      bvec[kz][ix] = 1.0;
	      cvec[kz][ix] = -exp(-1.0*sqrt(mesh->g33[ix][jy]/mesh->g11[ix][jy])*kwave*mesh->dx[ix][jy]);
	    }
	  }else {
	    // Zero value at inner boundary or INVERT_IN_SET
	    for (ix=0;ix<xbndry;ix++){
	      avec[kz][ix]=dcomplex(0.,0.);
	      bvec[kz][ix]=dcomplex(1.,0.);
	      cvec[kz][ix]=dcomplex(0.,0.);
	    }
	  }
	}
      }else if(mesh->lastX()) {
	// OUTER BOUNDARY
      
	if(!(flags & (INVERT_OUT_RHS | INVERT_OUT_SET))) {
	  for (ix=0;ix<xbndry;ix++)
	    bk[kz][ncx-ix] = 0.;
	}

	if(kz == 0) {
	  // DC
	
	  if(flags & INVERT_DC_OUT_GRAD) {
	    // Zero gradient at outer boundary
	    for (ix=0;ix<xbndry;ix++){
	      cvec[kz][ncx-ix]=dcomplex(0.,0.);
	      bvec[kz][ncx-ix]=dcomplex(1.,0.);
	      avec[kz][ncx-ix]=dcomplex(-1.,0.);
	    }
	  }else {
	    // Zero value at outer boundary or INVERT_OUT_SET
	    for (ix=0;ix<xbndry;ix++){
	      cvec[kz][ncx-ix]=dcomplex(0.,0.);
	      bvec[kz][ncx-ix]=dcomplex(1.,0.);
	      avec[kz][ncx-ix]=dcomplex(0.,0.);
	    }
	  }
	}else {
	  // AC
	
	  if(flags & INVERT_AC_OUT_GRAD) {
	    // Zero gradient at outer boundary
	    for (ix=0;ix<xbndry;ix++){
	      cvec[kz][ncx-ix]=dcomplex(0.,0.);
	      bvec[kz][ncx-ix]=dcomplex(1.,0.);
	      avec[kz][ncx-ix]=dcomplex(-1.,0.);
	    }
	  }else if(flags & INVERT_AC_OUT_LAP) {
	    // Use decaying zero-Laplacian solution in the boundary
	    BoutReal kwave=kz*2.0*PI/mesh->zlength; // wave number is 1/[rad]
	    for (ix=0;ix<xbndry;ix++) {
	      avec[kz][ncx-ix] = -exp(-1.0*sqrt(mesh->g33[ncx-ix][jy]/mesh->g11[ncx-ix][jy])*kwave*mesh->dx[ncx-ix][jy]);;
	      bvec[kz][ncx-ix] = 1.0;
	      cvec[kz][ncx-ix] = 0.0;
	    }
	  }else {
	    // Zero value at outer boundary or LAPLACE_OUT_SET
	    for (ix=0;ix<xbndry;ix++){
	      cvec[kz][ncx-ix]=dcomplex(0.,0.);
	      bvec[kz][ncx-ix]=dcomplex(1.,0.);
	      avec[kz][ncx-ix]=dcomplex(0.,0.);
	    }
	  }
	}
      }
    }
  }  
}

/**********************************************************************************
 *                              LEGACY INTERFACE
 *
 * These functions are depreciated, and will be removed in future
 **********************************************************************************/

/// Returns the coefficients for a tridiagonal matrix for laplace. Used by Delp2 too
void laplace_tridag_coefs(int jx, int jy, int jz, dcomplex &a, dcomplex &b, dcomplex &c, 
                          const Field2D *ccoef, const Field2D *d) {
  Laplacian::defaultInstance()->tridagCoefs(jx,jy, jz, a, b, c, ccoef, d);
}

int invert_laplace(const FieldPerp &b, FieldPerp &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  
  Laplacian *lap = Laplacian::defaultInstance();
  
  if(a != NULL) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);
  
  if(c != NULL) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);
  
  if(d != NULL) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);
  
  lap->setFlags(flags);
  
  x = lap->solve(b);
  
  return 0;
}

int invert_laplace(const Field3D &b, Field3D &x, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  
  Timer timer("invert"); ///< Start timer
  
  Laplacian *lap = Laplacian::defaultInstance();
  
  if(a != NULL) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);
  
  if(c != NULL) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);
  
  if(d != NULL) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);
  
  lap->setFlags(flags);
  
  x.allocate(); // Make sure x is allocated

  x = lap->solve(b, x);
  
}
const Field3D invert_laplace(const Field3D &b, int flags, const Field2D *a, const Field2D *c, const Field2D *d) {
  
  Timer timer("invert"); ///< Start timer 
  
  Laplacian *lap = Laplacian::defaultInstance();
  
  if(a != NULL) {
    lap->setCoefA(*a);
  }else
    lap->setCoefA(0.0);
  
  if(c != NULL) {
    lap->setCoefC(*c);
  }else
    lap->setCoefC(1.0);
  
  if(d != NULL) {
    lap->setCoefD(*d);
  }else
    lap->setCoefD(1.0);
  
  lap->setFlags(flags);
  
  Field3D x = lap->solve(b);
  
  return x;
}

