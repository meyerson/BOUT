/**************************************************************
 * Smoothing operators
 *
 *
 * 2010-05-17 Ben Dudson <bd512@york.ac.uk>
 *    * Added nonlinear filter
 * 
 **************************************************************
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
 **************************************************************/

#include <math.h>

#include <globals.hxx>
#include <smoothing.hxx>
#include <bout_types.hxx>
#include <fft.hxx>
#include <dcomplex.hxx>

// Smooth using simple 1-2-1 filter
const Field3D smooth_x(const Field3D &f, bool BoutRealspace) {
  Field3D fs, result;

  if(BoutRealspace) {
    fs = f.shiftZ(true); // Shift into BoutReal space
  }else
    fs = f; 

  result.allocate();
  
  // Copy boundary region
  for(int jy=0;jy<mesh->ngy;jy++)
    for(int jz=0;jz<mesh->ngz;jz++) {
      result[0][jy][jz] = fs[0][jy][jz];
      result[mesh->ngx-1][jy][jz] = fs[mesh->ngx-1][jy][jz];
    }

  // Smooth using simple 1-2-1 filter

  for(int jx=1;jx<mesh->ngx-1;jx++)
    for(int jy=0;jy<mesh->ngy;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	result[jx][jy][jz] = 0.5*fs[jx][jy][jz] + 0.25*( fs[jx-1][jy][jz] + fs[jx+1][jy][jz] );
      }

  if(BoutRealspace)
    result = result.shiftZ(false); // Shift back

  // Need to communicate boundaries
  mesh->communicate(result);

  return result;
}


const Field3D smooth_y(const Field3D &f) {
  Field3D result;

  result.allocate();
  
  // Copy boundary region
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz;jz++) {
      result[jx][0][jz] = f[jx][0][jz];
      result[jx][mesh->ngy-1][jz] = f[jx][mesh->ngy-1][jz];
    }
  
  // Smooth using simple 1-2-1 filter

  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=1;jy<mesh->ngy-1;jy++)
      for(int jz=0;jz<mesh->ngz;jz++) {
	result[jx][jy][jz] = 0.5*f[jx][jy][jz] + 0.25*( f[jx][jy-1][jz] + f[jx][jy+1][jz] );
	//result[jx][jy][jz] = 0.5*f[jx][jy][jz] + 0.25*( f[jx][jy-1][jz] + f[jx][jy+1][jz] );
      }

  // Need to communicate boundaries
  mesh->communicate(result);
  
  return result;
}

const Field2D averageY(const Field2D &f) {
  return mesh->averageY(f);
}

const Field3D smoothXY(const Field3D &f) {
  Field3D result;
  result.allocate();

  for(int x=2;x<mesh->ngx-2;x++)
    for(int y=2;y<mesh->ngy-2;y++)
      for(int z=0;z<mesh->ngz;z++) {
        result[x][y][z] = 0.5*f[x][y][z] + 0.125*( 0.5*f[x+1][y][z] + 0.125*(f[x+2][y][z] + f[x][y][z] + f[x+1][y-1][z] + f[x+1][y+1][z]) +
                                                   0.5*f[x-1][y][z] + 0.125*(f[x][y][z] + f[x-2][y][z] + f[x-1][y-1][z] + f[x-1][y+1][z]) +
                                                   0.5*f[x][y-1][z] + 0.125*(f[x+1][y-1][z] + f[x-1][y-1][z] + f[x][y-2][z] + f[x][y][z]) +
                                                   0.5*f[x][y+1][z] + 0.125*(f[x+1][y+1][z] + f[x-1][y+1][z] + f[x][y][z] + f[x][y+2][z]));
      }
  
  return result;
}

/// Nonlinear filtering to remove grid-scale noise
/*!
  From a paper:

  W.Shyy et. al. JCP 102 (1) September 1992 page 49

  "On the Suppression of Numerical Oscillations Using a Non-Linear Filter"
  
 */
void nl_filter(rvec &f, BoutReal w) {
  for(size_t i=1; i<f.size()-1; i++) {
    
    BoutReal dp = f[i+1] - f[i];
    BoutReal dm = f[i-1] - f[i];
    if(dp*dm > 0.) {
      // Local extrema - adjust
      BoutReal ep, em, e; // Amount to adjust by
      if(fabs(dp) > fabs(dm)) {
	ep = w*0.5*dp;
	em = w*dm;
	e = (fabs(ep) < fabs(em)) ? ep : em; // Pick smallest absolute
	// Adjust
	f[i+1] -= e;
	f[i] += e;
      }else {
	ep = w*0.5*dm;
	em = w*dp;
	e = (fabs(ep) < fabs(em)) ? ep : em; // Pick smallest absolute
	// Adjust
	f[i-1] -= e;
	f[i] += e;
      }
    }
  }
}

const Field3D nl_filter_x(const Field3D &f, BoutReal w) {
#ifdef CHECK
  msg_stack.push("nl_filter_x( Field3D )");
#endif
  
  Field3D fs;
  fs = f.shiftZ(true); // Shift into BoutReal space
  Field3D result;
  rvec v;
  
  for(int jy=0;jy<mesh->ngy;jy++)
    for(int jz=0;jz<mesh->ngz-1;jz++) {
      fs.getXArray(jy, jz, v);
      nl_filter(v, w);
      result.setXArray(jy, jz, v);
    }
  
  result = result.shiftZ(false); // Shift back

#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D nl_filter_y(const Field3D &fs, BoutReal w) {
#ifdef CHECK
  msg_stack.push("nl_filter_x( Field3D )");
#endif
  
  Field3D result;
  rvec v;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz-1;jz++) {
      fs.getYArray(jx, jz, v);
      nl_filter(v, w);
      result.setYArray(jx, jz, v);
    }
  
#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D nl_filter_z(const Field3D &fs, BoutReal w) {
#ifdef CHECK
  msg_stack.push("nl_filter_z( Field3D )");
#endif
  
  Field3D result;
  rvec v;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jy=0;jy<mesh->ngy;jy++) {
      fs.getZArray(jx, jy, v);
      nl_filter(v, w);
      result.setZArray(jx, jy, v);
    }

#ifdef CHECK
  msg_stack.pop();
#endif
  return result;
}

const Field3D nl_filter(const Field3D &f, BoutReal w)
{
  Field3D result;
  /// Perform filtering in Z, Y then X
  result = nl_filter_x(nl_filter_y(nl_filter_z(f, w), w), w);
  /// Communicate boundaries
  mesh->communicate(result);
  return result;
}

const Field3D lowPass_Y(const Field3D &var, int ymax)
{
  Field3D result;
  static dcomplex *f = NULL;
  int jx, jy, jz;
  //static double *g = NULL;
  //static BoutReal *g;  //an 1d array given x and z coords, need to allocate the correct size
  Field2D g; // nx x ny data container
  BoutReal **d; //an array of doubl, again nx x ny
  
  
#ifdef CHECK
  msg_stack.push("lowPass_Y(Field3D, %d)", ymax);
#endif
  int ncy = mesh->ngy;
 
  if(!var.isAllocated())
    return var;
  
  if(f == NULL)
    f = new dcomplex[ncy/2 + 1];
 
  if((ymax >= ncy/2) || (ymax < 0)) {
    // Removing nothing
    return var;
  }
  
  result.allocate();
  
  g = 0.0;
  d = g.getData();
  BoutReal rptr =0.0;
  int stat = g.getData(2,8,11,&rptr);
  //output.write("(%e),%d \n ",rptr,stat);

  //SurfaceIter * iterateSurfaces ( ) ;
  // Create an iterator over surfaces
  SurfaceIter* surf = mesh->iterateSurfaces();
  //int ysize = surf->ysize();
  //output.write("%d \n ",ysize);
  // for(surf->first(); !surf->isDone(); surf->next()) {
  //   int ysize = surf->ysize(); // Size of this surface
  //   //output.write("%d \n ",ysize);
  // }

  // output.write("sizeof(d) = %d,%d,%d,%d,\n %d,%d,%d \n ",sizeof(d),sizeof(d[0]),
  // 	       sizeof(*d),sizeof(*d[0]),
  // 	       mesh->ngz,mesh->ngy,mesh->ngx);
  // output.write("sizeof data[jx] = %d\n ",sizeof(var[0]));



  //return var;
  //g = var[0][0]; //crappy for now
  // #pragma omp parallel for
  // for(int j=0;j<mesh->ngx*mesh->ngy;j++) {
  //   for(int jz=0;jz<(mesh->ngz-1);jz++)
  //     d[0][j] = var[0][j][jz];  
  // }
  //return var;
  
  for(int jx=0;jx<mesh->ngx;jx++)
    for(int jz=0;jz<mesh->ngz-1;jz++){
      for(jy=0;jy<mesh->ngy;jy++){
	//output.write("(%d,%d,%d) \n ",jz,mesh->ngy,jy);
	BoutReal ycoord = mesh->GlobalY(jy);
	BoutReal xcoord = mesh->GlobalX(jx);
	
	//output.write("(%e,%e) \n ",xcoord,ycoord);
	d[0][jy] = var[jx][jy][jz];  
      }
      //x and z fixed at this point
      rfft(d[0], ncy, f);
      //for(jy=ymax+1;jy<=ncy/2;jy++)
      for(jy=ymax+1;jy<=ncy/2;jy++)
	f[jy] = 0.0;
	//f[jy] = f[jy];///1000.0;
      //f[1] = 1.0;
      int size = sizeof(f) / sizeof(f[0] );
      //output.write("%d, %d, %d \n ",size,sizeof(f),sizeof(d[0]));
      //f[0] = 1.0;

      irfft(f, ncy, d[0]);
      
      for(jy=0;jy<mesh->ngy;jy++)
	result[jx][jy][jz] = d[0][jy];
      //result[jx][ncy][jz] = result[jx][jy][0];
    }
  

  
  mesh->communicate(result);
  
  return result;

  //copy the data in an array g[x][z][y]

  //fft,filter, ifft g
  
}
