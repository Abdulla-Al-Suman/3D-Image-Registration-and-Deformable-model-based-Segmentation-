
  
#include <mex.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <iostream>

inline 
static
int access(int M, int N, int O, int x, int y, int z) {
  if (x<0) x=0; else if (x>=N) x=N-1;
  if (y<0) y=0; else if (y>=M) y=M-1;
  if (z<0) z=0; else if (z>=O) z=O-1;
  return y + M*(x + N*z);
}

inline 
static
int access_unchecked(int M, int N, int O, int x, int y, int z) {
  return y + M*(x + N*z);
}

inline
static
void indices_linear(
    int &f000_i,
    int &f100_i,
    int &f010_i,
    int &f110_i,
    int &f001_i,
    int &f101_i,
    int &f011_i,
    int &f111_i,
    const int x, const int y, const int z,
    const mwSize &M, const mwSize &N, const mwSize &O) {
  if (x<=1 || y<=1 || z<=1 || x>=N-2 || y>=M-2 || z>=O-2) {
    f000_i = access(M,N,O, x,   y  , z);
    f100_i = access(M,N,O, x+1, y  , z);

    f010_i = access(M,N,O, x,   y+1, z);
    f110_i = access(M,N,O, x+1, y+1, z);

    f001_i = access(M,N,O, x,   y  , z+1);
    f101_i = access(M,N,O, x+1, y  , z+1);

    f011_i = access(M,N,O, x,   y+1, z+1);
    f111_i = access(M,N,O, x+1, y+1, z+1);
  } else {
    f000_i = access_unchecked(M,N,O, x,   y  , z);
    f100_i = access_unchecked(M,N,O, x+1, y  , z);

    f010_i = access_unchecked(M,N,O, x,   y+1, z);
    f110_i = access_unchecked(M,N,O, x+1, y+1, z);

    f001_i = access_unchecked(M,N,O, x,   y  , z+1);
    f101_i = access_unchecked(M,N,O, x+1, y  , z+1);

    f011_i = access_unchecked(M,N,O, x,   y+1, z+1);
    f111_i = access_unchecked(M,N,O, x+1, y+1, z+1);
  }
}



static
void interpolate_linear(double *pO, const double *pF, 
  const double *pX, const double *pY, const double *pZ,
  const mwSize ND, const mwSize M, const mwSize N, const mwSize O,
  const double s_x, const double o_x,
  const double s_y, const double o_y,
  const double s_z, const double o_z) {
  const mwSize LO = M*N*O;
  for (mwSize i=0; i<ND; ++i) {
    const double &x_ = pX[i];
    const double &y_ = pY[i];
    const double &z_ = pZ[i];
    
    const double x = s_x*x_+o_x;
    const double y = s_y*y_+o_y;
    const double z = s_z*z_+o_z;

    const double x_floor = floor(x);
    const double y_floor = floor(y);
    const double z_floor = floor(z);

    const double dx = x-x_floor;
    const double dy = y-y_floor;
    const double dz = z-z_floor;

    const double wx0 = 1.0-dx;
    const double wx1 = dx;

    const double wy0 = 1.0-dy;
    const double wy1 = dy;

    const double wz0 = 1.0-dz;
    const double wz1 = dz;

    int f000_i, f100_i, f010_i, f110_i;
    int f001_i, f101_i, f011_i, f111_i;

    // TODO: Use openmp

    indices_linear(
        f000_i, f100_i, f010_i, f110_i, 
        f001_i, f101_i, f011_i, f111_i, 
        int(x_floor-1), int(y_floor-1), int(z_floor-1), M, N, O);



      pO[i] =
        wz0*(
            wy0*(wx0 * pF[f000_i] + wx1 * pF[f100_i]) +
            wy1*(wx0 * pF[f010_i] + wx1 * pF[f110_i])
            )+
        wz1*(
            wy0*(wx0 * pF[f001_i] + wx1 * pF[f101_i]) +
            wy1*(wx0 * pF[f011_i] + wx1 * pF[f111_i])
            );

  }
}



void mexFunction(int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]) {

  if (nlhs>1)
    mexErrMsgTxt("Wrong number of output arguments for Z = interp3_mex(Fx, Fy, Fz, F, X, Y, Z, method)");

  const mxArray *Fx = NULL;
  const mxArray *Fy = NULL;
  const mxArray *Fz = NULL;
  const mxArray *F = NULL;
  const mxArray *X = NULL;
  const mxArray *Y = NULL;
  const mxArray *Z = NULL;
  const mxArray *method = NULL;

  if (nrhs==7) {
//     // ba_interp(F, X, Y, Z);
//     F = prhs[0];
//     X = prhs[1];
//     Y = prhs[2];
//     Z = prhs[3];
//   } else if (nrhs==5) {
//     // ba_interp(F, X, Y, Z, 'method');
//     F = prhs[0];
//     X = prhs[1];
//     Y = prhs[2];
//     Z = prhs[3];
//     method = prhs[4];
//   } else if (nrhs==7) {
    // ba_interp(Fx, Fy, Fz, F, X, Y, Z);
    Fx= prhs[0];
    Fy= prhs[1];
    Fz= prhs[2];
    F = prhs[3];
    X = prhs[4];
    Y = prhs[5];
    Z = prhs[6];

  } else if (nrhs==8) {
    // ba_interp(Fx, Fy, Fz, F, X, Y, Z, 'method');
    Fx= prhs[0];
    Fy= prhs[1];
    Fz= prhs[2];
    F = prhs[3];
    X = prhs[4];
    Y = prhs[5];
    Z = prhs[6];
    method = prhs[7];
  } else {
    mexErrMsgTxt("Wrong number of input arguments for Z = ba_interp3(Fx, Fy, Fz, F, X, Y, Z, method)");
  }
  if ((Fx && !mxIsDouble(Fx)) ||(Fy && !mxIsDouble(Fy)) ||(Fz && !mxIsDouble(Fz)) ||
      (F && !mxIsDouble(F)) ||
      (X && !mxIsDouble(X)) || (Y && !mxIsDouble(Y)) || (Z && !mxIsDouble(Z)))
    mexErrMsgTxt("ba_interp3 takes only double arguments for Fx,Fy,Fz,F,X,Y,Z");

  const mwSize *F_dims = mxGetDimensions(F);
  const mwSize *X_dims = mxGetDimensions(X);
  const mwSize *Y_dims = mxGetDimensions(Y);
  const mwSize *Z_dims = mxGetDimensions(Z);

  const mwSize M=F_dims[0];
  const mwSize N=F_dims[1];
  const mwSize O=F_dims[2];

  if (Fx && mxGetNumberOfElements(Fx)<2) mexErrMsgTxt("Fx needs at least two elements.");
  if (Fy && mxGetNumberOfElements(Fy)<2) mexErrMsgTxt("Fy needs at least two elements.");
  if (Fz && mxGetNumberOfElements(Fz)<2) mexErrMsgTxt("Fz needs at least two elements.");

  if ((mxGetNumberOfDimensions(X) != mxGetNumberOfDimensions(Y)) ||
      (mxGetNumberOfDimensions(X) != mxGetNumberOfDimensions(Z)) ||
      (mxGetNumberOfElements(X) != mxGetNumberOfElements(Y)) ||
      (mxGetNumberOfElements(X) != mxGetNumberOfElements(Z)))
    mexErrMsgTxt("X, Y, Z should have the same size");

  mwSize P=1;

  mwSize outDims[50];
  if (mxGetNumberOfDimensions(X) + mxGetNumberOfDimensions(F) - 3 > 50)
    mexErrMsgTxt("Can't have that many dimensions in interpolated data.");

  for (mwSize i=0; i<mxGetNumberOfDimensions(X); ++i) outDims[i] = X_dims[i];

  for (mwSize i=3; i<mxGetNumberOfDimensions(F); ++i) {
    outDims[mxGetNumberOfDimensions(X)+i-3] = F_dims[i];
    P *= F_dims[i];
  }


  plhs[0] = mxCreateNumericArray(mxGetNumberOfDimensions(X) + mxGetNumberOfDimensions(F) - 3, outDims, mxDOUBLE_CLASS, mxREAL);

  const mwSize ND = mxGetNumberOfElements(X);

  const double *pF = mxGetPr(F);
  const double *pX = mxGetPr(X);
  const double *pY = mxGetPr(Y);
  const double *pZ = mxGetPr(Z);
  double       *pO = mxGetPr(plhs[0]);  
  
    if (Fx) { 
    
    const double x_low = mxGetPr(Fx)[0]; const double x_high = mxGetPr(Fx)[mxGetNumberOfElements(Fx)-1];
    const double y_low = mxGetPr(Fy)[0]; const double y_high = mxGetPr(Fy)[mxGetNumberOfElements(Fy)-1];
    const double z_low = mxGetPr(Fz)[0]; const double z_high = mxGetPr(Fz)[mxGetNumberOfElements(Fz)-1];
    
    const double s_x = double(1-N)/(x_low - x_high);
    const double s_y = double(1-M)/(y_low - y_high);
    const double s_z = double(1-O)/(z_low - z_high);

    const double o_x = double(1)-x_low*s_x;
    const double o_y = double(1)-y_low*s_y;
    const double o_z = double(1)-z_low*s_z;
    

     interpolate_linear(pO, pF, pX, pY, pZ, ND, M, N, O, s_x, o_x, s_y, o_y, s_z, o_z);
  
    }
    else
        interpolate_linear(pO, pF, pX, pY, pZ, ND, M, N, O, double(1), double(0), double(1), double(0), double(1), double(0));
    
      

}
