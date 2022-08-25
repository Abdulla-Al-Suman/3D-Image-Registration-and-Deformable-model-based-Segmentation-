// Hessian matrix calculation
#include "mex.h"


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

mwSize  ndims, ndims1;
const mwSize *dims, *dims1;
double *in1,*in2,*out;
int ite, ne;
// Check number of arguments

    in1 = mxGetPr(prhs[0]); // get a pointers to the first input argument
    in2 = mxGetPr(prhs[1]); // ...
    
    ne=int(mxGetNumberOfElements(prhs[0]));

    
        // create the output
    ndims = mxGetNumberOfDimensions(prhs[2]);
    dims = mxGetDimensions(prhs[2]);
    plhs[0] = mxCreateNumericArray(ndims,dims,mxDOUBLE_CLASS,mxREAL);
    // get a pointer to the beginning
    out = mxGetPr(plhs[0]);
    
    *out=0.0;
    

    for (ite=0 ; ite < ne ; ++ite)
        *out=*out+in1[ite]*in2[ite];

} ;