#include "mex.h"
#include "math.h"
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
    double *Dout;
    double a=0;
    int i,n;
  
    n = mxGetScalar(prhs[0]);
    for(i=1;i<=n;i++)
        a=a+i; 
    
    plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
    Dout=mxGetData(plhs[0]);
    
    Dout[0]=a;
}

