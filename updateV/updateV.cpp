// updateV.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"


#include "mex.h" 

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>


void doV(double *XTU, double *V, double lamda1, double tol, int n, int k, int maxinner, double *Vout)
{
    int idx;
	for (idx=0; idx< k*n ; idx++ )
   {		
		// get s*
		double temp = XTU[idx]/lamda1;
		double s = (temp <0 ? 0:temp)-V[idx];

		if (0==V[idx] && s<0)
		{
			Vout[idx]=0;
		}
		else{
			Vout[idx] = V[idx]+s;
		}

   }

}

void usage()
{
	printf("Error calling vThreads.\n");
	printf("Usage:  V_res = vThreads(XTU,V,lambda1,tol, k^2)\n");
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	double *XTU, *V, *values, *Vout, *U;
	double tol, lamda1; 
	int maxinner, k, n, i;
	
	// Input Arguments

	XTU = mxGetPr(prhs[0]);
	n = mxGetM(prhs[0]); // rows
	k = mxGetN(prhs[0]); // columns

	//printf("xtu: is a %d by %d matrix\n",  n, k);
	
	V = mxGetPr(prhs[1]);
	if ( (mxGetM(prhs[1]) != n) || (mxGetN(prhs[1])!=k) ) {
		//usage();
		printf("Error: V should be a %d by %d matrix\n",  n, k);
	}


	values = mxGetPr(prhs[2]);
	lamda1 = values[0];

	values = mxGetPr(prhs[3]);
	tol = values[0];

	values = mxGetPr(prhs[4]);
	maxinner = values[0];


	/// Output arguments
	plhs[0] = mxCreateDoubleMatrix(n,k,mxREAL);	
	Vout = mxGetPr(plhs[0]);
	
	doV(XTU, V, lamda1, tol, n, k, maxinner, Vout);

}


