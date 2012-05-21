#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

void R_init_waveTiling( DllInfo *info );
void R_unload_waveTiling( DllInfo *info );

void MAPMARG(double *Din,int *K, double *vareps,double *Bout,double *varBout,double *Phiout, double *X, double *XtX, int *q,int *N,int *ends)
{
int i,j,m,n,Jcount=0;
double phi;
for (i=0;i<*K;i++) 
{
	if (i>=ends[Jcount]) Jcount++;
	for (m=0;m<*q;m++) 
	{	
		phi=0;
		for (j=0;j<*N;j++) 
		{	
			phi+=X[j+m*(*N)]*Din[j+i*(*N)];		
		}
		Bout[m*(*K)+i]=phi;
		phi=phi*phi/XtX[m]/XtX[m]/(vareps[Jcount])-1/XtX[m];
		if (phi<=0) 
			{
				Bout[m*(*K)+i]=0;	
				varBout[m*(*K)+i]=0;
				Phiout[m*(*K)+i]=0;
			}
		else 	
			{
				Bout[m*(*K)+i]=Bout[m*(*K)+i]/(XtX[m]+1/phi);				    
				varBout[m*(*K)+i]=(vareps[Jcount])/(XtX[m]+1/phi);
				Phiout[m*(*K)+i]=phi;
			};
	
	}
}
} 

void MAPMARGEQSMOOTH(double *Din,int *K, double *vareps,double *Bout,double *varBout,double *Phiout, double *X, double *XtX, int *q,int *N,int *ends)
{
int i,j,m,n,Jcount=0;
double phi,XDm,XD;
for (i=0;i<*K;i++) 
{
	if (i>=ends[Jcount]) Jcount++;
	XD=0;
	for (m=0;m<*q;m++) 
	{
		XDm=0;
		for (j=0;j<*N;j++) 
		{	
			XDm+=X[j+m*(*N)]*Din[j+i*(*N)];
		}
		Bout[m*(*K)+i]=XDm;
		XD+=XDm;
	}
	phi=XD*XD/(*q)/(vareps[Jcount])-1;
	for (m=0;m<*q;m++)
	{
		if (phi<=0) 
			{
				Bout[m*(*K)+i]=0;	
				varBout[m*(*K)+i]=0;
				Phiout[m*(*K)+i]=0;
			}
		else 	
			{
				Bout[m*(*K)+i]=Bout[m*(*K)+i]/(XtX[m]+1/phi);
				varBout[m*(*K)+i]=(vareps[Jcount])/(XtX[m]+1/phi);
				Phiout[m*(*K)+i]=phi;
			};
	
	}
}
}

void MAPMARGIMP(double *Din,int *K, double *vareps,double *Bout,double *varBout,double *Phiout, double *X, double *XtX, int *q,int *N,int *ends)
{
int i,j,m,n,Jcount=0;
double phi;
for (i=0;i<*K;i++) 
{
	if (i>=ends[Jcount]) Jcount++;
	for (m=0;m<*q;m++) 
	{	
		phi=0;
		for (j=0;j<*N;j++) 
		{	
			phi+=X[j+m*(*N)]*Din[j+i*(*N)];		
		}
		Bout[m*(*K)+i]=phi;
		phi=phi*phi/XtX[m]/XtX[m]/(vareps[Jcount])/3-1/XtX[m];
		if (phi<=0) 
			{
				Bout[m*(*K)+i]=0;	
				varBout[m*(*K)+i]=0;
				Phiout[m*(*K)+i]=0;
			}
		else 	
			{
				Bout[m*(*K)+i]=Bout[m*(*K)+i]/(XtX[m]+1/phi);				    
				varBout[m*(*K)+i]=(vareps[Jcount])/(XtX[m]+1/phi);
				Phiout[m*(*K)+i]=phi;
			};
	
	}
}
}

/* Registration information for DLL */
R_CMethodDef cMethods[] = {
    { "MAPMARG", ( DL_FUNC ) &MAPMARG, 11, /*{ REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP }*/ },
    { NULL, NULL, 0 },
    { "MAPMARGEQSMOOTH", ( DL_FUNC ) &MAPMARGEQSMOOTH, 11, /*{ REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP }*/  },
    { NULL, NULL, 0 },
    { "MAPMARGIMP", ( DL_FUNC ) &MAPMARGIMP, 11, /*{ REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP }*/  },
    { NULL, NULL, 0 }
};


void R_init_waveTiling( DllInfo *info ) {
    R_registerRoutines( info, cMethods, NULL, NULL, NULL );
}

void R_unload_waveTiling( DllInfo *info ) {
  /* Release resources. */
}



