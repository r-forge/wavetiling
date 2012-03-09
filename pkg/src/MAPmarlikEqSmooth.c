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
