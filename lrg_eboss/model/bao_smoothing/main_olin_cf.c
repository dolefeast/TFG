#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define Pi (4.*atan(1.))
#include "gsl/gsl_multimin.h"
double lowlim,highlim;

//#define Pi (3.14159265358979323846264338)
double P_interpol(double k0, double *k, double *P, long int N);

double Xi_interpol(double k0, double *k, double *P, long int N);

double chi2 (const gsl_vector *A, void *params);

double poly(double r, double cm3, double cm2, double cm1, double c0, double c1);

double *XI,*R;
int NN;

int main(int argc, char *argv[])
{
double kf;

long int i,N,j;
double *k,*Pk;
double *kgrid,*Pgrid;
long int Ngrid;
//double lowlim,highlim;
int power=atoi(argv[1]);//15, 16, 17
Ngrid=pow(2,power);//32768;
kgrid=malloc(sizeof(double)*Ngrid*2);
Pgrid=malloc(sizeof(double)*Ngrid*2);
FILE *f;
char name[200],name2[200];

N=576;
//N=617;
//N=543;
//N=13775;// Pk_fiducial_eboss_z1.51_hr.txt
//N=13807;
k=malloc(sizeof(double)*N);
Pk=malloc(sizeof(double)*N);

//sprintf(name2,"Pk_linear_z0");lowlim=86.;highlim=150.;
//sprintf(name2,"mc_pk");;lowlim=86.;highlim=150.;
sprintf(name2,"eboss_comb_z070_matterpower");lowlim=86.;highlim=150.;
//sprintf(name2,"eboss_comb_z072_higherprecision_matterpower");lowlim=86.;highlim=150.;
//sprintf(name2,"Nseries_cmass_z055_higherprecision_matterpower");lowlim=86.;highlim=150.;
//sprintf(name2,"Nseries_cmass_z055_matterpower");lowlim=86.;highlim=150.;
//sprintf(name2,"OmC_cosmo_z055_matterpower");lowlim=86.;highlim=150.;
//sprintf(name2,"OmX_cosmo_z070_matterpower");lowlim=76.;highlim=140.;
//sprintf(name2,"Pk_fiducial_eboss_z1.51");
sprintf(name,"%s.dat",name2);

f=fopen(name,"r");
if(f==NULL){printf("No es pot obrir %s\n",name);return 0;}

for(i=0;i<N;i++)
{
fscanf(f,"%lf %lf\n",&k[i],&Pk[i]);
}
fclose(f);

kf=k[0]*1;
kgrid[0]=0;
Pgrid[0]=0;


sprintf(name,"%s_interpol.txt",name2);
//f=fopen(name,"w");
if(f==NULL){printf("No es pot obrir %s\n",name);return 0;}

for(i=1;i<Ngrid*2;i++)
{
if(i<=Ngrid)
{
kgrid[i]=i*kf;
Pgrid[i]=P_interpol(kgrid[i],k,Pk,N);
//fprintf(f,"%lf %.16lf\n",kgrid[i],Pgrid[i]);
}
else
{
kgrid[i]=-(2*Ngrid-i)*kf;
Pgrid[i]=P_interpol(fabs(kgrid[i]),k,Pk,N);

}

}
//fclose(f);

//grid for CF
double *r,*xi;
r=malloc(sizeof(double)*Ngrid);
xi=malloc(sizeof(double)*Ngrid);

for(i=0;i<Ngrid;i++)
{
xi[i]=0;
}

//for(j=0;j<Ngrid*2;j++)//sum over kmodes
for(j=0;j<Ngrid;j++)//sum over positive kmodes only
{
for(i=0;i<Ngrid;i++)//sum over R
{
r[i]=2.*Pi*i*0.2;
//r[i]=2.*3.*i*0.2*5;
xi[i]=xi[i]+kf*Pgrid[j]*kgrid[j]*sin(r[i]*kgrid[j]);
}
}
R=malloc(sizeof(double)*Ngrid);
XI=malloc(sizeof(double)*Ngrid);

//printf("First loop done\n");
sprintf(name,"Xi_%s_%d.txt",name2,power);
f=fopen(name,"w");
if(f==NULL){printf("No es pot obrir %s\n",name);return 0;}

for(i=0;i<Ngrid;i++)
{
if(i==0)
{
xi[i]=0;

}
else
{
xi[i]=xi[i]/r[i];

}
fprintf(f,"%lf %lf\n",r[i],xi[i]);
XI[i]=xi[i];
R[i]=r[i];
}
fclose(f);
//Here we cut xi

NN=Ngrid;


const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *lam, *x;
  gsl_multimin_function minex_func;
  size_t iter = 0;
  int status;
  double size;

double cm3,cm2,cm1,c0,c1;
cm3=1;
cm2=1;
cm1=1;
c0=1;
c1=1;
  x = gsl_vector_alloc (5);
  gsl_vector_set (x, 0, cm3);
  gsl_vector_set (x, 1, cm2);
  gsl_vector_set (x, 2, cm1);
  gsl_vector_set (x, 3, c0);
  gsl_vector_set (x, 4, c1);
 
  lam = gsl_vector_alloc (5);
  gsl_vector_set (lam, 0, 10);
  gsl_vector_set (lam, 1, 10);
  gsl_vector_set (lam, 2, 10);
  gsl_vector_set (lam, 3, 10);
  gsl_vector_set (lam, 4, 10);

  minex_func.n = 5;
  minex_func.f = chi2;
    s = gsl_multimin_fminimizer_alloc (T, 5);

    gsl_multimin_fminimizer_set (s, &minex_func, x, lam);

do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);

      if (status)
        break;

      size = gsl_multimin_fminimizer_size (s);
      status = gsl_multimin_test_size (size, 1e-4);

      if (status == GSL_SUCCESS)
        {

        }
    }
  while (status == GSL_CONTINUE && iter < 9000);



cm3=gsl_vector_get (s->x, 0);
cm2=gsl_vector_get (s->x, 1);
cm1=gsl_vector_get (s->x, 2);
c0=gsl_vector_get (s->x, 3);
c1=gsl_vector_get (s->x, 4);


printf("%lf %lf %lf %lf %lf %d\n",cm3,cm2,cm1,c0,c1,iter);
  gsl_vector_free(x);
  gsl_vector_free(lam);
  gsl_multimin_fminimizer_free (s);

//Refill those xi[i] between r=86 and 150 with the bestfit polynomial extracted from 50-86 and 150-190.

sprintf(name,"Xi_sm_%s_%d.txt",name2,power);
f=fopen(name,"w");
if(f==NULL){printf("No es pot obrir %s\n",name);return 0;}

for(i=0;i<Ngrid;i++)
{

//if(r[i]>86. && r[i]<150.)
if(r[i]>lowlim && r[i]<highlim)
{
xi[i]=poly(r[i],cm3,cm2,cm1,c0,c1);
}
fprintf(f,"%lf %lf\n",r[i],xi[i]);
}
fclose(f);

double *Pgrid_fin,*kgrid_fin;
Pgrid_fin=malloc(sizeof(double)*Ngrid*2);
kgrid_fin=malloc(sizeof(double)*Ngrid*2);

for(i=0;i<Ngrid*2;i++)
{
Pgrid_fin[i]=0;
if(i<=Ngrid){kgrid_fin[i]=i*kf;}
if(i>Ngrid){kgrid_fin[i]=-(Ngrid*2-i)*kf;}

}



//for(i=0;i<Ngrid*2;i++)//over k
for(i=0;i<Ngrid;i++)//over positive k
{
for(j=0;j<Ngrid;j++)//over R
{

Pgrid_fin[i]=Pgrid_fin[i]+2.*Pi*xi[j]*r[j]*sin(r[j]*kgrid_fin[i]);

}
}

FILE *g,*h;

sprintf(name,"Olinkirkby_%s_%d.txt",name2,power);
g=fopen(name,"w");
//f=fopen("Pk_fin.txt","w");
sprintf(name,"Pk_smkirkby_%s_%d.txt",name2,power);
h=fopen(name,"w");
for(i=0;i<Ngrid;i++)
{
if(i==0)
{
Pgrid_fin[i]=0;
}
else
{
Pgrid_fin[i]=(Pgrid_fin[i]/(kgrid_fin[i])+Pgrid_fin[2*Ngrid-i]/(kgrid_fin[2*Ngrid-i])    )/(2.*Pi)*4./5.;

}
//fprintf(f,"%lf %.16lf\n",kgrid_fin[i],Pgrid_fin[i]);
if(kgrid_fin[i]>=0.003 && kgrid_fin[i]<=0.75){
fprintf(g,"%e %e\n",kgrid_fin[i],Pgrid[i]/Pgrid_fin[i]);
fprintf(h,"%e %e\n",kgrid_fin[i],Pgrid_fin[i]);
}
}
//fclose(f);
fclose(g);
fclose(h);

return 0;
}

double chi2 (const gsl_vector *A, void *params)
{
int points=0;
double cm3=gsl_vector_get(A, 0);
double cm2=gsl_vector_get(A, 1);
double cm1=gsl_vector_get(A, 2);
double c0=gsl_vector_get(A, 3);
double c1=gsl_vector_get(A, 4);

int i=0;

double chi2=0;
for(i=0;i<NN;i++)
{
//if(R[i]>50. && R[i]<86)
if(R[i]>50. && R[i]<lowlim)
{
chi2=chi2+1000*pow(poly(R[i],cm3,cm2,cm1,c0,c1)-XI[i],2);
points++;
}

//if(R[i]>150. && R[i]<190)
if(R[i]>highlim && R[i]<190)
{
chi2=chi2+1000*pow(poly(R[i],cm3,cm2,cm1,c0,c1)-XI[i],2);
points++;
}

}

printf("%lf (%d) %lf %lf %lf %lf %lf\n",chi2,points,cm3,cm2,cm1,c0,c1);

return  chi2;
}
double poly(double r, double cm3, double cm2, double cm1, double c0, double c1)
{
double xi;

xi=cm3*pow(r,-3)+cm2*pow(r,-2)+cm1*pow(r,-1)+c0+c1*r;
return xi;
}

double P_interpol(double k0, double *k, double *P, long int N)
{
double P0,m,n;
long int i;
if(k0<=k[0] || k0>=k[N-1])
{

if(k0>=k[N-1]){P0=0;printf("Warning\n");}

if(k0<=k[0])
{
P0=pow(k0,1)*P[0]/k[0];
}
}
else
{
i=-1;
do
{
i=i+1;
}while(k0>k[i]);
if(i==0)
{
m=( (P[i]) - (P[i+1]) )/( (k[i]) - (k[i+1]) );
n=(P[i])-m*(k[i]);
}
else
{
m=( (P[i]) - (P[i-1]) )/( ( k[i]) - (k[i-1]) );
n=(P[i])-m*(k[i]);
}
P0=m*(k0)+n;
}

return P0;
}

double Plog_interpol(double k0, double *k, double *P, long int N)
{
double P0,m,n;
long int i;
if(k0<=k[0] || k0>=k[N-1])
{

if(k0>=k[N-1]){P0=pow(k0,-2.7)*P[N-1]/k[N-1];}

if(k0<=k[0])
{
P0=pow(k0,1)*P[0]/k[0];
}
}
else
{
i=-1;
do
{
i=i+1;
}while(k0>k[i]);
if(i==0)
{
m=( log10(P[i]) - log10(P[i+1]) )/( log10(k[i]) - log10(k[i+1]) );
n=log10(P[i])-m*log10(k[i]);
}
else
{
m=( log10(P[i]) - log10(P[i-1]) )/( log10(k[i]) - log10(k[i-1]) );
n=log10(P[i])-m*log10(k[i]);
}
P0=pow(10,m*(k0)+n);
}

return P0;
}


double Xi_interpol(double k0, double *k, double *P, long int N)
{
double P0,m,n;
long int i;
if(k0<=k[0] || k0>=k[N-1])
{

P0=0;
}
else
{
i=-1;
do
{
i=i+1;
}while(k0>k[i]);
if(i==0)
{
m=( (P[i]) - (P[i+1]) )/( (k[i]) - (k[i+1]) );
n=(P[i])-m*(k[i]);
}
else
{
m=( (P[i]) - (P[i-1]) )/( ( k[i]) - (k[i-1]) );
n=(P[i])-m*(k[i]);
}
P0=m*(k0)+n;
}

return P0;
}


