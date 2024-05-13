#include<stdlib.h>
#include<stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>
#include<math.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define PI 3.14159265

float sinc(float x)
{ 
 if(x==0)
  return 1;
 else
  return(sin(x)/x);
}


int main (void)
{
  int i;
  const int n = 1024;
  float xmin=-500.0,xmax=500.0,dx=(xmax-xmin)/(n-1),cons=(dx/sqrt(2*PI));
  double data[2*n];

  gsl_fft_complex_wavetable * wavetable;
  gsl_fft_complex_workspace * workspace;

  for (i = 0; i < n; i++)
    {
      REAL(data,i) = sinc(xmin + i*dx);
      IMAG(data,i) = 0.0;
    }

  

  wavetable = gsl_fft_complex_wavetable_alloc (n);
  workspace = gsl_fft_complex_workspace_alloc (n);
  float *realfactor=(float*)malloc(n*sizeof(float));
  float *imagfactor=(float*)malloc(n*sizeof(float));
  float *k=(float*)malloc(n*sizeof(float));
  for (i = 0; i < (int) wavetable->nf; i++)
    {
       printf ("# factor %d: %zu\n", i,
               wavetable->factor[i]);
    }

  gsl_fft_complex_forward (data, 1, n,
                           wavetable, workspace);

  

  
  
   for (i = 0; i < n; i += 1)
 {  
   if(i<(n/2))
    k[i]=(i/(n*dx))*(2*PI);
    
   else k[i]=((i-n)/(n*dx))*(2*PI);
   
    realfactor[i]=cos(k[i]*xmin);
    imagfactor[i]=sin(-k[i]*xmin);
 }
  
 char fname[50];
 sprintf(fname,"fftgsl.txt");
 FILE *fp=fopen(fname,"w");
  for (i = 0; i < n; i += 1)
 {  
   REAL(data,i)=cons*(realfactor[i]*REAL(data,i)-imagfactor[i]*IMAG(data,i));
   IMAG(data,i)=cons*(realfactor[i]*IMAG(data,i)+imagfactor[i]*REAL(data,i));
   
   fprintf(fp,"%f  %f \n",k[i],REAL(data,i));
 }
  
  
  fclose(fp);
  gsl_fft_complex_wavetable_free (wavetable);
  gsl_fft_complex_workspace_free (workspace);
  free(k);
  free(realfactor);
  free(imagfactor);
  return 0;
}