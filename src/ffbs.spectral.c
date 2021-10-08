#include <R.h>
#include <Rmath.h>
#include <stdlib.h>
#include <fftw3.h>



/* printf("%d",j); */
/* printf("%s","\n"); */
/* printf("%f",out[j-1][1]); */
/* printf("%s","\n"); */


/* void fftc(int *n, double wr[], double wc[], int *inverse){ */
/*   fftw_complex *in, *out; */
/*   fftw_plan p; */
/*   int i; */
/*   in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * *n * *n); */
/*   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * *n * *n); */
/*   for(i=0; i<(*n * *n);i++){ */
/*     in[i][0]=wr[i]; */
/*     in[i][1]=wc[i]; */
/*   } */
/*   if(*inverse==1){ */
/*     p = fftw_plan_dft_2d(*n,*n, in, out, FFTW_FORWARD, FFTW_ESTIMATE); */
/*   } */
/*   else{ */
/*     p = fftw_plan_dft_2d(*n,*n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); */
/*   } */
/*   fftw_execute(p); */
/*   fftw_destroy_plan(p); */
/*   for(i=0; i<(*n * *n);i++){ */
/*     wr[i]=out[i][0]; */
/*     wc[i]=out[i][1]; */
/*   } */
/*   fftw_free(in); fftw_free(out); */
/* } */

void real_fft(int *n, double yh[], int *inverse, int indCos[], int indW[], int indWCon[], int *NFc){
  fftw_complex *in, *out;
  fftw_plan p;
  int i;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * *n * *n);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * *n * *n);
  if(*inverse==1){
    for(i=0; i<(*n * *n);i++){
      in[i][0]=yh[i];
      in[i][1]=0.0;
    }
    p = fftw_plan_dft_2d(*n,*n, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<4;i++){
      yh[i]=out[indW[i]-1][0]/ *n;
    }
    for(i=0; i<*NFc;i++){
      yh[indCos[i]-1]=sqrt(2) * out[indW[indCos[i]-1]-1][0]/ *n;
      yh[indCos[i]]=-sqrt(2) * out[indW[indCos[i]]-1][1]/ *n;
    }
  }else{
    /* for(i=0; i<(*n * *n);i++){ */
    /*   in[i][0]=0; */
    /*   in[i][1]=0; */
    /* } */
    for(i=0; i<4;i++){
       in[indW[i]-1][0]=yh[i];
       in[indW[i]-1][1]=0;
    }
    for(i=0; i<*NFc;i++){
      in[indW[indCos[i]-1]-1][0]=yh[indCos[i]-1]/sqrt(2);
      in[indW[indCos[i]-1]-1][1]=-yh[indCos[i]]/sqrt(2);
      in[indWCon[i]-1][0]=yh[indCos[i]-1]/sqrt(2);
      in[indWCon[i]-1][1]=yh[indCos[i]]/sqrt(2);
    }
    p = fftw_plan_dft_2d(*n,*n, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(p);
    for(i=0; i<(*n * *n);i++){
      yh[i]=out[i][0]/ *n;
    }
  }
/* printf("      *********         test 1       *********          "); */
  if (NULL != in) fftw_free(in);
  if (NULL != out) fftw_free(out);
  if (NULL != p) fftw_destroy_plan(p);


}

void TSreal_fft(int *n, int *T, double yh[], int *inverse, int indCos[], int indW[], int indWCon[], int *NFc){
  int t;
  for(t=0; t<*T;t++){
    real_fft(n, &yh[t * *n * *n], inverse, indCos, indW, indWCon, NFc);
  }
}


/* void kf2x2diagC(double res[], double spec[], double G11[], double G12[], double *tau2, int *T, int *NFc){ */
/*   int t,j; */
/*   double rtt,rtt1; */
/*   for(j=1; j<=*NFc;j++){ */
/*     rtt=spec[j-1]; */
/*     rtt1=0; */
/*     for(t=1; t<=*T;t++){ */
/*       rtt1=rtt*(pow(G11[j-1],2)+pow(G12[j-1],2))+spec[j-1]; */
/*       res[(j-1)* *T + t-1]=rtt1; */
/*       rtt=rtt1*(1-rtt1/(*tau2+rtt1)); */
/*       /\* Rprintf(&t); *\/ */
/*     } */
/*    } */
/* } */

void propagate_spectral(double xtp1[], double xt[], double G11C[], double G11[], double G12[], int *NFc, int *ns){
  int j;
  for(j=0; j<*ns;j++){
    xtp1[j]=G11C[j] * xt[j];
  }
  for(j=0; j<*NFc;j++){
    xtp1[2 * j+*ns]=G11[j] * xt[2 * j+*ns] + G12[j] * xt[2 * j+*ns+1];
    xtp1[2 * j+*ns + 1]=G11[j] * xt[2 * j+*ns + 1] - G12[j] * xt[2 * j+*ns];
  }
}

void kf_spectral(double wFT[], double mtt1[], double mtt[], double rtt1[], double rtt[], double specCosOnly[], double G11C[], double specCosSine[], double G11[], double G12[], double *tau2, int *T, int *NFc, int *ns){
  int t,j, NF;
  NF=(2 * *NFc + *ns);
  double rtt1Temp,rttTemp;
  for(j=0; j<*ns;j++){
    rttTemp=specCosOnly[j];
    rtt1Temp=0;
    rtt[j]=rttTemp;
    for(t=0; t<*T;t++){
      rtt1Temp=rttTemp*(pow(G11C[j],2))+specCosOnly[j];
      rtt1[j + t*NF]=rtt1Temp;
      rttTemp=rtt1Temp*(1-rtt1Temp/(*tau2+rtt1Temp));
      rtt[j + (t+1)*NF]=rttTemp;
    }
  }
  for(j=0; j<*NFc;j++){
    rttTemp=specCosSine[j];
    rtt1Temp=0;
    rtt[*ns + 2*j]=rttTemp;
    rtt[*ns + 2*j+1]=rttTemp;
    for(t=0; t<*T;t++){
      rtt1Temp=rttTemp*(pow(G11[j],2)+pow(G12[j],2))+specCosSine[j];
      rtt1[*ns + 2*j + t*NF]=rtt1Temp;
      rtt1[*ns + 2*j + 1 +t*NF]=rtt1Temp;
      rttTemp=rtt1Temp*(1-rtt1Temp/(*tau2+rtt1Temp));
      rtt[*ns + 2*j + (t+1)*NF]=rttTemp; 
      rtt[*ns + 2*j + 1 +(t+1)*NF]=rttTemp; 
    }
  }
  for(t=0; t<*T;t++){
    propagate_spectral(&mtt1[t * NF], &mtt[t * NF], G11C, G11, G12, NFc, ns);
    for(j=0; j<NF;j++){
      mtt[(t+1) * NF+j] = mtt1[t * NF+j]+(rtt1[t * NF+j] * (wFT[t * NF+j]-mtt1[t * NF+j])/(*tau2+rtt1[t * NF+j]));
    }
  }
}

void bs_spectral(double simAlpha[], double mtt[], double mtt1[], double rtt[], double rtt1[], double spec[], double G11C[], double G11[], double G12[], int *T, int *NFc, int *ns){
  int NF=(2 * *NFc + *ns);
  int j,t;
  double *AlMinusMtt1 = (double *) malloc(NF * sizeof(*AlMinusMtt1));
  double *Prop = (double *) malloc(NF * sizeof(*Prop));
  double simTemp, mt, RttBar;
  double *G12t = (double *) malloc(NF * sizeof(*G12t));
  for(j=0; j<*NFc;j++){
    G12t[j]=-G12[j];
  }
  for(j=0; j<NF;j++){
    simTemp=mtt[*T * NF +j]+sqrt(rtt[*T * NF +j])*rnorm(0, 1);
    simAlpha[(*T-1) * NF +j] = simTemp;
    AlMinusMtt1[j]=simTemp-mtt1[(*T-1) * NF +j];
  }
  for(t=(*T-1); t>0;t--){
    propagate_spectral(Prop, AlMinusMtt1, G11C, G11, G12t, NFc, ns);
    for(j=0; j<NF;j++){
      mt = mtt[t * NF +j]+rtt[t * NF +j]/rtt1[t * NF +j]*Prop[j];
      RttBar = rtt[t * NF +j] * (1-(rtt1[t * NF +j]-spec[j])/rtt1[t * NF +j]);
      simTemp=mt+sqrt(RttBar)*rnorm(0, 1);
      simAlpha[(t-1) * NF +j] = simTemp;
      AlMinusMtt1[j]=simTemp-mtt1[(t-1) * NF +j];
    }
  }
}

void ll_spectral(double *ll, double wFT[], double mtt1[], double rtt1[],  int *T, int *NF, double *tau2){
  *ll=0;
  int t,j;
  for(t=0; t<*T;t++){
    for(j=0; j<*NF;j++){
      *ll = *ll-log(*tau2+rtt1[t* *NF +j]) - (wFT[t* *NF +j]-mtt1[t* *NF +j])*(wFT[t* *NF +j]-mtt1[t* *NF +j])/(*tau2+rtt1[t* *NF +j]);
	}
      }
  *ll = *ll/2-*T * *NF * log(2*M_PI)/2;
}

void ffbs_spectral(double wFT[], double *bw, double *ll, double specCosOnly[], double G11C[], double specCosSine[], double G11[], double G12[], double specAll[], double *tau2, int *T, int *NFc, int *ns){/* , double simAlpha[] */
  int NF=(2 * *NFc + *ns);
  int j;
  double *rtt = (double *) malloc(NF * (*T +1)*sizeof(*rtt));
  double *rtt1 = (double *) malloc(NF * (*T)*sizeof(*rtt1));
  double *mtt = (double *) malloc(NF * (*T+1)*sizeof(*mtt));
  double *mtt1 = (double *) malloc(NF * (*T)*sizeof(*mtt1));
  for(j=0; j<NF;j++){
    mtt[j]=0.0;
  }
  kf_spectral(wFT, mtt1, mtt, rtt1, rtt, specCosOnly, G11C, specCosSine, G11, G12, tau2, T, NFc, ns);
  if(*ll==1){
    ll_spectral(ll, wFT, mtt1, rtt1, T, &NF, tau2);
  }
  if(*bw==1){
    bs_spectral(wFT, mtt, mtt1, rtt, rtt1, specAll, G11C, G11, G12, T, NFc, ns);
  }
  /* wFT=simAlpha; */
  /* printf("%f",*ll); , double rtt[], double rtt1[], double mtt[], double mtt1[] */
 free(rtt);
 free(rtt1);
 free(mtt);
 free(mtt1);
 rtt=NULL;
 rtt1=NULL;
 mtt=NULL;
 mtt1=NULL;
}




/* void ffbsspectralC(double rtt1[], double specCosOnly[], double G11C[], double spec[], double G11[], double G12[], double *tau2, int *T, int *NFc, int *ns){ */
/*   int t,j, NF; */
/*   NF=(2 * *NFc + *ns); */
/*   double rtt1Temp,rttTemp; */
/*   /\* double *rtt1 = (double *) malloc(NF * (*T+1)*sizeof(*rtt1)); *\/ */
/*   double *rtt = (double *) malloc(NF * (*T+1)*sizeof(*rtt)); */
/*   for(j=1; j<=*ns;j++){ */
/*     for(t=1; t<=*T;t++){ */
/*       rttTemp=specCosOnly[j-1]; */
/*       rtt1Temp=0; */
/*       for(t=1; t<=*T;t++){ */
/*   	rtt1Temp=rttTemp*(pow(G11C[j-1],2))+specCosOnly[j-1]; */
/*   	rtt1[(j-1) * *T +t-1]=rtt1Temp; */
/*   	rttTemp=rtt1Temp*(1-rtt1Temp/(*tau2+rtt1Temp)); */
/*   	rtt[(j-1) * *T +t-1]=rttTemp; */
/*       } */
/*     } */
/*   } */
/*   for(j=1; j<=*NFc;j++){ */
/*     rttTemp=spec[j-1]; */
/*     rtt1Temp=0; */
/*     for(t=1; t<=*T;t++){ */
/*       rtt1Temp=rttTemp*(pow(G11[j-1],2)+pow(G12[j-1],2))+spec[j-1]; */
/*       rtt1[*ns * *T + 2*(j-1)* *T+t-1]=rtt1Temp; */
/*       rtt1[*ns * *T + 2*(j-1)* *T+*T+t-1]=rtt1Temp; */
/*       rttTemp=rtt1Temp*(1-rtt1Temp/(*tau2+rtt1Temp)); */
/*       rtt[*ns * *T + 2*(j-1)* *T+t-1]=rttTemp;  */
/*       rtt[*ns * *T + 2*(j-1)* *T+*T+t-1]=rttTemp;  */
/*     } */
/*   } */
/* } */
