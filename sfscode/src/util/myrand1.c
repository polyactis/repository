#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define IA 16807					/* from ran1 */
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-07
#define RNMX (1.0-EPS)
#define PI 3.141592654
#define E 2.718281828

float ran1(long *idum);
float poidev(float xm, long *idum);
int geomdev(float p, long *idum);
float bnldev(float pp, int n, long *idum) ;
float gammln(float xx);
double rexp(double lambda,long *idum);
double rgamma1(double alpha, long *idum);
double rgamma2(double alpha, long *idum);
double rgamma(double alpha, double beta, long *idum);
double rstd_normal(long *idum);
double rnormal(double mean, double sd, long *idum);
void moment(double data[], int n, float *ave, float *var);

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  do{
    if (*idum <= 0 || !iy){
      if (-(*idum) < 1) 
	*idum=1;
      else 
	*idum = -(*idum);
      for (j=NTAB+7;j>=0;j--){
	k=(*idum)/IQ;
	*idum=IA*(*idum-k*IQ)-IR*k; 
	if (*idum < 0) *idum += IM; 
	if (j < NTAB) 
	  iv[j] = *idum;
      }
      iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if (*idum < 0) 
      *idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    temp = AM*iy;
  }while(temp <= EPS || temp >= RNMX);
  return temp;
}

float poidev(float xm, long *idum)
{
  static float sq,alxm,g,oldm=(-1.0);
  float em,t,y;

  if (xm < 12.0) 
    {
      if (fabsf(xm - oldm)>EPS) 
	{
	  oldm=xm;
	  g=(float)exp(-xm);
	}
      em = -1;
      t=1.0;
      do {
	++em;
	t *= ran1(idum);
      } while (t > g);
    } 
  else 
    {
      if (fabsf(xm - oldm)>EPS) 
	{
	  oldm=xm;
	  sq=(float)sqrt(2.0*xm);
	  alxm=(float)log(xm);
	  g=xm*alxm-gammln((float)(xm+1.0));
	}
      do {
	do {
	  y=(float)tan(PI*ran1(idum));
	  em=sq*y+xm;
	} while (em < 0.0);
	em=(float)floor(em);
	t=(float)(0.9*(1.0+y*y)*exp(em*alxm-gammln((float)(em+1.0))-g));
      } while (ran1(idum) > t);
    }
  return em;
}

int geomdev(float p, long *idum)
/* this follows from simple inverse-CDF method outlined in Ross 2001 */
/* has mean 1/p, support={1,2,...} */
{
  double u;
  u = log(ran1(idum));
  u /= log(1-p);
  return(1+(int)u);
}

float gammln(float xx)
{
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
			24.01409824083091,-1.231739572450155,
			0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return (float)(-tmp+log(2.5066282746310005*ser/x));
}

double rexp(double lambda,long *idum)
/*
 * Generates from an exponential distribution
 */
{
  double random, uniform;
  uniform = ran1(idum);
  random = - (1/lambda) * log(uniform);
  return random;
}


double rgamma1(double alpha, long *idum)
/* 
 * Generates from a gamma distribution with alpha < 1
 */
{
  double uniform0, uniform1;
  double random, x;
  uniform0 = ran1(idum);
  uniform1 = ran1(idum);
  if (uniform0 > E/(alpha + E))
    {
      random = -log((alpha + E)*(1-uniform0)/(alpha*E));
      if ( uniform1 > pow(random,alpha - 1))
	return -1;
      else 
	return random;
    }
  else
    {
      x = (alpha + E) * uniform0 / E;
      random = pow(x,1/alpha);
      if ( uniform1 > exp(-random))
	return -1;
      else
	return random;
    } 
}

double rgamma2(double alpha, long *idum)
/*
 * Generates from a gamma distribution with alpha > 1
 */
{
  double uniform1,uniform2;
  double c1,c2,c3,c4,c5,w,foo;
  double random;
  int done = 1;
  c1 = alpha - 1;
  c2 = (alpha - 1/(6 * alpha))/c1;
  c3 = 2 / c1;
  c4 = c3 + 2;
  c5 = 1 / sqrt(alpha);
  do
    {
      uniform1 = ran1(idum);
      uniform2 = ran1(idum);
      if (alpha > 2.5)
	{
          uniform1 = uniform2 + c5 * (1 - 1.86 * uniform1);
	}
    }
  while ((uniform1 >= 1) || (uniform1 <= 0));
  w = c2 * uniform2 / uniform1;
  if ((c3 * uniform1 + w + 1/w) > c4)
    {
      foo=c3 * log(uniform1);
      foo+=w;
      foo-=log(w);
      if (foo  >= 1)
	{
	  done = 0;
	}
    }
  if (done == 0)
    return -1;
  random = c1 * w; 
  return random;
} 

double rgamma(double alpha, double beta, long *idum)
/*
 * Generates from a general gamma(alpha,beta) distribution
 */   
{
  double random = 0.0;
  if (alpha < 1)
    do {
      random = rgamma1(alpha,idum)/beta; 
    } while (random < 0 );
  if (fabs(alpha - 1)<EPS)
    random = rexp(1,idum)/beta; 
  if (alpha > 1)
    do {
      random = rgamma2(alpha,idum)/beta; 
    } while (random < 0);
  return random;
}

double rstd_normal(long *idum)
/*
 * Generates from a standard normal(0,1) distribution
 */ 
{
  double uniform1,uniform2;
  double theta,r;
  double random;
  uniform1 = ran1(idum);
  uniform2 = ran1(idum);
  theta = 2 * PI * uniform1;
  r = sqrt(2 * ( - log(uniform2)));
  random = r * cos(theta);
  return random;
}
   
double rnormal(double mean, double sd, long *idum)
/*
 * Generates from a general normal(mu,sigma) distribution
 */
{
  double random;
  random = mean + sd * rstd_normal(idum);
  return random;
}

void moment(double data[], int n, float *ave, float *var)
{
  int j;
  float s;
  
  if (n <= 1) 
    {    
      printf("n must be at least 2 in moment");
      exit(1);
    }
  s=0.0;
  (*var) = 0.0;
  for (j=1;j<=n;j++){
    s += (float)data[j];
    *var += (float)(data[j]*data[j]/n);
  }
  printf("\n");
  *ave=s/n;
  *var=(*var-(*ave)*(*ave))*n/(n-1);
}

float bnldev(float pp, int n, long *idum) 
/* Returns as a floating-point number an integer value that is a 
   random deviate drawn from a binomial distribution of n trials 
   each of probability pp, using ran1(idum) as asource of 
   uniform random deviates. */
{ 
  int j; 
  static int nold=(-1); 
  float am,em,g,angle,p,bnl,sq,t,y; 
  static float pold=(-1.0),pc,plog,pclog,en,oldg; 
  p=(float)(pp <= 0.5 ? pp : 1.0-pp); 
  /* The binomial distribution is invariant under changing pp to 1-pp, 
     if we also change the answer to n minus itself; we'll remember 
     to do this below. */
  am=n*p; /* This is the mean of the deviate to be produced. */
  if (n < 25) { /* Use the direct method while n is not too large. 
		   This can require up to 25 calls to ran1. */
    bnl=0.0; 
    for (j=1; j<=n; j++) 
      if(ran1(idum) < p) ++bnl; 
  }else if (am <1.0) { /* If fewer than one event is expected out of 25 
			  or more trials, then the distribution is quite 
			  accurately Poisson. Use direct Poisson method. */
    g=expf(-am); 
    t=1.0; 
    for (j=0;j<n;j++) { 
      t*= ran1(idum); 
      if(t < g) break; 
    } 
    bnl = j + .0; 
  }else { /* Use the rejection method.*/
    if (n != nold) { /*If n has changed, then compute useful quantities. */
      en = n; 
      oldg = gammln(en+1.0); 
      nold=n; 
    }
    if (fabs(p - pold) > EPS) { /*If p has changed, then compute useful
				  quantities.*/
      pc=1.0-p; 
      plog=logf(p); 
      pclog=logf(pc); 
      pold=p; 
    } 
    sq=sqrtf(2.0*am*pc); /* rejection method with a Lorentzian comparison
			    function. */
    do { 
      do{ 
	angle=PI*ran1(idum); 
	y=tanf(angle); 
	em=sq*y+am; 
      }while (em < 0.0 || em >= (en+1.0)); /* Reject. */
      em=floor(em); /* Trick for integer-valued distribution. */
      t=1.2*sq*(1.0+y*y)*expf(oldg-gammln(em+1.0)
			     -gammln(en-em+1.0)+em*plog+(en-em)*pclog);
    }while (ran1(idum) > t); /* Reject. This happens about 1.5 times 
				per deviate, on average. */
    bnl=em; 
  } 
  if (fabs(p - pp) > EPS) bnl=n-bnl; 
  /*Remember to undo the symmetry transformation.*/
  return bnl; 
} 

#undef PI
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
