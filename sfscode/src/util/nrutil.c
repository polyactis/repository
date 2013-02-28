/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.c.  Do not confuse this file with the same-named
   file nrutil.c that is supplied in the 'misc' subdirectory.
   *That* file is the one from the book, and contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include "nrutil.h"
#define NR_END 0
#define FREE_ARG char*
static double temp;
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp; 
static long ltemp;
#define SWAPL(a,b) ltemp=(a);(a)=(b);(b)=ltemp; 
#define M 7 
#define NSTACK 50 

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s \n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  assert(0);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
  float *v;
  
  v=(float *)calloc((nl-nl+nh+1),sizeof(float));
  if (!v) nrerror("allocation failure in vector()");
  
  return v;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;
  
  v=(int *)calloc((nl-nl+nh+1),sizeof(int));
  if (!v) nrerror("allocation failure in ivector()");
  
  return v;
}

char *cvector(long nl, long nh)
/* allocate an char vector with subscript range v[nl..nh] */
{
  char *v;
  
  v=(char *)calloc((nl-nl+nh+1),sizeof(char));
  
  if (!v) nrerror("allocation failure in cvector()");
  
  return v;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
  unsigned long *v;
  
  v=(unsigned long *)calloc((nl-nl+nh+1),sizeof(long));
  if (!v) nrerror("allocation failure in lvector()");
  
  return v;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;
  v=(double *)calloc((nl-nl+nh+1),sizeof(double));
  if (!v) nrerror("allocation failure in dvector()");
  
  return v;
  
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;
  
  /* allocate pointers to rows */
  m=(float **) calloc((nrow+NR_END),sizeof(float*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(float *) calloc((nrow*ncol+NR_END),sizeof(float));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

char **cmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate an char matrix with subscript range v[nl..nh] */
{
  char **m;
  int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;;
  
  m=(char **)calloc((nrow+NR_END),sizeof(char*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(char *) calloc((nrow*ncol+NR_END),sizeof(char));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;
  
  /* allocate pointers to rows */
  m=(double **) calloc((nrow+NR_END),sizeof(double*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  /* allocate rows and set pointers to them */
  m[nrl]=(double *) calloc((nrow*ncol+NR_END),sizeof(double));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
  int **m;
  
  /* allocate pointers to rows */
  m=(int **) calloc((nrow+NR_END),sizeof(int*));
  if (!m) nrerror("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;
  
  
  /* allocate rows and set pointers to them */
  m[nrl]=(int *) calloc((nrow*ncol+NR_END),sizeof(int));
  if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;
  
  for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

int ***icube(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate an int cube with subscript range p[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j;
  int ***p;
  
  /* //allocate for the 1d array */
  p = (int ***)calloc((nrh-nrl+1),sizeof(int**));
  if(p == NULL) nrerror("allocation failure 1 in icube()");
  
  /*   //create 2d array */
  for(i = 0; i < nrh-nrl+1; i++)
    {
      p[i] = (int**)calloc((nch-ncl+1),sizeof(int*));
      if(p[i] == NULL) nrerror("allocation failure 2 in icube()");
    }
/*   //expand that to a 3d array */
  for(i = 0; i < nrh-nrl+1; i++)
    {
      for(j = 0; j < nch-ncl+1; j++)
	{
	  p[i][j] = (int*)calloc((ndh-ndl+1),sizeof(int));
	  if(p[i][j] == NULL) nrerror("allocation failure 3 in icube()");
	}
    }
  return p;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
		  long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
  long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
  float **m;
  
  /* allocate array of pointers to rows */
  m=(float **) calloc((nrow+oldch-oldch+NR_END),sizeof(float*));
  if (!m) nrerror("allocation failure in submatrix()");
  m += NR_END;
  m -= newrl;
  
  /* set pointers to rows */
  for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;
  
  /* return pointer to array of pointers to rows */
  return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
   declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
   and ncol=nch-ncl+1. The routine should be called with the address
   &a[0][0] as the first argument. */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  float **m;
  
  /* allocate pointers to rows */
  m=(float **) calloc((nrow+NR_END),sizeof(float*));
  if (!m) nrerror("allocation failure in convert_matrix()");
  m += NR_END;
  m -= nrl;
  
  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  float ***t;
  
  /* allocate pointers to pointers to rows */
  t=(float ***) calloc((nrow+NR_END),sizeof(float**));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(float **) calloc((nrow*ncol+NR_END),sizeof(float*));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(float *) calloc((nrow*ncol*ndep+NR_END),sizeof(float));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
  /* return pointer to array of pointers to rows */
  return t;
}

int ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate an int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  int ***t;
  
  /* allocate pointers to pointers to rows */
  t = malloc((nrow+NR_END)*sizeof(*t));
  if (!t) nrerror("allocation failure 1 in i3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  for(i=nrl; i<=nrh; i++){
    t[i] = malloc((ncol+NR_END)*sizeof(**t));
    if (!t[i]) nrerror("allocation failure 2 in i3tensor()");
    t[i] += NR_END;
    t[i] -= ncl;
  }

  /* allocate rows and set pointers to them */
  for(i=nrl; i<=nrh; i++){
    for(j=ncl; j<=nch; j++){
      t[i][j] = malloc((ndep+NR_END)*sizeof(***t));
      if (!t[i][j]) nrerror("allocation failure 3 in i3tensor()");
      t[i][j] += NR_END;
      t[i][j] -= ndl;
    }  
  }
  /* return pointer to array of pointers to rows */
  return t;
}

double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  double ***t;
  
  /* allocate pointers to pointers to rows */
  t = malloc((nrow+NR_END)*sizeof(*t));
  if (!t) nrerror("allocation failure 1 in i3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  for(i=nrl; i<=nrh; i++){
    t[i] = malloc((ncol+NR_END)*sizeof(**t));
    if (!t[i]) nrerror("allocation failure 2 in i3tensor()");
    t[i] += NR_END;
    t[i] -= ncl;
  }

  /* allocate rows and set pointers to them */
  for(i=nrl; i<=nrh; i++){
    for(j=ncl; j<=nch; j++){
      t[i][j] = malloc((ndep+NR_END)*sizeof(***t));
      if (!t[i][j]) nrerror("allocation failure 3 in i3tensor()");
      t[i][j] += NR_END;
      t[i][j] -= ndl;
    }  
  }
  /* return pointer to array of pointers to rows */
  return t;
}

char ***c3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a char 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  char ***t;
  
  /* allocate pointers to pointers to rows */
  t=(char ***) calloc((nrow+NR_END),sizeof(char**));
  if (!t) nrerror("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;
  
  /* allocate pointers to rows and set pointers to them */
  t[nrl]=(char **) calloc((nrow*ncol+NR_END),sizeof(char*));
  if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  /* allocate rows and set pointers to them */
  t[nrl][ncl]=(char *) calloc((nrow*ncol*ndep+NR_END),sizeof(char));
  if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
  /* return pointer to array of pointers to rows */
  return t;
}

double ****d4tensor(long nrh, long nch, long ndh, long nwh)
/* allocate a double 4tensor with range t[0..nrh][0..nch][0..ndh][0..nwh] */
{
  long i,j,k;
  double ****t=NULL;
  
  /* allocate pointers to pointers to rows */
  assert(t = malloc((nrh+1)*sizeof(*t)));
  for(i=0; i<=nrh; i++){
    assert(t[i] = malloc((nch+1)*sizeof(**t)));
    for(j=0; j<=nch; j++){
      assert(t[i][j] = malloc((ndh+1)*sizeof(***t)));
      for(k=0; k<=ndh; k++){
	assert(t[i][j][k] = malloc((nwh+1)*sizeof(****t)));
      }
    }
  }
  /* return pointer to 4tensor of pointers to rows */
  return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
  free((FREE_ARG) (v+nl+nh-nh-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free((FREE_ARG) (v+nl+nh-nh-NR_END));
}

void free_cvector(char *v, long nl, long nh)
/* free an char vector allocated with cvector() */
{
  free((FREE_ARG) (v+nl+nh-nh-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
  free((FREE_ARG) (v+nl+nh-nh-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free((FREE_ARG) (v+nl+nh-nh-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
  free((FREE_ARG) (m[nrl]+ncl+nrh-nrh-NR_END));
  free((FREE_ARG) (m+nrl+nch-nch-NR_END));
}

void free_cmatrix(char **m, long nrl, long nrh, long ncl, long nch)
/* free a char matrix allocated by cmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl+nch-nch-NR_END));
  free((FREE_ARG) (m+nrl+nrh-nrh-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl+nch-nch-NR_END));
  free((FREE_ARG) (m+nrl+nrh-nrh-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free((FREE_ARG) (m[nrl]+ncl+nch-nch-NR_END));
  free((FREE_ARG) (m+nrl+nrh-nrh-NR_END));
}

void free_icube(int ***p,long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
{
  int i,j;
/*   //free last dimension first, if successfully allocated */
  i=ndl+ndh;
  for(i = 0; i < nrh-nrl+1; i++)
    {
      for(j = 0; j < nch-ncl+1; j++) free(p[i][j]);
    }
  
/*   //then 2d */
  for(i = 0; i < nrh-nrl+1; i++) free(p[i]);
  
/*   //then first */
  free(p);
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
  free((FREE_ARG) (b+nrl+nrh-nrh+ncl-ncl+nch-nch-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
  free((FREE_ARG) (b+nrl+nrh-nrh+nch-ncl+nch-nch-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl+ndh-ndh-NR_END));
  free((FREE_ARG) (t[nrl]+ncl+nch-nch-NR_END));
  free((FREE_ARG) (t+nrl+nrh-nrh-NR_END));
}

void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free a float i3tensor allocated by i3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl+ndh-ndh-NR_END));
  free((FREE_ARG) (t[nrl]+ncl+nch-nch-NR_END));
  free((FREE_ARG) (t+nrl+nrh-nrh-NR_END));
}

void free_d3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free a float d3tensor allocated by f3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl+ndh-ndh-NR_END));
  free((FREE_ARG) (t[nrl]+ncl+nch-nch-NR_END));
  free((FREE_ARG) (t+nrl+nrh-nrh-NR_END));
}

void free_c3tensor(char ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh)
/* free a char c3tensor allocated by f3tensor() */
{
  free((FREE_ARG) (t[nrl][ncl]+ndl+ndh-ndh-NR_END));
  free((FREE_ARG) (t[nrl]+ncl+nch-nch-NR_END));
  free((FREE_ARG) (t+nrl+nrh-nrh-NR_END));
}

void sort(unsigned long n, double arr[]) 
/*Sorts an array arr[0..n] into ascending numerical order using the 
  Quicksort algorithm. n is input; arr is replaced on output by its 
  sorted rearrangement. */
{ 
  long i,ir=n,j,k,l=0;
  unsigned long *istack; 
  int jstack=0; 
  double a; 
  
  istack=lvector(0,NSTACK); 
  
  for (;;) { 
    if (ir-l < M) { /*Insertion sort when subarray small enough. */
      for (j=l+1;j<=ir;j++) { 
	a=arr[j]; 
	for (i=j-1; i>=l; i--) { 
	  if (arr[i] <= a) break; 
	  arr[i+1]=arr[i];
	} 
	arr[i+1]=a; 
      } 
      if(jstack == 0) break; 
      ir=istack[jstack--]; /*Pop stack and begin a new round of partitioning. */
      l=istack[jstack--]; 
    }
    else { 
      k=(l+ir) >> 1;   /* Choose median of left, center, and right elements */
                       /* as partitioning element a. Also rearrange so that */
                       /* a[l] <= a[l+1] <= a [ir]. */ 
      SWAP(arr[k],arr[l+1]); 
      if(arr[l] >arr[ir]) { 
	SWAP(arr[l],arr[ir]);
      }
      if(arr[l+1] > arr[ir]) { 
	SWAP(arr[l+1],arr[ir]);
      }
      if(arr[l] > arr[l+1]) { 
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1; /* Initialize pointers for partitioning. */
      j=ir; 
      a=arr[l+1]; /* Partitioning element. */
      for (;;) {  /* Beginning of innermost loop. */
	do i++; while (arr[i] <a); /* Scan up to find element >a. */
	do j--; while (arr[j] >a); /* Scan down to find element <a. */
	if (j < i) break; /* Pointers crossed. Partitioning complete. */
	SWAP(arr[i],arr[j]); /* Exchange elements. */
      } /* End of inner most loop. */
      arr[l+1]=arr[j]; /*Insert partitioning element. */
      arr[j]=a; 
      jstack += 2; 
      /* Push pointers to larger subarray on stack, process smaller subarray 
	 immediately. */
      if(jstack >NSTACK) nrerror("NSTACK too small in sort."); 
      if(ir-i+1 >=j-l) { 
	istack[jstack]=ir; 
	istack[jstack-1]=i; 
	ir=j-1; 
      }
      else { 
	istack[jstack]=j-1; 
	istack[jstack-1]=l; 
	l=i; 
      } 
    } 
  } 
  free_lvector(istack, 0, NSTACK); 
}

void sortlong(unsigned long n, long arr[]) 
/*Sorts an array arr[0..n] into ascending numerical order using the 
  Quicksort algorithm. n is input; arr is replaced on output by its 
  sorted rearrangement. */
{ 
  long i,ir=n,j,k,l=0;
  unsigned long *istack; 
  int jstack=0; 
  double a; 
  
  istack=lvector(0,NSTACK); 
  
  for (;;) { 
    if (ir-l < M) { /*Insertion sort when subarray small enough. */
      for (j=l+1;j<=ir;j++) { 
	a=arr[j]; 
	for (i=j-1; i>=l; i--) { 
	  if (arr[i] <= a) break; 
	  arr[i+1]=arr[i];
	} 
	arr[i+1]=a; 
      } 
      if(jstack == 0) break; 
      ir=istack[jstack--]; /*Pop stack and begin a new round of partitioning. */
      l=istack[jstack--]; 
    }
    else { 
      k=(l+ir) >> 1;   /* Choose median of left, center, and right elements */
                       /* as partitioning element a. Also rearrange so that */
                       /* a[l] <= a[l+1] <= a [ir]. */ 
      SWAP(arr[k],arr[l+1]); 
      if(arr[l] >arr[ir]) { 
	SWAP(arr[l],arr[ir]);
      }
      if(arr[l+1] > arr[ir]) { 
	SWAP(arr[l+1],arr[ir]);
      }
      if(arr[l] > arr[l+1]) { 
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1; /* Initialize pointers for partitioning. */
      j=ir; 
      a=arr[l+1]; /* Partitioning element. */
      for (;;) {  /* Beginning of innermost loop. */
	do i++; while (arr[i] <a); /* Scan up to find element >a. */
	do j--; while (arr[j] >a); /* Scan down to find element <a. */
	if (j < i) break; /* Pointers crossed. Partitioning complete. */
	SWAP(arr[i],arr[j]); /* Exchange elements. */
      } /* End of inner most loop. */
      arr[l+1]=arr[j]; /*Insert partitioning element. */
      arr[j]=a; 
      jstack += 2; 
      /* Push pointers to larger subarray on stack, process smaller subarray 
	 immediately. */
      if(jstack >NSTACK) nrerror("NSTACK too small in sort."); 
      if(ir-i+1 >=j-l) { 
	istack[jstack]=ir; 
	istack[jstack-1]=i; 
	ir=j-1; 
      }
      else { 
	istack[jstack]=j-1; 
	istack[jstack-1]=l; 
	l=i; 
      } 
    } 
  } 
  free_lvector(istack, 0, NSTACK); 
}


void partialsortlong(unsigned long b, unsigned long n, long arr[]) 
/*Sorts an array arr[0..n] into ascending numerical order using the 
  Quicksort algorithm. n is input; arr is replaced on output by its 
  sorted rearrangement. */
{ 
  long i,ir=n,j,k,l=b;
  unsigned long *istack; 
  int jstack=0; 
  double a; 
  
  istack=lvector(0,NSTACK); 
  
  for (;;) { 
    if (ir-l < M) { /*Insertion sort when subarray small enough. */
      for (j=l+1;j<=ir;j++) { 
	a=arr[j]; 
	for (i=j-1; i>=l; i--) { 
	  if (arr[i] <= a) break; 
	  arr[i+1]=arr[i];
	} 
	arr[i+1]=a; 
      } 
      if(jstack == 0) break; 
      ir=istack[jstack--]; /*Pop stack and begin a new round of partitioning. */
      l=istack[jstack--]; 
    }
    else { 
      k=(l+ir) >> 1;   /* Choose median of left, center, and right elements */
                       /* as partitioning element a. Also rearrange so that */
                       /* a[l] <= a[l+1] <= a [ir]. */ 
      SWAP(arr[k],arr[l+1]); 
      if(arr[l] >arr[ir]) { 
	SWAP(arr[l],arr[ir]);
      }
      if(arr[l+1] > arr[ir]) { 
	SWAP(arr[l+1],arr[ir]);
      }
      if(arr[l] > arr[l+1]) { 
	SWAP(arr[l],arr[l+1]);
      }
      i=l+1; /* Initialize pointers for partitioning. */
      j=ir; 
      a=arr[l+1]; /* Partitioning element. */
      for (;;) {  /* Beginning of innermost loop. */
	do i++; while (arr[i] <a); /* Scan up to find element >a. */
	do j--; while (arr[j] >a); /* Scan down to find element <a. */
	if (j < i) break; /* Pointers crossed. Partitioning complete. */
	SWAP(arr[i],arr[j]); /* Exchange elements. */
      } /* End of inner most loop. */
      arr[l+1]=arr[j]; /*Insert partitioning element. */
      arr[j]=a; 
      jstack += 2; 
      /* Push pointers to larger subarray on stack, process smaller subarray 
	 immediately. */
      if(jstack >NSTACK) nrerror("NSTACK too small in sort."); 
      if(ir-i+1 >=j-l) { 
	istack[jstack]=ir; 
	istack[jstack-1]=i; 
	ir=j-1; 
      }
      else { 
	istack[jstack]=j-1; 
	istack[jstack-1]=l; 
	l=i; 
      } 
    } 
  } 
  free_lvector(istack, 0, NSTACK); 
}

void sort2(unsigned long n, double arr[], long brr[]) 
/* Sorts an array arr[0..n] into ascending order using Quicksort, */
/* while making the corresponding rearrangement of the array brr[1..n]. */
{ 
  long i, ir=n, j, k, l=0, *istack=NULL; 
  int jstack=0; 
  float a, b, temp;
  
  assert(istack=calloc(NSTACK+1, sizeof(long))); 
  for (;;) { /* Insertion sort when subarray small enough. */
    if (ir-l < M) { 
      for (j=l+1;j<=ir;j++) { 
	a=arr[j]; 
	b=brr[j]; 
	for (i=j-1;i>=l;i--) { 
	  if (arr[i] <= a) break; 
	  arr[i+1]=arr[i]; 
	  brr[i+1]=brr[i]; 
	} 
	arr[i+1]=a; 
	brr[i+1]=b; 
      } 
      if(jstack == 0)  break;
      ir=istack[jstack--];/* Pop stack and begin a new round of partitioning. */
      l=istack[jstack--]; 
    }
    else { 
      k=(l+ir) >> 1; /* Choose median of left, center and right elements 
			as partitioning element a. Also rearrange so that 
			a[l] <= a[l+1] <= a[ir]. */
      SWAP(arr[k],arr[l+1]);
      SWAPL(brr[k],brr[l+1]);
      if(arr[l] >arr[ir]) { 
	SWAP(arr[l],arr[ir]);
	SWAPL(brr[l],brr[ir]);
      }
      if(arr[l+1] > arr[ir]) { 
	SWAP(arr[l+1],arr[ir]);
	SWAPL(brr[l+1],brr[ir]);
      }
      if(arr[l] >arr[l+1]) { 
	SWAP(arr[l],arr[l+1]);
	SWAPL(brr[l],brr[l+1]);
      } 
      i=l+1; /* Initialize pointers for partitioning. */
      j=ir; 
      a=arr[l+1]; /* Partitioning element. */
      b=brr[l+1]; 
      for (;;) { /* Beginning of innermost loop. */
	do i++; while (arr[i] <a); /*Scan up to find element >a. */
	do j--; while (arr[j] >a); /* Scan down to find element <a. */
	if (j < i) break; /*Pointers crossed. Partitioning complete. */
	SWAP(arr[i],arr[j]); /* Exchange elements of both arrays. */
	SWAPL(brr[i],brr[j]);
      } /* End of innermost loop. */
      arr[l+1]=arr[j]; /* Insert partitioning element in both arrays. */
      arr[j]=a; 
      brr[l+1]=brr[j]; 
      brr[j]=b; 
      jstack += 2; 
      /* Push pointers to larger subarray on stack, process smaller 
	 subarray immediately. */
      if(jstack >NSTACK) nrerror("NSTACK too small in sort2."); 
      if(ir-i+1 >=j-l) { 
	istack[jstack]=ir; 
	istack[jstack-1]=i; 
	ir=j-1; 
      }else { 
	istack[jstack]=j-1; 
	istack[jstack-1]=l; 
	l=i; 
      }
    }
  }
  free(istack);
}

void sort2long(unsigned long n, long arr[], long brr[]) 
/* Sorts an array arr[0..n] into ascending order using Quicksort, */
/* while making the corresponding rearrangement of the array brr[1..n]. */
{ 
  long i, ir=n, j, k, l=0, *istack=NULL; 
  int jstack=0; 
  float a, b, temp;
  
  assert(istack=calloc(NSTACK+1, sizeof(long))); 
  for (;;) { /* Insertion sort when subarray small enough. */
    if (ir-l < M) { 
      for (j=l+1;j<=ir;j++) { 
	a=arr[j]; 
	b=brr[j]; 
	for (i=j-1;i>=l;i--) { 
	  if (arr[i] <= a) break; 
	  arr[i+1]=arr[i]; 
	  brr[i+1]=brr[i]; 
	} 
	arr[i+1]=a; 
	brr[i+1]=b; 
      } 
      if(jstack == 0)  break;
      ir=istack[jstack--];/* Pop stack and begin a new round of partitioning. */
      l=istack[jstack--]; 
    }
    else { 
      k=(l+ir) >> 1; /* Choose median of left, center and right elements 
			as partitioning element a. Also rearrange so that 
			a[l] <= a[l+1] <= a[ir]. */
      SWAP(arr[k],arr[l+1]);
      SWAPL(brr[k],brr[l+1]);
      if(arr[l] >arr[ir]) { 
	SWAP(arr[l],arr[ir]);
	SWAPL(brr[l],brr[ir]);
      }
      if(arr[l+1] > arr[ir]) { 
	SWAP(arr[l+1],arr[ir]);
	SWAPL(brr[l+1],brr[ir]);
      }
      if(arr[l] >arr[l+1]) { 
	SWAP(arr[l],arr[l+1]);
	SWAPL(brr[l],brr[l+1]);
      } 
      i=l+1; /* Initialize pointers for partitioning. */
      j=ir; 
      a=arr[l+1]; /* Partitioning element. */
      b=brr[l+1]; 
      for (;;) { /* Beginning of innermost loop. */
	do i++; while (arr[i] <a); /*Scan up to find element >a. */
	do j--; while (arr[j] >a); /* Scan down to find element <a. */
	if (j < i) break; /*Pointers crossed. Partitioning complete. */
	SWAP(arr[i],arr[j]); /* Exchange elements of both arrays. */
	SWAPL(brr[i],brr[j]);
      } /* End of innermost loop. */
      arr[l+1]=arr[j]; /* Insert partitioning element in both arrays. */
      arr[j]=a; 
      brr[l+1]=brr[j]; 
      brr[j]=b; 
      jstack += 2; 
      /* Push pointers to larger subarray on stack, process smaller 
	 subarray immediately. */
      if(jstack >NSTACK) nrerror("NSTACK too small in sort2."); 
      if(ir-i+1 >=j-l) { 
	istack[jstack]=ir; 
	istack[jstack-1]=i; 
	ir=j-1; 
      }else { 
	istack[jstack]=j-1; 
	istack[jstack-1]=l; 
	l=i; 
      }
    }
  }
  free(istack);
}

/* void doSWAP(double a, double b)
{
  SWAP(a, b);
  }*/
