/* #define WINDOWS // uncomment this if using windows intel compiler */

/*
  TO DO:

  * allow option of only printing subset of simulations.  can do this by printing to a tmp file named using initial seed.

  * CNVs (keep list of all occurrances w/frequency and impose stepping stone model for increase/decrease frequency, and allow for novel inserts at some rate per generation

  * insertion/deletion from motif files (e.g., alus)

  * use openMP to take advantage of multi-core systems

*/


/* PROGRAM:  sfs_code
   AUTHOR:   Ryan D. Hernandez
   
   DESCRIPTION:    Simulate sequence data under a realistic mutation model
   with recombination, selection, and a few other bells and whistles   

   Copyright (C) 2007  Ryan D. Hernandez
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/* #define UNIT_TEST        /\* Many debugging 'features' added *\/ */
/* #define VERBOSE_DEBUG    /\* print out many details at each generation *\/ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#ifndef WINDOWS
#include <sys/time.h>
#endif
#include <time.h>
#include <assert.h>
#include <float.h>

#include "../util/myrand1.c"
#include "../util/nrutil.c"
#include "../util/nrutil.h"
#include "sfs_code.h"

/* make sure that small numbers defined */
#ifndef FLT_EPSILON  /* prefer to use system default */
#define FLT_EPSILON 1.2e-07
#endif
#ifndef DBL_EPSILON  /* prefer to use system default */
#define DBL_EPSILON 2.3e-16
#endif
#ifndef LONG_MAX
#define LONG_MAX 2147483647
#endif

/* DEFINE GLOBAL VARIABLES */
/* FILE *TMP; */
FILE *outfile;
FILE *errfile;
FILE *tmpOUT=NULL;
char tmpOUTfile[100];
FILE *tmpTRAJ=NULL;
char tmpTRAJfile[100];
FILE *freqFile=NULL;
FILE *ancFile=NULL;
FILE *trajFile=NULL;
struct POPparams *ppars;  /* population specific parameters */
struct GLOBparams gpars;  /* global parameters */
struct EvolEvents *devents=NULL;  /* list of demographic events */
struct mutStorage *storage=NULL;  /* store ancestral information for 
				     successive iterations */
struct population *popn;      /* information about each population */
struct mutArray **mutAr=NULL; /* global structure containing all mutations */
struct history *trash=NULL;   /* freed nodes; reduce calls to malloc */
struct deadList *dead=NULL;   /* store dead lineages during gen, free at end */
int ITERFAILED=0;
int INFREQRANGE=1;
int *selfed=NULL;
int *selfing=NULL;
int *compR=NULL;
long **parStore=NULL;
char **ancestralSEQ=NULL;  /* ancestral seq stored for successive iterations */
double **Q;  /* context-dependent mutation rates */
long INITIALSEED = 0;

/* these inputs require global variables above */
#include "gencontextrate.h"
#include "btfuncs.c"

int main(int argc, char *argv[])
{
  long i, j, k, r;
  long numMutSeg=0, trials=0;
#ifndef WINDOWS
  struct timeval time_now;
#endif
  double *ST;            /* absolute mut rates, transition probs, 
			    stat freqs, rel. mut rates, hit probs */
  double *fitQuant;
  double **mutClass_site;      /* array indicating rate class for each site */
  double *mutClass_locus;      /*  */
  int *selClass;               /* indicates selected loci (1) or not (0) */
  struct EvolEvents *tmpStore = NULL, *tmp;

/*   assert(TMP = fopen("rec.txt","w")); */
  outfile = stdout;
  errfile = stderr;
  if(argc < 3 || !isdigit(argv[1][0]) || !isdigit(argv[2][0])){
    if(argc>1 && (strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0))
      helpMenu();
    else{
      fprintf(stderr,"usage:  sfs_code <NPOP> <ITERATIONS> [options]\n");
      exit(1);
    }
  }
  
#ifdef UNIT_TEST
  fprintf(errfile,"********** RUNNING IN DEBUG MODE **********\n");
#endif /* UNIT_TEST */
    
  /* must specify seed at command-line in WINDOWS mode */
#ifndef WINDOWS
  /*generate default seed from system time (using millionth of a second inc.)*/
  gettimeofday(&time_now, NULL);
  gpars.seed = -(time_now.tv_usec*747);
  if (gpars.seed>=0)
    gpars.seed = -gpars.seed-67;
#else
  gpars.seed = 0;
#endif

  /* GENERATE ARRAY OF CONTEXT-DEPENDENT MUTATION RATES */
  Q=dmatrix(0,63,0,3);
  assert(ST = calloc(64,sizeof(*ST)));
  assert(fitQuant = calloc(100,sizeof(*fitQuant)));
  GenContextRate(ST,fitQuant);

  /* 1st command-line argument is the number of populations */
  gpars.NPOP = atoi(argv[1]);
  /* 2nd command-line argument is the total number of iterations*/
  gpars.ITER = atoi(argv[2]);
  
  for(gpars.iter=0;gpars.iter<gpars.ITER;gpars.iter++){
    trials = 0;
    /* re-initialize parameters every iteration */
    assert(ppars = malloc(gpars.NPOP*sizeof(*ppars)));
    setDefaultGPARS();
    for(i=0; i<gpars.NPOP; i++)
      setDefaultPPARS(i);

    popn[0].conSeq = NULL;

    /* parse command line options */
    if(argc > 3){
      parseCommandSFSCODE(argc, argv, gpars.iter);
#ifdef WINDOWS
      if(gpars.seed == 0){
	fprintf(errfile,"must specify seed at command-line in WINDOWS mode.  \
If you don't want to be in WINDOWS mode, comment out the first line of \
sfs_code.c.\n");
	abort();
      }
#endif
#ifdef UNIT_TEST
      printEvolEvents(devents);
      checkEvolEventTimes(devents, -1.0, gpars.P*ppars[0].N);
/*       getchar(); */
#endif
      if(INITIALSEED == 0)
	INITIALSEED = gpars.seed;
      if(gpars.trajMaxFreq<1-FLT_EPSILON || gpars.trajMinFreq>FLT_EPSILON){
	sprintf(tmpOUTfile,"tmpOUT%ld",INITIALSEED);
	if(!(tmpOUT = fopen(tmpOUTfile,"w"))){
	  fprintf(errfile,"cannot write to %s\n",tmpOUTfile);
	  abort();
	}
	if(trajFile != NULL){
	  sprintf(tmpTRAJfile,"tmpTRAJ%ld",INITIALSEED);
	  if(!(tmpTRAJ = fopen(tmpTRAJfile,"w"))){
	    fprintf(errfile,"cannot write to %s\n",tmpTRAJfile);
	    abort();
	  }
	}
      }
    }
    
    /* adjust sequence lengths of coding regions for full codons */
    assert(gpars.sumL = malloc(gpars.R*sizeof(*gpars.sumL)));
    assert(gpars.sumLL = malloc(gpars.R*sizeof(*gpars.sumLL)));
    gpars.sumL[0] = 0.0;
    gpars.sumLL[0] = 0.0;
    gpars.totalNucs = 0.0;
    k = 0;
    compR = ivector(0, gpars.R-1);
    for(r=0; r<gpars.R; r++){
      compR[r] = 0;
      if(ppars[0].selDistType[r] != 0 || ppars[0].f0[r] < 1-FLT_EPSILON)
	k = 1;
      if(r>0){
	gpars.sumL[r] = gpars.sumL[r-1];
	gpars.sumLL[r] = gpars.sumLL[r-1];
      }
      if(gpars.ANNOTATE[r] == 'C'){
	if(gpars.L[r]%3 == 1){
	  gpars.L[r]--;
	  if(gpars.SKIPSITES != NULL){
	    assert(gpars.SKIPSITES[r] =
		   realloc(gpars.SKIPSITES[r], gpars.L[r]*sizeof(int)));
	  }
	}
	else if(gpars.L[r]%3 == 2){
	  gpars.L[r]++;
	  if(gpars.SKIPSITES != NULL){
	    assert(gpars.SKIPSITES[r] =
		   realloc(gpars.SKIPSITES[r], gpars.L[r]*sizeof(int)));
	    gpars.SKIPSITES[r][gpars.L[r]-1] = 0;
	  }
	}
      }
      gpars.sumL[r] += gpars.L[r];
      gpars.sumLL[r] += gpars.L[r];
      if(gpars.R > 1 && gpars.LinkDistType == 'p' && r < gpars.R-1){
	if(gpars.LINK[r] >= 0)
	  gpars.sumLL[r] += gpars.LINK[r];
	else{
	  if(gpars.recMap != NULL){
	    fprintf(errfile,"SFS_CODE error: cannot have independent loci \
while using the recombination map (loci %ld-%ld seem to be independent.\n",
		    r, r+1);
	    abort();
	  }
	}
      }
      if(r > 0)
	if(gpars.L[r] != gpars.L[r-1])
	  gpars.EQL = '0';
      if(popn[0].conSeq != NULL && popn[0].conSeq[r] != NULL){
	for(j=0; j<gpars.L[r]; j++){
	  if(popn[0].conSeq[r][j] != '4')
	    gpars.totalNucs++;
	}
      }
      else{
	gpars.totalNucs = gpars.sumL[r];
      }
    }

    if(k == 0)
      ppars[0].neutpop = 1;
    if(gpars.recMap != NULL && 
       gpars.sumLL[gpars.R-1] != gpars.recMapPos[gpars.mapPoints-1]){
      fprintf(errfile,"error, rec. map does not match sequence length!\n");
      fprintf(errfile,"map = %ld basepairs, sum(L) = %.0lf\n",
	      gpars.recMapPos[gpars.mapPoints-1], gpars.sumLL[gpars.R-1]);
      abort();
    }
    /* print command line to screen on first iteration only */
    if(gpars.iter == 0){
      if(tmpOUT != NULL){
	for(j=0; j<argc; j++)
	  fprintf(tmpOUT,"%s ",argv[j]);
	fprintf(tmpOUT,"\n");
	fprintf(tmpOUT,"SEED = %ld\n",INITIALSEED);
	fflush(tmpOUT);
      }
      else{
	for(j=0; j<argc; j++)
	  fprintf(outfile,"%s ",argv[j]);
	fprintf(outfile,"\n");
	fprintf(outfile,"SEED = %ld\n",INITIALSEED);
	fflush(outfile);
      }
    }
    
    gpars.maxN = ppars[0].N; /* largest population size */
    assert(mutClass_site = malloc(gpars.R*sizeof(*mutClass_site)));
    for(r=0; r<gpars.R; r++)
      assert(mutClass_site[r] = malloc(gpars.L[r]*sizeof(**mutClass_site)));

    assert(mutClass_locus = malloc(gpars.R*sizeof(*mutClass_locus)));
    assert(selClass=calloc(gpars.R,sizeof(*selClass)));

    /*initialize structures, note that Rtree points to beginning of information,
     * root defined here is just a pointer to it. */
    assert(popn[0].BigHead = malloc(sizeof(struct history**)));
    assert(popn[0].BigHead[0] = malloc(gpars.R*sizeof(struct history*)));
    for(r=0; r<gpars.R; r++){
      popn[0].BigHead[0][r] = popTrash();
    }

    if(tmpOUT != NULL){
      fprintf(tmpOUT,"//iteration:%i/%i\n",gpars.iter+1,gpars.ITER);
      fflush(tmpOUT);
    }
    else{
      fprintf(outfile,"//iteration:%i/%i\n",gpars.iter+1,gpars.ITER);
      fflush(outfile);
    }
    
    if(outfile != stdout){
      fprintf(errfile,"//iteration:%i/%i",gpars.iter+1,gpars.ITER);
      if(gpars.iter == 0)
	fprintf(errfile," (seed = %ld)",INITIALSEED);
      fprintf(errfile,"\n");
      fflush(errfile);
    }
    
    /* CALL FUNCTIONS TO RUN PROGRAM */
    if(gpars.iter == 0)
      Initialize(ST);
    else{
      gpars.BURN = gpars.BURN2; /* burn-in for successive iterations may differ
				   from initial burn-in */
      j = 0;
      if(popn[0].conSeq == NULL){
	assert(popn[0].conSeq = malloc(gpars.R*sizeof(char*)));
	j = 1;
      }
      for(r=0; r<gpars.R; r++){
	if(j == 1)
	  assert(popn[0].conSeq[r] = malloc((gpars.L[r]+1)*sizeof(char)));
#ifdef WINDOWS
	strcpy_s(popn[0].conSeq[r], gpars.L[r]+1, ancestralSEQ[r]);
#else
	strcpy(popn[0].conSeq[r], ancestralSEQ[r]);
#endif
      }
    }

#ifdef UNIT_TEST
    fprintf(errfile,"sequence initialized\n");
    fflush(errfile);
#endif

    /* now do all the work... */
    bign(mutClass_site, mutClass_locus, selClass, fitQuant, &numMutSeg);
    
#ifdef UNIT_TEST
    printf("\nmutations spaces taken:\n");
    for(r=0; r<gpars.R; r++){ /* all should be cleared */
      printf("locus %ld : %ld\n", r, mutAr[r]->numMuts);
    }
    fflush(stdout);
    printf("\nmaxchr:\n");
    for(r=0; r<gpars.NPOP; r++){ /* all should be cleared */
      printf("population %ld : %ld\n", r, popn[0].maxchr);
    }
    fflush(stdout);
#endif
    /* free everything but ancestral information used to seed next iteration */
    while(tmpStore != NULL){
      tmp = tmpStore->nextEvent;
      free(tmpStore);
      tmpStore = tmp;
    }
    for(r = 0; r<gpars.R; r++){
      pushTrash(&popn[0].BigHead[0][r]->Rtree);
      popn[0].BigHead[0][r]->Rtree = NULL;
      free(popn[0].BigHead[0][r]);
      free(mutClass_site[r]);
      free(popn[0].conSeq[r]);
      if(gpars.SKIPSITES != NULL)
	free(gpars.SKIPSITES[r]);
      if(gpars.TRACKANC)
	free(popn[0].ancestry[0][r]);
    }
    free(popn[0].BigHead[0]);
    free(popn[0].BigHead);
    popn[0].BigHead = NULL;
    free(popn[0].conSeq);
    popn[0].conSeq = NULL;
    if(gpars.TRACKANC){
      free(popn[0].ancestry[0]);
      free(popn[0].ancestry);
    }
    if(gpars.SKIPSITES != NULL)
      free(gpars.SKIPSITES);
    free(popn);
    popn = NULL;
    free(mutClass_site);
    mutClass_site = NULL;
    
    for(i=0; i<gpars.NPOP; i++){
      free(gpars.mig_mat[i]);
      free(ppars[i].GAMMA);
      free(ppars[i].ProbPos);
      free(ppars[i].ProbDel);
      free(ppars[i].ProbNeut);
      free(ppars[i].selDistType);
      free(ppars[i].ProbNegGamma);
      free(ppars[i].alphaN);
      free(ppars[i].alphaP);
      free(ppars[i].lambdaN);
      free(ppars[i].lambdaP);
      free(ppars[i].normMean);
      free(ppars[i].normVar);
      free(ppars[i].RateParamLoci);
      free(ppars[i].f0);
    }
    free(gpars.mig_mat);
    free(mutClass_locus);
    free(gpars.L);
    free(gpars.LINK);
    free(gpars.ANNOTATE);
    free(gpars.SEX);
    freeEvolEvents(devents);
    devents = NULL;
    free(ppars);
    free(selClass);
    free(gpars.sumL);
    free(gpars.sumLL);
    free(gpars.recMap);
    free(gpars.recMapPos);
    free(compR);
    freeMutArray();
    mutAr = NULL;
    while(dead->next != NULL){
      if(dead->next->next == NULL){
	free(dead->next);
	dead->next = NULL;
	break;
      }
      else{
	dead->next = dead->next->next;
	free(dead->next->prev);
	dead->next->prev = NULL;
      }
    }

    if(ITERFAILED || !INFREQRANGE){
      gpars.iter--;
      ITERFAILED = 0;
      trials++;
      if(trials >= gpars.trajMaxReps){
	fprintf(errfile,"maximum iterations reached.  could not achieve sufficient number of replicates.\n");
	fflush(errfile);
	abort();
      }
    }
    else{
      if(tmpOUT != NULL){
	fclose(tmpOUT);
	if(!(tmpOUT = fopen(tmpOUTfile,"r"))){
	  fprintf(errfile,"cannot read %s after simulation\n",tmpOUTfile);
	  abort();
	}
	else{
	  char c;
	  while((c=fgetc(tmpOUT)) != EOF){
	    fprintf(outfile, "%c",c);
	  }
	  fclose(tmpOUT);
	  tmpOUT = NULL;
	  unlink(tmpOUTfile);
	}
      }
      if(tmpTRAJ != NULL){
	fclose(tmpTRAJ);
	if(!(tmpTRAJ = fopen(tmpTRAJfile,"r"))){
	  fprintf(errfile,"cannot read %s after simulation\n",tmpTRAJfile);
	  abort();
	}
	else{
	  char c;
	  while((c=fgetc(tmpTRAJ)) != EOF){
	    fprintf(trajFile, "%c",c);
	  }
	  fclose(tmpTRAJ);
	  fflush(trajFile);
	  tmpTRAJ = NULL;
	  unlink(tmpTRAJfile);
	}
      }
    }
    
    if(tmpTRAJ != NULL){
      fclose(tmpTRAJ);
      unlink(tmpTRAJfile);
      tmpTRAJ = NULL;
    }
    if(tmpOUT != NULL){
      fclose(tmpOUT);
      unlink(tmpOUTfile);
      tmpOUT = NULL;
    }
  }
  freeTrash(trash->Rtree);
  free(trash);

  free(dead);

  if(outfile != stdout)  fclose(outfile);
  if(errfile != stderr)  fclose(errfile);
  if(freqFile != NULL)   fclose(freqFile);
  if(ancFile != NULL)    fclose(ancFile);
  if(trajFile != NULL){
    fclose(trajFile);
    trajFile = NULL;
  }
  for(r=0; r<gpars.R; r++)  free(ancestralSEQ[r]);
  free(ancestralSEQ);
  freeStorage(storage);
  free_dmatrix(Q,0,63,0,3);
  free(ST);
  free(fitQuant);
/*   fclose(TMP); */
  return(0);
}

/* ------------------------------------------------------------------------- */

void setDefaultPPARS(const int pop)
{
  int i, j;

  ppars[pop].ALIVE = 0;
  ppars[pop].N = 500;
  ppars[pop].Nt = 500.0;
  ppars[pop].pFEMALES = 0.5;
  ppars[pop].MALES = (long)(ppars[pop].N*ppars[pop].pFEMALES+.5);
  ppars[pop].pMaleMig = 1-ppars[pop].pFEMALES;
  ppars[pop].pMaleRec = 0.5;
  ppars[pop].popAlpha = 0.0;
  ppars[pop].K = 0;
  ppars[pop].tauAlpha = -1;
  ppars[pop].GenEffect = 1;
  ppars[pop].SS = 6;
  ppars[pop].THETA = 0.001;
  ppars[pop].INSRATE = 0.0;
  ppars[pop].DELRATE = 0.0;
  ppars[pop].INDELlength = 1.0;
  ppars[pop].longINSRATE = 0.0;
  ppars[pop].longDELRATE = 0.0;
  ppars[pop].longINDELlength = 1.0;
  ppars[pop].INVRATE = 0.0;
  ppars[pop].INVlength = 1.0;
  ppars[pop].RHO = 0.0;
  ppars[pop].fGC = -1.0;
  ppars[pop].GCtract = 0;
  ppars[pop].BGC = 0.5;
  ppars[pop].SELF = 0.0;
  ppars[pop].neutpop = 0;
  assert(ppars[pop].selDistType = malloc(sizeof(*ppars[pop].selDistType)));
  ppars[pop].selDistType[0] = 0;
  ppars[pop].propSelLoci = 0.0;
  assert(ppars[pop].GAMMA = malloc(sizeof(*ppars[pop].GAMMA)));
  ppars[pop].GAMMA[0] = 0;
  assert(ppars[pop].ProbPos = malloc(sizeof(*ppars[pop].ProbPos)));
  ppars[pop].ProbPos[0] = 0.0;
  assert(ppars[pop].ProbDel = malloc(sizeof(*ppars[pop].ProbDel)));
  ppars[pop].ProbDel[0] = 0.0;
  assert(ppars[pop].ProbNeut = malloc(sizeof(*ppars[pop].ProbNeut)));
  ppars[pop].ProbNeut[0] = 1.0;
  assert(ppars[pop].ProbNegGamma = malloc(sizeof(*ppars[pop].ProbNegGamma)));
  ppars[pop].ProbNegGamma[0] = 0.0;
  assert(ppars[pop].alphaN = malloc(sizeof(*ppars[pop].alphaN)));
  ppars[pop].alphaN[0] = 0.17;
  assert(ppars[pop].lambdaN = malloc(sizeof(*ppars[pop].lambdaN)));
  ppars[pop].lambdaN[0] = 1/7800;
  assert(ppars[pop].alphaP = malloc(sizeof(*ppars[pop].alphaP)));
  ppars[pop].alphaP[0] = 25.0;
  assert(ppars[pop].lambdaP = malloc(sizeof(*ppars[pop].lambdaP)));
  ppars[pop].lambdaP[0] = 1.0;
  assert(ppars[pop].normMean = malloc(sizeof(*ppars[pop].normMean)));
  ppars[pop].normMean[0] = 0.0;
  assert(ppars[pop].normVar = malloc(sizeof(*ppars[pop].normVar)));
  ppars[pop].normVar[0] = 1.0;
  assert(ppars[pop].f0 = malloc(sizeof(*ppars[pop].f0)));
  ppars[pop].f0[0] = 1.0;
  ppars[pop].RateClassSites = 1;
  ppars[pop].RateParamSites = 0.17;
  ppars[pop].RateClassLoci = 1;
  ppars[pop].RateParamLoci = malloc(sizeof(double));
  ppars[pop].RateParamLoci[0] = 0.17;
  ppars[pop].KAPPA = 1.0;
  ppars[pop].KAPPACpG = 1.0;
  ppars[pop].PSI = 0.9;
  for(i=0; i<4; i++){
    ppars[pop].baseFreq[i] = 0.25;
    for(j=0; j<4; j++){
      ppars[pop].transProb[i][j] = (i==j ? 0 : 0.25/0.75);
      if(j>0)
	ppars[pop].transProb[i][j] += ppars[pop].transProb[i][j-1];
    }
  }
}

void setDefaultGPARS()
{
  long i, j, r;

  if(trash == NULL){
    assert(trash = malloc(sizeof(*trash)));
    trash->event = NULL;
    trash->Rtree = NULL;
    trash->Ltree = NULL;
    trash->Parent = NULL;
  }

  gpars.NPOPDEF = 1;
  gpars.INFSITES = 0;
  gpars.P = 2;
  gpars.tetraType = 0;
  gpars.R = 1;
  gpars.EQL = '1';
  assert(gpars.L = malloc(sizeof(*gpars.L)));
  gpars.L[0] = 5000;

  assert(gpars.LINK = malloc(sizeof(*gpars.LINK)));
  gpars.LINK[0] = 0.0;
  gpars.recMap = NULL;
  gpars.recMapPos = NULL;
  gpars.LinkDistType = 'p';
  assert(gpars.ANNOTATE = malloc(sizeof(*gpars.ANNOTATE)));
  gpars.ANNOTATE[0] = 'C';
  assert(gpars.SEX = malloc(sizeof(*gpars.SEX)));
  gpars.SEX[0] = '0';
  gpars.substMod = 0;
  assert(gpars.mig_mat = malloc(gpars.NPOP*sizeof(*gpars.mig_mat)));
  gpars.BURN = 5;
  gpars.BURN2 = 2;
  gpars.PRINTSEQ = 1;
  gpars.KEEPFIXED = 0;
  gpars.FIXSTOP = 1;
  gpars.ADDITIVE = 0;
  gpars.FITTEST = 0;
  gpars.TRACKANC = 0;
  gpars.TRACKTRAJ = 0;
  gpars.trajMinFreq = 0;
  gpars.trajMaxFreq = 1;
  gpars.trajPop = -1;
  gpars.trajLoc = -1;
  gpars.trajSite = -1;
  if(INITIALSEED == 0)
    gpars.trajMaxReps = 10000;
  gpars.SKIPSITES = NULL; /* only allocate if using --mutation */
  gpars.PRINTGEN = 0;

  assert(popn = malloc(gpars.NPOP*sizeof(*popn)));
  for(i=0; i<gpars.NPOP; i++){
    assert(gpars.mig_mat[i] = malloc(gpars.NPOP*sizeof(**gpars.mig_mat)));
    for(j=0; j<gpars.NPOP; j++)  gpars.mig_mat[i][j] = 0.0;
  }

  if(mutAr == NULL){
    /* need to create mutation array for each locus */
    assert(mutAr = malloc(gpars.R*sizeof(*mutAr)));
    for(r=0; r<gpars.R; r++){
      assert(mutAr[r] = malloc(sizeof(**mutAr)));
      /* create first mutation as default */
      mutAr[r]->numMuts = 1;
      mutAr[r]->mutIndex = 0;
      assert(mutAr[r]->muts = malloc(sizeof(struct history*)));
      mutAr[r]->muts[0] = popTrash();
      assert(mutAr[r]->muts[0]->event=malloc(sizeof(struct event)));
      mutAr[r]->muts[0]->event->site = 0;
      assert(mutAr[r]->muts[0]->event->genFix =malloc(gpars.NPOP*sizeof(long)));
      assert(mutAr[r]->muts[0]->event->genDead=malloc(gpars.NPOP*sizeof(long)));
      mutAr[r]->muts[0]->event->nSites = 0;
      mutAr[r]->muts[0]->event->nucs = NULL;
      assert(mutAr[r]->muts[0]->event->fixed = malloc(gpars.NPOP*sizeof(char)));
      mutAr[r]->muts[0]->event->free = '1'; /* free by default */
      mutAr[r]->muts[0]->event->index = 0;
      mutAr[r]->muts[0]->event->checkGen = -1;
      mutAr[r]->muts[0]->event->checkStep = -1;
      mutAr[r]->muts[0]->event->checkPop = -1;
      assert(mutAr[r]->muts[0]->event->nextHap=malloc(gpars.NPOP*sizeof(long)));
      assert(mutAr[r]->muts[0]->event->numCarriers = 
	     malloc(gpars.NPOP*sizeof(long)));
      assert(mutAr[r]->muts[0]->event->maxHaps=malloc(gpars.NPOP*sizeof(long)));
      assert(mutAr[r]->muts[0]->event->hapFreq=
	     malloc(gpars.NPOP*sizeof(long**)));
      for(j=0; j<gpars.NPOP; j++){
	mutAr[r]->muts[0]->event->nextHap[j] = 0;
	mutAr[r]->muts[0]->event->genFix[j] = 0;
	mutAr[r]->muts[0]->event->genDead[j] = -LONG_MAX;
	mutAr[r]->muts[0]->event->fixed[j] = '0'; /* not fixed by default */
	mutAr[r]->muts[0]->event->numCarriers[j] = 0;
	mutAr[r]->muts[0]->event->maxHaps[j] = 1;
	assert(mutAr[r]->muts[0]->event->hapFreq[j] = malloc(sizeof(long*)));
	mutAr[r]->muts[0]->event->hapFreq[j][0] = NULL; /* no carriers yet */
      }
    }
  }
  if(dead == NULL){
    assert(dead = malloc(sizeof(*dead)));
    dead->chr = 0;
    dead->loc = 0;
    dead->next = NULL;
    dead->prev = NULL;
  }
}

/* ------------------------------------------------------------------------- */

void Initialize(const double *ST)
{
  int k, SEQBURN=100;
  long i, r, acc, reg;
  float rn;
  char *nucs;
  assert(nucs = malloc(4*sizeof(char)));
  nucs[3] = '\0';

  if(popn[0].conSeq == NULL){
    assert(popn[0].conSeq = malloc(gpars.R*sizeof(char*)));
    for(r=0; r<gpars.R; r++)
      popn[0].conSeq[r] = NULL;
  }

#ifdef VERBOSE_DEBUG
  printf("generating random sequence\n");
#endif
  /* GENERATE INITIAL SEQUENCE */
  for(reg=0; reg<gpars.R; reg++){
    if(popn[0].conSeq[reg] != NULL){
      if(gpars.ANNOTATE[reg] == 'C' && gpars.FIXSTOP == 1){
	for(i=0; i<gpars.L[reg]-2; i+=3){
	  if((k = 16*(popn[0].conSeq[reg][i]-'0') +
	      4*(popn[0].conSeq[reg][i+1]-'0') + 
	      popn[0].conSeq[reg][i+2]-'0') == 39 || k == 45 || k == 47){
	    /* fix in frame stop codons */
	    rn=ran1(&gpars.seed);
	    if(rn<=ppars[0].baseFreq[0])
	      popn[0].conSeq[reg][i+(int)(rn*3)]='0';
	    else if(rn<=ppars[0].baseFreq[0]+ppars[0].baseFreq[1])
	      popn[0].conSeq[reg][i+(int)(rn*3)]='1';
	    else if(rn<=(ppars[0].baseFreq[0]+ppars[0].baseFreq[1]+
			 ppars[0].baseFreq[2]))
	      popn[0].conSeq[reg][i+(int)(rn*3)]='2';
	    else popn[0].conSeq[reg][i+(int)(rn*3)]='3';
	    i -= 3;
	  }
	}
	continue;
      }
      continue;
    }
    assert(popn[0].conSeq[reg] = malloc((gpars.L[reg]+1)*sizeof(char)));
    if(reg > 0 && gpars.L[reg] <= gpars.L[reg-1] && gpars.substMod != 0){
#ifdef WINDOWS
      strncpy_s(popn[0].conSeq[reg], gpars.L[reg],
		popn[0].conSeq[reg-1], gpars.L[reg]);
#else
      strncpy(popn[0].conSeq[reg], popn[0].conSeq[reg-1], gpars.L[reg]);
#endif
      SEQBURN = 100;
    }
    else{
      SEQBURN=100;
      for(i=0;i<gpars.L[reg];i++){
	rn=ran1(&gpars.seed);
	if(rn<=ppars[0].baseFreq[0])
	  popn[0].conSeq[reg][i]='0';
	else if(rn<=ppars[0].baseFreq[0]+ppars[0].baseFreq[1])
	  popn[0].conSeq[reg][i]='1';
	else if(rn<=(ppars[0].baseFreq[0]+ppars[0].baseFreq[1]+
		     ppars[0].baseFreq[2]))
	  popn[0].conSeq[reg][i]='2';
	else popn[0].conSeq[reg][i]='3';
	
	if(gpars.ANNOTATE[reg]=='C' && i%3 == 2 && i>0)
	  if((k = 16*(popn[0].conSeq[reg][i-2]-'0') +
	      4*(popn[0].conSeq[reg][i-1]-'0') + 
	      popn[0].conSeq[reg][i]-'0') == 39 || k == 45 || k == 47)
	    i-=3;  /* stop codons */
      }
    }
    popn[0].conSeq[reg][gpars.L[reg]] = '\0';
    if(gpars.substMod > 0){ /* more complicated than random sequence */
      for(acc=0; acc<SEQBURN*gpars.L[reg]; acc++){
	if(acc%gpars.L[reg] > 0)
	  nucs[0] = popn[0].conSeq[reg][acc%gpars.L[reg]-1];
	else /* loop back to end of sequence */
	  nucs[0] = popn[0].conSeq[reg][gpars.L[reg]-1]; 
	nucs[1] = popn[0].conSeq[reg][acc%gpars.L[reg]];
	if(acc%gpars.L[reg] < gpars.L[reg]-1)
	  nucs[2] = popn[0].conSeq[reg][acc%gpars.L[reg]+1];
	else /* loop to beginning of sequence */
	  nucs[2] = popn[0].conSeq[reg][0];
	i = 0;
	if((nucs[0]=='0' && nucs[1]=='1') || (nucs[1]=='0' && nucs[2]=='1'))
	  i = 1;
	popn[0].conSeq[reg][acc%gpars.L[reg]] = MutantNuc(nucs, i, 0);
	if(gpars.ANNOTATE[reg]=='C' && (acc%gpars.L[reg])%3 == 2){
	  if((k = (16*(popn[0].conSeq[reg][acc%gpars.L[reg]-2]-'0') + 
		   4*(popn[0].conSeq[reg][acc%gpars.L[reg]-1]-'0') + 
		   popn[0].conSeq[reg][acc%gpars.L[reg]]-'0')) == 39 ||
	     k == 45 || k == 47) /* stop codons */
	    acc-=3;
	}
      }
    }
  }
  free(nucs);
  nucs = NULL;
}

/* ------------------------------------------------------------------------- */

int AA(const char *c) /* takes codon, returns AA */
{
  if(c[0] == '4' || c[1] == '4' || c[2] == '4')
    return 21;  /* nothing! */
  if(c[0]=='0'){
    if(c[1]=='0')
      return 1;      /* PROLINE */
    else if(c[1]=='1')
      return 2;      /* ARGININE */
    else if(c[1]=='2')
      return 3;      /* LEUCINE */
    else if(c[1]=='3'){
      if(c[2]=='0' || c[2]=='2')
	return 4;  /* HISTIDINE */
      else if(c[2]=='1' || c[2]=='3')
	return 5;  /* GLUTAMINE */
      else{
	fprintf(errfile,"%c%c%c not an AA (%ld)\n",c[0],c[1],c[2],INITIALSEED);
	exitNOW("\n");
	return(-1);
      }
    }
    else{
      fprintf(errfile,"%c%c%c not an AA (%ld)\n",c[0],c[1],c[2],INITIALSEED);
      exitNOW("\n");
      return(-1);
    }
  }
  else if(c[0]=='1'){
    if(c[1]=='0')
      return 6;      /* ALANINE */
    else if(c[1]=='1')
      return 7;      /* GLYCINE */
    else if(c[1]=='2')
      return 8;      /* VALINE */
    else if(c[1]=='3'){
      if(c[2]=='0' || c[2]=='2')
	return 9;  /* ASPARTIC ACID */
      else
	return 10; /* GLUTAMIC ACID */
    }
    else{
      fprintf(errfile,"%c%c%c not an AA (%ld)\n",c[0],c[1],c[2],INITIALSEED);
      exitNOW("\n");
      return(-1);
    }
  }
  else if(c[0]=='2'){
    if(c[1]=='0')
      return 11;     /* SERINE */
    else if(c[1]=='1'){
      if(c[2]=='0' || c[2]=='2')
	return 12; /* CYSTEINE */
      else if(c[2]=='1')
	return 13; /* TRYPTOPHAN */
      else if(c[2]=='3')
	return 0;  /* STOP */
      else{
	fprintf(errfile,"%c%c%c not an AA (%ld)\n",c[0],c[1],c[2],INITIALSEED);
	exitNOW("\n");
	return(-1);
      }
    }
    else if(c[1]=='2'){
      if(c[2]=='0' || c[2]=='2')
	return 14;  /* PHENYLALANINE */
      else
	return 3;   /* LEUCINE */
    }
    else if(c[2]=='0' || c[2]=='2')
      return 15;      /* TYROSINE */
    else
      return 0;       /* STOP */
  }
  else if(c[1]=='0')
    return 16;          /* THREONINE */
  else if(c[1]=='1'){
    if(c[2]=='0' || c[2]=='2')
      return 11;      /* SERINE */
    else
      return 2;       /* ARGININE */
  }
  else if(c[1]=='2'){
    if(c[2]=='1')
      return 17;      /* METHIONINE (START) */
    else
      return 18;      /* ISOLEUCINE */
  }
  else if(c[1]=='3'){
    if(c[2]=='0' || c[2]=='2')
      return 19;      /* ASPARAGINE */
    else
      return 20;      /* LYSINE */
  }
  else{
    printf("error in AA (%ld): %c%c%c\n",INITIALSEED,c[0],c[1],c[2]);
    abort();
    return(-1);
  }
  
  /*************************************************************************
   * 000 = CCC P pro | 100 = GCC A ala | 200 = TCC S ser | 300 = ACC T thr *
   * 001 = CCG P pro | 101 = GCG A ala | 201 = TCG S ser | 301 = ACG T thr *
   * 002 = CCT P pro | 102 = GCT A ala | 202 = TCT S ser | 302 = ACT T thr *
   * 003 = CCA P pro | 103 = GCA A ala | 203 = TCA S ser | 303 = ACA T thr *
   *************************************************************************
   * 010 = CGC R arg | 101 = GGC G gly | 210 = TGC C cys | 310 = AGC S ser *
   * 011 = CGG R arg | 111 = GGG G gly | 211 = TGG W try | 311 = AGG R arg *
   * 012 = CGT R arg | 112 = GGT G gly | 212 = TGT C cys | 312 = AGT S ser *
   * 013 = CGA R arg | 113 = GGA G gly | 213 = TGA * OPA | 313 = AGA R arg *
   *************************************************************************
   * 020 = CTC L leu | 120 = GTC V val | 220 = TTC F phe | 320 = ATC I iso *
   * 021 = CTG L leu | 121 = GTG V val | 221 = TTG L leu | 321 = ATG M met *
   * 022 = CTT L leu | 122 = GTT V val | 222 = TTT F phe | 322 = ATT I iso *
   * 023 = CTA L leu | 123 = GTA V val | 223 = TTA L leu | 323 = ATA I iso *
   *************************************************************************
   * 030 = CAC H his | 130 = GAC D asp | 230 = TAC Y tyr | 330 = AAC N asn *
   * 031 = CAG Q gln | 131 = GAG E glu | 231 = TAG * AMB | 331 = AAG K lys *
   * 032 = CAT H his | 132 = GAT D asp | 232 = TAT Y tyr | 332 = AAT N asn *
   * 033 = CAA Q gln | 133 = GAA E glu | 233 = TAA * OCH | 333 = AAA K lys *
   *************************************************************************/
}

/* ------------------------------------------------------------------------- */

void bign(double **mutClass_site, double *mutClass_locus, int *selClass,
	  const double *fitQuant, long *numMutSeg)
{
  int BURNING=1, DOMEST=0, genEff, pop=0, ITERSUCCESS=0;
  long i, j, k, r, offset;
  double *part, *relFit;      /* Partition array for each population*/
  long *Pars;           /* Parents to found next generation */
  long popSS, count;
  int MIG=0;
  int **DomestAllele, *survive=NULL;
  unsigned long t, nextDevent, CHANGED=0;
  long PNanc;
  double *rates_site, *rates_loci;
  double sum;
  double ***popFreqs=NULL, **mutLocs=NULL;
  long **orderMutLocs=NULL;
  struct EvolEvents *tmp;
  struct mutStorage *sampStore = NULL;
  
  /* initialize loads of things */
  PNanc = gpars.P*ppars[0].N; /* ancestral population size */
  popn[0].maxchr = 0;
  popn[0].gen = (long)(-gpars.BURN*PNanc);
  assert(popn[0].nextchr = malloc(gpars.R*sizeof(long)));
  for(i=0; i<gpars.R; i++)
    popn[0].nextchr[i] = 0;
  assert(popn[0].parentLoc = malloc(PNanc*sizeof(long*)));
  for(i=0; i<PNanc; i++)
    assert(popn[0].parentLoc[i] = malloc(gpars.R*sizeof(long)));

  assert(popn[0].numCopies = malloc(sizeof(long*))); /* parentLoc */
  assert(popn[0].chrGen = malloc(sizeof(long*))); 
  assert(popn[0].numCopies[0] = malloc(gpars.R*sizeof(long))); /* locus */
  assert(popn[0].chrGen[0] = malloc(gpars.R*sizeof(long))); /* locus */

  assert(popn[0].extremeMuts = malloc(sizeof(long**))); /* chr */
  assert(popn[0].extremeMuts[0] = malloc(gpars.R*sizeof(long*)));/*locus*/

  if(gpars.substMod == 5){
    assert(popn[0].hpSUM = malloc(sizeof(double**)));
    assert(popn[0].hpSUM[0] = malloc(gpars.R*sizeof(double*))); /* locus */
    for(r=0; r<gpars.R; r++){
      assert(popn[0].hpSUM[0][r] = malloc(gpars.L[r]*sizeof(double)));/* site */
      popn[0].hpSUM[0][r][0] = 0.0;
    }
  }

  assert(selfed = calloc(ppars[0].N, sizeof(*selfed)));
  assert(selfing = calloc(ppars[0].N, sizeof(*selfing)));
  assert(parStore = malloc(gpars.P*ppars[0].N*sizeof(long*)));
  for(i=0; i<gpars.P*ppars[0].N; i++){
    assert(parStore[i] = malloc(gpars.R*sizeof(long)));
    if(i<ppars[0].N){
      selfed[i] = 0;
      selfing[i] = 0;
    }
  }
  assert(popn[0].indfit = malloc(sizeof(double*)));
  assert(popn[0].indfit[0] = malloc(gpars.R*sizeof(double)));

  assert(popn[0].fixed = malloc(gpars.R*sizeof(struct history*)));

  assert(popn[0].polySites=malloc(gpars.R*sizeof(int*)));
  part = dvector(0,ppars[0].N-1); /* fitness partition */
  relFit = dvector(0,ppars[0].N-1); /* relative fitnesses */
  assert(DomestAllele = malloc(gpars.R*sizeof(*DomestAllele)));
  assert(Pars = malloc(PNanc*sizeof(*Pars)));
  
  if(gpars.TRACKANC){
    assert(popn[0].ancestry = malloc(sizeof(int**))); /* chrms */
    assert(popn[0].ancestry[0] = malloc(gpars.R*sizeof(int*))); /* loci */
  }
  
  for(r=0;r<gpars.R;r++){
    assert(popn[0].polySites[r] = malloc(gpars.L[r]*sizeof(int)));
    popn[0].fixed[r] = popTrash();
    assert(popn[0].extremeMuts[0][r] = malloc(2*sizeof(long)));
    popn[0].extremeMuts[0][r][0] = gpars.L[r]; /* min at end of seq */
    popn[0].extremeMuts[0][r][1] = 0; /* max at beg of seq */
    popn[0].numCopies[0][r] = PNanc;
    popn[0].chrGen[0][r] = 0;
    popn[0].indfit[0][r] = (gpars.ADDITIVE ? 0.0 : 1.0);
    selClass[r] = 0;
    assert(DomestAllele[r] = malloc(gpars.L[r]*sizeof(**DomestAllele)));
    if(gpars.TRACKANC)
      assert(popn[0].ancestry[0][r] = malloc(gpars.L[r]*sizeof(int)));
    for(i=0; i<PNanc; i++)  popn[0].parentLoc[i][r] = 0;
    for(j=0; j<gpars.L[r]; j++){
      popn[0].polySites[r][j] = 0;
      DomestAllele[r][j] = 0;
      if(gpars.TRACKANC)
	popn[0].ancestry[0][r][j] = 0;
    }
  }
  /* done with most of initialization... */
  
  if(ancestralSEQ != NULL){ /* piggy-back off previous burn-in */
    regeneratePOP();
    renewConSeq();
    if(gpars.TRACKANC){
      for(j=0; j<=popn[0].maxchr; j++){
	for(r=0; r<gpars.R; r++){
	  for(i=0; i<gpars.L[r]; i++){
	    popn[0].ancestry[j][r][i] = 0;
	  }
	}
      }
    }
#ifdef UNIT_TEST
    checkForErrors(0, "After regeneratePOP");
#endif
  }
  
  /* only necessary if using context-dependent mutation model */
  if(gpars.substMod == 5){
    for(j=0; j<=popn[0].maxchr; j++){
      for(r=0; r<gpars.R; r++){
	if(popn[0].numCopies[j][r] > 0){
	  getNewHitProbsReg(&popn[0].hpSUM[j][r], &popn[0].BigHead[j][r], r, 0,
			    0, gpars.L[r]-1);
#ifdef UNIT_TEST
	  checkHitProbs(popn[0].hpSUM[j][r], "error at beginning", 0, j, r);
	  if(popn[0].hpSUM[j][r][gpars.L[r]-1] < DBL_EPSILON){
	    fprintf(errfile,"hpSUM[%d][%ld][%ld][%lu] == 0!! (%ld)\n",
		    0, j, r, gpars.L[r]-1,INITIALSEED);
	    abort();
	  }
#endif
	}
      }
    }
  }
  
  /* determine rate classes for each site */
  rates_site = dvector(0,ppars[0].RateClassSites-1);
  if(ppars[0].RateClassLoci == 0)
    rates_loci = dvector(0,gpars.R-1);
  else
    rates_loci = dvector(0,ppars[0].RateClassLoci-1);

  if(ppars[0].RateClassSites == 1){ /* uniform */
    for(r=0;r<gpars.R;r++)
      for(i=0;i<gpars.L[r];i++)
	mutClass_site[r][i]=(i+1)/((double)gpars.L[r]);
  }
  else{
    for(i=0;i<ppars[0].RateClassSites;i++) 
      rates_site[i]=  /* gamma, mean = 1 */
	rgamma(ppars[0].RateParamSites,ppars[0].RateParamSites,&gpars.seed);
    for(r=0;r<gpars.R;r++){
      sum = 0.0;
      for(i=0;i<gpars.L[r];i++){
	mutClass_site[r][i]=
	  rates_site[(int)(ppars[0].RateClassSites*ran1(&gpars.seed))];
	sum += mutClass_site[r][i];
      }
      mutClass_site[r][0]=mutClass_site[r][0]/sum;
      for(i=1;i<gpars.L[r];i++)
	mutClass_site[r][i]=mutClass_site[r][i-1]+mutClass_site[r][i]/sum;
    }
  }
  
  /* determine rate classes for each locus */
  if(ppars[0].RateClassLoci == 0){ /* defined at command-line */
    sum = 0.0;
    for(r=0; r<gpars.R; r++)
      sum += ppars[0].RateParamLoci[r];
    mutClass_locus[0] = ppars[0].RateParamLoci[0]/sum;
    for(r=1; r<gpars.R; r++)
      mutClass_locus[r] = mutClass_locus[r-1] + ppars[0].RateParamLoci[r]/sum;
    if(sum < FLT_EPSILON){
      fprintf(errfile, "must have some non-zero rate classes\n");
      abort();
    }
  }
  else if(ppars[0].RateClassLoci == 1){
    sum = 0.0;
    for(r=0; r<gpars.R; r++)
      sum += gpars.L[r]+.0;
    mutClass_locus[0] = gpars.L[0]/sum;
    for(r=1; r<gpars.R; r++)
      mutClass_locus[r] = mutClass_locus[r-1]+gpars.L[r]/sum;
  }
  else{
    for(r=0; r<ppars[0].RateClassLoci; r++)  /* gamma, mean = 1 */
      rates_loci[r]=rgamma(ppars[0].RateParamLoci[0],ppars[0].RateParamLoci[0],
			   &gpars.seed);
    sum = 0.0;
    for(r=0; r<gpars.R; r++){
      mutClass_locus[r]=
	rates_loci[(int)(ppars[0].RateClassLoci*ran1(&gpars.seed))];
      sum += mutClass_locus[r];
    }
    mutClass_locus[0] /= sum;
    for(r=1; r<gpars.R; r++)
      mutClass_locus[r]=mutClass_locus[r-1]+mutClass_locus[r]/sum;
  }
  
  /* randomly choose selected loci */
  if(fabs(ppars[0].propSelLoci - 1.0) < DBL_EPSILON)
    for(r=0; r<gpars.R; r++)
      selClass[r] = 1;
  else{
    count=(long)(gpars.R*ppars[0].propSelLoci+.5); /* +.5 to round up */
    for(r=0; r<count; r++){
      i=(long)(gpars.R*ran1(&gpars.seed));
      if(selClass[i] != 1)
	selClass[i] = 1;
      else
	r--;
    }
  }
  
  /* determine if migration is necessary */
  for(i=0; i<gpars.NPOP; i++){
    ppars[i].P0 = ppars[0].N; /* ancestral population size */
    if(MIG == 0)
      for(j=0; j<gpars.NPOP; j++)
	if(gpars.mig_mat[i][j] > DBL_EPSILON){
	  MIG = 1;
	  break;
	}
  }
  /* I GIVE YOU LIFE; commence generations */
#ifdef UNIT_TEST
  fprintf(errfile,"COMMENCE BURNIN\n");
  fflush(errfile);
#endif
  t = 0;
  BURNING = 1;
  ppars[0].ALIVE = 1;
  nextDevent = gpars.BURN*PNanc;
  
  while(1){
    if(BURNING == 1 && t == nextDevent){
#ifdef UNIT_TEST
      fprintf(errfile,"BURNIN COMPLETED\n");
      fflush(errfile);
#endif
      BURNING = 0;
      t = 0;
      popn[0].gen = 0;
      if(devents != NULL)
	nextDevent = devents->tau*PNanc;
      if(gpars.KEEPFIXED == 0){
	for(r=0; r<gpars.R; r++){ /* not keeping substitutions during BURNIN */
	  if(popn[0].fixed[r]->Rtree != NULL){
	    freeFixedHistory(popn[0].fixed[r]->Rtree, 0);
	    popn[0].fixed[r]->Rtree = NULL;
	  }
	}
      }
      
      /* print ancestral sequence */
      if(gpars.PRINTSEQ && gpars.trajMinFreq<FLT_EPSILON &&
	 gpars.trajMaxFreq>1-FLT_EPSILON){
	if(tmpOUT == NULL){
	  for(r=0; r<gpars.R; r++){
	    fprintf(outfile,">locus_%ld\n",r);
	    for(i=0; i<gpars.L[r]; i++)
	      fprintf(outfile,"%c",ToCGTA(popn[0].conSeq[r][i]));
	    fprintf(outfile,"\n");
	  }
	  fflush(outfile);
	}
	else{
	  for(r=0; r<gpars.R; r++){
	    fprintf(tmpOUT,">locus_%ld\n",r);
	    for(i=0; i<gpars.L[r]; i++)
	      fprintf(tmpOUT,"%c",ToCGTA(popn[0].conSeq[r][i]));
	    fprintf(tmpOUT,"\n");
	  }
	  fflush(tmpOUT);
	}
      }
      
      /* store population at end of burnin to be used next iteration */
#ifdef UNIT_TEST
      (*numMutSeg) = 0;
#endif
      offset = gpars.P*ppars[0].N*gpars.BURN;
      freeStorage(storage);
      storage = NULL;
      for(j=0; j<gpars.P*ppars[0].N; j++){
	for(r=0; r<gpars.R; r++){
#ifdef UNIT_TEST
	  countSeg(popn[0].BigHead[popn[0].parentLoc[j][r]][r]->Rtree,
		   numMutSeg);
#endif
	  storeInfo(&storage,popn[0].BigHead[popn[0].parentLoc[j][r]][r]->Rtree,
		    0, j, r, offset, 0);
	}
      }

#ifdef UNIT_TEST
      {
	long foo = 0;
	countSegStorage(storage, &foo);
	if(foo != (*numMutSeg)){
	  fprintf(errfile,"sfs_code error(%ld):  did not store mutations \
correctly! Expecting %ld mutations, stored %ld.\n", INITIALSEED,*numMutSeg,foo);
	  abort();
	}
      }
#endif
      if(ancestralSEQ == NULL){
	assert(ancestralSEQ = malloc(gpars.R*sizeof(*ancestralSEQ)));
	for(r=0; r<gpars.R; r++){
	  assert(ancestralSEQ[r] = 
		 malloc((gpars.L[r]+1)*sizeof(**ancestralSEQ)));
	}
      }
      for(r=0; r<gpars.R; r++){
#ifdef WINDOWS
	strcpy_s(ancestralSEQ[r], gpars.L[r]+1, popn[0].conSeq[r]);
#else
	strcpy(ancestralSEQ[r], popn[0].conSeq[r]);
#endif
      }
      continue;
    }
    else if(BURNING == 0 && devents == NULL)  break; /* done! */
    else if(devents != NULL && t >= nextDevent){ /* evolutionary event! */
#ifdef UNIT_TEST
      printf("eventType=%d; parIndex=%d\n",devents->eventType,
	     devents->parIndex);
#endif
      if(devents->eventType == 0){ /* split population into 2 */
#ifdef UNIT_TEST
	long foo = 0;
	fprintf(errfile,"eventType 0: t=%ld, pop%d->%d\n",t,devents->popi,
		devents->popj);
	fflush(errfile);
	for(j=0; j<gpars.P*ppars[devents->popi].N; j++)
	  for(r=0; r<gpars.R; r++)
	    countSeg(popn[devents->popi].BigHead
		     [popn[devents->popi].parentLoc[j][r]][r]->Rtree,
		     &foo);
	checkForErrors(devents->popi, "Before split");
#endif
	Split(devents->popi, devents->popj);
	ppars[devents->popj].ALIVE = 1;
	if(gpars.TRACKANC){
	  for(j=0; j<=popn[devents->popj].maxchr; j++){
	    for(r=0; r<gpars.R; r++){
	      for(i=0; i<gpars.L[r]; i++){
		popn[devents->popj].ancestry[j][r][i] = devents->popj;
	      }
	    }
	  }
	}
	
#ifdef UNIT_TEST
	{
	  long foo2 = 0;
	  for(j=0; j<gpars.P*ppars[devents->popj].N; j++)
	    for(r=0; r<gpars.R; r++)
	      countSeg(popn[devents->popj].BigHead[popn[devents->popj].
						   parentLoc[j][r]][r]->Rtree,
		       &foo2);
	  if(foo != foo2){
	    fprintf(errfile,"error in split!  not all mutations copied.\n");
	    fprintf(errfile,"pop %d has %ld; pop %d has %ld\n",devents->popi,
		    foo, devents->popj, foo2);
	    abort();
	  }
	  checkForErrors(devents->popj, "After split");
	}
#endif
      }
      else if(devents->eventType == 1){ /* domestication event */
#ifdef UNIT_TEST
	fprintf(errfile,"eventType 1 from pop %d->%d with allele freq %f and \
N = %1.0f\n",devents->popi,devents->popj, devents->freq, devents->newP.Nt);
	fflush(errfile);
#endif
	
	Split(devents->popi, devents->popj);
	domesticate(DomestAllele, survive, devents->popj);
	ppars[devents->popj].ALIVE = 1;
	DOMEST = 1;
      }
      else if(devents->eventType == 7){ /* change global parameter */
	if(devents->parIndex == 0){ /* migration matrix */
#ifdef UNIT_TEST
	  fprintf(errfile,"changing migration parameters to %d from %d\n",
		  devents->popi, devents->popj);
	  fflush(errfile);
#endif
	  if(devents->popi == -1 && devents->popj == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      for(j=0; j<gpars.NPOP; j++)
		gpars.mig_mat[i][j] = devents->newG.mig_mat[i][j];
	  }
	  else{
	    gpars.mig_mat[devents->popi][devents->popj] =
	      devents->newG.mig_mat[devents->popi][devents->popj];
	  }
	  for(i=0; i<gpars.NPOP; i++)
	    free(devents->newG.mig_mat[i]);
	  free(devents->newG.mig_mat);
	  
	  MIG = 0; /* default: assume no migration */
	  for(i=0; i<gpars.NPOP; i++){
	    if(MIG == 0)
	      for(j=0; j<gpars.NPOP; j++)
		if(i != j)
		  if(gpars.mig_mat[i][j] > DBL_EPSILON){
		    MIG = 1; /* non-zero migration rate! */
		    break;
		  }
	  }
	}
	else if(devents->parIndex == 1){
	  gpars.TRACKTRAJ = 1;
	  gpars.trajPop = devents->popi;
	  gpars.trajLoc = devents->locus;
	  gpars.trajSite = devents->site;
#ifdef UNIT_TEST
	  printf("start tracking site %ld!\n",gpars.trajSite);
	  fflush(stdout);
#endif
	}
      }
      else if(devents->eventType == 8){
	admixture(devents, PNanc, selClass);
	freeDead(devents->popi, 0);
#ifdef UNIT_TEST
	printf("admixed!\n");
	fflush(stdout);
	checkForErrors(devents->popi,"after admixture\n");
#endif
	free(devents->ancPops);
	free(devents->maleFreqs);
	free(devents->femaleFreqs);
      }
      else if(devents->eventType == 9){
	int mutpop = 0, CpG=0;
	long mutind, mutindPL, mutloc, mutreg;
	char n[5], newnuc, oc[3], nc[3];
	struct event *event = NULL;
	assert(event = malloc(sizeof(struct event)));
	assert(event->genFix = malloc(gpars.NPOP*sizeof(long)));
	assert(event->genDead = malloc(gpars.NPOP*sizeof(long)));
	for(i=0; i<gpars.NPOP; i++){
	  event->genFix[i] = 0;
	  event->genDead[i] = -LONG_MAX;
	}
	if(devents->popi >= 0)
	  mutpop = devents->popi;
	else{ /* choose population randomly */
	  pop=0;
	  for(i=0; i<gpars.NPOP; i++){
	    if(ppars[i].ALIVE){
	      pop++;
	    }
	  }
	  j=(long)(ran1(&gpars.seed)*pop);
	  for(mutpop=0; mutpop<gpars.NPOP; mutpop++){
	    if(ppars[mutpop].ALIVE){
	      if(j==0){
		break;
	      }
	      else{
		j--;
	      }
	    }
	  }
	}
	mutreg = devents->locus;
	mutloc = devents->site;
	mutind=(long)(gpars.P*ppars[mutpop].N*ran1(&gpars.seed));
	createNewHistory(mutind, mutreg, mutpop, 1, 1);
	mutindPL = popn[mutpop].parentLoc[mutind][mutreg];
	if(mutloc > 1){
	  n[3] = retNucHistory(mutloc-2,&popn[mutpop].BigHead[mutindPL][mutreg],
			       mutreg, mutpop, 1);
	}
	else
	  n[3] = 'N';
	if(mutloc > 0)
	  n[0] = retNucHistory(mutloc-1,&popn[mutpop].BigHead[mutindPL][mutreg],
			       mutreg, mutpop, 1);
	else
	  n[0] = 'N';
	n[1] = retNucHistory(mutloc,&popn[mutpop].BigHead[mutindPL][mutreg],
			     mutreg, mutpop, 1);
	if(mutloc < gpars.L[mutreg]-1)
	  n[2] = retNucHistory(mutloc+1,&popn[mutpop].BigHead[mutindPL][mutreg],
			       mutreg, mutpop, 1);
	else
	  n[2] = 'N';
	if(mutloc<gpars.L[mutreg]-2)
	  n[4] = retNucHistory(mutloc+2,&popn[mutpop].BigHead[mutindPL][mutreg],
			       mutreg, mutpop, 1);
	else
	  n[4] = 'N';
	if((n[0]=='0' && n[1]=='1') || (n[1]=='0' && n[2]=='1'))
	  CpG=1;
	else
	  CpG=0;
	do{
	  newnuc=MutantNuc(n, CpG, mutpop);
	  if((mutloc%3)==0){
	    nc[0] = newnuc;
	    nc[1] = n[2];
	    nc[2] = n[4];
	    oc[0] = n[1];
	    oc[1] = nc[1];
	    oc[2] = nc[2];
	  }
	  else if((mutloc%3)==1){
	    nc[0] = n[0];
	    nc[1] = newnuc;
	    nc[2] = n[2];
	    oc[0] = nc[0];
	    oc[1] = n[1];
	    oc[2] = nc[2];
	  }
	  else{
	    nc[0] = n[3];
	    nc[1] = n[0];
	    nc[2] = newnuc;
	    oc[0] = nc[0];
	    oc[1] = nc[1];
	    oc[2] = n[1];
	  }
	}while(gpars.ANNOTATE[mutreg] == 'C' && AA(nc)==0);
	event->site = mutloc;
	event->gen = popn[mutpop].gen;
	event->genFix[mutpop] = 0;
	event->genDead[mutpop] = -LONG_MAX;
	event->ancNuc = n[1];
	event->derNuc = newnuc;
	event->nonsyn =
	  (char)(gpars.ANNOTATE[mutreg]=='C'?(AA(oc)!=AA(nc))+'0':'0');
	event->ancAA = (gpars.ANNOTATE[mutreg]=='C' ? AA(oc) : 21);
	event->derAA = (gpars.ANNOTATE[mutreg]=='C' ? AA(nc) : 21);
	event->fit = devents->gamma/PNanc;
	event->CpG = (char)(CpG+'0');
	event->fiveP = n[0];
	event->threeP= n[2];
	if(gpars.SEX[mutreg] == '0')
	  event->axy = 'A'; /* autosome */
	else{ /* sex chromosome */
	  if(mutind < gpars.P*ppars[pop].MALES || mutind%gpars.P < gpars.P/2)
	    event->axy = 'X';  /* X-linked */
	  else
	    event->axy = 'Y';  /* Y-linked */
	}
	addHistoryNode(event, &popn[mutpop].BigHead[mutindPL][mutreg]->Rtree,
		       NULL, mutpop, mutindPL, mutreg, 0);
	popn[mutpop].polySites[mutreg][mutloc]++;
	if(gpars.ADDITIVE)
	  popn[mutpop].indfit[mutindPL][mutreg] += event->fit;
	else
	  popn[mutpop].indfit[mutindPL][mutreg] *= 1.0+event->fit;
	if(fabs(event->fit) > 0)
	  ppars[mutpop].neutpop = 0;

	if(popn[mutpop].extremeMuts[mutindPL][mutreg][0] > mutloc)
	  popn[mutpop].extremeMuts[mutindPL][mutreg][0] = mutloc;
	if(popn[mutpop].extremeMuts[mutindPL][mutreg][1] < mutloc)
	  popn[mutpop].extremeMuts[mutindPL][mutreg][1] = mutloc;
	free(event->genFix);
	free(event->genDead);
	free(event);
#ifdef UNIT_TEST
	printf("mutation injected\n");
#endif
      }
      else if(devents->eventType == 3){ /* change pop size */
#ifdef UNIT_TEST
	fprintf(errfile,"eventType 3 to pop %d with strength %f\n",
		devents->popi,devents->nu);
	fflush(errfile);
#endif
	if(devents->popi == -1){ /* change all population sizes */
	  for(i=0; i<gpars.NPOP; i++){
	    ppars[i].Nt *= devents->nu;
	    ppars[i].THETA *= devents->nu;
	    ppars[i].INSRATE *= devents->nu;
	    ppars[i].DELRATE *= devents->nu;
	    if(ppars[i].RHO > DBL_EPSILON)
	      ppars[i].RHO *= devents->nu;
	    else if(ppars[i].fGC > DBL_EPSILON)
	      ppars[i].fGC *= devents->nu;
	    for(r=0; r<gpars.R; r++){
	      ppars[i].GAMMA[r] *= devents->nu;
	      ppars[i].lambdaN[r] /= devents->nu;
	      ppars[i].lambdaP[r] /= devents->nu;
	      ppars[i].normMean[r] *= devents->nu;
	      ppars[i].normVar[r] *= (devents->nu*devents->nu);
	    }
	    ppars[i].popAlpha = 0.0; /* growth rate to zero! */
	    for(j=0; j<gpars.NPOP; j++)
	      if(i != j && gpars.mig_mat[i][j] > DBL_EPSILON)
		gpars.mig_mat[i][j] *= devents->nu;
	  }
	}
	else{
	  ppars[devents->popi].Nt *= devents->nu;
	  ppars[devents->popi].THETA *= devents->nu;
	  ppars[devents->popi].INSRATE *= devents->nu;
	  ppars[devents->popi].DELRATE *= devents->nu;
	  if(ppars[devents->popi].RHO > DBL_EPSILON)
	    ppars[devents->popi].RHO *= devents->nu;
	  else if(ppars[devents->popi].fGC > DBL_EPSILON)
	    ppars[devents->popi].fGC *= devents->nu;
	  for(r=0; r<gpars.R; r++){
	    ppars[devents->popi].GAMMA[r] *= devents->nu;
	    ppars[devents->popi].lambdaN[r] /= devents->nu;
	    ppars[devents->popi].lambdaP[r] /= devents->nu;
	    ppars[devents->popi].normMean[r] *= devents->nu;
	    ppars[devents->popi].normVar[r] *= (devents->nu*devents->nu);
	  }
	  ppars[devents->popi].popAlpha = 0; /* growth rate to zero! */
	  for(i=0; i<gpars.NPOP; i++)
	    if(devents->popi != i && 
	       gpars.mig_mat[devents->popi][i] > DBL_EPSILON) 
	      gpars.mig_mat[devents->popi][i] *= devents->nu;
	}
      }
      else if(devents->eventType == 5){ /* kill population */
	if(devents->popi == -1){
	  for(i=0; i<gpars.NPOP; i++)  ppars[i].ALIVE = 0;	
	}
	else{
	  ppars[devents->popi].ALIVE = 0;
	}
      }
      else if(devents->eventType == 6){ /* change population parameter */
#ifdef UNIT_TEST
	fprintf(errfile,"eventType 6: t=%ld, ",t);
	fflush(errfile);
#endif
	if(devents->parIndex == 0){ /* specific population size */
	  double foo;
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].Nt = devents->newP.Nt;
	      foo = devents->newP.Nt/ppars[i].N;
	      ppars[i].THETA *= foo;
	      ppars[i].INSRATE *= foo;
	      ppars[i].DELRATE *= foo;
	      if(ppars[i].RHO > DBL_EPSILON)
		ppars[i].RHO *= foo;
	      else if(ppars[i].fGC > DBL_EPSILON)
		ppars[i].fGC *= foo;
	      for(r=0; r<gpars.R; r++){
		ppars[i].GAMMA[r] *= foo;
		ppars[i].lambdaN[r] /= foo;
		ppars[i].lambdaP[r] /= foo;
		ppars[i].normMean[r] *= foo;
		ppars[i].normVar[r] *= (foo*foo);
	      }
	      ppars[i].popAlpha = 0.0;
	      for(j=0; j<gpars.NPOP; j++)
		if(i != j && gpars.mig_mat[i][j] > DBL_EPSILON)
		  gpars.mig_mat[i][j] *= foo;
	    }
	  }
	  else{
	    ppars[devents->popi].Nt = devents->newP.Nt;
	    foo = devents->newP.Nt/ppars[devents->popi].N;
	    ppars[devents->popi].THETA *= foo;
	    ppars[devents->popi].INSRATE *= foo;
	    ppars[devents->popi].DELRATE *= foo;
	    if(ppars[devents->popi].RHO > DBL_EPSILON)
	      ppars[devents->popi].RHO *= foo;
	    else if(ppars[devents->popi].fGC > DBL_EPSILON)
	      ppars[devents->popi].fGC *= foo;
	    for(r=0; r<gpars.R; r++){
	      ppars[devents->popi].GAMMA[r] *= foo;
	      ppars[devents->popi].lambdaN[r] /= foo;
	      ppars[devents->popi].lambdaP[r] /= foo;
	      ppars[devents->popi].normMean[r] *= foo;
	      ppars[devents->popi].normVar[r] *= (foo*foo);
	    }
	    ppars[devents->popi].popAlpha = 0.0;
	    for(i=0; i<gpars.NPOP; i++)
	      if(i != devents->popi && 
		 gpars.mig_mat[devents->popi][i] > DBL_EPSILON)
		gpars.mig_mat[devents->popi][i] *= foo;
	  }
	}
	else if(devents->parIndex == 1){ /* change exponential growth rate */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].popAlpha = devents->newP.popAlpha;
	    }
	  }
	  else{
	    ppars[devents->popi].popAlpha = devents->newP.popAlpha;
	  }
	}
	else if(devents->parIndex == 2){ /* logistic growth rate */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].popAlpha = devents->newP.popAlpha;
	      ppars[i].K = devents->newP.K;
	      ppars[i].tauAlpha = t;
	      ppars[i].P0 = ppars[i].N;
	    }
	  }
	  else{
	    ppars[devents->popi].popAlpha = devents->newP.popAlpha;
	    ppars[devents->popi].K = devents->newP.K;
	    ppars[devents->popi].tauAlpha = t;
	    ppars[devents->popi].P0 = ppars[devents->popi].N;
	  }
	}
	else if(devents->parIndex == 3){ /* THETA */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      ppars[i].THETA = devents->newP.THETA;
	  }
	  else
	    ppars[devents->popi].THETA = devents->newP.THETA;
	}
	else if(devents->parIndex == 4){ /* RHO */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      ppars[i].RHO = devents->newP.RHO;
	  }
	  else
	    ppars[devents->popi].RHO = devents->newP.RHO;
	}
	else if(devents->parIndex == 5){ /* SELF */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      ppars[i].SELF = devents->newP.SELF;
	  }
	  else
	    ppars[devents->popi].SELF = devents->newP.SELF;
	}
	else if(devents->parIndex == 7){ /* distribution selective effects */
	  i = (devents->popi == -1 ? 0 : devents->popi);
	  for(; (devents->popi == -1 && i<gpars.NPOP) || i==devents->popi; i++){
	    if(devents->locus == -1){
	      for(r=0; r<gpars.R; r++){
		ppars[i].selDistType[r] = devents->newP.selDistType[r];
		if(ppars[i].selDistType[r] != 0)
		  ppars[i].neutpop = 0;
	      }
	      j=0;
	    }
	    else{
	      j=devents->locus;
	      ppars[i].selDistType[j] = devents->newP.selDistType[j];
	      if(ppars[i].selDistType[j] != 0)	
		ppars[i].neutpop = 0;
	    }
	    switch(ppars[i].selDistType[j]){
	    case 0:
	      break;
	    case 1:
	      if(devents->locus == -1){
		for(r=0; r<gpars.R; r++){
		  ppars[i].GAMMA[r] = devents->newP.GAMMA[r];
		  ppars[i].ProbPos[r] = devents->newP.ProbPos[r];
		  ppars[i].ProbDel[r] = devents->newP.ProbDel[r];
		  ppars[i].ProbNeut[r] = devents->newP.ProbNeut[r];
		}
	      }
	      else{
		ppars[i].GAMMA[devents->locus] = 
		  devents->newP.GAMMA[devents->locus];
		ppars[i].ProbPos[devents->locus] = 
		  devents->newP.ProbPos[devents->locus];
		ppars[i].ProbDel[devents->locus] = 
		  devents->newP.ProbDel[devents->locus];
		ppars[i].ProbNeut[devents->locus] = 
		  devents->newP.ProbNeut[devents->locus];
	      }
	      break;
	    case 2:
	      if(devents->locus == -1){
		for(r=0; r<gpars.R; r++){
		  ppars[i].ProbNegGamma[r] = devents->newP.ProbNegGamma[r];
		  ppars[i].alphaN[r] = devents->newP.alphaN[r];
		  ppars[i].alphaP[r] = devents->newP.alphaP[r];
		  ppars[i].lambdaN[r] = devents->newP.lambdaN[r];
		  ppars[i].lambdaP[r] = devents->newP.lambdaP[r];
		}
	      }
	      else{
		j = devents->locus;
		ppars[i].ProbNegGamma[j] = devents->newP.ProbNegGamma[j];
		ppars[i].alphaN[j] = devents->newP.alphaN[j];
		ppars[i].alphaP[j] = devents->newP.alphaP[j];
		ppars[i].lambdaN[j] = devents->newP.lambdaN[j];
		ppars[i].lambdaP[j] = devents->newP.lambdaP[j];
	      }
	      break;
	    case 3:
	      if(devents->locus == -1){
		for(r=0; r<gpars.R; r++){
		  ppars[i].normMean[r] = devents->newP.normMean[r];
		  ppars[i].normVar[r] = devents->newP.normVar[r];
		}
	      }
	      else{
		j = devents->locus;
		ppars[i].normMean[j] = devents->newP.normMean[j];
		ppars[i].normVar[j] = devents->newP.normVar[j];
	      }
	      break;
	    case 4:
	      break;
	    default:
	      if(devents->locus == -1)
		fprintf(errfile,"sfs_code error(%ld): bad selDistType = %d\n",
			INITIALSEED,ppars[i].selDistType[0]);
	      else
		fprintf(errfile,"sfs_code error(%ld): bad selDistType[%ld].[%ld] = %d\n",
			INITIALSEED,i,devents->locus,ppars[i].selDistType[devents->locus]);
	      abort();
	      break;
	    }
	  }
	}
	else if(devents->parIndex == 8){ /* KAPPA */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].KAPPA = devents->newP.KAPPA;
	      for(j=0; j<4; j++){
		for(k=0; k<4; k++){
		  ppars[i].transProb[j][k] = 
		    (j==k ? 0 : (j==((k+2)%4) ? 
				 ppars[i].KAPPA*ppars[i].baseFreq[k] :
				 ppars[i].baseFreq[k]));
		  if(k>0)
		    ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
		}
		for(k=0; k<4; k++){
		  ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
		}
	      }
	    }
	  }
	  else{
	    ppars[devents->popi].KAPPA = devents->newP.KAPPA;
	    i = devents->popi;
	    for(j=0; j<4; j++){
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] = 
		  (j==k ? 0 : (j==((k+2)%4) ? 
			       ppars[i].KAPPA*ppars[i].baseFreq[k] :
			       ppars[i].baseFreq[k]));
		if(k>0)
		  ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
	      }
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	      }
	    }
	  }
	}
	else if(devents->parIndex == 9){ /* PSI */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      ppars[i].PSI = devents->newP.PSI;
	  }
	  else
	    ppars[devents->popi].PSI = devents->newP.PSI;
	}
	else if(devents->parIndex == 10){ /* f0, selective constraint */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      if(devents->locus == -1)
		for(r=0; r<gpars.R; r++){
		  ppars[i].f0[r] = devents->newP.f0[r];
		  if(ppars[i].f0[r] <= 1-FLT_EPSILON)
		    ppars[i].neutpop = 0;
		}
	      else{
		ppars[i].f0[devents->locus] = devents->newP.f0[devents->locus];
		if(ppars[i].f0[devents->locus] <= 1-FLT_EPSILON)
		  ppars[i].neutpop = 0;
	      }
	    }
	  }
	  else{
	    if(devents->locus == -1){
	      for(r=0; r<gpars.R; r++){
		ppars[devents->popi].f0[r] = devents->newP.f0[r];
		if(ppars[devents->popi].f0[r] <= 1-FLT_EPSILON)
		  ppars[devents->popi].neutpop = 0;
	      }
	    }
	    else{
	      j=devents->locus;
	      ppars[devents->popi].f0[j] = devents->newP.f0[j];
	      if(ppars[devents->popi].f0[j] <= 1-FLT_EPSILON)
		ppars[devents->popi].neutpop = 0;
	    }
	  }
	}
	else if(devents->parIndex == 11){ /* pMaleMig */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      ppars[i].pMaleMig = devents->newP.pMaleMig;
	  }
	  else
	    ppars[devents->popi].pMaleMig = devents->newP.pMaleMig;
	}
	else if(devents->parIndex == 12){ /* GenEffect */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      ppars[i].GenEffect = devents->newP.GenEffect;
	  }
	  else
	    ppars[devents->popi].GenEffect = devents->newP.GenEffect;
	}
	else if(devents->parIndex == 13){ /* indel */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].INSRATE = devents->newP.INSRATE;
	      ppars[i].DELRATE = devents->newP.DELRATE;
	      ppars[i].INDELlength = devents->newP.INDELlength;
	    }
	  }
	  else{
	    ppars[devents->popi].INSRATE = devents->newP.INSRATE;
	    ppars[devents->popi].DELRATE = devents->newP.DELRATE;
	    ppars[devents->popi].INDELlength = devents->newP.INDELlength;
	  }
	}
	else if(devents->parIndex == 14){ /* longIndel */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].longINSRATE = devents->newP.longINSRATE;
	      ppars[i].longDELRATE = devents->newP.longDELRATE;
	      ppars[i].longINDELlength = devents->newP.longINDELlength;
	    }
	  }
	  else{
	    ppars[devents->popi].longINSRATE = devents->newP.longINSRATE;
	    ppars[devents->popi].longDELRATE = devents->newP.longDELRATE;
	    ppars[devents->popi].longINDELlength =devents->newP.longINDELlength;
	  }
	}
	else if(devents->parIndex == 15){ /* inversion */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].INVRATE = devents->newP.INVRATE;
	      ppars[i].INVlength = devents->newP.INVlength;
	    }
	  }
	  else{
	    ppars[devents->popi].INVRATE = devents->newP.INVRATE;
	    ppars[devents->popi].INVlength =devents->newP.INVlength;
	  }
	}
	else if(devents->parIndex == 16){ /* RHO */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      ppars[i].BGC = devents->newP.BGC;
	      ppars[i].fGC = devents->newP.fGC;
	      ppars[i].GCtract = devents->newP.GCtract;
	    }
	  }
	  else
	    ppars[devents->popi].RHO = devents->newP.RHO;
	}
	else if(devents->parIndex == 17){ /* pMaleRec */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++)
	      ppars[i].pMaleRec = devents->newP.pMaleRec;
	  }
	  else
	    ppars[devents->popi].pMaleRec = devents->newP.pMaleRec;
	}
	else if(devents->parIndex == 18){ /* baseFreq */
	  if(devents->popi == -1){
	    for(i=0; i<gpars.NPOP; i++){
	      for(j=0; j<4; j++){
		ppars[i].baseFreq[j] = devents->newP.baseFreq[j];
	      }
	      for(j=0; j<4; j++){
		for(k=0; k<4; k++){
		  if(gpars.substMod == 6){
		    int t = j+k-(j==0||k==0 ? 1 : 0);
		    ppars[i].transProb[j][k] =
		      (j==k ? 0 : ppars[i].GTR[t]*ppars[i].baseFreq[k]);
		    if(k>0)
		      ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
		  }
		  else{
		    ppars[i].transProb[j][k] = 
		      (j==k ? 0 : (j==((k+2)%4) ? 
				   ppars[i].KAPPA*ppars[i].baseFreq[k] :
				   ppars[i].baseFreq[k]));
		    if(k>0)
		      ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
		  }
		}
		for(k=0; k<4; k++){
		  ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
		}
	      }
	    }
	  }
	  else{
	    for(j=0; j<4; j++)
	      ppars[devents->popi].baseFreq[j] = devents->newP.baseFreq[j];
	    i = devents->popi;
	    for(j=0; j<4; j++){
	      for(k=0; k<4; k++){
		if(gpars.substMod == 6){
		  int t = j+k-(j==0||k==0 ? 1 : 0);
		  ppars[i].transProb[j][k] =
		    (j==k ? 0 : ppars[i].GTR[t]*ppars[i].baseFreq[k]);
		  if(k>0)
		    ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
		}
		else{
		  ppars[i].transProb[j][k] = 
		    (j==k ? 0 : (j==((k+2)%4) ? 
				 ppars[i].KAPPA*ppars[i].baseFreq[k] :
				 ppars[i].baseFreq[k]));
		  if(k>0)
		    ppars[i].transProb[j][k] += ppars[i].transProb[j][k-1];
		}
	      }
	      for(k=0; k<4; k++){
		ppars[i].transProb[j][k] /= ppars[i].transProb[j][3];
	      }
	    }
	  }
	}
	else{
	  fprintf(errfile,"unfortunately, parIndex %d is not implemented\n",
		  devents->parIndex);
	  abort();
	}
      }
      else{
	fprintf(errfile,"sfs_code error:  Unfortunately, eventType %d has not \
been implemented yet. Please contact author.\n",devents->eventType);
	abort();
      }
      tmp = devents->nextEvent;
      free(devents);
      devents = tmp;
      if(devents != NULL)
	nextDevent = devents->tau*PNanc;
      continue;
    }
    
    CHANGED = 0;
    if(gpars.autoRestart && (ITERFAILED || ITERSUCCESS)){
      t++;
      continue; /* rapidly finish sim */
    }
    for(pop=0; pop<gpars.NPOP; pop++){
      if(!ppars[pop].ALIVE)  continue;
      else if(ppars[pop].GenEffect<0 && t%abs(ppars[pop].GenEffect)!=0)
	continue;
      for(genEff=0; ((ppars[pop].GenEffect<0 && genEff==0) || 
		     genEff<ppars[pop].GenEffect); genEff++){
	popn[pop].gen++;
	if(ppars[pop].N <= 1){
	  fprintf(errfile,"error, population %d size decreased to <=1 at generation %ld, which is doom for sexually reproducing species.  If you are simulating a selfing population, contact the author of the program for assistance in enabling this feature\n", pop, popn[pop].gen);
	  abort();
	}
	if(gpars.PRINTGEN == 1){
	  printf("***   GEN=%ld; POP %d; GE %d; N=%ld; MIG=%d; \
ISEED=%ld(%ld); TRACK=%d  ******\n",popn[pop].gen, pop, genEff, ppars[pop].N, MIG,
		 INITIALSEED, gpars.seed, gpars.TRACKTRAJ);
	  fflush(stdout);
	}
#ifdef UNIT_TEST
	fprintf(errfile,"********   GENERATION  %ld; POP %d; GE %d; N=%ld; \
MIG=%d; ISEED=%ld(%ld); TRACK=%d  ******\n",t, pop, genEff, ppars[pop].N, MIG,
		INITIALSEED, gpars.seed, gpars.TRACKTRAJ);
	fflush(errfile);
#endif
	if(ppars[pop].N > gpars.maxN){
	  assert(selfed = realloc(selfed, ppars[pop].N*sizeof(*selfed)));
	  assert(selfing = realloc(selfing, ppars[pop].N*sizeof(*selfing)));
	  assert(parStore=realloc(parStore,gpars.P*ppars[pop].N*sizeof(long*)));
	  for(i=gpars.P*gpars.maxN; i<gpars.P*ppars[pop].N; i++)
	    assert(parStore[i] = malloc(gpars.R*sizeof(long)));
	  for(i=0; i<ppars[pop].N; i++){
	    selfed[i] = selfing[i] = 0;
	  }
	  gpars.maxN = ppars[pop].N;
	}
#ifdef UNIT_TEST
	checkForErrors(pop, "Beginning of generation");
#endif /* UNIT_TEST */
      
	/* Selecting parents (with replacement) for the next generation */
	partition(part, relFit, selClass, pop);
/* 	if(ppars[pop].neutpop==0){ */
/* 	  for(i=0; i<ppars[pop].N; i++){ */
/* 	    printf("%ld (%f %f) %.10f\n",i,popn[pop].indfit[popn[pop].parentLoc[gpars.P*i][0]][0],popn[pop].indfit[popn[pop].parentLoc[gpars.P*i+1][0]][0],part[i]); */
/* 	    fflush(stdout); */
/* 	  } */
/* 	  getchar(); */
/* 	} */
	/* update population size if necessary */
	if(fabs(ppars[pop].popAlpha) > DBL_EPSILON && CHANGED <= genEff){
	  double foo = 1.0;
	  if(ppars[pop].tauAlpha == -1){
	    foo = exp(ppars[pop].popAlpha/PNanc);
	  }
	  else if(t > ppars[pop].tauAlpha){
	    foo = ppars[pop].popAlpha;
	    foo += log(ppars[pop].K+ppars[pop].P0*
		       (exp(ppars[pop].popAlpha*
			    (t-ppars[pop].tauAlpha-1.0))-1));
	    foo -= log(ppars[pop].K+ppars[pop].P0*
		       (exp(ppars[pop].popAlpha*(t-ppars[pop].tauAlpha))-1.));
	    foo = exp(foo);
	    if(fabs(1.0-foo) < DBL_EPSILON && ppars[pop].Nt > ppars[pop].K){
	      foo = 1.0;
	      ppars[pop].popAlpha = 0.0;
	      ppars[pop].tauAlpha = -1;
	    }
	  }
	  ppars[pop].Nt *= foo;
	  ppars[pop].THETA *= foo;
	  ppars[pop].INSRATE *= foo;
	  ppars[pop].DELRATE *= foo;
	  if(ppars[pop].RHO > DBL_EPSILON)
	    ppars[pop].RHO *= foo;
	  else if(ppars[pop].fGC > DBL_EPSILON)
	    ppars[pop].fGC *= foo;
	  for(r=0; r<gpars.R; r++){
	    ppars[pop].GAMMA[r] *= foo;
	    ppars[pop].lambdaN[r] /= foo;
	    ppars[pop].lambdaP[r] /= foo;
	    ppars[pop].normMean[r] *= foo;
	    ppars[pop].normVar[r] *= (foo*foo);
	  }
	  for(i=0; i<gpars.NPOP; i++)
	    if(i != pop && gpars.mig_mat[pop][i] > DBL_EPSILON)
	      gpars.mig_mat[pop][i] *= foo;
	  CHANGED++;
	}
	if(fabs(ppars[pop].Nt - ppars[pop].N) > 1.0){ /* pop size changed */
#ifdef UNIT_TEST
	  fprintf(errfile,"N[%d](%ld) = %ld(%f; %ld; %ld)\n",pop,t,ppars[pop].N,
		  ppars[pop].Nt, ppars[pop].tauAlpha, ppars[pop].P0);
	  fflush(errfile);
#endif
	  if(floor(ppars[pop].Nt) > gpars.maxN){ /* largest population yet! */
	    /* must expand part, relFit, and pars */
	    assert(part = realloc(part, ((long)ppars[pop].Nt)*sizeof(*part)));
	    assert(relFit = realloc(relFit, ((long)ppars[pop].Nt)*
				    sizeof(*relFit)));
	    assert(Pars = realloc(Pars, gpars.P*((long)ppars[pop].Nt)*
				  sizeof(*Pars)));
	    assert(selfed = realloc(selfed,
				    ((long)ppars[pop].Nt)*sizeof(*selfed)));
	    assert(selfing = realloc(selfing,
				     ((long)ppars[pop].Nt)*sizeof(*selfing)));
	    assert(parStore=realloc(parStore,gpars.P*((long)ppars[pop].Nt)*
				    sizeof(long*)));
	    for(i=gpars.P*gpars.maxN; i<gpars.P*((long)ppars[pop].Nt); i++)
	      assert(parStore[i] = malloc(gpars.R*sizeof(long)));
	    for(i=0; i<((long)ppars[pop].Nt); i++){
	      selfed[i] = selfing[i] = 0;
	    }
	    gpars.maxN = (long)ppars[pop].Nt;
	  }
	  if(gpars.P*ppars[pop].Nt >= 2147483647){
	    fprintf(errfile, "sfs_code error(%ld): population size has grown \
too large.  Maximum allowed size is a long (2,147,483,647 chromosomes).	\
Please contact author if this is insufficient.\n",INITIALSEED);
	    abort();
	  }
#ifdef UNIT_TEST
	  printf("pop %d size %ld -> %ld\n",pop, ppars[pop].N, (long)(ppars[pop].Nt));
#endif
	  ChangePopSize(relFit, (long)(ppars[pop].Nt), pop, NULL, gpars.FITTEST);
	  partition(part, relFit, selClass, pop);
#ifdef UNIT_TEST
	  freeDead(pop, 1);
	  checkForErrors(pop, "After Size Change");
#endif /* UNIT_TEST */
	}
#ifdef VERBOSE_DEBUG
	{
	  fprintf(errfile,"\nNCbef%3ld: ",t);
	  for(i=0; i<=popn[pop].maxchr; i++){
	    fprintf(errfile,"%2ld ",popn[pop].numCopies[i][0]);
	  }
	  fprintf(errfile,"\n");
	  fprintf(errfile,"\nTREES BEFORE\n");
	  for(i=0; i<=popn[pop].maxchr; i++){
	    for(r=0; r<gpars.R; r++){
	      int foo=0;
	      carryMut(popn[pop].BigHead[i][r]->Rtree, 4, &foo);
	      if(foo>0){
		fprintf(errfile,"pop %d chr %ld locus %ld freq %ld carries mut\n",pop, i, r, popn[pop].numCopies[i][r]);
		long foo2=0;
		PrintHistTree(popn[pop].BigHead[i][r]->Rtree, &foo2, pop);
	      }
	    }
	  }
	}
#endif /* VERBOSE_DEBUG */
	
	/* select parents to populate new generation */
	Selecting(relFit, Pars, pop);

#ifdef VERBOSE_DEBUG
	{
	  fprintf(errfile,"\nTREES AFTER\n");
	  for(i=0; i<=popn[pop].maxchr; i++){
	    for(r=0; r<gpars.R; r++){
	      int foo=0;
	      carryMut(popn[pop].BigHead[i][r]->Rtree, 4, &foo);
	      if(foo>0){
		fprintf(errfile,"pop %d chr %ld locus %ld freq %ld carries mut\n",pop, i, r, popn[pop].numCopies[i][r]);
		long foo2=0;
		PrintHistTree(popn[pop].BigHead[i][r]->Rtree, &foo2, pop);
	      }
	    }
	  }
	  fprintf(errfile,"\nNCaft%3ld: ",t);
	  for(i=0; i<=popn[pop].maxchr; i++){
	    fprintf(errfile,"%2ld ",popn[pop].numCopies[i][0]);
	  }
	  fprintf(errfile,"\n");
	}
#endif /* VERBOSE_DEBUG */
#ifdef UNIT_TEST
/* 	updateMutFreqs(); */
/* 	checkForErrors(pop, "After Selecting"); */
/* 	for(r=0; r<gpars.R; r++){ */
/* 	  for(i=0; i<=popn[pop].maxchr; i++){ */
/* 	    if(popn[pop].numCopies[i][r] == 0 && */
/* 	       popn[pop].BigHead[i][r]->Rtree != NULL){ */
/* 	      fprintf(errfile, "error after NextGen: BH[%d][%ld][%ld]->Rtree != NULL\n", pop, i, r); */
/* 	      fprintf(errfile, "but numCopies = %ld\n", */
/* 		      popn[pop].numCopies[i][r]); */
/* 	      abort(); */
/* 	    } */
/* 	  } */
/* 	} */
#endif
#ifdef VERBOSE_DEBUG
	{
	  for(i=0; i<mutAr[0]->numMuts; i++){
	    if(mutAr[0]->muts[i]->event->free == '0'){
	      printf("-- %ld = site=%ld, ps=%d, carriers: ", i,
		     mutAr[0]->muts[i]->event->site,
		     popn[pop].polySites[0][mutAr[0]->muts[i]->event->site]);
	      for(j=0; j<mutAr[0]->muts[i]->event->maxHaps[pop]; j++)
		if(mutAr[0]->muts[i]->event->hapFreq[pop][j] != NULL)
		  printf("%ld ",*mutAr[0]->muts[i]->event->hapFreq[pop][j]);
	      printf("\n");
	    }
	  }
	}
#endif /* VERBOSE_DEBUG */

	/* perform migration into pop if necessary */
	if(MIG && !BURNING){
	  freeDead(pop, 2);
#ifdef UNIT_TEST
	  checkForErrors(pop, "before migrate");
#endif
	  Migrate(pop);

#ifdef VERBOSE_DEBUG
	  {
	    fprintf(errfile,"\nNCaft%3ld: ",t);
	    for(i=0; i<=popn[pop].maxchr; i++){
	      fprintf(errfile,"%2ld ",popn[pop].numCopies[i][0]);
	    }
	    fprintf(errfile,"\n");
	    fprintf(errfile,"PLaft%3ld: ",t);
	    for(i=0; i<gpars.P*ppars[pop].N; i++){
	      fprintf(errfile,"%2ld ",popn[pop].parentLoc[i][0]);
	    }
	    fprintf(errfile,"\n");
	  }
#endif /* VERBOSE_DEBUG */
#ifdef UNIT_TEST
/* 	  updateMutFreqs(); */
/* 	  checkForErrors(pop, "After Migrate"); */
#endif /* UNIT_TEST */
#ifdef VERBOSE_DEBUG
	  for(i=0; i<gpars.L[0]; i++)
	    fprintf(errfile,"%d;",popn[pop].polySites[0][i]);
	  fprintf(errfile,"\n");
#endif
	}
	
	/* Mutation */
	mutate(mutClass_site, mutClass_locus, selClass, fitQuant, pop);

	/* indels */
	indels(pop, fitQuant);
#ifdef UNIT_TEST /* check that everything went ok */
/* 	checkForErrors(pop, "After indels"); */
#endif

	/* inversions */
	inversions(pop, fitQuant);
#ifdef UNIT_TEST /* check that everything went ok */
/* 	checkForErrors(pop, "After inversions"); */
#endif

	freeDead(pop, 3);
#ifdef UNIT_TEST
	checkForErrors(pop, "end of pop gen");
#endif
      }
      CHANGED = 0;
    }
    if(gpars.TRACKTRAJ){
      long nmuts = 0;
      INFREQRANGE = 0;
      for(r=0; r<gpars.R; r++){
	if(gpars.trajLoc == -1 || gpars.trajLoc == r){
	  for(i=0; i<mutAr[r]->numMuts; i++){
	    if(mutAr[r]->muts[i]->event->free == '0'){
	      if(gpars.trajSite==-1 ||
		 gpars.trajSite==mutAr[r]->muts[i]->event->site){
		long cnt = 0;
		long nfixed = 0;
		nmuts++;
		for(pop=0; pop<gpars.NPOP; pop++){
		  if(gpars.trajPop == -1 || gpars.trajPop == pop){
		    if(mutAr[r]->muts[i]->event->fixed[pop] == '1'){
		      nfixed++;
		      if(gpars.trajMaxFreq>=1-FLT_EPSILON)
			INFREQRANGE = 1;
		    }
		    else{
		      cnt = 0;
		      for(k=0; k<mutAr[r]->muts[i]->event->maxHaps[pop]; k++){
			if(mutAr[r]->muts[i]->event->hapFreq[pop][k] != NULL){
			  cnt += *(mutAr[r]->muts[i]->event->hapFreq[pop][k]);
			}
		      }
		      mutAr[r]->muts[i]->event->numCarriers[pop] = cnt;
		      if(cnt == gpars.P*ppars[pop].N){
			nfixed++;
			freeDead(pop,7);
		      }
		      if(cnt >=
			 gpars.trajMinFreq*gpars.P*ppars[pop].N-FLT_EPSILON &&
			 cnt <=
			 gpars.trajMaxFreq*gpars.P*ppars[pop].N+FLT_EPSILON){
			INFREQRANGE = 1;
		      }
		    }
		  }
		}
		if(trajFile!=NULL && nfixed<gpars.NPOP){
		  if(tmpTRAJ == NULL){
		    fprintf(trajFile, "%ld %ld %ld %ld %c %c %1.4f ", 
			    (BURNING?(long)(-gpars.BURN*PNanc)+t:t), r, 
			    mutAr[r]->muts[i]->event->site,
			    mutAr[r]->muts[i]->event->gen,
			    mutAr[r]->muts[i]->event->ancNuc,
			    mutAr[r]->muts[i]->event->derNuc,
			    mutAr[r]->muts[i]->event->fit);
		    for(pop=0; pop<gpars.NPOP; pop++){
		      if(gpars.trajPop == -1 || gpars.trajPop == pop){
			if(mutAr[r]->muts[i]->event->fixed[pop] == '1'){
			  fprintf(trajFile, "%ld %ld ", gpars.P*ppars[pop].N,
				  gpars.P*ppars[pop].N);
			}
			else{
			  fprintf(trajFile, "%ld %ld ",
				  mutAr[r]->muts[i]->event->numCarriers[pop],
				  gpars.P*ppars[pop].N);
			}
		      }
		    }
		    fprintf(trajFile,"\n");
		  }
		  else{
		    fprintf(tmpTRAJ, "%ld %ld %ld %ld %c %c %1.4f ", 
			    (BURNING?(long)(-gpars.BURN*PNanc)+t:t), r, 
			    mutAr[r]->muts[i]->event->site,
			    mutAr[r]->muts[i]->event->gen,
			    mutAr[r]->muts[i]->event->ancNuc,
			    mutAr[r]->muts[i]->event->derNuc,
			    mutAr[r]->muts[i]->event->fit);
		    for(pop=0; pop<gpars.NPOP; pop++){
		      if(gpars.trajPop == -1 || gpars.trajPop == pop){
			if(mutAr[r]->muts[i]->event->fixed[pop] == '1'){
			  fprintf(tmpTRAJ, "%ld %ld ", gpars.P*ppars[pop].N,
				  gpars.P*ppars[pop].N);
			}
			else{
			  fprintf(tmpTRAJ, "%ld %ld ",
				  mutAr[r]->muts[i]->event->numCarriers[pop],
				  gpars.P*ppars[pop].N);
			}
		      }
		    }
		    fprintf(tmpTRAJ,"\n");
		  }
		}
		if(nfixed == gpars.NPOP){
		  if(gpars.trajSite != -1){
		    if(gpars.trajMaxFreq<1-FLT_EPSILON){
		      ITERFAILED = 1;
		    }
		    else{
		      ITERSUCCESS = 1;
		      INFREQRANGE = 1;
		    }
		    break;
		  }
		  else
		    continue;
		}
	      }
	    }
	    else if(gpars.trajSite != -1){
	      if(gpars.trajMinFreq > FLT_EPSILON){
		ITERFAILED = 1;
	      }
	      else
		INFREQRANGE = 1;
	      break;
	    }
	  }
	}
      }
      if(nmuts == 0){
	if(gpars.trajMinFreq>FLT_EPSILON){
	  ITERFAILED = 1;
	}
	else{
	  ITERSUCCESS = 1;
	  INFREQRANGE = 1;
	}
      }
    }
    t++;
#ifdef UNIT_TEST /* check that everything went ok */
    {
      fflush(stdout);
      fflush(outfile);
      fflush(errfile);
      for(pop=0; pop<gpars.NPOP; pop++){
	if(ppars[pop].ALIVE == 0)
	  continue;
	freeDead(pop,5);
	checkForErrors(pop, "end of gen");
	for(r=0; r<gpars.R; r++){
	  for(i=0; i<=popn[pop].maxchr; i++){
	    if(popn[pop].numCopies[i][r] == 0 &&
	       popn[pop].BigHead[i][r]->Rtree != NULL){
	      fprintf(errfile, "error end gen: BH[%d][%ld][%ld]->Rtree != NULL\n",
		      pop, i, r);
	      fprintf(errfile, "but numCopies = %ld\n",
		      popn[pop].numCopies[i][r]);
	      abort();
	    }
	  }
	}
      }
      fflush(stdout);
      fflush(outfile);
      fflush(errfile);
    }
#endif
  }

  if(!ITERFAILED){
    if(DOMEST) fprintf(errfile,"\n"); /* need spacer */
    if(tmpOUT == NULL){
      fprintf(outfile,"Nc:");
      for(pop=0; pop<gpars.NPOP; pop++){
	if(pop < gpars.NPOP-1)
	  fprintf(outfile,"%ld,",ppars[pop].N);
	else
	  fprintf(outfile,"%ld;\n",ppars[pop].N);
      }
      fflush(outfile);
    }
    else{
      fprintf(tmpOUT,"Nc:");
      for(pop=0; pop<gpars.NPOP; pop++){
	if(pop < gpars.NPOP-1)
	  fprintf(tmpOUT,"%ld,",ppars[pop].N);
	else
	  fprintf(tmpOUT,"%ld;\n",ppars[pop].N);
      }
      fflush(tmpOUT);
    }
    
    if(freqFile != NULL){ /* report population frequencies in ancestral pop! */
      assert(popFreqs = malloc(gpars.R*sizeof(double**)));
      assert(mutLocs = malloc(gpars.R*sizeof(double*)));
      assert(orderMutLocs = malloc(gpars.R*sizeof(long*)));
      for(r=0; r<gpars.R; r++){
	assert(popFreqs[r] = malloc(mutAr[r]->numMuts*sizeof(double*)));
	assert(mutLocs[r] = malloc(mutAr[r]->numMuts*sizeof(double)));
	assert(orderMutLocs[r] = malloc(mutAr[r]->numMuts*sizeof(long)));
	for(i=0; i<mutAr[r]->numMuts; i++){
	  assert(popFreqs[r][i] = malloc(gpars.NPOP*sizeof(double)));
	  mutLocs[r][i] = mutAr[r]->muts[i]->event->site;
	  orderMutLocs[r][i] = i;
	
	  for(pop=0; pop<gpars.NPOP; pop++){
	    if(mutAr[r]->muts[i]->event->free == '1')
	      popFreqs[r][i][pop] = 0.0;
	    else if(mutAr[r]->muts[i]->event->fixed[pop] == '1')
	      popFreqs[r][i][pop] = 1.0;
	    else{
	      mutAr[r]->muts[i]->event->numCarriers[pop] = 0.0;
	      for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[pop]; j++){
		if(mutAr[r]->muts[i]->event->hapFreq[pop][j] != NULL)
		  mutAr[r]->muts[i]->event->numCarriers[pop] +=
		    (*mutAr[r]->muts[i]->event->hapFreq[pop][j]);
	      }
	      popFreqs[r][i][pop] = (mutAr[r]->muts[i]->event->numCarriers[pop]/
				     (gpars.P*ppars[pop].N+.0));
	    }
	  }
	}
	sort2(mutAr[r]->numMuts-1, &mutLocs[r][0], &orderMutLocs[r][0]);
      }
    }
    
    /* Now sampling each pop akin to population contraction!! */
    if(tmpOUT == NULL)
      fprintf(outfile,"MALES:");
    else
      fprintf(tmpOUT,"MALES:");
#ifdef UNIT_TEST
    (*numMutSeg) = 0;
#endif
    for(pop=0; pop<gpars.NPOP; pop++){
      if(ppars[pop].SS == -1)
	popSS = ppars[pop].N;
      else if(ppars[pop].SS == 0){
	if(pop < gpars.NPOP-1){
	  if(tmpOUT == NULL)
	    fprintf(outfile,"0,");
	  else
	    fprintf(tmpOUT, "0,");
	}
	else{
	  if(tmpOUT == NULL)
	    fprintf(outfile,"0;\n");
	  else
	    fprintf(tmpOUT,"0;\n");
	}
	continue;
      }
      else
	popSS = ppars[pop].SS;
      
#ifdef UNIT_TEST
      for(i=0; i<=popn[pop].maxchr; i++)
	if(popn[pop].numCopies[i][0] < 0){
	  fprintf(errfile,"nc[%ld] = %ld (%ld)\n",i, popn[pop].numCopies[i][0],
		  INITIALSEED);
	  abort();
	}
#endif
      ChangePopSize(NULL, popSS, pop, survive, 0);
      updateMutFreqs();
      adjustCarriers(pop);
      if(tmpOUT == NULL){
	if(pop < gpars.NPOP-1)
	  fprintf(outfile,"%ld,",ppars[pop].MALES);
	else
	  fprintf(outfile,"%ld;\n",ppars[pop].MALES);
      }
      else{
	if(pop < gpars.NPOP-1)
	  fprintf(tmpOUT,"%ld,",ppars[pop].MALES);
	else
	  fprintf(tmpOUT,"%ld;\n",ppars[pop].MALES);
      }

      for(r=0; r<gpars.R; r++){
#ifdef UNIT_TEST
	countSeg(popn[pop].fixed[r]->Rtree, numMutSeg);
#endif
	storeInfo(&sampStore, popn[pop].fixed[r]->Rtree, pop, -1, r, 0,
		  popn[pop].gen);
#ifdef UNIT_TEST
	j=0;
	printf("fixed history pop=%d, locus=%ld:\n",pop, r);
	fflush(stdout);
	PrintHistTree(popn[pop].fixed[r]->Rtree, &j, pop);
	fflush(stdout);
#endif
	for(i=0; i<gpars.P*ppars[pop].SS; i++){
#ifdef UNIT_TEST
	  countSeg(popn[pop].BigHead[popn[pop].parentLoc[i][r]][r]->Rtree,
		   numMutSeg);
	  j=0;
#ifdef VERBOSE_DEBUG
	  printf("seg history pop=%d, locus=%ld; ind=%ld(%ld):\n",pop, r,
		 i, popn[pop].parentLoc[i][r]);
	  fflush(stdout);
	  PrintHistTree(popn[pop].BigHead[popn[pop].parentLoc[i][r]][r]->Rtree,
			&j, pop);
	  fflush(stdout);	
#endif
#endif
	  storeInfo(&sampStore,
		    popn[pop].BigHead[popn[pop].parentLoc[i][r]][r]->Rtree,
		    pop, i, r, 0, t);
	}
      }
#ifdef UNIT_TEST
      {
	long foo = 0;
	countSegStorage(sampStore, &foo);
	if(foo != (*numMutSeg)){
	  fprintf(errfile,"sfs_code error(%ld):  did not store mutations \
correctly! Expecting %ld mutations, stored %ld.\n",INITIALSEED,*numMutSeg, foo);
	  abort();
	}
	else{
	  printf("%ld derived alleles stored correctly...\n", *numMutSeg);
	}
      }
#endif
    }
    
    if(freqFile != NULL){
      long index = 0;
      for(r=0; r<gpars.R; r++){
	for(i=0; i<mutAr[r]->numMuts; i++){
	  index = orderMutLocs[r][i];
	  j = 0;
	  for(pop=0; pop<gpars.NPOP; pop++){
	    if(popFreqs[r][index][pop] > DBL_EPSILON){
	      j++;
	      break;
	    }
	  }
	  if(j == 0)
	    continue;

	  fprintf(freqFile,"%d\t%ld\t%.0f\t%ld\t",gpars.iter,r,mutLocs[r][i],
		  mutAr[r]->muts[index]->event->gen);
	  if(mutAr[r]->muts[index]->event->nonsyn == '0' || 
	     mutAr[r]->muts[index]->event->nonsyn == '1'){
	    fprintf(freqFile,"S\t%c\t%c\t",mutAr[r]->muts[index]->event->ancNuc,
		    mutAr[r]->muts[index]->event->derNuc);
	  }
	  else if(mutAr[r]->muts[index]->event->nonsyn == 'i'){
	    fprintf(freqFile,"I\t-\t+\t");
	  }
	  else if(mutAr[r]->muts[index]->event->nonsyn == 'd'){
	    fprintf(freqFile,"D\t+\t-\t");
	  }
	  else if(mutAr[r]->muts[index]->event->nonsyn == 'v'){
	    fprintf(freqFile,"D\tf\tr\t");
	  }
	  for(pop=0; pop<gpars.NPOP; pop++){
	    if(mutAr[r]->muts[index]->event->fixed[pop] == '0'){
	      mutAr[r]->muts[index]->event->numCarriers[pop] = 0;
	      for(j=0; j<mutAr[r]->muts[index]->event->maxHaps[pop]; j++){
		if(mutAr[r]->muts[index]->event->hapFreq[pop][j] != NULL){
		  mutAr[r]->muts[index]->event->numCarriers[pop] +=
		    (*mutAr[r]->muts[index]->event->hapFreq[pop][j]);
		}
	      }
	      if(mutAr[r]->muts[index]->event->numCarriers[pop] == 0){
		fprintf(freqFile,"0.0\t");
	      }
	      else if(mutAr[r]->muts[index]->event->numCarriers[pop] == 
		      gpars.P*ppars[pop].N){
		fprintf(freqFile,"1.0\t");
	      }
	      else{
		fprintf(freqFile,"%e\t",mutAr[r]->muts[index]->event->numCarriers[pop]/(gpars.P*ppars[pop].N+.0));
	      }
	    }
	    else{
	      fprintf(freqFile,"1.0\t");
	    }
	    if(popFreqs[r][index][pop] > 1-DBL_EPSILON)
	      fprintf(freqFile,"1.0\t");
	    else if(popFreqs[r][index][pop] < DBL_EPSILON)
	      fprintf(freqFile,"0.0\t");
	    else
	      fprintf(freqFile,"%e\t", popFreqs[r][index][pop]);
	  }
	  fprintf(freqFile,"\n");
	}
      }
      for(r=0; r<gpars.R; r++){
	free(mutLocs[r]);
	free(orderMutLocs[r]);
	for(i=0; i<mutAr[r]->numMuts; i++){
	  free(popFreqs[r][i]);
	}
	free(popFreqs[r]);
      }
      free(popFreqs);
      free(mutLocs);
      free(orderMutLocs);
      fflush(freqFile);
    }

    if(gpars.TRACKANC){
      for(pop=0; pop<gpars.NPOP; pop++){
	for(i=0; i<ppars[pop].N; i++){
	  for(k=0; k<gpars.P; k++){
	    fprintf(ancFile, "pop%dind%ldchr%ld ", pop, i, k);
	    for(r=0; r<gpars.R; r++){
	      for(j=0; j<gpars.L[r]; j++){
		fprintf(ancFile, "%d ",
			popn[pop].ancestry[popn[pop].parentLoc[gpars.P*i+k]
					   [r]][r][j]);
	      }
	    }
	    fprintf(ancFile, "\n");
	  }
	}
      }
    }
    
    count=0;
    if(tmpOUT == NULL){
      printStorage(sampStore, &count, outfile);
      fprintf(outfile,"\n");
      fflush(outfile);
    }
    else{
      printStorage(sampStore, &count, tmpOUT);
      fprintf(tmpOUT,"\n");
      fflush(tmpOUT);
    }    
#ifdef UNIT_TEST
    {
      long *hapFreq;
      int numHap = 0;
      int indicator = 0, max;
      char **haps = cmatrix(0,ppars[0].N,0,12*gpars.R);
      char *temp = cvector(0,12*gpars.R);
      char *tshort = cvector(0,100);
    
      assert(hapFreq = calloc(ppars[0].N, sizeof(*hapFreq)));
      for(i=0; i<ppars[0].N; i++){
	hapFreq[i] = 0;
	temp[0] = '\0';
	for(r=0; r<gpars.R; r++){
	  if(popn[0].parentLoc[2*i][r] <= popn[0].parentLoc[2*i+1][r]){
#ifdef WINDOWS
	    sprintf_s(tshort,100,"%ld %ld ",popn[0].parentLoc[2*i][r],
		      popn[0].parentLoc[2*i+1][r]);
#else
	    sprintf(tshort,"%ld %ld ",popn[0].parentLoc[2*i][r],
		    popn[0].parentLoc[2*i+1][r]);
#endif
	  }
	  else{
#ifdef WINDOWS
	    sprintf_s(tshort,100,"%ld %ld ",popn[0].parentLoc[2*i+1][r],
		      popn[0].parentLoc[2*i][r]);
#else
	    sprintf(tshort,"%ld %ld ",popn[0].parentLoc[2*i+1][r],
		    popn[0].parentLoc[2*i][r]);
#endif
	  }
#ifdef WINDOWS
	  strcat_s(temp, 12*gpars.R*sizeof(char), tshort);
#else
	  strcat(temp,tshort);
#endif
	}
	indicator = 0;
	for(j=0; j<=i; j++){
	  if(strcmp(haps[j], temp) == 0){
	    indicator = 1;
	    hapFreq[j]++;
	    break;
	  }
	}
	if(indicator == 0){
#ifdef WINDOWS
	  sprintf_s(haps[numHap], 12*gpars.R+1, temp);
#else
	  sprintf(haps[numHap], "%s", temp);
#endif
	  hapFreq[numHap]++;
	  numHap++;
	}
      }
    
      fprintf(errfile,"\nthere were %d diplotypes\n",numHap);
      fprintf(errfile,"diplotype frequencies:\n");
      max = 0;
      for(i=0; i<numHap; i++){
	max = 0;
	indicator = 0;
	for(j=0; j<numHap; j++){
	  if(hapFreq[j] >= max){
	    max = hapFreq[j];
	    indicator = j;
	  }
	}
	fprintf(errfile,"%d ",max);
	hapFreq[indicator] = -1;
      }
      fprintf(errfile,"\n");
    
      /* FREE ALLOCATED MEMORY */
      free_cvector(tshort,0,100);
      free_cvector(temp,0,12*gpars.R);
      free_cmatrix(haps,0,ppars[0].N,0,4*gpars.R);
      free(hapFreq);
    }
#endif
    
    freeStorage(sampStore);
  }
  for(r=0; r<gpars.R; r++){
    free(DomestAllele[r]);
  }
  free(DomestAllele);
  free(survive);
  for(i=1; i<gpars.NPOP; i++){
    for(r=0; r<gpars.R; r++){
      free(popn[i].conSeq[r]);
    }
    free(popn[i].conSeq);
  }
  for(i=0; i<gpars.NPOP; i++){
    for(j=0; j<=popn[i].maxchr; j++){
      for(r=0; r<gpars.R; r++){
	if(i == 0 && j == 0){
	  pushTrash(&popn[i].BigHead[j][r]->Rtree);
	  popn[i].BigHead[j][r]->Rtree = NULL;
	}
	else{
	  pushTrash(&popn[i].BigHead[j][r]->Rtree);
	  popn[i].BigHead[j][r]->Rtree = NULL;
	  free(popn[i].BigHead[j][r]);
	  if(gpars.TRACKANC)
	    free(popn[i].ancestry[j][r]);
	}
	free(popn[i].extremeMuts[j][r]);
	if(gpars.substMod == 5)  free(popn[i].hpSUM[j][r]);
	if(j == 0){
	  free(popn[i].polySites[r]);
	  pushTrash(&popn[i].fixed[r]->Rtree);
	  popn[i].fixed[r]->Rtree = NULL;
	  free(popn[i].fixed[r]);
	  popn[i].fixed[r] = NULL;
	}
      }
      if(i > 0 || j > 0){ 
	free(popn[i].BigHead[j]);
	popn[i].BigHead[j] = NULL;
	if(gpars.TRACKANC){
	  free(popn[i].ancestry[j]);
	}
      }
      free(popn[i].numCopies[j]);
      free(popn[i].chrGen[j]);
      free(popn[i].indfit[j]);
      free(popn[i].extremeMuts[j]);
      if(gpars.substMod == 5)  free(popn[i].hpSUM[j]);
    }
    if(i > 0){
      free(popn[i].BigHead);
      popn[i].BigHead = NULL;
      if(gpars.TRACKANC)
	free(popn[i].ancestry);
    }
    free(popn[i].nextchr);
    free(popn[i].numCopies);
    free(popn[i].chrGen);
    free(popn[i].indfit);
    free(popn[i].polySites);
    free(popn[i].fixed);
    free(popn[i].extremeMuts);
    if(gpars.substMod == 5)  free(popn[i].hpSUM);
    for(j=0; j<gpars.P*ppars[i].N; j++)
      free(popn[i].parentLoc[j]);
    free(popn[i].parentLoc);
  }
  /*   clearMutArray(); */
  free_dvector(part,0,gpars.maxN-1);
  free_dvector(relFit,0,gpars.maxN-1);
  free(Pars);
  free_dvector(rates_site,0,ppars[0].RateClassSites-1);
  free_dvector(rates_loci,0,ppars[0].RateClassLoci-1);
  free(selfed);
  free(selfing);
  for(i=0; i<gpars.P*gpars.maxN; i++)
    free(parStore[i]);
  free(parStore);
}

/* ------------------------------------------------------------------------- */

void partition(double *part, double *relFit, const int *selClass,const int pop)
{
  int p;
  long i, r;
  double sumfit = 0.0;	/* The sum of all fitness */

#ifdef UNIT_TEST
  double tfit = 1.0;
#endif

  if(ppars[pop].neutpop){
    part[0] = relFit[0] = 1/(double)ppars[pop].N;
    for(i=1; i<ppars[pop].N; i++){ /* complete neutrality */
      part[i] = 1.0/ppars[pop].N;
      relFit[i] = relFit[i-1]+1.0/ppars[pop].N;
    }
    return;
  }

#ifdef UNIT_TEST
  {
    long popsize = gpars.P*ppars[pop].N;
    for(r=0; r<gpars.R; r++){
      char *check = cvector(0,popsize-1);
      for(i=0; i<popsize; i++)
	check[i] = '0';
      for(i=0; i<popsize; i++){
	if(check[popn[pop].parentLoc[i][r]] == '0'){
	  double tfit2 = (gpars.ADDITIVE ? 0.0 : 1.0);
	  check[popn[pop].parentLoc[i][r]] = '1';
	  getFitTree(popn[pop].BigHead[popn[pop].parentLoc[i][r]][r]->Rtree,
		     &tfit2, pop, r);
	  if(tfit2 > FLT_EPSILON && 
	     fabs(tfit2 - popn[pop].indfit[popn[pop].parentLoc[i][r]][r]) >
	     FLT_EPSILON){
	    fprintf(errfile,"error in partition, before getting started(%ld): ",
		    INITIALSEED);
	    fprintf(errfile,"tfit2 = %f; fits[%ld][%ld] = %f\n",
		    tfit2,i,r,popn[pop].indfit[popn[pop].parentLoc[i][r]][r]);
	    abort();
	  }
	}
      }
      free(check);
    }
  }
#endif
  
  for(i=0; i<ppars[pop].N; i++){
    part[i] = 1.0;
#ifdef UNIT_TEST
    tfit = 1.0;
#endif
    
    for(p=0; p<gpars.P; p++){
      for(r=0; r<gpars.R; r++){
	if(ppars[pop].selDistType[r] == 0 && ppars[pop].f0[r] >= 1-FLT_EPSILON
	   && gpars.SKIPSITES==NULL)
	  continue;
	if(gpars.ADDITIVE)
	  part[i] += popn[pop].indfit[popn[pop].parentLoc[gpars.P*i+p][r]][r];
	else
	  part[i] *= popn[pop].indfit[popn[pop].parentLoc[gpars.P*i+p][r]][r];
	if(part[i] < 0){
	  part[i] = 0.0;/* if fitness < EPS, lethal mutation! */
	  p = gpars.P;  /* no need to continue */
	  break;
	}
#ifdef UNIT_TEST
	long foo;
	getFitTree(popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i+p][r]]
		   [r]->Rtree, &tfit, pop, r);
	if(tfit < DBL_EPSILON)  tfit = 0.0;
	if(fabs(tfit - part[i]) > FLT_EPSILON){
	  fprintf(errfile, "fitnesses not right for %ld (%ld) locus %ld:\n",
		  popn[pop].parentLoc[gpars.P*i+p][r], gpars.P*i+p, r);
	  fprintf(errfile,"expecting %f, obs %f (ADD=%d)\n",tfit, part[i],
		  gpars.ADDITIVE);
	  fflush(errfile);
	  foo = 0;
	  PrintHistTree(popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i+p][r]]
			[r]->Rtree, &foo, pop);
	  fflush(stdout);
	  abort();
	}
#endif
      }
    }
    sumfit += part[i];
#ifdef UNIT_TEST
    if(part[i] > FLT_EPSILON && fabs(part[i] - tfit) > FLT_EPSILON){
      fprintf(errfile,"error in partition(%ld): part[%ld] = %lf; tfit = %lf\n",
	      INITIALSEED,i,part[i],tfit);
      abort();
    }
#endif
  }
  
  if(sumfit > FLT_EPSILON){
    part[0] /= sumfit;
    relFit[0] = part[0];
    for(i=1; i<ppars[pop].N; i++){
      part[i] /= sumfit;
      relFit[i] = relFit[i-1]+part[i];
    }
  }
  else{
    fprintf(errfile,"the entire population is effectively dead (%ld)...  \
generation %ld; too much negative selection! sumfit = %lf\n",
	    INITIALSEED, popn[pop].gen, sumfit);
    sumfit = 0;
    for(i=0; i<=popn[pop].maxchr; i++){
      sumfit += part[i];
      printf("i=%ld; fit=%lf (%f)\n",i, part[i], sumfit);
      fflush(stdout);
      for(r=0; r<gpars.R; r++){
	long foo = 0;
	if(popn[pop].numCopies[i][r] == 0 || 
	   popn[pop].BigHead[i][r]->Rtree == NULL ||
	   fabs(1-popn[pop].indfit[i][r]) <= FLT_EPSILON)
	  continue;
	printf("r=%5ld, n=%5ld, fit=%5.5lf\n",r, popn[pop].numCopies[i][r],
	       popn[pop].indfit[i][r]);
	PrintHistTree(popn[pop].BigHead[i][r]->Rtree, &foo, pop);
      }
    }
    abort();
  }
}

/* ------------------------------------------------------------------------- */

void Selecting(double *relFit, long *Pars, const int pop)
{
  long i, n, r, NSELF, NMSELF, NFSELF, iW=0, iW2=0, iB=0, iB2=0;
  float s;
  long p, rem, remS, remMale, remFemale, numXW=0, numXB=0;
  struct locRec *within=NULL;
  struct interRec *between=NULL;

  /* generate recombination positions among offspring */
  distribRec(pop, &within, &numXW, &between, &numXB);
  
#ifdef UNIT_TEST
#ifdef VERBOSE_DEBUG
  printf("recombinants: %ld\n", numXW);
  for(i=0; i<numXW; i++){
    printf("IND=%ld:\n",within[i].id);
    for(n=0; n<within[i].numEvents; n++){
      printf("\t%ld.%ld(%d)\n",within[i].loci[n], within[i].pos[n],
	     within[i].type[n]);
    }
  }
  fflush(stdout);
#endif
#endif
  
  /* Note that for haploid, only 1 chromosome.
     for diploid mating, each parent leaves 1 chromosome to their child.
     for tetraploid mating, each parent leaves 2 chromosomes to child. */
  NSELF = (long)(ppars[pop].SELF*ppars[pop].N); /* selfed individuals */
#ifdef UNIT_TEST
  printf("selfing %ld individuals\n",NSELF);
  fflush(stdout);
#endif
  if(NSELF > 0){
    NMSELF = NFSELF = 0; /* males versus females */
    rem = NSELF;
    while(rem > 0){
      n = invCDF(ran1(&gpars.seed), relFit, 0, ppars[pop].N-1);
      /* accept a female selfing up to the total number of females */
      if(n < ppars[pop].MALES && NFSELF < ppars[pop].MALES){
	selfed[n]++;
	rem--;
	NFSELF++;
      }
      else if(n >= ppars[pop].MALES &&
	      NMSELF < ppars[pop].N-ppars[pop].MALES){
	selfed[n]++;
	NMSELF++;
	rem--;
      }
      else{
	continue;
      }
      /* find offspring to be produced by selfing */
      if(selfing[n] == 0)
	selfing[n] = 1;
      else{
	if(n<ppars[pop].MALES){
	  while(1){
	    n = (long)(ran1(&gpars.seed)*ppars[pop].MALES);
	    if(selfing[n] == 0){
	      selfing[n] = 1;
	      break;
	    }
	  }
	}
	else{
	  while(1){
	    n = (long)(ran1(&gpars.seed)*(ppars[pop].N-ppars[pop].MALES)+
		       ppars[pop].MALES);
	    if(selfing[n] == 0){
	      selfing[n] = 1;
	      break;
	    }
	  }
	}
      }
    }
#ifdef UNIT_TEST
    {
      fprintf(errfile,"NSELF = %ld\n",NSELF);
      fflush(errfile);
      assert(NSELF <= ppars[pop].N);
      rem = 0;
      for(i=0; i<ppars[pop].N; i++){
	rem += selfing[i];
      }
      assert(rem == NSELF);
    }
#endif /* UNIT_TEST */
  }

  /* now loop over individuals and identify parents */
  remS = 0;
  for(i=0; i<ppars[pop].N; i++){
    if(selfing[i]){  /* this individual to be selfed */
      while(selfed[remS] == 0){
	if(i<ppars[pop].MALES)
	  remS = (remS+1)%ppars[pop].MALES;
	else
	  remS = (remS+1<ppars[pop].N ? remS+1 : ppars[pop].MALES);
      }
      if(gpars.P == 2){
	if(i >= ppars[pop].MALES && gpars.SEX[0] == '1'){
	  Pars[2*i] = 2*remS;
	  Pars[2*i+1] = 2*remS+1;
	} /* males must give X and Y */
	else if((s = ran1(&gpars.seed)) < 1/3.0){ /* homozygous for 0 */
	  Pars[2*i] = Pars[2*i+1] = 2*remS;
	}
	else if(s > 2/3.0){ /* homozygous for 1 */
	  Pars[2*i] = Pars[2*i+1] = 2*remS+1;
	}
	else{ /* heterozygous */
	  Pars[2*i] = 2*remS;
	  Pars[2*i+1] = 2*remS+1;
	}
	selfed[remS]--;
	NSELF--;
      }
      else if(gpars.P == 4 && gpars.tetraType == 0){ /* auto */
	/* choose two chromosomes w/o replacement, twice */
	for(p=0; p<4; p+=2){ /* p=0,2 */
	  s = ran1(&gpars.seed);
	  if(s < 0.5){
	    Pars[4*i+p] = 4*remS;
	    if(s < 1.0/6.0)  Pars[4*i+1+p] = 4*remS+1;
	    else if(s < 2.0/6.0)  Pars[4*i+1+p] = 4*remS+2;
	    else  Pars[4*i+1+p] = 3;
	  }
	  else if(s < 5.0/6.0){
	    Pars[4*i+p] = 4*remS+1;
	    if(s < 4.0/6.0)  Pars[4*i+1+p] = 4*remS+2;
	    else  Pars[4*i+1+p] = 4*remS+3;
	  }
	  else{
	    Pars[4*i+p] = 4*remS+2;
	    Pars[4*i+1+p] = 4*remS+3;
	  }
	}
	selfed[remS]--;
	NSELF--;
      }
      else if(gpars.P == 4 && gpars.tetraType == 1){ /* allo */
	for(p=0; p<4; p+=2){ /* p=0,2 */
	  s = ran1(&gpars.seed); /* pairs of chrs are indep */
	  if(s < 0.3333333) /* both get 0 */
	    Pars[4*i+p] = Pars[4*i+1+p] = 4*remS+p;
	  else if(s > 0.6666666) /* both get 1 */
	    Pars[4*i+p] = Pars[4*i+1+p] = 4*remS+1+p;
	  else{ /* one of each */
	    Pars[4*i+p] = 4*remS+p;
	    Pars[4*i+1+p] = 4*remS+1+p;
	  }
	}
	selfed[remS]--;
	NSELF--;
      }
      else{
	fprintf(errfile,"sfs_code error(%ld):selfing must have P=2/4\n",
		INITIALSEED);
	fprintf(errfile,"if P==4, must set tetraType = 0 or 1\n");
	abort();
      }
      selfing[i] = 0;
    }
    else{ /* individual outcrossed */
      if(gpars.P == 2){
	s = ran1(&gpars.seed);
	if(ppars[pop].neutpop){
	  remFemale = (long)(s*ppars[pop].MALES);
	  s = ran1(&gpars.seed);
	  remMale=(long)(ppars[pop].N-s*(ppars[pop].N-ppars[pop].MALES));
	}
	else{
	  remFemale = invCDF(s*relFit[ppars[pop].MALES-1], relFit, 0,
			     ppars[pop].MALES-1);
	  s = ran1(&gpars.seed);
	  remMale = invCDF(1.0-(s*(1.0-relFit[ppars[pop].MALES-1])),
			   relFit, ppars[pop].MALES, ppars[pop].N-1);
	}
	/* choose even/odd parental chrm by evenness of ID numbers */
	if(gpars.SEX[0] == '1'){
	  if((i+remFemale+popn[pop].parentLoc[gpars.P*remFemale][0])%2==0)
	    Pars[2*i] = 2*remFemale;
	  else
	    Pars[2*i] = 2*remFemale + 1;
	  Pars[2*i+1] = 2*remMale + (i >= ppars[pop].MALES);
	}
	else{
	  s = ran1(&gpars.seed);
	  if(s<0.25){
	    Pars[2*i] = 2*remFemale;
	    Pars[2*i+1] = 2*remMale;
	  }
	  else if(s<0.5){
	    Pars[2*i] = 2*remFemale;
	    Pars[2*i+1] = 2*remMale+1;
	  }
	  else if(s<0.75){
	    Pars[2*i] = 2*remFemale + 1;
	    Pars[2*i+1] = 2*remMale;
	  }
	  else{
	    Pars[2*i] = 2*remFemale + 1;
	    Pars[2*i+1] = 2*remMale + 1;
	  }
	}
      }
      else{
	if(gpars.P == 1) 
	  Pars[i] = invCDF(ran1(&gpars.seed), relFit, 0, ppars[pop].N-1);
	else if(gpars.P == 4 && gpars.tetraType == 0){ /* autotetraploid */
	  s = ran1(&gpars.seed);
	  remFemale = invCDF(s*relFit[ppars[pop].MALES-1], relFit, 0,
			     ppars[pop].MALES-1);
	  s = ran1(&gpars.seed);
	  remMale = invCDF(1.0-(s*(1.0-relFit[ppars[pop].MALES-1])),
			   relFit, ppars[pop].MALES, ppars[pop].N-1);
	  for(p=0; p<4; p+=2){ /* p=0,2 */
	    rem = (p==0 ? remFemale : remMale);
	    s = ran1(&gpars.seed); /* two chromosome from one indiv. */
	    if(s < 0.5){
	      Pars[4*i+p] = 4*rem;
	      if(s < 1.0/6.0)  Pars[4*i+1+p] = 4*rem+1;
	      else if(s < 2.0/6.0)  Pars[4*i+1+p] = 4*rem+2;
	      else  Pars[4*i+1+p] = 4*rem+3;
	    }
	    else if(s < 5.0/6.0){
	      Pars[4*i+p] = 4*rem+1;
	      if(s < 4.0/6.0)  Pars[4*i+1+p] = 4*rem+2;
	      else  Pars[4*i+1+p] = 4*rem+3;
	    }
	    else{
	      Pars[4*i+p] = 4*rem+2;
	      Pars[4*i+1+p] = 4*rem+3;
	    }
	  }
	} /* end autotetraploid case */
	else if(gpars.P == 4 && gpars.tetraType == 1){ /* allotetraploid */
	  s = ran1(&gpars.seed);
	  remFemale = invCDF(s*relFit[ppars[pop].MALES-1], relFit, 0,
			     ppars[pop].MALES-1);
	  s = ran1(&gpars.seed);
	  remMale = invCDF(1.0-(s*(1.0-relFit[ppars[pop].MALES-1])),
			   relFit, ppars[pop].MALES, ppars[pop].N-1);
	  for(p=0; p<2; p++){
	    rem = (p==1 ? remFemale : remMale);
	    s=ran1(&gpars.seed);
	    if(s < 0.25){ /* get first of each pair */
	      Pars[4*i+p] = 4*rem;
	      Pars[4*i+2+p] = 4*rem+2;
	    }
	    else if(s < 0.5){ /* get second of each pair */
	      Pars[4*i+p] = 4*rem+1;
	      Pars[4*i+2+p] = 4*rem+3;
	    }
	    else if (s < 0.75){ /* 0::1 */
	      Pars[4*i+p] = 4*rem;
	      Pars[4*i+2+p] = 4*rem+3;
	    }
	    else{ /* 1::0 */
	      Pars[4*i+p] = 4*rem+1;
	      Pars[4*i+2+p] = 4*rem+2;
	    }
	  }
	} /* end allotetraploid case */
	else{
	  fprintf(errfile,"sfs_code error(%ld):  P must be 1,2, or 4.\n",
		  INITIALSEED);
	  fprintf(errfile,"if P==4, tetraType must be set to 0, 1\n");
	  abort();
	}
      }
    } /* end non-selfing case */
    
    /* update allele counts and perform inter-locus recombination */
    for(p=0; p<gpars.P; p++){
      for(rem=0; rem<gpars.R; rem++){
	/* store parentLoc for future reference */
	parStore[gpars.P*i+p][rem] = popn[pop].parentLoc[gpars.P*i+p][rem];
	
	/* first swap parental chrms for inter-locus recombination */
	if(rem > 0 && between != NULL){
	  while(iB<numXB && between[iB].id==gpars.P*i+p &&
		iB2<between[iB].numEvents && between[iB].loci[iB2]==rem-1){
	    iB2++;
	    if(gpars.P==2)
	      Pars[gpars.P*i+p] += (Pars[gpars.P*i+p]%2==0 ? 1 : -1);
	    else if(gpars.P==4 && gpars.tetraType==0){
	      /* autotetraploid: can choose any other chromosome;*/
	      r = Pars[gpars.P*i+p];
	      Pars[gpars.P*i+p] += (int)(ran1(&gpars.seed)*4);
	      if(Pars[gpars.P*i+p]%4 < r%4)/*only if spill over to next ind*/
		Pars[gpars.P*i+p] -= 4;
	    }
	    else if(gpars.P==4 && gpars.tetraType==0){
	      /* allotetraploid: choose (0,1) or (2,3) */
	      r = Pars[gpars.P*i+p];
	      Pars[gpars.P*i+p] += (int)(ran1(&gpars.seed)*2);
	      if(Pars[gpars.P*i+p]%2 < r%2)/*only if spill over to next ind*/
		Pars[gpars.P*i+p] -= 2;
	    }
	    if(iB2>=between[iB].numEvents){
	      free(between[iB].loci);
	      between[iB].loci = NULL;
	      iB++;
	      iB2 = 0;
	      break;
	    }
	  }
	}
	else if(rem>0 && gpars.recMap==NULL &&
		fabs(gpars.LINK[rem-1])>FLT_EPSILON &&
		(i<ppars[pop].MALES || gpars.SEX[rem] == '0')){
	  s = 0.0; /* swap probability */
	  if(((gpars.LinkDistType=='g' && gpars.LINK[rem-1]>=0.5-FLT_EPSILON) 
	      ||
	      (gpars.LinkDistType=='p' && gpars.LINK[rem-1] < -FLT_EPSILON)))
	    s = 0.5; /* independence */
	  else if(gpars.LinkDistType=='g' && gpars.LINK[rem-1] > DBL_EPSILON)
	    s = gpars.LINK[rem-1]; /* genetic distance, use directly */
	  else if(gpars.LinkDistType=='p' && gpars.LINK[rem-1] > 0)
	    s = gpars.LINK[rem-1]*ppars[pop].RHO/(2.0*ppars[pop].N*gpars.P);
	  s = (s>0.5 ? 0.5 : s); /* 0.5 is the maximal value */
	  
	  if(s>DBL_EPSILON){
	    n=0;
	    /* think of case pMaleRec=0.5, then the prob of recombination 
	       in either sex is s = 2*s*pMaleRec */
	    if(p < gpars.P/2.0){
	      if(ran1(&gpars.seed) < 2*s*(1-ppars[pop].pMaleRec))
		n = 1;
	    }
	    else{
	      if(ran1(&gpars.seed) < 2*s*(ppars[pop].pMaleRec))
		n = 1;
	    }
	    if(n==1){ /* this chromosome is recombinant */
	      if(gpars.P==2)
		Pars[gpars.P*i+p] += (Pars[gpars.P*i+p]%2==0 ? 1 : -1);
	      else if(gpars.P==4 && gpars.tetraType==0){
		/* autotetraploid: can choose any other chromosome;*/
		r = Pars[gpars.P*i+p];
		Pars[gpars.P*i+p] += (int)(ran1(&gpars.seed)*4);
		if(Pars[gpars.P*i+p]%4 < r%4)/*only if goes to next ind*/
		  Pars[gpars.P*i+p] -= 4;
	      }
	      else if(gpars.P==4 && gpars.tetraType==0){
		/* allotetraploid: choose (0,1) or (2,3) */
		r = Pars[gpars.P*i+p];
		Pars[gpars.P*i+p] += (int)(ran1(&gpars.seed)*2);
		if(Pars[gpars.P*i+p]%2 < r%2)/*only if goes to next ind*/
		  Pars[gpars.P*i+p] -= 2;
	      }
	    }
	  }
	}
	
	/* perform recombination/Gene Conversion within locus */
	if(within != NULL && iW<numXW && within[iW].id==gpars.P*i+p &&
	   iW2<within[iW].numEvents && within[iW].loci[iW2]==rem){
	  getRecombinant(gpars.P*i+p, pop, Pars[gpars.P*i+p],
			 &Pars[gpars.P*i+p], within, iW, &iW2);
	}
	else{
	  if(gpars.P*i+p >= Pars[gpars.P*i+p] &&
	     popn[pop].parentLoc[gpars.P*i+p][rem] !=
	     parStore[Pars[gpars.P*i+p]][rem]){
	    popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][rem]][rem]--;
	    if(popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][rem]][rem]
	       == 0){
	      popn[pop].chrGen[popn[pop].parentLoc[gpars.P*i+p][rem]][rem]
		= popn[pop].gen;
	      newDead(popn[pop].parentLoc[gpars.P*i+p][rem], rem);
	    }
	    popn[pop].numCopies[parStore[Pars[gpars.P*i+p]][rem]][rem]++;
	    popn[pop].parentLoc[gpars.P*i+p][rem] =
	      parStore[Pars[gpars.P*i+p]][rem];
	  }
	  else if(gpars.P*i+p < Pars[gpars.P*i+p] &&
		  popn[pop].parentLoc[gpars.P*i+p][rem] !=
		  popn[pop].parentLoc[Pars[gpars.P*i+p]][rem]){
	    popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][rem]][rem]--;
	    if(popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][rem]][rem]
	       == 0){
	      popn[pop].chrGen[popn[pop].parentLoc[gpars.P*i+p][rem]][rem] =
		popn[pop].gen;
	      newDead(popn[pop].parentLoc[gpars.P*i+p][rem], rem);
	    }
	    popn[pop].numCopies[popn[pop].parentLoc[Pars[gpars.P*i+p]]
				[rem]][rem]++;
	    popn[pop].parentLoc[gpars.P*i+p][rem] =
	      popn[pop].parentLoc[Pars[gpars.P*i+p]][rem];
	  }
	}
      }
      if(within != NULL && iW<numXW && within[iW].id==gpars.P*i+p){
#ifdef UNIT_TEST
	if(iW2 != within[iW].numEvents){
	  fprintf(errfile,"not all recomb events found in selecting()\n");
	  fprintf(errfile,"ind=%ld;  iW=%ld/%ld; seed=%ld\n",gpars.P*i+p, iW,
		  numXW, INITIALSEED);
	  for(i=0; i<within[iW].numEvents; i++)
	    fprintf(errfile,"id=%ld: %ld.%ld\n",within[iW].id, 
		    within[iW].loci[i], within[iW].pos[i]);
	  fflush(errfile);
	  abort();
	}
#endif
	free(within[iW].loci);
	free(within[iW].pos);
	free(within[iW].type);
	within[iW].loci = NULL;
	within[iW].pos = NULL;
	within[iW].type = NULL;
	iW++;
	iW2 = 0;
      }
    } /* end iteration over ploidy */
  } /* end iteration of individuals */
  
#ifdef UNIT_TEST
  if(between != NULL){
    for(i=0; i<numXB; i++){
      if(between[i].loci != NULL){
	fprintf(errfile,"error: between[%ld].id=%ld (out of %ld) not freed!\n",
		i,between[i].id, numXB);
	fprintf(errfile,"events:\n");
	for(iW2=0; iW2<between[i].numEvents; iW2++){
	  fprintf(errfile,"%ld\n",between[i].loci[iW2]);
	}
	fflush(errfile);
	abort();
      }
    }
  }

  if(within != NULL){
    for(i=0; i<numXW; i++){
      if(within[i].loci != NULL || within[i].pos != NULL || 
	 within[i].type != NULL){
	fprintf(errfile,"error: within[%ld].id=%ld (out of %ld) not freed!\n",i,
		within[i].id, numXW);
	fprintf(errfile,"events:\n");
	for(iW2=0; iW2<within[i].numEvents; iW2++){
	  fprintf(errfile,"%ld.%ld\n",within[i].loci[iW2],within[i].pos[iW2]);
	}
	fflush(errfile);
	abort();
      }
    }
  }
#endif
  if(within != NULL)
    free(within);
  if(between != NULL)
    free(between);
#ifdef VERBOSE_DEBUG
  {
    printf("parents:\n");
    for(i=0; i<gpars.P*ppars[pop].N; i++){
      printf("%ld:%ld[--%ld(%ld)--] ",i,Pars[i], popn[pop].parentLoc[i][0],
	     popn[pop].numCopies[popn[pop].parentLoc[i][0]][0]);
    }
    printf("\n");
    fflush(stdout);
    printf("PL:\n");
    for(i=0; i<gpars.P*ppars[pop].N; i++){
      printf("%ld ",popn[pop].parentLoc[i][0]);
    }
    printf("\n");
    fflush(stdout);
  }
#endif
}


/* ------------------------------------------------------------------------- */

void distribRec(const int pop, struct locRec **within, long *numX1,
		struct interRec **between, long *numX2)
{
  int i1, i2;
  long i, k, ind, loc, pos;
  long numX=0, x;
  float rn;

  if(ppars[pop].RHO<FLT_EPSILON && ppars[pop].fGC<FLT_EPSILON)
    return;
  
  /* 2 cases:  uniform rate across loci, or use recombination map */
  if(gpars.recMap == NULL){ /* uniform recomb rate */
    if(ppars[pop].RHO < FLT_EPSILON) /* GC only */
      numX = (long)poidev(ppars[pop].fGC*gpars.totalNucs/2.0, &gpars.seed);
    else
      numX = (long)poidev(ppars[pop].RHO*gpars.totalNucs/2.0, &gpars.seed);
  }
  else{ /* use recombination map */
    if(ppars[pop].RHO < FLT_EPSILON) /* GC only */
      numX = (long)poidev(ppars[pop].fGC*gpars.sumLL[gpars.R-1]/2.0,
			  &gpars.seed);
    else
      numX = (long)poidev(ppars[pop].RHO*gpars.sumLL[gpars.R-1]/2.0,
			  &gpars.seed);
  }

#ifdef UNIT_TEST
  fprintf(errfile,"numX = %ld\n",numX);
  fflush(errfile);
#endif

  if(numX == 0)
    return; /* nothing to do */
  
  else if(numX < 0.1*ppars[pop].N){
/*     printf("few: ");fflush(stdout); */
    /* few recombinants: use quicker method */
    /* split population into numX bins, and choose 1 individual randomly
       from each bin. */
    assert(*within = calloc(numX, sizeof(**within)));/*final size mostly known*/
    for(i=0; i<numX; i++){
      /* draw individual uniformly from group (1st chr of indiv) */
      ind = (long)((i+ran1(&gpars.seed))*floor(ppars[pop].N/(numX+.0)))*gpars.P;
      
      /* now get locus and position */
      if(gpars.recMap == NULL){ /* uniform recomb rate */
	/* position just a linear transformation of a random number */
	pos = (long)(ran1(&gpars.seed)*(gpars.totalNucs-1));
	if(gpars.R == 0)
	  loc = 0;
	else{
	  loc = invCDF(pos, gpars.sumL, 0, gpars.R-1);
	  if(loc>0)
	    pos -= gpars.sumL[loc-1];
	  if(pos == gpars.L[loc]){
	    pos = 0;
	    loc++;
	  }
	}
	if(gpars.SEX[loc] == '1' && ind >= ppars[pop].MALES*gpars.P)
	  continue; /* males do not recomb on sex chrms, nothing to do */
      }
      else{ /* use recombination map */
	rn = ran1(&gpars.seed);
	k = invCDF(rn, gpars.recMap, 0, gpars.mapPoints-1);
	
	/* determine site by inverse transformation of rn */
	if(k == 0)
	  pos = (long)(rn/gpars.recMap[0]*gpars.recMapPos[0]);
	else
	  pos = (long)(gpars.recMapPos[k-1] + (rn-gpars.recMap[k-1])/
		       (gpars.recMap[k]-gpars.recMap[k-1])*
		       (gpars.recMapPos[k]-gpars.recMapPos[k-1]));
	if(gpars.R == 0)
	  loc = 0;
	else{
	  loc = invCDF(pos, gpars.sumLL, 0, gpars.R-1);
	  if(loc > 0)
	    pos -= gpars.sumLL[loc-1]; /* correct for locus boundaries */
	}
	pos--;

	if(gpars.SEX[loc] == '1' && ind >= ppars[pop].MALES*gpars.P)
	  continue; /* males do not recomb on sex chrms, nothing to do */
	
	if(pos >= gpars.L[loc]){ /* interlocus recombination! */
	  /* check to make sure this is a recombination event! */
	  if(ppars[pop].fGC<FLT_EPSILON || ppars[pop].fGC>1-FLT_EPSILON ||
	     ran1(&gpars.seed)<ppars[pop].fGC){
	    if(*between == NULL)
	      assert(*between = calloc(1, sizeof(**between)));
	    else
	      assert(*between=realloc(*between,((*numX2)+1)*sizeof(**between)));
	    (*between)[*numX2].id = ind;
	    (*between)[*numX2].numEvents = 1;
	    assert((*between)[*numX2].loci = malloc(sizeof(long)));
	    (*between)[*numX2].loci[0] = loc;
	    (*numX2)++;
	  }
	  continue;
	}
      }
      
      /* if we've made it this far, then there is a recombination event within
	 a locus, so set values of structure within */
      (*within)[*numX1].id = ind;
      (*within)[*numX1].numEvents = 1;
      assert((*within)[*numX1].loci = malloc(sizeof(long)));
      (*within)[*numX1].loci[0] = loc;
      assert((*within)[*numX1].pos = malloc(sizeof(long)));
      (*within)[*numX1].pos[0] = pos;
      assert((*within)[*numX1].type = malloc(sizeof(int)));

      /* pick which chromosome receives event with pMaleRec bias */
      if(gpars.P==2){
	if(ran1(&gpars.seed) < ppars[pop].pMaleRec)
	  (*within)[*numX1].id++; /* paternal recombination */
      }
      else if(gpars.P==4){
	if(ran1(&gpars.seed) < ppars[pop].pMaleRec){/*choose 1 of 2 male chrms*/
	  if(ran1(&gpars.seed) < 0.5) 
	    (*within)[*numX1].id += 3;
	  else
	    (*within)[*numX1].id += 2;
	}
	else{ /* choose 1 of 2 female chrms */
	  if(ran1(&gpars.seed) < 0.5)
	    (*within)[*numX1].id += 1;
	}
      }
      else{
	fprintf(errfile,"error, for recombination must have ploidy==2 or 4\n");
	abort();
      }

      /* determine whether each event is GC, crossover, or both */
      i1=i2=0;
      if(ppars[pop].fGC<-0.5) /* only modeling cross-overs */
	i2 = 1;
      else if(ppars[pop].RHO<FLT_EPSILON) /* only modeling gene conversion */
	i1 = 1;
      else if(ppars[pop].fGC<=1+FLT_EPSILON){
	i1 = 1; /* definitely a gene conversion event */
	if(ran1(&gpars.seed)<ppars[pop].fGC)
	  i2 = 1; /* also a cross-over */
      }
      else if(ppars[pop].fGC>1+FLT_EPSILON){
	i2 = 1; /* definitely a cross-over event */
	if(ran1(&gpars.seed)<1.0/ppars[pop].fGC)
	  i1 = 1; /* also gene conversion */
      }
      else /* GC & CO at equal rates */
	i1 = i2 = 1;
      
      if(i1 == 0) /* crossover with no GC */
	(*within)[*numX1].type[0] = 0;
      else if(i2 == 0) /* GC without crossover */
	(*within)[*numX1].type[0] = 1;
      else /* GC with crossover */
	(*within)[*numX1].type[0] = 2;
      (*numX1)++;
    }
/*     printf("\n");fflush(stdout); */
  }
  else{
/*     printf("many: ");fflush(stdout); */
    /* many recombinants: use slower method */
    /* for each chromosome, draw a Binomial number of GC/crossovers
       until there are no more events left */
    ind = -1; /* loop over all chromosomes */
    /* final size of within unknown, generate dynamically */
    assert(*within = calloc(*numX1+1, sizeof(**within)));
    (*within)[*numX1].loci = NULL;
    (*within)[*numX1].pos = NULL;
    (*within)[*numX1].type = NULL;
    
    /* only need between if using recombination map with multiple loci */
    if(gpars.R>1 && gpars.recMap!=NULL){ 
      /* final size of between unknown, generate dynamically */
      assert(*between = malloc(sizeof(**between)));
      (*between)[0].loci = NULL;
      (*between)[0].numEvents = 0;
    }
    while(numX > 0){
      ind++;
      if(ind%gpars.P<gpars.P/2.0) /* descends from mother */
	x = (long)(bnldev(2*(1-ppars[pop].pMaleRec)/(gpars.P*ppars[pop].N-ind),
			  numX, &gpars.seed));
      else /* descends from father */
	x = (long)(bnldev(2*ppars[pop].pMaleRec/(gpars.P*ppars[pop].N-ind),
			  numX, &gpars.seed));
      /* note that if Nm=Nf=P*N/2 (num male chrms=num female chrms) then
	 gpars.P*Nm/2*2*pMaleRec/gpars.P/N+gpars.P*Nf/2*2*(1-pMaleRec)/P/N = 1! 
	 so the probability of randomly choosing a maternal/paternal chromosome
	 is as written above */
      if(x == 0)
	continue;
      else
	numX -= x;

/*       printf("ind=%ld; x=%ld; numX1=%ld; numX2=%ld\t",ind,x,*numX1,*numX2);fflush(stdout); */

      (*within)[*numX1].id = ind;
      (*within)[*numX1].numEvents = 0;
      if((*within)[*numX1].loci==NULL){
	assert((*within)[*numX1].loci = malloc(x*sizeof(long)));
	assert((*within)[*numX1].pos = malloc(x*sizeof(long)));
	assert((*within)[*numX1].type = malloc(x*sizeof(int)));
      }
      else{
	assert((*within)[*numX1].loci = realloc((*within)[*numX1].loci,
						x*sizeof(long)));
	assert((*within)[*numX1].pos = realloc((*within)[*numX1].pos,
					       x*sizeof(long)));
	assert((*within)[*numX1].type = realloc((*within)[*numX1].type,
						x*sizeof(int)));
      }
      for(i=0; i<x; i++){
	/* get locus & position */
	if(gpars.recMap == NULL){ /* uniform recomb rate */
	  pos = ran1(&gpars.seed)*(gpars.sumL[gpars.R-1]-1);
	  if(gpars.R == 0)
	    loc = 0;
	  else
	    loc = invCDF(pos, gpars.sumL, 0, gpars.R-1);
	  if(gpars.SEX[loc] == '1' && ind >= ppars[pop].MALES*gpars.P){
	    /* males do not recombine on sex chromosomes, nothing to do */
	    continue;
	  }
	}
	else{ /* use recombination map */
	  rn = ran1(&gpars.seed);
	  k = invCDF(rn, gpars.recMap, 0, gpars.mapPoints-1);
	  
	  /* determine site by inverse transformation of rn */
	  if(k == 0)
	    pos = (long)(rn/gpars.recMap[0]*gpars.recMapPos[0]);
	  else
	    pos = (long)(gpars.recMapPos[k-1] + (rn-gpars.recMap[k-1])/
			 (gpars.recMap[k]-gpars.recMap[k-1])*
			 (gpars.recMapPos[k]-gpars.recMapPos[k-1]));
	  if(gpars.R==0)
	    loc = 0;
	  else
	    loc = invCDF(pos, gpars.sumLL, 0, gpars.R-1);
	  if(gpars.SEX[loc] == '1' && ind >= ppars[pop].MALES*gpars.P){
	    /* males do not recombine on sex chromosomes, nothing to do */
	    continue;
	  }

	  if((loc==0 && pos>gpars.L[0]) || 
	     (loc>0 && pos>gpars.sumLL[loc-1]+gpars.L[loc])){
	    /* interlocus recombination! */
	    if(ppars[pop].fGC<FLT_EPSILON || ppars[pop].fGC>1-FLT_EPSILON ||
	       ran1(&gpars.seed)<ppars[pop].fGC){
	      /* make sure it is a cross-over event */
	      if((*between)[*numX2].loci == NULL)
		assert((*between)[*numX2].loci = malloc(sizeof(long)));
	      else
		assert((*between)[*numX2].loci = 
		       realloc((*between)[*numX2].loci, 
			       ((*between)[*numX2].numEvents+1)*sizeof(long)));
	      (*between)[*numX2].id = ind;
	      (*between)[*numX2].loci[(*between)[*numX2].numEvents] = loc;
	      (*between)[*numX2].numEvents++;
	    }
/* 	    printf("BnumX2=%ld(%ld)\n",*numX2,(*between)[*numX2].numEvents);fflush(stdout); */
	    continue;
	  }
	}
	/* if we've made it this far, then there is a recombination event,
	   so set values of structure within */
	(*within)[*numX1].loci[(*within)[*numX1].numEvents] = loc;
	(*within)[*numX1].pos[(*within)[*numX1].numEvents] = pos;
	/* determine whether each event is GC, crossover, or both */
	i1=i2=0;
	if(ppars[pop].fGC<-0.5) /* only modeling cross-overs */
	  i2 = 1;
	else if(ppars[pop].RHO<FLT_EPSILON) /* only modeling gene conversion */
	  i1 = 1;
	else if(ppars[pop].fGC<1-FLT_EPSILON){
	  i1 = 1; /* definitely a gene conversion event */
	  if(ran1(&gpars.seed)<ppars[pop].fGC)
	    i2 = 1; /* also a cross-over */
	}
	else if(ppars[pop].fGC>1+FLT_EPSILON){
	  i2 = 1; /* definitely a cross-over event */
	  if(ran1(&gpars.seed)<1.0/ppars[pop].fGC)
	    i1 = 1; /* also gene conversion */
	}
	else /* GC & CO at equal rates */
	  i1 = i2 = 1;
	
	if(i1 == 0) /* crossover with no GC */
	  (*within)[*numX1].type[(*within)[*numX1].numEvents] = 0;
	else if(i2 == 0) /* GC without crossover */
	  (*within)[*numX1].type[(*within)[*numX1].numEvents] = 1;
	else /* GC with crossover */
	  (*within)[*numX1].type[(*within)[*numX1].numEvents] = 2;
	(*within)[*numX1].numEvents++;
      }
      if((*within)[*numX1].numEvents>0){
	if((*within)[*numX1].numEvents>1){
	  sort2long((unsigned long)((*within)[*numX1].numEvents-1),
		    &(*within)[*numX1].pos[0],
		    &(*within)[*numX1].loci[0]);
	}
	for(i=0; i<(*within)[*numX1].numEvents; i++){
	  if((*within)[*numX1].loci[i] > 0){
	    if(gpars.recMap != NULL)
	      (*within)[*numX1].pos[i] -=
		gpars.sumLL[(*within)[*numX1].loci[i]-1];
	    else
	      (*within)[*numX1].pos[i] -=
		gpars.sumL[(*within)[*numX1].loci[i]-1];
	  }
/* 	  printf("L=%ld; P=%ld;\t",(*within)[*numX1].loci[i],(*within)[*numX1].pos[i]);fflush(stdout); */
	}
	(*numX1)++;
	assert((*within) = realloc(*within, (*numX1+1)*sizeof(**within)));
	(*within)[*numX1].loci = NULL;
	(*within)[*numX1].pos = NULL;
	(*within)[*numX1].type = NULL;
	(*within)[*numX1].numEvents = 0;
      }
      if(gpars.R>1 && gpars.recMap!=NULL){
	if((*between)[*numX2].numEvents > 0){
	  if((*between)[*numX2].numEvents > 1){
	    sortlong((unsigned long)((*between)[*numX2].numEvents-1),
		     &(*between)[*numX2].loci[0]);
	  }
	  (*numX2)++;
	  assert((*between) = realloc(*between, (*numX2+1)*sizeof(**between)));
	  (*between)[*numX2].loci = NULL;
	  (*between)[*numX2].numEvents = 0;
	}
      }
      /*       printf("\n");fflush(stdout); */
    }
    if(*numX1==0){
      if((*within)[*numX1].loci != NULL){
	free((*within)[*numX1].loci);
	free((*within)[*numX1].pos);
	free((*within)[*numX1].type);
      }
      free(*within);
      (*within) = NULL;
    }
    else if((*within)[*numX1].numEvents==0){
      if((*within)[*numX1].loci != NULL){
	free((*within)[*numX1].loci);
	free((*within)[*numX1].pos);
	free((*within)[*numX1].type);
      }
      assert((*within) = realloc(*within, (*numX1)*sizeof(**within)));
    }
    
    if(gpars.R>1 && gpars.recMap!=NULL){
      if(*numX2==0){
	free(*between);
	(*between) = NULL;
      }
      else if((*between)[*numX2].numEvents==0){
	assert((*between) = realloc(*between, (*numX2)*sizeof(**between)));
/* 	(*numX2)--; */
      }
    }
  }
#ifdef UNIT_TEST
  for(i=0; i<*numX1; i++){
    long nloc = -1;
    loc = -1;
    if(i>0){
      if((*within)[i].id <= (*within)[i-1].id){
	fprintf(errfile,"error in distribRec() %ld: individuals not sorted!\n",
		INITIALSEED);
	fprintf(errfile,"within[%ld].id=%ld;  within[%ld].id=%ld! tot=%ld\n",
		i-1,(*within)[i-1].id,i,(*within)[i].id, *numX1);
	abort();
      }
    }
    for(k=0; k<(*within)[i].numEvents; k++){
      if((*within)[i].loci[k] == loc){
	if((*within)[i].pos[k] < nloc){
	  fprintf(errfile,"error in distribRec() %ld: positions not sorted\n",
		  INITIALSEED);
	  abort();
	}
	loc = (*within)[i].pos[k];
      }
      else{
	if((*within)[i].loci[k] == loc){
	  fprintf(errfile,"error in distribRec() %ld: loci not sorted\n",
		  INITIALSEED);
	  abort();
	}
	loc = (*within)[i].loci[k];
	nloc = (*within)[i].pos[k];
      }
    }
  }
#endif
  return;
}

/* ------------------------------------------------------------------------- */

void getRecombinant(const long ind, const int pop, long oldPar, long *newPar,
		    struct locRec *within, const long iW, long *iW2)
{
  int TR1=0, TR2=0, D, F, j;/* tract len, donor strand{0,1}, current par chrm */
  long x[2], xPL[2], txPL[2], PL, tPL;
  long i, k, min, max, numx=0, bTR, bTR1, first=-1, last=0;
  long loc=within[iW].loci[*iW2];
  float rn;
  double delta;
  struct history *ptr=NULL;
  
  /* choose sister-chromosome to swap with */
  PL = popn[pop].parentLoc[ind][loc];
  x[0] = oldPar;
  if(gpars.P == 2) /* obvious in diploid case */
    x[1] = x[0] + (x[0]%2==0 ? 1 : -1);
  else if(gpars.P == 4 && gpars.tetraType == 0){ /* autotetraploid */
    rn = ran1(&gpars.seed);
    if(rn < 0.3333)
      x[1] = x[0]+1;
    else if(rn < 0.6667)
      x[1] = x[0]+2;
    else
      x[1] = x[0]+3;
    if(x[0]%2 > x[1]%2)
      x[1] -= 4;
  }
  else
    x[1] = x[0] + (x[0]%2==0 ? 1 : -1);

  if(ind >= x[0])
    xPL[0] = parStore[x[0]][loc];
  else
    xPL[0] = popn[pop].parentLoc[x[0]][loc];
  if(ind >= x[1])
    xPL[1] = parStore[x[1]][loc];
  else
    xPL[1] = popn[pop].parentLoc[x[1]][loc];

#ifdef VERBOSE_DEBUG
  printf("parent 0 (%ld.%ld:%ld)\n", x[0], xPL[0],
	 popn[pop].numCopies[xPL[0]][loc]);
  k=0;
  PrintHistTree(popn[pop].BigHead[xPL[0]][loc]->Rtree, &k, pop);
  printf("\n\n");
  fflush(stdout);
  printf("parent 1 (%ld.%ld:%ld)\n", x[1], xPL[1],
	 popn[pop].numCopies[xPL[1]][loc]);
  k=0;
  PrintHistTree(popn[pop].BigHead[xPL[1]][loc]->Rtree, &k, pop);
  printf("\n\n");
  fflush(stdout);
#endif

  /* check to make sure that we actually need to do something...  
     if there is no biased gene conversion, it suffices to check if the parental
     chromosomes point to the same location or if the haplotypes are empty
     or only carry the same mutation.  If there is gene conversion, need to
     loop over all events to make sure none spill onto the previous locus */
  j=0;
  if(!gpars.TRACKANC && fabs(ppars[pop].BGC-0.5)<FLT_EPSILON &&
     (gpars.R==1 || ppars[pop].fGC<FLT_EPSILON) &&
     (xPL[0] == xPL[1] || (popn[pop].BigHead[xPL[0]][loc]->Rtree == NULL &&
			   popn[pop].BigHead[xPL[1]][loc]->Rtree == NULL) ||
      (popn[pop].extremeMuts[xPL[0]][loc][0] ==
       popn[pop].extremeMuts[xPL[0]][loc][1] &&
       popn[pop].extremeMuts[xPL[0]][loc][0] ==
       popn[pop].extremeMuts[xPL[1]][loc][0] &&
       popn[pop].extremeMuts[xPL[0]][loc][0] ==
       popn[pop].extremeMuts[xPL[1]][loc][1]))){
    while((*iW2)<within[iW].numEvents && within[iW].loci[*iW2]==loc){
      if(within[iW].type[*iW2]%2 == 0) /* recomb; swap parental chromosome */
	j = (j+1)%2;
      (*iW2)++;
    }
    (*newPar) = x[j];
    popn[pop].numCopies[PL][loc]--;
    if(popn[pop].numCopies[PL][loc]==0){
      popn[pop].chrGen[PL][loc] = popn[pop].gen;
      newDead(PL, loc);
    }
    popn[pop].numCopies[xPL[j]][loc]++;
    popn[pop].parentLoc[ind][loc] = xPL[j];
    return;
  }
  
#ifdef UNIT_TEST
  long numMuts1=0, numMuts2=0;
  if(gpars.R>1 || within[iW].numEvents!=1 || within[iW].type[0]!=0){
    printf("not checking for correct mut copying in getRecombinant()...\n");
    printf("gpars.R=%ld; numEvents=%ld; type[0]=%d\n",gpars.R,
	   within[iW].numEvents, within[iW].type[0]);
    fflush(stdout);
  }
  else{
    countMuts(popn[pop].BigHead[xPL[0]][0]->Rtree, &numMuts1, 0,
	      within[iW].pos[0]);
    countMuts(popn[pop].BigHead[xPL[1]][0]->Rtree, &numMuts1, 
	      within[iW].pos[0]+1, gpars.L[0]-1);
  }
#endif
#ifdef VERBOSE_DEBUG
  printf("par0=%ld(%ld:%ld); par1=%ld(%ld:%ld)\n",x[0],xPL[0],
	 popn[pop].numCopies[xPL[0]][loc], x[1],xPL[1],
	 popn[pop].numCopies[xPL[1]][loc]);
  printf("tree 0 [%ld, %ld]:\n", popn[pop].extremeMuts[xPL[0]][loc][0],
	 popn[pop].extremeMuts[xPL[0]][loc][1]);
  k=0;
  PrintHistTree(popn[pop].BigHead[xPL[0]][loc]->Rtree, &k, pop);
  printf("tree 1 [%ld, %ld]:\n", popn[pop].extremeMuts[xPL[1]][loc][0],
	 popn[pop].extremeMuts[xPL[1]][loc][1]);
  k=0;
  PrintHistTree(popn[pop].BigHead[xPL[1]][loc]->Rtree, &k, pop);
#endif
  
  /* create blank haplotype, and copy necessary mutations to it */
#ifdef VERBOSE_DEBUG
  printf("offsp = %ld (%ld ->",ind,PL);
#endif
  
  createNewHistory(ind, loc, pop, 0, 0);
  PL = popn[pop].parentLoc[ind][loc];
  if(gpars.substMod == 5)
    memcpy(popn[pop].hpSUM[PL][loc], popn[pop].hpSUM[xPL[0]][loc],
	   gpars.L[loc]*sizeof(double));

#ifdef VERBOSE_DEBUG
  printf("%ld)\n",PL);
  fflush(stdout);
#endif

  max = -1;
  i = *iW2;
  j = 0;
  numx = 0;
  while((i<within[iW].numEvents && within[iW].loci[i]==loc) || 
	(i<=within[iW].numEvents && within[iW].loci[i-1]==loc)){
    numx++;
    min = max+1;
    
    if(i<within[iW].numEvents && within[iW].loci[i]==loc){
      if(first < 0)
	first = within[iW].pos[i];
      last = within[iW].pos[i];

      if(within[iW].type[i]>0){ /* gene conversion event */
	F = (ran1(&gpars.seed)<0.5 ? 0 : 1); /* convert before copy */
	TR1 = geomdev(1/ppars[pop].GCtract, &gpars.seed);
	TR2 = geomdev(1/ppars[pop].GCtract, &gpars.seed);
	bTR = within[iW].pos[i]-TR1-TR2+1;
	if(first<0 || bTR<first)
	  first = (bTR<0 ? 0:bTR);
	D = (j+1)%2; /* DSB always on inherited chromosome */

#ifdef VERBOSE_DEBUG
	printf("gene conversion event:TR(%ld)=%ld-%ld(donor=%d,j=%d,type=%d)\n",
	       within[iW].loci[i], bTR, within[iW].pos[i], D, j,
	       within[iW].type[i]);
	fflush(stdout);
#endif
	
	/* note that a conversion tract can overlap a previous event/locus */
	if(bTR > min){/*GC does not include previous locus or event */
	  /* copy [min,bTR) */
	  copyPartialHistory(&popn[pop].BigHead[PL][loc]->Rtree,
			     &popn[pop].BigHead[xPL[j]][loc]->Rtree,
			     min, bTR-1, pop, PL, loc);
	  if(gpars.TRACKANC){
	    memcpy(&(popn[pop].ancestry[PL][loc][min]),
		   &(popn[pop].ancestry[xPL[j]][loc][min]),
		   (bTR-min)*sizeof(int));
	  }
	}
	else{
	  if(min>0){ /* overlaps previous event */
	    /* reconvert [bTR/0,min] */
#ifdef VERBOSE_DEBUG
	    printf("before freeHistoryNodeRange(%ld,%ld):\n",bTR, min-1);
	    k=0;
	    PrintHistTree(popn[pop].BigHead[PL][loc]->Rtree, &k,pop);
#endif
	    popn[pop].BigHead[PL][loc]->Rtree =
	      freeHistoryNodeRange((bTR<0?0:bTR), min,
				   popn[pop].BigHead[PL][loc]->Rtree, loc, pop,
				   PL);
#ifdef VERBOSE_DEBUG
	    printf("after freeHistoryNodeRange:\n");
	    k=0;
	    PrintHistTree(popn[pop].BigHead[PL][loc]->Rtree, &k,pop);
	    fflush(stdout);
#endif
	  }
	  
	  if(loc>0 && bTR<0 && gpars.LinkDistType=='p' && gpars.LINK[loc-1]>=0&&
	     bTR+gpars.LINK[loc-1]<0){ /* overlaps previous locus! */
	    if(ind >= x[0])
	      txPL[0] = parStore[x[0]][loc-1];
	    else
	      txPL[0] = popn[pop].parentLoc[x[0]][loc-1];
	    if(ind >= x[1])
	      txPL[1] = parStore[x[1]][loc-1];
	    else
	      txPL[1] = popn[pop].parentLoc[x[1]][loc-1];
	    
#ifdef VERBOSE_DEBUG
	    printf("x[0]=%ld; x[1]=%ld; ",x[0], x[1]);
	    fflush(stdout);
	    printf("parStore[0][%ld]=%ld; ", loc-1, txPL[0]);
	    fflush(stdout);
	    printf("parStore[1][%ld]=%ld\n", loc-1, txPL[1]);
	    fflush(stdout);
#endif
	    
	    if(txPL[0] != txPL[1] &&
	       (popn[pop].BigHead[txPL[0]][loc-1]->Rtree!=NULL ||
		popn[pop].BigHead[txPL[1]][loc-1]->Rtree!=NULL)){
	      /*overlaps prev locus: create new history there, delete events 
		in range, then add back in converted tract */
	      bTR1 = gpars.L[loc-1]+(within[iW].pos[i]+gpars.LINK[loc-1]-
				     TR1-TR2);
	      tPL = popn[pop].parentLoc[ind][loc-1]; /*BEFORE CREATING NEW!*/
	      
#ifdef VERBOSE_DEBUG
	      k=0;
	      printf("par0(%ld) [%ld, %ld]; loc=%ld\n",txPL[0],
		     popn[pop].extremeMuts[txPL[0]][loc-1][0],
		     popn[pop].extremeMuts[txPL[0]][loc-1][1],
		     loc-1);
	      PrintHistTree(popn[pop].BigHead[txPL[0]][loc-1]->Rtree, &k, pop);
	      k=0;
	      printf("par1(%ld) [%ld, %ld]; loc=%ld\n",txPL[1],
		     popn[pop].extremeMuts[txPL[1]][loc-1][0],
		     popn[pop].extremeMuts[txPL[1]][loc-1][1],
		     loc-1);
	      PrintHistTree(popn[pop].BigHead[txPL[1]][loc-1]->Rtree, &k,
			    pop);
	      k=0;
	      printf("pre-offsp (%ld) [%ld,%ld]; loc=%ld\n", 
		     tPL, popn[pop].extremeMuts[tPL][loc-1][0],
		     popn[pop].extremeMuts[tPL][loc-1][1], loc-1);
	      PrintHistTree(popn[pop].BigHead[tPL][loc-1]->Rtree, &k,pop);
	      fflush(stdout);
#endif
	      
	      createNewHistory(ind, loc-1, pop, 0, 0);
	      tPL = popn[pop].parentLoc[ind][loc-1]; /*AFTER CREATING NEW*/
	      if(bTR1 <= 0) /* completely covers previous locus! */
		bTR1 = 0; /* GC only extends to previous locus */
	      else if(bTR1>0){ /*copy back mutations carried by inhereted chrm*/
		copyPartialHistory(&popn[pop].BigHead[tPL][loc-1]->Rtree,
				   &popn[pop].BigHead[txPL[j]][loc-1]->Rtree,
				   0, bTR1-1, pop, tPL, loc-1);
		if(gpars.TRACKANC){
		  memcpy(popn[pop].ancestry[tPL][loc-1], 
			 popn[pop].ancestry[txPL[j]][loc-1], bTR1*sizeof(int));
		}
	      }
	      
#ifdef VERBOSE_DEBUG
	      k=0;
	      printf("med-offsp (%ld) [%ld,%ld]; loc=%ld\n", 
		     tPL, popn[pop].extremeMuts[tPL][loc-1][0],
		     popn[pop].extremeMuts[tPL][loc-1][1], loc-1);
	      PrintHistTree(popn[pop].BigHead[tPL][loc-1]->Rtree, &k,pop);
	      printf("convertTract(%ld: %ld-%ld) D=%d\n", loc-1, bTR1,
		     gpars.L[loc-1]-1, D);
	      k=0;
	      PrintHistTree(popn[pop].BigHead[txPL[D]][loc-1]->Rtree,&k,pop);
	      fflush(stdout);
#endif
	      if(F == 0){ /* strand synthesis from donor first */
		copyPartialHistory(&popn[pop].BigHead[tPL][loc-1]->Rtree,
				   &popn[pop].BigHead[txPL[D]][loc-1]->Rtree,
				   bTR1, bTR1+TR1-1, pop, tPL, loc-1);
		if(gpars.TRACKANC){
		  memcpy(&(popn[pop].ancestry[tPL][loc-1][bTR1]),
			 &(popn[pop].ancestry[txPL[D]][loc-1][bTR1]),
			 TR1*sizeof(int));
		}
		if(bTR1+TR1 <= gpars.L[loc-1]-1){
		  convertTract(&popn[pop].BigHead[tPL][loc-1]->Rtree, 
			       &popn[pop].BigHead[txPL[D]][loc-1]->Rtree, 
			       bTR1+TR1, gpars.L[loc-1]-1, pop, tPL, loc-1,
			       txPL[(D+1)%2]);
		  if(fabs(ppars[pop].BGC-0.5) > FLT_EPSILON)
		    addToTract(&popn[pop].BigHead[tPL][loc-1]->Rtree, 
			       &popn[pop].BigHead[txPL[(D+1)%2]][loc-1]->Rtree,
			       bTR1+TR1, gpars.L[loc-1]-1, pop, tPL, loc-1,
			       txPL[D]);
		  if(gpars.TRACKANC){
		    memcpy(&(popn[pop].ancestry[tPL][loc-1][bTR1+TR1]),
			   &(popn[pop].ancestry[txPL[D]][loc-1][bTR1+TR1]),
			   (gpars.L[loc-1]-(bTR1+TR1))*sizeof(int));
		  }
		}
	      }
	      else{ /* strand conversion from donor first */
		convertTract(&popn[pop].BigHead[tPL][loc-1]->Rtree, 
			     &popn[pop].BigHead[txPL[D]][loc-1]->Rtree, 
			     bTR1, bTR1+TR1-1, pop, tPL, loc-1, txPL[(D+1)%2]);
		if(fabs(ppars[pop].BGC-0.5) > FLT_EPSILON)
		  addToTract(&popn[pop].BigHead[tPL][loc-1]->Rtree, 
			     &popn[pop].BigHead[txPL[(D+1)%2]][loc-1]->Rtree,
			     bTR1, bTR1+TR1-1, pop, tPL, loc-1, txPL[D]);
		if(gpars.TRACKANC){
		  memcpy(&(popn[pop].ancestry[tPL][loc-1][bTR1]), 
			 &(popn[pop].ancestry[txPL[D]][loc-1][bTR1]),
			 TR1*sizeof(int));
		}
		if(bTR1+TR1 <= gpars.L[loc-1]-1){
		  copyPartialHistory(&popn[pop].BigHead[tPL][loc-1]->Rtree,
				     &popn[pop].BigHead[txPL[D]][loc-1]->Rtree,
				     bTR1+TR1, gpars.L[loc-1]-1,pop,tPL,loc-1);
		  if(gpars.TRACKANC){
		    memcpy(&(popn[pop].ancestry[tPL][loc-1][bTR1+TR1]),
			   &(popn[pop].ancestry[txPL[D]][loc-1][bTR1+TR1]),
			   (gpars.L[loc-1]-bTR1+TR1)*sizeof(int));
		  }
		}
	      }
	      if(gpars.substMod == 5){
		/* recalculate hit probabilities  */
		memcpy(popn[pop].hpSUM[tPL][loc-1],
		       popn[pop].hpSUM[txPL[j]][loc-1],
		       gpars.L[loc-1]*sizeof(double));
		getNewHitProbsReg(&popn[pop].hpSUM[tPL][loc-1],
				  &popn[pop].BigHead[tPL][loc-1],
				  loc-1, pop, bTR1, gpars.L[loc-1]-1);
		
#ifdef UNIT_TEST
		checkHitProbs(popn[pop].hpSUM[tPL][loc-1],
			      "error in getRecombinant", pop, tPL, loc-1);
#endif
	      }
	      if(popn[pop].BigHead[tPL][loc-1]->Rtree != NULL){
		if(ptr==NULL)
		  ptr = popTrash();
		ptr->Rtree = NULL;
		getSmallestHist(popn[pop].BigHead[tPL][loc-1]->Rtree, &ptr);
		if(ptr->Rtree != NULL){
		  popn[pop].extremeMuts[tPL][loc-1][0] =
		    ptr->Rtree->event->site;
		  getLargestHist(popn[pop].BigHead[tPL][loc-1]->Rtree, &ptr);
		  popn[pop].extremeMuts[tPL][loc-1][1] =
		    ptr->Rtree->event->site;
		  ptr->Rtree = NULL;
		}
		else{ /* carries no mutations */
		  popn[pop].extremeMuts[tPL][loc-1][0] = gpars.L[loc-1];
		  popn[pop].extremeMuts[tPL][loc-1][1] = 0;
		}
	      }
#ifdef VERBOSE_DEBUG
	      k=0;
	      printf("post-offsp (%ld) [%ld,%ld]; loc=%ld; bTR1=%ld\n", tPL,
		     popn[pop].extremeMuts[tPL][loc-1][0],
		     popn[pop].extremeMuts[tPL][loc-1][1],loc-1, bTR1);
	      PrintHistTree(popn[pop].BigHead[tPL][loc-1]->Rtree, &k,pop);
	      fflush(stdout);
#endif
	    }
	  }
	}
	
	max = within[iW].pos[i];
	
#ifdef VERBOSE_DEBUG
	printf("convertTract(%ld: %ld-%ld)\n", loc, min, max);
	fflush(stdout);
#endif
	if(F == 0){ /* strand synthesis from donor first */
	  if(bTR+TR1-1 >= 0){
	    min = (bTR>=0 ? bTR : 0);
	    copyPartialHistory(&popn[pop].BigHead[PL][loc]->Rtree,
			       &popn[pop].BigHead[xPL[D]][loc]->Rtree,
			       min, bTR+TR1-1, pop, PL, loc);
	    if(gpars.TRACKANC){
	      memcpy(&(popn[pop].ancestry[PL][loc][min]),
		     &(popn[pop].ancestry[xPL[D]][loc][min]),
		     (bTR+TR1-min)*sizeof(int));
	    }
	  }
	  min = (bTR+TR1>=0 ? bTR+TR1 : 0);
	  convertTract(&popn[pop].BigHead[PL][loc]->Rtree, 
		       &popn[pop].BigHead[xPL[D]][loc]->Rtree, min, max,
		       pop, PL, loc, xPL[(D+1)%2]);
	  if(fabs(ppars[pop].BGC-0.5) > FLT_EPSILON)
	    addToTract(&popn[pop].BigHead[PL][loc]->Rtree, 
		       &popn[pop].BigHead[xPL[(D+1)%2]][loc]->Rtree, min,
		       max, pop, PL, loc, xPL[D]);
	  if(gpars.TRACKANC){
	    memcpy(&(popn[pop].ancestry[PL][loc][min]), 
		   &(popn[pop].ancestry[xPL[D]][loc][min]),
		   (max-min+1)*sizeof(int));
	  }
	}
	else{ /* strand synthesis from donor first */
	  if(bTR+TR1-1 >= 0){
	    min = (bTR>=0 ? bTR : 0);
	    convertTract(&popn[pop].BigHead[PL][loc]->Rtree, 
			 &popn[pop].BigHead[xPL[D]][loc]->Rtree, min, bTR+TR1-1,
			 pop, PL, loc, xPL[(D+1)%2]);
	    if(fabs(ppars[pop].BGC-0.5) > FLT_EPSILON)
	      addToTract(&popn[pop].BigHead[PL][loc]->Rtree, 
			 &popn[pop].BigHead[xPL[(D+1)%2]][loc]->Rtree, min,
			 bTR+TR1-1, pop, PL, loc, xPL[D]);
	    if(gpars.TRACKANC){
	      memcpy(&(popn[pop].ancestry[PL][loc][min]), 
		     &(popn[pop].ancestry[xPL[(D+1)%2]][loc][min]),
		     (bTR+TR1-min)*sizeof(int));
	    }
	  }
	  min = (bTR+TR1>=0 ? bTR+TR1 : 0);
	  copyPartialHistory(&popn[pop].BigHead[PL][loc]->Rtree,
			     &popn[pop].BigHead[xPL[D]][loc]->Rtree,
			     min, max, pop, PL, loc);
	  if(gpars.TRACKANC){
	    memcpy(&(popn[pop].ancestry[PL][loc][min]),
		   &(popn[pop].ancestry[xPL[D]][loc][min]),
		   (max-min+1)*sizeof(int));
	  }
	}
      }
      else
	max = within[iW].pos[i];
    }
    else
      max = gpars.L[loc]-1;
    
    if(i == within[iW].numEvents ||
       (i < within[iW].numEvents && 
	(within[iW].loci[i] != loc || within[iW].type[i]==0))){
#ifdef VERBOSE_DEBUG
      printf("min=%ld; max=%ld; j=%d [%ld, %ld]\n",min,max,j,
	     popn[pop].extremeMuts[xPL[j]][loc][0], 
	     popn[pop].extremeMuts[xPL[j]][loc][1]);
      fflush(stdout);
#endif    
      
      if(popn[pop].extremeMuts[xPL[j]][loc][0] <= max &&
	 popn[pop].extremeMuts[xPL[j]][loc][1] >= min){ /* copy when necessary*/
	
#ifdef VERBOSE_DEBUG
	printf("copying %ld [%ld, %ld]\n",xPL[j],min,max);
	k=0;
	PrintHistTree(popn[pop].BigHead[xPL[j]][loc]->Rtree, &k,pop);
#endif
	
	copyPartialHistory(&popn[pop].BigHead[PL][loc]->Rtree,
			   &popn[pop].BigHead[xPL[j]][loc]->Rtree,
			   min, max, pop, PL, loc);
#ifdef VERBOSE_DEBUG
	printf("result %ld\n",PL);
	k=0;
	PrintHistTree(popn[pop].BigHead[PL][loc]->Rtree, &k,pop);
	fflush(stdout);
#endif
      }
      if(gpars.TRACKANC){
	memcpy(&(popn[pop].ancestry[PL][loc][min]),
	       &(popn[pop].ancestry[xPL[j]][loc][min]),
	       (max-min+1)*sizeof(int));
      }
    }
    if(i<within[iW].numEvents && within[iW].loci[i]==loc &&
       within[iW].type[i]%2==0)
      j = (j+1)%2;
    i++;
  }

  if(gpars.substMod == 5){
    /* recalculate hit probabilities between first and last event */
    getNewHitProbsReg(&popn[pop].hpSUM[PL][loc], &popn[pop].BigHead[PL][loc],
		      loc, pop, (first>=0?first:0),
		      (last<gpars.L[loc]-1?last:gpars.L[loc]-1));
    last += 2;
    
    if(last<gpars.L[loc]-1){ /* more sites to update */
      delta = -popn[pop].hpSUM[xPL[j]][loc][last];
      delta += popn[pop].hpSUM[PL][loc][last];
      if(j!=0){ /* cross-over, memcpy HP from other parent */
	memcpy(&popn[pop].hpSUM[PL][loc][last+1],
	       &popn[pop].hpSUM[xPL[j]][loc][last+1],
	       (gpars.L[loc]-last-1)*sizeof(double));
      }
      
      
      for(k=last+1; k<gpars.L[loc]; k++)
	popn[pop].hpSUM[PL][loc][k] += delta;
      popn[pop].hpSUM[PL][loc][gpars.L[loc]-1] =
	popn[pop].hpSUM[PL][loc][gpars.L[loc]-2];
    }
#ifdef UNIT_TEST
    checkHitProbs(popn[pop].hpSUM[PL][loc], "error in getRecombinant", pop, 
		  PL, loc);
    
#endif
  }
  
#ifdef UNIT_TEST
  if(gpars.R==1 && within[iW].numEvents==1 && ppars[pop].fGC<FLT_EPSILON){
    countMuts(popn[pop].BigHead[PL][0]->Rtree, &numMuts2, 0,
	      gpars.L[0]-1);
    if(numMuts1 != numMuts2){
      fflush(stdout);
      fflush(errfile);
      fprintf(errfile,"error in copying mutations from getRecombinant()\n");
      fprintf(errfile,"seed=%ld\n",INITIALSEED);
      fprintf(errfile,"parents: %ld(%ld) %ld(%ld);  rec site=%ld\n",x[0], 
	      xPL[0], x[1], xPL[1], within[iW].pos[0]);
      fprintf(errfile,"parent0:\n");
      fflush(errfile);
      long foo=0;
      PrintHistTree(popn[pop].BigHead[xPL[0]][0]->Rtree, &foo, pop);
      fprintf(errfile,"\nparent1:\n");
      fflush(errfile);
      foo=0;
      PrintHistTree(popn[pop].BigHead[xPL[1]][0]->Rtree, &foo, pop);
      fprintf(errfile,"\noffspring:\n");
      fflush(errfile);
      foo=0;
      PrintHistTree(popn[pop].BigHead[PL][0]->Rtree, &foo, pop);
      abort();
    }
  }
#endif
  
  numx--;/* off by 1 error: n intervals => n-1 events */
#ifdef VERBOSE_DEBUG
  printf("offspring (%ld.%ld)\n", loc, within[iW].pos[*iW2]);
  k=0;
  PrintHistTree(popn[pop].BigHead[PL][loc]->Rtree, &k, pop);
  printf("\n\n\n");
  fflush(stdout);
#endif
  (*newPar) = x[j];
  (*iW2) += numx;
  
  if(!ppars[pop].neutpop){
    popn[pop].indfit[PL][loc] = (gpars.ADDITIVE ? 0.0 : 1.0);
    getFitTree(popn[pop].BigHead[PL][loc]->Rtree, &popn[pop].indfit[PL][loc],
	       pop, loc);
  }
  
  if(popn[pop].BigHead[PL][loc]->Rtree != NULL){
    if(ptr==NULL)/* get smallest mutation */
      ptr = popTrash();
    getSmallestHist(popn[pop].BigHead[PL][loc]->Rtree, &ptr);
    if(ptr->Rtree != NULL){
      popn[pop].extremeMuts[PL][loc][0] = ptr->Rtree->event->site;
      if(popn[pop].extremeMuts[PL][loc][0] > popn[pop].extremeMuts[PL][loc][1])
	popn[pop].extremeMuts[PL][loc][1] = popn[pop].extremeMuts[PL][loc][0];
      /* get largest mutation site */
      ptr->Rtree = NULL;
      getLargestHist(popn[pop].BigHead[PL][loc]->Rtree, &ptr);
      popn[pop].extremeMuts[PL][loc][1] = ptr->Rtree->event->site;
     }
    else{ /* all mutations removed */
      popn[pop].extremeMuts[PL][loc][0] = gpars.L[loc];
      popn[pop].extremeMuts[PL][loc][1] = 0;
    }
  }
  else{
    popn[pop].extremeMuts[PL][loc][0] = gpars.L[loc];
    popn[pop].extremeMuts[PL][loc][1] = 0;
  }
  
  if(ptr!=NULL){
    ptr->Rtree = NULL;
    free(ptr);
    ptr = NULL;
  }
  return;
}

/* ------------------------------------------------------------------------- */

void mutate(double **mutClass_site, double *mutClass_locus, const int *selClass,
	    const double *fitQuant, const int pop)
{
  long nummutes;        /* Number of mutations each generation,
			  drawn from a Poisson w/ mean theta/2 */
  long mutreg;
  long mutind, i, r, mutcount = 0, mutindPL;
  unsigned long mutloc;
  /* mutind: chromosome in which a mutation is going to occur
     mutloc: site at which a mutation is to occur  */
  char  oc[3], nc[3], n[5];   /* old and new codons (before and after mutation 
				 (resp.)) */
  long m;
  char newnuc='0';     /* the nucleotide before and after it mutates (resp.) */
  double step;         /* equivalent to a counter */
  int CpG=0;           /* note whether site to mutate is CpG, default=0 */
  double fit=0.0;      /* record the value of the fitness for new mutation */
  struct event *event; /* hold data for linked list */
  int **tmpSITE;
  unsigned long *tmpL;
  double delta = 0;
  
  tmpL = lvector(0,gpars.R-1);
  assert(tmpSITE = malloc(gpars.R*sizeof(*tmpSITE)));
  for(r=0; r<gpars.R; r++)  tmpSITE[r] = NULL;
  
  assert(event = malloc(sizeof(struct event)));
  assert(event->genFix = malloc(gpars.NPOP*sizeof(long)));
  assert(event->genDead = malloc(gpars.NPOP*sizeof(long)));
  for(i=0; i<gpars.NPOP; i++){
    event->genFix[i] = 0;
    event->genDead[i] = -LONG_MAX;
  }
  
  /* DRAW NUMBER OF MUTATIONS */
  nummutes=(long)poidev(ppars[pop].THETA*gpars.totalNucs/gpars.P,&gpars.seed);
  
#ifdef UNIT_TEST
  fprintf(errfile,"nummutes (%ld) = %ld\n",popn[pop].gen, nummutes);
  fflush(errfile);
#endif
  
  /* DETERMINE MUTATION SITE AND TYPE */
  for(m=1; m<=nummutes; m++){
#ifdef VERBOSE_DEBUG
    printf("   mutation %ld\n",m);
    fflush(stdout);
#endif    
    if(++mutcount >= gpars.totalNucs){
      fprintf(errfile,"Too many mutation attempts in single generation (%ld)",
	      INITIALSEED);
      fprintf(errfile,"(mutcount = %lu; MAX = %.0f)\n", mutcount,
	      gpars.totalNucs);
      fprintf(errfile,"Try reducing the mutation rate.\n");
      abort();
    }
    
    /* RANDOMLY SELECT CHROMOSOME (1,...,PN)*/
    mutind=(long)(gpars.P*ppars[pop].N*ran1(&gpars.seed));
    
    if(gpars.R == 1)
      mutreg = 0;
    else{
      float rn = ran1(&gpars.seed);
      mutreg = invCDF(rn, mutClass_locus, 0, gpars.R-1);
    }
    
#ifdef VERBOSE_DEBUG
    printf("mutind=%ld(%ld) ",mutind, popn[pop].parentLoc[mutind][mutreg]);
    fflush(stdout);
    long foo=0;
#endif
    
    createNewHistory(mutind, mutreg, pop, 1, 1);
    mutindPL = popn[pop].parentLoc[mutind][mutreg];
    
#ifdef VERBOSE_DEBUG
    {
      printf("mutate ind=%ld(%ld) reg=%ld gen=%ld:\n",
	     mutind, mutindPL, mutreg, popn[pop].gen);
      printf("popn[%d].BigHead[%ld][%ld]:\n",pop,mutindPL,mutreg);
      fflush(stdout);
      foo = 0;
      PrintHistTree(popn[pop].BigHead[mutindPL][mutreg]->Rtree, &foo, pop);
      fflush(stdout);
    }
#endif
      
    /* SELECT MUTATION SITE */
    if(gpars.substMod == 5){
      do{
	getMutSite(&mutloc, popn[pop].hpSUM[mutindPL][mutreg], mutreg);
	if(popn[pop].conSeq[mutreg][mutloc] == 'N'){ /* don't mutate site... */
	  if(gpars.R == 1)
	    mutreg = 0;
	  else{
	    float rn = ran1(&gpars.seed);
	    mutreg = invCDF(rn, mutClass_locus, 0, gpars.R-1);
	  }
	  continue;
	}
	
#ifdef VERBOSE_DEBUG
	fprintf(errfile,"hpSUM[%ld][%ld][%ld] = %f;  ", mutind, mutreg, 
		gpars.L[mutreg]-1, 
		popn[pop].hpSUM[mutindPL][mutreg][gpars.L[mutreg]-1]);
	fflush(errfile);
#endif
	
	/* SAVE NUCLEOTIDE THAT IS TO MUTATE FOR COMPARISON */
	if(mutloc>1)
	  n[3] = retNucHistory(mutloc-2, &popn[pop].BigHead[mutindPL][mutreg],
			       mutreg, pop, 1);
	n[0] = retNucHistory(mutloc-1, &popn[pop].BigHead[mutindPL][mutreg],
			     mutreg, pop, 1);
	n[1] = retNucHistory(mutloc, &popn[pop].BigHead[mutindPL][mutreg],
			     mutreg, pop, 1);
	n[2] = retNucHistory(mutloc+1, &popn[pop].BigHead[mutindPL][mutreg],
			     mutreg, pop, 1);
	if(mutloc<gpars.L[mutreg]-2)
	  n[4] = retNucHistory(mutloc+2, &popn[pop].BigHead[mutindPL][mutreg],
			       mutreg, pop, 1);
	
	if((n[1]=='0' && n[2]=='1') || (n[0]=='0' && n[1]=='1'))
	  CpG=1;
	else
	  CpG=0;
	
	/* CHOOSE MUTATION TYPE */
	newnuc=MutantNuc(n, CpG, pop);
	/* SAVE NEW CODON FOR COMPARISON */
	if((mutloc%3)==0 && mutloc<gpars.L[mutreg]-2){
	  nc[0] = newnuc;
	  nc[1] = n[2];
	  nc[2] = n[4];
	  oc[0] = n[1];
	  oc[1] = nc[1];
	  oc[2] = nc[2];
	}
	else if((mutloc%3)==1){
	  nc[0] = n[0];
	  nc[1] = newnuc;
	  nc[2] = n[2];
	  oc[0] = n[0];
	  oc[1] = n[1];
	  oc[2] = n[2];
	}
	else if(mutloc>1){
	  nc[0] = n[3];
	  nc[1] = n[0];
	  nc[2] = newnuc;
	  oc[0] = nc[0];
	  oc[1] = nc[1];
	  oc[2] = n[1];
	}
	else{
	  fprintf(errfile,"error (%ld)\n",INITIALSEED);
	  exitNOW("error in mutate! (site out of bounds)\n");
	}
      }while((gpars.ANNOTATE[mutreg]=='C' && AA(nc) == 0) ||
	     (gpars.SKIPSITES != NULL && gpars.SKIPSITES[mutreg][mutloc]>0)); 
      /* don't allow stop codons */
      if((newnuc=='0' && n[2]=='1') || (n[0]=='0' && newnuc=='1'))
	CpG=1;
    }
    else{
      do{
	do{
	  if(ppars[pop].RateClassSites == 1){ /* site in [0,L-1] */
	    do{
	      mutloc = (long)(ran1(&gpars.seed)*gpars.L[mutreg]);
	    }while(mutloc >= gpars.L[mutreg]);
	  }
	  else{
	    step=ran1(&gpars.seed);  /* select randomly according to site 
					specific mutation rate */
	    for(i=0; i<gpars.L[mutreg]; i++){
	      if(step <= mutClass_site[mutreg][i])
		break;
	    }
	    mutloc = i;
	  }
	  if(mutloc > 1){
	    n[3] = retNucHistory(mutloc-2, &popn[pop].BigHead[mutindPL][mutreg],
				 mutreg, pop, 1);
	  }
	  else
	    n[3] = 'N';
	  if(mutloc > 0)
	    n[0] = retNucHistory(mutloc-1, &popn[pop].BigHead[mutindPL][mutreg],
				 mutreg, pop, 1);
	  else
	    n[0] = 'N';
	  n[1] = retNucHistory(mutloc, &popn[pop].BigHead[mutindPL][mutreg],
			       mutreg, pop, 1);
	  if(mutloc < gpars.L[mutreg]-1)
	    n[2] = retNucHistory(mutloc+1, &popn[pop].BigHead[mutindPL][mutreg],
				 mutreg, pop, 1);
	  else
	    n[2] = 'N';
	  if(mutloc<gpars.L[mutreg]-2)
	    n[4] = retNucHistory(mutloc+2, &popn[pop].BigHead[mutindPL][mutreg],
				 mutreg, pop, 1);
	  else
	    n[4] = 'N';
	  if((n[0]=='0' && n[1]=='1') || (n[1]=='0' && n[2]=='1'))
	    CpG=1;
	  else
	    CpG=0;
	}while((gpars.SKIPSITES != NULL && gpars.SKIPSITES[mutreg][mutloc]>0) ||
	       ((gpars.substMod == 1 || gpars.substMod == 3) && 
		CpG == 0 && ran1(&gpars.seed) < ppars[pop].PSI));
	
	/* CHOOSE MUTATION TYPE */
	newnuc=MutantNuc(n, CpG, pop);
	/* SAVE NEW CODON FOR COMPARISON */
	if((mutloc%3)==0){
	  nc[0] = newnuc;
	  nc[1] = n[2];
	  nc[2] = n[4];
	  oc[0] = n[1];
	  oc[1] = nc[1];
	  oc[2] = nc[2];
	}
	else if((mutloc%3)==1){
	  nc[0] = n[0];
	  nc[1] = newnuc;
	  nc[2] = n[2];
	  oc[0] = nc[0];
	  oc[1] = n[1];
	  oc[2] = nc[2];
	}
	else{
	  nc[0] = n[3];
	  nc[1] = n[0];
	  nc[2] = newnuc;
	  oc[0] = nc[0];
	  oc[1] = nc[1];
	  oc[2] = n[1];
	}
	if((newnuc=='0' && n[2]=='1') || (n[0]=='0' && newnuc=='1'))
	  CpG=1;
      }while(gpars.ANNOTATE[mutreg]=='C' && AA(nc) == 0); 
      /* don't allow stop codons */
    }
    if(gpars.INFSITES == 1){
      int cnt = 0;
      for(i=0; i<gpars.NPOP; i++){
	if(ppars[i].ALIVE && popn[i].polySites[mutreg][mutloc] > 0){
	  cnt++;
	  m--;
	  break;
	}
      }
      if(cnt>0){
#ifdef VERBOSE_DEBUG
	fprintf(errfile,"mutation rejected\n");
	fflush(errfile);
#endif
	continue;
      }
    }
    /* don't allow more than 1 mutation at site in a given generation */
    if(tmpL[mutreg] <= mutloc){
      if(tmpL[mutreg] == 0)
	assert(tmpSITE[mutreg] = malloc((mutloc+1)*sizeof(**tmpSITE)));
      else
	assert(tmpSITE[mutreg] = 
	       realloc(tmpSITE[mutreg],(mutloc+1)*sizeof(**tmpSITE)));
      for(i=tmpL[mutreg]; i<=mutloc; i++)
	tmpSITE[mutreg][i] = 0;
      tmpSITE[mutreg][mutloc] = 1;
      tmpL[mutreg] = mutloc+1;
    }
    else{
      if(tmpSITE[mutreg][mutloc] == 1){
	m--;
#ifdef VERBOSE_DEBUG
	fprintf(errfile,"mutation rejected\n");
	fflush(errfile);
#endif
	continue;
      }
      else
	tmpSITE[mutreg][mutloc] = 1;
    }
    if(ppars[pop].f0[mutreg] < 1.0-DBL_EPSILON && 
       ((gpars.ANNOTATE[mutreg] == 'C' && AA(oc)!=AA(nc)) || 
	gpars.ANNOTATE[mutreg] == 'N') && 
       (ran1(&gpars.seed) > ppars[pop].f0[mutreg]))
      fit = -1.0;
    else if(ppars[pop].neutpop ||
	    (gpars.ANNOTATE[mutreg]=='C' && (AA(oc)==AA(nc))) || 
	    (ppars[pop].propSelLoci>FLT_EPSILON && selClass[mutreg]==0))
      fit=0;
    else{
      fit=NewFits(fitQuant, pop, mutreg); 
    }
    event->site = mutloc;
    event->gen = popn[pop].gen;
    event->genFix[pop] = 0;
    event->genDead[pop] = -LONG_MAX;
    event->ancNuc = n[1];
    event->derNuc = newnuc;
    event->nonsyn =(char)(gpars.ANNOTATE[mutreg]=='C'?(AA(oc)!=AA(nc))+'0':'0');
    event->ancAA = (gpars.ANNOTATE[mutreg]=='C' ? AA(oc) : 21);
    event->derAA = (gpars.ANNOTATE[mutreg]=='C' ? AA(nc) : 21);
    event->fit = fit;
    event->CpG = (char)(CpG+'0');
    event->fiveP = n[0];
    event->threeP= n[2];
    if(gpars.SEX[mutreg] == '0')
      event->axy = 'A'; /* autosome */
    else{ /* sex chromosome */
      if(mutind < gpars.P*ppars[pop].MALES || mutind%gpars.P < gpars.P/2)
	event->axy = 'X';  /* X-linked */
      else
	event->axy = 'Y';  /* Y-linked */
    }

#ifdef UNIT_TEST
    checkHistTreeOrder(popn[pop].BigHead[mutindPL][mutreg]->Rtree, 0,
		       gpars.L[mutreg], pop);
    
    if(popn[pop].BigHead[mutindPL][mutreg] == NULL){
      fprintf(errfile,"(*BigHead)[%lu][%ld] == NULL!!! (%ld)\n",mutindPL,mutreg,
	      INITIALSEED);
      abort();
    }
#endif

    addHistoryNode(event, &popn[pop].BigHead[mutindPL][mutreg]->Rtree, NULL,
		   pop, mutindPL, mutreg, 0);
    popn[pop].polySites[mutreg][mutloc]++;
#ifdef UNIT_TEST
    if(popn[pop].polySites[mutreg][mutloc]>1){
      fprintf(errfile,"ps[%d][%ld][%ld] -> %d!!!\n",pop, mutreg,
	      mutloc,
	      popn[pop].polySites[mutreg][mutloc]);
      fflush(errfile);
#ifdef VERBOSE_DEBUG
      //getchar();
#endif
    }
#endif

#ifdef VERBOSE_DEBUG
    {
      printf(">>r=%ld,site=%ld(%ld),gen=%ld(%ld) = ",mutreg,mutloc,
	     event->site, popn[pop].gen, event->gen);
      if(event->gen != popn[pop].gen){
	printf("gen not right!!!\n");
	abort();
      }
      i = mutAr[mutreg]->mutIndex;
      printf("%3ld: %4lu %5ld %5ld %c %c %c %2i %2i %1.2f %c %d %c %c[%c%c]\n",
	     i, mutAr[mutreg]->muts[i]->event->site,
	     mutAr[mutreg]->muts[i]->event->gen,
	     mutAr[mutreg]->muts[i]->event->genFix[pop],
	     mutAr[mutreg]->muts[i]->event->ancNuc,
	     mutAr[mutreg]->muts[i]->event->derNuc,
	     mutAr[mutreg]->muts[i]->event->nonsyn,
	     mutAr[mutreg]->muts[i]->event->ancAA, 
	     mutAr[mutreg]->muts[i]->event->derAA,
	     mutAr[mutreg]->muts[i]->event->fit,
	     mutAr[mutreg]->muts[i]->event->CpG,
	     popn[pop].polySites[mutreg][mutAr[mutreg]->muts[i]->event->site],
	     mutAr[mutreg]->muts[i]->event->fiveP,
	     mutAr[mutreg]->muts[i]->event->threeP,
	     mutAr[mutreg]->muts[i]->event->free,
	     mutAr[mutreg]->muts[i]->event->fixed[pop]);
      /* 	abort(); */
      fflush(stdout);
    }
#endif

    tmpSITE[mutreg][mutloc] = 1;
    if(gpars.ADDITIVE)
      popn[pop].indfit[mutindPL][mutreg] += event->fit;
    else
      popn[pop].indfit[mutindPL][mutreg] *= 1.0+event->fit;
    /* check for lethal mutation! */
    if(gpars.ADDITIVE && popn[pop].indfit[mutindPL][mutreg] < -1.0)
      popn[pop].indfit[mutindPL][mutreg] = -1.0;
    else if(!gpars.ADDITIVE && popn[pop].indfit[mutindPL][mutreg] < 0) 
      popn[pop].indfit[mutindPL][mutreg] = 0.0;
    if(popn[pop].extremeMuts[mutindPL][mutreg][0] > mutloc)
      popn[pop].extremeMuts[mutindPL][mutreg][0] = mutloc;
    if(popn[pop].extremeMuts[mutindPL][mutreg][1] < mutloc)
      popn[pop].extremeMuts[mutindPL][mutreg][1] = mutloc;

#ifdef UNIT_TEST
#ifdef VERBOSE_DEBUG
    fprintf(errfile,"mut at mutindPL %ld, site %ld, gen %ld\n",mutindPL, mutloc,
	    popn[pop].gen);
    fflush(errfile);
    foo=0;
    PrintHistTree(popn[pop].BigHead[mutindPL][mutreg]->Rtree, &foo, pop);
    fflush(stdout);
#endif /* VERBOSE_DEBUG */

    checkHistTreeOrder(popn[pop].BigHead[mutindPL][mutreg]->Rtree, 0,
		       gpars.L[mutreg], pop);
#endif /* UNIT_TEST */
    
#ifdef VERBOSE_DEBUG
    foo=0;
    printf("end of mutate:\n");
    fflush(stdout);
    PrintHistTree(popn[pop].BigHead[mutindPL][mutreg]->Rtree, &foo, pop);
    fflush(stdout);
    if(gpars.substMod == 5)
      fprintf(errfile,"mutate:  %ld %ld %c%c%c%c%c  %f %f %f\n", 
	      mutind, mutloc, n[3], n[0], n[1], n[2], n[4],
	      popn[pop].hpSUM[mutindPL][mutreg][mutloc-1], 
	      popn[pop].hpSUM[mutindPL][mutreg][mutloc], 
	      popn[pop].hpSUM[mutindPL][mutreg][mutloc+1]);
#endif /* VERBOSE_DEBUG */
    if(gpars.substMod == 5){
      delta = 0;
      if(mutloc>1){
	for(i=0;i<4;i++)
	  if(i != n[0]-'0')
	    delta += Q[16*(n[3]-'0')+4*(n[0]-'0')+(newnuc-'0')][i] - 
	      Q[16*(n[3]-'0')+4*(n[0]-'0')+(n[1]-'0')][i];
	popn[pop].hpSUM[mutindPL][mutreg][mutloc-1] += delta;
      }
      
      for(i=0;i<4;i++){
	if(i != n[1]-'0')
	  delta -= Q[16*(n[0]-'0')+4*(n[1]-'0')+(n[2]-'0')][i];
	if(i != newnuc-'0')
	  delta += Q[16*(n[0]-'0')+4*(newnuc-'0')+(n[2]-'0')][i];
      }
      popn[pop].hpSUM[mutindPL][mutreg][mutloc] += delta;
      
      if(mutloc<gpars.L[mutreg]-2){
	for(i=0;i<4;i++)
	  if(i != n[2]-'0')
	    delta += Q[16*(newnuc-'0')+4*(n[2]-'0')+(n[4]-'0')][i] - 
	      Q[16*(n[1]-'0')+4*(n[2]-'0')+(n[4]-'0')][i];
	popn[pop].hpSUM[mutindPL][mutreg][mutloc+1] += delta;
      }
      else
	popn[pop].hpSUM[mutindPL][mutreg][gpars.L[mutreg]-1] = 
	  popn[pop].hpSUM[mutindPL][mutreg][gpars.L[mutreg]-2];
      
      for(i=mutloc+2; i<gpars.L[mutreg]; i++)
	popn[pop].hpSUM[mutindPL][mutreg][i] += delta;
    }   
#ifdef VERBOSE_DEBUG
    if(mutloc > 0 && mutloc < gpars.L[mutreg]-1){
      fprintf(errfile,"getting nucs:  %c%c%c\n\n", 
	      retNucHistory(mutloc-1, &popn[pop].BigHead[mutindPL][mutreg],
			    mutreg, pop, 0), 
	      retNucHistory(mutloc, &popn[pop].BigHead[mutindPL][mutreg],
			    mutreg, pop, 0), 
	      retNucHistory(mutloc+1, &popn[pop].BigHead[mutindPL][mutreg],
			    mutreg, pop, 0));
    }
#endif
  }
  free(event->genFix);
  free(event->genDead);
  free(event);
  for(i=0; i<gpars.R; i++){
    if(tmpSITE[i] != NULL)
      free(tmpSITE[i]);
  }
  free(tmpSITE);
  free_lvector(tmpL, 0, gpars.R);
}

/* ------------------------------------------------------------------------- */

void indels(const int pop, const double *fitQuant)
{
  long i, j, ii, ix, ir, ipl, step;
  long numID, preLen=0;
  float rn;
  double sumRates;
  int IDlen, maxlen=100, type;
  struct event *indel;
  
  sumRates = (ppars[pop].INSRATE + ppars[pop].DELRATE +
	      ppars[pop].longINSRATE + ppars[pop].longDELRATE);
  if(sumRates < FLT_EPSILON)
    return;
  /* determine number of events as sum of poissons */
  numID = poidev(sumRates*gpars.totalNucs/gpars.P, &gpars.seed);
  
  if(numID == 0)
    return;

  assert(indel = malloc(sizeof(struct event)));
  assert(indel->genFix = malloc(gpars.NPOP*sizeof(long)));
  assert(indel->genDead = malloc(gpars.NPOP*sizeof(long)));
  assert(indel->nucs = malloc((maxlen+1)*sizeof(char)));;
  indel->genFix[pop] = 0;
  indel->genDead[pop] = -LONG_MAX;
  
  for(i=0; i<numID; i++){
    /* determine type by competing rates */
    rn = ran1(&gpars.seed);
    if(rn < ppars[pop].INSRATE/sumRates)
      type = 0;
    else if(rn < (ppars[pop].INSRATE+ppars[pop].DELRATE)/sumRates)
      type = 1;
    else if(rn < (ppars[pop].INSRATE+ppars[pop].DELRATE+
		  ppars[pop].longINSRATE)/sumRates)
      type = 2;
    else
      type = 3;
    
    /* randomly choose chromosome */
    ii = (long)(ran1(&gpars.seed)*gpars.P*ppars[pop].N);

    do{
      /* length of insertion ~geom(1/mean_length) */
      if(type < 2)
	IDlen = geomdev(1/ppars[pop].INDELlength, &gpars.seed);
      else
	IDlen = poidev(ppars[pop].longINDELlength, &gpars.seed);
      
      if(IDlen > maxlen && type%2 == 0){
	assert(indel->nucs = realloc(indel->nucs, (IDlen+1)*sizeof(char)));
	maxlen = IDlen;
      }
      
      /* randomly choose locus wp given by relative lengths */
      if(gpars.R == 1){
	ir = 0;
	if(gpars.ANNOTATE[0] == 'C'){
	  printf("insertion/deletions only available for non-coding regions\n");
	  printf("use option --annotate (-a) to identify non-coding regions\n");
	  abort();
	}
      }
      else{
	step = 0;
	do{
	  rn = ran1(&gpars.seed);
	  ir = invCDF(rn*gpars.sumL[gpars.R-1], gpars.sumL, 0, gpars.R-1);
	  step++;
	}while(gpars.ANNOTATE[ir] == 'C' && step<14*gpars.R);
	if(step==14*gpars.R){
	  printf("no non-coding regions found for indels after 14R steps\n");
	  printf("use option --annotate (-a) to identify non-coding regions\n");
	  abort();
	}
      }
      /* randomly choose site in locus */
      ix = (long)(ran1(&gpars.seed)*(gpars.L[ir]-2))+1;
      if(popn[pop].conSeq[ir][ix] == 'N'){
	i--;
	continue;
      }
      if(ir == 0)
	preLen = ix;
      else
	preLen = gpars.sumL[ir-1]+ix;
    }while(type%2==1 && gpars.sumL[ir]-preLen < IDlen);
    if(type == 0 || type == 2){
      indel->ancNuc = '-'; /* ancestor had nothing */
      indel->derNuc = '+'; /* bases added */
      indel->nonsyn = 'i'; /* 'i'nsertion */
    }
    else{
      indel->ancNuc = '+'; /* ancestor had something */
      indel->derNuc = '-'; /* bases deleted */
      indel->nonsyn = 'd'; /* 'd'eletion */
    }
    
    if(gpars.SEX[ir] == '0')
      indel->axy = 'A'; /* autosome */
    else{ /* sex chromosome */
      if(ii < gpars.P*ppars[pop].MALES || ii%gpars.P < gpars.P/2)
	indel->axy = 'X';  /* X-linked */
      else
	indel->axy = 'Y';  /* Y-linked */
    }

    indel->CpG = '0'; /* doesn't matter */
    if(ix > 0)
      indel->fiveP = popn[pop].conSeq[ir][ix-1];
    else if(ir > 0)
      indel->fiveP = popn[pop].conSeq[ir-1][gpars.L[ir-1]-1];
    else
      indel->fiveP = '-';
    if(ix < gpars.L[ir]-1)
      indel->threeP = popn[pop].conSeq[ir][ix+1];
    else if(ir < gpars.R-1)
      indel->threeP = popn[pop].conSeq[ir+1][0];
    else
      indel->threeP = '-';

    indel->ancAA = 21; /* default for nothing */
    indel->derAA = 21; /* default for nothing */
    indel->gen = popn[pop].gen;
    indel->site = ix;
    indel->fit = NewFits(fitQuant, pop, ir);
    indel->nSites = IDlen;
    if(type%2 == 0){
      for(j=0; j<IDlen; j++){
	rn = ran1(&gpars.seed);
	if(rn < 0.25)
	  indel->nucs[j] = '0';
	else if(rn < 0.5)
	  indel->nucs[j] = '1';
	else if(rn < 0.75)
	  indel->nucs[j] = '2';
	else
	  indel->nucs[j] = '3';
      }
      indel->nucs[j] = '\0';
    }
    else{
      indel->nucs[0] = '-';
      indel->nucs[1] = '\0';
    }

    createNewHistory(ii, ir, pop, 1, 1);
    ipl = popn[pop].parentLoc[ii][ir];

    addHistoryNode(indel,&popn[pop].BigHead[ipl][ir]->Rtree,NULL,pop,ipl,ir,0);
    
    if(gpars.ADDITIVE)
      popn[pop].indfit[ipl][ir] += indel->fit;
    else
      popn[pop].indfit[ipl][ir] *= 1.0+indel->fit;
    if(popn[pop].indfit[ipl][ir] < 0) /* lethal mutation! */
      popn[pop].indfit[ipl][ir] = 0.0;
    if(popn[pop].extremeMuts[ipl][ir][0] > ix)
      popn[pop].extremeMuts[ipl][ir][0] = ix;
    if(popn[pop].extremeMuts[ipl][ir][1] < ix)
      popn[pop].extremeMuts[ipl][ir][1] = ix;
  }

  free(indel->nucs);
  free(indel->genFix);
  free(indel->genDead);
  free(indel);
}

/* ------------------------------------------------------------------------- */

void inversions(const int pop, const double *fitQuant)
{
  long i, ii, ix, ir, ipl, step;
  long numID, preLen;
  float rn;
  int IDlen, maxlen=1;
  struct event *inv;
  
  /* determine number of events as sum of poissons */
  numID = poidev(ppars[pop].INVRATE*gpars.totalNucs/gpars.P, &gpars.seed);
  if(numID == 0)
    return;

  assert(inv = malloc(sizeof(struct event)));
  assert(inv->genFix = malloc(gpars.NPOP*sizeof(long)));
  assert(inv->genDead = malloc(gpars.NPOP*sizeof(long)));
  assert(inv->nucs = malloc((maxlen+1)*sizeof(char)));;
  inv->genFix[pop] = 0;
  inv->genDead[pop] = -LONG_MAX;

  for(i=0; i<numID; i++){
    /* randomly choose chromosome */
    ii = (long)(ran1(&gpars.seed)*gpars.P*ppars[pop].N);

    do{
      /* length of inversion ~pois(mean_length) */
      IDlen = poidev(ppars[pop].INVlength, &gpars.seed);
      
      /* randomly choose locus wp given by relative lengths */
      if(gpars.R == 1){
	ir = 0;
	if(gpars.ANNOTATE[0] == 'C'){
	  printf("insertion/deletions only available for non-coding regions\n");
	  printf("use option --annotate (-a) to identify non-coding regions\n");
	  abort();
	}
      }
      else{
	step = 0;
	do{
	  rn = ran1(&gpars.seed);
	  ir = invCDF(rn*gpars.sumL[gpars.R-1], gpars.sumL, 0, gpars.R-1);
	  step++;
	}while(gpars.ANNOTATE[ir] == 'C' && step<14*gpars.R);
	if(step==14*gpars.R){
	  printf("no non-coding regions found for inversions\n");
	  printf("use option --annotate (-a) to identify non-coding regions\n");
	  abort();
	}
      }
      /* randomly choose site in locus */
      ix = (long)(ran1(&gpars.seed)*(gpars.L[ir]-2))+1;
      if(ir == 0)
	preLen = ix;
      else
	preLen = gpars.sumL[ir-1]+ix;
    }while(gpars.sumL[ir]-preLen < IDlen);
    
    if(popn[pop].conSeq[ir][ix] == 'N'){
      i--;
      continue;
    }
    
    inv->ancNuc = 'f'; /* ancestor was forward */
    inv->derNuc = 'r'; /* bases swapped (reversed)*/
    inv->nonsyn = 'v'; /* in'v'ersion */
    
    if(gpars.SEX[ir] == '0')
      inv->axy = 'A'; /* autosome */
    else{ /* sex chromosome */
      if(ii < gpars.P*ppars[pop].MALES || ii%gpars.P < gpars.P/2)
	inv->axy = 'X';  /* X-linked */
      else
	inv->axy = 'Y';  /* Y-linked */
    }
    
    inv->CpG = '0'; /* doesn't matter */
    if(ix > 0)
      inv->fiveP = popn[pop].conSeq[ir][ix-1];
    else if(ir > 0)
      inv->fiveP = popn[pop].conSeq[ir-1][gpars.L[ir-1]-1];
    else
      inv->fiveP = '-';
    if(ix < gpars.L[ir]-1)
      inv->threeP = popn[pop].conSeq[ir][ix+1];
    else if(ir < gpars.R-1)
      inv->threeP = popn[pop].conSeq[ir+1][0];
    else
      inv->threeP = '-';
    
    inv->ancAA = 21; /* default for nothing */
    inv->derAA = 21; /* default for nothing */
    inv->gen = popn[pop].gen;
    inv->site = ix;
    inv->fit = NewFits(fitQuant, pop, ir);
    inv->nSites = IDlen;
    inv->nucs[0] = '\0';

    createNewHistory(ii, ir, pop, 1, 1);
    ipl = popn[pop].parentLoc[ii][ir];

    addHistoryNode(inv,&popn[pop].BigHead[ipl][ir]->Rtree,NULL,pop,ipl,ir,0);
    
    if(gpars.ADDITIVE)
      popn[pop].indfit[ipl][ir] += inv->fit;
    else
      popn[pop].indfit[ipl][ir] *= 1.0+inv->fit;
    if(popn[pop].indfit[ipl][ir] < 0) /* lethal mutation! */
      popn[pop].indfit[ipl][ir] = 0.0;
    if(popn[pop].extremeMuts[ipl][ir][0] > ix)
      popn[pop].extremeMuts[ipl][ir][0] = ix;
    if(popn[pop].extremeMuts[ipl][ir][1] < ix)
      popn[pop].extremeMuts[ipl][ir][1] = ix;
  }

  free(inv->nucs);
  free(inv->genFix);
  free(inv->genDead);
  free(inv);
}

/* ------------------------------------------------------------------------- */

char MutantNuc(const char *nuc, const int CpG, const int pop)
{
  int i,j,n;
  float table[6][4];
  float rn, sum=0, sum1=0;
  char c0,c1,c2, ret;
  /*********************************************************
   *  Enter the transition table such that the (i,j) entry *
   *  is the probability of mutation from i to j given     *
   *  a mutation at i (i=4,5 indicates CpG)                *
   *********************************************************/
  
  rn=ran1(&gpars.seed);
  c0=nuc[0];
  c1=nuc[1];
  c2=nuc[2];
  if(gpars.substMod < 3 || (gpars.substMod==3 && CpG==0) ||
     gpars.substMod == 6){
    /* nucleotide model, draw from transProb */
    return((rn<ppars[pop].transProb[c1-'0'][0] ? '0' :
	    (rn<ppars[pop].transProb[c1-'0'][1] ? '1' :
	     (rn<ppars[pop].transProb[c1-'0'][2] ? '2' : '3'))));
  }
  else if(gpars.substMod==3 && CpG==0){
    /* let a=Pr(transversion), b=Pr(transition), solve:
     * 2*a+b = 1; b = KAPPA*a */
    if(rn<=ppars[pop].KAPPACpG/(ppars[pop].KAPPACpG+2))
      return(TsNuc(c1));  /* this is a transition */
    else
      return(TvNuc(c1));   /* this is one of the transversions */
  }
  else if(gpars.substMod == 4){
    table[0][1]=.134;  /* C->G   rates from Zhang & Gerstein, 2003, 
			  directly from probability, assuming K2P */
    table[0][2]=.732;  /* C->T */
    table[0][3]=.134;  /* C->A */
    table[1][0]=.120;  /* G->C */
    table[1][2]=.120;  /* G->T */
    table[1][3]=.760;  /* G->A */
    table[2][0]=.646;  /* T->C */
    table[2][1]=.177;  /* T->G */
    table[2][3]=.177;  /* T->A */
    table[3][0]=.1965; /* A->C */
    table[3][1]=.607;  /* A->G */
    table[3][2]=.1965; /* A->T */
    table[4][1]=.091;  /* C->G and C is CpG */
    table[4][2]=.828;  /* C->T and C is CpG */
    table[4][3]=.091;  /* C->A and C is CpG */
    table[5][0]=.0775; /* G->C and G is CpG */
    table[5][2]=.0775; /* G->T and G is CpG */
    table[5][3]=.845;  /* G->A and G is CpG */
    
    n=c1-'0';
    if(rn<=table[n+4*CpG][(n+1)%4])
      return (NumToNuc((n+1)%4));
    else if(rn<=table[n+4*CpG][(n+1)%4]+table[n+4*CpG][(n+2)%4])
      return (NumToNuc((n+2)%4));
    else
      return (NumToNuc((n+3)%4));
  }
  else if(gpars.substMod == 5){
#ifdef UNIT_TEST
    if(c0-'0' < 0 || c0-'0' > 3 ||
       c1-'0' < 0 || c1-'0' > 3 ||
       c2-'0' < 0 || c2-'0' > 3){
      fprintf(errfile,"bad nucleotide (%ld): %c%c%c\n",INITIALSEED,c0,c1,c2);
      exitNOW("error in MutantNuc\n");
    }
#endif /* UNIT_TEST */
    sum=0.0;
    for(i=1;i<=3;i++)
      sum+=Q[16*(c0-'0')+4*(c1-'0')+(c2-'0')][((c1-'0')+i)%4];
    j=0;
    sum1=0.0;
    while(rn>sum1 && ++j<3)
      sum1+=Q[16*(c0-'0')+4*(c1-'0')+(c2-'0')][((c1-'0')+j)%4]/sum;
    ret=NumToNuc(((c1-'0')+j)%4);
    if(ret-'0'>4){
      fprintf(errfile,"%c (%ld)",ret,INITIALSEED);
      exitNOW("Crap!  got bad nuc...");
      return (-1);
    }
    return(NumToNuc(((c1-'0')+j)%4));
  }
  else{
    fprintf(errfile,"error, did not find substMod = %d\n",gpars.substMod);
    abort();
  }
}

/* ------------------------------------------------------------------------- */

double NewFits(const double *fitQuant, const int pop, const long locus)
{
  float rn;
  double fit=0.0;

  if(ppars[pop].selDistType[locus] == 0 || 
     (ppars[pop].selDistType[locus] == 1 && 
      ppars[pop].ProbNeut[locus] >= 1.0-DBL_EPSILON))
    fit = 0;
  else if(ppars[pop].selDistType[locus] == 1){
    /* DISCRETE DISTRIBUTION */
    rn=ran1(&gpars.seed);
    if(rn<=ppars[pop].ProbDel[locus])
      fit = -ppars[pop].GAMMA[locus]/gpars.P/ppars[pop].N;
    else if(rn<=(ppars[pop].ProbDel[locus]+ppars[pop].ProbNeut[locus]))
      fit = 0.0;
    else
      fit = ppars[pop].GAMMA[locus]/gpars.P/ppars[pop].N;
  }
  else if(ppars[pop].selDistType[locus] == 2){  /* MIXED-GAMMA DISTRIBUTION */
    rn=ran1(&gpars.seed);
    if(rn<=ppars[pop].ProbNegGamma[locus])
      fit = -rgamma(ppars[pop].alphaN[locus],ppars[pop].lambdaN[locus],
		    &gpars.seed);
    else
      fit = rgamma(ppars[pop].alphaP[locus],ppars[pop].lambdaP[locus],
		   &gpars.seed);
    fit /= (gpars.P*ppars[pop].N);
  }
  else if(ppars[pop].selDistType[locus] == 3)     /* NORMAL */
    fit=(rnormal(ppars[pop].normMean[locus],ppars[pop].normVar[locus],
		 &gpars.seed)/gpars.P/ppars[pop].N);
  else if(ppars[pop].selDistType[locus]==4){      /* some other creation... */
    rn = 100*ran1(&gpars.seed);
    fit= -fitQuant[(int)rn]/gpars.P/ppars[pop].N;
  }
  return(fit);
}

/* ------------------------------------------------------------------------- */

void Split(const int fromPop, const int toPop)
{
  long i, k, r;
  
  /* INITIALIZE NEW POPULATION */
  if(toPop >= gpars.NPOPDEF){
    assert(popn = realloc(popn, (toPop+1)*sizeof(struct population)));
    gpars.NPOPDEF = toPop+1;
  }

  /* copy all ppars except certain ones */
  ppars[toPop].ALIVE = 1;
  ppars[toPop].N = ppars[fromPop].N;
  ppars[toPop].Nt = ppars[fromPop].Nt;
  /* don't copy pFEMALES */
  ppars[toPop].MALES = ppars[fromPop].MALES;
  /* don't copy pMaleMig or popAlpha, tauAlpha, K, P0, GenEffect, SS */
  ppars[toPop].THETA = ppars[fromPop].THETA;
  ppars[toPop].INSRATE = ppars[fromPop].INSRATE;
  ppars[toPop].DELRATE = ppars[fromPop].DELRATE;
  ppars[toPop].RHO = ppars[fromPop].RHO;
  /* don't copy SELF */
  if(ppars[toPop].neutpop == 1){
    for(i=0; i<gpars.R; i++)
      ppars[toPop].selDistType[i] = 0;
  }
  else{
    for(r=0; r<gpars.R; r++){
      ppars[toPop].selDistType[r] = ppars[fromPop].selDistType[r];
      /* don't copy selDistType propSelLoci */
      ppars[toPop].GAMMA[r] = ppars[fromPop].GAMMA[r];
      ppars[toPop].ProbPos[r] = ppars[fromPop].ProbPos[r];
      ppars[toPop].ProbDel[r] = ppars[fromPop].ProbDel[r];
      ppars[toPop].ProbNeut[r] = ppars[fromPop].ProbNeut[r];
      ppars[toPop].ProbNegGamma[r] = ppars[fromPop].ProbNegGamma[r];
      ppars[toPop].alphaN[r] = ppars[fromPop].alphaN[r];
      ppars[toPop].lambdaN[r] = ppars[fromPop].lambdaN[r];
      ppars[toPop].alphaP[r] = ppars[fromPop].alphaP[r];
      ppars[toPop].lambdaP[r] = ppars[fromPop].lambdaP[r];
      ppars[toPop].normMean[r] = ppars[fromPop].normMean[r];
      ppars[toPop].normVar[r] = ppars[fromPop].normVar[r];
      /* don't copy f0 */
    }
  }
  if(ppars[fromPop].neutpop != 0) /* non-neutral segregating variation */
    ppars[toPop].neutpop = ppars[fromPop].neutpop;
  /* don't copy RateClassSites, RateClassSites/loci, RateParamSites/loci */
  /* don't copy KAPPA or PSI */
  
  popn[toPop].maxchr = popn[fromPop].maxchr;
  popn[toPop].gen = popn[fromPop].gen;

  assert(popn[toPop].nextchr = malloc(gpars.R*sizeof(long)));
  for(i=0; i<gpars.R; i++)
    popn[toPop].nextchr[i] = 0;
  assert(popn[toPop].parentLoc = malloc(gpars.P*ppars[toPop].N*sizeof(long*)));
  for(i=0; i<gpars.P*ppars[toPop].N; i++){
    assert(popn[toPop].parentLoc[i] = malloc(gpars.R*sizeof(long)));
    memcpy(popn[toPop].parentLoc[i], popn[fromPop].parentLoc[i],
	   gpars.R*sizeof(long));
  }
  
  assert(popn[toPop].numCopies = malloc((popn[toPop].maxchr+1)*sizeof(long*)));
  assert(popn[toPop].chrGen = malloc((popn[toPop].maxchr+1)*sizeof(long*)));
  assert(popn[toPop].extremeMuts=malloc((popn[toPop].maxchr+1)*sizeof(long**)));
  assert(popn[toPop].polySites = malloc(gpars.R*sizeof(int*)));
  assert(popn[toPop].indfit = malloc((popn[toPop].maxchr+1)*sizeof(double*)));
  assert(popn[toPop].fixed = malloc(gpars.R*sizeof(struct history*)));
  assert(popn[toPop].conSeq = malloc(gpars.R*sizeof(char*))); /* [locus] */
  assert(popn[toPop].BigHead =
	 malloc((popn[toPop].maxchr+1)*sizeof(struct history**)));
  if(gpars.TRACKANC)
    assert(popn[toPop].ancestry = malloc((popn[toPop].maxchr+1)*sizeof(int**)));
  for(i=0; i<=popn[toPop].maxchr; i++){
    assert(popn[toPop].numCopies[i] = malloc(gpars.R*sizeof(long)));
    assert(popn[toPop].chrGen[i] = malloc(gpars.R*sizeof(long)));
    memcpy(popn[toPop].numCopies[i], popn[fromPop].numCopies[i],
	   gpars.R*sizeof(long));
    memcpy(popn[toPop].chrGen[i], popn[fromPop].chrGen[i],
	   gpars.R*sizeof(long));
    assert(popn[toPop].extremeMuts[i] = malloc(gpars.R*sizeof(long*)));
    assert(popn[toPop].indfit[i] = malloc(gpars.R*sizeof(double)));
    memcpy(popn[toPop].indfit[i], popn[fromPop].indfit[i], 
	   gpars.R*sizeof(double));
    assert(popn[toPop].BigHead[i] = malloc(gpars.R*sizeof(struct history*)));
    if(gpars.TRACKANC)
      assert(popn[toPop].ancestry[i] = malloc(gpars.R*sizeof(int*)));
    for(r=0; r<gpars.R; r++){
      if(i == 0){
	assert(popn[toPop].polySites[r] = malloc(gpars.L[r]*sizeof(int)));
	memcpy(popn[toPop].polySites[r], popn[fromPop].polySites[r],
	       gpars.L[r]*sizeof(int));
	assert(popn[toPop].conSeq[r] = malloc((gpars.L[r]+1)*sizeof(char)));
#ifdef WINDOWS
	strcpy_s(popn[toPop].conSeq[r], gpars.L[r]+1, popn[fromPop].conSeq[r]);
#else
	strcpy(popn[toPop].conSeq[r], popn[fromPop].conSeq[r]);
#endif
	popn[toPop].fixed[r] = popTrash();
	copyHistory(&popn[toPop].fixed[r], &popn[fromPop].fixed[r], r, toPop,
		    -1, fromPop, 0, 1, 1);
      }
      assert(popn[toPop].extremeMuts[i][r] = malloc(2*sizeof(long)));
      popn[toPop].extremeMuts[i][r][0] = popn[fromPop].extremeMuts[i][r][0];
      popn[toPop].extremeMuts[i][r][1] = popn[fromPop].extremeMuts[i][r][1];
      popn[toPop].BigHead[i][r] = popTrash();
      copyHistory(&popn[toPop].BigHead[i][r], &popn[fromPop].BigHead[i][r], r,
		  toPop, i, fromPop, 0, 0, 1);
      if(gpars.TRACKANC){
	assert(popn[toPop].ancestry[i][r] = malloc(gpars.L[r]*sizeof(int)));
	memcpy(popn[toPop].ancestry[i][r], popn[fromPop].ancestry[i][r],
	       gpars.L[r]*sizeof(int));
      }
    }
  }
  k = 0;
  for(r=0; r<gpars.R; r++){
    if(ppars[toPop].selDistType[r] != 0 || ppars[toPop].f0[r] < 1-FLT_EPSILON)
      k = 1;
  }
  if(k == 0){
    ppars[toPop].neutpop = 1;
  }
  if(gpars.substMod == 5){
    assert(popn[toPop].hpSUM = malloc((popn[toPop].maxchr+1)*sizeof(double**)));
    for(i=0; i<=popn[toPop].maxchr; i++){
      assert(popn[toPop].hpSUM[i] = malloc(gpars.R*sizeof(double*)));
      for(r=0; r<gpars.R; r++){
	assert(popn[toPop].hpSUM[i][r] = malloc(gpars.L[r]*sizeof(double)));
	memcpy(popn[toPop].hpSUM[i][r], popn[fromPop].hpSUM[i][r],
	       gpars.L[r]*sizeof(double));
      }
    }
  }
  
  if(ppars[toPop].pFEMALES != ppars[fromPop].pFEMALES){
    /* increase population size to twice difference in females */
    ChangePopSize(NULL, 2*fabs(ppars[toPop].pFEMALES-ppars[fromPop].pFEMALES)*
		  ppars[toPop].N+ppars[toPop].N, toPop, NULL, gpars.FITTEST);
    /* reduce size back down */
    ChangePopSize(NULL, ppars[fromPop].N, toPop, NULL, gpars.FITTEST);
  }
}

/* ------------------------------------------------------------------------- */

void admixture(struct EvolEvents *devents, const long PNanc, 
	       const int *selClass){
  long i, malesIN, femalesIN, cntF, cntM, newN;
  int toPop=devents->popi, fromPop=devents->ancPops[0];
  double *part=NULL, *relFit=NULL;
  
  /* first split ancPops[0] -> new population */
  Split(fromPop, toPop);
  
  /* subsample/expand as necessary to get right population size */
  if(devents->newP.Nt <= 0)
    newN = ppars[(int)(-devents->newP.Nt+.1)].N;
  else
    newN = (long)(devents->newP.Nt);
  
  if(gpars.FITTEST){
    part = dvector(0,ppars[toPop].N-1); /* fitness partition */
    relFit = dvector(0,ppars[toPop].N-1); /* relative fitnesses */
    partition(part, relFit, selClass, toPop);
  }
  ChangePopSize(relFit, newN, toPop, NULL, gpars.FITTEST);
  
  /* use migration technique to incorporate other populations */
  cntF = (long)(devents->femaleFreqs[0]*(gpars.P*ppars[toPop].MALES));
  cntM = (long)(devents->maleFreqs[0]*(gpars.P*(ppars[toPop].N-ppars[toPop].MALES)));
  for(i=1; i<devents->nAncPops; i++){
    femalesIN = (long)(devents->femaleFreqs[i]*(gpars.P*ppars[toPop].MALES));
    malesIN = (long)(devents->maleFreqs[i]*
		     (gpars.P*(ppars[toPop].N-ppars[toPop].MALES)));
    copyNchrs(femalesIN, malesIN, toPop, devents->ancPops[i], cntF, cntM);
    cntF += femalesIN;
    cntM += malesIN;
    if(ppars[i].neutpop == 0){
      ppars[toPop].neutpop = 0;
    }
  }
}

/* ------------------------------------------------------------------------- */

void Migrate(const int pop)
{
  long i, numIN;
  long malesIN;
  /* first determine total number of migrants entering from each population */
  for(i=0; i<gpars.NPOP; i++){
    if(i == pop || ppars[i].ALIVE == 0 || gpars.mig_mat[pop][i]<FLT_EPSILON)
      continue;
    /* draw actual number of migrants from poisson distribution */
    numIN = (int)(poidev(gpars.mig_mat[pop][i]/gpars.P,&gpars.seed));
    if(numIN < 1)  continue;
    malesIN = (int)(bnldev(ppars[i].pMaleMig, numIN, &gpars.seed)+.5);
    copyNchrs(numIN-malesIN, malesIN, pop, i, -1, -1);
  }
}

/* ------------------------------------------------------------------------- */

void ChangePopSize(double *relFit, const long newSize, const int pop, 
		   int *survive, const int fittest)
{
  int n;
  long i, j, p, r, newMALES, oldMALES, oldSize, rem=-1;
  float rn;

#ifdef UNIT_TEST
  long COUNT=0;
#endif
  if(newSize == ppars[pop].N)  return; /* nothing to do */

  /* set index of first male in population */
  newMALES = (long)(newSize*ppars[pop].pFEMALES);
  oldSize = ppars[pop].N;
  oldMALES = ppars[pop].MALES;
  ppars[pop].N = newSize;
  ppars[pop].MALES = newMALES;
  
#ifdef UNIT_TEST
  {
    fprintf(errfile,"pop %d: size change %ld -> %ld (#FEMALES %ld->%ld) with \
FITTEST=%d\n", pop,oldSize,newSize,oldMALES,newMALES,fittest);
    fflush(stdout);
  }
#endif /* UNIT_TEST */

  if(newSize == 0){
    if(ppars[pop].SS != 0){
      fprintf(errfile,"sfs_code error (%ld):  population %d size has shrunk \
to zero!\n",INITIALSEED,pop);
      abort();
    }
    else{
      ppars[pop].ALIVE = 0;
      if(survive != NULL)
	free(survive);
    }
  }
  if(newSize != ppars[pop].SS && (newMALES == 0 || newMALES >= newSize)){
    if(ppars[pop].SS != 0){
      fprintf(errfile,"sfs_code error(%ld): Male/female pop size went to zero!",
	      INITIALSEED);
      fprintf(errfile,"There are now %ld females and %ld males.\n", newMALES,
	      newSize-newMALES);
      abort();
    }
    else{
      ppars[pop].ALIVE = 0;
      if(survive != NULL)
	free(survive);
      return;
    }
  }
  
  if(newSize > oldSize){ /* deal with growth first */
    /* realloc parentLoc to accommodate extra individuals */
    assert(popn[pop].parentLoc = realloc(popn[pop].parentLoc,
					 gpars.P*newSize*sizeof(long*)));
    for(i=gpars.P*oldSize; i<gpars.P*newSize; i++)
      assert(popn[pop].parentLoc[i] = malloc(gpars.R*sizeof(long)));
    /* move first group of males to end of population to open up middle */
    /* just requires changing pointers, no frequency changes */
    for(i=0; i<newSize-oldSize && i<oldSize-oldMALES; i++){
      for(p=0; p<gpars.P; p++){
	for(r=0; r<gpars.R; r++){
	  popn[pop].parentLoc[gpars.P*(newSize-i-1)+p][r] =
	    popn[pop].parentLoc[gpars.P*(oldMALES+i)+p][r];
	}
      }
    }
    /* first expand females */
    for(i=oldMALES; i<newMALES; i++){
      for(p=0; p<gpars.P; p++){
	rn = ran1(&gpars.seed);
	if(fittest == 0 || ppars[pop].neutpop)
	  j = (long)(rn*oldMALES);
	else{
	  j = invCDF(rn*relFit[oldMALES-1], relFit, 0, oldMALES-1);
	}
	rn = ran1(&gpars.seed); /* choose chromosomes */
	n = gpars.P-1;
	while(n > 0 && rn <= (n+.0)/gpars.P)
	  n--;
	
	for(r=0; r<gpars.R; r++){
	  popn[pop].parentLoc[gpars.P*i+p][r] =
	    popn[pop].parentLoc[gpars.P*j+n][r];
	  popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][r]][r]++;
	  
#ifdef UNIT_TEST
	  if(p==0 && r==0)
	    COUNT++;
#endif
	}
      }
    }
    /* now expand males */
    for(i=newMALES; i<newMALES+((newSize-newMALES)-(oldSize-oldMALES)); i++){
      for(p=0; p<gpars.P; p++){
	rn = ran1(&gpars.seed);
	if(fittest == 0 || ppars[pop].neutpop)
	  j = (long)(oldSize-rn*(oldSize-oldMALES));
	else
	  j = invCDF(1.0-rn*(1.0-relFit[oldMALES-1]), relFit, oldMALES,
		     oldSize-1);
	rn = ran1(&gpars.seed); /* choose chromosomes */
	for(r=0; r<gpars.R; r++){
	  popn[pop].parentLoc[gpars.P*i+p][r] = 
	    popn[pop].parentLoc[gpars.P*j+p][r];
	  popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][r]][r]++;
	  
#ifdef UNIT_TEST
	  if(p==0 && r==0)
	    COUNT++;
#endif
	}
      }
    }
#ifdef UNIT_TEST
    if(COUNT != newSize-oldSize){
      fprintf(errfile,"sfs_code error (%ld):  COUNT != size difference \
(%ld != %ld)\n", INITIALSEED,COUNT,gpars.P*(newSize-oldSize));
      abort();
    }
#endif
  }
  else{ /* decreasing population size */
    if(survive == NULL){
      assert(survive = malloc(oldSize*sizeof(*survive)));
      for(i=0; i<oldSize; i++)  survive[i] = 0;
      /* determine which females survives */
      for(i=0; i<newMALES; i++){
	rn = ran1(&gpars.seed);
	if(fittest == 0 || ppars[pop].neutpop)
	  j = (long)(rn*oldMALES);
	else
	  j = invCDF(rn*relFit[oldMALES-1], relFit,0, oldMALES-1);
	if(survive[j] == 1){ /* choose w/o replacement */
	  i--;
	  continue;
	}
	survive[j] = 1;
      }
      /* determine which males survive */
      for(i=newMALES; i<newSize; i++){
	rn = ran1(&gpars.seed);
	if(fittest == 0 || ppars[pop].neutpop)
	  j = (long)(oldSize - rn*(oldSize-oldMALES));
	else
	  j = invCDF(1-rn*(1-relFit[oldMALES-1]), relFit, oldMALES, oldSize-1);
	if(survive[j] == 1){ /* choose w/o replacement */
	  i--;
	  continue;
	}
	survive[j] = 1;
      }
    }
    
#ifdef UNIT_TEST
    {
      j = 0;
      for(i=0; i<oldSize; i++)
	j += survive[i];
      if(j != newSize){
	fprintf(errfile,"error: survive does not sum up!!  %ld!=%ld\n",
		j,newSize);
	fflush(errfile);
      }
    }
#endif
    
    /* kill others */
    for(i=0; i<oldSize; i++){
      if(survive[i] == 0){
	if(rem == -1)
	  rem = i; /* first dead chrm */
	for(p=0; p<gpars.P; p++){
	  for(r=0; r<gpars.R; r++){
	    popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][r]][r]--;
	    if(popn[pop].numCopies[popn[pop].parentLoc[gpars.P*i+p][r]][r]==0){
	      popn[pop].chrGen[popn[pop].parentLoc[gpars.P*i+p][r]][r] =
		popn[pop].gen;
	      newDead(popn[pop].parentLoc[gpars.P*i+p][r], r);
	    }
	  }
	}
      }
      else if(i >= newMALES && i < oldMALES){ 
	/* move female to lower spot */
#ifdef UNIT_TEST
	if(rem == -1){
	  fprintf(errfile,"sfs_code internal error (%ld): rem == -1\n",
		  INITIALSEED);
	  abort();
	}
#endif
	if(rem >= oldMALES || survive[rem] == 1){
	  do{ /* find first available spot for female */
	    rem++;
	    if(rem >= oldMALES) rem = 0;
	  }while(survive[rem] == 1);
	}
	for(p=0; p<gpars.P; p++)
	  for(r=0; r<gpars.R; r++)
	    popn[pop].parentLoc[gpars.P*rem+p][r] = 
	      popn[pop].parentLoc[gpars.P*i+p][r];
	survive[rem] = 1; /* do not overwrite an individual */
	survive[i] = 0;   /* this spot is now free */
      }
      else if(i >= oldMALES){
	/* move males into the first spots after newMALES (no freq change) */
	do{ /* find first available spot for male */
	  rem++;
	  if(rem >= oldSize || rem < newMALES)
	    rem = newMALES;
	}while(survive[rem] == 1);
	for(p=0; p<gpars.P; p++)
	  for(r=0; r<gpars.R; r++)
	    popn[pop].parentLoc[gpars.P*rem+p][r] =
	      popn[pop].parentLoc[gpars.P*i+p][r];
	survive[rem] = 1; /* do not overwrite an individual */
	survive[i] = 0;   /* this spot is now free */
      }
    }
    
#ifdef UNIT_TEST
    printf("maxchr = %ld\n",popn[pop].maxchr);
    for(r=0; r<gpars.R; r++){
      j=0;
      for(i=0; i<=popn[pop].maxchr; i++){
#ifdef VERBOSE_DEBUG
	printf("nc[%ld] = %ld:  ",i,popn[pop].numCopies[i][r]);
	for(rem=0; rem<gpars.P*oldSize; rem++){
	  if(popn[pop].parentLoc[rem][r] == i && survive[rem/gpars.P] == 1)
	    printf("%ld ",rem/gpars.P);
	}
	printf("\n");
#endif
	if(popn[pop].numCopies[i][r] < 0){
	  printf("numCopies[%ld][%ld] = %ld < 0!!!\n",i,r,
		 popn[pop].numCopies[i][r]);
	  abort();
	}
	j += popn[pop].numCopies[i][r];
      }
      if(j != gpars.P*newSize){
	fprintf(errfile,"problem in ChangePopSize (%ld): newSize[%ld] = %ld, \
not %ld\n", INITIALSEED,r,j,gpars.P*newSize);
	abort();
      }
    }
    fflush(stdout);
#endif
    
    for(i=gpars.P*newSize; i<gpars.P*oldSize; i++){
      free(popn[pop].parentLoc[i]);
    }
    assert(popn[pop].parentLoc = 
	   realloc(popn[pop].parentLoc, gpars.P*newSize*sizeof(long*)));
    free(survive);
  }
}

/* ------------------------------------------------------------------------- */

void domesticate(int **DomestAllele, int *survive, int pop)
{
  int bestR=-1, bestSite=-1, GotIt=0;
  long i, j, k, m, r, numAlive=0, site, reg=0;
  float bestFreq=-1.0, tfreq;
  survive = ivector(0,ppars[pop].N-1);
  
  if(gpars.P != 2){
    fprintf(errfile,"sfs_code error(%ld): unfortunately, domestication has only\
 been implemented for the diploid case.\n",INITIALSEED);
    abort();
  }
  
#ifdef UNIT_TEST
  {
    int errstate=0;
    for(r=0; r<gpars.R; r++){
      unsigned long *count = lvector(0,gpars.P*ppars[pop].N);
      for(i=0; i<gpars.P*ppars[pop].N; i++){
	count[popn[pop].parentLoc[i][r]]++;
      }
      for(i=0; i<=popn[pop].maxchr; i++){
	if(count[i] != popn[pop].numCopies[i][r]){
	  fprintf(errfile, 
		  "error:  count[%ld] != numCopies[%d][%ld][%ld] \
(%ld vs %ld)\n", i, pop, i, r, count[i], popn[pop].numCopies[i][r]);
	  fflush(errfile);
	  errstate=1;
	}
      }
      free_lvector(count,0,gpars.P*ppars[pop].N);
    }
    if(errstate == 1){
      fprintf(errfile,"error at beginning (%ld)!!!\n",INITIALSEED);
      fflush(errfile);
      abort();
    }
  }
#endif /*UNIT_TEST*/
  
  /* pick the domestication allele */
  if(devents->locus == -1)
    devents->locus = (int)(gpars.R/2.0+.1);
  for(r=0; r<gpars.R; r++){
    for(i=0; i<2; i++){ /* reg = R/2 -/+ r */
      if((i == 0 && devents->locus < r) ||
	 (i == 1 && devents->locus + r >= gpars.R))
	continue;
      else if(i == 0)
	reg = devents->locus - r;
      else
	reg = devents->locus + r;
      
      for(j=0; j<mutAr[reg]->numMuts; j++){
	if(mutAr[reg]->muts[j]->event->free == '0'){
	  site = mutAr[reg]->muts[j]->event->site;
	  if(mutAr[reg]->muts[j]->event->numCarriers[pop] == 0 ||
	     DomestAllele[reg][site] == 1)
	    continue;
	  else{
	    tfreq = (float)mutAr[reg]->muts[j]->event->numCarriers[pop];
	    tfreq /= (gpars.P*ppars[pop].N + 0.0);
	    if(tfreq > devents->freq - devents->freq/10 &&
	       tfreq <= devents->freq + devents->freq/10){
	      bestR = reg;
	      bestSite = site;
	      bestFreq = tfreq;
	      GotIt = 1;
	      break;
	    }
	    else{
	      if(bestR == -1 || 
		 fabs(tfreq-devents->freq)<fabs(bestFreq-devents->freq)){
		bestR = reg;
		bestSite = site;
		bestFreq = tfreq;
	      }
	    }
	  }
	}
	if(GotIt) break;
      }
      if(GotIt) break;
    }
    if(GotIt) break;
  }
  
  if(!GotIt && bestR==-1){
    fprintf(errfile,"DID NOT FIND ADEQUATE MUTATION FOR DOMESTICATION! S=%ld\n",
	    INITIALSEED);
    for(r=0; r<gpars.R; r++){
      fprintf(errfile,"ERROR STATS FOR LOCUS %ld\n",r);
      for(i=0; i<mutAr[r]->numMuts; i++){
	fprintf(errfile,
		"%4ld: %4lu %5ld %5ld %c %c %c %2i %2i %1.2f %c %c %c [%c%c] %d"
		,i, mutAr[r]->muts[i]->event->site,
		mutAr[r]->muts[i]->event->gen,
		mutAr[r]->muts[i]->event->genFix[pop],
		mutAr[r]->muts[i]->event->ancNuc,
		mutAr[r]->muts[i]->event->derNuc,
		mutAr[r]->muts[i]->event->nonsyn,
		mutAr[r]->muts[i]->event->ancAA, 
		mutAr[r]->muts[i]->event->derAA,
		mutAr[r]->muts[i]->event->fit, mutAr[r]->muts[i]->event->CpG,
		mutAr[r]->muts[i]->event->fiveP,
		mutAr[r]->muts[i]->event->threeP,
		mutAr[r]->muts[i]->event->free,
		mutAr[r]->muts[i]->event->fixed[pop],
		popn[pop].polySites[r][mutAr[r]->muts[i]->event->site]);
	for(k=0; k<gpars.NPOP; k++){
	  fprintf(errfile," (p=%ld; mH=%ld); %ld: ", k,
		  mutAr[r]->muts[i]->event->maxHaps[k],
		  mutAr[r]->muts[i]->event->numCarriers[k]);
	  m = 0;
	  for(j=0; j<mutAr[r]->muts[i]->event->maxHaps[k]; j++){
	    if(mutAr[r]->muts[i]->event->hapFreq[k][j] != NULL){
	      fprintf(errfile,
		      "%ld.%ld ",j,*mutAr[r]->muts[i]->event->hapFreq[k][j]);
	      m++;
	    }
	  }
	}
	fprintf(errfile, "\n");
      }
      fprintf(errfile,"\n");
      fflush(errfile);
    }
    abort();
  }
  
#ifdef UNIT_TEST
  fprintf(errfile,"GotIt = %d\nbestR=%d; bestSite=%d; bestFreq=%f\n",GotIt,
	  bestR, bestSite, bestFreq);
#endif
  DomestAllele[bestR][bestSite] = 1;
  if(tmpOUT == NULL){
    fprintf(outfile,"DOM:%d,%d,%d,%d,%0.4f;",devents->popi,devents->popj,bestR,
	    bestSite,bestFreq);
    fflush(outfile);
  }
  else{
    fprintf(tmpOUT,"DOM:%d,%d,%d,%d,%0.4f;",devents->popi,devents->popj,bestR,
	    bestSite,bestFreq);
    fflush(tmpOUT);
  }
  /* determine who survives */
  for(i=0; i<ppars[pop].N; i++)
    survive[i]=0;
  numAlive = 0;
  /* first choose homozygous female carriers */
  for(i=0; i<ppars[pop].MALES; i++){
    if((retNucHistory(bestSite, 
		      &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i]
					 [bestR]][bestR], bestR, pop, 1) !=
	popn[pop].conSeq[bestR][bestSite]) &&
       (retNucHistory(bestSite,
		      &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i+1]
					 [bestR]][bestR], bestR, pop, 1) != 
	popn[pop].conSeq[bestR][bestSite])){
      survive[i] = 1;
      if((++numAlive) == (long)(ppars[pop].pFEMALES*devents->newP.Nt))
	break;
    }
  }
  /* then choose heterozygous female carriers */
  if(numAlive < (long)(ppars[pop].pFEMALES*devents->newP.Nt)){
    for(i=0; i<ppars[pop].MALES; i++){
      if(((retNucHistory(bestSite,
			 &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i]
					    [bestR]][bestR], bestR, pop, 1) != 
	   popn[pop].conSeq[bestR][bestSite]) ||
	  (retNucHistory(bestSite, 
			 &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i+1]
					    [bestR]][bestR], bestR, pop, 1) != 
	   popn[pop].conSeq[bestR][bestSite])) && survive[i] == 0){
	survive[i] = 1;
	if(++numAlive == (long)(ppars[pop].pFEMALES*devents->newP.Nt))
	  break;
      }
    }
    /* if insufficient carriers, choose other randomly */
    while(numAlive < (long)(ppars[pop].pFEMALES*devents->newP.Nt)){
      j = (long)(ran1(&gpars.seed)*ppars[pop].MALES);
      if(survive[j] == 0){ /* choose w/o replacement */
	survive[j] = 1;
	numAlive++;
      }
    }
  }
  /* now choose homozygous male carriers */
  for(i=ppars[pop].MALES; i<ppars[pop].N; i++){
    if((retNucHistory(bestSite, 
		      &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i]
					 [bestR]][bestR], bestR, pop, 1) != 
	popn[pop].conSeq[bestR][bestSite]) &&
       (retNucHistory(bestSite,
		      &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i+1]
					 [bestR]][bestR], bestR, pop, 1) != 
	popn[pop].conSeq[bestR][bestSite])){
      survive[i] = 1;
      if(++numAlive == (long)(devents->newP.Nt+.5))
	break;
    }
  }
  /* then choose heterozygous male carriers */
  if(numAlive < (long)(devents->newP.Nt+.5)){
    for(i=ppars[pop].MALES; i<ppars[pop].N; i++){
      if(((retNucHistory(bestSite,
			 &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i]
					    [bestR]][bestR], bestR, pop, 1) != 
	   popn[pop].conSeq[bestR][bestSite]) ||
	  (retNucHistory(bestSite,
			 &popn[pop].BigHead[popn[pop].parentLoc[gpars.P*i+1]
					    [bestR]][bestR], bestR, pop, 1) != 
	   popn[pop].conSeq[bestR][bestSite])) && survive[i] == 0){
	survive[i] = 1;
	if(++numAlive == (long)(devents->newP.Nt+.5))
	  break;
      }
    }
    /* if insufficient carriers, choose other randomly */
    while(numAlive < (long)(devents->newP.Nt+.5)){
      j = (long)((1.0-ran1(&gpars.seed)*
		  (1.-ppars[pop].pFEMALES))*ppars[pop].N);
      if(survive[j] == 0){ /* choose w/o replacement */
	survive[j] = 1;
	numAlive++;
      }
    }
  }
  /* use ChangePopSize to shrink population */  
  ChangePopSize(NULL, (long)(devents->newP.Nt+.5), pop, survive, gpars.FITTEST);
}

/* ------------------------------------------------------------------------- */

int NucToNum(const char *c)
{
  switch(*c){
  case '0': return(0);
  case '1': return(1);
  case '2': return(2);
  case '3': return(3);
  default:
    printf("error in NucToNum (%ld): %c\n",INITIALSEED,*c);
    abort();
    return(-1);
  }
}

/* ------------------------------------------------------------------------- */

char NumToNuc(const int c)
{
  switch(c){
  case 0: return('0');
  case 1: return('1');
  case 2: return('2');
  case 3: return('3');
  default:
    printf("error in NumToNuc (%ld): %i\n",INITIALSEED,c);
    abort();
    return(-1);
  }
}

/* ------------------------------------------------------------------------- */

char NumToAA(const int c)
{
  switch(c){
  case 0: return '*';  /* STOP CODON */
  case 1: return 'P';  /* PROLINE */
  case 2: return 'R';  /* ARGININE */
  case 3: return 'L';  /* LEUCINE */
  case 4: return 'H';  /* HISTIDINE */
  case 5: return 'Q';  /* GLUTAMINE */
  case 6: return 'A';  /* ALANINE */
  case 7: return 'G';  /* GLYCINE */
  case 8: return 'V';  /* VALINE */
  case 9: return 'D';  /* ASPARTIC ACID */
  case 10: return 'E'; /* GLUTAMIC ACID */
  case 11: return 'S'; /* SERINE */
  case 12: return 'C'; /* CYSTEINE */
  case 13: return 'W'; /* TRYPTOPHAN */
  case 14: return 'F'; /* PHENYLALANINE */
  case 15: return 'Y'; /* TYROSINE */
  case 16: return 'T'; /* THREONINE */
  case 17: return 'M'; /* METHIONINE (START) */
  case 18: return 'I'; /* ISOLEUCINE */
  case 19: return 'N'; /* ASPARAGINE */
  case 20: return 'K'; /* LYSINE */
  case 21: return 'X'; /* NON-CODING */
  default: 
    fprintf(stderr, "sfs_code error(%ld): %d is not a valid amino acid...\n",
	    INITIALSEED,c);
    abort();
    return 'Z';
  }
}

/* ------------------------------------------------------------------------- */

int AAToNum(const char c)
{
  switch(c){
  case '*': return 0;   /* STOP CODON */
  case 'P': return 1;   /* PROLINE */
  case 'R': return 2;   /* ARGININE */
  case 'L': return 3;   /* LEUCINE */
  case 'H': return 4;   /* HISTIDINE */
  case 'Q': return 5;   /* GLUTAMINE */
  case 'A': return 6;   /* ALANINE */
  case 'G': return 7;   /* GLYCINE */
  case 'V': return 8;   /* VALINE */
  case 'D': return 9;   /* ASPARTIC ACID */
  case 'E': return 10;  /* GLUTAMIC ACID */
  case 'S': return 11;  /* SERINE */
  case 'C': return 12;  /* CYSTEINE */
  case 'W': return 13;  /* TRYPTOPHAN */
  case 'F': return 14;  /* PHENYLALANINE */
  case 'Y': return 15;  /* TYROSINE */
  case 'T': return 16;  /* THREONINE */
  case 'M': return 17;  /* METHIONINE (START) */
  case 'I': return 18;  /* ISOLEUCINE */
  case 'N': return 19;  /* ASPARAGINE */
  case 'K': return 20;  /* LYSINE */
  case 'X': return 21;  /* NOTHING!! */
  default: 
    fprintf(stderr, "sfs_code error (%ld): %c is not a valid amino acid...\n",
	    INITIALSEED,c);
    abort();
    return -1;
  }
}

/* ------------------------------------------------------------------------- */

char ToCGTA(const char c)
{
  switch(c){
  case '0': return('C');
  case '1': return('G');
  case '2': return('T');
  case '3': return('A');
  case 'N': return('N');
  default:
    fprintf(errfile,"error in ToCGTA (%ld): %c\n",INITIALSEED,c);
    abort();
    return('N');
  }
}

/* ------------------------------------------------------------------------- */

char FromCGTA(const char c)
{
  switch(c){
  case 'C': return('0');
  case 'G': return('1');
  case 'T': return('2');
  case 'A': return('3');
  case 'c': return('0');
  case 'g': return('1');
  case 't': return('2');
  case 'a': return('3');
  default : return('N'); /* any non-nucleotide value */
  }
}

/* ------------------------------------------------------------------------- */

char TsNuc(const char c)
{
  if(c=='0')
    return('2');
  if(c=='1')
    return('3');
  if(c=='2')
    return('0');
  if(c=='3')
    return('1');
  else{
    fprintf(errfile,"error in TsNuc (%ld): %c\n",INITIALSEED,c);
    exit(-1);
  }
    
}

/* ------------------------------------------------------------------------- */

char TvNuc(const char c)
{
  float rn=ran1(&gpars.seed);

  if(c=='0'){
    if(rn<0.5)
      return('1');
    else return('3');
  }
  if(c=='1'){
    if(rn<0.5)
      return('0');
    else return('2');
  }
  if(c=='2'){
    if(rn<0.5)
      return('1');
    else return('3');
  }
  if(c=='3'){
    if(rn<0.5)
      return('0');
    else return('2');
  }
  else{
    fprintf(errfile,"error in TvNuc (%ld): %c\n",INITIALSEED,c);
    exit(-1);
  }
}

/* ------------------------------------------------------------------------- */

