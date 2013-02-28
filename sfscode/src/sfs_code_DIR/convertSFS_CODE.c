#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <assert.h>
#include <float.h>

#include "../util/myrand1.c"
#include "../util/nrutil.c"
#include "../util/nrutil.h"
/* #include "sfs_code.h" */
#include "convertSFS_CODE.h"

#define EPSILON 1e-8
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?	\
		   (iminarg1) : (iminarg2))

int main(int argc, char *argv[])
{
  int arg, j, ind;                    /* loop counters */
  long i, comLen, totLOCI=0;   /* data counters */
  int TOTPOPS=0, *INDperPOP=NULL, narg, ploid=2, linkTYPE=0;
  int *MALES=NULL, *NC, *SEX=NULL, numSEX=0;
  long seed, *maxSeqLen=NULL, numITS=0, numStored=0, Na=500, numLINES=0, totSS=0;
  double *LINK=NULL, maxGen=0.0;
  struct mutStorage **data=NULL, *tmpDATA;  /* data for all simulations */
  char c,c1,c2,c3,c4,c5,c6, args[100];
  char ***ancSeq=NULL;            /* all ancestral sequences */
  char **SFSCODE;                 /* command line to SFS_CODE */
  FILE *infile;

  if(argc == 1){
    fprintf(stderr,"usage:   %s <infile> [options]\nUse %s -h or %s --help\n",
	    argv[0],argv[0],argv[0]);
    exit(1);
  }
  else if(strcmp(argv[1],"--help") == 0 || strcmp(argv[1],"-h") == 0){
    helpConvert();
    return 0;
  }
  else if(!(infile = fopen(argv[1],"r"))){
    fprintf(stderr, "convertSFS_CODE error:  cannot read file %s\n",
	    argv[1]);
    exit(1);
  }
  
  assert(tmpDATA = malloc(sizeof(*tmpDATA)));
  assert(tmpDATA->pops = malloc(sizeof(*tmpDATA->pops)));
  assert(tmpDATA->chrs = malloc(sizeof(*tmpDATA->chrs)));
  assert(tmpDATA->event = malloc(sizeof(*tmpDATA->event)));
  tmpDATA->event->nSites = 0;
  tmpDATA->event->nucs = NULL;
  tmpDATA->Rtree = NULL;
  tmpDATA->Ltree = NULL;
  
  /* dynamically update the length of each sequence */
  comLen = 3; /* min. # words in SFS_CODE string */
  assert(SFSCODE = malloc(comLen*sizeof(*SFSCODE)));
  for(i=0; i<comLen; i++){
    assert(SFSCODE[i] = malloc(100*sizeof(**SFSCODE)));
  }
  /* read command-line to SFS_CODE */
  narg=0;
  i=0;
  while((c=fgetc(infile)) != EOF){
    if(c == '\n'){
      SFSCODE[narg][i] = '\0';
      if(strlen(SFSCODE[narg]) > 0)
	narg++;
      break;
    }
    else if(isspace(c)){
      SFSCODE[narg][i] = '\0';
      narg++;
      i=0;
    }
    else
      SFSCODE[narg][i++] = c;
    
    if(narg >= comLen){
      assert(SFSCODE = realloc(SFSCODE, (narg+1)*sizeof(*SFSCODE)));
      for(j=comLen; j<=narg; j++)
	assert(SFSCODE[j] = malloc(50*sizeof(**SFSCODE)));
      comLen = narg+1;
    }
  }

  /* parse sfs_code command line to get necessary parameters */
  TOTPOPS = atoi(SFSCODE[1]);
  /* initialize with default values */
  assert(INDperPOP = malloc(TOTPOPS*sizeof(*INDperPOP)));
  assert(MALES = malloc(TOTPOPS*sizeof(*MALES)));
  assert(NC = malloc(TOTPOPS*sizeof(*NC)));
  for(i=0; i<TOTPOPS; i++){
    INDperPOP[i] = 6;
    totSS += 6;
    MALES[i] = 0;
    NC[i] = 0;
  }
  totLOCI = 1;
  assert(maxSeqLen = malloc(totLOCI*sizeof(*maxSeqLen)));
  maxSeqLen[0] = 5002; 
  assert(ancSeq = malloc(sizeof(*ancSeq))); /* iterations */
  assert(ancSeq[0] = malloc(totLOCI*sizeof(**ancSeq))); /* 1st iteration*/
  for(i=0; i<totLOCI; i++)
    assert(ancSeq[0][i] = malloc((maxSeqLen[i]+1)*sizeof(***ancSeq)));
  
  arg = 3;
  while(arg < narg){
    if(SFSCODE[arg][0] == '-'){
      if(strcmp(SFSCODE[arg], "--annotate") == 0 || 
	 strcmp(SFSCODE[arg], "-a") == 0){ /* potentially get annotation file */
	arg++;
	if(strcmp(SFSCODE[arg], "F") == 0){
	  char c, *line;
	  FILE *annfile;
	  long maxLen = 1000;
	  arg++;
	  assert(annfile=fopen(SFSCODE[arg++],"r"));
	  assert(line = malloc(maxLen*sizeof(*line)));
	  if(fscanf(annfile, "%ld", &totLOCI) != 1){
	    printf("first line of annotation file must have the \
number of loci\n");
	    abort();
	  }
	  if(totLOCI > 1){
	    assert(maxSeqLen = realloc(maxSeqLen,totLOCI*sizeof(*maxSeqLen)));
	    assert(ancSeq[0] = realloc(ancSeq[0],totLOCI*sizeof(**ancSeq)));
	    assert(LINK = malloc((totLOCI-1)*sizeof(*LINK)));
	    LINK[0] = 0.0;
	    if(SEX == NULL)
	      assert(SEX = malloc((totLOCI)*sizeof(*SEX)));
	    SEX[0] = 0;
	  }
	  while((c=fgetc(annfile)) != EOF && (isspace(c) || c == ';'));
	  ungetc(c, annfile);
	  for(i=0; i<totLOCI; i++){
	    j = 0;
	    while((c=fgetc(annfile)) != EOF && c != ','){ /* length */
	      line[j] = c;
	      j++;
	    }
	    line[j] = '\0';
	    if(sscanf(line, "%ld", &maxSeqLen[i]) != 1){
	      printf("did not read length of locus %ld (%s)\n",i,line);
	      abort();
	    }
	    if(i == 0)
	      assert(ancSeq[0][0] =
		     realloc(ancSeq[0][0], (maxSeqLen[0]+1)*sizeof(***ancSeq)));
	    else
	      assert(ancSeq[0][i] = malloc((maxSeqLen[i]+1)*sizeof(***ancSeq)));
	    if(i<totLOCI-1)
	      LINK[i] = 0.0;
	    SEX[i] = 0;
	    while((c=fgetc(annfile)) != EOF && c != '\n');
	    /* skip to next line */
	  }
	  free(line);
	  fclose(annfile);
	}
      }
      else if(strcmp(SFSCODE[arg], "--ploidy") == 0 ||
	      strcmp(SFSCODE[arg], "-P") == 0){
	arg++;
	ploid = atoi(SFSCODE[arg++]);
      }
      else if(strcmp(SFSCODE[arg], "--length") == 0 ||
	      strcmp(SFSCODE[arg], "-L") == 0){
	arg++;
	totLOCI = atoi(SFSCODE[arg++]);
	if(totLOCI > 1){
	  assert(maxSeqLen = realloc(maxSeqLen,totLOCI*sizeof(*maxSeqLen)));
	  assert(ancSeq[0] = realloc(ancSeq[0],totLOCI*sizeof(**ancSeq)));
	  assert(LINK = malloc((totLOCI-1)*sizeof(*LINK)));
	  LINK[0] = 0.0;
	  assert(SEX = malloc((totLOCI)*sizeof(*SEX)));
	  SEX[0] = 0;
	}
	maxSeqLen[0] = atoi(SFSCODE[arg++]);
	if(maxSeqLen[0] != 5000){
	  assert(ancSeq[0][0] = 
		 realloc(ancSeq[0][0], (maxSeqLen[0]+3)*sizeof(***ancSeq)));
	}
	maxSeqLen[0] += 3;
	if(SFSCODE[arg][0] == '-'){
	  for(i=1; i<totLOCI; i++){
	    maxSeqLen[i] = maxSeqLen[0];
	    assert(ancSeq[0][i] = malloc(maxSeqLen[i]*sizeof(***ancSeq)));
	    SEX[i] = 0;
	  }
	}
	else{
	  for(i=1; i<totLOCI; i++){
	    if(SFSCODE[arg][0] == '-' || SFSCODE[arg][0] == 'R'){
	      for(j=i; j<totLOCI; j++){
		maxSeqLen[j] = maxSeqLen[j%i];
		assert(ancSeq[0][j] = malloc(maxSeqLen[j]*sizeof(***ancSeq)));
		if(j<totLOCI-1)
		  LINK[j] = 0.0;
		SEX[j] = 0;
	      }
	      if(SFSCODE[arg][0] == 'R')
		arg++;
	      break;
	    }
	    else{
	      maxSeqLen[i] = atoi(SFSCODE[arg++])+3;
	      assert(ancSeq[0][i] = malloc(maxSeqLen[i]*sizeof(***ancSeq)));
	      if(i<totLOCI-1)
		LINK[i] = 0.0;
	      SEX[i] = 0;
	    }
	  }
	}
      }
      else if(strcmp(SFSCODE[arg], "--LINK") == 0 ||
	      strcmp(SFSCODE[arg], "-l") == 0){
	arg++;
	if(SFSCODE[arg][0] == 'r')
	  linkTYPE = 1;
	else if(SFSCODE[arg][0] == 'p')
	  linkTYPE = 0;
	arg++;
	if(strcmp(SFSCODE[arg],"-1") != 0 && SFSCODE[arg][0] == '-'){
	  fprintf(stderr,"convertSFS_CODE error: sfs_code command line not \
read properly...  expecting linkage distance\n");
	  exit(1);
	}
	LINK[0] = atof(SFSCODE[arg++]);
	for(i=1; i<totLOCI-1; i++){
	  if(strcmp(SFSCODE[arg], "R") == 0 || 
	     (SFSCODE[arg][0] == '-' && strcmp(SFSCODE[arg],"-1") != 0)){
	    /* start repeating */
	    for(j=i; j<totLOCI-1; j++)
	      LINK[j] = LINK[j%i];
	    if(SFSCODE[arg][0] == 'R')
	      arg++;
	    break;;
	  }
	  else
	    LINK[i] = atof(argv[arg++]);
	}
      }
      else if(strcmp(SFSCODE[arg], "--popSize") == 0 ||
	      strcmp(SFSCODE[arg], "-N") == 0){
	arg++;
	ind = 1;
	if(strcmp(SFSCODE[arg], "P") == 0){
	  arg++;
	  if(atoi(SFSCODE[arg]) > 0){
	    arg += 2; /* skip non pop 0 declarations */
	    ind = 0;
	  }
	  else arg++;
	}
	if(ind == 1){
	  Na = atol(SFSCODE[arg]);
	  arg++;
	}
      }
      else if(strcmp(SFSCODE[arg], "--sex") == 0 || 
	      strcmp(SFSCODE[arg], "-x") == 0){
	arg++;
	printf("arg=%d: %s\n",arg, SFSCODE[arg]);
	if(SEX == NULL)
	  assert(SEX = malloc((totLOCI)*sizeof(*SEX)));
	SEX[0] = atoi(SFSCODE[arg]);
	numSEX += SEX[0];
	arg++;
	if(arg >= narg || SFSCODE[arg][0] == '-'){
	  for(i=1; i<totLOCI; i++){
	    SEX[i] = SEX[0];
	    numSEX += SEX[i];
	  }
	}
	else{
	  for(i=1; i<totLOCI; i++){
	    if(arg >= narg || SFSCODE[arg][0] == '-'){
	      fprintf(stderr,"convertSFS_CODE error:  did not read SFS_CODE \
command line properly.  Expecting %ld arguments to option --sex (-x)\n", 
		      totLOCI);
	    }
	    SEX[i] = atoi(SFSCODE[arg++]);
	    numSEX += SEX[i];
	  }
	}
      }
      else if(strcmp(SFSCODE[arg],"--sampSize") == 0 ||
	      strcmp(SFSCODE[arg], "-n") == 0){
	arg++;
	INDperPOP[0] = atoi(SFSCODE[arg++]);
	totSS = INDperPOP[0];
	if(arg >= narg || SFSCODE[arg][0] == '-'){
	  for(i=1; i<TOTPOPS; i++){
	    INDperPOP[i] = INDperPOP[0];
	    totSS += INDperPOP[i];
	  }
	}
	else{
	  for(i=1; i<TOTPOPS; i++){
	    INDperPOP[i] = atoi(SFSCODE[arg++]);
	    totSS += INDperPOP[i];
	  }
	}
      }
      else if(strcmp(SFSCODE[arg], "-TE") == 0){
	arg++;
	if(SFSCODE[arg][0] == 'P'){
	  arg += 2;
	}
	if(atof(SFSCODE[arg]) > maxGen){
	  maxGen = atof(SFSCODE[arg]);
	}
	arg++;
      }
      else{ /* skip unnecessary arguments */
	arg++;
      }
    }
    else{ /* skip unnecessary arguments */
      arg++;
    }
  }
  maxGen *= Na;
  totSS = 0;
  for(i=0; i<TOTPOPS; i++){
    INDperPOP[i] *= ploid;
    totSS += INDperPOP[i];
  }

  /* read seed */
  if(fscanf(infile, "SEED = %ld\n", &seed) != 1){
    fprintf(stderr,"convertSFS_CODE error: did not read seed on line %d\n", 2);
    exit(1);
  }
  
  /* first read in all data */
  numLINES = 2;
  while(1){ /* get all iterations in a loop */
    while((c=fgetc(infile)) != EOF){
      if(c == '/'){
	if((c=fgetc(infile)) == '/' && c != EOF){
	  numLINES++;
/* 	  printf("reading iteration %d\n",numITS+1); */
	  numITS++;
	  if(numITS == 1)
	    assert(data = malloc(sizeof(*data))); /* index for each iteration */
	  else{
	    assert(data = realloc(data, numITS*sizeof(*data)));
	    assert(ancSeq = realloc(ancSeq, numITS*sizeof(*ancSeq)));
	    assert(ancSeq[numITS-1] = malloc(totLOCI*sizeof(**ancSeq)));
	    for(i=0; i<totLOCI; i++){
	      assert(ancSeq[numITS-1][i] = 
		     malloc((maxSeqLen[i]+1)*sizeof(***ancSeq)));
	    }
	  }
	  assert(data[numITS-1] = malloc(sizeof(**data)));
	  data[numITS-1]->pops=NULL;
	  data[numITS-1]->chrs=NULL;
	  data[numITS-1]->event=NULL;
	  data[numITS-1]->Rtree=NULL;
	  data[numITS-1]->Ltree=NULL;
	  data[numITS-1]->Parent=NULL;
	  break;
	}
      }
    }
    if(c == EOF)
      break;
    if(fscanf(infile, "iteration:%ld/%d\n",&i, &j) != 2){
      fprintf(stderr, "convertSFS_CODE error:iteration %ld not read properly \
on line %ld\n", numITS, numLINES);
      exit(1);
    }
    while(1){ /* get ancestral sequences */
      if(fscanf(infile, ">locus_%d\n", &j) == 1){
	numLINES++;
	if(j >= totLOCI){
	  fprintf(stderr,"convertSFS_CODE error: expected %ld loci, read %d on \
line %ld\n", totLOCI, j+1, numLINES);
	  exit(1);
	}
	i=0;
	while((c = fgetc(infile)) != '\n' && c != EOF){
	  if(i >= maxSeqLen[j]){
	    fprintf(stderr, "convertSFS_CODE error: expected sequence length \
%ld, got something larger on line %ld!\n",maxSeqLen[j], numLINES);
	    exit(1);
	  }
	  ancSeq[numITS-1][j][i++] = c;
	}
	numLINES++;
	ancSeq[numITS-1][j][i] = '\0';
	/* keep allocated space as small as possible */
	maxSeqLen[j] = i;
      }
      else{
	break;
      }
    }

    /* skip over domestication allele info */
    int t1, t2, t3, t4;
    float t5;
    while(fscanf(infile, "DOM:%d,%d,%d,%d,%f;",&t1, &t2, &t3, &t4, &t5) == 5);
    i=0;
    if(TOTPOPS == 1){
      if(fscanf(infile,"Nc:%d;\n",&NC[i]) != 1){
	fprintf(stderr,"converSFS_CODE error: did not read Nc correctly on \
line %ld\n", numLINES);
	exit(1);
      }
    }
    else{
      if(fscanf(infile,"Nc:%d,",&NC[i]) != 1){
	fprintf(stderr,"converSFS_CODE error: did not read Nc correctly on \
line %ld\n", numLINES);
	exit(1);
      }
      i = 1;
      for(; i<TOTPOPS-1; i++){
	if(fscanf(infile,"%d,",&NC[i]) != 1){
	  fprintf(stderr,"converSFS_CODE error: read Nc incorrectly on line %ld\
.\n", numLINES);
	  exit(1);
	}
      }
      if(i < TOTPOPS){
	if(fscanf(infile,"%d;\n",&NC[i]) != 1){
	  fprintf(stderr,"converSFS_CODE error: read Nc incorrectly on line %ld\
.\n", numLINES);
	  exit(1);
	}
      }
    }
    numLINES++;
    i=0;
    if(TOTPOPS == 1){
      if(fscanf(infile,"MALES:%d;\n",&MALES[i]) != 1){
	fprintf(stderr,"converSFS_CODE error: did not read MALES correctly on \
line %ld\n", numLINES);
	exit(1);
      }
    }
    else{
      if(fscanf(infile,"MALES:%d,",&MALES[i]) != 1){
	fprintf(stderr,"converSFS_CODE error: did not read MALES correctly on \
line %ld\n", numLINES);
	exit(1);
      }
      i = 1;
      for(; i<TOTPOPS-1; i++){
	if(fscanf(infile,"%d,",&MALES[i]) != 1){
	  fprintf(stderr,"converSFS_CODE error: read MALES incorrectly on line \
%ld.\n", numLINES);
	  exit(1);
	}
      }
      if(i < TOTPOPS){
	if(fscanf(infile,"%d;\n",&MALES[i]) != 1){
	  fprintf(stderr,"converSFS_CODE error: read MALES incorrectly on line \
%ld.\n", numLINES);
	  exit(1);
	}
      }
    }
    numLINES++;
    for(i=0; i<TOTPOPS; i++)  MALES[i] *= ploid; /*index of first male chrom*/
    /* read stored mutation events */
    numStored = 0;
    while(1){
      if(fscanf(infile, "%ld,%c,%ld,%ld,%ld,%c%c%c,%c,%c",
		&tmpDATA->locus, &tmpDATA->event->axy, &tmpDATA->event->site, 
		&tmpDATA->event->gen, &tmpDATA->event->genFix, &c1, &c2, &c3,
		&c4, &tmpDATA->event->nonsyn) == 10){
	tmpDATA->event->fiveP = c1;
	tmpDATA->event->ancNuc = c2;
	tmpDATA->event->threeP = c3;
	tmpDATA->event->derNuc = c4;
      }
      else{
	break;
      }
      if(tmpDATA->event->nonsyn == '0' || tmpDATA->event->nonsyn == '1'){
	if(fscanf(infile, ",%c,%c,%lf,%ld", &c5, &c6, &tmpDATA->event->fit, 
		  &tmpDATA->numCarry) == 4){
	  tmpDATA->event->ancAA = AAToNum(c5);
	  tmpDATA->event->derAA = AAToNum(c6);
	}
	else{
	  fprintf(stderr,"error reading non-indel portion of event %ld\n",
		  numStored);
	  abort();
	}
      }
      else{
	if(fscanf(infile, ",%ld,", &tmpDATA->event->nSites) != 1){
	  fprintf(stderr,"error reading nSites for event %ld\n",numStored);
	  abort();
	}
	if(tmpDATA->event->nonsyn == 'i'){
	  assert(tmpDATA->event->nucs =
		 malloc((tmpDATA->event->nSites+1)*sizeof(char)));
	  for(i=0; i<tmpDATA->event->nSites; i++){
	    if((c=fgetc(infile)) != EOF && c != ',')
	      tmpDATA->event->nucs[i] = c;
	    else{
	      fprintf(stderr,"error reading nucs of event %ld\n",numStored);
	      abort();
	    }
	  }
	  tmpDATA->event->nucs[i] = '\0';
	}
	if(tmpDATA->event->nonsyn == 'i' && ((c=fgetc(infile)) == EOF ||
					     c != ',')){
	  fprintf(stderr,"error reading after nucs of event %ld\n",numStored);
	  abort();
	}
	if(fscanf(infile, "%lf,%ld", &tmpDATA->event->fit, 
		  &tmpDATA->numCarry) != 2){
	  fprintf(stderr,"error reading fitness/numCarry of event %ld\n",
		  numStored);
	  abort();
	}
      }
      assert(tmpDATA->pops =
	     realloc(tmpDATA->pops,
		     tmpDATA->numCarry*sizeof(*tmpDATA->pops)));
      assert(tmpDATA->chrs =
	     realloc(tmpDATA->chrs,
		     tmpDATA->numCarry*sizeof(*tmpDATA->chrs)));
      for(i=0; i<tmpDATA->numCarry; i++){
	if(fscanf(infile, ",%d.%ld",&tmpDATA->pops[i], &tmpDATA->chrs[i])!=2){
	  fprintf(stderr,"convertSFS_CODE error:  cannot read data iteration \
%ld after storing %ld items on line %ld\n", numITS, numStored, numLINES);
	  exit(1);
	}
	if(tmpDATA->pops[i] >= TOTPOPS){
	  fprintf(stderr,"convertSFS_CODE error: incorrect number of \
populations!  Expecting %d, found at least %d on line %ld\n",TOTPOPS,
		  tmpDATA->pops[i]+1, numLINES);
	  exit(1);
	}
	if(tmpDATA->chrs[i] != -1 &&
	   tmpDATA->chrs[i] >= INDperPOP[tmpDATA->pops[i]]){
	  fprintf(stderr,"convertSFS_CODE error:  incorrect number of \
individuals simulated in population %d.  Expecting %d, got at least %ld on \
line %ld\n", tmpDATA->pops[i], INDperPOP[tmpDATA->pops[i]], tmpDATA->chrs[i]+1,
		  numLINES);
	  exit(1);
	}
      }
      if((c=fgetc(infile)) != ';'){
	fprintf(stderr,"convertSFS_CODE error: cannot read data iteration \
%ld after storing %ld items on line %ld\n",numITS, numStored, numLINES);
      }
      buildStorage(&data[numITS-1]->Rtree, tmpDATA, NULL);
      numStored++;
      do{
	c=fgetc(infile);
	if(c == '\n')
	  numLINES++;
      }while(c != EOF && isspace(c));
      if(c == EOF)
	break;
      ungetc(c, infile);
    }
    if(c == EOF)
      break;
    /*     if(numITS == 100)  break; */
  }
  
  /* now parse command line and perform desired analysis */
  arg = 2;
  while(arg < argc){
    if(argv[arg][0] == '-'){
      strcpy(args, &argv[arg][1]);
      if(strcmp(args,"-help") == 0 || strcmp(args,"h") == 0){
	helpConvert();
	return 0;
      }
      else if(strcmp(args,"-alignment") == 0 || strcmp(args,"a") == 0){
	FILE *out=stdout;
	int numPOPS=-1, numINDIV=-1, numPI=-1, numLOCI=-1, printANC=0;
	int *pops=NULL, *indivs=NULL, *loci=NULL, **keep, *keepLOCI;
	int printITS=numITS;
	assert(keep = malloc(TOTPOPS*sizeof(*keep)));
	for(i=0; i<TOTPOPS; i++){
	  assert(keep[i] = malloc(INDperPOP[i]*sizeof(**keep)));
	  for(j=0; j<INDperPOP[i]; j++){
	    keep[i][j] = 0;
	  }
	}
	assert(keepLOCI = malloc(totLOCI*sizeof(*keepLOCI)));
	for(i=0; i<totLOCI; i++)
	  keepLOCI[i] = 0;
	arg++;
	while(arg < argc){
	  if(arg >= argc || argv[arg][0] == '-') break;
	  else{
	    if(strcmp(argv[arg], "A") == 0){
	      arg++;
	      printANC = 1; /* print ancestral seq */
	    }
	    else if(strcmp(argv[arg], "F") == 0){ /* define output file */
	      arg++;
	      assert(out = fopen(argv[arg],"w"));
	      arg++;
	    }
	    else if(strcmp(argv[arg], "L") == 0){ /* only print certain loci */
	      arg++;
	      numLOCI = atoi(argv[arg]);
	      arg++;
	      assert(loci = malloc(numLOCI*sizeof(*loci)));
	      for(i=0; i<numLOCI; i++){
		if(arg >= argc || argv[arg][0] == '-'){
		  fprintf(stderr,"converSFS_CODE error: not enough loci \
specified.  Expected %d, only got %ld.\n",numLOCI,i);
		  exit(1);
		}
		loci[i] = atoi(argv[arg]);
		arg++;
	      }
	    }
	    else if(strcmp(argv[arg], "I") == 0){ /* only print certain indivs*/
	      if(numPOPS > 0 || numPI > 0){
		fprintf(stderr, "convertSFS_CODE error in alignment options:\
use only one of \"P\", \"I\", or \"P.I\".\n");
		exit(1);
	      }
	      arg++;
	      numINDIV = atoi(argv[arg]);
	      arg++;
	      assert(indivs = malloc(numINDIV*sizeof(*indivs)));
	      for(i=0; i<numINDIV; i++){
		if(arg >= argc || argv[arg][0] == '-'){
		  fprintf(stderr,"converSFS_CODE error: not enough individuals\
specified.  Expected %d, only got %ld.\n",numINDIV,i);
		  exit(1);
		}
		indivs[i] = atoi(argv[arg]);
		for(j=0; j<TOTPOPS; j++){
		  if(INDperPOP[j] > indivs[i])
		    keep[j][indivs[i]] = 1;
		}
		arg++;
	      }
	    }
	    else if(strcmp(argv[arg], "P") == 0){ /* only print certain pops */
	      if(numINDIV > 0 || numPI > 0){
		fprintf(stderr, "convertSFS_CODE error in alignment options:\
use only one of \"P\", \"I\", or \"P.I\".\n");
		exit(1);
	      }
	      arg++;
	      numPOPS = atoi(argv[arg]);
	      arg++;
	      assert(pops = malloc(numPOPS*sizeof(*pops)));
	      for(i=0; i<numPOPS; i++){
		if(arg >= argc || argv[arg][0] == '-'){
		  fprintf(stderr,"converSFS_CODE error: not enough populations\
specified.  Expected %d, only got %ld.\n",numPOPS,i);
		  exit(1);
		}
		pops[i] = atoi(argv[arg]);
		for(j=0; j<INDperPOP[pops[i]]; j++)
		  keep[pops[i]][j] = 1;
		arg++;
	      }
	    }
	    else if(strcmp(argv[arg], "P.I") == 0){ /* spec. indivs from pops */
	      if(numPOPS > 0 || numINDIV > 0){
		fprintf(stderr, "convertSFS_CODE error in alignment options:\
use only one of \"P\", \"I\", or \"P.I\".\n");
		exit(1);
	      }
	      arg++;
	      numPI = atoi(argv[arg]);
	      arg++;
	      assert(pops = malloc(numPI*sizeof(*pops)));
	      assert(indivs = malloc(numPI*sizeof(*indivs)));
	      for(i=0; i<numPI; i++){
		if(arg >= argc || argv[arg][0] == '-'){
		  fprintf(stderr,"converSFS_CODE error: not enough population.\
individuals specified.  Expected %d, only got %ld.\n",numPI,i);
		  exit(1);
		}
		if(sscanf(argv[arg],"%d.%d",&pops[i], &indivs[i]) != 2){
		  fprintf(stderr,"convertSFS_CODE error: improper population.\
individual specified: %s\n",argv[arg]);
		  exit(1);
		}
		if(pops[i] < TOTPOPS && indivs[i] < INDperPOP[pops[i]])
		  keep[pops[i]][indivs[i]] = 1;
		else{
		  fprintf(stderr,"converSFS_CODE error: in alignment option P.I\
 combination %d.%d does not exist\n",pops[i],indivs[i]);
		  exit(1);
		}
		arg++;
	      }
	    }
	    else if(strcmp(argv[arg], "ITS") == 0){/* subset of alignments */
	      arg++;
	      printITS = atoi(argv[arg++]);
	      if(printITS > numITS)
		printITS = numITS;
	    }
	    else{
	      fprintf(stderr,"convertSFS_CODE error: %s not a defined option \
for printing alignments\n",argv[arg]);
	      exit(1);
	    }
	  }
	}
	printAlign(out, numPOPS, pops, numINDIV, indivs, numPI, numLOCI, loci,
		   keep, keepLOCI, printITS, maxSeqLen, numITS, totLOCI, 
		   TOTPOPS, INDperPOP, printANC, ancSeq, data);
	if(out != stdout) fclose(out);
      }
      else if(strcmp(args, "-AGE") == 0 || strcmp(args,"A") == 0){ /* AGE */
	int i, its2print = -1;
	char TYPE = '2';
	int pop=-1; /* default to all populations (-1) */
	FILE *out=stdout;

	/* read arguments to this option */
	arg++;
	while(arg < argc){
	  if(argv[arg][0] == '-')  break;
	  else if(strcmp(argv[arg], "F") == 0){ /* print to file */
	    arg++;
	    if(strcmp(argv[arg],"a") == 0){
	      arg++;
	      assert(out = fopen(argv[arg++],"a"));
	    }
	    else{
	      assert(out = fopen(argv[arg++],"w"));
	    }
	  }
	  else if(strcmp(argv[arg], "ITS") == 0){/*print subset of iterations*/
	    arg++;
	    if(arg < argc){
	      if(argv[arg][0] == '-'){
		its2print = numITS;
		if(strcmp(argv[arg], "-1") == 0)
		  arg++;
	      }
	      else{
		its2print = atoi(argv[arg]);
		arg++;
	      }
	    }
	    else
	      its2print = numITS;
	  }
	  else if(strcmp(argv[arg], "P") == 0){/*indicate single population*/
	    arg++;
	    if(arg >= argc){
	      fprintf(stderr,"error in AGE option, give population\n");
	      abort();
	    }
	    pop = atoi(argv[arg]);
	    if(pop < 0 || pop >= TOTPOPS){
	      fprintf(stderr,"invalid population specified: AGE option\n");
	      abort();
	    }
	    arg++;
	  }
	  else if(strcmp(argv[arg], "T") == 0){ /* mutant type */
	    arg++;
	    TYPE = argv[arg++][0];
	    if(TYPE != '0' && TYPE != '1' && TYPE != '2' &&
	       TYPE != 'i' && TYPE != 'd'){
	      fprintf(stderr,"convertSFS_CODE error: For --SFS, type must be \
0 (synonymous only), 1 (nonsynonymous only), or 2 (all mutants). Read %c.\n",
		      TYPE);
	      exit(1);
	    }
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error: %s is not a valid option \
for --SFS\n",argv[arg]);
	    exit(1);
	  }
	}
	
	for(i=0; i<numITS; i++){
	  if(its2print > 0 && i >= its2print)
	    break;
	  printAges(data[i]->Rtree, TYPE, pop, ploid*maxGen, ploid*Na, out, i);
	}
	if(out != stdout)
	  fclose(out);
      }
      else if(strcmp(args, "-derivedAlleles") == 0 || strcmp(args,"D") == 0){
	/* get number of derived alleles carried by each individual */
	int **cnt, i, j, k, *ss=NULL;
	char TYPE='2'; /* 0=>synon, 1=>nonsynon, 2=>all */
	int its2print=-1;
	int HH = 2; /* 0=heterozygous, 1=homozygous, 2=both */
	FILE *out=stdout;

	/* read arguments to this option */
	arg++;
	while(arg < argc){
	  if(argv[arg][0] == '-')  break;
	  else if(strcmp(argv[arg], "F") == 0){ /* print to file */
	    arg++;
	    if(strcmp(argv[arg],"a") == 0){
	      arg++;
	      assert(out = fopen(argv[arg++],"a"));
	    }
	    else{
	      assert(out = fopen(argv[arg++],"w"));
	    }
	  }
	  else if(strcmp(argv[arg], "H") == 0){ /* het/homo SNPs only */
	    arg++;
	    if(arg >= argc || (argv[arg][0] != '0' && argv[arg][0] != '1' &&
			       argv[arg][0] != '2')){
	      fprintf(stderr,"error, for hetero/homo SNPs in derivedAlleles \
option, must indicate 0, 1, or 2 after flag: H\n");
	      abort();
	    }
	    HH = atoi(argv[arg]);
	    if(HH != 2 && ploid != 2){
	      fprintf(stderr, "error, can only use het/homo SNPs in diploid \
case.\n");
	      abort();
	    }
	    arg++;
	  }
	  else if(strcmp(argv[arg], "ITS") == 0){/*print subset of iterations*/
	    arg++;
	    if(arg < argc){
	      if(argv[arg][0] == '-'){
		its2print = numITS;
		if(strcmp(argv[arg], "-1") == 0)
		  arg++;
	      }
	      else{
		its2print = atoi(argv[arg]);
		arg++;
	      }
	    }
	    else
	      its2print = numITS;
	  }
	  else if(strcmp(argv[arg], "T") == 0){ /* mutant type */
	    arg++;
	    TYPE = argv[arg++][0];
	    if(TYPE != '0' && TYPE != '1' && TYPE != '2' &&
	       TYPE != 'i' && TYPE != 'd'){
	      fprintf(stderr,"convertSFS_CODE error: For --SFS, type must be \
0 (synonymous only), 1 (nonsynonymous only), or 2 (all mutants). Read %c.\n",
		      TYPE);
	      exit(1);
	    }
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error: %s is not a valid option \
for --SFS\n",argv[arg]);
	    exit(1);
	  }
	}
	
	assert(cnt = malloc(TOTPOPS*sizeof(*cnt)));
	assert(ss = malloc(TOTPOPS*sizeof(*ss)));
	for(i=0; i<TOTPOPS; i++){
	  ss[i] = (int)(INDperPOP[i]/(ploid+0.0));
	  assert(cnt[i] = malloc(ss[i]*sizeof(**cnt)));
	  for(j=0; j<ss[i]; j++)
	    cnt[i][j] = 0;
	}
	for(i=0; i<numITS; i++){
	  if(its2print > 0 && i >= its2print)
	    break;
	  updateDerivedAlleles(&cnt, data[i]->Rtree, ss, TYPE, ploid, HH);
	  for(j=0; j<TOTPOPS; j++){
	    for(k = 0; k<ss[j]; k++){
	      fprintf(out,"%d ",cnt[j][k]);
	      cnt[j][k] = 0;
	    }
	    fprintf(out, "\t");
	  }
	  fprintf(out, "\n");
	}
	if(out != stdout)
	  fclose(out);
	for(i=0; i<TOTPOPS; i++)
	  free(cnt[i]);
	free(cnt);
	free(ss);
      }
      else if(strcmp(args, "-fitness") == 0 || strcmp(args,"f") == 0){
	/* get fitness of mutations subject to some criteria */
	int i;
	int its2print=-1;
	int pop = 0; /* for single population */
	int *share = NULL, nShare=0; /* only mutations shared across list */
	int private = 0;  /* only mutations private to given pop(s) */
	int fixed = 0;  /* 0=exclude fixed; 1=include fixed */
	FILE *out=stdout;
	
	/* read arguments to this option */
	arg++;
	while(arg < argc){
	  if(argv[arg][0] == '-')  break;
	  else if(strcmp(argv[arg], "f") == 0){
	    arg++;
	    fixed = 1;
	  }
	  else if(strcmp(argv[arg], "F") == 0){ /* print to file */
	    arg++;
	    if(strcmp(argv[arg],"a") == 0){
	      arg++;
	      assert(out = fopen(argv[arg++],"a"));
	    }
	    else{
	      assert(out = fopen(argv[arg++],"w"));
	    }
	  }
	  else if(strcmp(argv[arg], "ITS") == 0){/*print subset of iterations*/
	    arg++;
	    if(arg < argc){
	      if(argv[arg][0] == '-'){
		its2print = numITS;
		if(strcmp(argv[arg], "-1") == 0)
		  arg++;
	      }
	      else{
		its2print = atoi(argv[arg]);
		arg++;
	      }
	    }
	    else
	      its2print = numITS;
	  }
	  else if(strcmp(argv[arg], "p") == 0){ /* only private mutations*/
	    arg++;
	    private = 1;
	    if(nShare != 0){
	      fprintf(stderr,"cannot use private and shared option\n");
	      abort();
	    }
	  }
	  else if(strcmp(argv[arg], "P") == 0){/*indicate single population*/
	    arg++;
	    if(arg >= argc){
	      fprintf(stderr,"error in fitness option, give population\n");
	      abort();
	    }
	    pop = atoi(argv[arg]);
	    if(pop < 0 || pop >= TOTPOPS){
	      fprintf(stderr,"invalid population specified: fitness option\n");
	      abort();
	    }
	    arg++;
	  }
	  else if(strcmp(argv[arg], "S") == 0){/*indicate shared populations*/
	    pop = -1;
	    if(private == 1){
	      fprintf(stderr,"cannot use private and shared option\n");
	      abort();
	    }
	    arg++;
	    if(arg >= argc){
	      fprintf(stderr,"error in fitness option, give population\n");
	      abort();
	    }
	    nShare = atoi(argv[arg]); /* number of populations indicated */
	    arg++;
	    assert(share = malloc(nShare*sizeof(*share)));
	    for(i=0; i<nShare; i++){
	      if(arg >= argc){
		fprintf(stderr,"not enough populations specified!\n");
		abort();
	      }
	      share[i] = atoi(argv[arg]);
	      arg++;
	      if(share[i] < 0 || share[i] >= TOTPOPS){
		fprintf(stderr,"invalid pop specified: fitness option\n");
		abort();
	      }
	    }
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error: %s is not a valid option \
for --SFS\n",argv[arg]);
	    exit(1);
	  }
	}
	
	for(i=0; i<numITS; i++){
	  if(numITS > 0 && its2print > 0 && i>=its2print)
	    break;
	  printFitness(data[i]->Rtree, pop, nShare, share, private, fixed, out);
	}
	free(share);
	if(out != stdout)
	  fclose(out);
      }
      else if(strcmp(args,"-MK") == 0 || strcmp(args,"M") == 0){
	/* generate McDonald-Kreitman tables */
	int i, TRUE=1, ING, OUTG, numLOCI=0, *loci=NULL, its2print=-1, OGss=1;
	FILE *out=stdout;

	if(TOTPOPS == 1){
	  fprintf(stderr,"convertSFS_CODE error: cannot generate MK table from \
single population.\n");
	  exit(1);
	}
	arg++;
	if(arg > argc-2 || argv[arg][0] == '-' || argv[arg+1][0] == '-'){
	  fprintf(stderr,"convertSFS_CODE error: insufficient arguments to \
--MK (-m) option.  Please identify both ingoup and outgroup\n");
	  exit(1);
	}
	ING = atoi(argv[arg++]);
	OUTG = atoi(argv[arg++]);
	while(arg < argc){
	  if(argv[arg][0] == '-') break;
	  else if(strcmp(argv[arg], "F") == 0){
	    arg++;
	    if(strcmp(argv[arg],"a") == 0){
	      arg++;
	      assert(out = fopen(argv[arg++],"a"));
	    }
	    else{
	      assert(out = fopen(argv[arg++],"w"));
	    }
	  }
	  else if(strcmp(argv[arg], "ITS") == 0){
	    arg++;
	    if(arg < argc){
	      if(argv[arg][0] == '-'){
		its2print = numITS;
		if(strcmp(argv[arg], "-1") == 0)
		  arg++;
	      }
	      else{
		its2print = atoi(argv[arg]);
		if(its2print > 0)
		  arg++;
		else
		  its2print = numITS;
	      }
	    }
	    else
	      its2print = numITS;
	  }
	  else if(strcmp(argv[arg], "L") == 0){
	    arg++;
	    if(argv[arg][0] == '-' && strcmp(argv[arg], "-1") != 0){
	      fprintf(stderr,"convertSFS_CODE error: For option \'L\' in --MK \
(-m), must specify number of loci (or -1 for printing all independently)\n");
	      exit(1);
	    }
	    numLOCI = atoi(argv[arg++]);
	    if(numLOCI == -1){ /* print each locus independently */
	      numLOCI = totLOCI;
	      assert(loci = malloc(numLOCI*sizeof(*loci)));
	      for(i=0; i<numLOCI; i++)
		loci[i] = i;
	    }
	    else{
	      assert(loci = malloc(numLOCI*sizeof(*loci)));
	      for(i=0; i<numLOCI; i++){
		if(arg >= argc || argv[arg][0] == '-'){
		  if(i > 0){
		    fprintf(stderr,"convertSFS_CODE error: expecting %d loci \
only read %d\n",numLOCI,i);
		    exit(1);
		  }
		  else if(i == 0){
		    for(i=0; i<numLOCI; i++)
		      loci[i] = i;
		    break;
		  }
		}
		else
		  loci[i] = atoi(argv[arg++]);
	      }
	    }
	  }
	  else if(strcmp(argv[arg], "OGSS") == 0){
	    arg++;
	    OGss = atoi(argv[arg++]);
	    if(OGss == -1)
	      OGss = INDperPOP[OUTG];
	    else if(OGss > INDperPOP[OUTG]){
	      fprintf(stderr,"convertSFS_CODE error:  in --MK (-m), outgroup \
size exceeds outgroup sample size.  Max <OG_size> = %d\n",INDperPOP[OUTG]);
	      exit(1);
	    }
	  }
	  else if(strcmp(argv[arg], "OBS") == 0){
	    arg++;
	    TRUE = 0;
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error: for --MK (-m), %s option not\
 implemented.\n",argv[arg]);
	    exit(1);
	  }
	}
	
	printMK(TRUE, ING, OUTG, numLOCI, loci, its2print, out, data, numITS, 
		totLOCI, INDperPOP, ancSeq, OGss);

	if(out != stdout) fclose(out);
	if(loci != NULL)  free(loci);
      }
      else if(strcmp(argv[arg], "--ms") == 0 || 
	      strcmp(argv[arg], "-m") == 0){ /* produce ms-style ouput file */
	char TYPE = '2';
	char *filename=NULL;
	FILE *outfile=stdout;
	arg++;
	
	while(arg < argc){
	  if(strcmp(argv[arg], "F") == 0){
	    arg++;
	    filename = cvector(0,strlen(argv[arg]));
	    strcpy(filename, argv[arg++]);
	    assert(outfile = fopen(filename,"w"));
	  }
	  else if (strcmp(argv[arg], "T") == 0){
	    arg++;
	    TYPE = argv[arg++][0];
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error:  %s is not a valid option \
for --structure (-S)\n",argv[arg]);
	    exit(1);
	  }
	}
	fprintf(outfile, "sfscode2ms %ld %ld\n%ld %ld %ld\n\n", totSS, numITS, -seed, -seed,
		-seed);
	printMS(data, outfile, numITS, TOTPOPS, totLOCI, maxSeqLen, INDperPOP,
		seed, TYPE);
	if(filename != NULL)  free(filename);
	if(outfile != stdout) fclose(outfile);
      }
      else if(strcmp(args, "-privateShared") == 0 || strcmp(args,"p") == 0){
	/* generate private/shared polymorphism table */
	int i, *pops=NULL;
	int its2print=-1;
	char TYPE = '2';
	int *ps=NULL;
	FILE *out=stdout;
	
	if(TOTPOPS < 2){
	  fprintf(stderr,"error, --privateShared requires at least 2 \
populations.  Sorry!\n");
	  abort();
	}
	
	/* read arguments to this option */
	arg++;
	while(arg < argc){
	  if(argv[arg][0] == '-')  break;
	  else if(strcmp(argv[arg], "F") == 0){ /* print to file */
	    arg++;
	    if(strcmp(argv[arg],"a") == 0){
	      arg++;
	      assert(out = fopen(argv[arg++],"a"));
	    }
	    else{
	      assert(out = fopen(argv[arg++],"w"));
	    }
	  }
	  else if(strcmp(argv[arg], "ITS") == 0){/*print subset of iterations*/
	    arg++;
	    if(arg < argc){
	      if(argv[arg][0] == '-'){
		its2print = numITS;
		if(strcmp(argv[arg], "-1") == 0)
		  arg++;
	      }
	      else{
		its2print = atoi(argv[arg]);
		arg++;
	      }
	    }
	    else
	      its2print = numITS;
	  }
	  else if(strcmp(argv[arg], "P") == 0){ /* define populations */
	    assert(pops = malloc(2*sizeof(int)));
	    arg++;
	    pops[0] = atoi(argv[arg++]);
	    pops[1] = atoi(argv[arg++]);
	  }
	  else if(strcmp(argv[arg], "T") == 0){ /* mutant type */
	    arg++;
	    TYPE = argv[arg++][0];
	    if(TYPE != '0' && TYPE != '1' && TYPE != '2' &&
	       TYPE != 'i' && TYPE != 'd'){
	      fprintf(stderr,"convertSFS_CODE error: For --privateShared, \
type must be 0 (synonymous only), 1 (nonsynonymous only), or 2 (all mutants). \
Read %c.\n", TYPE);
	      abort();
	    }
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error: %s is not a valid option \
for --SFS\n",argv[arg]);
	    exit(1);
	  }
	}
	
	if(pops == NULL){
	  assert(pops = malloc(2*sizeof(*pops)));
	  pops[0] = 0;
	  pops[1] = 1;
	}

	assert(ps = malloc(3*sizeof(*ps)));
	for(i=0; i<3; i++) 
	  ps[i] = 0; 
	for(i=0; i<numITS; i++){
	  if(its2print > 0 && numITS > 0 && i >= its2print)
	    break;
	  updatePrivShared(data[i]->Rtree, &ps, TYPE, pops);
	}
	
	fprintf(out,"%d\t%d\t%d\n",ps[0], ps[1], ps[2]);
	free(ps);
	free(pops);
	if(out != stdout)
	  fclose(out);
      }
      else if(strcmp(args, "-SFS") == 0 || strcmp(args,"S") == 0){ /* SFS */
	int i, j, PAUTO=0, PSEX=0, MONLY=0;
	char TYPE = '2'; /* 0=>synon, 1=>nonsynon, 2=>all */
	int TRUE = 1, its2print=-1, numPOPS=-1, numLOCI=0, *numOG=NULL;
	int *pops=NULL, *loci=NULL, **og=NULL, OGsize=1, **keep = NULL, tmp,tp1;
	FILE *out=stdout;

	assert(keep = malloc(TOTPOPS*sizeof(*keep)));
	for(i=0; i<TOTPOPS; i++){
	  assert(keep[i] = malloc(INDperPOP[i]*sizeof(**keep)));
	  for(j=0; j<INDperPOP[i]; j++)
	    keep[i][j] = 1;
	}

	arg++;
	/* read arguments to this option */
	while(arg < argc){
	  if(argv[arg][0] == '-')  break;
	  else if(strcmp(argv[arg], "A") == 0){
	    /* print only autosomal loci */
	    arg++;
	    PAUTO = 1;
	    if(PSEX == 1){
	      fprintf(stderr,"convertSFS_CODE error:  Use either \'A\' or \'X\'\
, but not both.\n");
	      exit(1);
	    }
	    numLOCI = 0;
	    for(i=0; i<totLOCI; i++){
	      if(SEX == NULL || !SEX[i]){
		assert(loci = realloc(loci, (numLOCI+1)*sizeof(*loci)));
		loci[numLOCI++] = i;
	      }
	    }
	    if(numLOCI == 0){
	      fprintf(stderr,"convertSFS_CODE error: no autosomal loci.\n");
	      exit(1);
	    }
	  }
	  else if(strcmp(argv[arg], "F") == 0){ /* print to file */
	    arg++;
	    if(strcmp(argv[arg],"a") == 0){
	      arg++;
	      assert(out = fopen(argv[arg++],"a"));
	    }
	    else{
	      assert(out = fopen(argv[arg++],"w"));
	    }
	  }
	  else if(strcmp(argv[arg], "ITS") == 0){/*print subset of iterations*/
	    arg++;
	    if(arg < argc){
	      if(argv[arg][0] == '-'){
		its2print = numITS;
		if(strcmp(argv[arg], "-1") == 0)
		  arg++;
	      }
	      else{
		its2print = atoi(argv[arg]);
		arg++;
	      }
	    }
	    else
	      its2print = numITS;
	  }
	  else if(strcmp(argv[arg],"L") == 0){
	    arg++;
	    numLOCI = atoi(argv[arg++]);
	    if(numLOCI == -1){
	      numLOCI = totLOCI;
	      assert(loci = malloc(numLOCI*sizeof(*loci)));
	      for(i=0; i<totLOCI; i++)
		loci[i] = i;
	    }
	    else{
	      assert(loci = malloc(numLOCI*sizeof(*loci)));
	      for(i=0; i<numLOCI; i++){
		if(argv[arg][0] == '-'){
		  fprintf(stderr, "convertSFS_CODE error: not enough arguments \
to \'L\' option of --SFS (-s).  Expected %d, got %d\n",numLOCI,i);
		  exit(1);
		}
		loci[i] = atoi(argv[arg++]);
	      }
	    }
	  }
	  else if(strcmp(argv[arg], "N") == 0){ /* random subsample */
	    arg++;
	    tmp = 0;
	    for(i=0; i<TOTPOPS; i++){
	      if(isdigit(argv[arg][0]))
		tmp = atoi(argv[arg++]);
	      for(j=0; j<INDperPOP[i]; j++)  keep[i][j] = 0;
	      j = tmp;
	      if(j > INDperPOP[i]){
		fprintf(stderr,"convertSFS_CODE error: requensted number of \
individuals for population %d too large.  Maximum sample size = %d, requested \
%d.\n",i,INDperPOP[i],tmp);
		exit(1);
	      }
	      while(j > 0){
		tp1 = (int)(INDperPOP[i]*ran1(&seed));
		if(keep[i][tp1] == 0){
		  keep[i][tp1] = 1;
		  j--;
		}
	      }
	    }
	  }
	  else if(strcmp(argv[arg], "OBS") == 0){ /* print observed SFS */
	    if(TRUE && numPOPS != -1){
	      fprintf(stderr,"convertSFS_CODE error: for --SFS (-s), define \
populations using \'OBS\' option, rather than \'P\' option.");
	      exit(1);
	    }
	    arg++;
	    TRUE = 0;
	    if(numPOPS == -1){
	      numPOPS = 1;
	      assert(pops = malloc(sizeof(*pops)));
	      assert(og = malloc(TOTPOPS*sizeof(*og)));
	      assert(numOG = malloc(TOTPOPS*sizeof(*numOG)));
	      for(i=0; i<TOTPOPS; i++){
		og[i] = NULL;
		numOG[i] = 0;
	      }
	    }
	    else{
	      numPOPS++;
	      assert(pops = realloc(pops, numPOPS*sizeof(*pops)));
	    }
	    pops[numPOPS-1] = atoi(argv[arg++]);
	    numOG[pops[numPOPS-1]]++;
	    if(numOG[pops[numPOPS-1]] == 1){
	      assert(og[pops[numPOPS-1]] = malloc(sizeof(**og)));
	    }
	    else{
	      assert(og[pops[numPOPS-1]] = 
		     realloc(og[pops[numPOPS-1]], 
			     numOG[pops[numPOPS-1]]*sizeof(**og)));
	    }
	    og[pops[numPOPS-1]][numOG[pops[numPOPS-1]]-1] = atoi(argv[arg++]);
	    if(arg < argc){
	      if(argv[arg][0] != '-'){
		OGsize = atoi(argv[arg]);
		if(OGsize > 0 && 
		   OGsize <= og[pops[numPOPS-1]][numOG[pops[numPOPS-1]]-1])
		  arg++;
		else
		  OGsize = 1;
	      }
	      else if(strcmp(argv[arg], "-1") == 0){
		OGsize=INDperPOP[og[pops[numPOPS-1]][numOG[pops[numPOPS-1]]-1]];
		arg++;
	      }
	    }
	  }
	  else if(strcmp(argv[arg],"P") == 0){
	    arg++;
	    if(!TRUE){
	      fprintf(stderr,"convertSFS_CODE error: for --SFS (-s), define \
populations using \'OBS\' option, rather than \'P\' option.");
	      exit(1);
	    }
	    numPOPS = atoi(argv[arg++]);
	    if(numPOPS <= 0 || numPOPS > TOTPOPS){
	      fprintf(stderr,"convertSFS_CODE error: for --SFS (-s), number \
of populations must be >0 and <=%d (total number of populations simulated).\n",
		      TOTPOPS);
	      exit(1);		      
	    }
	    assert(pops = malloc(numPOPS*sizeof(*pops)));
	    for(i=0; i<numPOPS; i++){
	      if(argv[arg][0] == '-'){
		fprintf(stderr, "convertSFS_CODE error: not enough arguments \
to \'P\' option of --SFS (-s).  Expected %d, got %d\n",numPOPS,i);
		exit(1);
	      }
	      pops[i] = atoi(argv[arg++]);
	    }
	  }
	  else if(strcmp(argv[arg], "T") == 0){ /* mutant type */
	    arg++;
	    TYPE = argv[arg++][0];
	    if(TYPE != '0' && TYPE != '1' && TYPE != '2' &&
	       TYPE != 'i' && TYPE != 'd'){
	      fprintf(stderr,"convertSFS_CODE error: For --SFS, type must be \
0 (synonymous only), 1 (nonsynonymous only), 2 (all mutants), \'i\' (insertion)\
 or \'d\' (deletion). Read %c.\n",
		      TYPE);
	      exit(1);
	    }
	  }
	  else if(strcmp(argv[arg], "X") == 0){
	    /* print only sex loci */
	    arg++;
	    PSEX = 1;
	    if(PAUTO){
	      fprintf(stderr,"convertSFS_CODE error:  Use either \'A\' or \'X\'\
, not both.\n");
	      exit(1);
	    }
	    numLOCI = 0;
	    for(i=0; i<totLOCI; i++){
	      if(SEX[i]){
		assert(loci = realloc(loci, (numLOCI+1)*sizeof(*loci)));
		loci[numLOCI++] = i;
	      }
	    }
	    if(numLOCI == 0){
	      fprintf(stderr,"convertSFS_CODE error:  no sex loci found\n");
	      exit(1);
	    }
	  }
	  else if(strcmp(argv[arg], "Y") == 0){
	    /* print only sex loci */
	    arg++;
	    PSEX = 1;
	    MONLY = 1;
	    if(PAUTO){
	      fprintf(stderr,"convertSFS_CODE error:  Use either \'A\' or \'X\'\
, not both.\n");
	      exit(1);
	    }
	    numLOCI = 0;
	    for(i=0; i<totLOCI; i++){
	      if(SEX[i]){
		assert(loci = realloc(loci, (numLOCI+1)*sizeof(*loci)));
		loci[numLOCI++] = i;
	      }
	    }
	    if(numLOCI == 0){
	      fprintf(stderr,"convertSFS_CODE error:  no sex loci found\n");
	      exit(1);
	    }
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error: %s is not a valid option \
for --SFS\n",argv[arg]);
	    exit(1);
	  }
	}
	
	if(numPOPS == -1){
	  numPOPS = TOTPOPS;
	  assert(pops = malloc(numPOPS*sizeof(*pops)));
	  for(i=0; i<numPOPS; i++)
	    pops[i] = i;
	}
	if(numLOCI == 0){
	  assert(loci = malloc(sizeof(*loci)));
	  assert(loci = realloc(loci, totLOCI*sizeof(*loci)));
	  for(i=0; i<totLOCI; i++)
	    loci[i] = i;
	}
	if(numOG == NULL){
	  assert(numOG = malloc(TOTPOPS*sizeof(*numOG)));
	  for(i=0; i<TOTPOPS; i++)
	    numOG[i] = 0;
	}

	printSFS(TYPE, TRUE, its2print, numPOPS, numLOCI, numOG, pops, loci, og,
		 out, data, TOTPOPS, totLOCI, INDperPOP, numITS, ancSeq, OGsize,
		 PSEX, MALES, MONLY, ploid, keep);
	if(out != stdout) fclose(out);
	for(i=0; i<TOTPOPS; i++)
	  free(keep[i]);
	free(keep);
      }
      else if(strcmp(argv[arg], "--structure") == 0 || 
	      strcmp(argv[arg], "-s") == 0){ /* produce structure input file */
	int numPOPS=0, *pops=NULL, numINDIV=0, *indivs=NULL, numPI=0;
	int **keep;
	char *filename=NULL;
	double CMMB = 1.0;
	arg++;

	assert(keep = malloc(TOTPOPS*sizeof(*keep)));
	for(i=0; i<TOTPOPS; i++){
	  assert(keep[i] = malloc(INDperPOP[i]*sizeof(**keep)));
	  for(j=0; j<INDperPOP[i]; j++){
	    keep[i][j] = 0;
	  }
	}

	while(arg < argc){
	  if(strcmp(argv[arg], "F") == 0){
	    arg++;
	    filename = cvector(0,strlen(argv[arg]));
	    strcpy(filename, argv[arg++]);
	  }
	  else if(strcmp(argv[arg], "CMMB") == 0){
	    arg++;
	    CMMB = atof(argv[arg++]);
	    if(CMMB < 0){
	      fprintf(stderr,"convertSFS_CODE error:  centiMorgans/Megabase \
< 0!!\n");
	      exit(1);
	    }
	  }
	  else if(strcmp(argv[arg], "I") == 0){ /* only print certain indivs*/
	    if(numPOPS > 0 || numPI > 0){
	      fprintf(stderr, "convertSFS_CODE error in structure options:\
use only one of \"P\", \"I\", or \"P.I\".\n");
	      exit(1);
	    }
	    arg++;
	    numINDIV = atoi(argv[arg]);
	    arg++;
	    assert(indivs = malloc(numINDIV*sizeof(*indivs)));
	    for(i=0; i<numINDIV; i++){
	      if(arg >= argc || argv[arg][0] == '-'){
		fprintf(stderr,"converSFS_CODE error: not enough individuals\
specified.  Expected %d, only got %ld.\n",numINDIV,i);
		exit(1);
	      }
	      indivs[i] = atoi(argv[arg]);
	      for(j=0; j<TOTPOPS; j++){
		if(INDperPOP[j] > indivs[i])
		  keep[j][indivs[i]] = 1;
	      }
	      arg++;
	    }
	  }
	  else if(strcmp(argv[arg], "P") == 0){ /* only print certain pops */
	    if(numINDIV > 0 || numPI > 0){
	      fprintf(stderr, "convertSFS_CODE error in structure options:\
use only one of \"P\", \"I\", or \"P.I\".\n");
	      exit(1);
	    }
	    arg++;
	    numPOPS = atoi(argv[arg]);
	    arg++;
	    assert(pops = malloc(numPOPS*sizeof(*pops)));
	    for(i=0; i<numPOPS; i++){
	      if(arg >= argc || argv[arg][0] == '-'){
		fprintf(stderr,"converSFS_CODE error: not enough populations\
specified.  Expected %d, only got %ld.\n",numPOPS,i);
		exit(1);
	      }
	      pops[i] = atoi(argv[arg]);
	      for(j=0; j<INDperPOP[pops[i]]; j++)
		keep[pops[i]][j] = 1;
	      arg++;
	    }
	  }
	  else if(strcmp(argv[arg], "P.I") == 0){ /* spec. indivs from pops */
	    if(numPOPS > 0 || numINDIV > 0){
	      fprintf(stderr, "convertSFS_CODE error in structure options:\
use only one of \"P\", \"I\", or \"P.I\".\n");
	      exit(1);
	    }
	    arg++;
	    numPI = atoi(argv[arg]);
	    arg++;
	    assert(pops = malloc(numPI*sizeof(*pops)));
	    assert(indivs = malloc(numPI*sizeof(*indivs)));
	    for(i=0; i<numPI; i++){
	      if(arg >= argc || argv[arg][0] == '-'){
		fprintf(stderr,"converSFS_CODE error: not enough population.\
individuals specified.  Expected %d, only got %ld.\n",numPI,i);
		exit(1);
	      }
	      if(sscanf(argv[arg],"%d.%d",&pops[i], &indivs[i]) != 2){
		fprintf(stderr,"convertSFS_CODE error: improper population.\
individual specified: %s\n",argv[arg]);
		exit(1);
	      }
	      if(pops[i] < TOTPOPS && indivs[i] < INDperPOP[pops[i]])
		keep[pops[i]][indivs[i]] = 1;
	      else{
		fprintf(stderr,"converSFS_CODE error: in alignment option P.I\
 combination %d.%d does not exist\n",pops[i],indivs[i]);
		exit(1);
	      }
	      arg++;
	    }
	  }
	  else{
	    fprintf(stderr,"convertSFS_CODE error:  %s is not a valid option \
for --structure (-S)\n",argv[arg]);
	    exit(1);
	  }
	}

       if(numPOPS == 0 && numINDIV == 0 && numPI == 0){
	  for(i=0; i<TOTPOPS; i++){
	    for(j=0; j<INDperPOP[i]; j++)
	      keep[i][j] = 1;
	  }
	}
	printSTRUCTURE(data, ancSeq, filename, numITS, TOTPOPS, totLOCI, LINK,
		       CMMB, linkTYPE, maxSeqLen, keep, INDperPOP);
	if(filename != NULL)  free(filename);
      }
      else{
	fprintf(stderr, "convertSFS_CODE error:  %s has not yet been \
implemented\n",argv[arg]);
 	exit(1);
      }
    }
    else{
      fprintf(stderr, "convertSFS_CODE error at argument %d: %s\n",
 	      arg,argv[arg]);
      exit(1);
    }
  }
 
  free(tmpDATA->pops); 
  free(tmpDATA->chrs); 
  if(tmpDATA->event->nucs != NULL)
    free(tmpDATA->event->nucs);
  free(tmpDATA->event); 
  free(tmpDATA); 
  free(MALES);

  for(i=0; i<numITS; i++){ 
    for(j=0; j<totLOCI; j++){
      free(ancSeq[i][j]);
    }
    freeStorage(data[i]);
    free(ancSeq[i]);
  }
  free(data);
  free(ancSeq);
  for(i=0; i<comLen; i++)
    free(SFSCODE[i]);
  free(SFSCODE);
  free(maxSeqLen);
  free(INDperPOP);
  fclose(infile);
  free(NC);
  free(SEX);
  free(LINK);
  return(0);
}

/* ------------------------------------------------------------------------- */

void printAlign(FILE *out, int numPOPS, int *pops, int numINDIV, int *indivs,
		int numPI, int numLOCI, int *loci, int **keep, int *keepLOCI,
		int printITS, long *maxSeqLen, long numITS, int totLOCI,
		int TOTPOPS, int *INDperPOP, int printANC, char ***ancSeq,
		struct mutStorage **data)
{
  int i, j, k, it;
  char ****seqs=NULL, **ancPopSeq=NULL;
  /* allocate space for sequences */
  int maxLOCUS = 0, maxPOP=0, *maxIND;
  long numInv=0, **invBK=NULL;
  assert(maxIND = malloc(TOTPOPS*sizeof(*maxIND)));
  for(i=0; i<TOTPOPS; i++)
    maxIND[i] = -1;
  
  if(numPOPS == -1 && numINDIV == -1 && numPI == -1){
    maxPOP = TOTPOPS-1;
    for(i=0; i<TOTPOPS; i++){
      maxIND[i] = INDperPOP[i]-1;
      for(j=0; j<INDperPOP[i]; j++)
	keep[i][j] = 1;
    }
  }
  else{
    for(i=0; i<TOTPOPS; i++){
      maxIND[i] = 0;
      for(j=0; j<INDperPOP[i]; j++){
	if(keep[i][j] == 1){
	  if(i > maxPOP)
	    maxPOP = i;
	  if(j > maxIND[i])
	    maxIND[i] = j;
	}
      }
    }
  }
  if(numLOCI > 0){
    for(i=0; i<numLOCI; i++){
      if(loci[i] > maxLOCUS)
	maxLOCUS = loci[i];
      keepLOCI[loci[i]] = 1;
    }
  }
  else{
    for(i=0; i<totLOCI; i++)
      keepLOCI[i] = 1;
  }
  assert(seqs = malloc(TOTPOPS*sizeof(*seqs)));
  for(i=0; i<TOTPOPS; i++)
    seqs[i] = NULL;
  for(i=0; i<=maxPOP; i++){
    if(maxIND[i] >= 0){
      assert(seqs[i] = malloc(INDperPOP[i]*sizeof(**seqs)));
      for(j=0; j<INDperPOP[i]; j++){
	seqs[i][j] = NULL;
	if(keep[i][j] == 1){
	  assert(seqs[i][j] = malloc(totLOCI*sizeof(***seqs)));
	  for(k=0; k<totLOCI; k++){
	    seqs[i][j][k] = NULL;
	    if(keepLOCI[k] == 1){
	      assert(seqs[i][j][k] = 
		     malloc((maxSeqLen[k]+1)*sizeof(****seqs)));
	    }
	  }
	}
      }
    }
  }
  if(printANC){
    assert(ancPopSeq = malloc(totLOCI*sizeof(*ancPopSeq)));
    for(i=0; i<totLOCI; i++){
      ancPopSeq[i] = NULL;
      if(keepLOCI[i] == 1)
	assert(ancPopSeq[i] = 
	       malloc((maxSeqLen[i]+1)*sizeof(**ancPopSeq)));
    }
  }
  for(it=0; it<printITS; it++){
    for(k=0; k<totLOCI; k++){
      if(keepLOCI[k] == 1){
	for(i=0; i<TOTPOPS; i++){
	  int printIT = 1;
	  for(j=0; j<INDperPOP[i]; j++){
	    if(keep[i][j] == 1){
	      if(printIT && printANC){ /* take care of ancestral seq */
		strcpy(ancPopSeq[k], ancSeq[it][k]);
		numInv = 0;
		updateSeq(&ancPopSeq[k], data[it]->Rtree, i, -1, k,
			  &numInv, &invBK);
		if(invBK != NULL){
		  free(invBK);
		  invBK = NULL;
		}
		fprintf(out,">it%dpop%dindAlocus%d\n",it,i,k);
		fprintf(out,"%s\n",ancPopSeq[k]);
		printIT = 0;
	      }
	      strcpy(seqs[i][j][k], ancSeq[it][k]);
	      numInv = 0;
	      updateSeq(&seqs[i][j][k], data[it]->Rtree, i, j, k,
			&numInv, &invBK);
	      if(invBK != NULL){
		free(invBK);
		invBK = NULL;
	      }
	      fprintf(out,">it%dpop%dind%dlocus%d\n",it,i,j,k);
	      fprintf(out,"%s\n",seqs[i][j][k]);
	    }
	  }
	}
      }
    }
  }
  
  if(printANC){
    for(i=0; i<totLOCI; i++)
      if(ancPopSeq[i] != NULL)
	free(ancPopSeq[i]);
    free(ancPopSeq);
  }
  for(i=0; i<TOTPOPS; i++){
    if(seqs[i] != NULL){
      for(j=0; j<INDperPOP[i]; j++)
	if(seqs[i][j] != NULL){
	  for(k=0; k<totLOCI; k++)
	    if(seqs[i][j][k] != NULL)
	      free(seqs[i][j][k]);
	  free(seqs[i][j]);
	}
      free(seqs[i]);
    }
    free(keep[i]);
  }
  free(seqs);
  free(keep);
  free(keepLOCI);
  free(maxIND);
  
  if(pops != NULL)  free(pops);
  if(indivs != NULL) free(indivs);
  if(loci != NULL)  free(loci);
}

/* ------------------------------------------------------------------------- */

void printAges(struct mutStorage *storage, char TYPE, int pop, double maxGen,
	       long PNa, FILE *out, long IT)
{
  int i, k = 0;
  if(storage == NULL)
    return;
  else if(TYPE == '2' || TYPE == storage->event->nonsyn){
    for(i=0; i<storage->numCarry; i++){
      if(pop == -1 && storage->chrs[i] != -1){
	k = 1;
	break;
      }
      else if(pop != -1 && storage->pops[i] == pop && storage->chrs[i] != -1){
	k = 1;
	break;
      }
    }
    if(k == 1){
      if(fabs(storage->event->fit) < 1e-9)
	fprintf(out,"%ld\t%ld\t0.0\t%f\n", IT, storage->event->site,
		(maxGen - storage->event->gen)/PNa);
      else
	fprintf(out,"%ld\t%ld\t%1.9f\t%f\n",IT, storage->event->site,storage->event->fit,
		(maxGen - storage->event->gen)/PNa);
    }
  }
  
  printAges(storage->Rtree, TYPE, pop, maxGen, PNa, out, IT);
  printAges(storage->Ltree, TYPE, pop, maxGen, PNa, out, IT);
}

/* ------------------------------------------------------------------------- */

void printSFS(char TYPE, int TRUE, int its2print, int numPOPS, int numLOCI,
	      int *numOG, int *pops, int *loci, int **og, FILE *out, 
	      struct mutStorage **data, int TOTPOPS, int totLOCI, 
	      int *INDperPOP, long numITS, char ***ancSeq, int OGsize, int PSEX,
	      int *MALES, int MONLY, int ploid, int **keep)
{
  int i, j, k, m, n, *ind;
  int ****sfs=NULL, *tmpSFS=NULL, maxSS=0;
  
  /* allocate necessary memory */
  /* create 2 sets of SFS: tmp (refreshed for every iter/locus/pop/OG) 
     and a permanent one that collects the temporary SFS as desired */
  assert(sfs = malloc(numPOPS*sizeof(*sfs)));/* [pop][loci][OG][freq] */
  assert(ind = malloc(numPOPS*sizeof(*ind)));
  for(i=0; i<numPOPS; i++){
    ind[i] = 0;
    for(j=0; j<INDperPOP[pops[i]]; j++)
      ind[i] += keep[pops[i]][j];
    if(numLOCI == 0)
      assert(sfs[i] = malloc(sizeof(**sfs))); /* sum over loci */
    else
      assert(sfs[i] = malloc(numLOCI*sizeof(**sfs)));
    for(j=0; j==0 || j<numLOCI; j++){
      if(numOG[i] == 0)
	assert(sfs[i][j] = malloc(sizeof(***sfs))); /* OG irrelevant */
      else
	assert(sfs[i][j] = malloc(numOG[i]*sizeof(***sfs)));
      for(k=0; k==0 || k<numOG[pops[i]]; k++){
	assert(sfs[i][j][k] = malloc(ind[i]*sizeof(****sfs)));
	for(m=0; m<ind[i]; m++)
	  sfs[i][j][k][m] = 0;
      }
      if(tmpSFS == NULL){
	maxSS = ind[i];
	assert(tmpSFS = malloc(maxSS*sizeof(*tmpSFS)));
      }
      else if(ind[i] > maxSS){
	maxSS = ind[i];
	assert(tmpSFS = realloc(tmpSFS, maxSS*sizeof(*tmpSFS)));
      }
    }
  }

  for(i=0; i<numITS; i++){
    if(its2print > 0 && i >= its2print)
      break;
    for(j=0; j<numPOPS; j++){
      for(m=0; m==0 || m<numOG[pops[j]]; m++){
	for(k=0; k<totLOCI; k++){
	  if(numLOCI != 0 && k >= numLOCI)  break;
	  for(n=1; n<ind[j]; n++)
	    tmpSFS[n] = 0;
	  if(TRUE)
	    updateTrueSFS(tmpSFS, data[i]->Rtree, pops[j], ind[j],
			  loci[k], TYPE, PSEX, MALES, MONLY, ploid, keep);
	  else
	    updateObsSFS(tmpSFS, &data[i], pops[j], ind[j],
			 loci[k], TYPE, og[pops[j]][m], ancSeq[i][loci[k]],
			 OGsize, PSEX, MALES, MONLY, ploid, keep, INDperPOP[j]);
	  if(numLOCI == 0){ /* sum over loci */
	    for(n=1; n<ind[j]; n++){
	      sfs[j][0][m][n] += tmpSFS[n];
	      tmpSFS[n] = 0;
	    }
	  }
	  else if(its2print == -1){
	    for(n=1; n<ind[j]; n++){
	      sfs[j][k][m][n] += tmpSFS[n];
	      tmpSFS[n] = 0;
	    }
	  }
	  else{
	    for(n=1; n<ind[j]; n++){
	      fprintf(out,"%d ",tmpSFS[n]);
	      tmpSFS[n] = 0;
	    }
	    fprintf(out,"\t");
	  }
	}
	if(numLOCI == 0 && its2print != -1){ /*sum over loci, print each iter */
	  for(n=1; n<ind[j]; n++){
	    fprintf(out,"%d ",sfs[j][0][m][n]);
	    sfs[j][0][m][n] = 0;
	  }
	  fprintf(out,"\t");
	}
      }
    }
    if(its2print != -1)
      fprintf(out,"\n");
  }
  if(its2print == -1){
    for(j=0; j<numPOPS; j++){
      for(m=0; m==0 || m<numOG[pops[j]]; m++){
	for(k=0; k==0 || k<numLOCI; k++){
	  for(n=1; n<ind[j]; n++)
	    fprintf(out,"%d ",sfs[j][k][0][n]);
	  fprintf(out,"\t");
	}
	fprintf(out,"\t");
      }
    }
    fprintf(out,"\n");
  }

  /* free memory */
  for(j=0; j<numPOPS; j++){
    for(k=0; k==0 || k<numLOCI; k++){
      for(m=0; m==0 || m<numOG[j]; m++)
	free(sfs[j][k][m]);
      free(sfs[j][k]);
    }
    free(sfs[j]);
  }
  free(sfs);
  free(tmpSFS);
  free(ind);
  if(loci != NULL)  free(loci);
  if(pops != NULL)  free(pops);
  for(i=0; i<TOTPOPS; i++)
    if(numOG[i] > 0)  free(og[i]);
  free(numOG);
  free(og);
}

/* ------------------------------------------------------------------------- */

void printMK(int TRUE, int ING, int OUTG, int numLOCI, int *loci, int its2print,
	     FILE *out, struct mutStorage **data, long numITS, int totLOCI,
	     int *INDperPOP, char ***ancSeq, int OGss)
{
  int i, j, k, **MK=NULL, *tmpMK=NULL;;
  if(numLOCI == 0){ /* sum across loci */
    assert(loci = malloc(totLOCI*sizeof(*loci)));
    for(i=0; i<totLOCI; i++)
      loci[i] = i;
    assert(MK = malloc(sizeof(*MK)));
    assert(MK[0] = malloc(4*sizeof(**MK)));
    for(i=0; i<4; i++)
      MK[0][i] = 0;
  }
  else{
    assert(MK = malloc(numLOCI*sizeof(*MK)));
    for(i=0; i<numLOCI; i++){
      assert(MK[i] = malloc(4*sizeof(**MK)));
      for(j=0; j<4; j++)
	MK[i][j] = 0;
    }
  }
  assert(tmpMK = malloc(4*sizeof(*tmpMK)));
  for(i=0; i<4; i++) tmpMK[i] = 0;
  for(i=0; i<numITS; i++){
/*     if(i < 2) continue; */
/*     printf("calculating table for iteration %d\n",i); */
/*     fflush(stdout); */
    for(j=0; j<totLOCI; j++){
      if(numLOCI != 0 && j >= numLOCI) break;
      if(TRUE)
	getTrueMK(tmpMK, ING, OUTG, data[i]->Rtree, loci[j], data[i]->Rtree,
		  INDperPOP[ING], OGss);
      else /* generate observed MK table */
	getObsMK(tmpMK, ING, OUTG, &data[i], loci[j], INDperPOP[ING], 
		 OGss, ancSeq[i][loci[j]]);
      
      if(numLOCI == 0){ /* sum across loci */
	for(k=0; k<4; k++){
	  MK[0][k] += tmpMK[k];
	  tmpMK[k] = 0;
	}
      }
      else{
	for(k=0; k<4; k++){
	  MK[loci[j]][k] += tmpMK[k];
	  tmpMK[k] = 0;
	}
      }
    }
    if(its2print != -1){
      for(j=0; j==0 || j<numLOCI; j++){
	for(k=0; k<4; k++){
	  fprintf(out,"%d ",MK[j][k]);
	}
	fprintf(out, "%1.4f\t",checkSig(MK[j][0]+.0, MK[j][1]+.0, 
					   MK[j][2]+.0,MK[j][3]+.0));
	MK[j][0] = MK[j][1] = MK[j][2] = MK[j][3] = 0;
      }
      fprintf(out, "\n");
    }
  }
  if(its2print == -1){
    for(j=0; j==0 || j<numLOCI; j++){
      for(k=0; k<4; k++){
	fprintf(out,"%d ",MK[j][k]);
      }
      fprintf(out, "%1.4f\t",checkSig(MK[j][0]+.0, MK[j][1]+.0, 
					 MK[j][2]+.0,MK[j][3]+.0));
      MK[j][0] = MK[j][1] = MK[j][2] = MK[j][3] = 0;
    }
    fprintf(out, "\n");
  }
  for(j=0; j==0 || j<numLOCI; j++)
    free(MK[j]);
  free(MK);
  free(tmpMK);
}

/* ------------------------------------------------------------------------- */

void printMS(struct mutStorage **data, FILE *outfile, long numITS, int TOTPOPS,
	     int totLOCI, long *maxSeqLen, int *INDperPOP, long seed, char TYPE)
{
  int i, chrs=0, *sumChrs=NULL;
  long Nsites=-1, it, *totLen=NULL, sumLen=0, maxMuts=0;
  char **ZeroOne=NULL;
  double *positions=NULL;
  
  assert(sumChrs = malloc(TOTPOPS*sizeof(*sumChrs)));
  sumChrs[0] = 0;
  chrs = INDperPOP[0];
  for(i=1; i<TOTPOPS; i++){
    sumChrs[i] = sumChrs[i-1]+INDperPOP[i-1];
    chrs += INDperPOP[i];
  }

  assert(positions=malloc(sizeof(*positions)));
  assert(ZeroOne = malloc(chrs*sizeof(*ZeroOne)));
  for(i=0; i<chrs; i++){
    assert(ZeroOne[i] = malloc(sizeof(**ZeroOne)));
    ZeroOne[i][0] = '\0';
  }
  
  assert(totLen = malloc(totLOCI*sizeof(*totLen)));
  for(i=0; i<totLOCI; i++){
    totLen[i] = sumLen;
    sumLen += maxSeqLen[i];
  }

  for(it=0; it<numITS; it++){
    for(i=0; i<chrs; i++)
      ZeroOne[i][0] = '\0';
    Nsites = -1;
    
    getZeroOne(data[it]->Rtree, &ZeroOne, &positions, &Nsites, sumChrs, chrs,
	       totLen, sumLen, TYPE, &maxMuts, INDperPOP);
    Nsites++; /* was zero based, now is 1 based */
    fprintf(outfile,"//\nsegsites: %ld\n",Nsites);
    if(Nsites > 0){
      fprintf(outfile,"positions: ");
      for(i=0; i<Nsites; i++)
	fprintf(outfile,"%f ",positions[i]);
      fprintf(outfile,"\n");
      for(i=0; i<chrs; i++)
	fprintf(outfile,"%s\n",ZeroOne[i]);
      fprintf(outfile, "\n");
    }
  }

  free(positions);
  free(sumChrs);
  free(totLen);
  for(i=0; i<chrs; i++)
    free(ZeroOne[i]);
  free(ZeroOne);
}

/* ------------------------------------------------------------------------- */

void getZeroOne(struct mutStorage *data, char ***ZeroOne, double **positions,
		long *Nsites, int *popInds, int totchrs, long *totLen, 
		long sumLen, char TYPE, long *maxMuts, int *INDperPOP)
{
  long i, j, pos=0;
  if(data != NULL){
    getZeroOne(data->Ltree, ZeroOne, positions, Nsites, popInds, totchrs, 
	       totLen, sumLen, TYPE, maxMuts, INDperPOP);
    if(TYPE=='2' || data->event->nonsyn==TYPE){
      for(i=0; i<(*Nsites); i++){
	if((*positions)[i] == data->event->site){
	  pos = i;
	  break;
	}
      }
      if(i >= *Nsites){
	(*Nsites)++;
	pos = *Nsites;
	if((*Nsites) > (*maxMuts)){
	  (*maxMuts)++;
	  assert((*positions) = 
		 realloc(*positions, (*Nsites+1)*sizeof(**positions)));
	  for(i=0; i<totchrs; i++){
	    assert((*ZeroOne)[i] = 
		   realloc((*ZeroOne)[i], (*Nsites+2)*sizeof(***ZeroOne)));
	    (*ZeroOne)[i][*Nsites] = '0';
	    (*ZeroOne)[i][*Nsites+1] = '\0';
	  }
	}
	(*positions)[*Nsites] = 
	  (totLen[data->locus]+data->event->site)/(sumLen+.0);
	for(i=0; i<totchrs; i++){
	  (*ZeroOne)[i][*Nsites] = '0';
	  (*ZeroOne)[i][*Nsites+1] = '\0';
	}
      }
      for(i=0; i<data->numCarry; i++){
	if(data->chrs[i] == -1){
	  for(j=0; j<INDperPOP[data->pops[i]]; j++){
	    (*ZeroOne)[popInds[data->pops[i]]+j][pos] = '1';
	  }
	}
	else
	  (*ZeroOne)[popInds[data->pops[i]]+data->chrs[i]][pos] = '1';
      }
    }
    getZeroOne(data->Rtree, ZeroOne, positions, Nsites, popInds, totchrs, 
	       totLen, sumLen, TYPE, maxMuts, INDperPOP);
  }
}

/* ------------------------------------------------------------------------- */

void printSTRUCTURE(struct mutStorage **data, char ***ancSeq, char *filename,
		    long numITS, int TOTPOPS, int totLOCI, double *LINK, 
		    double CMMB, int linkTYPE, long *maxSeqLen, int **keep,
		    int *INDperPOP)
{
  int i, j, k, m, n, *Nsites=NULL;
  long **sites=NULL, psite=0;
  char str[100], c='x';
  FILE *out=stdout;

  for(i=0; i<numITS; i++){
    if(filename != NULL){
      if(numITS > 1){
	sprintf(str,"%d.%s",i,filename);
	assert(out = fopen(str,"w"));
      }
      else
	assert(out = fopen(filename, "w"));
    }
    
    assert(sites = malloc(totLOCI*sizeof(*sites)));
    assert(Nsites = malloc(totLOCI*sizeof(*Nsites)));
    for(j=0; j<totLOCI; j++){
      sites[j] = NULL;
      Nsites[j] = 0;
    }
    for(j=0; j<totLOCI; j++){
      getPolyFixedSites(data[i]->Rtree, j, &Nsites[j], &sites[j]);
    }
    for(j=0; j<totLOCI; j++){
      for(k=0; k<Nsites[j]; k++)
	fprintf(out,"s.%d.%d ",j,k);
    }
    fprintf(out, "\n");
    for(j=0; j<totLOCI; j++){
      for(k=0; k<Nsites[j]; k++){
	if(j == 0 && k == 0){
	  fprintf(out, "-1 ");
	  psite = 0;
	}
	else if(k > 0){
	  psite += (sites[j][k]-sites[j][k-1]);
	  fprintf(out, "%ld ",psite);
	  if(k == Nsites[j]-1)
	    psite += (maxSeqLen[j]-sites[j][k]);
	}
	else{ /* j > 0 && k == 0 */
	  if(fabs(LINK[j-1]) < EPSILON){ /* completely linked loci */
	    psite += sites[j][k];
	    fprintf(out, "%ld ", psite);
	  }
	  else{
	    if(linkTYPE == 0){ /* physical distance */
	      if(LINK[j-1] < 0){ /* independence */
		psite = 0;
		fprintf(out, "-1 ");
	      }
	      else{
		psite += LINK[j-1]+sites[j][k];
		fprintf(out, "%ld ",psite);
	      }
	    }
	    else{
	      /* convert recombination distance using Haldane's mapping func */
	      psite += (long)((-log(1-2*LINK[j-1])/2/CMMB*1e6)+sites[j][k]);
	      fprintf(out, "%ld ",psite);
	    }
	  }
	}
      }
    }
    fprintf(out,"\n");
    for(j=0; j<TOTPOPS; j++){
      for(k=0; k<INDperPOP[j]; k++){
	if(keep[j][k] == 1){
	  fprintf(out, "I.%d.%d ",j,(int)(k/2));
	  fprintf(out, "%d ",j);
	  for(m=0; m<totLOCI; m++){
	    for(n=0; n<Nsites[m]; n++){
	      c = ancSeq[i][m][sites[m][n]];
	      getNucStorage(j, k, m, sites[m][n], data[i]->Rtree, &c, &data[i],
			    1);
	      fprintf(out, "%c ",c);
	    }
	  }
	  fprintf(out, "\n");
	}
	else{
	  printf("fooey... %d %d\n",j, k);
	}
      }
    }
    
    if(out != stdout)  fclose(out);
    for(j=0; j<totLOCI; j++)
      free(sites[j]);
    free(sites);
    sites = NULL;
    free(Nsites);
    Nsites = NULL;
  }
}

/* ------------------------------------------------------------------------- */


void getSynNsSites(char *anc, char **cod, int Mc, int N, int *V, int **type)
{
  int i, j, k, *tmpType, *tmpType2;
  char c1[3], c2[3];
  tmpType = ivector(0,2);
  tmpType2 = ivector(0,2);
  tmpType[0] = tmpType[1] = tmpType[2] = 10; /* large value */
  tmpType2[0] = tmpType2[1] = tmpType2[2] = 10; /* large value */
  
  if(N==0)  return;
  else if(Mc == 1){ /* just compare ancestral to derived */
    synVSns(anc, cod[0], &tmpType);
    memcpy((*type), tmpType, 3*sizeof(int));
  }
  else if(Mc == 2){ /* 2 codons */
    if(normCod(anc, cod[0]) == normCod(anc, cod[1])){
      synVSns(anc, cod[0], &tmpType); /* codons equidistant from anc */
      synVSns(anc, cod[1], &tmpType2);
      for(i=0; i<3; i++){
	if(tmpType[i] == 0 || tmpType2[i] == 0)
	  (*type)[i] = tmpType[i]+tmpType2[i];
	else{ /* check to make sure only minimal path chosen */
	  if(cod[0][i] == cod[1][i])
	    (*type)[i] = IMIN(tmpType[i], tmpType2[i]);
	  else{
	    if(tmpType[i] == 1 && tmpType2[i] == 1)
	      (*type)[i] = 3;
	    if(tmpType[i] == 2 && tmpType2[i] == 2)
	      (*type)[i] = 4;
	    if((tmpType[i] == 1 && tmpType2[i] == 2) ||
	       (tmpType[i] == 2 && tmpType2[i] == 1))
	      (*type)[i] = 5;
	  }
	}
      }
    }
    else{
      /* temp codon c1: update all seg. sites fixed for derived allele  */
      memcpy(c1, anc, 3*sizeof(char));
      for(i=0; i<3; i++){
	if(cod[0][i] != anc[i] && cod[0][i] == cod[1][i]) /* both derived */
	  c1[i] = cod[0][i];
      }
      synVSns(anc, c1, type);
      /* all descending codons carry c1 + potentially other mutations */
      memcpy(c2, cod[0], 3*sizeof(char));
      synVSns(c1, c2, &tmpType2);
      for(i=0; i<3; i++){
	(*type)[i] += tmpType2[i];
	tmpType2[i] = 10;
      }
      memcpy(c2, cod[1], 3*sizeof(char));
      synVSns(c1, c2, &tmpType2);
      for(i=0; i<3; i++){
	(*type)[i] += tmpType2[i];
      }
    }
  }
  else if(Mc == 3){ /* 3 codons... */
    /* count up # ways to get synon vs non-synon change at each position */
    int Syn[3]={0, 0, 0}, Ns[3]={0, 0, 0}, AAnum=0;
    char der[3];
    memcpy(der, anc, 3*sizeof(char));
    /* first count observed number of amino acids */
    for(i=0; i<Mc; i++){
      int keep=1;
      if(AA(cod[i]) != AA(anc)){
	for(j=0; j<i; j++)
	  if(AA(cod[i]) == AA(cod[j])){
	    keep = 0;
	    break;
	  }
	if(keep)
	  AAnum++;
      }
    }
    for(i=0; i<3; i++){ /* get derived allele at each position */ 
      for(j=0; j<3; j++)
	if(cod[i][j] != anc[i]){
	  der[i] = cod[i][j];
	  break;
	}
    }
    for(i=0; i<3; i++){
      for(j=0; j<3; j++){
	if(j==i)  continue;
	k = 3 - i - j;
	memcpy(c1, anc, 3*sizeof(char));
	c1[i] = der[i];
	if(AA(c1) == AA(der))  Syn[i]++;
	else  Ns[i]++;
	c1[j] = der[j];
	if(AA(c1) == AA(der))  Syn[j]++;
	else  Ns[j]++;
	c1[k] = der[k];
	if(AA(c1) == AA(der))  Syn[k]++;
	else  Ns[k]++;
      }
    }
    (*type)[1] = 2; /* 2nd position always non-synon */
    if(Syn[0] == 0){
      (*type)[0] = 2; 
      if(AAnum == 2 && Syn[2] != 0)
	(*type)[2] = 1; /* likely had a synonymous mutation */
      else  (*type)[2] = 2; /* all 3 positions are nonsynonymous */
    }
    else if(Syn[2] == 0){
      (*type)[2] = 2;
      if(AAnum == 2 && Syn[0] != 0)
	(*type)[0] = 1;
      else  (*type)[0] = 2;
    }
    else{
      if(AAnum == 2){ /* might have been synonymous change */
	if(Syn[0] < Syn[2]){
	  (*type)[0] = 2;
	  (*type)[2] = 1;
	}
	else{
	  (*type)[0] = 1;
	  (*type)[2] = 2;
	}
      }
      else{ /* all mutations were nonsynonymous */
	(*type)[0] = (*type)[2] = 2;
      }
    }
  }
  else{
    fprintf(stderr,"convertSFS_CODE error: complex codon:  Mc=%d, N=%d:\n",
	    Mc, N);
    fprintf(stderr,"%c%c%c\n",anc[0],anc[1],anc[2]);
    for(i=0; i<Mc; i++){
      for(j=0; j<3; j++)
	fprintf(stderr,"%c",cod[i][j]);
      fprintf(stderr,"\n");
    }
  }
  free_ivector(tmpType, 0, 2);
  free_ivector(tmpType2, 0, 2);
}

/* ------------------------------------------------------------------------- */

int normCod(char *c1, char *c2)
{
  int i, x=0;
  for(i=0; i<3; i++)
    if(c1[i] != c2[i])  x++;
  return(x);
}

/* ------------------------------------------------------------------------- */

void synVSns(char *anc, char *der, int **type)
{
  int i, j, k, tmpTYPE[3]={10,10,10};
  char c[3], d[3];

  for(i=2; i>=0; i--){ /* mutate 3' end first, but loop over all positions */
    for(j=0; j<3; j++){
      if(j==i) continue;
      k = 3 - i -j; /* i + j + k = 3 */
      memcpy(c,anc,3*sizeof(char));
      memcpy(d,anc,3*sizeof(char));
      if(anc[i] == der[i])  tmpTYPE[i] = 0;
      else{
	c[i] = der[i];
	if(AA(c) == 0)  tmpTYPE[i] = 100; /* stop codon, disallow */
	else tmpTYPE[i] = 1+(AA(d)!=AA(c));
	d[i] = der[i];
      }
      if(anc[j] == der[j])  tmpTYPE[j] = 0;
      else{
	c[j] = der[j];
	if(AA(c) == 0)  tmpTYPE[j] = 100; /* stop codon, disallow */
	else  tmpTYPE[j] = 1+(AA(d)!=AA(c));
	d[j] = der[j];
      }
      if(anc[k] == der[k])  tmpTYPE[k] = 0;
      else{
	c[k] = der[k];
	if(AA(c) == 0)  tmpTYPE[k] = 100; /* stop codon, disallow */
	else  tmpTYPE[k] = 1+(AA(d)!=AA(c));
	d[k] = der[k];
      }
      if((tmpTYPE[0]+tmpTYPE[1]+tmpTYPE[2]) <
	 ((*type)[0]+(*type)[1]+(*type)[2])){
	(*type)[0] = tmpTYPE[0];
	(*type)[1] = tmpTYPE[1];
	(*type)[2] = tmpTYPE[2];
	if((*type)[1] == 1){
	  printf("ERROR:  %c%c%c > %c%c%c [%d%d%d] = %d%d%d! (synon. pos 2!)\n",
		 anc[0],anc[1],anc[2],der[0],der[1],der[2],i,j,k,(*type)[0],
		 (*type)[1],(*type)[2]);
	  exit(1);
	}
      }
    }
  }
}

/* ------------------------------------------------------------------------- */

int makeRecomb(char *anc, char *c1, char*c2, char *res)
{ /* does recombining c1 and c2 produce res? */
  int i, j;
  char t1[3], t2[3];

  for(i=0; i<3; i++){
    if(res[i] != c1[i] && res[i] != c2[i])
      return 0; /* not possible */
  }
  
  for(i=1; i<3; i++){
    for(j=0; j<3; j++){
      if(j<i){
	t1[j] = c1[j];
	t2[j] = c2[j];
      }
      else{
	t1[j] = c2[j];
	t2[j] = c1[j];
      }
    }
    if(normCod(t1, res) == 0)
      return 1;
    if(normCod(t2, res) == 0)
      return 1;
  }
  return 0;
}

/* ------------------------------------------------------------------------- */

void reduceCodons(char *anc, char **cods, int numCods, int *Mc, char ***MinCods)
{
  int i, j, k, m, keep;

  /* first copy all unique codons to MinCods */
  for(i=0; i<numCods; i++){
    if(normCod(anc, cods[i]) > 0){
      keep = 1;
      for(j=0; j<i; j++)
	if(normCod(cods[i], cods[j]) == 0){
	  keep = 0;
	  break;
	}
      if(keep){ /* codon does not match previous */
	(*Mc)++;
	if(*Mc == 1){
	  assert((*MinCods) = malloc(sizeof(**MinCods)));
	  assert((*MinCods)[0] = malloc(3*sizeof(***MinCods)));
	}
	else{
	  assert((*MinCods) = realloc(*MinCods, (*Mc)*sizeof(**MinCods)));
	  assert((*MinCods)[(*Mc)-1] = malloc(3*sizeof(***MinCods)));
	}
	memcpy((*MinCods)[(*Mc)-1], cods[i], 3*sizeof(char));
      }
    }
  }
  if((*Mc) == 1)  return;
  else{ /* try to eliminate codons from recombinant pairs */
    for(i=0; i<(*Mc); i++){
      for(j=i+1; j<(*Mc); j++){
	if(makeRecomb(anc, (*MinCods)[i], anc, (*MinCods)[j])){
	  /* anc X i = j, so eliminate j */
	  for(k=j; k<(*Mc)-1; k++){
	    for(m=0; m<3; m++)
	      (*MinCods)[k][m] = (*MinCods)[k+1][m];
	  }
	  j--;
	  (*Mc)--;
	  break;
	}
      }
    }
  }
}

/* ------------------------------------------------------------------------- */

void getPolySitesPop(struct mutStorage *storage, int loc, int *Nsites, 
		     long **sites, int pop, int *keep)
{
  int i;
  if(storage != NULL){
    if(storage->locus >= loc)
      getPolySitesPop(storage->Ltree, loc, Nsites, sites, pop, keep);
    if(storage->locus == loc){
      for(i=0; i<storage->numCarry; i++){
	if((pop == -1 || storage->pops[i] == pop) && storage->chrs[i] != -1 &&
	   (keep == NULL || keep[storage->chrs[i]])){
	  if((*Nsites) == 0){
	    assert((*sites) = malloc(sizeof(**sites)));
	    (*sites)[0] = storage->event->site;
	    (*Nsites) = 1;
	  }
	  else{
	    if(storage->event->site != (*sites)[(*Nsites)-1]){
	      (*Nsites)++;
	      assert((*sites) = realloc((*sites), (*Nsites)*sizeof(**sites)));
	      (*sites)[(*Nsites)-1] = storage->event->site;
	    }
	  }
	  break;
	}
      }
    }
    if(storage->locus <= loc)
      getPolySitesPop(storage->Rtree, loc, Nsites, sites, pop, keep);
  }
}

/* ------------------------------------------------------------------------- */

void getPolyFixedSites(struct mutStorage *storage, int loc, int *Nsites, 
		       long **sites)
{
  if(storage != NULL){
    if(storage->locus >= loc)
      getPolyFixedSites(storage->Ltree, loc, Nsites, sites);
    if(storage->locus == loc){
      if((*Nsites) == 0){
	assert((*sites) = malloc(sizeof(**sites)));
	(*sites)[0] = storage->event->site;
	(*Nsites) = 1;
      }
      else{
	if(storage->event->site != (*sites)[(*Nsites)-1]){
	  (*Nsites)++;
	  assert((*sites) = realloc((*sites), (*Nsites)*sizeof(**sites)));
	  (*sites)[(*Nsites)-1] = storage->event->site;
	}
      }
    }
    if(storage->locus <= loc)
      getPolyFixedSites(storage->Rtree, loc, Nsites, sites);
  }
}

/* ------------------------------------------------------------------------- */

void getObsFixedSites(struct mutStorage *storage, int loc, int *Nsites, 
		      long **sites, int pop1, int pop2, int ss1, int ss2,
		      struct mutStorage **rootStor, char *ancSeq)
{
  int i, j, numAl, keep;
  char c[4], tc;
  
  if(storage != NULL){
    if(storage->locus >= loc)
       getObsFixedSites(storage->Ltree, loc, Nsites, sites, pop1, pop2, ss1, 
			ss2, rootStor, ancSeq);
    if(storage->locus == loc){
      keep = 1;
      /* ensure site isn't already declared fixed */
      for(j=(*Nsites)-1; j>=0; j--)  
	if((*sites)[j] == storage->event->site){
	  keep = 0;
	  break;
	}
      if(keep){
	keep = 0;
	for(i=0; i<storage->numCarry; i++){
	  if(storage->pops[i] == pop1 || storage->pops[i] == pop2){
	    keep = 1;
	    break;
	  }
	}
	if(keep){
	  /* get all alleles from pop2, and compare to pop1 */
	  numAl = 0;
	  for(i=0; i<ss2; i++){
	    tc = ancSeq[storage->event->site];
	    getNucStorage(pop2, i, loc, storage->event->site,(*rootStor)->Rtree,
			  &tc, rootStor, 0);
	    if(numAl == 0){
	      numAl++;
	      c[0] = tc;
	    }
	    else if(numAl == 1 && tc != c[0]){
	      numAl++;
	      c[1] = tc;
	    }
	    else if(numAl == 2 && tc != c[0] && tc != c[1]){
	      numAl++;
	      c[2] = tc;
	    }
	    else if(numAl == 3 && tc != c[0] && tc != c[1] && tc != c[2]){
	      numAl++;
	      c[3] = tc;
	    }
	    if(tc != 'A' && tc != 'C' && tc != 'G' && tc != 'T'){
	      printf("og allele = %c!!!\n",tc);
	    }
	  }
	  for(i=0; i<ss1; i++){
	    if(!keep)  break;
	    tc = ancSeq[storage->event->site];
	    getNucStorage(pop1, i, loc, storage->event->site,(*rootStor)->Rtree,
			  &tc, rootStor, 0);
	    for(j=0; j<numAl; j++)
	      if(tc == c[j]){ /* alleles match, no fixed difference */
		keep = 0;
		break;
	      }
	  }
	  if(keep){  /* fixed difference! */
	    if((*Nsites) == 0){
	      assert((*sites) = malloc(sizeof(**sites)));
	      (*sites)[0] = storage->event->site;
	      (*Nsites) = 1;
	    }
	    else{
	      (*Nsites)++;
	      assert((*sites) = realloc((*sites), (*Nsites)*sizeof(**sites)));
	      (*sites)[(*Nsites)-1] = storage->event->site;
	    }
	  }
	}
      }
    }
    if(storage->locus <= loc)
      getObsFixedSites(storage->Rtree, loc, Nsites, sites, pop1, pop2, ss1,
		       ss2, rootStor, ancSeq);
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
    fprintf(stderr,"error in ToCGTA: %c\n",c);
    exit(-1);
  }
}

/* ------------------------------------------------------------------------- */

char FromCGTA(char c)
{
  switch(c){
  case('C'):
    return('0');
  case('c'):
    return('0');
  case('G'):
    return('1');
  case('g'):
    return('1');
  case('T'):
    return('2');
  case('t'):
    return('2');
  case('A'):
    return('3');
  case('a'):
    return('3');
  default:
    fprintf(stderr,"error:  %c is not a valid entry\n",c);
    abort();
    break;
  }
}

/* ------------------------------------------------------------------------- */

void memcpyEVENT(struct event *to, const struct event *from)
{
  if(to == NULL)
    exitNOW("NULL event in memcpyEVENT\n");
  to->site = from->site;
  to->gen = from->gen;
  to->genFix = from->genFix;
  to->ancNuc = from->ancNuc;
  to->derNuc = from->derNuc;
  to->nonsyn = from->nonsyn;
  to->ancAA = from->ancAA;
  to->derAA = from->derAA;
  to->fit = from->fit;
  to->CpG = from->CpG;
  to->fiveP = from->fiveP;
  to->threeP = from->threeP;
  if(to->nonsyn == 'i' || to->nonsyn == 'd' || to->nonsyn == 'v'){
    to->nSites = from->nSites;
    if(from->nucs != NULL){
      if(to->nucs == NULL)
	assert(to->nucs = malloc((to->nSites+1)*sizeof(char)));
      else
	assert(to->nucs = realloc(to->nucs, (to->nSites+1)*sizeof(char)));
      strcpy(to->nucs, from->nucs);
    }
  }
}

/* ------------------------------------------------------------------------- */

void exitNOW(const char *s)
{
  fprintf(stderr,"%s\n",s);
  abort();
}

/* ------------------------------------------------------------------------- */

void buildStorage(struct mutStorage **storage, struct mutStorage *tmp,
		  struct mutStorage **parent)
{
  int i;
  if((*storage) == NULL){ /* reached leaf node */
    assert((*storage) = malloc(sizeof(struct mutStorage)));
    assert((*storage)->event = malloc(sizeof(struct event)));
    (*storage)->event->nucs = NULL;
    (*storage)->numCarry = tmp->numCarry;
    (*storage)->pops = ivector(0,tmp->numCarry);
    assert((*storage)->chrs = malloc(tmp->numCarry*sizeof(*(*storage)->chrs)));
    (*storage)->locus = tmp->locus;
    for(i=0; i<tmp->numCarry; i++){
      (*storage)->pops[i] = tmp->pops[i];
      (*storage)->chrs[i] = tmp->chrs[i];
    }
    memcpyEVENT((*storage)->event,tmp->event);
    (*storage)->Rtree = NULL;
    (*storage)->Ltree = NULL;
    if(parent == NULL)  (*storage)->Parent = NULL;
    else  (*storage)->Parent = (*parent);
  }
  else{
    if((*storage)->locus < tmp->locus)
      buildStorage(&(*storage)->Rtree, tmp, storage);
    else if((*storage)->locus > tmp->locus)
      buildStorage(&(*storage)->Ltree, tmp, storage);
    else if((*storage)->event->site < tmp->event->site)
      buildStorage(&(*storage)->Rtree, tmp, storage);
    else if((*storage)->event->site > tmp->event->site)
      buildStorage(&(*storage)->Ltree, tmp, storage);
    else if((*storage)->event->gen < tmp->event->gen)
      buildStorage(&(*storage)->Rtree, tmp, storage);
    else if((*storage)->event->gen > tmp->event->gen)
      buildStorage(&(*storage)->Ltree, tmp, storage);
    else if((*storage)->event->genFix < tmp->event->genFix)
      buildStorage(&(*storage)->Rtree, tmp, storage);
    else if((*storage)->event->genFix > tmp->event->genFix)
      buildStorage(&(*storage)->Ltree, tmp, storage);
    else if((*storage)->event->fit < tmp->event->fit)
      buildStorage(&(*storage)->Rtree, tmp, storage);
    else /* new tmp into Ltree */
      buildStorage(&(*storage)->Ltree, tmp, storage);
  }
}

/* ------------------------------------------------------------------------- */

void printStorage(struct mutStorage *storage, int *x)
{
  int i;
  if(storage != NULL){
    printStorage(storage->Ltree, x);
    printf("%ld,%ld,%ld,%ld,%c%c%c,%c,%c,%d,%d,%1.4f,%ld",
	   storage->locus, storage->event->site, storage->event->gen,
	   storage->event->genFix, storage->event->fiveP, 
	   storage->event->ancNuc, storage->event->threeP, 
	   storage->event->derNuc, storage->event->nonsyn,
	   storage->event->ancAA, storage->event->derAA, 
	   storage->event->fit, storage->numCarry);
    for(i=0; i<storage->numCarry; i++)
      printf(",%d.%ld",storage->pops[i], storage->chrs[i]);
    printf(";");
    if((++(*x))%20 == 0)
      printf("\n");
    printStorage(storage->Rtree, x);
  }
}

/* ------------------------------------------------------------------------- */

void freeStorage(struct mutStorage *storage)
{
  if(storage != NULL){
    freeStorage(storage->Ltree);
    freeStorage(storage->Rtree);
    free(storage->pops);
    free(storage->chrs);
    free(storage->event);
    free(storage);
    storage = NULL;
  }
}

/* ------------------------------------------------------------------------- */

void updateSeq(char **seq,struct mutStorage *storage, int pop, int chr, int loc,
	       long *numInv, long ***invBK)
{
  long i, j, k;
  char *invSeq=NULL;
  if(storage != NULL){
    updateSeq(seq, storage->Rtree, pop, chr, loc, numInv, invBK);
    if(storage->locus == loc){
/*       printf("site=%ld; numCarry = %ld; type=%c\n",storage->event->site, storage->numCarry,storage->event->nonsyn); */
/*       fflush(stdout); */
      for(i=0; i<storage->numCarry; i++){
/* 	printf("carrier: %ld (%d)\n",storage->chrs[i],chr); */
	if(storage->pops[i] == pop && 
	   (storage->chrs[i] == -1 || storage->chrs[i] == chr)){
/* 	  printf("match chr\n"); */
/* 	  fflush(stdout); */
	  if(storage->event->nonsyn == '0' || storage->event->nonsyn == '1'){
/* 	    printf("mutation\n"); */
/* 	    fflush(stdout); */
	    k = 0;
	    for(j=0; j<(*numInv); j++){
	      if(storage->event->site >= (*invBK)[j][0] &&
		 storage->event->site < (*invBK)[j][1]){
		k = 1;
		(*seq)[(*invBK)[j][0]+(*invBK)[j][0]-
		       storage->event->site] =
		  storage->event->derNuc;
		break;
	      }
	    }
	    if(k == 0)
	      (*seq)[storage->event->site] = storage->event->derNuc;
	  }
	  else if(storage->event->nonsyn == 'i'){
/* 	    printf("insertion\n"); */
/* 	    fflush(stdout); */
	    j = strlen(*seq);
	    assert((*seq) = 
		   realloc(*seq, (j+storage->event->nSites+1)*
			   sizeof(char)));
	    j += storage->event->nSites;
	    (*seq)[j] = '\0';
	    for(; j>storage->event->site+storage->event->nSites; j--){
	      (*seq)[j] = (*seq)[j-storage->event->nSites];
	    }
	    for(k=0; k<storage->event->nSites; k++){
	      (*seq)[storage->event->site+k+1] = storage->event->nucs[k];
	    }
	  }
	  else if(storage->event->nonsyn == 'd'){
/* 	    printf("deletion\n"); */
/* 	    fflush(stdout); */
	    j = strlen(*seq);
	    for(j=storage->event->site+1;
		j<=strlen(*seq)-storage->event->nSites; j++){
	      (*seq)[j] = (*seq)[j+storage->event->nSites];
	    }
	    (*seq)[j] = '\0';
	    assert((*seq) = realloc(*seq, (j+1)*sizeof(char)));
	  }
	  else if(storage->event->nonsyn == 'v'){
/* 	    printf("inversion: nSites = %ld\n",storage->event->nSites); */
/* 	    fflush(stdout); */
	    (*numInv)++;
	    if((*numInv) == 1){
	      assert((*invBK) = malloc(sizeof(*invBK)));
	      assert((*invBK)[0] = malloc(2*sizeof(**invBK)));
	    }
	    else{
	      assert((*invBK) = realloc(*invBK, (*numInv)*sizeof(**invBK)));
	      assert((*invBK)[(*numInv)-1] = malloc(2*sizeof(***invBK)));
	    }
	    (*invBK)[(*numInv)-1][0] = storage->event->site+1;
	    (*invBK)[(*numInv)-1][1] = (storage->event->site + 
					storage->event->nSites+1);
	    assert(invSeq = malloc(storage->event->nSites*sizeof(char)));
	    for(j=0; j<storage->event->nSites; j++){
	      invSeq[j] = (*seq)[storage->event->site+1+j];
/* 	      printf("%ld:  %c\n",j,invSeq[j]); */
	    }
	    for(j=0; j<storage->event->nSites; j++){
	      (*seq)[storage->event->site+j+1] =
		invSeq[storage->event->nSites-j-1];
/* 	      printf("-> %ld:  %c\n",storage->event->site+storage->event->nSites-j,(*seq)[storage->event->site+storage->event->nSites-j]); */
	    }
	    free(invSeq);
	  }
	  else{
	    fprintf(stderr,"error, event->nonsyn = %c\n",
		    storage->event->nonsyn);
	    abort();
	  }
	}
      }
    }
    updateSeq(seq, storage->Ltree, pop, chr, loc, numInv, invBK);
  }
}

/* ------------------------------------------------------------------------- */

void rotateStorageNode(struct mutStorage *node, struct mutStorage **tree)
{
  struct mutStorage *parent, *grandparent;
  int x=0;
  
  if (node->Parent == NULL) return;
  
  parent = node->Parent;
  grandparent = parent->Parent;
  
  if (parent->Ltree == node){
    parent->Ltree = node->Rtree;
    if (parent->Ltree != NULL) 
      parent->Ltree->Parent = parent;
    node->Rtree = parent;
  } 
  else if (parent->Rtree == node){
    parent->Rtree = node->Ltree;
    if (parent->Rtree != NULL) 
      parent->Rtree->Parent = parent;
    node->Ltree = parent;
  } 
  else{
    fprintf(stderr,"rotateStorageNode error: parent's children not right\n");
    printStorage((*tree)->Rtree,&x);
    printf("\n");
    x=0;
    printStorage(node,&x);
    printf("\n");
    assert(0);
  }
  
  parent->Parent = node;
  node->Parent = grandparent;
  
  if (grandparent == NULL){
    (*tree)->Rtree = node;
    node->Parent = NULL;
  }
  else if (grandparent->Ltree == parent) 
    grandparent->Ltree = node;
  else if (grandparent->Rtree == parent) 
    grandparent->Rtree = node;
  else{
    fprintf(stderr,"rotateStorageNode error: grandparent's kids not right\n");
    x=0;
    printStorage((*tree)->Rtree,&x);
    printf("\n");
    x=0;
    printStorage(node,&x);
    printf("\n");
    assert(0);
  }
}

/* ------------------------------------------------------------------------- */

void splayStorage(struct mutStorage *node, struct mutStorage **tree)
{
  struct mutStorage *parent, *grandparent;
  
  if (node == NULL) return;
  while(1){
    if (node->Parent == NULL) break;
    
    parent = node->Parent;
    grandparent = parent->Parent;
    
    /* If the node's parent is the root of the tree, do one rotation */
    
    if (grandparent == NULL) 
      rotateStorageNode(node, tree);
    
    /* If we have a zig-zig, then rotate my parent, then rotate me */
    else if ((parent->Ltree  == node && grandparent->Ltree  == parent) ||
	     (parent->Rtree == node && grandparent->Rtree == parent)){
      rotateStorageNode(parent, tree);
      rotateStorageNode(node, tree);

      /* If we have a zig-zag, then rotate me twice */
    }
    else{
      rotateStorageNode(node, tree);
      rotateStorageNode(node, tree);
    }
  }
}

/* ------------------------------------------------------------------------- */

void getNucStorage(int pop, int ind, int locus, long site, 
		   struct mutStorage *storage, char *nuc, 
		   struct mutStorage **root, int toSplay)
{
  int i;
  if(storage != NULL){
    if(storage->locus > locus || 
       (storage->locus == locus && storage->event->site >= site))
      getNucStorage(pop, ind, locus, site, storage->Ltree, nuc, root, toSplay);
    if(storage->locus == locus && storage->event->site == site){
      for(i=0; i<storage->numCarry; i++){
	if(storage->pops[i] == pop && 
	   (storage->chrs[i] == -1 || storage->chrs[i]  == ind))
	  (*nuc) = storage->event->derNuc;
      }
      if(toSplay)
	splayStorage(storage, root); /* move node to root */
    }
    if(storage->locus < locus ||
       (storage->locus == locus && storage->event->site <= site)){
      getNucStorage(pop, ind, locus, site, storage->Rtree, nuc, root, toSplay);
    }
  }
}

/* ------------------------------------------------------------------------- */

void updateObsSFS(int *sfs, struct mutStorage **storage, int pop, int ss,
		  int locus, char TYPE, int og, char *ancSeq, int ogSS, int PSEX,
		  int *MALES, int MONLY, int ploid, int **keep, int max)
{
  int i, j, k, site, id[3], Mc, *type, N, *V, keepIT, f, ind;
  int numSegSites=0, numAl;
  long *SegSites=NULL;
  char **cods=NULL, *anc=NULL, **MinCods=NULL, c1='x';

  V=ivector(0,2);
  type = ivector(0,2);
  /* get list of polymorphic sites for each locus */
  getPolySitesPop((*storage)->Rtree, locus, &numSegSites, &SegSites, pop,
		  keep[pop]);
  
  assert(cods = malloc(ss*sizeof(*cods)));
  for(i=0; i<ss; i++)
    assert(cods[i] = malloc(3*sizeof(**cods)));
  assert(anc = malloc(3*sizeof(*anc)));
  for(i=0; i<numSegSites; i++){
    V[0] = V[1] = V[2] = 0;
    keepIT = 1;
    /* get observed codons */
    site = SegSites[i];
    for(j=0; j<3; j++){
      anc[j] = ancSeq[site-(site%3)+j];
      getNucStorage(og, 0, locus, (site+(j-(site%3))), (*storage)->Rtree,
		    &anc[j], storage, 1);
      for(k=1; k<ogSS; k++){
	c1=anc[j];
	getNucStorage(og, k, locus, (site+(j-(site%3))), (*storage)->Rtree, &c1,
		      storage, 1);
	if(c1 != anc[j]){
	  V[j] = 1; /* mark sites with polymorphic outgroup */
	}
      }
    }
    id[0] = id[1] = id[2] = 0;
    ind = j = 0;
    while(ind < ss && ind < max){
      if(!keep[j]){
	j++;
	continue;
      }
      if(PSEX && ((!MONLY && (j < MALES[pop] || j%ploid < ploid/2)) ||
		  (MONLY && j >= MALES[pop] && j%ploid >= ploid/2))){
	j++;
	continue;
      }
      for(k=0; k<3; k++){
	cods[ind][k] = ancSeq[site+(k-(site%3))];
	getNucStorage(pop, j, locus, (site+(k-(site%3))), (*storage)->Rtree,
		      &cods[ind][k], storage, 1);
	if(ind>0 && cods[ind][k] == cods[ind-1][k])
	  id[k]++;
      }
      j++;
      ind++;
    }
    for(k=0; k<3; k++){
      if(id[k] == ss-1){ /* nucleotide fixed in sample */
	anc[k] = cods[0][k];
	V[k] = 0; /* ignore potential outgroup polymorphism at fixed sites */
      }
      if(V[k] == 1)  keepIT = 0; /* polymorphic in both populations: skip codon */
    }
    if(keepIT){
      Mc=0;
      reduceCodons(anc, cods, ss, &Mc, &MinCods);
      N=0;
      for(j=0; j<3; j++){
	type[j] = 10; /* default large value */
	numAl=0;
	for(k=0; k<Mc; k++){
	  if(anc[j] != MinCods[k][j]){
	    if(numAl == 0){
	      c1 = MinCods[k][j];
	      numAl++;
	      V[N++] = j;
	    }
	    else if(MinCods[k][j] != c1 && numAl == 1){
	      keepIT = 0; /* discard codon if site has 3 alleles */
	      break;
	    }
	  }
	}
      }
      if(keepIT){ /* ensures at most 2 nucleotides at each site */
	getSynNsSites(anc, MinCods, Mc, N, V, &type);
	for(j=0; j<N; j++){
	  if(TYPE == '2' || TYPE-'0' == (type[V[j]]-1)){
	    f=0;
	    ind = 0;
	    for(k=0; k<max; k++){
	      if(keep[k]){
		if(cods[ind][V[j]] != anc[V[j]])
		  f++;
		ind++;
	      }
	    }
	    sfs[f]++;
	  }
	}
      }
    }
    /* skip sites in same codon... */
    if(site%3 < 2){
      if(i<numSegSites-1){ 
	if(SegSites[i+1] <= site+(2-site%3))
	  i++;
	if(i<numSegSites-1){ /* can have 3 segsites in single codon */
	  if(SegSites[i+1] <= site+(2-site%3))
	    i++;
	}
      }
    }
    for(j=0; j<Mc; j++)
      if(MinCods != NULL){
	free(MinCods[j]);
	MinCods[j] = NULL;
      }
    if(MinCods != NULL){
      free(MinCods);
      MinCods = NULL;
    }
  }
  free(SegSites);
  free(anc);
  for(i=0; i<ss; i++)
    free(cods[i]);
  free(cods);
  free_ivector(V, 0, 2);
  free_ivector(type,0,2);
}

/* ------------------------------------------------------------------------- */

void getObsMK(int *MK, int ING, int OUTG, struct mutStorage **storage, int loc,
	      int ss1, int ss2, char *ancSeq)
{
  int PS=0, PN=0, FS=0, FN=0;
  int i, j, k, m, id[3], Mc, McOG, *type, *tmpTYPE, N, *V, keep;
  int numSegSites=0, numFixedSites=0, numAl;
  long *SegSites=NULL, *FixedSites=NULL, site; 
  char **cods=NULL, *anc=NULL, **ancCods=NULL, **MinCods=NULL, c1, c2, c3;
  
  V=ivector(0,2);
  type = ivector(0,2);
  tmpTYPE = ivector(0,2);

  /* first work on polymorphisms */
  /* get list of polymorphic sites for each locus */
  getPolySitesPop((*storage)->Rtree, loc, &numSegSites, &SegSites, ING, NULL);
  assert(cods = malloc(ss1*sizeof(*cods)));
  for(i=0; i<ss1; i++)
    assert(cods[i] = malloc(3*sizeof(**cods)));
  assert(anc = malloc(3*sizeof(*anc)));

  for(i=0; i<numSegSites; i++){
    V[0] = V[1] = V[2] = 0;
    /* get observed codons */
    site = SegSites[i];
    id[0] = id[1] = id[2] = 0;
    for(j=0; j<ss1; j++){
      for(k=0; k<3; k++){
	cods[j][k] = ancSeq[site+(k-(site%3))];
	getNucStorage(ING, j, loc, (site+(k-(site%3))), (*storage)->Rtree,
		      &cods[j][k], storage, 1);
	if(j>0 && cods[j][k] == cods[j-1][k])
	  id[k]++;
      }
    }
    memcpy(anc, cods[0], 3*sizeof(char));
    Mc=0;
    N=0;
    reduceCodons(anc, cods, ss1, &Mc, &MinCods);
    keep = 1;
    for(j=0; j<3; j++){
      type[j] = 10; /* default large value */
      numAl=0;
      for(k=0; k<Mc; k++){
	if(anc[j] != MinCods[k][j]){
	  if(numAl == 0){
	    c1 = MinCods[k][j];
	    numAl++;
	    V[N++] = j;
	  }
/* 	  else if(MinCods[k][j] != c1 && numAl == 1){ */
/* 	    keep = 0; /\* discard codon if site has 3 alleles *\/ */
/* 	    break; */
/* 	  } */
	}
      }
    }
    if(keep){ /* ensures at most 2 nucleotides at each site */
      getSynNsSites(anc, MinCods, Mc, N, V, &type);
      for(j=0; j<N; j++){
	if(type[V[j]] == 1)
	  PS++;
	else if(type[V[j]] == 2)
	  PN++;
	else if(type[V[j]] == 3)
	  PS += 2;
	else if(type[V[j]] == 4)
	  PN += 2;
	else if(type[V[j]] == 5){
	  PS++;
	  PN++;
	}
      }
    }
    /* skip sites in same codon... */
    if(site%3 < 2){
      if(i<numSegSites-1){ 
	if(SegSites[i+1] <= site+(2-site%3))
	  i++;
	if(i<numSegSites-1){ /* can have 3 segsites in single codon */
	  if(SegSites[i+1] <= site+(2-site%3))
	    i++;
	}
      }
    }
    for(j=0; j<Mc; j++)  free(MinCods[j]);
    free(MinCods);
    MinCods = NULL;
  }

  /* Now work on fixed differences */
  assert(ancCods = malloc(6*sizeof(*ancCods))); /* at most 6 ancestral codons */
  for(i=0; i<6; i++)  assert(ancCods[i] = malloc(3*sizeof(**ancCods)));
  getObsFixedSites((*storage)->Rtree, loc, &numFixedSites, &FixedSites, ING,
		   OUTG, ss1, ss2, storage, ancSeq);
/*   printf("observed fixed sites = %d\n",numFixedSites); */
  for(i=0; i<numFixedSites; i++){
    site = FixedSites[i];
    
    /* get observed codons in ingroup */
    Mc = 0;
    for(j=0; j<ss1; j++){
      c1 = ancSeq[site-(site%3)];
      c2 = ancSeq[site+1-(site%3)];
      c3 = ancSeq[site+2-(site%3)];
      getNucStorage(ING, j,loc,site-(site%3),(*storage)->Rtree, &c1,storage,1);
      getNucStorage(ING, j,loc,site+1-(site%3),(*storage)->Rtree,&c2,storage,1);
      getNucStorage(ING, j,loc,site+2-(site%3),(*storage)->Rtree,&c3,storage,1);
      if(Mc == 0){
	Mc++;
	cods[0][0] = c1;
	cods[0][1] = c2;
	cods[0][2] = c3;
      }
      else{
	keep = 1;
	for(k=0; k<Mc; k++){
	  if(cods[k][0] == c1 && cods[k][1] == c2 && cods[k][2] == c3){
	    keep = 0;
	    break;
	  }
	}
	if(keep){
	  cods[Mc][0] = c1;
	  cods[Mc][1] = c2;
	  cods[Mc][2] = c3;
	  Mc++;
	}
      }
    }
    
    /* get observed codons in outgroup */
    McOG=0;
    for(j=0; j<ss2; j++){
      c1 = ancSeq[site-(site%3)];
      c2 = ancSeq[site+1-(site%3)];
      c3 = ancSeq[site+2-(site%3)];
      getNucStorage(OUTG,j,loc,site-(site%3),(*storage)->Rtree,&c1,storage,1);
      getNucStorage(OUTG,j,loc,site+1-(site%3),(*storage)->Rtree,&c2,storage,1);
      getNucStorage(OUTG,j,loc,site+2-(site%3),(*storage)->Rtree,&c3,storage,1);
      if(McOG == 0){
	McOG++;
	ancCods[0][0] = c1;
	ancCods[0][1] = c2;
	ancCods[0][2] = c3;
      }
      else{
	keep = 1;
	for(k=0; k<McOG; k++){
	  if(ancCods[k][0] == c1 && ancCods[k][1] == c2 && ancCods[k][2] == c3){
	    keep = 0;
	    break;
	  }
	}
	if(keep){
	  if(McOG >= 6){
	    assert(ancCods = realloc(ancCods, (McOG+1)*sizeof(*ancCods)));
	    assert(ancCods[McOG] = malloc(3*sizeof(**ancCods)));
	  }
	  ancCods[McOG][0] = c1;
	  ancCods[McOG][1] = c2;
	  ancCods[McOG][2] = c3;
	  McOG++;
	}
      }
    }

    /* eliminate polymorphisms with matching alleles */
    for(m=0; m<3; m++){
      keep = 1;
      for(j=0; j<Mc; j++){
	if(!keep)  break;
	for(k=0; k<McOG; k++){
	  if(cods[j][m] == ancCods[k][m]){
	    keep = 0;
	    c1 = cods[j][m];
	    break;
	  }
	}
      }
      if(!keep){
	for(j=0; j<Mc; j++)  cods[j][m] = c1;
	for(j=0; j<McOG; j++)  ancCods[j][m] = c1;
      }
    }
    for(j=1; j<Mc; j++){
      if(cods[j][0] == cods[j-1][0] && cods[j][1] == cods[j-1][1] && 
	 cods[j][2] == cods[j-1][2]){ /* exact match, eliminate */
	for(k=j+1; k<Mc; k++){
	  cods[k-1][0] = cods[k][0];
	  cods[k-1][1] = cods[k][1];
	  cods[k-1][2] = cods[k][2];
	}
	Mc--;
      }
    }
    for(j=1; j<McOG; j++){
      if(ancCods[j][0] == ancCods[j-1][0] && 
	 ancCods[j][1] == ancCods[j-1][1] && 
	 ancCods[j][2] == ancCods[j-1][2]){ /* exact match, eliminate */
	for(k=j+1; k<McOG; k++){
	  ancCods[k-1][0] = ancCods[k][0];
	  ancCods[k-1][1] = ancCods[k][1];
	  ancCods[k-1][2] = ancCods[k][2];
	}
	McOG--;
      }
    }
    
    if(McOG == 1 && Mc == 1){ /* simplest case */
      type[0] = type[1] = type[2] = 10;
      tmpTYPE[0] = tmpTYPE[1] = tmpTYPE[2] = 10;
      synVSns(ancCods[0], cods[0], &type);
      synVSns(cods[0], ancCods[0], &type);
      for(j=0; j<3; j++){
	type[j] = (tmpTYPE[j] < type[j] ? tmpTYPE[j] : type[j]);
	if(type[j] == 1)  FS++;
	else if(type[j] == 2)  FN++;
      }
    }
    else{ 
      if(McOG == 1){
	type[0] = type[1] = type[2] = 10;
	tmpTYPE[0] = tmpTYPE[1] = tmpTYPE[2] = 10;
	
	synVSns(ancCods[0], cods[0], &type);
	for(j=1; j<Mc; j++){
	  synVSns(ancCods[0], cods[j], &tmpTYPE);
	  for(k=0; k<3; k++){
	    if(tmpTYPE[k] != 0 && (type[k] == 0 || type[k] > tmpTYPE[k]))
	      type[k] = tmpTYPE[k];
	  }
	}
      }
      else if(Mc == 1){
	type[0] = type[1] = type[2] = 10;
	tmpTYPE[0] = tmpTYPE[1] = tmpTYPE[2] = 10;
	
	synVSns(ancCods[0], cods[0], &type);
	for(j=1; j<McOG; j++){
	  synVSns(ancCods[j], cods[0], &tmpTYPE);
	  for(k=0; k<3; k++){
	    if(tmpTYPE[k] != 0 && (type[k] == 0 || type[k] > tmpTYPE[k]))
	      type[k] = tmpTYPE[k];
	  }
	}
      }
      else{
	fprintf(stderr,"complex situation:  codons in outgroup = {");
	for(j=0; j<McOG; j++){
	  fprintf(stderr,"%c%c%c, ",ancCods[j][0],ancCods[j][1],ancCods[j][2]);
	}
	fprintf(stderr,"}\ncodons in ingroup = {");
	for(j=0; j<Mc; j++){
	  fprintf(stderr,"%c%c%c, ",cods[j][0],cods[j][1],cods[j][2]);
	}
	fprintf(stderr,"}\n");
      }
      for(j=0; j<3; j++){
	if(type[j] == 1)  FS++;
	else if(type[j] == 2)  FN++;
      }
    }
    
    /* skip sites in same codon... */
    if(site%3 < 2){
      if(i<numFixedSites-1){ 
	if(FixedSites[i+1] <= site+(2-site%3))
	  i++;
	if(i<numFixedSites-1){ /* can have 3 segsites in single codon */
	  if(FixedSites[i+1] <= site+(2-site%3))
	    i++;
	}
      }
    }
/*     printf("%d (%ld):  FS=%d  FN=%d  FN+FS=%d;  %d%d%d\n", */
/* 	   i,FixedSites[i],FS,FN,FN+FS,type[0],type[1],type[2]); */
/*     if(type[0] == 1) printf("%ld\n",FixedSites[i] - FixedSites[i]%3); */
/*     if(type[1] == 1) printf("%ld\n",FixedSites[i] - FixedSites[i]%3+1); */
/*     if(type[2] == 1) printf("%ld\n",FixedSites[i] - FixedSites[i]%3+2); */
/*     if(FixedSites[i] == -1){ */
/*       printf("ING:"); */
/*       for(i=0; i<Mc; i++){ */
/* 	printf("\t%c%c%c\n",cods[i][0],cods[i][1],cods[i][2]); */
/*       } */
/*       printf("OUTG:"); */
/*       for(i=0; i<McOG; i++){ */
/* 	printf("\t%c%c%c\n",ancCods[i][0],ancCods[i][1],ancCods[i][2]); */
/*       } */
/*       printf("TYPE:\t%d%d%d\n",type[0],type[1],type[2]); */
/*       printf("\n\n"); */
/*       printf("FULL:"); */
/*       for(j=0; j<ss1; j++){ */
/* 	c1 = ancSeq[site-(site%3)]; */
/* 	c2 = ancSeq[site+1-(site%3)]; */
/* 	c3 = ancSeq[site+2-(site%3)]; */
/* 	getNucStorage(ING, j,loc,site-(site%3),(*storage)->Rtree, */
/* 		      &c1,storage,1); */
/* 	getNucStorage(ING, j,loc,site+1-(site%3),(*storage)->Rtree, */
/* 		      &c2,storage,1); */
/* 	getNucStorage(ING, j,loc,site+2-(site%3),(*storage)->Rtree, */
/* 		      &c3,storage,1); */
/* 	printf("\t%c%c%c\n",c1,c2,c3); */
/*       } */
/*       printf("OG:"); */
/*       for(j=0; j<ss2; j++){ */
/* 	c1 = ancSeq[site-(site%3)]; */
/* 	c2 = ancSeq[site+1-(site%3)]; */
/* 	c3 = ancSeq[site+2-(site%3)]; */
/* 	getNucStorage(OUTG, j,loc,site-(site%3),(*storage)->Rtree, */
/* 		      &c1,storage,1); */
/* 	getNucStorage(OUTG, j,loc,site+1-(site%3),(*storage)->Rtree, */
/* 		      &c2,storage,1); */
/* 	getNucStorage(OUTG, j,loc,site+2-(site%3),(*storage)->Rtree, */
/* 		      &c3,storage,1); */
/* 	printf("\t%c%c%c\n",c1,c2,c3); */
/*       } */
/*       exit(1); */
/*     } */
  }
  
  /*   exit(1); */
  
  /* MK = <PS> <PN> <FS> <FN> */
  MK[0] += PS;
  MK[1] += PN;
  MK[2] += FS;
  MK[3] += FN;
  free(SegSites);
  free(FixedSites);
  free(anc);
  for(i=0; i<ss1; i++)
    free(cods[i]);
  free(cods);
  free_ivector(V, 0, 2);
  free_ivector(type,0,2);
  free_ivector(tmpTYPE,0,2);
  for(i=0; i<6; i++)  free(ancCods[i]);
  free(ancCods);
}

/* ------------------------------------------------------------------------- */

void updateTrueSFS(int *sfs, struct mutStorage *storage, int pop, int ss,
		   int loc,char TYPE, int PSEX, int *MALES, int MONLY, int ploid,
		   int **keep)
{
  int i, f;
  if(storage != NULL){
    if(storage->locus == loc){
      if(TYPE == '2' || TYPE == storage->event->nonsyn){
	f = 0;
	for(i=0; i<storage->numCarry; i++){
	  if(storage->pops[i] == pop && storage->chrs[i] != -1 &&
	     keep[pop][storage->chrs[i]] &&
	     (!PSEX ||
	      (PSEX && ((!MONLY && (storage->chrs[i] < MALES[pop] || 
				    storage->chrs[i]%ploid < ploid/2))
			||(MONLY && storage->chrs[i] >= MALES[pop] && 
			   storage->chrs[i]%ploid == 0)))))
	    f++;
	}
	if(f > 0 && ((!PSEX && f < ss) || 
		     (PSEX && ((!MONLY && f < 3.0/4.0*ss)||
			       (MONLY && f < 1.0/4.0*ss)))))
	  sfs[f]++;
      }
    }
    if(storage->locus >= loc)
      updateTrueSFS(sfs, storage->Ltree, pop, ss, loc, TYPE, PSEX, MALES, MONLY,
		    ploid, keep);
    if(storage->locus <= loc)
      updateTrueSFS(sfs, storage->Rtree, pop, ss, loc, TYPE, PSEX, MALES, MONLY,
		    ploid, keep);
  }
}

/* ------------------------------------------------------------------------- */

void updateDerivedAlleles(int ***cnt, struct mutStorage *storage, int *ss,
			  char TYPE, int ploid, int HH)
{
  int i, j;
  if(storage != NULL){
    if(TYPE == '2' || TYPE == storage->event->nonsyn){
      if(HH == 2){
	for(i=0; i<storage->numCarry; i++){
	  if(storage->chrs[i] == -1){
	    for(j=0; j<ss[storage->pops[i]]; j++)
	      (*cnt)[storage->pops[i]][(int)(floor(j/(ploid+0.0)))]++;
	  }
	  else{
	    (*cnt)[storage->pops[i]][(int)(floor(storage->chrs[i]/
						 (ploid+.0)))]++;
	  }
	}
      }
      else{
	for(i=0; i<storage->numCarry; i++){
	  if(storage->chrs[i] == -1){ /* all indiv homozygous */
	    if(HH == 1)
	      for(j=0; j<ss[storage->pops[i]]; j++)
		(*cnt)[storage->pops[i]][j]++;
	  }
	  else if(i < storage->numCarry-1 &&
		  storage->pops[i] == storage->pops[i+1] &&
		  (int)(storage->chrs[i]/2.0) == (int)(storage->chrs[i+1]/2.0)){
	    /* found homozygous individual */
	    if(HH == 1)
	      (*cnt)[storage->pops[i]][(int)(storage->chrs[i]/2.0)]++;
	    i++;
	  }
	  else if(i == storage->numCarry-1 ||
		  (i < storage->numCarry-1 && 
		   ((int)(storage->chrs[i]/2.0) != (int)(storage->chrs[i+1]/2.0)
		    || storage->pops[i] != storage->pops[i+1]))){
	    /* heterozygous individual */
	    if(HH == 0)
	      (*cnt)[storage->pops[i]][(int)(storage->chrs[i]/2.0)]++;
	  }
	}
      }
    }
    updateDerivedAlleles(cnt, storage->Ltree, ss, TYPE, ploid, HH);
    updateDerivedAlleles(cnt, storage->Rtree, ss, TYPE, ploid, HH);
  }
}

/* ------------------------------------------------------------------------- */

void printFitness(struct mutStorage *storage, int pop, int nShare, int *share,
		  int private, int fixed, FILE *out)
{
  int i, j, k, *p;

  if(storage == NULL)
    return;
  if(fabs(storage->event->fit) > 1e-9){ /* worth looking at... */
    if(pop == -1){ /* shared */
      k = 0; /* populations found so far */
      assert(p = malloc(nShare*sizeof(*p)));
      for(i=0; i<nShare; i++)
	p[i] = 0;
      for(i=0; i<storage->numCarry; i++){
	if(k == nShare)
	  break;
	for(j=0; j<nShare; j++){
	  if(storage->pops[i] == share[j] && p[j] == 0){
	    p[j] = 1;
	    k++;
	    break;
	  }
	}
      }
      if(k == nShare){
	fprintf(out,"%1.9f\n",storage->event->fit);
      }
      free(p);
    }
    else{
      k = 0;
      for(i=0; i<storage->numCarry; i++){
	if(storage->pops[i] == pop)
	  k = 1; /* keep it! */
	if(private == 1 && storage->pops[i] != pop){
	  /* reject because its not a private mutation */
	  k = 0;
	  break;
	}
	if(fixed == 0 && storage->pops[i] == pop && storage->chrs[i] == -1){
	  /* reject because it is fixed in population */
	  k = 0;
	  break;
	}
      }
      if(k == 1)
	fprintf(out,"%1.9f\n",storage->event->fit);
    }
  }
  printFitness(storage->Ltree, pop, nShare, share, private, fixed, out);
  printFitness(storage->Rtree, pop, nShare, share, private, fixed, out);
}

/* ------------------------------------------------------------------------- */

void updatePrivShared(struct mutStorage *storage, int **ps, char TYPE,int *pops)
{
  int i, s0=0, s1=0;

  if(storage == NULL)
    return;
  
  if(TYPE == '2' || TYPE == storage->event->nonsyn){
    for(i=0; i<storage->numCarry; i++){
      if(storage->pops[i] == pops[0] && s0 == 0){
/* 	if(storage->chrs[i] == -1) */
/* 	  s0 = -1; */
/* 	else */
	  s0 = 1;
	if(s0 == 1 && s1 == 1){
	  break;
	}
      }
      else if(storage->pops[i] == pops[1] && s1 == 0){
/* 	if(storage->chrs[i] == -1) */
/* 	  s1 = -1; */
/* 	else */
	  s1 = 1;
	if(s0==1 && s1==1){
	  break;
	}
      }
    }
    if(s0 == 1 && s1 == 1)
      (*ps)[2]++;
    else if(s0 == 1 && s1 == 0)
      (*ps)[0]++;
    else if(s1 == 1 && s0 == 0)
      (*ps)[1]++;
  }
  updatePrivShared(storage->Ltree, ps, TYPE, pops);
  updatePrivShared(storage->Rtree, ps, TYPE, pops);
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
  default: 
    fprintf(stderr, "convertSFS_CODE error: %d is not a valid amino acid...\n",
	    c);
    abort();
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
  default:
    return 0;
  }
}

/* ------------------------------------------------------------------------- */

int AA(const char *c)
{
  char s[60];
  
  if(c[0]=='C'){
    if(c[1]=='C')
      return 1;      /* PROLINE */
    else if(c[1]=='G')
      return 2;      /* ARGININE */
    else if(c[1]=='T')
      return 3;      /* LEUCINE */
    else if(c[1]=='A'){
      if(c[2]=='C' || c[2]=='T')
	return 4;  /* HISTIDINE */
      else if(c[2]=='G' || c[2]=='A')
	return 5;  /* GLUTAMINE */
      else{
	fprintf(stderr,"%c%c%c not an AA\n",c[0],c[1],c[2]);
	exitNOW("\n");
	return(-1);
      }
    }
    else{
      fprintf(stderr,"%c%c%c not an AA\n",c[0],c[1],c[2]);
      exitNOW("\n");
      return(-1);
    }
  }
  else if(c[0]=='G'){
    if(c[1]=='C')
      return 6;      /* ALANINE */
    else if(c[1]=='G')
      return 7;      /* GLYCINE */
    else if(c[1]=='T')
      return 8;      /* VALINE */
    else if(c[1]=='A'){
      if(c[2]=='C' || c[2]=='T')
	return 9;  /* ASPARTIC ACID */
      else
	return 10; /* GLUTAMIC ACID */
    }
    else{
      fprintf(stderr,"%c%c%c not an AA\n",c[0],c[1],c[2]);
      exitNOW("\n");
      return(-1);
    }
  }
  else if(c[0]=='T'){
    if(c[1]=='C')
      return 11;     /* SERINE */
    else if(c[1]=='G'){
      if(c[2]=='C' || c[2]=='T')
	return 12; /* CYSTEINE */
      else if(c[2]=='G')
	return 13; /* TRYPTOPHAN */
      else if(c[2]=='A')
	return 0;  /* STOP */
      else{
	fprintf(stderr,"%c%c%c not an AA\n",c[0],c[1],c[2]);
	exitNOW("\n");
	return(-1);
      }
    }
    else if(c[1]=='T'){
      if(c[2]=='C' || c[2]=='T')
	return 14;  /* PHENYLALANINE */
      else
	return 3;   /* LEUCINE */
    }
    else if(c[2]=='C' || c[2]=='T')
      return 15;      /* TYROSINE */
    else
      return 0;       /* STOP */
  }
  else if(c[1]=='C')
    return 16;          /* THREONINE */
  else if(c[1]=='G'){
    if(c[2]=='C' || c[2]=='T')
      return 11;      /* SERINE */
    else
      return 2;       /* ARGININE */
  }
  else if(c[1]=='T'){
    if(c[2]=='G')
      return 17;      /* METHIONINE (START) */
    else
      return 18;      /* ISOLEUCINE */
  }
  else if(c[1]=='A'){
    if(c[2]=='C' || c[2]=='T')
      return 19;      /* ASPARAGINE */
    else
      return 20;      /* LYSINE */
  }
  else{
    sprintf(s,"error in AA: %c%c%c\n",c[0],c[1],c[2]);
    exitNOW(s);
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

void isFixedSame(struct mutStorage *storage, int pop, int locus, 
		 struct event *ev, int ss, int *fixed)
{
  int i, f;
  if((*fixed) == 1)  return;
  if(storage != NULL){
    if(storage->locus == locus && storage->event->site == ev->site &&
       storage->event->gen == ev->gen){ /* check for same exact mutation */
      f = 0;
      for(i=0; i<storage->numCarry; i++){
	if(storage->pops[i] == pop){
	  if(storage->chrs[i] == -1){
	    (*fixed) = 1;
	    return;
	  }
	  else if(storage->chrs[i] < ss)
	    f++;
	}
      }
      if(f == ss){
	(*fixed) = 1;
	return;
      }
    }
    if(storage->locus >= locus)
      isFixedSame(storage->Ltree, pop, locus, ev, ss, fixed);
    if(storage->locus <= locus)
      isFixedSame(storage->Rtree, pop, locus, ev, ss, fixed);
  }
}

/* ------------------------------------------------------------------------- */

void getTrueMK(int *MK, int ING, int OUTG, struct mutStorage *storage, int loc,
	       struct mutStorage *rootStor, int ss1, int ss2)
{ /* MK = <PS> <PN> <FS> <FN> */
  int i, j, f1, f2, PS=0, PN=0, FS=0, FN=0;
  if(storage != NULL){
    if(storage->locus >= loc)
      getTrueMK(MK, ING, OUTG, storage->Ltree, loc, rootStor, ss1, ss2);
    
    if(storage->locus == loc){
      f1 = f2 = 0;
      for(i=0; i<storage->numCarry; i++){
	if(storage->pops[i] == ING && storage->chrs[i] == -1){ 
	  /* fixed difference? */
	  j = 0;
	  isFixedSame(rootStor, OUTG, loc, storage->event, ss2, &j);
	  if(j) /* mutation fixed in both, ignore */
	    break;
	  else{
	    if(storage->event->nonsyn == '0')  FS = 1;
	    else  FN = 1;
	  }
	}
	else if(storage->pops[i] == OUTG && storage->chrs[i] == -1){ 
	  /* fixed difference? */
	  j = 0;
	  isFixedSame(rootStor, ING, loc, storage->event, ss1, &j);
	  if(j) /* mutation fixed in both, ignore */
	    break;
	  else{
	    if(storage->event->nonsyn == '0')  FS = 1;
	    else  FN = 1;
	  }
	}
	else if(storage->pops[i] == ING && storage->chrs[i] < ss1)
	  f1++;
	else if(storage->pops[i] == OUTG && storage->chrs[i] < ss2)
	  f2++;
      }
      if((f1 == ss1 && f2 == 0) || (f1 == 0 && f2 == ss2)){ /* fixed! */
	if(storage->event->nonsyn == '0')  FS = 1;
	else  FN = 1;
      }
      else if((f2 == 0 || f2 == ss2) && f1 != 0){ /* polymorphic in ING */
	if(storage->event->nonsyn == '0')  PS = 1;
	else  PN = 1;
      }
      /* MK = <PS> <PN> <FS> <FN> */
      MK[0] += PS;
      MK[1] += PN;
      MK[2] += FS;
      MK[3] += FN;
    }
    if(storage->locus <= loc)
      getTrueMK(MK, ING, OUTG, storage->Rtree, loc, rootStor, ss1, ss2);
  }
}

/* ------------------------------------------------------------------------- */

double LnNchooseK(int n, int k)
{
  double v=0;

  v = lgamma(n+1.0) - lgamma(k+1.0) - lgamma(n-k+1.0);
  return v;
}

/* ------------------------------------------------------------------------- */

double HypergeometricPMF(int a, int r1, int r2, int c1) 
{
  return exp(LnNchooseK(r1, a)+LnNchooseK(r2, c1-a)-LnNchooseK(r1+r2, c1));
}

/* ------------------------------------------------------------------------- */

double FisherExact(float a, float b, float c, float d)
{
  double r1, r2, c1, minA, maxA;
  double pvalue = 0, p0 = 0, pfinal = 0;
  int i;

  r1 = (a + b);
  r2 = (c + d);
  c1 = (a + c);
  
  minA = (c1-r2 > 0 ? c1-r2 : 0);
  maxA = (r1 > c1 ? c1 : r1);
  
  pfinal = p0 = HypergeometricPMF(a,r1,r2,c1);
  i = minA;
  while((pvalue = HypergeometricPMF(i,r1,r2,c1)) <= p0 && i < a){
    pfinal += pvalue;
    i++;
  }
  i=maxA;
  while((pvalue = HypergeometricPMF(i,r1,r2,c1)) <= p0 && i > a){
    pfinal += pvalue;
    i--;
  }
  return (pfinal);
}

/* ------------------------------------------------------------------------- */

double chisq(float a, float b, float c, float d)
{
  double Ea, Eb, Ec, Ed, r1, r2, c1, c2, n, t;

  r1 = a+b;
  r2 = c+d;
  c1 = a+c;
  c2 = b+d;
  n = r1+r2;
  Ea = r1*c1/n;
  Eb = r1*c2/n;
  Ec = r2*c1/n;
  Ed = r2*c2/n;
  
  t = powf(a-Ea,2)/Ea;
  t += powf(b-Eb,2)/Eb;
  t += powf(c-Ec,2)/Ec;
  t += powf(d-Ed,2)/Ed;
  return(chisqPval(t));
}

/* ------------------------------------------------------------------------- */

float chisqPval(double t)
{
  int i;
  double q[100], p[100];
  q[0]=0.0001570879; p[0]=0.99; q[1]=0.0006284502; p[1]=0.98; q[2]=0.001414383; p[2]=0.97; 
  q[3]=0.002515382; p[3]=0.96; q[4]=0.00393214; p[4]=0.95; q[5]=0.005665552; p[5]=0.94; 
  q[6]=0.007716716; p[6]=0.93; q[7]=0.01008693; p[7]=0.92; q[8]=0.01277771; p[8]=0.91; 
  q[9]=0.01579077; p[9]=0.9; q[10]=0.01912805; p[10]=0.89; q[11]=0.02279170; p[11]=0.88; 
  q[12]=0.0267841; p[12]=0.87; q[13]=0.03110785; p[13]=0.86; q[14]=0.03576578; p[14]=0.85; 
  q[15]=0.04076098; p[15]=0.84; q[16]=0.04609676; p[16]=0.83; q[17]=0.05177672; p[17]=0.82; 
  q[18]=0.05780468; p[18]=0.81; q[19]=0.06418475; p[19]=0.8; q[20]=0.07092134; p[20]=0.79; 
  q[21]=0.07801912; p[21]=0.78; q[22]=0.08548308; p[22]=0.77; q[23]=0.09331851; p[23]=0.76; 
  q[24]=0.1015310; p[24]=0.75; q[25]=0.1101266; p[25]=0.74; q[26]=0.1191116; p[26]=0.73; 
  q[27]=0.1284927; p[27]=0.72; q[28]=0.1382770; p[28]=0.71; q[29]=0.1484719; p[29]=0.7; 
  q[30]=0.1590854; p[30]=0.69; q[31]=0.1701258; p[31]=0.68; q[32]=0.1816021; p[32]=0.67; 
  q[33]=0.1935236; p[33]=0.66; q[34]=0.2059001; p[34]=0.65; q[35]=0.2187422; p[35]=0.64; 
  q[36]=0.2320608; p[36]=0.63; q[37]=0.2458676; p[37]=0.62; q[38]=0.2601749; p[38]=0.61; 
  q[39]=0.2749959; p[39]=0.6; q[40]=0.2903443; p[40]=0.59; q[41]=0.3062346; p[41]=0.58; 
  q[42]=0.3226825; p[42]=0.57; q[43]=0.3397042; p[43]=0.56; q[44]=0.3573172; p[44]=0.55; 
  q[45]=0.3755398; p[45]=0.54; q[46]=0.3943916; p[46]=0.53; q[47]=0.4138933; p[47]=0.52; 
  q[48]=0.4340671; p[48]=0.51; q[49]=0.4549364; p[49]=0.5; q[50]=0.4765263; p[50]=0.49; 
  q[51]=0.4988633; p[51]=0.48; q[52]=0.521976; p[52]=0.47; q[53]=0.5458947; p[53]=0.46; 
  q[54]=0.5706519; p[54]=0.45; q[55]=0.5962824; p[55]=0.44; q[56]=0.6228235; p[56]=0.43; 
  q[57]=0.6503152; p[57]=0.42; q[58]=0.6788007; p[58]=0.41; q[59]=0.7083263; p[59]=0.4; 
  q[60]=0.738942; p[60]=0.39; q[61]=0.7707019; p[61]=0.38; q[62]=0.8036645; p[62]=0.37; 
  q[63]=0.8378932; p[63]=0.36; q[64]=0.8734571; p[64]=0.35; q[65]=0.9104313; p[65]=0.34; 
  q[66]=0.9488978; p[66]=0.33; q[67]=0.9889465; p[67]=0.32; q[68]=1.030676; p[68]=0.31; 
  q[69]=1.074194; p[69]=0.3; q[70]=1.119621; p[70]=0.29; q[71]=1.167090; p[71]=0.28; 
  q[72]=1.216747; p[72]=0.27; q[73]=1.268757; p[73]=0.26; q[74]=1.323304; p[74]=0.25; 
  q[75]=1.380594; p[75]=0.24; q[76]=1.440861; p[76]=0.23; q[77]=1.504371; p[77]=0.22; 
  q[78]=1.571426; p[78]=0.21; q[79]=1.642374; p[79]=0.2; q[80]=1.717618; p[80]=0.19; 
  q[81]=1.797624; p[81]=0.18; q[82]=1.882943; p[82]=0.17; q[83]=1.974226; p[83]=0.16; 
  q[84]=2.072251; p[84]=0.15; q[85]=2.177959; p[85]=0.14; q[86]=2.292505; p[86]=0.13; 
  q[87]=2.417321; p[87]=0.12; q[88]=2.554221; p[88]=0.11; q[89]=2.705543; p[89]=0.1; 
  q[90]=2.874373; p[90]=0.09; q[91]=3.064902; p[91]=0.08; q[92]=3.28302; p[92]=0.07; 
  q[93]=3.537385; p[93]=0.06; q[94]=3.841459; p[94]=0.05; q[95]=4.217885; p[95]=0.04; 
  q[96]=4.709292; p[96]=0.03; q[97]=5.411894; p[97]=0.02; q[98]=6.634897; p[98]=0.01;  

  i = 0;
  while(i <= 98 && q[i++] < t);
  return(p[i-1]);
}

/* ------------------------------------------------------------------------- */

double checkSig(float a, float b, float c, float d)
{
  float r1, r2, c1, c2, n;
  r1 = a+b;
  r2 = c+d;
  c1 = a+c;
  c2 = b+d;
  n = r1+r2;
  if(((r1 < r2 ? r1 : r2) * (c1 < c2 ? c1 : c2) / n) < 5000)
    return(FisherExact(a, b, c, d));
  else
    return(chisq(a, b, c, d));
}

/* ------------------------------------------------------------------------- */

void helpConvert()
{
  printf("PROGRAM:       convertSFS_CODE\n");
  printf("DEVELOPER:     Ryan D. Hernandez\n");
  printf("RELEASE DATE:  ???\n\n");
  printf("DESCRIPTION:\n");
  printf("\tConvert simulated data output from sfs_code to a more useful\n");
  printf("\tform.  Examples are the overall site-frequency spectrum, or\n");
  printf("\tthe input into another program such as structure, or the\n");
  printf("\tfasta sequence alignment.\n\n");
  printf("USAGE:  ./convertSFS_CODE <input_file> [options]\n\n");
  
  printf("OPTIONS:\n");
  
  printf("\t--help (-h)\n");
  printf("\t\tPrint this help menu\n");
  
  printf("\t--alignment (-a) [A] [F <filename>] [P <#pops> <P1>..<P#pops>] \n");
  printf("\t\t  [I <#indiv> <I1>..<I#indiv>] [P.I <#comb> <P1.I1>..<Pn.In>]\n");
  printf("\t\t  [L <#loci> <L1>..<L#loci>] [ITS <#its>]\n");
  printf("\t\tPrint the alignment in fasta format.  Use \"A\" option to print");
  printf("\n\t\tancestral sequences.  Use \"F\" option\n");
  printf("\t\tto print the output to the file <filename>.  Use \"P\" option\n");
  printf("\t\tto print only specific populations, or \"I\" option for\n");
  printf("\t\tspecific individuals.  Use \"P.I\" option to print only\n");
  printf("\t\tspcific individuals from specific populations, and use \"L\"\n");
  printf("\t\tto print output from specific loci.  Use \"T\" option to\n");
  printf("\t\tprint only first <#its> iterations.  Note: cannot use\n");
  printf("\t\tboth \"P\" and \"I\" options (use \"P.I\" instead).  By\n");
  printf("\t\tdefault, a fasta style alignment of all individuals across\n");
  printf("\t\tall loci and iterations will be printed to the screen.\n\n");

  printf("\t--MK (-m) <ingroup> <outgroup> [OBS] [F [a] <file>]\n");
  printf("\t\t   [L <#loci> <L1>..<L#loc>] [ITS [#its]] [OGSS <OG_size>\n");
  printf("\t\tExtract McDonald-Kreitman (MK) tables for comparing <ingroup>\n");
  printf("\t\tto <outgroup>.  By default, will print the MK table for the\n");
  printf("\t\ttrue data (saved during simulation) summed over all loci &\n");
  printf("\t\titerations, with each entry space delimited.  Using \"OBS\"\n");
  printf("\t\toption will print the observed MK table based on parsimony\n");
  printf("\t\t(by default a single sequence from the outgroup will be used\n");
  printf("\t\tfor divergence, but specifying OG_size can change this).\n");
  printf("\t\tUse \"F\" to print output to the file <file> (include \'a\'\n");
  printf("\t\tto append to file instead of over-writing it).  Use \"L\" to\n");
  printf("\t\tprint specific loci (<#loci> = -1 prints seperate MK tables at");
  printf("\n\t\teach locus).  To print each iteration indepdently, use \n");
  printf("\t\t\'ITS\' option (to print only a subset of iterations, \n");
  printf("\t\tinclude the number of iterations [#its].  Using the \"OGSS\"\n");
  printf("\t\toption sets the outgroup size to <OG_size>, such that any\n");
  printf("\t\tmutation carried by the first <OG_size> chromosomes in the\n");
  printf("\t\toutgroup will be considered a fixed difference (setting\n");
  printf("\t\t<OG_size> = -1 is uses the entire sample size).\n");
  printf("\t\tEach MK table is printed in the following order:\n");
  printf("\t\t<PS> <PN> <FS> <FN>, where \'P\' refers to polymorphic,\n");
  printf("\t\t\'F\' fixed, \'S\' synonymous and \'N\' is nonsynonymous.\n\n");

  printf("\t--SFS (-s) [OBS <P1> <P2> [OG_size]] [T <0/1/2>] [F [a] <file>]\n");
  printf("\t\t   [L <#loci> <L1>..<L#loc>] [P <#pops> <P1>..<P#pops>]\n");
  printf("\t\t   [ITS [#its]]\n");
  printf("\t\tExtract Site-Frequency Spectra.  By default, will print the\n");
  printf("\t\tSFS for the true data (saved during simulation) summed\n");
  printf("\t\tover all loci & iterations, each entry space delimited and\n");
  printf("\t\tpopulations separated by a tab.  Using \"OBS\" option will\n");
  printf("\t\tprint the observed SFS for population <P1> using <P2> as an \n");
  printf("\t\toutgroup to parsimoniously determine the ancestral state of\n");
  printf("\t\teach SNP (arbitrarily chooses first sequence from outgroup as\n");
  printf("\t\treference unless OG_size is specified).  Use \"T\" option to\n");
  printf("\t\tchange the type of SFS extracted:\n");
  printf("\t\t0=synonymous, 1=nonsynonymous, 2=all mutants (default).\n");
  printf("\t\tUse \"F\" to print output to the file <file> (include \'a\'\n");
  printf("\t\tto append to file instead of over-writing it).  Use \"L\" to\n");
  printf("\t\tprint specific loci (<#loci> = -1 prints seperate SFS for\n");
  printf("\t\teach locus).  Use \"P\" to print the SFS only from specific\n");
  printf("\t\tpopulations (only valid for true data, i.e. not when using\n");
  printf("\t\tOBS flag).  To print each iteration indepdently, use \'ITS\'\n");
  printf("\t\tflag (to print only a subset of iterations, include the\n");
  printf("\t\tnumber of iterations [#its].\n");  
}

/* ------------------------------------------------------------------------- */

