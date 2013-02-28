struct mutArray{
  struct history **muts; /* array of events */
  long numMuts;          /* total number of events */
  long mutIndex;         /* next putatively free mutation */
};

struct event{
  char ancNuc;           /* ancestral nucleotide */
  char derNuc;           /* derived nucleotide */
  char nonsyn;           /* indicator whether or not a nonsynonymous change */
  char axy;              /* chromosome type (autosome, X, or Y) */
  char CpG;
  char fiveP;
  char threeP;
  int ancAA;
  int derAA;
  long gen;              /* generation mutation arose */
  long *genFix;          /* generation mutation fixed in each population */
  long *genDead;         /* generation mutation lost from each population */
  unsigned long site;
  double fit;
  int nSites;           /* number of sites inserted/deleted */
  char *nucs;            /* nucleotides inserted */

  char free;             /* indicates whether information free */
  char *fixed;           /* indicates whether mutation is fixed in 1+ pops */
  long index;            /* index in mutation array */
  long *nextHap;         /* next potentially free haplotype ptr in hapFreq */
  long *numCarriers;     /* number of haplotypes carrying mut in each pop */
  long *maxHaps;         /* maximum number haplotypes allocated for hapFreq */
  long ***hapFreq;       /* array of pointers to hap. freq. in each pop */
  int checkPop;          /* only check frequency once per generation */
  long checkGen;         /* different pops require new check */
  int checkStep;         /* differnt steps require different counters */
};

struct history{
  struct event *event;
  struct history* Ltree;
  struct history* Rtree;
  struct history* Parent;
};

struct POPparams{
  int ALIVE;         /* indicator (0/1), whether this population is alive */
  long N;            /* population size */
  double Nt;         /* actual population size (nec. for exp. growth) */
  float pFEMALES;    /* proportion of females in population (ladies first) */
  long MALES;        /* the first male in the population */
  float pMaleMig;    /* probability that a male migrates out */
  float pMaleRec;    /* probability that the paternal chrm is recombinant */
  double popAlpha;   /* exponential population growth rate */
  long tauAlpha;     /* time logistic growth started */
  long K;            /* Final population size for logistic growth */
  long P0;           /* population size when logistic/exp. growth begins */
  int GenEffect;     /* generation time effect */
  int SS;            /* size of sample for pop */
  double THETA;      /* population scaled mutation rate */
  double INSRATE;    /* rate of insertions */
  double DELRATE;    /* rate of insertions */
  double INDELlength;/* mean length of indel */
  double longINSRATE;    /* rate of insertions */
  double longDELRATE;    /* rate of insertions */
  double longINDELlength;/* mean length of indel */
  double INVRATE;        /* rate of invertions */
  double INVlength;      /* mean length of inversions */
  double RHO;        /* population scaled recombination rate */
  double fGC;        /* rate of gene conversion / recombination */
  double GCtract;    /* mean length of GC tract ~Geom(1/GCtract) */
  double BGC;        /* biased gene conversion parameter */
  double SELF;       /* selfing rate (0,1) */
  int neutpop;       /* {0,1} if allowing selection or not */
  int *selDistType;  /* neutral=0; discrete=1; mixed_gamma=2; normal=3; else=4*/
  double propSelLoci;/* proportion of loci under selection model */
  double *GAMMA;     /* pop. scaled selection coef. (discrete case) */
  double *ProbPos;   /* probability of advantageous mutation */
  double *ProbDel;   /* probability of deleterious mutation */
  double *ProbNeut;  /* probability of neutral mutation */
  double *ProbNegGamma; /* probability of negative gamma distribution */
  double *alphaN;       /* parameter for Gamma distn (deleterious case) */
  double *lambdaN;      /* parameter for Gamma distn (deleterious case) */
  double *alphaP;       /* parameter for Gamma distn (advantageous case) */
  double *lambdaP;      /* parameter for Gamma distn (advantageous case) */
  double *normMean;     /* mean of normal distribution of sel. effects */
  double *normVar;      /* variance of normal distribution of sel. effects */
  float *f0;            /* proportion of non-lethal nonsynonymous mutations */
  unsigned int RateClassSites; /* number of mutation rate classes w/i locus */
  double RateParamSites;       /* Gamma rate param for site classes */
  unsigned int RateClassLoci;  /* number of mutation rate classes among loci */
  double *RateParamLoci;       /* Gamma rate parm for loci classes */

  double KAPPA;      /* transition/transversion bias, if substMod=2|3 */
  double KAPPACpG;   /* transition/transversion bias, for CpGs */
  double PSI;        /* non-CpG rejection rate, if substMod=1|3 */
  double GTR[6];     /* GTR model parameters */
  double baseFreq[4];  /* nucleotide base frequencies */
  double transProb[4][4]; /* nucleotide transition probabilities */
};

struct GLOBparams{
  int iter;          /* current iteration */
  int ITER;          /* number of iterations to simulate */
  long int seed;     /* random number seed */
  int NPOP;          /* total number of populations to be simulated */
  int NPOPDEF;       /* total populations defined */
  long maxN;         /* maximum size of any population */
  int INFSITES;      /*indicates whether (1) or not (0) to use inf. sites mod.*/
  int P;             /* ploidy */
  int tetraType;     /* 0=>autotetraploid; 1=>allotetraploid */
  long R;            /* number of loci */
  char EQL;          /* indicates that the loci are all the same length */
  unsigned long *L;  /* length of each locus */
  double *sumL;      /* sum of lengths of loci */
  double *sumLL;     /* sum of lengths of loci with linkage */
  double totalNucs;  /* total number of nucleotides (excluding possible Ns) */
  double *LINK;      /* linkage between loci */
  long mapPoints;    /* total number of points in recombination map */
  double *recMap;    /* recombination map */
  long *recMapPos;   /* corresponding base position of rec. rate point */
  char LinkDistType; /* either 0 or 1 for physical versus recom. frac. */
  char *ANNOTATE;    /* annotation of each locus: C or N for coding or non */
  char *SEX;         /* whether or not locus is sex chromosome */
  int substMod;      /* 0:JC69, 1:JC69+CpG, 2:K2P, 3:K2P+CpG, 4:ZG03, 5:C-D */
  double **mig_mat;  /* migration rate matrix */
  double BURN;       /* burn-in time for simulation */
  double BURN2;      /* burn-in time for simulation */
  int PRINTSEQ;      
  int KEEPFIXED;     /* whether or not to keep substitutions during burn-in */
  int FIXSTOP;       /* fix in frame stop codons in imported sequence */
  int ADDITIVE;      /* Additive (1) versus multiplicative (0) */
  int FITTEST;       /* 0=>uniform; 1=>demographic events use rel fit. */
  int TRACKANC;      /* track ancestry if 1, else do not */
  int TRACKTRAJ;     /* track trajectory if 1, else do not */
  int autoRestart;   /* automatically restart if mutation is fixed/lost */
  long trajPop;      /* population to track trajectory of */
  long trajLoc;      /* locus to track trajectory of */
  long trajSite;     /* site to track trajectory of */
  float trajMinFreq; /* minimum frequency at time of sampling to keep output */
  float trajMaxFreq; /* maximum frequency at time of sampling to keep output */
  long trajMaxReps;  /* maximum replicate iterations to attempt */
  int **SKIPSITES;   /* only used with --mutation, prevents unwanted mutants */
  int PRINTGEN;      /* print the generation to screen throughout simulation */
};

struct EvolEvents{
  double tau;        /* population scaled time that event takes place */
  int popi;          /* default/original population */
  int popj;          /* secondary/new population */
  int nAncPops;      /* number of ancestral populations for admixture */
  int *ancPops;      /* ancestral populations for admixture */
  float *maleFreqs;  /* proportion of males from each ancestral population */
  float *femaleFreqs;/* proportion of females from each ancestral population */
  float nu;          /* magnitude of parameter size change */
  float freq;        /* frequency of domestication allele */
  long locus;        /* desired locus for event */
  long site;         /* desired site for event */
  double gamma;      /* desired selection coefficient (default 0) */
  int eventType;     /* types of events:
			(0) split population i to found j (i->j)
			(1) create domesticated species
			(3) change population size
			(4) change ploidy
			(5) kill population popi (popi = -1 kills all)
			(6) change 1 or more POPparams
			(7) change 1 or more GLOBparams (excluding P)
			(8) join populations (admixture)
			(9) add mutation event
		     */
  int parIndex;      /* parameter index for changing POPparams:
			(0) specific population size
			(1) population growth(decay) rate (exponential)
			(2) logistic growth rate
			(3) THETA
			(4) RHO
			(5) selfing rate
			(6) neutPop  #implemented as (7) by changing selDistType
			(7) selDistType
			(8) KAPPA
			(9) PSI
			(10) f0
			(11) pMaleMig
			(12) genEffect
			(13) indel
			(14) longIndel
			(15) inversion
			(16) gene conversion
			(17) pMaleRec
			(18) baseFreq

			parameter index for changing GLOBparams:
			(0) migration matrix
			(1) track trajectory
		     */
  struct POPparams newP;
  struct GLOBparams newG;
  struct EvolEvents *nextEvent;
};

struct mutStorage{          /* a binary tree to store mutation info */
  int *pops;                /* populations that carry mutation */
  long int numCarry;        /* total number of carriers. */
  long int *chrs;           /* the specific chromosome carrying mut */
  long locus;               /* locus on which mutation lies */
  struct event *event;      /* mutation information */
  struct mutStorage *Rtree; /* muts with larger locus and/or site */
  struct mutStorage *Ltree; /* muts with smaller locus and/or site */
  struct mutStorage *Parent;/* parent, only used in convertSFS_CODE.c */
}; /* for fixed differences, chrs will be -1 */

struct population{
  long maxchr;              /* highest haplotype number across loci */
  long *nextchr;            /* first free chromosome in pop at each locus */
  long gen;                 /* generations experienced by population */
  long **parentLoc;         /* location of haplotype in memory */
  long **numCopies;         /* number of carriers of haplotype/locus */
  long **chrGen;            /* generation haplotype went to freq=0 */
  long ***extremeMuts;      /* min/max site of mutation carried by hap/locus */
  int **polySites;          /* indicates #muts seg at given site/locus */
  double ***hpSUM;          /* sum of hit probabilities across loci/haplotype */
  double **indfit;          /* fitness of each haplotype */
  char **conSeq;            /* consensus sequence/locus */
  int ***ancestry;           /* a matrix of ancestries of each site */
  struct history **fixed;   /* tree of fixed differences */
  struct history ***BigHead;/* splay tree of mutations carried by each hap */
};

struct locRec{
  long id;        /* individual ID: [0,N) */
  long numEvents; /* total number of GC+rec events */
  long *loci;     /* array of loci containing events */
  long *pos;      /* array of positions of events */
  int  *type;     /* 0=rec, 1=GC, 2=both */
};

struct interRec{
  long id;
  long numEvents;
  long *loci;
};

struct deadList{
  long chr;    /* chr (i.e., in parentLoc) */
  long loc;    /* locus */
  struct deadList *next; /* doubly-linked list */
  struct deadList *prev; /*  */
};

void setDefaultPPARS(const int pop);

void setDefaultGPARS();
/* set all parameters to their default values */

void Initialize(const double *ST);
/* generates sequence from stationary distribution of mutation model */

int AA(const char *c);
/*  Takes in an array of integers and returns an integer: */
/*  0=stop codon, 1-20 = amino acids. */

void bign(double **mutClass_site, double *mutClass_locus, int *selClass,
	  const double *fitQuant, long *numMutSeg);
/*  Carries out one realization of the evolutionary process */

void partition(double *part, double *relFit, const int *selClass,const int pop);
/*  Partitions (0,1), so that each individual gets a segment */
/*  proportional to their relative fitness. */

void Selecting(double *relFit, long *pars, const int pop);
/*  Draws parental chromosomes from the preceeding generation */
/*  to found the new generation */

void NextGen(long **pars, const long popsize, int pop);
/*  Generates new set of linked lists after a round of random mating */

void distribRec(const int pop, struct locRec **within, long *numX1,
		struct interRec **between, long *numX2);

void getRecombinant(long ind, int p, long oldPar, long *newPar,
		    struct locRec *within, long iW, long *iW2);

void recombine(const int pop);
/*  Draws a poisson number (mean PNr) of recombination events in the entire */
/*  population each generation, then randomly distributes those events across */
/*  individuals */

void getLargestHist(struct history *tree, struct history **node);

void getSmallestHist(struct history *tree, struct history **node);

void recombineHists(const unsigned long site, struct history **tree1,
		    struct history **tree2, const long xreg, int pop,
		    const long x1, const long x2);

void swapHapFreq(struct history *tree, const int pop, const long xreg,
		 const long chrTo, const long chrFrom);

void mutate(double **mutClass_site, double *mutClass_locus, const int *selClass,
	    const double *fitQuant, const int pop);
/*  Draws a poisson number (mean theta/2) of mutations to enter the population*/
/*  each generation, then randomly distributes mutations across individuals */

void indels(const int pop, const double *fitQuant);

void inversions(const int pop, const double *fitQuant);

char MutantNuc(const char *nuc, const int CpG, const int pop);
/*  Takes in as argument the nucleotide chosen to mutate, and returns */
/*  the chosen mutant nucleotide. */

double NewFits(const double *fitQuant, const int pop,
	       const long locus);
/*  returns a fitness value corresponding to discrete, mixed gamma, and normal*/
/*  distributions. */

void Split(const int fromPop, const int toPop);
/*  Generates divergent population. */

void admixture(struct EvolEvents *devents, const long PNanc, 
	       const int *selClass);

void Migrate(const int pop);

void copyNchrs(long femalesIN, long malesIN, const int popTo, const int popFrom,
	       long femaleToStart, long maleToStart);

void ChangePopSize(double *relFit, const long newSize, const int pop, 
		   int *survive, const int fittest);

void domesticate(int **DomestAllele, int *survive, int pop);

void checkForErrors(int pop, char *errMessIn);

void PrintErrorStats(const int pop, const long r);

int NucToNum(const char *c);

char NumToNuc(const int c);

int AAToNum(const char c);

char NumToAA(const int c);

char ToCGTA(const char c);

char FromCGTA(const char c);

char TsNuc(const char c);

char TvNuc(const char c);

void GenContextRate(double *ST, double *fitQuant);

void GenInitFreqMat(float *f);

void exitNOW(const char *s);

char retNucHistory(const unsigned long site, struct history **BigHead, 
		   const long locus, int pop, int SPLAY);

void getNucHistory(const unsigned long site, char *retNuc, 
		   struct history *BigHead, struct history **root, 
		   const long locus, int pop, int SPLAY);

void addHistoryNode(const struct event *data, struct history **BigHead, 
		    struct history **MOM, const int pop, const long chr,
		    const long locus,  const int k);

void copyHistoryNode(struct history **to, struct history **from,
		     struct history **MOM, const int pop, const long chr,
		     const long locus);

void checkCarryNum(struct history *tree, long **checkCarry, long numCopies,
		   long *badSite);

void emptyMut(const long mutreg, const int pop);

void updateMutFreqs();

void adjustCarriers(int pop);

void freeDead(const int p, const int step);

struct history* popTrash();

void pushTrash(struct history **tree);

void freeTrash(struct history *tree);

void freeHistory(struct history *tree, const int pop, const long chr,
		 const long locus, int step, const int inc);

void checkForFixed(const int pop, long chr, long locus, struct history *tree,
		   int inc);

void carryMut(struct history *tree, long index, int *foo);

void freeFixedHistory(struct history *tree, const int pop);

void checkStillFixed(struct history *tree, int popTo, int popFrom, long locus,
		     long chr);

void clearMutList(long locus);

void addToStorage(struct mutStorage **storage, struct event *data,
		  const int pop, const long chr, const long locus,
		  const long offset, const long genFix);

void storeInfo(struct mutStorage **storage, struct history *data, const int pop,
	       const long chr, const long locus, const long offset,
	       const long genFix);

void regeneratePOP();

void getPolySites(struct mutStorage *storage, long loc, long *Nsites, 
		  long **sites);

struct history* freeHistoryNode(struct event *event, struct history *tree,
				const long locus);

void setHistoryParents(struct history *tree);

void rotateHistoryNode(struct history *node, struct history **tree, int pop);

void splayHistory(struct history *node, struct history **tree, int pop);

void PrintHistTree(struct history *BigHead, long *i, int pop);

void getNewHitProbsReg(double **hpTMP, struct history **BigHead,
		       const long locus, const int pop, const long first,
		       const long last);

void getMutSite(unsigned long *site, double *hpSUM, const long locus);

long invCDF(float rn, double *dist, long beg, long end);

void getFitTree(struct history *BigHead, double *fit, const int pop,
		const long locus);

void copyHistory(struct history **to, struct history **from, const long locus,
		 const int popTo, const long chrTo, const int popFrom, 
		 const int updatePS, const int fixed, const int cpNC);

void copyPartialHistory(struct history **to, struct history **from,
			const long min, const long max, const int pop,
			const long chr, const long locus);

void convertTract(struct history **to, struct history **from, const long min,
		  const long max, const int pop, const long chr,
		  const long locus, const long altPar);

void addToTract(struct history **to, struct history **from, const long min,
		const long max, const int pop, const long chr,
		const long locus, const long altPar);

void reconvertTract(struct history **to, const long donor, const long altPar,
		    const long min, const long max, const int pop,
		    const long chr, const long locus);
 
void copyFixed(struct history *tree, int pop, long chr, long locus);

void memcpyEVENT(struct event *to, const struct event *from, const int pop);

void cpyToFixed(struct history **tree, struct history **fixed, const int pop);

void createNewHistory(const long chr, const long locus,int pop,int IBD, int FZ);

void checkHistTreeOrder(struct history *head, long min, long max, int pop);

void checkHitProbs(const double *hpTMP, const char *err, const int pop,
		   const long chr, const long locus);

void parseCommandSFSCODE(int argc, char *argv[], int it);

void helpMenu(){
  printf("PROGRAM:  \tSelection on Finite Sites under COmplex Demographic \
Effects\n\t\t(SFS_CODE)\n");
  printf("DEVELOPER: \tRyan D. Hernandez\n");
  printf("COPYRIGHT (C):\t2007\n");
  printf("DESCRIPTION:\n");
  printf("\tSFS_CODE simulates population genetic data under a wide range of\n\
\tdemographic and selective models.\n\n");

  printf("USAGE: ./sfs_code <NPOP> <ITER> [--option <argument(s)>]...\n");
  printf("\twhere <NPOP> is the total number of populations to be simulated\n");
  printf("\tand <ITER> is the number of iterations you would like to run.\n\n");
  
  printf("OPTIONS:\n");

  printf("\t(Note: arguments in '<..>' are required, in '[..]' are optional.\n\
\tIncluding the optional argument '[P <pop>]' will set the\n\
\tparameters just for population <pop>.  Both short name [single -,\n\
\tsingle letter] and long name [double --, full name] are provided,\n\
\tbut you only need to enter one of them.)\n\n");
  
  printf("-A  --noSeq\n\
\tdon't print ancestral sequence\n\n\
-a  --annotate  [F <filename>] [<a1> [<a2>...<aR>] [R]]\n\
\tindicate whether each locus is coding/non-coding\n\n\
-B  --BURN  <burn>\n\
\tset initial burn-in length, generations/P/N\n\n\
-b  --BURN2  <burn>\n\
\tset burn-in length of subsequent iterations > 1\n\n\
-C  --PSI  [P <pop>] <psi>\n\
\tset the CpG mutation bias\n\n\
-c  --constraint  [P <pop>] [L <locus>] <f0>\n\
\tset the non-lethal mutation rate\n\n\
-TD  <t> <i> <j> <allele_freq> <N> [locus]\n\
\tonly used with -T to create a domesticated population\n\n\
-Td  <t> [P <pop>] <v>\n\
\tonly used with -T to set the demographic effects\n\n\
-TE  <t> [pop]\n\
\tonly used with -T to end the simulation for a population at time t\n\n\
-e  --errfile  [a] <file> \n\
\tprint error to specific file\n\n\
-F  --popFreq  [a] <file>\n\
\tcreate file with population & sample frequencies for each event\n\n\
-f  --propFemale [P <pop>] <pf>\n\
\tset the proportion of females in a population\n\n\
-G  --GenEffect  <pop> <G>\n\
\tset the generation effect for a population\n\n\
-Tg  [P <pop>] <alpha>\n\
\tonly used with -T to set the exponential population growth rate\n\n\
-H  --geneConversion  [P <pop>] [B <BGC>] <f> <lambda> \n\
\tset parameters for gene conversion\n\n\
-h  --help \n\
\thelp menu\n\n\
-I  --INF_SITES \n\
\tturn on infinite-sites model [only at single time!!]\n\n\
-i  --self  [P <pop>] <s>\n\
\tset the selfing [not really inbreeding] rate\n\n\
-J  {NOT USED...}\n\n\
-j  {NOT USED...}\n\n\
-K  --KAPPA  [P <pop>] <kappa> \n\
\tset transition transversion rate\n\n\
-Tk  [P <pop>] <K> <r>\n\
\tused with -T to implement logistic growth rate\n\n\
-L  --length  <nloci> <L1> [<L2>...<Ln>] [R]\n\
\tset sequence lengths and number of loci\n\n\
-l  --linkage  <p/g> <d1> [<d2>...<dn-1>] [R]\n\
\tset linkage between adjacent loci\n\n\
-M  --substMod  <mod> [args]\n\
\tset the substitution model\n\n\
-m  --mig_mat {see documentation...}\n\
\tset the migration rates to and from populations\n\n\
-N  --popSize  [P <pop>] <size>\n\
\tset the effective size of a population\n\n\
-n  --sampSize  [P <pop>] <SS1> [<SS2>...<SSNpops>]\n\
\tset the number of individuals sampled from a population\n\n\
-O  {NOT USED...}\n\n\
-o  --outfile  [a] <file>\n\
\tprint output to specific file\n\n\
-P  --ploidy  <P>\n\
\tset the ploidy, P=1 is haploid, P=2 is diploid, P=4 is tetraploid\n\n\
-p  --tetraType  <0/1>\n\
\tset the type of tetraploid [allo- or auto-]\n\n\
-Q  {NOT USED...}\n\n\
-q  {NOT USED...}\n\n\
-R  --ReportBurnFixed\n\
\treport mutations that fixed during burn-in\n\n\
-r  --rho  [P <pop>] [F <filename>] <rho>\n\
\tset the recombination rate for a population\n\n\
-TS  <t> <i> <j>\n\
\tused with -T to split the population at a given time\n\n\
-s  --seed  <value>\n\
\tset seed for random number generator\n\n\
-T*\n\
\tused in combination with other parameters for timed effects\n\n\
-t  [P <pop>] <theta>\n\
\tset the value of theta for a population\n\n\
-U  --longIndel  [P <pop>] <INSrate> <DELrate> <mean_length>\n\
\trate of insertion/deletions with Poisson length distribution\n\n\
-u  --indel  [P <pop>] <INSrate> <DELrate> <mean_length> \n\
\trate of insertion/deletions with geometric length distribution\n\n\
-V  --rateClassSites  [P <pop>] <nclasses> <alpha> \n\
\tmutation rate variation among SITES\n\n\
-v  --rateClassLoci  [P <pop>] <nclasses> <alpha>\n\
\tmutation rate variation among LOCI\n\n\
-W  --selDistType  [P <pop>] [L <locus>] <type> [args] \n\
\tset distribution of selective effects (often denoted w in pop.gen.)\n\n\
-w  --neutPop  <pop>\n\
\tallow a population to be neutral\n\n\
-X  {NOT USED...}\n\n\
-x  --sex  <x1> [<x2 >..<xR>] [R]\n\
\tindicate whether each locus is autosomal/sex\n\n\
-Y  --pMaleRec  [P <pop>] <p>\n\
\tproportion of recombinants occuring in male meioses vs female.\n\n\
-y  --pMaleMig  [P <pop>] <pmale>\n\
\tproportion of migrants that are male\n\n\
-Z  --additive \n\
\tspecify an additive model of selective effects instead of multiplicative\n\n\
-z  --inversions  [P <pop>] <INVrate> <mean length>\n\
\trate of inversions with Poisson length distribution\n\n");

  exit(0);
}

