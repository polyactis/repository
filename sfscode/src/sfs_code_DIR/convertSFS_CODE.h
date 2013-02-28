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
  long genFix;          /* generation mutation fixed in each population */
  unsigned long site;
  double fit;
  long nSites;            /* number of sites inserted/deleted */
  char *nucs;            /* nucleotides inserted */
};

struct mutStorage{ /* a binary tree to store mutation info */
  int *pops;                /* populations that carry mutation */
  long int numCarry;    /* total number of carriers. */
  long int *chrs;       /* the specific chromosome carrying mut */
  long locus;       /* locus on which mutation lies */
  struct event *event;      /* mutation information */
  struct mutStorage *Rtree; /* muts with larger locus and/or site */
  struct mutStorage *Ltree; /* muts with smaller locus and/or site */
  struct mutStorage *Parent;/* parent, only used in convertSFS_CODE.c */
}; /* for fixed differences, chrs will be -1 */


void printAlign(FILE *out, int numPOPS, int *pops, int numINDIV, int *indivs,
		int numPI, int numLOCI, int *loci, int **keep, int *keepLOCI,
		int printITS, long *maxSeqLen, long numITS, int totLOCI,
		int TOTPOPS, int *INDperPOP, int printANC, char ***ancSeq,
		struct mutStorage **data);

void printAges(struct mutStorage *storage, char TYPE, int pop, double maxGen,
	       long PNa, FILE *out, long IT);

void printSFS(char TYPE, int TRUE, int its2print, int numPOPS, int numLOCI,
	      int *numOG, int *pops, int *loci, int **og, FILE *out, 
	      struct mutStorage **data, int TOTPOPS, int totLOCI, 
	      int *INDperPOP, long numITS, char ***ancSeq, int OGsize, int PSEX,
	      int *MALES, int MONLY, int ploid, int **keep);

void printMK(int TRUE, int ING, int OUTG, int numLOCI, int *loci, int its2print,
	     FILE *out, struct mutStorage **data, long numITS, int totLOCI,
	     int *INDperPOP, char ***ancSeq, int OGss);

void printMS(struct mutStorage **data, FILE *outfile, long numITS, int TOTPOPS,
	     int totLOCI, long *maxSeqLen, int *INDperPOP, long seed, char TYPE);

void getZeroOne(struct mutStorage *data, char ***ZeroOne, double **positions,
		long *Nsites, int *popInds, int totchrs, long *totLen, 
		long sumLen, char TYPE, long *maxMuts, int *INDperPOP);

void printSTRUCTURE(struct mutStorage **data, char ***ancSeq, char *filename,
		    long numITS, int TOTPOPS, int totLOCI, double *LINK, 
		    double CMMB, int linkTYPE, long *maxSeqLen, int **keep,
		    int *INDperPOP);

void getSynNsSites(char *anc, char **cod, int Mc, int N, int *V, int **type);

int normCod(char *c1, char *c2);

void synVSns(char *anc, char *der, int **type);

int makeRecomb(char *anc, char *c1, char*c2, char *res);

void reduceCodons(char *anc, char **cods, int numCods, int *Mc,char ***MinCods);

void getPolySitesPop(struct mutStorage *storage, int loc, int *Nsites, 
		     long **sites, int pop, int *keep);

void getPolyFixedSites(struct mutStorage *storage, int loc, int *Nsites, 
		       long **sites);

void getObsFixedSites(struct mutStorage *storage, int loc, int *Nsites, 
		      long **sites, int pop1, int pop2, int ss1, int ss2,
		      struct mutStorage **rootStor, char *ancSeq);

char ToCGTA(const char c);

char FromCGTA(char c);

void memcpyEVENT(struct event *to, const struct event *from);

void exitNOW(const char *s);

void buildStorage(struct mutStorage **storage, struct mutStorage *tmp,
		  struct mutStorage **parent);

void printStorage(struct mutStorage *storage, int *x);

void freeStorage(struct mutStorage *storage);

void updateSeq(char **seq,struct mutStorage *storage, int pop, int chr, int loc,
	       long *numInv, long ***invBK);

void rotateStorageNode(struct mutStorage *node, struct mutStorage **tree);

void splayStorage(struct mutStorage *node, struct mutStorage **tree);

void getNucStorage(int pop, int ind, int locus, long site, 
		   struct mutStorage *storage, char *nuc, 
		   struct mutStorage **root, int toSplay);

void updateObsSFS(int *sfs, struct mutStorage **storage, int pop, int ss,
		  int locus, char TYPE, int og, char *ancSeq, int ogSS, int PSEX,
		  int *MALES, int MONLY, int ploid, int **keep, int max);

void getObsMK(int *MK, int ING, int OUTG, struct mutStorage **storage, int loc,
	      int ss1, int ss2, char *ancSeq);

void updateTrueSFS(int *sfs, struct mutStorage *storage, int pop, int ss,
		   int loc,char TYPE, int PSEX, int *MALES, int MONLY, int ploid,
		   int **keep);

void updateDerivedAlleles(int ***cnt, struct mutStorage *storage, int *ss,
			  char TYPE, int ploid, int HH);

void printFitness(struct mutStorage *storage, int pop, int nShare, int *share,
		  int private, int fixed, FILE *out);

void updatePrivShared(struct mutStorage *storage, int **ps, char TYPE,int *pops);

char NumToAA(const int c);

int AAToNum(const char c);

int AA(const char *c);

void isFixedSame(struct mutStorage *storage, int pop, int locus, 
		 struct event *ev, int ss, int *fixed);

void getTrueMK(int *MK, int ING, int OUTG, struct mutStorage *storage, int loc,
	       struct mutStorage *rootStor, int ss1, int ss2);

double LnNchooseK(int n, int k);

double HypergeometricPMF(int a, int r1, int r2, int c1);

double FisherExact(float a, float b, float c, float d);

double chisq(float a, float b, float c, float d);

float chisqPval(double t);

double checkSig(float a, float b, float c, float d);

void helpConvert();
