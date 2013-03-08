/****************************************************************************************
SimPed: A simulation program to generate haplotype and genotype data for pedigree structures.

The SimPed program generates haplotype and/or genotype data for pedigrees as follows.  All of 
the founders within the pedigree are assigned haplotypes and/or alleles conditional on the 
user specified frequencies. Once assignment is completed each founder has two haplotypes.  
Starting at the top of the pedigree structure the first offspring of the founder is randomly 
assigned one of the founder's haplotypes. The allele at the first marker from this haplotype 
is assigned to the offspring.  It is then determined based upon the genetic map whether a 
recombination event has occurred between the first and second marker loci. If with probability, 
a recombination event has occurred then at the second marker locus the allele from the founder's 
other haplotype is assigned to the offspring.  If recombination event has not occurred with 
probability (1- ) then the allele from the founder's same haplotype is assigned at the offspring's 
second marker locus.  This procedure is repeated until alleles for all markers loci have been 
assigned from one founder to their offspring.  The process is then repeated this time assigning 
alleles to the offspring from their other parent.  In this manner the haplotypes flow down the 
pedigree tree as all nonfounders are assigned haplotypes conditional on parental haplotypes.
Once all individuals within the pedigree have been assigned haplotypes, for those individuals/
marker loci for which it was specified that they are unavailable the genotypes are made unknown 
(i.e. 0 0).
********************************************************************************************/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#define MAX_PED 100		//maximum number of pedigrees
#define MAX_IND 3000		//maximum number of individuals
#define MAX_LOCUS  10000	//maximum number of loci
#define MAX_HAPLOTYPES  1000	//maximum number of haplotypes
#define MAX_ALLELE      20      //maximum number of alleles at a locus

struct genotype{
	int geno_chars[MAX_LOCUS][2];
};

typedef struct person{
	char pedId[20];
	char id[20];
	char father[20];
	char mother[20];
	int  sex;
	int  filled;
	struct genotype *geno;
}Person;

Person  *indi[MAX_IND];
int     membersOfFamily[MAX_PED];
int     seeds[3];
float   recomb_fractions[MAX_LOCUS];
int     totalInd, totalPed, totalLoci;
int     extra_col, replicates, inter_dist_type;
int     haplotype[MAX_HAPLOTYPES][MAX_LOCUS];
float   possibilities[MAX_HAPLOTYPES], allele_freq[MAX_ALLELE];
long    offset;		//the current value of the file-position indicator for parameter file
static int index;	//be reset in main() after each replicate; be increased in 
			//buildFounderGeno(), which is called by readMarkerInfo().
			//it's used to track the number of markers be processed when the program
			//assigns genotypes/haplotypes to founders in each replicate
static int familyNo;	//be initialized in main() and get increased in writeResults()
			//it's used to re-number the result pedigrees after each replicate
int 	iteration;      //the last number in the 5th line of parameter file


/************************************************/
/*                                              */
/*              randomnum	                */
/*                                              */
/************************************************/
float randomnum()  //random number generator
{
  	/*Global variables used: seeds[3]*/
  	float r;
	seeds[0] = seeds[0] % 177 * 171 - seeds[0] / 177 * 2;
	seeds[1] = seeds[1] % 176 * 172 - seeds[1] / 176 * 35; 
	seeds[2] = seeds[2] % 178 * 170 - seeds[2] / 178 * 63;
	if (seeds[0] < 0)
	    seeds[0] += 30269;
  	if (seeds[1] < 0)
    	    seeds[1] += 30307;
	if (seeds[2] < 0)
            seeds[2] += 30323;
  	r = seeds[0]/ 30269.00 + seeds[1] / 30307.00 + seeds[2] / 30323.00;
  	
	return (r - (long)r);
}

/************************************************/
/*                                              */
/*              buildFounderGeno                */
/*                                              */
/************************************************/
/*this subfunction is called by "readMarkerInfo". It determines 
the genotypes for all the founders by using random generator.*/
void buildFounderGeno(int type, int loci, int haplo_or_allele)
{
	float randnum;
	int   i, j, t, tmp;
	int   count;

	for (i=0; i<totalInd; i++)
	{
		if (strcmp(indi[i]->father, "0") == 0 &&
		    strcmp(indi[i]->mother, "0") == 0)
		{
			count = 0;
			while (count < 2)
			{
				if (type == 1)   //haplotype information
				{
					randnum = randomnum();
					for (j=0; j<haplo_or_allele; j++)
						if(randnum < possibilities[j])
						{	tmp = j;  break;      }    
					for (t=0; t<loci; t++)
					    indi[i]->geno->geno_chars[index+t][count] = haplotype[tmp][t];
				}
				if (type == 2)   //genotype information
				{
					randnum = randomnum();
					for (j=0; j<haplo_or_allele; j++)
						if (randnum < allele_freq[j])
						{	tmp = j+1;  break;    }   
					indi[i]->geno->geno_chars[index][count] = tmp; 
				}
				count ++;
			}
			indi[i]->filled = 1;
		}		
	}//end of assigning genotypes to founders
	index = index + loci;    //track the markers 'simped' has processed for founders
}	

/************************************************/
/*                                              */
/*              assignGenotype                  */
/*                                              */
/************************************************/
void assignGenotype(int thisPerson, int parent, char c)
{
	int j, k, side;
	float randnum;

	if (c == 'f')	k = 0;		//comes from father
	if (c == 'm')	k = 1; 	    	//comes from mother
	
	//for 1st locus, to determine which haplotype should be select from father or mother
	randnum = randomnum();
	if (randnum < 0.5 ) {
		indi[thisPerson]->geno->geno_chars[0][k] = indi[parent]->geno->geno_chars[0][0];
		side = 0;
	} else {
		indi[thisPerson]->geno->geno_chars[0][k] = indi[parent]->geno->geno_chars[0][1];
		side = 1;
	}

	//select haplotypes for all the other loci according to the recombination fractions
	for (j=1; j<totalLoci; j++)
	{
		randnum = randomnum();

		//if recombination ocurrs, pick the next haplotype in the other side;
		//otherwise keep picking the next haplotype from the same side
        	if (randnum <= recomb_fractions[j-1]) 
		   side = 1 - side;

		indi[thisPerson]->geno->geno_chars[j][k] = indi[parent]->geno->geno_chars[j][side];
	}
}

/************************************************/
/*						*/
/*		readPed				*/
/*						*/
/************************************************/
void readPed(char* pedfile)
{ 
	int  i;
	char currentFamily[20];
	FILE *fp;

	i = 0;	
	totalPed = 0;
	indi[i] = (struct person *) malloc(sizeof(struct person)); //dynamic allocation

	if( (fp = fopen(pedfile, "r")) == NULL) 
        { 
                printf("ERROR: cannot open file \'%s\'\n\n", pedfile); 
                exit(0);
        }

	fscanf(fp, "%s", indi[i]->pedId);
	strcpy(currentFamily, indi[i]->pedId);
    	while (!feof(fp))	
	{	
		if (strcmp(currentFamily, indi[i]->pedId) == 0)  //individual from same family
		{
			membersOfFamily[totalPed]++;
		}  else  { //individual from a new family
			totalPed ++;
			membersOfFamily[totalPed] = 1;
			strcpy(currentFamily, indi[i]->pedId);
		}
		fscanf(fp, "%s%s%s%d%*[^\n]", indi[i]->id, indi[i]->father, indi[i]->mother, &indi[i]->sex);

		i++;
		indi[i] = (struct person *) malloc(sizeof(struct person));  //dynamic allocation
		fscanf(fp, "%s", indi[i]->pedId);
	}
	fclose(fp);
	
	totalInd = i;
	if (totalInd > MAX_IND)
	{
		printf("\nERROR: maximum number %d of individuals exceeded.\n\n", MAX_IND);
		exit(1);
	}
	totalPed += 1;
	if (totalPed > MAX_PED)
	{
		printf("\nERROR: maximum number %d of pedigrees exceeded.\n\n", MAX_PED);
		exit(1);
	}
}

/************************************************/
/*                                              */
/*              readParam                       */
/*                                              */
/************************************************/
void readParam(char* paramfile)
{
        int   i, j;
	int   pattern, lp, mod;
        float f;
        FILE  *fp;
        char  file1[50], file2[50];
	float tmp_array[MAX_LOCUS];
	char  c;

        fp = fopen(paramfile, "r");
                
        //read the 1st line of input, which are the file names
        fscanf(fp, "%s%s%*[^\n]", file1, file2);

        //read the 2nd line of input, which are the three random seeds
        fscanf(fp, "%d%d%d%*[^\n]", &seeds[0], &seeds[1], &seeds[2]);
                
        /*read in the 3rd line of input, which is # of columns between the sex column
        and the first marker genotype column. */
        fscanf(fp, "%d%*[^\n]", &extra_col);
        
        //read in the 4th line of input, which is the number of replicates
        fscanf(fp, "%d%*[^\n]", &replicates);

        //read in the 5th line of input 
        fscanf(fp, "%d%d%*[^\n]", &totalLoci, &iteration);
        if (totalLoci > MAX_LOCUS)
        {
                printf("\nERROR: maximum number %d of markers exceeded.\n\n", MAX_LOCUS);
                exit(1);
        }
        if ((totalLoci % iteration) != 0)
        {
                printf("ERROR in parameter file: the total number of loci divided by");
                printf(" the number of \ntimes the pattern should be repeated");
                printf(" must be a whole number.\n\n");
                exit(0);
        }

        //read in the 6th line of input
        fscanf(fp, "%d%*[^\n]", &inter_dist_type);
        
        /*read in the 7th line of input, which are a set of recombination fractions or inter-
        marker distances. Convert inter-marker distances to recombination fractions.*/
	fscanf(fp, "%d", &pattern);
        for (i=0; i<pattern; i++)
        {
                fscanf(fp, "%f", &f);
                if (inter_dist_type == 1)
                {
                        tmp_array[i] = f;
                } else if (inter_dist_type == 2) { //convert Kosambi dist. to theta
                        tmp_array[i] = 0.5 * (exp(4*f)-1) / (exp(4*f)+1);
                } else if (inter_dist_type == 3) { //convert Haldane dist. to theta
                        tmp_array[i] = 0.5 * (1- exp(-2 * f));
                }
        }
        fscanf(fp, "%*[^\n]");
	
	//get the recombination fractions according to the pattern
	mod = (totalLoci-1) % pattern;
	lp = (int)(totalLoci-1) / pattern;
	for (i=0; i<lp; i++)
		for (j=0; j<pattern; j++)
			recomb_fractions[i*pattern+j] = tmp_array[j];
	if (mod != 0)
		for (j=0; j<mod; j++)
			recomb_fractions[i*pattern+j] = tmp_array[j];

	offset = ftell(fp);	//record the current value of file-position indicator
	fclose(fp);
}

/************************************************/
/*                           		        */
/*              readMarkerInfo  	        */
/*                                              */
/************************************************/
void readMarkerInfo(char* paramfile)
{
	int   t1, t2;
	int   repeat_h, repeat_g;
	FILE  *fp;
	int   ii, jj, type;
	int   no_of_haplo, no_of_loci, no_of_allele;
	float sum, freq;
	float tmp_matrix[totalLoci][MAX_ALLELE];
	int   ary_no_of_allele[totalLoci];
	
	int  tmp_loci_of_set = 0;

	fp = fopen(paramfile, "r");

	//skip the information of the first seven lines.
	if (fseek(fp, offset, 0) != 0)
	{
		printf("ERROR: cannot move file pointer in parameter file.\n\n");
		exit(0);
	}

	//keep read in the following lines for genotype info. or haplotype info.
	fscanf(fp, "%d", &type);  //the 8th line of input parameter file

	while (!feof(fp))
	{
		if (type == 1) //prepare to read in a set of information for haplotypes
		{
			fscanf(fp, "%d%d%d%*[^\n]", &no_of_haplo, &no_of_loci, &repeat_h);
			tmp_loci_of_set += no_of_loci * repeat_h;	
		
		        if (no_of_haplo > MAX_HAPLOTYPES)
        		{
                		printf("\nERROR: maximum number %d of haplotypes exceeded.\n\n", MAX_HAPLOTYPES);
                		exit(1);
        		}

			//read in the frequency for each haplotype
			for (ii=0; ii<no_of_haplo; ii++)
			{
				fscanf(fp, "%f", &possibilities[ii]);
			}
			fscanf(fp, "%*[^\n]");
			
			for (ii=0; ii<no_of_haplo; ii++)
			{
				//make some adjustments to the array for possibilities
				possibilities[ii] = possibilities[ii-1] + possibilities[ii];
				
				//read in observed alleles for each haplotype
				for (jj=0; jj<no_of_loci; jj++)
				{
					fscanf(fp, "%d", &haplotype[ii][jj]);
				}
				fscanf(fp, "%*[^\n]");
			}

			for (ii=0; ii<repeat_h; ii++)
				buildFounderGeno(type, no_of_loci, no_of_haplo);
				
		} //end of type==1
		
		if (type == 2)  //prepare to read in a set of information for genotype
		{
			fscanf(fp, "%d%d%*[^\n]", &no_of_loci, &repeat_g);
			tmp_loci_of_set += no_of_loci * repeat_g;
			
			//read in the lines for allele frequencies and store the information
			for (ii=0; ii<no_of_loci; ii++)
			{
				fscanf(fp, "%d", &t1);
				ary_no_of_allele[ii] = t1;
			        if (ary_no_of_allele[ii] > MAX_ALLELE)
			        {
                			printf("\nERROR: maximum number %d of alleles exceeded.\n\n", MAX_ALLELE);
                			exit(1);
        			}
				sum = 0;
				for (jj=0; jj<t1-1; jj++)
				{
					fscanf(fp, "%f", &freq);
					sum += freq;
					tmp_matrix[ii][jj] = sum;   //store the allele freq in each line into a matrix
				}
				tmp_matrix[ii][jj] = 1;   
				fscanf(fp, "%*[^\n]");				
			}
        		for (ii=0; ii<repeat_g; ii++)
                		for (jj=0; jj<no_of_loci; jj++)
				{
					for (t2=0; t2<ary_no_of_allele[jj]; t2++) 
                        			allele_freq[t2] = tmp_matrix[jj][t2];
					buildFounderGeno(type, 1, ary_no_of_allele[jj]);
				}
		} //end of type == 2

		fscanf(fp, "%d", &type);
	}
	
	if (tmp_loci_of_set != totalLoci/iteration)
	{
		printf("ERROR in parameter file: the total number of loci on line 5 does not match\n");
		printf("the sum of the number of loci in proceeding lines.\n\n");
		exit(0);
	} 	 	
	fclose(fp);
}

/************************************************/
/*                                              */
/*              simulate                        */
/*                                              */
/************************************************/
void simulate()
{
        int ii, jj, k;
	int tmp, counter;
        int currentOne, currentPa, currentMa;

        tmp = 0;
        for (jj=0; jj<totalPed; jj++)
        {       //process the simulation family by family
                while (1) 
                {
                        counter = 0;
                        for (ii=0; ii<membersOfFamily[jj]; ii++)
                        {       //for every member in a pedigree
                                
                                //convert the index inside family to the index for all the individuals
                                currentOne = tmp + ii;
                                //check if the individual's genotypes are available
                                if (indi[currentOne]->filled == 1)    counter++;
                                else { //begin simulation  
					for (k=0; k<membersOfFamily[jj]; k++)
					{
					    if (strcmp(indi[tmp+k]->id, indi[currentOne]->father) == 0)
						currentPa = tmp + k;      //find father's real index
					    if (strcmp(indi[tmp+k]->id, indi[currentOne]->mother) == 0)
                                                currentMa = tmp + k;	  //find mother's real index
					}
					if (indi[currentPa]->filled == 1 && indi[currentMa]->filled == 1)
					{	
						assignGenotype(currentOne, currentPa, 'f');
						assignGenotype(currentOne, currentMa, 'm');
						indi[currentOne]->filled = 1;
 					}      
                                }
                        }
			//check to see if everybody of this family has genotypes assigned
			//if yes, break the infinite loop and process the next family
			if (counter == membersOfFamily[jj]) 	break; 
                }
                tmp += membersOfFamily[jj];	//set the start index for the next family
        }
}

/************************************************/
/*                                              */
/*              writeResult                     */
/*                                              */
/************************************************/
writeResult(char* pedfile, FILE *fp2)
{
	int i, j, lineNo;
	int  sex;
	char sped[20], currentFamily[20]; 
	char sperson[20], fa[20], ma[20];
	char sequence[totalLoci];
	char str[extra_col][5];
	char c;
	FILE *fp1;

        fp1 = fopen(pedfile, "r");
	
	fscanf(fp1, "%s", sped);
	strcpy(currentFamily, sped);
	lineNo = 0;
        while (!feof(fp1))
        {
		if (strcmp(currentFamily, sped) != 0)  //read in a new pedigree
		{
			familyNo ++;
			strcpy(currentFamily, sped); 
		} 
		fprintf(fp2, "%d  ", familyNo+1);
		fscanf(fp1, "%s%s%s%d", sperson, fa, ma, &sex);
		fprintf(fp2, "%s  %s  %s  %d", sperson, fa, ma, sex);
		for (i=0; i<extra_col; i++)
		{
			fscanf(fp1, "%s", str[i]);			
			fprintf(fp2, "  %s", str[i]);
		}
		fprintf(fp2, "  ");
		
		//read in a sequence of 0's and 1's, which indicates a marker should be simulated or not
		j = 0;
		while ((c=fgetc(fp1)) != '\n')
		{	
			if (c!=' ' && c!='\t')
			{
				sequence[j] = c;
				j++;
			}
		}
		if (j == 1)  //in pedigree file, a sigle 0 or 1 is used to denote all markers should be generated or not
			if (sequence[0] == '1')
			    for (i=0; i<totalLoci; i++)
				fprintf(fp2, " %d %d", indi[lineNo]->geno->geno_chars[i][0], indi[lineNo]->geno->geno_chars[i][1]);
			else //sequence[0]=='0'
			    for (i=0; i<totalLoci; i++)
				fprintf(fp2, " 0 0");
		else
			for (i=0; i<totalLoci; i++)
			    if (sequence[i] == '1')
				fprintf(fp2, " %d %d", indi[lineNo]->geno->geno_chars[i][0], indi[lineNo]->geno->geno_chars[i][1]);
			    else
				fprintf(fp2, " 0 0");

		fprintf(fp2, "\n");

		fscanf(fp1, "%s", sped);  //read in the next line
		lineNo ++;
	}
	fclose(fp1);
}

/************************************************/
/*                                              */
/*              main function                   */
/*                                              */
/************************************************/
main(int argc, char *argv[])
{
	char parameter_file[50], pedigree_file[50];	
	char output[50];
	FILE *fp;
	int i, j;

	printf("\n  ********************************************************\n");
	printf("  *                                                      *\n");
	printf("  *               Program   SimPed                       *\n");
	printf("  *                                                      *\n");
	printf("  ********************************************************\n");
	printf("\nUsage: simped parameter_file\n");
		
	if (argc > 2)
	{
		printf("ERROR: Too many command line arguments!\n\n");
		exit(0);
	} 
	if (argc == 2)   //get parameter file if on command line
		strcpy(parameter_file, argv[1]); 
	else {		 //file name is not on command line
		printf("\nParameter file -> ");
		scanf("%s", parameter_file);
	}

        if( (fp = fopen(parameter_file, "r")) == NULL)
        {
                printf("ERROR: cannot open file \'%s\'\n\n", parameter_file);
                exit(0);
        }                
        //read the 1st line of input, which are the file names
        fscanf(fp, "%s%s%*[^\n]", pedigree_file, output);
	fclose(fp);

	printf("\nConstants in effect:\n");
	printf("  Maximum number of pedigrees\t\t\t%d\n", MAX_PED);
	printf("  Maximum number of individuals\t\t\t%d\n", MAX_IND);
	printf("  Maximum number of loci\t\t\t%d\n", MAX_LOCUS);
	printf("  Maximum number of haplotypes\t\t\t%d\n", MAX_HAPLOTYPES);
	printf("  Maximum number of alleles at a locus\t\t%d\n", MAX_ALLELE);
	printf("\nInput pedigree file: %s\n", pedigree_file);
	printf("Output file: %s\n\n", output);

	if( (fp = fopen(output, "w")) == NULL)
	{
		printf("ERROR: cannot open file \'%s\'\n\n", output);
		exit(0);
	}
	
	//get the pedigree structure for each pedigree
	readPed(pedigree_file);
	//get the common parameters -- the first 7 lines in paramter file
	readParam(parameter_file);

	//dynamic memory allocation for storage of everybody's genotypes
	for (i=0; i<totalInd; i++)
	{
		indi[i]->geno = (struct genotype *) malloc(sizeof(struct genotype));
		if (indi[i]->geno == NULL)   
		{
			printf("ERROR: Cannot allocate memory\n\n");
			exit(1);
		}
	}

	//initialize the static variable, which is increased in function writeResults.
	familyNo = 0;
	for (i=0; i<replicates; i++)
	{
		//reset the static variable, which is increased in function 
		//buildFounderGeno called by readMarkerInfo.
		index = 0;

		//1. get the haplotype or/and genotype information from the parameter file, 
		//2. assign genotypes for all founders by using the random number generator
		for (j=1; j<=iteration; j++)
			readMarkerInfo(parameter_file);

		//determine genotypes for non-founders according to their parents' genotypes
		//and the recombination fractions
		simulate(); 

		//write the results of this replicate
		writeResult(pedigree_file, fp);

		//prepare for the next replicate
		familyNo ++;
		for (j=0; j<totalInd; j++)
			indi[j]->filled = 0;
	}
	fclose(fp);
}	


